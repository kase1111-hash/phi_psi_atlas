#!/usr/bin/env python3
"""
Test 1: pLDDT Confidence at Mechanical Singularity Hotspots
============================================================

Reads AlphaFold CIF files directly from your alphafold_cache and checks
whether the hotspot residue for each disease singularity falls in a
high-confidence (pLDDT >= 70) or low-confidence region.

If hotspots cluster in low-confidence disordered loops, the geometric
signal may be an AlphaFold prediction artifact rather than real biology.

Usage:
  python test1_plddt_validation.py \
    --singularities human_annotated_singularities.csv \
    --cache G:/Mechanical Taxonomy/alphafold_cache/human
    --output test1_plddt_report.txt

The cache directory should contain files named like: P04062.cif
(just UniProt ID + .cif extension)
"""

import os
import sys
import csv
import gzip
import argparse
import statistics
from collections import Counter


def extract_plddt_from_cif(filepath):
    """Extract per-residue pLDDT from an AlphaFold mmCIF file.
    
    Tries multiple strategies:
    1. _ma_qa_metric_local loop (AlphaFold v3/v4 CIF)
    2. _atom_site B_iso_or_equiv on CA atoms (fallback)
    
    Returns dict: {residue_number(int): plddt_value(float)}
    """
    opener = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'
    
    try:
        with opener(filepath, mode) as f:
            lines = f.readlines()
    except Exception as e:
        return None, f"Read error: {e}"
    
    # Strategy 1: _ma_qa_metric_local
    plddt = _parse_ma_qa_metric(lines)
    if plddt:
        return plddt, "ma_qa_metric_local"
    
    # Strategy 2: _atom_site B-factor on CA atoms
    plddt = _parse_atom_site_bfactor(lines)
    if plddt:
        return plddt, "atom_site_bfactor"
    
    return None, "no_plddt_found"


def _parse_ma_qa_metric(lines):
    """Parse _ma_qa_metric_local loop for pLDDT values."""
    plddt = {}
    in_loop = False
    col_names = []
    col_map = {}
    reading_data = False
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('_ma_qa_metric_local.'):
            in_loop = True
            reading_data = False
            col_name = line.split('.')[1].split()[0]
            col_map[col_name] = len(col_names)
            col_names.append(col_name)
            continue
        
        if in_loop and not line.startswith('_'):
            if not line or line.startswith('#') or line.startswith('loop_') or line.startswith('data_'):
                if plddt:
                    break
                in_loop = False
                continue
            
            reading_data = True
            parts = line.split()
            if len(parts) >= len(col_names):
                try:
                    # Find residue number column
                    seq_col = None
                    for key in ['label_seq_id', 'seq_id']:
                        if key in col_map:
                            seq_col = col_map[key]
                            break
                    
                    # Find pLDDT value column
                    val_col = None
                    for key in ['metric_value']:
                        if key in col_map:
                            val_col = col_map[key]
                            break
                    
                    if seq_col is not None and val_col is not None:
                        res_num = int(parts[seq_col])
                        val = float(parts[val_col])
                        plddt[res_num] = val
                except (ValueError, IndexError):
                    pass
    
    return plddt if plddt else None


def _parse_atom_site_bfactor(lines):
    """Parse _atom_site loop, extract B-factor from CA atoms as pLDDT."""
    plddt = {}
    in_loop = False
    col_names = []
    col_map = {}
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('_atom_site.'):
            in_loop = True
            col_name = line.split('.')[1].split()[0]
            col_map[col_name] = len(col_names)
            col_names.append(col_name)
            continue
        
        if in_loop and not line.startswith('_'):
            if not line or line.startswith('#') or line.startswith('loop_') or line.startswith('data_'):
                if plddt:
                    break
                if line.startswith('#'):
                    in_loop = False
                continue
            
            parts = line.split()
            if len(parts) >= len(col_names):
                try:
                    # Only CA atoms
                    atom_col = col_map.get('label_atom_id', col_map.get('auth_atom_id'))
                    if atom_col is not None and parts[atom_col] == 'CA':
                        seq_col = col_map.get('label_seq_id', col_map.get('auth_seq_id'))
                        bfac_col = col_map.get('B_iso_or_equiv')
                        
                        if seq_col is not None and bfac_col is not None:
                            res_num = int(parts[seq_col])
                            val = float(parts[bfac_col])
                            plddt[res_num] = val
                except (ValueError, IndexError):
                    pass
    
    return plddt if plddt else None


def run_test1(singularities_csv, cache_dir, output_file):
    """Main Test 1 execution."""
    
    # Load singularities
    proteins = []
    with open(singularities_csv, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                proteins.append({
                    'uid': r['uid'].strip(),
                    'gene': r['gene'].strip(),
                    'hotspot': int(r['hotspot_residue']),
                    'n_residues': int(r['n_residues']),
                    'gini': float(r['gini']),
                    'concentration': float(r['concentration_ratio']),
                    'disease': r.get('disease_names', ''),
                    'disease_count': int(r.get('disease_count', 0)),
                })
            except (ValueError, KeyError) as e:
                pass
    
    total = len(proteins)
    disease_proteins = [p for p in proteins if p['disease_count'] > 0]
    
    print(f"Loaded {total} singularities ({len(disease_proteins)} with disease annotations)")
    print(f"Searching for structures in: {cache_dir}")
    print()
    
    # Process each protein
    results = []
    found = 0
    not_found = 0
    parse_fail = 0
    
    for p in proteins:
        # Look for CIF file — try multiple naming patterns
        uid = p['uid']
        candidates = [
            os.path.join(cache_dir, f"AF-{uid}-F1-model_v6.cif"),
            os.path.join(cache_dir, f"AF-{uid}-F1-model_v5.cif"),
            os.path.join(cache_dir, f"AF-{uid}-F1-model_v4.cif"),
            os.path.join(cache_dir, f"AF-{uid}-F1-model_v3.cif"),
            os.path.join(cache_dir, f"AF-{uid}-F1-model_v2.cif"),
            os.path.join(cache_dir, f"{uid}.cif"),
            os.path.join(cache_dir, f"AF-{uid}-F1-model_v6.cif.gz"),
            os.path.join(cache_dir, f"{uid}.cif.gz"),
        ]
        
        filepath = None
        for c in candidates:
            if os.path.isfile(c):
                filepath = c
                break
        
        if not filepath:
            not_found += 1
            results.append({**p, 'status': 'no_structure',
                          'hotspot_plddt': None, 'window_mean': None,
                          'protein_mean': None, 'category': 'MISSING'})
            continue
        
        found += 1
        
        # Extract pLDDT
        plddt, method = extract_plddt_from_cif(filepath)
        
        if not plddt:
            parse_fail += 1
            results.append({**p, 'status': f'parse_fail ({method})',
                          'hotspot_plddt': None, 'window_mean': None,
                          'protein_mean': None, 'category': 'PARSE_FAIL'})
            continue
        
        # Get hotspot pLDDT
        hotspot_val = plddt.get(p['hotspot'])
        
        # Try ±1, ±2 if exact residue not found (numbering offsets)
        if hotspot_val is None:
            for offset in [1, -1, 2, -2]:
                hotspot_val = plddt.get(p['hotspot'] + offset)
                if hotspot_val is not None:
                    break
        
        # Window ±5 residues
        window = [plddt[r] for r in range(p['hotspot'] - 5, p['hotspot'] + 6) if r in plddt]
        window_mean = statistics.mean(window) if window else None
        
        # Whole protein stats
        all_vals = list(plddt.values())
        protein_mean = statistics.mean(all_vals)
        protein_min = min(all_vals)
        
        # Categorize hotspot confidence
        if hotspot_val is not None:
            if hotspot_val >= 90:
                cat = 'VERY_HIGH'
            elif hotspot_val >= 70:
                cat = 'CONFIDENT'
            elif hotspot_val >= 50:
                cat = 'LOW'
            else:
                cat = 'VERY_LOW'
        else:
            cat = 'RESIDUE_NOT_FOUND'
        
        results.append({
            **p,
            'status': 'ok',
            'hotspot_plddt': hotspot_val,
            'window_mean': window_mean,
            'protein_mean': protein_mean,
            'protein_min': protein_min,
            'n_plddt_residues': len(plddt),
            'category': cat,
            'parse_method': method,
        })
        
        if found % 50 == 0:
            print(f"  Processed {found} structures...")
    
    print(f"\nDone. Found={found}, Missing={not_found}, ParseFail={parse_fail}")
    
    # =====================================================
    # GENERATE REPORT
    # =====================================================
    R = []
    R.append("=" * 72)
    R.append("  TEST 1: pLDDT CONFIDENCE AT MECHANICAL SINGULARITY HOTSPOTS")
    R.append("  Question: Are singularity hotspots AlphaFold artifacts?")
    R.append("=" * 72)
    R.append("")
    R.append(f"  Input: {total} singularities ({len(disease_proteins)} disease-associated)")
    R.append(f"  Structures found:      {found}")
    R.append(f"  Structures missing:    {not_found}")
    R.append(f"  Parse failures:        {parse_fail}")
    R.append("")
    
    # Filter to successfully analyzed
    ok = [r for r in results if r['status'] == 'ok' and r['hotspot_plddt'] is not None]
    ok_disease = [r for r in ok if r['disease_count'] > 0]
    
    if not ok:
        R.append("  ** NO pLDDT DATA EXTRACTED — CHECK CACHE PATH AND FILE FORMAT **")
        report_text = "\n".join(R)
        print(report_text)
        with open(output_file, 'w') as f:
            f.write(report_text)
        return
    
    # ── ALL SINGULARITIES ──
    R.append(f"  ═══ ALL SINGULARITIES (n={len(ok)}) ═══")
    R.append("")
    _report_category(R, ok, "  ")
    
    # ── DISEASE SINGULARITIES ONLY ──
    if ok_disease:
        R.append("")
        R.append(f"  ═══ DISEASE SINGULARITIES ONLY (n={len(ok_disease)}) ═══")
        R.append("")
        _report_category(R, ok_disease, "  ")
    
    # ── HOTSPOT vs PROTEIN-WIDE COMPARISON ──
    R.append("")
    R.append("  ═══ HOTSPOT vs PROTEIN-WIDE pLDDT ═══")
    R.append("")
    diffs = [r['hotspot_plddt'] - r['protein_mean'] for r in ok]
    avg_diff = statistics.mean(diffs)
    R.append(f"  Average (hotspot - protein_mean): {avg_diff:+.1f}")
    if avg_diff >= 0:
        R.append(f"  Hotspots are in EQUAL or BETTER predicted regions than average")
    elif avg_diff > -10:
        R.append(f"  Hotspots are slightly lower confidence than average (small effect)")
    else:
        R.append(f"  WARNING: Hotspots are in notably lower confidence regions")
    
    # ── VERDICT ──
    R.append("")
    R.append("  " + "=" * 50)
    
    passed = len([r for r in ok if r['category'] in ('VERY_HIGH', 'CONFIDENT')])
    pass_rate = passed / len(ok)
    
    if pass_rate >= 0.90:
        R.append(f"  VERDICT: TEST 1 PASSED")
        R.append(f"  {passed}/{len(ok)} ({100*pass_rate:.1f}%) hotspots have pLDDT >= 70")
        R.append(f"  Mechanical singularities are NOT prediction artifacts.")
        R.append(f"  Safe to proceed.")
    elif pass_rate >= 0.70:
        R.append(f"  VERDICT: TEST 1 PARTIALLY PASSED")
        R.append(f"  {passed}/{len(ok)} ({100*pass_rate:.1f}%) hotspots have pLDDT >= 70")
        R.append(f"  Most hotspots are real, but filter low-confidence cases.")
    else:
        R.append(f"  VERDICT: TEST 1 FAILED")
        R.append(f"  Only {passed}/{len(ok)} ({100*pass_rate:.1f}%) hotspots have pLDDT >= 70")
        R.append(f"  Singularities may be AlphaFold artifacts. Investigate further.")
    R.append("  " + "=" * 50)
    
    # ── LOW CONFIDENCE HOTSPOTS ──
    low_conf = [r for r in ok if r['category'] in ('LOW', 'VERY_LOW')]
    if low_conf:
        R.append("")
        R.append(f"  LOW-CONFIDENCE HOTSPOTS ({len(low_conf)} proteins):")
        R.append(f"  {'Gene':<12s} {'UID':<12s} {'Hotspot':>7s} {'pLDDT':>6s} {'ProtAvg':>7s} {'Disease'}")
        R.append(f"  {'-'*12} {'-'*12} {'-'*7} {'-'*6} {'-'*7} {'-'*30}")
        for r in sorted(low_conf, key=lambda x: x['hotspot_plddt']):
            R.append(f"  {r['gene']:<12s} {r['uid']:<12s} {r['hotspot']:>7d} "
                     f"{r['hotspot_plddt']:>6.1f} {r['protein_mean']:>7.1f} "
                     f"{r['disease'][:40]}")
    
    # ── VALIDATION TRIO ──
    R.append("")
    R.append("  VALIDATION TRIO:")
    trio_uids = {'P04062': 'GBA1', 'Q9Y6N1': 'COX11', 'Q9BRN9': 'TM2D3'}
    for uid, gene in trio_uids.items():
        match = [r for r in results if r['uid'] == uid]
        if match:
            r = match[0]
            if r['hotspot_plddt'] is not None:
                R.append(f"  {gene:<10s} | res {r['hotspot']:>5d} | "
                        f"pLDDT {r['hotspot_plddt']:>5.1f} | {r['category']}")
            else:
                R.append(f"  {gene:<10s} | {r['status']}")
        else:
            R.append(f"  {gene:<10s} | not found in singularity list")
    
    # Print and save
    report_text = "\n".join(R)
    print("\n" + report_text)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(report_text)
    print(f"\nReport saved to {output_file}")
    
    # Save per-protein CSV
    csv_out = output_file.replace('.txt', '_results.csv')
    with open(csv_out, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['uid', 'gene', 'hotspot_residue', 'hotspot_plddt',
                     'window_plddt', 'protein_mean_plddt', 'plddt_category',
                     'gini', 'concentration_ratio', 'n_residues',
                     'disease_count', 'disease_names', 'status'])
        for r in results:
            w.writerow([
                r['uid'], r['gene'], r['hotspot'],
                f"{r['hotspot_plddt']:.1f}" if r.get('hotspot_plddt') is not None else '',
                f"{r['window_mean']:.1f}" if r.get('window_mean') is not None else '',
                f"{r['protein_mean']:.1f}" if r.get('protein_mean') is not None else '',
                r['category'], r['gini'], r['concentration'],
                r['n_residues'], r['disease_count'],
                r['disease'][:80], r['status']
            ])
    print(f"Per-protein CSV saved to {csv_out}")


def _report_category(R, ok, indent="  "):
    """Add category breakdown to report."""
    cats = Counter(r['category'] for r in ok)
    n = len(ok)
    
    very_high = cats.get('VERY_HIGH', 0)
    confident = cats.get('CONFIDENT', 0)
    low = cats.get('LOW', 0)
    very_low = cats.get('VERY_LOW', 0)
    
    R.append(f"{indent}pLDDT at hotspot residue:")
    R.append(f"{indent}  Very High (>= 90): {very_high:4d}  ({100*very_high/n:.1f}%)")
    R.append(f"{indent}  Confident (70-89): {confident:4d}  ({100*confident/n:.1f}%)")
    R.append(f"{indent}  Low (50-69):       {low:4d}  ({100*low/n:.1f}%)")
    R.append(f"{indent}  Very Low (< 50):   {very_low:4d}  ({100*very_low/n:.1f}%)")
    R.append(f"{indent}")
    
    hotspot_vals = [r['hotspot_plddt'] for r in ok]
    R.append(f"{indent}Hotspot pLDDT stats:")
    R.append(f"{indent}  Mean:   {statistics.mean(hotspot_vals):.1f}")
    R.append(f"{indent}  Median: {statistics.median(hotspot_vals):.1f}")
    R.append(f"{indent}  Min:    {min(hotspot_vals):.1f}")
    R.append(f"{indent}  Max:    {max(hotspot_vals):.1f}")
    if len(hotspot_vals) > 1:
        R.append(f"{indent}  StdDev: {statistics.stdev(hotspot_vals):.1f}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Test 1: pLDDT confidence at mechanical singularity hotspots")
    parser.add_argument('--singularities', required=True,
                       help='human_annotated_singularities.csv')
    parser.add_argument('--cache', required=True,
                       help='Path to alphafold_cache/human/ directory with {uid}.cif files')
    parser.add_argument('--output', default='test1_plddt_report.txt',
                       help='Output report filename (default: test1_plddt_report.txt)')
    
    args = parser.parse_args()
    
    if not os.path.isfile(args.singularities):
        print(f"Error: {args.singularities} not found")
        sys.exit(1)
    if not os.path.isdir(args.cache):
        print(f"Error: {args.cache} is not a directory")
        sys.exit(1)
    
    run_test1(args.singularities, args.cache, args.output)

#!/usr/bin/env python3
"""
Cross-Species Curvature Conservation Survey
=============================================
Runs the curvature survey across multiple proteomes and compares:
  1. Helix-energy relationship (is r = +0.63 universal?)
  2. Fold-class bending energy distributions
  3. Overall Σκ²/residue statistics

USAGE:
    # Point at your proteome database root directory.
    # Expects subdirectories per species, each containing .cif files:
    #   proteome_db/
    #     human/
    #       AF-P68871-F1-model_v6.cif
    #       ...
    #     mouse/
    #       AF-Q9QZ07-F1-model_v4.cif
    #       ...

    python cross_species_survey.py /path/to/proteome_db/ cross_species_results.csv

    # Limit per species (for testing):
    python cross_species_survey.py /path/to/proteome_db/ test.csv --limit 200

    # Process specific species only:
    python cross_species_survey.py /path/to/proteome_db/ results.csv --species human mouse zebrafish

    # Resume interrupted run:
    python cross_species_survey.py /path/to/proteome_db/ results.csv --resume

REQUIRES: numpy, pandas (for summary). Uses proteome_curvature_survey.py functions.
LICENSE: CC0 1.0
"""

import os
import sys
import csv
import time
import argparse
import numpy as np
from pathlib import Path


# ══════════════════════════════════════════════════════════════
# Import core functions from the survey script.
# If it's not importable, we embed the essentials.
# ══════════════════════════════════════════════════════════════

try:
    from proteome_curvature_survey import (
        parse_mmcif_backbone, compute_dihedrals, torus_curvature_standalone,
        assign_ss_from_dihedrals, dssp_assign, gini_coefficient,
        process_one_structure, extract_uniprot_id
    )
    HAVE_SURVEY = True
    print("Loaded functions from proteome_curvature_survey.py")
except ImportError:
    HAVE_SURVEY = False
    print("WARNING: proteome_curvature_survey.py not found on PYTHONPATH.")
    print("Make sure it's in the same directory or on your PYTHONPATH.")
    print("e.g.: set PYTHONPATH=G:\\Hemo-Mapping;%PYTHONPATH%")
    sys.exit(1)


def find_species_dirs(root_path):
    """Find subdirectories that contain .cif files."""
    root = Path(root_path)
    species = []
    
    # Check if root itself contains CIF files (single species)
    cifs_in_root = list(root.glob('*.cif'))
    if cifs_in_root:
        species.append(('root', root, len(cifs_in_root)))
        return species
    
    # Check subdirectories
    for d in sorted(root.iterdir()):
        if d.is_dir():
            cifs = list(d.glob('*.cif'))
            if cifs:
                species.append((d.name, d, len(cifs)))
    
    return species


def process_species(species_name, species_dir, output_csv, limit=0, use_dssp=True, resume_ids=None):
    """Process all CIF files for one species, append to CSV."""
    
    cif_files = sorted(species_dir.glob('*.cif'))
    if limit > 0:
        cif_files = cif_files[:limit]
    
    total = len(cif_files)
    n_success = 0
    n_fail = 0
    n_skip = 0
    t0 = time.time()
    
    # CSV fieldnames (same as proteome_curvature_survey.py + species column)
    fieldnames = [
        'species', 'uniprot_id', 'source_file', 'length', 'n_segments',
        'n_gauss_peaks', 'gauss_positions', 'total_kappa_sq', 'kappa_sq_per_res',
        'mean_kappa', 'max_kappa', 'std_kappa', 'kappa_gini',
        'ss_frac_H', 'ss_frac_E', 'ss_frac_C',
        'mean_plddt', 'min_plddt', 'plddt_below_70',
        'kappa_sq_q25', 'kappa_sq_q50', 'kappa_sq_q75', 'kappa_sq_q90',
        'top_peak_kappa_sq', 'top_peak_position', 'processing_time'
    ]
    
    # Check if file exists and has header
    file_exists = os.path.exists(output_csv)
    csvfile = open(output_csv, 'a', newline='')
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    if not file_exists:
        writer.writeheader()
    
    for i, cif_path in enumerate(cif_files):
        basename = os.path.basename(str(cif_path))
        uid = extract_uniprot_id(str(cif_path))
        
        # Resume support
        if resume_ids and f"{species_name}:{basename}" in resume_ids:
            n_skip += 1
            continue
        
        try:
            result = process_one_structure(str(cif_path), fast=True, use_dssp=use_dssp)
            
            if result is None:
                n_fail += 1
                continue
            
            result['species'] = species_name
            result['uniprot_id'] = uid
            result['source_file'] = basename
            
            writer.writerow(result)
            n_success += 1
            
            if n_success % 500 == 0 or n_success <= 5:
                elapsed = time.time() - t0
                rate = (n_success + n_skip) / max(elapsed, 1)
                remaining = total - i - 1
                eta = remaining / max(rate, 0.1)
                print(f"    [{i+1}/{total}] {basename}: L={result['length']}, "
                      f"Σκ²/r={result['kappa_sq_per_res']:.1f}, "
                      f"H={result['ss_frac_H']:.2f} "
                      f"({rate:.1f}/s, ETA {eta/60:.0f}min)")
        
        except Exception as e:
            n_fail += 1
            continue
    
    csvfile.flush()
    csvfile.close()
    
    elapsed = time.time() - t0
    print(f"    {species_name}: {n_success} processed, {n_fail} failed, "
          f"{n_skip} skipped, {elapsed/60:.1f}min")
    
    return n_success


def summarize_results(output_csv):
    """Print cross-species summary."""
    import pandas as pd
    
    df = pd.read_csv(output_csv)
    hq = df[(df['mean_plddt'] >= 70) & (df['length'] >= 50) & (df['length'] <= 2000)].copy()
    
    species_list = sorted(hq['species'].unique())
    
    print(f"\n{'='*80}")
    print(f"  CROSS-SPECIES CURVATURE SURVEY — {len(species_list)} species, {len(hq)} HQ proteins")
    print(f"{'='*80}")
    
    # Per-species summary
    print(f"\n  {'Species':>20s} {'n':>6s} {'Σκ²/r':>8s} {'±':>1s} {'SD':>7s} {'CV%':>6s} "
          f"{'med':>7s} {'H%':>5s} {'E%':>5s} {'r(H,κ²)':>8s} {'pLDDT':>6s}")
    print(f"  {'-'*20} {'-'*6} {'-'*8} {'-'*1} {'-'*7} {'-'*6} "
          f"{'-'*7} {'-'*5} {'-'*5} {'-'*8} {'-'*6}")
    
    all_r_values = []
    for sp in species_list:
        sub = hq[hq['species'] == sp]
        ksq = sub['kappa_sq_per_res']
        cv = ksq.std() / ksq.mean() * 100 if ksq.mean() > 0 else 0
        r_helix = sub['kappa_sq_per_res'].corr(sub['ss_frac_H'])
        all_r_values.append(r_helix)
        
        print(f"  {sp:>20s} {len(sub):>6d} {ksq.mean():>8.1f}   {ksq.std():>7.1f} {cv:>5.0f}% "
              f"{ksq.median():>7.1f} {sub['ss_frac_H'].mean()*100:>4.0f}% "
              f"{sub['ss_frac_E'].mean()*100:>4.0f}% {r_helix:>+7.3f} "
              f"{sub['mean_plddt'].mean():>5.1f}")
    
    # Test: is the helix-energy correlation universal?
    print(f"\n{'='*80}")
    print(f"  HELIX-ENERGY RELATIONSHIP ACROSS SPECIES")
    print(f"{'='*80}")
    print(f"  Mean r(helix, Σκ²/res):  {np.mean(all_r_values):+.3f}")
    print(f"  Std:                      {np.std(all_r_values):.3f}")
    print(f"  Min:                      {np.min(all_r_values):+.3f}")
    print(f"  Max:                      {np.max(all_r_values):+.3f}")
    print(f"  Human reference:          +0.631")
    
    if np.std(all_r_values) < 0.05 and np.mean(all_r_values) > 0.5:
        print(f"\n  CONCLUSION: Helix-energy relationship is UNIVERSAL (low variance, consistently positive)")
    elif np.mean(all_r_values) > 0.3:
        print(f"\n  CONCLUSION: Helix-energy relationship is GENERAL but variable across species")
    else:
        print(f"\n  CONCLUSION: Helix-energy relationship is NOT universal")
    
    # Fold-class comparison across species
    print(f"\n{'='*80}")
    print(f"  FOLD-CLASS BENDING ENERGY BY SPECIES")
    print(f"{'='*80}")
    
    hq['ss_class'] = np.where(hq['ss_frac_H'] > 0.40, 'α-rich',
                     np.where(hq['ss_frac_E'] > 0.30, 'β-rich',
                     np.where((hq['ss_frac_H'] > 0.15) & (hq['ss_frac_E'] > 0.15), 'α/β',
                     'coil')))
    
    print(f"\n  Mean Σκ²/residue by fold class:")
    print(f"  {'Species':>20s} {'α-rich':>10s} {'β-rich':>10s} {'α/β':>10s} {'coil':>10s}")
    for sp in species_list:
        sub = hq[hq['species'] == sp]
        vals = []
        for fc in ['α-rich', 'β-rich', 'α/β', 'coil']:
            fsub = sub[sub['ss_class'] == fc]
            if len(fsub) > 10:
                vals.append(f"{fsub['kappa_sq_per_res'].mean():>8.0f}")
            else:
                vals.append(f"{'n/a':>8s}")
        print(f"  {sp:>20s} {''.join(v + '  ' for v in vals)}")
    
    # Overall regression: does SS explain similar variance across species?
    print(f"\n  SS + Length → Σκ²/res regression R² by species:")
    for sp in species_list:
        sub = hq[hq['species'] == sp]
        if len(sub) < 50:
            continue
        X = sub[['length', 'ss_frac_H', 'ss_frac_E', 'ss_frac_C']].values
        y = sub['kappa_sq_per_res'].values
        X_aug = np.column_stack([np.ones(len(X)), X])
        try:
            beta, _, _, _ = np.linalg.lstsq(X_aug, y, rcond=None)
            y_pred = X_aug @ beta
            r2 = 1 - np.sum((y - y_pred)**2) / np.sum((y - y.mean())**2)
            print(f"    {sp:>20s}: R² = {r2:.3f} (n={len(sub)})")
        except:
            pass
    
    print(f"\n  Human reference: R² = 0.402")


def main():
    parser = argparse.ArgumentParser(description='Cross-species curvature survey')
    parser.add_argument('proteome_dir', help='Root directory containing species subdirectories with .cif files')
    parser.add_argument('output', help='Output CSV path')
    parser.add_argument('--limit', '-n', type=int, default=0,
                        help='Max files per species (0 = all)')
    parser.add_argument('--species', nargs='*', default=None,
                        help='Process only these species (directory names)')
    parser.add_argument('--no-dssp', action='store_true',
                        help='Skip DSSP, use dihedral SS only (faster)')
    parser.add_argument('--resume', action='store_true',
                        help='Skip already-processed entries')
    args = parser.parse_args()
    
    # Find species
    all_species = find_species_dirs(args.proteome_dir)
    if not all_species:
        print(f"ERROR: No .cif files found in {args.proteome_dir} or its subdirectories")
        sys.exit(1)
    
    # Filter species if requested
    if args.species:
        all_species = [(n, d, c) for n, d, c in all_species if n in args.species]
    
    print(f"Found {len(all_species)} species:")
    for name, path, count in all_species:
        print(f"  {name:>25s}: {count:>6d} CIF files")
    
    total_files = sum(min(c, args.limit) if args.limit > 0 else c for _, _, c in all_species)
    print(f"\nTotal: {total_files} files to process")
    
    # Resume support
    resume_ids = set()
    if args.resume and os.path.exists(args.output):
        import csv as csvmod
        with open(args.output, 'r') as f:
            for row in csvmod.DictReader(f):
                resume_ids.add(f"{row.get('species','')}:{row.get('source_file','')}")
        print(f"Resume: {len(resume_ids)} entries already processed")
    
    # Process each species
    t_total = time.time()
    for name, path, count in all_species:
        n_to_process = min(count, args.limit) if args.limit > 0 else count
        print(f"\n  Processing {name} ({n_to_process} files)...")
        
        process_species(
            species_name=name,
            species_dir=path,
            output_csv=args.output,
            limit=args.limit,
            use_dssp=not args.no_dssp,
            resume_ids=resume_ids if args.resume else None
        )
    
    elapsed = time.time() - t_total
    print(f"\n{'='*60}")
    print(f"  ALL SPECIES COMPLETE ({elapsed/60:.1f} min total)")
    print(f"{'='*60}")
    
    # Summary
    try:
        summarize_results(args.output)
    except Exception as e:
        print(f"Summary failed: {e}")
        print("Run summary manually with: python -c \"import pandas; ...\"")


if __name__ == '__main__':
    main()

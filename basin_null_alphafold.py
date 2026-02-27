#!/usr/bin/env python3
"""
BASIN-CONTROL NULL ANALYSIS — AlphaFold Edition

Runs the Markov-walk basin-control null on all cached AlphaFold structures.
This is the decisive test: are helices smoother than basin-constrained
random walks through the same Ramachandran region?

Usage:
    python basin_null_alphafold.py                  # Run on all cached structures
    python basin_null_alphafold.py --reps 10        # 10 null reps per segment (default: 5)
    python basin_null_alphafold.py --quick           # Fast mode: 2 reps, first 200 proteins

The script:
  1. Loads all cached AlphaFold structures (mmCIF/PDB from data/alphafold/)
  2. Extracts dihedrals and assigns SS
  3. Builds per-SS step pools (Δφ, Δψ) from real data
  4. For each real segment, generates N Markov-walk surrogates
  5. Classifies real and null segments identically
  6. Compares: helix non-constant fraction (real) vs (null)
  7. Computes bootstrap CI, LOO, effect sizes
  8. Saves results to results/basin_null_results.json

This answers the central question:
  Does the curvature suppression finding from n=23 hold at proteome scale?
"""

import sys
import os
import json
import time
import math
import numpy as np
from pathlib import Path
from collections import defaultdict

# Import the pipeline
sys.path.insert(0, str(Path(__file__).parent))
from cornu_spirals_on_T2 import (
    extract_dihedrals_mmcif, extract_dihedrals_manual,
    dssp_assign_pure_numpy, assign_ss_from_dihedrals,
    torus_curvature, classify_spiral,
    build_superpotential_grid,
    ALPHAFOLD_DIR, ALPHAFOLD_CACHE_DIR, RESULTS_DIR,
)

# ─── Configuration ───
N_REPS = 5          # null replicates per segment
SEED = 42
MIN_SEG_LEN = 4

MACRO_MAP = {
    'geodesic': 'constant_k', 'circular_arc': 'constant_k',
    'sinusoidal': 'oscillatory', 'damped_oscillation': 'oscillatory',
    'damped_osc': 'oscillatory',
    'linear': 'monotone', 'clothoid': 'monotone', 'fermat': 'monotone',
    'exponential': 'monotone',
    'sigmoid': 'transition', 'step': 'transition',
    'quadratic': 'polynomial',
    'gauss_peak': 'localized',
    'too_short': 'too_short',
}


def load_structure(fpath):
    """Extract dihedrals from a cached structure file."""
    fpath = str(fpath)
    if fpath.endswith('.pdb'):
        return extract_dihedrals_manual(fpath, chain_id='A')
    else:
        return extract_dihedrals_mmcif(fpath)


def segment_protein(phi, psi, backbone_coords, has_O):
    """Segment a protein and classify each segment. Returns segments list and SS array."""
    valid = ~(np.isnan(phi) | np.isnan(psi))
    phi_v = phi[valid]
    psi_v = psi[valid]
    bb_v = backbone_coords[valid]
    ho_v = has_O[valid]
    
    if len(phi_v) < 6:
        return [], None, phi_v, psi_v
    
    # SS assignment
    ss = None
    if np.all(ho_v):
        try:
            ss = dssp_assign_pure_numpy(bb_v)
        except Exception:
            pass
    if ss is None:
        ss = assign_ss_from_dihedrals(phi_v, psi_v)
    
    # Curvature
    kappa, s, _ = torus_curvature(phi_v, psi_v)
    if len(kappa) < 4:
        return [], ss, phi_v, psi_v
    
    # Segment by SS boundaries
    ss_changes = np.where(ss[1:] != ss[:-1])[0] + 1
    boundaries = np.concatenate([[0], ss_changes, [len(ss)]])
    
    segments = []
    for i in range(len(boundaries) - 1):
        start = int(boundaries[i])
        end = int(boundaries[i + 1]) - 1
        seg_ss = ss[start]
        seg_len = end - start + 1
        
        if seg_len < 2:
            continue
        
        seg_k = kappa[max(0, start):min(end + 1, len(kappa))]
        seg_s = s[max(0, start):min(end + 1, len(s))]
        
        if len(seg_k) < 2:
            continue
        
        seg_s_rel = seg_s - seg_s[0]
        
        # Classify
        if seg_len >= MIN_SEG_LEN and len(seg_k) >= MIN_SEG_LEN:
            slope, r_sq, class_name, fits, osc_info = classify_spiral(seg_k, seg_s_rel)
        else:
            class_name = 'too_short'
        
        segments.append({
            'start': start, 'end': end, 'length': seg_len,
            'ss_type': seg_ss, 'class': class_name,
            'macro': MACRO_MAP.get(class_name, 'constant_k'),
        })
    
    return segments, ss, phi_v, psi_v


def build_step_pools(all_protein_data):
    """Build per-SS (Δφ, Δψ) step pools and start positions from all proteins."""
    step_pools = {'H': [], 'E': [], 'C': []}
    start_pools = {'H': [], 'E': [], 'C': []}
    
    for phi_v, psi_v, segments in all_protein_data:
        for seg in segments:
            if seg['class'] == 'too_short' or seg['length'] < MIN_SEG_LEN:
                continue
            ss = seg['ss_type']
            if ss not in step_pools:
                continue
            
            s, e = seg['start'], seg['end']
            if e >= len(phi_v):
                continue
            
            seg_phi = phi_v[s:e+1]
            seg_psi = psi_v[s:e+1]
            
            # Periodic differences
            dphi = np.arctan2(np.sin(np.diff(seg_phi)), np.cos(np.diff(seg_phi)))
            dpsi = np.arctan2(np.sin(np.diff(seg_psi)), np.cos(np.diff(seg_psi)))
            steps = np.column_stack([dphi, dpsi])
            
            step_pools[ss].append(steps)
            start_pools[ss].append((seg_phi[0], seg_psi[0]))
    
    # Stack into arrays
    pooled = {}
    starts = {}
    for ss in ['H', 'E', 'C']:
        if step_pools[ss]:
            pooled[ss] = np.vstack(step_pools[ss])
            starts[ss] = start_pools[ss]
    
    return pooled, starts


def run_basin_null(all_protein_data, step_pools, start_pools, n_reps=5, seed=42):
    """Generate Markov-walk null segments and classify them.
    
    For each real segment, generate n_reps surrogates by:
      1. Pick a random start from the same SS pool
      2. Random-walk using steps drawn from the SS step pool
      3. Compute torus curvature and classify
    
    Returns: real_classes, null_classes (lists of (ss_type, macro_class))
    """
    rng = np.random.RandomState(seed)
    
    real_classes = []
    null_classes = []
    
    total_segs = sum(
        1 for _, _, segs in all_protein_data
        for seg in segs
        if seg['class'] != 'too_short' and seg['length'] >= MIN_SEG_LEN
        and seg['ss_type'] in step_pools
    )
    
    done = 0
    t0 = time.time()
    
    for phi_v, psi_v, segments in all_protein_data:
        for seg in segments:
            if seg['class'] == 'too_short' or seg['length'] < MIN_SEG_LEN:
                continue
            ss = seg['ss_type']
            if ss not in step_pools or len(step_pools[ss]) < 5:
                continue
            
            # Record real classification
            real_classes.append((ss, seg['macro']))
            
            pool = step_pools[ss]
            n_starts = len(start_pools[ss])
            seg_len = seg['length']
            
            # Generate null surrogates
            for rep in range(n_reps):
                si = rng.randint(n_starts)
                phi0, psi0 = start_pools[ss][si]
                
                step_indices = rng.randint(0, len(pool), size=seg_len - 1)
                steps = pool[step_indices]
                
                phi_path = np.empty(seg_len)
                psi_path = np.empty(seg_len)
                phi_path[0] = phi0
                psi_path[0] = psi0
                np.cumsum(steps[:, 0], out=phi_path[1:])
                phi_path[1:] += phi0
                np.cumsum(steps[:, 1], out=psi_path[1:])
                psi_path[1:] += psi0
                
                kappa_null, s_null, _ = torus_curvature(phi_path, psi_path)
                if len(kappa_null) >= MIN_SEG_LEN:
                    s_rel = s_null[:len(kappa_null)] - s_null[0] if len(s_null) > len(kappa_null) else \
                            np.arange(len(kappa_null), dtype=float)
                    _, _, cls_null, _, _ = classify_spiral(kappa_null, s_rel)
                    mc_null = MACRO_MAP.get(cls_null, 'constant_k')
                    null_classes.append((ss, mc_null))
            
            done += 1
            if done % 200 == 0:
                elapsed = time.time() - t0
                rate = done / elapsed if elapsed > 0 else 0
                eta = (total_segs - done) / rate if rate > 0 else 0
                print(f"    Basin null: {done}/{total_segs} segments "
                      f"({done/total_segs:.0%}, {rate:.0f} seg/s, ETA {eta/60:.0f} min)")
    
    return real_classes, null_classes


def compute_statistics(real_classes, null_classes, all_protein_data):
    """Compute per-SS comparison statistics: chi2, effect sizes, bootstrap, LOO."""
    from scipy import stats as sp_stats
    
    results = {}
    all_macros = sorted(set(mc for _, mc in real_classes + null_classes))
    
    for ss in ['H', 'E', 'C']:
        real_ss = [(s, mc) for s, mc in real_classes if s == ss]
        null_ss = [(s, mc) for s, mc in null_classes if s == ss]
        
        if not real_ss or not null_ss:
            continue
        
        real_total = len(real_ss)
        null_total = len(null_ss)
        
        # Count macro classes
        real_counts = {mc: sum(1 for _, m in real_ss if m == mc) for mc in all_macros}
        null_counts = {mc: sum(1 for _, m in null_ss if m == mc) for mc in all_macros}
        
        # Non-constant fractions
        real_nc = sum(1 for _, mc in real_ss if mc != 'constant_k')
        null_nc = sum(1 for _, mc in null_ss if mc != 'constant_k')
        real_nc_rate = real_nc / real_total
        null_nc_rate = null_nc / null_total
        
        # Chi-squared (full distribution)
        real_arr = np.array([real_counts.get(mc, 0) for mc in all_macros])
        null_arr = np.array([null_counts.get(mc, 0) for mc in all_macros])
        null_expected = null_arr * (real_total / null_total)
        mask = null_expected > 0
        if mask.sum() >= 2:
            chi2, p_val = sp_stats.chisquare(real_arr[mask], f_exp=null_expected[mask])
        else:
            chi2, p_val = 0.0, 1.0
        
        # Binary test (constant vs non-constant)
        a = real_total - real_nc   # real constant
        b = real_nc                # real non-constant
        c = null_total - null_nc   # null constant
        d = null_nc                # null non-constant
        
        table = np.array([[a, b], [c, d]])
        chi2_bin, p_bin = sp_stats.chi2_contingency(table)[:2]
        
        # Effect sizes
        arr = abs_risk_reduction = null_nc_rate - real_nc_rate
        odds_ratio = (a * d) / max(b * c, 1)
        n_total = table.sum()
        cramers_v = float(np.sqrt(chi2_bin / max(n_total, 1)))
        
        results[ss] = {
            'real_total': real_total,
            'null_total': null_total,
            'real_nonconst_frac': round(real_nc_rate, 4),
            'null_nonconst_frac': round(null_nc_rate, 4),
            'real_const_frac': round(1 - real_nc_rate, 4),
            'null_const_frac': round(1 - null_nc_rate, 4),
            'chi2_full': round(float(chi2), 2),
            'p_full': float(p_val),
            'chi2_binary': round(float(chi2_bin), 2),
            'p_binary': float(p_bin),
            'abs_risk_reduction': round(float(abs_risk_reduction), 4),
            'odds_ratio': round(float(odds_ratio), 3),
            'cramers_v': round(float(cramers_v), 4),
            'real_counts': {k: v for k, v in real_counts.items() if v > 0},
            'null_fracs': {k: round(v/null_total, 4) for k, v in null_counts.items() if v > 0},
            'direction': 'suppressed' if real_nc_rate < null_nc_rate else 'enriched',
        }
    
    # ── Bootstrap + LOO for helix (the key claim) ──
    if 'H' in results:
        print("\n    Helix robustness: bootstrap (1000 reps) + LOO...")
        
        # Map helix segments to proteins
        helix_by_protein = defaultdict(list)
        for idx, (phi_v, psi_v, segments) in enumerate(all_protein_data):
            for seg in segments:
                if seg['ss_type'] == 'H' and seg['class'] != 'too_short' and seg['length'] >= MIN_SEG_LEN:
                    helix_by_protein[idx].append(seg['macro'])
        
        protein_ids = list(helix_by_protein.keys())
        n_proteins = len(protein_ids)
        null_nc_rate = results['H']['null_nonconst_frac']
        
        if n_proteins >= 3:
            rng_boot = np.random.RandomState(999)
            
            # Bootstrap
            boot_fracs = []
            for _ in range(1000):
                boot_ids = rng_boot.choice(protein_ids, size=n_proteins, replace=True)
                boot_segs = []
                for pid in boot_ids:
                    boot_segs.extend(helix_by_protein[pid])
                if boot_segs:
                    nc = sum(1 for mc in boot_segs if mc != 'constant_k')
                    boot_fracs.append(nc / len(boot_segs))
            
            boot_arr = np.array(boot_fracs)
            boot_exceeds = float(np.mean(boot_arr >= null_nc_rate))
            
            results['H']['bootstrap'] = {
                'n_reps': 1000,
                'n_proteins': n_proteins,
                'mean_nonconst': round(float(np.mean(boot_arr)), 4),
                'ci_95_low': round(float(np.percentile(boot_arr, 2.5)), 4),
                'ci_95_high': round(float(np.percentile(boot_arr, 97.5)), 4),
                'frac_exceeds_null': round(boot_exceeds, 4),
            }
            
            # LOO
            loo_fracs = []
            for leave_out in protein_ids:
                loo_segs = []
                for pid in protein_ids:
                    if pid != leave_out:
                        loo_segs.extend(helix_by_protein[pid])
                if loo_segs:
                    nc = sum(1 for mc in loo_segs if mc != 'constant_k')
                    loo_fracs.append(nc / len(loo_segs))
            
            loo_arr = np.array(loo_fracs)
            results['H']['leave_one_out'] = {
                'n_proteins': n_proteins,
                'mean_nonconst': round(float(np.mean(loo_arr)), 4),
                'min_nonconst': round(float(np.min(loo_arr)), 4),
                'max_nonconst': round(float(np.max(loo_arr)), 4),
                'all_below_null': bool(np.max(loo_arr) < null_nc_rate),
            }
    
    return results


def print_results(results):
    """Print formatted basin null results."""
    print()
    print("=" * 72)
    print("  BASIN-CONTROL NULL RESULTS")
    print("=" * 72)
    
    for ss in ['H', 'E', 'C']:
        if ss not in results:
            continue
        r = results[ss]
        ss_name = {'H': 'HELIX', 'E': 'STRAND', 'C': 'COIL'}[ss]
        
        print(f"\n  ── {ss_name} ──")
        print(f"  Real segments:  {r['real_total']}")
        print(f"  Null segments:  {r['null_total']}")
        print(f"  Real non-const: {r['real_nonconst_frac']:.1%}")
        print(f"  Null non-const: {r['null_nonconst_frac']:.1%}")
        print(f"  Direction:      {r['direction'].upper()}")
        print(f"  ΔAR (null-real): {r['abs_risk_reduction']:+.1%}")
        print(f"  Odds ratio:     {r['odds_ratio']:.2f}")
        print(f"  χ² (binary):    {r['chi2_binary']:.2f}")
        print(f"  p (binary):     {r['p_binary']:.2e}")
        print(f"  Cramér's V:     {r['cramers_v']:.3f}")
        
        # Distribution
        print(f"\n  Macro-class distribution:")
        print(f"    {'Class':18s} {'Real':>8s} {'Null':>8s}")
        all_macros = sorted(set(list(r['real_counts'].keys()) + list(r['null_fracs'].keys())))
        for mc in all_macros:
            rc = r['real_counts'].get(mc, 0)
            nf = r['null_fracs'].get(mc, 0)
            rc_pct = rc / r['real_total'] if r['real_total'] > 0 else 0
            print(f"    {mc:18s} {rc_pct:7.1%} {nf:7.1%}")
        
        # Bootstrap / LOO
        if 'bootstrap' in r:
            b = r['bootstrap']
            print(f"\n  Bootstrap (1000 reps, {b['n_proteins']} proteins):")
            print(f"    Mean non-const: {b['mean_nonconst']:.1%}")
            print(f"    95% CI:         [{b['ci_95_low']:.1%}, {b['ci_95_high']:.1%}]")
            print(f"    Null rate:      {r['null_nonconst_frac']:.1%}")
            print(f"    Frac ≥ null:    {b['frac_exceeds_null']:.1%}")
        
        if 'leave_one_out' in r:
            l = r['leave_one_out']
            print(f"\n  Leave-one-out ({l['n_proteins']} proteins):")
            print(f"    Range: [{l['min_nonconst']:.1%}, {l['max_nonconst']:.1%}]")
            print(f"    All below null: {l['all_below_null']}")
    
    # ── THE VERDICT ──
    print()
    print("=" * 72)
    if 'H' in results:
        r = results['H']
        if r['direction'] == 'suppressed' and r['p_binary'] < 0.05:
            print("  VERDICT: Helix curvature SUPPRESSION confirmed at proteome scale")
            print(f"           (p = {r['p_binary']:.2e}, ARR = {r['abs_risk_reduction']:+.1%})")
            if 'bootstrap' in r:
                print(f"           Bootstrap {r['bootstrap']['frac_exceeds_null']:.0%} exceed null")
        elif r['direction'] == 'suppressed' and r['p_binary'] >= 0.05:
            print("  VERDICT: Helix curvature suppression TREND but not significant")
            print(f"           (p = {r['p_binary']:.3f}, ARR = {r['abs_risk_reduction']:+.1%})")
        else:
            print("  VERDICT: Helix curvature suppression NOT FOUND at proteome scale")
            print(f"           Real non-const = {r['real_nonconst_frac']:.1%}, "
                  f"Null non-const = {r['null_nonconst_frac']:.1%}")
            print(f"           Helices are curvature-ENRICHED relative to basin null")
    print("=" * 72)
    print()


# ─── Main ───
def main():
    import argparse
    parser = argparse.ArgumentParser(description="Basin-control null for AlphaFold barcodes")
    parser.add_argument('--reps', type=int, default=5, help='Null reps per segment (default: 5)')
    parser.add_argument('--quick', action='store_true', help='Quick mode: 2 reps, max 200 proteins')
    parser.add_argument('--max-proteins', type=int, default=0, help='Max proteins to analyze (0=all)')
    args = parser.parse_args()
    
    n_reps = 2 if args.quick else args.reps
    max_proteins = 200 if args.quick else args.max_proteins
    
    print()
    print("=" * 72)
    print("  BASIN-CONTROL NULL ANALYSIS — AlphaFold Edition")
    print(f"  Null reps per segment: {n_reps}")
    print("=" * 72)
    print()
    
    # Find cached structures — check both alphafold_cache/ and data/alphafold/
    search_dirs = []
    for d in [ALPHAFOLD_CACHE_DIR, ALPHAFOLD_DIR,
              Path("alphafold_cache"), Path("data/alphafold")]:
        if d.exists() and d not in search_dirs:
            search_dirs.append(d)
    
    if not search_dirs:
        print(f"  [ERROR] No AlphaFold cache found.")
        print(f"  Searched: alphafold_cache/, data/alphafold/")
        sys.exit(1)
    
    structure_files = []
    for d in search_dirs:
        structure_files.extend(d.glob("*.cif"))
        structure_files.extend(d.glob("*.pdb"))
    
    # Deduplicate by filename stem
    seen = set()
    unique_files = []
    for f in sorted(structure_files):
        if f.stem not in seen:
            seen.add(f.stem)
            unique_files.append(f)
    structure_files = unique_files
    
    if max_proteins > 0:
        structure_files = structure_files[:max_proteins]
    
    print(f"  Found {len(structure_files)} cached structures")
    print()
    
    # Phase 1: Load and segment all proteins
    print("  PHASE 1: Loading structures and classifying segments...")
    t0 = time.time()
    
    all_protein_data = []  # list of (phi_v, psi_v, segments)
    n_loaded = 0
    n_failed = 0
    
    for i, fpath in enumerate(structure_files):
        if (i + 1) % 100 == 0:
            print(f"    Loading {i+1}/{len(structure_files)}...")
        
        try:
            phi, psi, res_ids, backbone_coords, has_O = load_structure(fpath)
            segments, ss, phi_v, psi_v = segment_protein(phi, psi, backbone_coords, has_O)
            
            if segments:
                all_protein_data.append((phi_v, psi_v, segments))
                n_loaded += 1
        except Exception as e:
            n_failed += 1
    
    t1 = time.time()
    total_segs = sum(len(segs) for _, _, segs in all_protein_data)
    classifiable_segs = sum(
        1 for _, _, segs in all_protein_data
        for seg in segs if seg['class'] != 'too_short'
    )
    
    print(f"\n  Loaded: {n_loaded} proteins ({n_failed} failed)")
    print(f"  Total segments: {total_segs} ({classifiable_segs} classifiable)")
    print(f"  Time: {t1-t0:.0f}s")
    
    # Quick check: SS distribution of real segments
    ss_counts = {'H': 0, 'E': 0, 'C': 0}
    ss_const = {'H': 0, 'E': 0, 'C': 0}
    for _, _, segs in all_protein_data:
        for seg in segs:
            if seg['class'] != 'too_short' and seg['length'] >= MIN_SEG_LEN:
                ss = seg['ss_type']
                if ss in ss_counts:
                    ss_counts[ss] += 1
                    if seg['macro'] == 'constant_k':
                        ss_const[ss] += 1
    
    print(f"\n  Real segments (before null):")
    for ss in ['H', 'E', 'C']:
        if ss_counts[ss] > 0:
            geo_pct = ss_const[ss] / ss_counts[ss]
            print(f"    {ss}: {ss_counts[ss]} segments, {geo_pct:.1%} constant-κ (geodesic/arc)")
    
    # Phase 2: Build step pools
    print(f"\n  PHASE 2: Building step pools...")
    step_pools, start_pools = build_step_pools(all_protein_data)
    
    for ss in ['H', 'E', 'C']:
        if ss in step_pools:
            print(f"    {ss}: {len(step_pools[ss])} steps, {len(start_pools[ss])} start positions")
    
    # Phase 3: Run basin null
    print(f"\n  PHASE 3: Running basin-control null ({n_reps} reps/segment)...")
    print(f"    Expected null segments: ~{classifiable_segs * n_reps}")
    t2 = time.time()
    
    real_classes, null_classes = run_basin_null(
        all_protein_data, step_pools, start_pools, n_reps=n_reps, seed=SEED
    )
    
    t3 = time.time()
    print(f"\n  Null generation complete: {len(null_classes)} null segments in {t3-t2:.0f}s")
    
    # Phase 4: Statistics
    print(f"\n  PHASE 4: Computing statistics...")
    results = compute_statistics(real_classes, null_classes, all_protein_data)
    
    # Print
    print_results(results)
    
    # Save
    out_path = RESULTS_DIR / "basin_null_results.json"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump({
            'n_proteins': n_loaded,
            'n_reps': n_reps,
            'n_real_segments': len(real_classes),
            'n_null_segments': len(null_classes),
            'results': results,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        }, f, indent=2)
    print(f"  Results saved to: {out_path}")
    print()


if __name__ == "__main__":
    main()

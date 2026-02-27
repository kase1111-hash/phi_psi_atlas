#!/usr/bin/env python3
"""
ALLOSTERIC THRESHOLD SENSITIVITY
=================================

Tests whether the hemoglobin allosteric result is robust across
different significance thresholds (1.0σ to 2.5σ).

If the result holds across thresholds → robust finding.
If it only works at 1.5σ → fragile/anecdotal.

Usage:
    python allosteric_sensitivity.py
"""

import json
import sys
import numpy as np
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

# ═══════════════════════════════════════════════════════════════════════
#  Load the raw Δκ data from allosteric analysis
# ═══════════════════════════════════════════════════════════════════════

# We need to recompute from structures. But since we have the analysis
# script, let's extract curvature profiles directly.

from allosteric_barcodes import (
    extract_curvature_profile, compare_states, DATA_DIR
)

# Allosteric zones (same as permutation test)
ZONES = {
    'O2_binding': {
        'A': list(range(55, 62)),
        'B': list(range(60, 67)),
    },
    'proximal_strain': {
        'A': list(range(84, 91)),
        'B': list(range(70, 77)),
    },
    'alpha1_beta2_interface': {
        'A': list(range(36, 45)) + list(range(130, 142)),
        'B': list(range(99, 112)),
    },
    'DPG_effector': {
        'A': [],
        'B': list(range(79, 86)),
    },
    'Bohr_effect': {
        'A': [],
        'B': list(range(135, 147)),
    },
}


def get_all_zone_residues(chain):
    """Get set of all functional zone residues for a chain."""
    zone_res = set()
    for zone_def in ZONES.values():
        zone_res |= set(zone_def.get(chain, []))
    return zone_res


def analyze_threshold(comp_a, comp_b, n_sigma):
    """Analyze at a given sigma threshold."""
    results = {}
    
    for chain_id, comp in [('A', comp_a), ('B', comp_b)]:
        abs_delta = comp['abs_delta']
        residues = comp['residues']
        
        mean_d = np.mean(abs_delta)
        std_d = np.std(abs_delta)
        threshold = mean_d + n_sigma * std_d
        
        sig_mask = abs_delta > threshold
        sig_res = set(residues[sig_mask])
        n_sig = len(sig_res)
        
        zone_res = get_all_zone_residues(chain_id)
        overlap = sig_res & zone_res
        
        results[chain_id] = {
            'n_sig': n_sig,
            'sig_residues': sig_res,
            'n_overlap': len(overlap),
            'threshold': threshold,
        }
    
    # Combined
    total_sig = results['A']['n_sig'] + results['B']['n_sig']
    total_overlap = results['A']['n_overlap'] + results['B']['n_overlap']
    
    # Zone hits
    sig_a = results['A']['sig_residues']
    sig_b = results['B']['sig_residues']
    
    zone_hits = 0
    for zone_name, zone_def in ZONES.items():
        zone_a = set(zone_def.get('A', []))
        zone_b = set(zone_def.get('B', []))
        if (sig_a & zone_a) or (sig_b & zone_b):
            zone_hits += 1
    
    # Permutation test (quick, 2000 reps)
    rng = np.random.RandomState(42)
    n_perm = 5000
    perm_overlaps = []
    
    res_a_list = list(comp_a['residues'])
    res_b_list = list(comp_b['residues'])
    n_sig_a = results['A']['n_sig']
    n_sig_b = results['B']['n_sig']
    
    zone_a_all = get_all_zone_residues('A')
    zone_b_all = get_all_zone_residues('B')
    
    if n_sig_a > 0 and n_sig_b > 0:
        for _ in range(n_perm):
            pa = set(rng.choice(res_a_list, size=min(n_sig_a, len(res_a_list)), replace=False))
            pb = set(rng.choice(res_b_list, size=min(n_sig_b, len(res_b_list)), replace=False))
            perm_overlaps.append(len(pa & zone_a_all) + len(pb & zone_b_all))
        
        perm_overlaps = np.array(perm_overlaps)
        p_value = np.mean(perm_overlaps >= total_overlap)
        expected = np.mean(perm_overlaps)
    else:
        p_value = 1.0
        expected = 0
    
    return {
        'n_sigma': n_sigma,
        'total_sig': total_sig,
        'total_overlap': total_overlap,
        'zone_hits': zone_hits,
        'p_value': p_value,
        'expected': expected,
        'pct_overlap': total_overlap / total_sig * 100 if total_sig > 0 else 0,
    }


def main():
    print("\n  ALLOSTERIC THRESHOLD SENSITIVITY ANALYSIS")
    print("  " + "═" * 55)
    
    # Load structures and compute profiles
    pdb_t = DATA_DIR / "2dn2.pdb"
    pdb_r = DATA_DIR / "2dn1.pdb"
    
    if not pdb_t.exists() or not pdb_r.exists():
        print("  Missing PDB files. Run allosteric_barcodes.py download first.")
        return
    
    print("  Loading curvature profiles...")
    prof_a_t = extract_curvature_profile(pdb_t, 'A')
    prof_a_r = extract_curvature_profile(pdb_r, 'A')
    prof_b_t = extract_curvature_profile(pdb_t, 'B')
    prof_b_r = extract_curvature_profile(pdb_r, 'B')
    
    comp_a = compare_states(prof_a_t, prof_a_r, smooth_sigma=0)
    comp_b = compare_states(prof_b_t, prof_b_r, smooth_sigma=0)
    
    print("  Running threshold sweep...\n")
    
    thresholds = [0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5]
    
    print(f"  {'σ':>5s} {'n_sig':>6s} {'overlap':>8s} {'%':>6s} {'zones':>6s} "
          f"{'expected':>9s} {'p-value':>9s} {'verdict':>12s}")
    print(f"  {'-'*70}")
    
    all_results = []
    for ns in thresholds:
        r = analyze_threshold(comp_a, comp_b, ns)
        all_results.append(r)
        
        if r['p_value'] < 0.001:
            verdict = '★★★'
        elif r['p_value'] < 0.01:
            verdict = '★★'
        elif r['p_value'] < 0.05:
            verdict = '★'
        else:
            verdict = 'ns'
        
        print(f"  {ns:5.2f} {r['total_sig']:6d} {r['total_overlap']:8d} "
              f"{r['pct_overlap']:5.1f}% {r['zone_hits']:5d}/5 "
              f"{r['expected']:9.1f} {r['p_value']:9.4f} {verdict:>12s}")
    
    print()
    print("  ★ = p < 0.05, ★★ = p < 0.01, ★★★ = p < 0.001, ns = not significant")
    print()
    
    # Summary
    sig_count = sum(1 for r in all_results if r['p_value'] < 0.05)
    print(f"  ROBUSTNESS: {sig_count}/{len(thresholds)} thresholds significant at p < 0.05")
    
    if sig_count >= 5:
        print("  → HIGHLY ROBUST: Result holds across wide threshold range")
    elif sig_count >= 3:
        print("  → MODERATELY ROBUST: Result holds at multiple thresholds")
    elif sig_count >= 1:
        print("  → FRAGILE: Result threshold-dependent")
    else:
        print("  → NOT ROBUST: No threshold yields significance")
    
    # Save
    with open(SCRIPT_DIR / "results" / "allosteric_sensitivity.json", 'w') as f:
        json.dump(all_results, f, indent=2)


if __name__ == "__main__":
    main()

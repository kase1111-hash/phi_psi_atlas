#!/usr/bin/env python3
"""
ALLOSTERIC PERMUTATION TEST
============================

Tests whether the 5/5 zone detection in hemoglobin T→R curvature
differences is statistically significant by shuffling residue labels.

Null hypothesis: significant Δκ residues are randomly distributed
along the chain (no enrichment at allosteric sites).

Test: For 10,000 permutations, randomly place the same NUMBER of
significant residues along the chain and count how often all 5
allosteric zones are hit.

Usage:
    python allosteric_permutation.py
    python allosteric_permutation.py --reps 50000

Author: Kase Knochenhauer / True North Construction LLC
"""

import sys
import numpy as np
from collections import defaultdict

N_REPS = 10000
SEED = 42

# ═══════════════════════════════════════════════════════════════════════
#  OBSERVED DATA from allosteric_barcodes.py
# ═══════════════════════════════════════════════════════════════════════

# Chain lengths (residues with curvature data)
CHAIN_A_LEN = 136
CHAIN_B_LEN = 141

# Significant residues (observed)
SIG_A = [11, 26, 27, 39, 55, 57, 65, 132, 133, 134]
SIG_B = [8, 9, 21, 27, 61, 62, 73, 82, 104, 106, 108, 114, 116, 127, 138]

# Residue ranges for each chain
RES_A = list(range(1, CHAIN_A_LEN + 1))  # 1-136 (approximate)
RES_B = list(range(1, CHAIN_B_LEN + 1))  # 1-141 (approximate)

# ═══════════════════════════════════════════════════════════════════════
#  ALLOSTERIC ZONES (±3 residues around key functional sites)
# ═══════════════════════════════════════════════════════════════════════

# Zone definitions: each zone has residues in α AND/OR β chain
# A hit = at least one significant residue within ±3 of the key site

ZONES = {
    'O2_binding': {
        # Distal histidine: αHis58(E7), βHis63(E7)
        'A': list(range(55, 62)),   # α55-61
        'B': list(range(60, 67)),   # β60-66
    },
    'proximal_strain': {
        # Proximal histidine / F helix start: αHis87(F8), βHis92(F8)
        # β73 is F helix start, mechanically coupled
        'A': list(range(84, 91)),   # α84-90
        'B': list(range(70, 77)),   # β70-76 (F helix start)
    },
    'alpha1_beta2_interface': {
        # Quaternary switch: α38-44, α130-141, β99-108
        'A': list(range(36, 45)) + list(range(130, 142)),
        'B': list(range(99, 112)),
    },
    'DPG_effector': {
        # 2,3-DPG binding: βLys82, βHis2, βHis143
        'A': [],
        'B': list(range(79, 86)),   # β79-85
    },
    'Bohr_effect': {
        # His146 salt bridge, C-terminal region
        'A': [],
        'B': list(range(135, 147)), # β135-146
    },
}


def count_zone_hits(sig_a, sig_b, zones):
    """Count how many zones have at least one significant residue."""
    sig_a_set = set(sig_a)
    sig_b_set = set(sig_b)
    hits = 0
    for zone_name, zone_def in zones.items():
        zone_a = set(zone_def.get('A', []))
        zone_b = set(zone_def.get('B', []))
        if (sig_a_set & zone_a) or (sig_b_set & zone_b):
            hits += 1
    return hits


def count_residue_overlap(sig_a, sig_b, zones):
    """Count total residues overlapping with any zone."""
    sig_a_set = set(sig_a)
    sig_b_set = set(sig_b)
    all_zone_a = set()
    all_zone_b = set()
    for zone_def in zones.values():
        all_zone_a |= set(zone_def.get('A', []))
        all_zone_b |= set(zone_def.get('B', []))
    return len(sig_a_set & all_zone_a) + len(sig_b_set & all_zone_b)


def main():
    global N_REPS
    
    # Parse args
    if '--reps' in sys.argv:
        idx = sys.argv.index('--reps')
        N_REPS = int(sys.argv[idx + 1])
    
    rng = np.random.RandomState(SEED)
    
    # Observed statistics
    n_sig_a = len(SIG_A)
    n_sig_b = len(SIG_B)
    
    observed_zones = count_zone_hits(SIG_A, SIG_B, ZONES)
    observed_residue_overlap = count_residue_overlap(SIG_A, SIG_B, ZONES)
    
    print(f"\n  ALLOSTERIC PERMUTATION TEST")
    print(f"  {'═'*55}")
    print(f"\n  Observed:")
    print(f"    Chain A: {n_sig_a} significant residues out of {CHAIN_A_LEN}")
    print(f"    Chain B: {n_sig_b} significant residues out of {CHAIN_B_LEN}")
    print(f"    Zone hits: {observed_zones}/5")
    print(f"    Residue overlap with zones: {observed_residue_overlap}")
    print(f"\n  Running {N_REPS:,} permutations...\n")
    
    # Permutation test
    zone_hit_counts = np.zeros(6, dtype=int)  # index = number of zones hit (0-5)
    residue_overlap_dist = []
    
    for i in range(N_REPS):
        # Randomly place same number of significant residues
        perm_a = sorted(rng.choice(RES_A, size=n_sig_a, replace=False))
        perm_b = sorted(rng.choice(RES_B, size=n_sig_b, replace=False))
        
        n_zones = count_zone_hits(perm_a, perm_b, ZONES)
        zone_hit_counts[n_zones] += 1
        
        r_overlap = count_residue_overlap(perm_a, perm_b, ZONES)
        residue_overlap_dist.append(r_overlap)
        
        if (i + 1) % 2000 == 0:
            print(f"    {i+1:,}/{N_REPS:,}...")
    
    residue_overlap_dist = np.array(residue_overlap_dist)
    
    # Results
    print(f"\n  {'═'*55}")
    print(f"  RESULTS")
    print(f"  {'═'*55}")
    
    print(f"\n  Zone hit distribution ({N_REPS:,} permutations):")
    for n_zones in range(6):
        count = zone_hit_counts[n_zones]
        pct = count / N_REPS * 100
        bar = '█' * int(pct / 2)
        marker = ' ← OBSERVED' if n_zones == observed_zones else ''
        print(f"    {n_zones}/5 zones: {count:6,} ({pct:5.1f}%) {bar}{marker}")
    
    # P-value for zone test
    p_zones = np.sum(zone_hit_counts[observed_zones:]) / N_REPS
    print(f"\n  P(≥{observed_zones}/5 zones by chance) = {p_zones:.4f}")
    
    # Residue overlap statistics
    mean_overlap = np.mean(residue_overlap_dist)
    std_overlap = np.std(residue_overlap_dist)
    p_residue = np.mean(residue_overlap_dist >= observed_residue_overlap)
    
    print(f"\n  Residue overlap with zones:")
    print(f"    Observed: {observed_residue_overlap}")
    print(f"    Null distribution: {mean_overlap:.1f} ± {std_overlap:.1f}")
    print(f"    P(≥{observed_residue_overlap} overlap) = {p_residue:.4f}")
    
    # Summary
    print(f"\n  {'═'*55}")
    print(f"  SUMMARY")
    print(f"  {'═'*55}")
    print(f"\n  Zone-level test (5/5 zones hit):")
    if p_zones < 0.001:
        print(f"    p = {p_zones:.6f} → HIGHLY SIGNIFICANT (p < 0.001)")
    elif p_zones < 0.01:
        print(f"    p = {p_zones:.4f} → SIGNIFICANT (p < 0.01)")
    elif p_zones < 0.05:
        print(f"    p = {p_zones:.4f} → SIGNIFICANT (p < 0.05)")
    else:
        print(f"    p = {p_zones:.4f} → NOT SIGNIFICANT at p < 0.05")
    
    print(f"\n  Residue-level test ({observed_residue_overlap} overlaps):")
    if p_residue < 0.05:
        print(f"    p = {p_residue:.4f} → SIGNIFICANT")
    else:
        print(f"    p = {p_residue:.4f} → Not significant")
    
    print(f"\n  Interpretation:")
    print(f"    The probability of randomly placed significant")
    print(f"    residues hitting all 5 allosteric zones is {p_zones:.4f}.")
    if p_zones < 0.05:
        print(f"    This confirms that Δκ localizes to functional sites")
        print(f"    beyond chance expectation.")


if __name__ == "__main__":
    main()

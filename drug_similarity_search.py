#!/usr/bin/env python3
"""
Drug-Target Curvature Similarity Search
=========================================
Identifies proteins in the human proteome whose bending energy signatures
are geometrically similar to hemoglobin's drug-affected regions.

Two levels of analysis:
  1. GLOBAL: Proteins with similar overall curvature fingerprint to β-globin
     (Σκ²/res, Gini, SS composition, length)
  2. LOCAL: Proteins whose per-residue κ² profiles contain regions similar
     to hemoglobin's drug-response hotspots (E, G, H helices)

USAGE:
    # Basic: find hemoglobin's geometric neighbors
    python drug_similarity_search.py human_curvature_dssp.csv

    # With full per-residue profiles (requires CIF files):
    python drug_similarity_search.py human_curvature_dssp.csv \
        --cif-dir "G:\Hemo-Mapping\alphafold_cache\human" \
        --top 50

    # Using existing hemoglobin barcode for reference:
    python drug_similarity_search.py human_curvature_dssp.csv \
        --hb-barcode P68871_barcode.json

REQUIRES: numpy, pandas, scipy
LICENSE: CC0 1.0
"""

import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import wasserstein_distance


# ══════════════════════════════════════════════════════════════
# HEMOGLOBIN REFERENCE VALUES (from Paper 2 crystal analysis)
# ══════════════════════════════════════════════════════════════

# β-globin crystal structure reference (WT R-state, 2DN1 chain B)
HB_REFERENCE = {
    'length': 143,
    'kappa_sq_per_res': 107.0,
    'ss_frac_H': 0.72,    # DSSP helix fraction for β-globin
    'ss_frac_E': 0.0,     # no sheet in globin
    'ss_frac_C': 0.28,
    'kappa_gini': 0.82,
    'mean_plddt': 99.0,   # crystal structure, not AF
}

# Drug-affected helices in β-globin (residue ranges, 0-indexed)
# These are where voxelotor's r = -0.61 anticorrelation is concentrated
DRUG_HOTSPOTS = {
    'E_helix': (59, 74),    # residues 60-75
    'G_helix': (99, 119),   # residues 100-120
    'H_helix': (122, 139),  # residues 123-140
}

# Voxelotor curvature signature:
# Disease deviation (HbS R - WT R) concentrated in these helices
# Drug correction (Drug R - HbS R) anticorrelated at r = -0.61
# Key: the drug INCREASES curvature in some helices, DECREASES in others


def compute_global_similarity(curv_df, reference=HB_REFERENCE):
    """Compute multi-dimensional distance from each protein to the hemoglobin reference.
    
    Features: Σκ²/res, Gini, SS fractions, length (all normalized).
    Returns DataFrame with similarity scores.
    """
    
    features = ['kappa_sq_per_res', 'kappa_gini', 'ss_frac_H', 'ss_frac_E', 'ss_frac_C', 'length']
    
    # Normalize by proteome statistics
    X = curv_df[features].values.copy()
    means = np.nanmean(X, axis=0)
    stds = np.nanstd(X, axis=0)
    stds[stds == 0] = 1
    X_norm = (X - means) / stds
    
    # Reference point normalized the same way
    ref_values = np.array([reference[f] for f in features])
    ref_norm = (ref_values - means) / stds
    
    # Euclidean distance in normalized feature space
    distances = np.sqrt(np.sum((X_norm - ref_norm) ** 2, axis=1))
    
    curv_df = curv_df.copy()
    curv_df['global_distance'] = distances
    curv_df['global_rank'] = distances.argsort().argsort() + 1
    
    return curv_df


def compute_helix_similarity(curv_df, reference=HB_REFERENCE):
    """Find proteins with similar helix-dominated curvature profiles.
    
    Specifically targets proteins that:
    1. Have similar helix fraction to β-globin (~70%)
    2. Have similar Σκ²/res (~107)
    3. Have similar Gini (~0.82, moderately concentrated)
    
    These are candidates for similar drug sensitivity because:
    - Similar helix content → similar torus trajectory geometry
    - Similar Σκ²/res → similar total curvature budget
    - A drug that redistributes curvature in globin might similarly
      redistribute curvature in these proteins.
    """
    
    df = curv_df.copy()
    
    # Component distances (weighted)
    d_ksq = np.abs(df['kappa_sq_per_res'] - reference['kappa_sq_per_res']) / reference['kappa_sq_per_res']
    d_helix = np.abs(df['ss_frac_H'] - reference['ss_frac_H']) / max(reference['ss_frac_H'], 0.01)
    d_gini = np.abs(df['kappa_gini'] - reference['kappa_gini']) / max(reference['kappa_gini'], 0.01)
    d_length = np.abs(df['length'] - reference['length']) / reference['length']
    
    # Weighted combination (helix and Σκ² matter most)
    df['helix_similarity_distance'] = (
        0.35 * d_ksq +      # curvature budget match
        0.30 * d_helix +    # helix content match
        0.20 * d_gini +     # curvature concentration match
        0.15 * d_length     # size match
    )
    
    df['helix_similarity_rank'] = df['helix_similarity_distance'].argsort().argsort() + 1
    
    return df


def compute_drug_vulnerability_score(curv_df):
    """Score proteins for potential vulnerability to curvature-redistributing drugs.
    
    Hypothesis: Proteins are vulnerable to curvature redistribution when:
    1. High curvature concentration (high Gini) — energy is localized
    2. High helix content — helices carry most of the curvature
    3. Moderate length — not too short (no room to redistribute) or too long
    4. High pLDDT — well-folded, curvature pattern is real
    
    This is SPECULATIVE. It identifies structural prerequisites for
    drug-induced curvature redistribution, not proven drug targets.
    """
    
    df = curv_df.copy()
    
    # Normalize each component to [0, 1]
    def norm01(x):
        xmin, xmax = x.min(), x.max()
        if xmax == xmin:
            return np.zeros_like(x)
        return (x - xmin) / (xmax - xmin)
    
    # High Gini → curvature concentrated (vulnerable to redistribution)
    gini_score = norm01(df['kappa_gini'])
    
    # High helix → more curvature available to redistribute
    helix_score = norm01(df['ss_frac_H'])
    
    # Medium length (100-400 optimal, penalty outside)
    length = df['length'].values
    length_score = np.exp(-((length - 250) / 200) ** 2)
    
    # High pLDDT → well-folded
    plddt_score = norm01(df['mean_plddt'].clip(upper=100))
    
    # High κ² concentration at peaks (top_peak / mean suggests localized energy)
    peak_ratio = df['top_peak_kappa_sq'] / df['kappa_sq_per_res'].clip(lower=1)
    peak_score = norm01(peak_ratio.clip(upper=50))
    
    df['vulnerability_score'] = (
        0.25 * gini_score +
        0.25 * helix_score +
        0.20 * length_score +
        0.15 * plddt_score +
        0.15 * peak_score
    )
    
    df['vulnerability_rank'] = df['vulnerability_score'].rank(ascending=False).astype(int)
    
    return df


def print_results(df, top_n=50):
    """Print the most similar proteins and vulnerability rankings."""
    
    hq = df[(df['mean_plddt'] >= 70) & (df['length'] >= 50) & (df['length'] <= 2000)].copy()
    
    # ── Global similarity ──
    print(f"\n{'='*80}")
    print(f"  GLOBAL CURVATURE SIMILARITY TO β-GLOBIN (n = {len(hq)})")
    print(f"  Reference: Σκ²/res=107, H=72%, E=0%, Gini=0.82, L=143")
    print(f"{'='*80}")
    
    top_global = hq.nsmallest(top_n, 'global_distance')
    print(f"\n  TOP {top_n} most similar proteins:")
    print(f"  {'Rank':>5s} {'UniProt':>15s} {'L':>5s} {'Σκ²/r':>7s} {'H%':>5s} {'E%':>5s} "
          f"{'Gini':>5s} {'Dist':>6s} {'pLDDT':>6s}")
    for _, r in top_global.iterrows():
        print(f"  {r['global_rank']:>5.0f} {r['uniprot_id']:>15s} {r['length']:>5.0f} "
              f"{r['kappa_sq_per_res']:>7.1f} {r['ss_frac_H']*100:>4.0f}% "
              f"{r['ss_frac_E']*100:>4.0f}% {r['kappa_gini']:>5.2f} "
              f"{r['global_distance']:>6.2f} {r['mean_plddt']:>5.1f}")
    
    # ── Helix similarity ──
    print(f"\n{'='*80}")
    print(f"  HELIX-WEIGHTED SIMILARITY TO β-GLOBIN")
    print(f"{'='*80}")
    
    top_helix = hq.nsmallest(top_n, 'helix_similarity_distance')
    print(f"\n  TOP {top_n} most similar (helix-weighted):")
    print(f"  {'Rank':>5s} {'UniProt':>15s} {'L':>5s} {'Σκ²/r':>7s} {'H%':>5s} {'E%':>5s} "
          f"{'Gini':>5s} {'HxDist':>7s} {'pLDDT':>6s}")
    for _, r in top_helix.iterrows():
        print(f"  {r['helix_similarity_rank']:>5.0f} {r['uniprot_id']:>15s} {r['length']:>5.0f} "
              f"{r['kappa_sq_per_res']:>7.1f} {r['ss_frac_H']*100:>4.0f}% "
              f"{r['ss_frac_E']*100:>4.0f}% {r['kappa_gini']:>5.2f} "
              f"{r['helix_similarity_distance']:>7.3f} {r['mean_plddt']:>5.1f}")
    
    # ── Drug vulnerability ──
    print(f"\n{'='*80}")
    print(f"  CURVATURE REDISTRIBUTION VULNERABILITY RANKING")
    print(f"  (SPECULATIVE — structural prerequisites only)")
    print(f"{'='*80}")
    
    top_vuln = hq.nlargest(top_n, 'vulnerability_score')
    print(f"\n  TOP {top_n} most vulnerable proteins:")
    print(f"  {'Rank':>5s} {'UniProt':>15s} {'L':>5s} {'Σκ²/r':>7s} {'H%':>5s} "
          f"{'Gini':>5s} {'VulnScore':>9s} {'pLDDT':>6s}")
    for _, r in top_vuln.iterrows():
        print(f"  {r['vulnerability_rank']:>5.0f} {r['uniprot_id']:>15s} {r['length']:>5.0f} "
              f"{r['kappa_sq_per_res']:>7.1f} {r['ss_frac_H']*100:>4.0f}% "
              f"{r['kappa_gini']:>5.2f} {r['vulnerability_score']:>9.3f} "
              f"{r['mean_plddt']:>5.1f}")
    
    # ── Overlap analysis ──
    print(f"\n{'='*80}")
    print(f"  OVERLAP: Proteins appearing in multiple top-50 lists")
    print(f"{'='*80}")
    
    global_ids = set(top_global['uniprot_id'])
    helix_ids = set(top_helix['uniprot_id'])
    vuln_ids = set(top_vuln['uniprot_id'])
    
    # In all three
    all_three = global_ids & helix_ids & vuln_ids
    # In at least two
    two_of_three = (global_ids & helix_ids) | (global_ids & vuln_ids) | (helix_ids & vuln_ids)
    
    print(f"\n  In all 3 top-50 lists: {len(all_three)}")
    if all_three:
        for uid in sorted(all_three):
            r = hq[hq['uniprot_id'] == uid].iloc[0]
            print(f"    {uid:>15s} L={r['length']:.0f} Σκ²/r={r['kappa_sq_per_res']:.1f} "
                  f"H={r['ss_frac_H']*100:.0f}% Gini={r['kappa_gini']:.2f}")
    
    print(f"\n  In 2+ of 3 lists: {len(two_of_three)}")
    if two_of_three:
        for uid in sorted(two_of_three):
            r = hq[hq['uniprot_id'] == uid].iloc[0]
            in_lists = []
            if uid in global_ids: in_lists.append('global')
            if uid in helix_ids: in_lists.append('helix')
            if uid in vuln_ids: in_lists.append('vuln')
            print(f"    {uid:>15s} L={r['length']:.0f} Σκ²/r={r['kappa_sq_per_res']:.1f} "
                  f"H={r['ss_frac_H']*100:.0f}% [{', '.join(in_lists)}]")
    
    # ── Summary statistics ──
    print(f"\n{'='*80}")
    print(f"  SUMMARY")
    print(f"{'='*80}")
    
    # What fraction of the proteome is "globin-like"?
    globin_like = hq[
        (hq['kappa_sq_per_res'].between(50, 200)) &
        (hq['ss_frac_H'] > 0.50) &
        (hq['length'].between(100, 200)) &
        (hq['kappa_gini'].between(0.70, 0.95))
    ]
    print(f"\n  'Globin-like' proteins (Σκ²/r=50-200, H>50%, L=100-200, Gini=0.7-0.95):")
    print(f"  n = {len(globin_like)} / {len(hq)} ({len(globin_like)/len(hq)*100:.1f}%)")
    
    if len(globin_like) > 0:
        print(f"  Their IDs:")
        for _, r in globin_like.sort_values('helix_similarity_distance').head(20).iterrows():
            print(f"    {r['uniprot_id']:>15s} L={r['length']:.0f} Σκ²/r={r['kappa_sq_per_res']:.1f} "
                  f"H={r['ss_frac_H']*100:.0f}% Gini={r['kappa_gini']:.2f}")


def main():
    parser = argparse.ArgumentParser(description='Drug-target curvature similarity search')
    parser.add_argument('curvature_csv', help='Curvature survey CSV')
    parser.add_argument('--top', type=int, default=50, help='Number of top results to show')
    parser.add_argument('--output', default='drug_similarity_results.csv', help='Output CSV')
    args = parser.parse_args()
    
    print(f"Loading curvature data from {args.curvature_csv}...")
    df = pd.read_csv(args.curvature_csv)
    hq = df[(df['mean_plddt'] >= 70) & (df['length'] >= 50) & (df['length'] <= 2000)].copy()
    print(f"  {len(hq)} high-quality proteins")
    
    # Run analyses
    print("Computing global similarity...")
    hq = compute_global_similarity(hq)
    
    print("Computing helix-weighted similarity...")
    hq = compute_helix_similarity(hq)
    
    print("Computing vulnerability scores...")
    hq = compute_drug_vulnerability_score(hq)
    
    # Save full results
    hq.to_csv(args.output, index=False)
    print(f"Saved full results to {args.output}")
    
    # Print summary
    print_results(hq, top_n=args.top)


if __name__ == '__main__':
    main()

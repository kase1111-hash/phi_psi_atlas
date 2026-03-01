#!/usr/bin/env python3
"""
Curvature Profile Comparison
==============================
Computes per-residue κ² profiles for hemoglobin's top geometric neighbors
and measures pairwise Wasserstein distances.

Tests: Do curvature-similar proteins (by global metrics) also share
similar SPATIAL curvature patterns?

USAGE:
    python profile_comparison.py human_curvature_dssp.csv \
        --cif-dir "G:\Hemo-Mapping\alphafold_cache\human" \
        --top 30

REQUIRES: numpy, pandas, scipy. Uses proteome_curvature_survey.py functions.
LICENSE: CC0 1.0
"""

import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import wasserstein_distance, pearsonr
from scipy.interpolate import interp1d

try:
    from proteome_curvature_survey import (
        parse_mmcif_backbone, compute_dihedrals, torus_curvature_standalone
    )
except ImportError:
    print("ERROR: proteome_curvature_survey.py must be in same directory or on PYTHONPATH")
    sys.exit(1)


def compute_kappa_sq_profile(cif_path, smooth_window=5):
    """Compute smoothed κ² profile for a single protein.
    Returns (kappa_sq_smoothed, length) or (None, 0) on failure.
    """
    try:
        residues = parse_mmcif_backbone(str(cif_path))
        if len(residues) < 10:
            return None, 0
        
        phi, psi, plddt, n = compute_dihedrals(residues)
        if phi is None or n < 10:
            return None, 0
        
        kappa, _ = torus_curvature_standalone(phi, psi)
        if len(kappa) < 5:
            return None, 0
        
        kappa_sq = kappa ** 2
        
        # Smooth with uniform kernel
        if smooth_window > 1 and len(kappa_sq) > smooth_window:
            kernel = np.ones(smooth_window) / smooth_window
            kappa_sq_smooth = np.convolve(kappa_sq, kernel, mode='same')
        else:
            kappa_sq_smooth = kappa_sq
        
        return kappa_sq_smooth, len(kappa_sq_smooth)
    
    except Exception as e:
        return None, 0


def normalize_profile(profile, n_points=100):
    """Resample profile to fixed length for comparison.
    Returns normalized profile of length n_points.
    """
    if len(profile) < 3:
        return np.zeros(n_points)
    
    x = np.linspace(0, 1, len(profile))
    x_new = np.linspace(0, 1, n_points)
    f = interp1d(x, profile, kind='linear', fill_value='extrapolate')
    resampled = f(x_new)
    
    # Normalize to unit area (probability distribution)
    total = resampled.sum()
    if total > 0:
        resampled = resampled / total
    
    return resampled


def find_cif_for_uniprot(uniprot_id, cif_dir):
    """Find CIF file for a UniProt ID in the directory."""
    cif_dir = Path(cif_dir)
    
    # Try common AlphaFold naming patterns
    patterns = [
        f"AF-{uniprot_id}-F1-model_v*.cif",
        f"*{uniprot_id}*.cif",
    ]
    
    for pattern in patterns:
        matches = list(cif_dir.glob(pattern))
        if matches:
            return str(matches[0])
    
    return None


def main():
    parser = argparse.ArgumentParser(description='Curvature profile comparison')
    parser.add_argument('curvature_csv', help='Curvature survey CSV with similarity rankings')
    parser.add_argument('--cif-dir', required=True, help='Directory containing CIF files')
    parser.add_argument('--top', type=int, default=30, help='Number of top similar proteins to compare')
    parser.add_argument('--output', default='profile_comparison_results.csv', help='Output CSV')
    args = parser.parse_args()
    
    # Load ranked data
    df = pd.read_csv(args.curvature_csv)
    
    # Check if similarity columns exist; if not, compute them
    if 'helix_similarity_rank' not in df.columns:
        print("Similarity rankings not found. Run drug_similarity_search.py first,")
        print("or pass the drug_similarity_results.csv file.")
        sys.exit(1)
    
    # Get top N by helix similarity
    top_ids = df.nsmallest(args.top, 'helix_similarity_distance')['uniprot_id'].tolist()
    
    print(f"Computing κ² profiles for top {len(top_ids)} hemoglobin-similar proteins...")
    
    # Compute profiles
    profiles = {}
    raw_profiles = {}
    
    for i, uid in enumerate(top_ids):
        cif_path = find_cif_for_uniprot(uid, args.cif_dir)
        if cif_path is None:
            print(f"  [{i+1}/{len(top_ids)}] {uid}: CIF not found")
            continue
        
        profile, length = compute_kappa_sq_profile(cif_path)
        if profile is not None:
            raw_profiles[uid] = profile
            profiles[uid] = normalize_profile(profile)
            print(f"  [{i+1}/{len(top_ids)}] {uid}: L={length}, profile computed")
        else:
            print(f"  [{i+1}/{len(top_ids)}] {uid}: failed to compute")
    
    if len(profiles) < 3:
        print("Too few profiles computed. Check CIF directory path.")
        sys.exit(1)
    
    print(f"\nComputed {len(profiles)} profiles. Computing pairwise distances...")
    
    # Pairwise Wasserstein distances between normalized profiles
    uids = list(profiles.keys())
    n = len(uids)
    W_matrix = np.zeros((n, n))
    R_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            w = wasserstein_distance(profiles[uids[i]], profiles[uids[j]])
            W_matrix[i, j] = w
            W_matrix[j, i] = w
            
            r, _ = pearsonr(profiles[uids[i]], profiles[uids[j]])
            R_matrix[i, j] = r
            R_matrix[j, i] = r
    
    # Results
    print(f"\n{'='*70}")
    print(f"  PAIRWISE PROFILE COMPARISON ({n} proteins)")
    print(f"{'='*70}")
    
    # Mean pairwise Wasserstein
    W_upper = W_matrix[np.triu_indices(n, k=1)]
    R_upper = R_matrix[np.triu_indices(n, k=1)]
    
    print(f"\n  Pairwise Wasserstein distance (normalized profiles):")
    print(f"    Mean: {W_upper.mean():.4f}")
    print(f"    SD:   {W_upper.std():.4f}")
    print(f"    Min:  {W_upper.min():.4f}")
    print(f"    Max:  {W_upper.max():.4f}")
    
    print(f"\n  Pairwise Pearson correlation:")
    print(f"    Mean: {R_upper.mean():.3f}")
    print(f"    SD:   {R_upper.std():.3f}")
    print(f"    Min:  {R_upper.min():.3f}")
    print(f"    Max:  {R_upper.max():.3f}")
    
    # Find most similar pairs
    print(f"\n  TOP 10 MOST SIMILAR PROFILE PAIRS:")
    pairs = []
    for i in range(n):
        for j in range(i+1, n):
            pairs.append((uids[i], uids[j], W_matrix[i,j], R_matrix[i,j]))
    
    pairs.sort(key=lambda x: x[2])
    for uid1, uid2, w, r in pairs[:10]:
        l1 = len(raw_profiles.get(uid1, []))
        l2 = len(raw_profiles.get(uid2, []))
        print(f"    {uid1:>12s} - {uid2:<12s}: W={w:.4f}, r={r:+.3f} (L={l1},{l2})")
    
    # ── Curvature peak position analysis ──
    print(f"\n{'='*70}")
    print(f"  CURVATURE PEAK POSITIONS (fraction of chain length)")
    print(f"{'='*70}")
    
    for uid in uids[:15]:
        raw = raw_profiles[uid]
        # Find top 3 peaks
        if len(raw) < 10:
            continue
        peak_positions = []
        sorted_idx = np.argsort(raw)[::-1]
        used = set()
        for idx in sorted_idx:
            # Skip if too close to already-found peak
            if any(abs(idx - u) < 5 for u in used):
                continue
            peak_positions.append(idx / len(raw))
            used.add(idx)
            if len(peak_positions) >= 3:
                break
        
        row = df[df['uniprot_id'] == uid]
        L = row['length'].values[0] if len(row) > 0 else len(raw)
        H = row['ss_frac_H'].values[0] if len(row) > 0 else 0
        
        peaks_str = ', '.join(f'{p:.2f}' for p in sorted(peak_positions))
        print(f"  {uid:>12s} (L={L:.0f}, H={H*100:.0f}%): peaks at [{peaks_str}]")
    
    # ── Save results ──
    results = []
    for uid in uids:
        row = df[df['uniprot_id'] == uid].iloc[0]
        mean_W = W_matrix[uids.index(uid)].mean() if n > 1 else 0
        results.append({
            'uniprot_id': uid,
            'length': row['length'],
            'kappa_sq_per_res': row['kappa_sq_per_res'],
            'ss_frac_H': row['ss_frac_H'],
            'kappa_gini': row['kappa_gini'],
            'helix_similarity_rank': row['helix_similarity_rank'],
            'mean_W_to_group': round(mean_W, 4),
            'profile_length': len(raw_profiles[uid]),
        })
    
    pd.DataFrame(results).to_csv(args.output, index=False)
    print(f"\nSaved to {args.output}")


if __name__ == '__main__':
    main()

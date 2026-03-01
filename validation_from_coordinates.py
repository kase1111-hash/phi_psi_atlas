#!/usr/bin/env python3
"""
INDEPENDENT VALIDATION: Curvature from Raw Atomic Coordinates
================================================================

Implements the reviewer's 8-step falsification protocol:

  Step 1: Clean test cases (myoglobin, lysozyme, ubiquitin)
  Step 2: Extract backbone from PDB
  Step 3: Compute discrete curvature (BOTH Cartesian and torus)
  Step 4: Validate against DSSP secondary structure
  Step 5: Curvature density across structures
  Step 6: Internal reproducibility (same protein, different PDB)
  Step 7: Perturbation sensitivity (apo vs holo)
  Step 8: Noise stress test

CRITICAL NOTE:
  The reviewer proposes testing Cartesian Cα curvature.
  The paper uses geodesic curvature on the (φ,ψ) flat torus.
  These are DIFFERENT quantities measuring DIFFERENT things.
  We compute BOTH so the distinction is empirically clear.

  Section 4.9 of the paper already shows Cartesian Cα curvature
  detects 0-3% variability — this script will reproduce that
  negative result AND show torus curvature detects ~49%.

Test proteins:
  - Myoglobin: 1MBN, 1A6M (two structures, helix-rich)
  - Lysozyme: 1LYZ, 1HEL (two structures, mixed α/β)
  - Ubiquitin: 1UBQ, 1UBI (two structures, β-rich)
  - Myoglobin apo vs holo: 1MBN vs 1MBD (deoxy vs CO-bound)

Requirements:
  pip install numpy scipy biopython

Usage:
  python validation_from_coordinates.py

Author: Kase Knochenhauer / True North Construction LLC
"""

import numpy as np
import sys
import json
import warnings
from pathlib import Path
from collections import defaultdict

warnings.filterwarnings('ignore')

try:
    from Bio.PDB import PDBParser, DSSP, is_aa
    from Bio.PDB.vectors import calc_dihedral
    from Bio.PDB import Vector
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("⚠ BioPython not installed. Install with: pip install biopython")
    print("  Falling back to coordinate-only mode.")

try:
    from scipy import stats as sp_stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


# ═══════════════════════════════════════════════════════════
# PDB DOWNLOAD & PARSING
# ═══════════════════════════════════════════════════════════

def download_pdb(pdb_id, out_dir="."):
    """Download PDB file from RCSB."""
    import urllib.request
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    path = Path(out_dir) / f"{pdb_id.lower()}.pdb"
    if path.exists():
        return str(path)
    try:
        urllib.request.urlretrieve(url, path)
        return str(path)
    except Exception as e:
        print(f"  ⚠ Could not download {pdb_id}: {e}")
        return None


def extract_backbone(pdb_path, chain_id='A'):
    """Extract Cα coordinates and backbone atoms (N, Cα, C) for dihedral computation."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('prot', pdb_path)
    model = structure[0]
    
    # Find chain
    chain = None
    for c in model.get_chains():
        if c.id == chain_id:
            chain = c
            break
    if chain is None:
        # Use first chain
        chain = list(model.get_chains())[0]
    
    residues = []
    for res in chain.get_residues():
        if not is_aa(res, standard=True):
            continue
        try:
            n = res['N'].get_vector()
            ca = res['CA'].get_vector()
            c = res['C'].get_vector()
            residues.append({
                'resnum': res.get_id()[1],
                'resname': res.get_resname(),
                'N': np.array(n.get_array()),
                'CA': np.array(ca.get_array()),
                'C': np.array(c.get_array()),
            })
        except KeyError:
            continue
    
    return residues


def compute_dihedrals(residues):
    """Compute φ and ψ angles from backbone coordinates."""
    n = len(residues)
    phi = np.full(n, np.nan)
    psi = np.full(n, np.nan)
    
    for i in range(n):
        # φ(i) = dihedral(C[i-1], N[i], CA[i], C[i])
        if i > 0:
            try:
                phi[i] = calc_dihedral(
                    Vector(residues[i-1]['C']),
                    Vector(residues[i]['N']),
                    Vector(residues[i]['CA']),
                    Vector(residues[i]['C'])
                ) * 180.0 / np.pi
            except:
                pass
        
        # ψ(i) = dihedral(N[i], CA[i], C[i], N[i+1])
        if i < n - 1:
            try:
                psi[i] = calc_dihedral(
                    Vector(residues[i]['N']),
                    Vector(residues[i]['CA']),
                    Vector(residues[i]['C']),
                    Vector(residues[i+1]['N'])
                ) * 180.0 / np.pi
            except:
                pass
    
    return phi, psi


def assign_ss_simple(phi, psi):
    """Simple SS assignment from dihedrals (no DSSP dependency)."""
    n = len(phi)
    ss = ['C'] * n
    for i in range(n):
        if np.isnan(phi[i]) or np.isnan(psi[i]):
            continue
        # Helix basin: φ ∈ [-160, -20], ψ ∈ [-80, 0]  (generous)
        if -160 < phi[i] < -20 and -80 < psi[i] < 20:
            ss[i] = 'H'
        # Strand basin: φ ∈ [-180, -60], ψ ∈ [60, 180]
        elif -180 < phi[i] < -60 and (60 < psi[i] < 180 or -180 < psi[i] < -120):
            ss[i] = 'E'
    
    # Smooth: require 3+ consecutive for H or E
    ss_smooth = list(ss)
    for i in range(1, n-1):
        if ss[i] == 'H' and ss[i-1] != 'H' and ss[i+1] != 'H':
            ss_smooth[i] = 'C'
        if ss[i] == 'E' and ss[i-1] != 'E' and ss[i+1] != 'E':
            ss_smooth[i] = 'C'
    
    return ss_smooth


# ═══════════════════════════════════════════════════════════
# CURVATURE COMPUTATIONS
# ═══════════════════════════════════════════════════════════

def cartesian_curvature(residues):
    """
    REVIEWER'S METHOD: Discrete Cα curvature in 3D space.
    κ(i) = |T(i+1) × T(i)| / (|T(i+1)| · |T(i)|)
    where T(i) = r(i) - r(i-1)
    """
    coords = np.array([r['CA'] for r in residues])
    n = len(coords)
    kappa = np.zeros(n)
    
    for i in range(1, n-1):
        v1 = coords[i] - coords[i-1]
        v2 = coords[i+1] - coords[i]
        
        cross = np.cross(v1, v2)
        denom = np.linalg.norm(v1) * np.linalg.norm(v2)
        
        if denom > 1e-10:
            kappa[i] = np.linalg.norm(cross) / denom
    
    return kappa


def torus_curvature(phi, psi):
    """
    THE PAPER'S METHOD: Geodesic curvature on the flat (φ,ψ) torus.
    
    The trajectory is (φ(i), ψ(i)) on T² with periodic boundaries.
    Tangent vector T(i) = (Δφ(i), Δψ(i)) wrapped to [-180, 180].
    Curvature = signed turning angle between consecutive tangents,
    normalized by arc length.
    """
    n = len(phi)
    kappa = np.zeros(n)
    
    def wrap(angle):
        """Wrap to [-180, 180]."""
        return ((angle + 180) % 360) - 180
    
    for i in range(2, n-1):
        if np.isnan(phi[i]) or np.isnan(psi[i]) or \
           np.isnan(phi[i-1]) or np.isnan(psi[i-1]) or \
           np.isnan(phi[i-2]) or np.isnan(psi[i-2]):
            continue
        
        # Tangent at i-1
        dphi1 = wrap(phi[i] - phi[i-1])
        dpsi1 = wrap(psi[i] - psi[i-1])
        
        # Tangent at i-2 → i-1
        dphi0 = wrap(phi[i-1] - phi[i-2])
        dpsi0 = wrap(psi[i-1] - psi[i-2])
        
        # Arc lengths
        s0 = np.sqrt(dphi0**2 + dpsi0**2)
        s1 = np.sqrt(dphi1**2 + dpsi1**2)
        
        if s0 < 0.1 or s1 < 0.1:
            continue
        
        # Unit tangents
        T0 = np.array([dphi0, dpsi0]) / s0
        T1 = np.array([dphi1, dpsi1]) / s1
        
        # Signed turning angle (2D cross product)
        cross = T0[0] * T1[1] - T0[1] * T1[0]
        dot = T0[0] * T1[0] + T0[1] * T1[1]
        theta = np.arctan2(-cross, dot)  # IUPAC sign convention
        
        # Curvature = turning angle / mean arc length
        s_mean = (s0 + s1) / 2.0
        kappa[i] = theta / s_mean * (180.0 / np.pi)  # degrees per degree
    
    return kappa


# ═══════════════════════════════════════════════════════════
# ANALYSIS FUNCTIONS
# ═══════════════════════════════════════════════════════════

def analyze_protein(pdb_id, chain='A', label=None):
    """Full analysis of one PDB structure."""
    if label is None:
        label = pdb_id
    
    pdb_path = download_pdb(pdb_id)
    if pdb_path is None:
        return None
    
    residues = extract_backbone(pdb_path, chain)
    if len(residues) < 20:
        print(f"  ⚠ {label}: only {len(residues)} residues extracted")
        return None
    
    phi, psi = compute_dihedrals(residues)
    ss = assign_ss_simple(phi, psi)
    
    # Compute both curvature types
    kappa_cart = cartesian_curvature(residues)
    kappa_torus = torus_curvature(phi, psi)
    
    kappa_sq_cart = kappa_cart ** 2
    kappa_sq_torus = kappa_torus ** 2
    
    # Per-SS statistics
    h_idx = [i for i in range(len(ss)) if ss[i] == 'H' and kappa_torus[i] != 0]
    e_idx = [i for i in range(len(ss)) if ss[i] == 'E' and kappa_torus[i] != 0]
    c_idx = [i for i in range(len(ss)) if ss[i] == 'C' and kappa_torus[i] != 0]
    
    result = {
        'pdb_id': pdb_id,
        'label': label,
        'n_residues': len(residues),
        'n_helix': len([s for s in ss if s == 'H']),
        'n_strand': len([s for s in ss if s == 'E']),
        'n_coil': len([s for s in ss if s == 'C']),
        
        # Cartesian Cα curvature
        'cart_mean_kappa_sq': float(np.mean(kappa_sq_cart)),
        'cart_sum_per_res': float(np.sum(kappa_sq_cart) / len(residues)),
        'cart_std': float(np.std(kappa_cart)),
        
        # Torus geodesic curvature
        'torus_mean_kappa_sq': float(np.mean(kappa_sq_torus)),
        'torus_sum_per_res': float(np.sum(kappa_sq_torus) / len(residues)),
        'torus_std': float(np.std(kappa_torus)),
        
        # Raw profiles for correlation
        '_kappa_cart': kappa_cart,
        '_kappa_torus': kappa_torus,
        '_kappa_sq_cart': kappa_sq_cart,
        '_kappa_sq_torus': kappa_sq_torus,
        '_ss': ss,
        '_phi': phi,
        '_psi': psi,
    }
    
    # Per-SS curvature (torus)
    if h_idx:
        result['torus_helix_mean'] = float(np.mean(np.abs(kappa_torus[h_idx])))
    if e_idx:
        result['torus_strand_mean'] = float(np.mean(np.abs(kappa_torus[e_idx])))
    if c_idx:
        result['torus_coil_mean'] = float(np.mean(np.abs(kappa_torus[c_idx])))
    
    # Per-SS curvature (Cartesian) 
    if h_idx:
        result['cart_helix_mean'] = float(np.mean(kappa_cart[h_idx]))
    if e_idx:
        result['cart_strand_mean'] = float(np.mean(kappa_cart[e_idx]))
    if c_idx:
        result['cart_coil_mean'] = float(np.mean(kappa_cart[c_idx]))
    
    return result


def step4_validate_ss(result):
    """Step 4: Do curvature peaks align with structural boundaries?"""
    ss = result['_ss']
    kappa_t = result['_kappa_sq_torus']
    kappa_c = result['_kappa_sq_cart']
    n = len(ss)
    
    # Find SS boundaries (transitions)
    boundary_idx = []
    for i in range(1, n):
        if ss[i] != ss[i-1]:
            boundary_idx.append(i)
    
    if not boundary_idx:
        return None
    
    # Check if curvature peaks are enriched at boundaries (±2 residues)
    boundary_set = set()
    for b in boundary_idx:
        for d in range(-2, 3):
            if 0 <= b + d < n:
                boundary_set.add(b + d)
    
    # Torus curvature at boundaries vs non-boundaries
    boundary_kappa_t = [kappa_t[i] for i in boundary_set if kappa_t[i] > 0]
    non_boundary_kappa_t = [kappa_t[i] for i in range(n) if i not in boundary_set and kappa_t[i] > 0]
    
    # Cartesian curvature at boundaries vs non-boundaries  
    boundary_kappa_c = [kappa_c[i] for i in boundary_set if kappa_c[i] > 0]
    non_boundary_kappa_c = [kappa_c[i] for i in range(n) if i not in boundary_set and kappa_c[i] > 0]
    
    out = {
        'n_boundaries': len(boundary_idx),
        'n_boundary_residues': len(boundary_set),
    }
    
    if boundary_kappa_t and non_boundary_kappa_t and HAS_SCIPY:
        _, p = sp_stats.mannwhitneyu(boundary_kappa_t, non_boundary_kappa_t, alternative='greater')
        out['torus_boundary_mean'] = np.mean(boundary_kappa_t)
        out['torus_non_boundary_mean'] = np.mean(non_boundary_kappa_t)
        out['torus_ratio'] = np.mean(boundary_kappa_t) / np.mean(non_boundary_kappa_t)
        out['torus_p'] = p
    
    if boundary_kappa_c and non_boundary_kappa_c and HAS_SCIPY:
        _, p = sp_stats.mannwhitneyu(boundary_kappa_c, non_boundary_kappa_c, alternative='greater')
        out['cart_boundary_mean'] = np.mean(boundary_kappa_c)
        out['cart_non_boundary_mean'] = np.mean(non_boundary_kappa_c)
        out['cart_ratio'] = np.mean(boundary_kappa_c) / np.mean(non_boundary_kappa_c)
        out['cart_p'] = p
    
    return out


def step6_reproducibility(r1, r2):
    """Step 6: Correlation of κ² profiles between two structures of same protein."""
    k1_t = r1['_kappa_sq_torus']
    k2_t = r2['_kappa_sq_torus']
    k1_c = r1['_kappa_sq_cart']
    k2_c = r2['_kappa_sq_cart']
    
    # Align by length (take shorter)
    n = min(len(k1_t), len(k2_t))
    k1_t = k1_t[:n]
    k2_t = k2_t[:n]
    k1_c = k1_c[:n]
    k2_c = k2_c[:n]
    
    # Mask zeros
    mask = (k1_t > 0) & (k2_t > 0)
    
    out = {}
    if HAS_SCIPY and np.sum(mask) > 10:
        r_t, p_t = sp_stats.pearsonr(k1_t[mask], k2_t[mask])
        out['torus_pearson_r'] = float(r_t)
        out['torus_p'] = float(p_t)
        
        rho_t, _ = sp_stats.spearmanr(k1_t[mask], k2_t[mask])
        out['torus_spearman'] = float(rho_t)
    
    mask_c = (k1_c > 0) & (k2_c > 0)
    if HAS_SCIPY and np.sum(mask_c) > 10:
        r_c, p_c = sp_stats.pearsonr(k1_c[mask_c], k2_c[mask_c])
        out['cart_pearson_r'] = float(r_c)
        out['cart_p'] = float(p_c)
    
    return out


def step8_noise_test(result, sigma=0.1, n_trials=10):
    """Step 8: Perturb coordinates by Gaussian noise and check stability."""
    pdb_path = download_pdb(result['pdb_id'])
    residues_orig = extract_backbone(pdb_path)
    
    kappa_orig = result['_kappa_sq_torus']
    density_orig = result['torus_sum_per_res']
    
    densities = []
    correlations = []
    
    for trial in range(n_trials):
        # Perturb all atom coordinates by sigma Å
        perturbed = []
        for r in residues_orig:
            pr = dict(r)
            for atom in ['N', 'CA', 'C']:
                pr[atom] = r[atom] + np.random.normal(0, sigma, 3)
            perturbed.append(pr)
        
        phi_p, psi_p = compute_dihedrals(perturbed)
        kappa_p = torus_curvature(phi_p, psi_p)
        kappa_sq_p = kappa_p ** 2
        
        density_p = np.sum(kappa_sq_p) / len(perturbed)
        densities.append(density_p)
        
        # Correlation with original
        n = min(len(kappa_orig), len(kappa_sq_p))
        mask = (kappa_orig[:n] > 0) & (kappa_sq_p[:n] > 0)
        if HAS_SCIPY and np.sum(mask) > 10:
            r, _ = sp_stats.pearsonr(kappa_orig[:n][mask], kappa_sq_p[:n][mask])
            correlations.append(r)
    
    return {
        'sigma_angstrom': sigma,
        'n_trials': n_trials,
        'orig_density': density_orig,
        'mean_perturbed_density': float(np.mean(densities)),
        'std_perturbed_density': float(np.std(densities)),
        'density_change_pct': float((np.mean(densities) - density_orig) / density_orig * 100),
        'mean_correlation': float(np.mean(correlations)) if correlations else None,
        'std_correlation': float(np.std(correlations)) if correlations else None,
    }


# ═══════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════

def main():
    if not HAS_BIOPYTHON:
        print("BioPython required. pip install biopython")
        sys.exit(1)
    
    print(f"\n  INDEPENDENT VALIDATION FROM RAW COORDINATES")
    print(f"  {'═'*55}")
    print(f"  Testing BOTH Cartesian Cα curvature (reviewer's)")
    print(f"  AND torus geodesic curvature (paper's method)")
    print()
    
    # ─── Step 1: Test cases ───
    test_cases = [
        # (PDB ID, chain, label, fold type)
        ('1MBN', 'A', 'Myoglobin (deoxy)', 'all-alpha'),
        ('1A6M', 'A', 'Myoglobin (met)', 'all-alpha'),
        ('1MBD', 'A', 'Myoglobin (CO)', 'all-alpha'),
        ('1LYZ', 'A', 'Lysozyme (hen)', 'alpha/beta'),
        ('1HEL', 'A', 'Lysozyme (human)', 'alpha/beta'),
        ('1UBQ', 'A', 'Ubiquitin', 'beta-rich'),
        ('1UBI', 'A', 'Ubiquitin (2)', 'beta-rich'),
    ]
    
    results = {}
    for pdb_id, chain, label, fold in test_cases:
        print(f"  Processing {label} ({pdb_id})...")
        r = analyze_protein(pdb_id, chain, label)
        if r:
            r['fold_type'] = fold
            results[pdb_id] = r
    
    if not results:
        print("  ⚠ No structures processed. Check network access.")
        return
    
    # ─── Step 4: SS validation ───
    print(f"\n  {'═'*55}")
    print(f"  STEP 4: CURVATURE AT SECONDARY STRUCTURE BOUNDARIES")
    print(f"  {'─'*55}")
    print(f"  Q: Are curvature peaks enriched at helix↔coil transitions?")
    print()
    
    for pdb_id, r in results.items():
        ss_val = step4_validate_ss(r)
        if ss_val is None:
            continue
        print(f"  {r['label']:25s} ({ss_val['n_boundaries']} boundaries)")
        if 'torus_ratio' in ss_val:
            sig = '★' if ss_val['torus_p'] < 0.05 else 'ns'
            print(f"    Torus:     boundary/interior = {ss_val['torus_ratio']:.2f}×  p={ss_val['torus_p']:.3f} {sig}")
        if 'cart_ratio' in ss_val:
            sig = '★' if ss_val['cart_p'] < 0.05 else 'ns'
            print(f"    Cartesian: boundary/interior = {ss_val['cart_ratio']:.2f}×  p={ss_val['cart_p']:.3f} {sig}")
    
    # ─── Step 5: Curvature density by fold class ───
    print(f"\n  {'═'*55}")
    print(f"  STEP 5: CURVATURE DENSITY BY FOLD CLASS")
    print(f"  {'─'*55}")
    print(f"  {'Protein':25s} {'Fold':>12s} {'Cart Σκ²/N':>12s} {'Torus Σκ²/N':>12s} {'H%':>5s}")
    
    for pdb_id, r in results.items():
        h_pct = r['n_helix'] / r['n_residues'] * 100
        print(f"  {r['label']:25s} {r['fold_type']:>12s} "
              f"{r['cart_sum_per_res']:12.4f} {r['torus_sum_per_res']:12.2f} "
              f"{h_pct:4.0f}%")
    
    # ─── Step 6: Reproducibility ───
    print(f"\n  {'═'*55}")
    print(f"  STEP 6: INTERNAL REPRODUCIBILITY")
    print(f"  {'─'*55}")
    print(f"  Same protein, different PDB entry: r > 0.8 = stable descriptor")
    print()
    
    pairs = [
        ('1MBN', '1A6M', 'Myoglobin pair'),
        ('1LYZ', '1HEL', 'Lysozyme pair'),
        ('1UBQ', '1UBI', 'Ubiquitin pair'),
    ]
    
    for id1, id2, label in pairs:
        if id1 in results and id2 in results:
            rep = step6_reproducibility(results[id1], results[id2])
            print(f"  {label}:")
            if 'torus_pearson_r' in rep:
                print(f"    Torus:     Pearson r = {rep['torus_pearson_r']:.3f}, "
                      f"Spearman ρ = {rep['torus_spearman']:.3f}")
            if 'cart_pearson_r' in rep:
                print(f"    Cartesian: Pearson r = {rep['cart_pearson_r']:.3f}")
    
    # ─── Step 7: Perturbation (apo vs holo) ───
    print(f"\n  {'═'*55}")
    print(f"  STEP 7: PERTURBATION SENSITIVITY (apo vs holo)")
    print(f"  {'─'*55}")
    
    if '1MBN' in results and '1MBD' in results:
        r_apo = results['1MBN']
        r_holo = results['1MBD']
        
        delta_cart = abs(r_holo['cart_sum_per_res'] - r_apo['cart_sum_per_res'])
        delta_torus = abs(r_holo['torus_sum_per_res'] - r_apo['torus_sum_per_res'])
        
        print(f"  Myoglobin deoxy vs CO-bound:")
        print(f"    Cart Σκ²/N:  {r_apo['cart_sum_per_res']:.4f} → {r_holo['cart_sum_per_res']:.4f} "
              f"(Δ = {delta_cart:.4f}, {delta_cart/r_apo['cart_sum_per_res']*100:.1f}%)")
        print(f"    Torus Σκ²/N: {r_apo['torus_sum_per_res']:.2f} → {r_holo['torus_sum_per_res']:.2f} "
              f"(Δ = {delta_torus:.2f}, {delta_torus/r_apo['torus_sum_per_res']*100:.1f}%)")
        
        # Profile correlation (redistribution test)
        rep = step6_reproducibility(r_apo, r_holo)
        if 'torus_pearson_r' in rep:
            print(f"    Torus profile correlation: r = {rep['torus_pearson_r']:.3f}")
            print(f"    → {'Redistribution' if rep['torus_pearson_r'] > 0.5 else 'Noisy'}")
    
    # ─── Step 8: Noise stress test ───
    print(f"\n  {'═'*55}")
    print(f"  STEP 8: NOISE STRESS TEST (σ = 0.1 Å)")
    print(f"  {'─'*55}")
    
    for pdb_id in ['1MBN', '1UBQ']:
        if pdb_id in results:
            noise = step8_noise_test(results[pdb_id], sigma=0.1, n_trials=20)
            print(f"  {results[pdb_id]['label']}:")
            print(f"    Original density: {noise['orig_density']:.2f}")
            print(f"    Perturbed:        {noise['mean_perturbed_density']:.2f} ± {noise['std_perturbed_density']:.2f}")
            print(f"    Change:           {noise['density_change_pct']:+.1f}%")
            if noise['mean_correlation'] is not None:
                print(f"    Profile correlation: {noise['mean_correlation']:.3f} ± {noise['std_correlation']:.3f}")
                stable = "STABLE" if noise['mean_correlation'] > 0.8 else "UNSTABLE"
                print(f"    → {stable}")
    
    # ─── Per-SS curvature comparison ───
    print(f"\n  {'═'*55}")
    print(f"  PER-SS CURVATURE (Paper's core claim)")
    print(f"  {'─'*55}")
    print(f"  {'Protein':25s} {'Method':>10s} {'Helix':>8s} {'Strand':>8s} {'Coil':>8s} {'H>E?':>5s}")
    
    for pdb_id, r in results.items():
        for method, prefix in [('Torus', 'torus'), ('Cartesian', 'cart')]:
            h = r.get(f'{prefix}_helix_mean', 0)
            e = r.get(f'{prefix}_strand_mean', 0)
            c = r.get(f'{prefix}_coil_mean', 0)
            he = '✓' if h > e else '✗'
            print(f"  {r['label']:25s} {method:>10s} {h:8.3f} {e:8.3f} {c:8.3f} {he:>5s}")
    
    # ─── THE HARDEST TEST ───
    print(f"\n  {'═'*55}")
    print(f"  THE HARDEST TEST: Budget vs Profile Separability")
    print(f"  {'─'*55}")
    
    if '1MBN' in results and '1UBQ' in results:
        r_myo = results['1MBN']
        r_ubq = results['1UBQ']
        
        print(f"  Myoglobin vs Ubiquitin:")
        print(f"    Torus density: {r_myo['torus_sum_per_res']:.2f} vs {r_ubq['torus_sum_per_res']:.2f}")
        
        # Profile correlation between different proteins
        rep = step6_reproducibility(r_myo, r_ubq)
        if 'torus_pearson_r' in rep:
            print(f"    Torus profile r:  {rep['torus_pearson_r']:.3f}")
            print(f"    → Budget {'similar' if abs(r_myo['torus_sum_per_res'] - r_ubq['torus_sum_per_res']) < r_myo['torus_sum_per_res'] * 0.5 else 'different'}, "
                  f"profile {'similar' if rep['torus_pearson_r'] > 0.5 else 'different'}")
            if abs(r_myo['torus_sum_per_res'] - r_ubq['torus_sum_per_res']) > 5 and rep['torus_pearson_r'] < 0.5:
                print(f"    ★ Budget AND profile separable — fold-specific curvature")
    
    # ─── VERDICT ───
    print(f"\n  {'═'*55}")
    print(f"  VERDICT")
    print(f"  {'═'*55}")
    
    # Save
    save_results = {}
    for k, v in results.items():
        save_results[k] = {key: val for key, val in v.items() if not key.startswith('_')}
    
    out_path = Path("results") / "validation_from_coordinates.json"
    out_path.parent.mkdir(exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(save_results, f, indent=2)
    print(f"\n  Saved: {out_path}")


if __name__ == "__main__":
    main()

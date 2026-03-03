#!/usr/bin/env python3
"""
Calmodulin Apo→Holo Gini Test
==============================
Computes kappa-squared curvature and Gini coefficient for calmodulin
in multiple conformational states to test the prediction:

    Apo (closed, rigid) → higher Gini (curvature concentrated at EF-hand hinges)
    Holo (open, flexible) → lower Gini (curvature redistributed as helices reorient)

This is a three-way validation test:
    1. Gini direction (apo > holo?)
    2. Consistency with NMR S² data (apo S²~0.85, holo S²~0.72)
    3. Consistency with atlas gradient (high-baseline proteins loosen under perturbation)

PDB structures:
    Apo:     1CFD (NMR, Kuboniwa/Tjandra 1995) — calcium-free CaM
    Holo:    1CLL (X-ray, 1.7 Å, Chattopadhyaya 1992) — Ca²⁺-saturated CaM
    Complex: 1CDL (X-ray, Meador 1992) — Ca²⁺-CaM + MLCK peptide

Methodology: Identical to atlas Papers 2/3 — discrete Frenet curvature on the
Ramachandran torus T², kappa-squared per residue, Gini coefficient.

Usage:
    python calmodulin_gini_test.py

    Requires: numpy, urllib (stdlib). No other dependencies.
    Downloads PDB files from RCSB automatically.

Author: Kase / True North Construction LLC
License: CC0 1.0 Universal (unpatentable prior art)
"""

import numpy as np
import os
import sys
import json
from urllib.request import urlretrieve
from pathlib import Path
from collections import OrderedDict

# ═══════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════

STRUCTURES = OrderedDict([
    ("1CFD", {
        "label": "Apo CaM (NMR, Ca²⁺-free)",
        "chain": "A",
        "is_nmr": True,
        "expected": "HIGH Gini (closed, rigid, curvature concentrated)",
    }),
    ("1CLL", {
        "label": "Holo CaM (X-ray, Ca²⁺-saturated)",
        "chain": "A",
        "is_nmr": False,
        "expected": "LOWER Gini (open, curvature redistributed)",
    }),
    ("1CDL", {
        "label": "CaM-MLCK Complex (X-ray, Ca²⁺ + target peptide)",
        "chain": "A",
        "is_nmr": False,
        "expected": "Variable — peptide may re-concentrate curvature",
    }),
])

# Additional structures for robustness check
EXTRA_STRUCTURES = OrderedDict([
    ("1QX5", {
        "label": "Apo CaM (X-ray, Ca²⁺-free crystal)",
        "chain": "A",
        "is_nmr": False,
        "expected": "HIGH Gini (cross-validates 1CFD with X-ray method)",
    }),
    ("1X02", {
        "label": "Holo CaM (NMR, Ca²⁺-saturated)",
        "chain": "A",
        "is_nmr": True,
        "expected": "LOWER Gini (cross-validates 1CLL with NMR method)",
    }),
])

DATA_DIR = Path("calmodulin_pdb")

# ═══════════════════════════════════════════════════════════════════
# PDB DOWNLOAD
# ═══════════════════════════════════════════════════════════════════

def download_pdb(pdb_id, data_dir=DATA_DIR):
    """Download PDB file from RCSB if not already present."""
    data_dir.mkdir(exist_ok=True)
    filepath = data_dir / f"{pdb_id.upper()}.pdb"
    if filepath.exists():
        print(f"  [cached] {filepath}")
        return filepath
    
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"  Downloading {url} ...")
    try:
        urlretrieve(url, filepath)
        print(f"  [ok] {filepath}")
    except Exception as e:
        print(f"  [FAIL] Could not download {pdb_id}: {e}")
        return None
    return filepath


# ═══════════════════════════════════════════════════════════════════
# PDB PARSING
# ═══════════════════════════════════════════════════════════════════

def parse_backbone_atoms(filepath, chain_id="A", model_num=1):
    """
    Parse backbone N, CA, C atoms from a PDB file.
    For NMR ensembles, takes the specified model (default: model 1).
    Returns dict: {resnum: {'N': [x,y,z], 'CA': [x,y,z], 'C': [x,y,z]}}
    """
    residues = {}
    current_model = 0
    in_target_model = False
    
    with open(filepath) as f:
        for line in f:
            if line.startswith("MODEL"):
                current_model = int(line[10:14].strip())
                in_target_model = (current_model == model_num)
                continue
            if line.startswith("ENDMDL"):
                if in_target_model:
                    break  # done with our model
                in_target_model = False
                continue
            
            # If no MODEL records (single-model X-ray), always parse
            if current_model == 0:
                in_target_model = True
            
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if not in_target_model:
                continue
            
            atom_name = line[12:16].strip()
            if atom_name not in ("N", "CA", "C"):
                continue
            
            chain = line[21].strip()
            if chain != chain_id:
                continue
            
            # Handle insertion codes
            resnum_str = line[22:27].strip()
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                continue
            
            icode = line[26].strip()
            if icode:
                continue  # skip insertion-code residues
            
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            
            if resnum not in residues:
                residues[resnum] = {}
            residues[resnum][atom_name] = np.array([x, y, z])
    
    return residues


# ═══════════════════════════════════════════════════════════════════
# DIHEDRAL ANGLE COMPUTATION
# ═══════════════════════════════════════════════════════════════════

def dihedral_angle(p1, p2, p3, p4):
    """
    Compute dihedral angle (in radians) for four points.
    Standard Praxeolitic method.
    """
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    
    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    
    if n1_norm < 1e-10 or n2_norm < 1e-10:
        return np.nan
    
    n1 = n1 / n1_norm
    n2 = n2 / n2_norm
    
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    
    return np.arctan2(y, x)


def compute_dihedrals(residues):
    """
    Compute phi and psi for each residue.
    phi(i) = dihedral(C[i-1], N[i], CA[i], C[i])
    psi(i) = dihedral(N[i], CA[i], C[i], N[i+1])
    Returns: {resnum: (phi, psi)} — angles in radians
    """
    sorted_resnums = sorted(residues.keys())
    dihedrals = {}
    
    for idx, rn in enumerate(sorted_resnums):
        r = residues[rn]
        if not all(k in r for k in ("N", "CA", "C")):
            continue
        
        # phi: C(i-1) - N(i) - CA(i) - C(i)
        phi = np.nan
        if idx > 0:
            prev_rn = sorted_resnums[idx - 1]
            prev = residues.get(prev_rn, {})
            if "C" in prev and (rn - prev_rn) == 1:  # sequential check
                phi = dihedral_angle(prev["C"], r["N"], r["CA"], r["C"])
        
        # psi: N(i) - CA(i) - C(i) - N(i+1)
        psi = np.nan
        if idx < len(sorted_resnums) - 1:
            next_rn = sorted_resnums[idx + 1]
            nxt = residues.get(next_rn, {})
            if "N" in nxt and (next_rn - rn) == 1:  # sequential check
                psi = dihedral_angle(r["N"], r["CA"], r["C"], nxt["N"])
        
        dihedrals[rn] = (phi, psi)
    
    return dihedrals


# ═══════════════════════════════════════════════════════════════════
# KAPPA-SQUARED (Discrete Frenet curvature on Ramachandran torus)
# ═══════════════════════════════════════════════════════════════════

def angular_diff(a, b):
    """Signed angular difference with proper wrapping on [-pi, pi]."""
    return np.arctan2(np.sin(a - b), np.cos(a - b))


def compute_kappa_squared(dihedrals):
    """
    Compute kappa-squared per residue using discrete Frenet frame on T².
    
    For residue i, using three consecutive residues (i-1, i, i+1):
        tangent₁ = (dphi₁, dpsi₁) / ||(dphi₁, dpsi₁)||
        tangent₂ = (dphi₂, dpsi₂) / ||(dphi₂, dpsi₂)||
        kappa = |tangent₂ - tangent₁| / mean_arc_length
        kappa² = kappa²
    
    This is curvature in CONFORMATIONAL space (phi, psi), 
    not spatial curvature of the 3D chain.
    """
    sorted_resnums = sorted(dihedrals.keys())
    kappa_sq = {}
    
    for idx in range(1, len(sorted_resnums) - 1):
        rn_prev = sorted_resnums[idx - 1]
        rn = sorted_resnums[idx]
        rn_next = sorted_resnums[idx + 1]
        
        # Sequential check
        if (rn - rn_prev) != 1 or (rn_next - rn) != 1:
            continue
        
        phi_prev, psi_prev = dihedrals[rn_prev]
        phi_curr, psi_curr = dihedrals[rn]
        phi_next, psi_next = dihedrals[rn_next]
        
        # Skip if any angle is nan
        if any(np.isnan(x) for x in [phi_prev, psi_prev, phi_curr, psi_curr, 
                                       phi_next, psi_next]):
            continue
        
        # Step 1: dihedral differences with angular wrapping
        dphi1 = angular_diff(phi_curr, phi_prev)
        dpsi1 = angular_diff(psi_curr, psi_prev)
        dphi2 = angular_diff(phi_next, phi_curr)
        dpsi2 = angular_diff(psi_next, psi_curr)
        
        # Step 2: arc lengths
        ds1 = np.sqrt(dphi1**2 + dpsi1**2)
        ds2 = np.sqrt(dphi2**2 + dpsi2**2)
        
        if ds1 < 1e-10 or ds2 < 1e-10:
            continue
        
        # Step 3: unit tangent vectors
        T1_phi = dphi1 / ds1
        T1_psi = dpsi1 / ds1
        T2_phi = dphi2 / ds2
        T2_psi = dpsi2 / ds2
        
        # Step 4: curvature = |dT/ds|
        dT_phi = T2_phi - T1_phi
        dT_psi = T2_psi - T1_psi
        ds_mid = 0.5 * (ds1 + ds2)
        
        kappa_sq[rn] = (dT_phi**2 + dT_psi**2) / (ds_mid**2)
    
    return kappa_sq


# ═══════════════════════════════════════════════════════════════════
# GINI COEFFICIENT
# ═══════════════════════════════════════════════════════════════════

def gini_coefficient(values):
    """
    Compute Gini coefficient of a distribution.
    Gini = 0: perfectly equal. Gini = 1: maximally concentrated.
    """
    v = np.array(values, dtype=float)
    v = v[~np.isnan(v)]
    if len(v) == 0:
        return np.nan
    v = np.sort(v)
    n = len(v)
    total = np.sum(v)
    if total < 1e-15:
        return 0.0
    cumsum = np.cumsum(v)
    gini = (2.0 * np.sum((np.arange(1, n + 1) * v)) / (n * total)) - (n + 1) / n
    return gini


# ═══════════════════════════════════════════════════════════════════
# RAMACHANDRAN CLASSIFICATION
# ═══════════════════════════════════════════════════════════════════

def rama_class(phi, psi):
    """Classify residue by Ramachandran region."""
    if np.isnan(phi) or np.isnan(psi):
        return "unk"
    phi_d = np.degrees(phi)
    psi_d = np.degrees(psi)
    if -180 <= phi_d <= -30 and -80 <= psi_d <= -10:
        return "alpha"
    if -180 <= phi_d <= -30 and 80 <= psi_d <= 180:
        return "beta"
    if -180 <= phi_d <= -30 and -10 <= psi_d <= 80:
        return "beta"  # extended beta
    if 30 <= phi_d <= 120 and -60 <= psi_d <= 60:
        return "alphaL"
    return "coil"


# ═══════════════════════════════════════════════════════════════════
# PER-RESIDUE ANALYSIS
# ═══════════════════════════════════════════════════════════════════

def analyze_structure(filepath, chain_id, label, is_nmr=False):
    """Full analysis pipeline for one structure."""
    
    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"  File: {filepath}  Chain: {chain_id}")
    print(f"{'='*70}")
    
    # For NMR ensembles, analyze model 1 and report ensemble stats
    models_to_analyze = [1]
    if is_nmr:
        # Check how many models exist
        model_count = 0
        with open(filepath) as f:
            for line in f:
                if line.startswith("MODEL"):
                    model_count += 1
        if model_count > 1:
            models_to_analyze = list(range(1, min(model_count + 1, 21)))  # up to 20 models
            print(f"  NMR ensemble: {model_count} models (analyzing {len(models_to_analyze)})")
    
    all_ginis = []
    all_kappa_sq_means = []
    primary_result = None
    
    for model_num in models_to_analyze:
        residues = parse_backbone_atoms(filepath, chain_id, model_num)
        
        if len(residues) < 10:
            print(f"  [WARN] Model {model_num}: only {len(residues)} residues parsed, skipping")
            continue
        
        dihedrals = compute_dihedrals(residues)
        kappa_sq = compute_kappa_squared(dihedrals)
        
        if len(kappa_sq) < 5:
            print(f"  [WARN] Model {model_num}: only {len(kappa_sq)} kappa² values, skipping")
            continue
        
        vals = list(kappa_sq.values())
        g = gini_coefficient(vals)
        all_ginis.append(g)
        all_kappa_sq_means.append(np.mean(vals))
        
        if model_num == 1:
            primary_result = {
                "dihedrals": dihedrals,
                "kappa_sq": kappa_sq,
                "gini": g,
                "n_residues": len(residues),
                "n_kappa": len(kappa_sq),
            }
    
    if primary_result is None:
        print("  [ERROR] No valid models found")
        return None
    
    # Report primary model
    ksq_vals = list(primary_result["kappa_sq"].values())
    print(f"\n  Residues parsed:  {primary_result['n_residues']}")
    print(f"  Kappa² computed:  {primary_result['n_kappa']}")
    print(f"  Mean kappa²:      {np.mean(ksq_vals):.3f}")
    print(f"  Median kappa²:    {np.median(ksq_vals):.3f}")
    print(f"  Max kappa²:       {np.max(ksq_vals):.3f}")
    print(f"  Sum kappa²:       {np.sum(ksq_vals):.1f}")
    print(f"\n  *** GINI = {primary_result['gini']:.4f} ***")
    
    # NMR ensemble statistics
    if len(all_ginis) > 1:
        print(f"\n  NMR Ensemble ({len(all_ginis)} models):")
        print(f"    Gini mean:   {np.mean(all_ginis):.4f}")
        print(f"    Gini SD:     {np.std(all_ginis):.4f}")
        print(f"    Gini range:  [{min(all_ginis):.4f}, {max(all_ginis):.4f}]")
        primary_result["ensemble_ginis"] = all_ginis
        primary_result["ensemble_gini_mean"] = np.mean(all_ginis)
        primary_result["ensemble_gini_sd"] = np.std(all_ginis)
    
    # Secondary structure breakdown
    dihedrals = primary_result["dihedrals"]
    kappa_sq = primary_result["kappa_sq"]
    
    ss_bins = {"alpha": [], "beta": [], "coil": [], "alphaL": []}
    for rn, ksq in kappa_sq.items():
        if rn in dihedrals:
            phi, psi = dihedrals[rn]
            ss = rama_class(phi, psi)
            if ss in ss_bins:
                ss_bins[ss].append(ksq)
    
    print(f"\n  Secondary structure breakdown:")
    for ss in ["alpha", "beta", "coil", "alphaL"]:
        vals = ss_bins[ss]
        if vals:
            print(f"    {ss:>6s}: n={len(vals):>3d}  mean_κ²={np.mean(vals):.3f}  "
                  f"sum_κ²={np.sum(vals):.1f}  ({100*len(vals)/len(kappa_sq):.0f}%)")
    
    # Top 10 curvature hotspots
    sorted_res = sorted(kappa_sq.items(), key=lambda x: -x[1])
    print(f"\n  Top 10 curvature hotspots:")
    for rn, ksq in sorted_res[:10]:
        phi, psi = dihedrals.get(rn, (np.nan, np.nan))
        ss = rama_class(phi, psi)
        print(f"    Res {rn:>4d} ({ss:>5s}): κ² = {ksq:.3f}")
    
    return primary_result


# ═══════════════════════════════════════════════════════════════════
# COMPARISON
# ═══════════════════════════════════════════════════════════════════

def compare_structures(result1, label1, result2, label2):
    """Compare two structures: delta-Gini and per-residue differences."""
    
    print(f"\n{'─'*70}")
    print(f"  COMPARISON: {label1}  →  {label2}")
    print(f"{'─'*70}")
    
    g1 = result1["gini"]
    g2 = result2["gini"]
    dg = g2 - g1
    
    print(f"\n  Gini({label1[:20]:>20s}) = {g1:.4f}")
    print(f"  Gini({label2[:20]:>20s}) = {g2:.4f}")
    print(f"  ΔGini = {dg:+.4f}  ({'STIFFENING' if dg > 0 else 'LOOSENING'})")
    
    # Per-residue delta-kappa-squared
    ksq1 = result1["kappa_sq"]
    ksq2 = result2["kappa_sq"]
    common = sorted(set(ksq1.keys()) & set(ksq2.keys()))
    
    if len(common) < 5:
        print(f"  [WARN] Only {len(common)} common residues — structures may not align")
        return
    
    deltas = {rn: ksq2[rn] - ksq1[rn] for rn in common}
    delta_vals = list(deltas.values())
    
    print(f"\n  Common residues: {len(common)}")
    print(f"  Mean |Δκ²|/res:  {np.mean(np.abs(delta_vals)):.3f}")
    print(f"  % perturbed (|Δκ²| > 5): {100*np.mean([abs(d)>5 for d in delta_vals]):.1f}%")
    
    # Biggest changes
    sorted_deltas = sorted(deltas.items(), key=lambda x: -abs(x[1]))
    print(f"\n  Top 10 largest |Δκ²| residues:")
    for rn, d in sorted_deltas[:10]:
        direction = "↑" if d > 0 else "↓"
        print(f"    Res {rn:>4d}: Δκ² = {d:+.3f} {direction}")
    
    # How many residues stiffened vs loosened
    n_stiff = sum(1 for d in delta_vals if d > 0.01)
    n_loose = sum(1 for d in delta_vals if d < -0.01)
    n_neutral = len(delta_vals) - n_stiff - n_loose
    print(f"\n  Residue-level: {n_stiff} stiffened, {n_loose} loosened, {n_neutral} neutral")
    
    # If we have ensemble data, compare distributions
    if "ensemble_ginis" in result1 and "ensemble_ginis" in result2:
        g1s = result1["ensemble_ginis"]
        g2s = result2["ensemble_ginis"]
        # Simple t-test
        from scipy import stats
        try:
            t_stat, p_val = stats.ttest_ind(g1s, g2s)
            print(f"\n  NMR Ensemble comparison:")
            print(f"    {label1}: Gini = {np.mean(g1s):.4f} ± {np.std(g1s):.4f} (n={len(g1s)})")
            print(f"    {label2}: Gini = {np.mean(g2s):.4f} ± {np.std(g2s):.4f} (n={len(g2s)})")
            print(f"    t = {t_stat:.2f}, p = {p_val:.4f}")
        except ImportError:
            print(f"\n  (scipy not available — skipping ensemble t-test)")


# ═══════════════════════════════════════════════════════════════════
# NMR HOTSPOT VALIDATION
# ═══════════════════════════════════════════════════════════════════

def validate_against_nmr(result, label):
    """
    Check if the curvature hotspot at Phe-92 (from Gemini's table)
    is real and where it ranks.
    """
    kappa_sq = result["kappa_sq"]
    
    print(f"\n  NMR Validation Check ({label}):")
    
    # Gemini claims Phe-92 is the κ²_max hotspot for calmodulin
    target_residue = 92
    if target_residue in kappa_sq:
        ksq_92 = kappa_sq[target_residue]
        rank = sum(1 for v in kappa_sq.values() if v > ksq_92) + 1
        percentile = 100 * (1 - rank / len(kappa_sq))
        print(f"    Res 92 (Phe): κ² = {ksq_92:.3f}  "
              f"(rank {rank}/{len(kappa_sq)}, {percentile:.0f}th percentile)")
    else:
        print(f"    Res 92 not found in kappa² data")
        # Find nearest
        nearby = [rn for rn in kappa_sq if abs(rn - 92) <= 3]
        for rn in nearby:
            print(f"    Nearby: Res {rn}: κ² = {kappa_sq[rn]:.3f}")


# ═══════════════════════════════════════════════════════════════════
# DOMAIN-LEVEL ANALYSIS
# ═══════════════════════════════════════════════════════════════════

def domain_analysis(result, label):
    """
    Calmodulin has two domains connected by a flexible linker:
        N-terminal domain: ~res 4-77
        Linker: ~res 77-81  
        C-terminal domain: ~res 82-148
    
    Compute Gini separately for each domain.
    """
    kappa_sq = result["kappa_sq"]
    
    n_term = {rn: v for rn, v in kappa_sq.items() if 4 <= rn <= 76}
    linker = {rn: v for rn, v in kappa_sq.items() if 77 <= rn <= 81}
    c_term = {rn: v for rn, v in kappa_sq.items() if 82 <= rn <= 148}
    
    print(f"\n  Domain-level analysis ({label}):")
    for name, domain in [("N-term (4-76)", n_term), 
                          ("Linker (77-81)", linker),
                          ("C-term (82-148)", c_term)]:
        if len(domain) < 3:
            print(f"    {name:>20s}: too few residues ({len(domain)})")
            continue
        vals = list(domain.values())
        g = gini_coefficient(vals)
        print(f"    {name:>20s}: n={len(domain):>3d}  Gini={g:.4f}  "
              f"mean_κ²={np.mean(vals):.3f}  sum_κ²={np.sum(vals):.1f}")


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║     CALMODULIN APO → HOLO GINI TEST                           ║")
    print("║     Section 18.5: Three-Way Dynamic Validation                 ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    
    print("\nPREDICTION:")
    print("  Apo CaM (closed, rigid, S²~0.85) → HIGHER Gini")
    print("  Holo CaM (open, flexible, S²~0.72) → LOWER Gini")
    print("  If confirmed: Gini captures the same physics as NMR order parameters")
    print("  Direction: this is a LOOSENING transition (ΔGini < 0)")
    
    # Download structures
    print("\n── Downloading PDB files ──")
    all_structures = {**STRUCTURES, **EXTRA_STRUCTURES}
    paths = {}
    for pdb_id in all_structures:
        path = download_pdb(pdb_id)
        if path:
            paths[pdb_id] = path
    
    if len(paths) < 2:
        print("\n[FATAL] Need at least 2 structures. Check network connection.")
        sys.exit(1)
    
    # Analyze each structure
    print("\n── Analyzing structures ──")
    results = {}
    for pdb_id, info in all_structures.items():
        if pdb_id not in paths:
            continue
        result = analyze_structure(
            paths[pdb_id], 
            info["chain"], 
            f"{pdb_id}: {info['label']}",
            info.get("is_nmr", False)
        )
        if result:
            results[pdb_id] = result
            validate_against_nmr(result, pdb_id)
            domain_analysis(result, pdb_id)
    
    # Comparisons
    print("\n\n" + "═"*70)
    print("  COMPARISONS")
    print("═"*70)
    
    # Primary comparison: Apo (1CFD) vs Holo (1CLL)
    if "1CFD" in results and "1CLL" in results:
        compare_structures(results["1CFD"], "Apo (1CFD)", 
                          results["1CLL"], "Holo (1CLL)")
    
    # Holo vs Complex
    if "1CLL" in results and "1CDL" in results:
        compare_structures(results["1CLL"], "Holo (1CLL)", 
                          results["1CDL"], "Complex (1CDL)")
    
    # Apo vs Complex (full transition)
    if "1CFD" in results and "1CDL" in results:
        compare_structures(results["1CFD"], "Apo (1CFD)", 
                          results["1CDL"], "Complex (1CDL)")
    
    # Cross-method validation
    if "1QX5" in results and "1CFD" in results:
        compare_structures(results["1QX5"], "Apo X-ray (1QX5)", 
                          results["1CFD"], "Apo NMR (1CFD)")
    
    if "1X02" in results and "1CLL" in results:
        compare_structures(results["1X02"], "Holo NMR (1X02)", 
                          results["1CLL"], "Holo X-ray (1CLL)")
    
    # Summary
    print("\n\n" + "╔"+"═"*68+"╗")
    print("║  SUMMARY" + " "*59 + "║")
    print("╠"+"═"*68+"╣")
    
    summary_order = ["1CFD", "1QX5", "1CLL", "1X02", "1CDL"]
    for pdb_id in summary_order:
        if pdb_id in results:
            g = results[pdb_id]["gini"]
            info = all_structures[pdb_id]
            ens_str = ""
            if "ensemble_gini_mean" in results[pdb_id]:
                ens_str = f" (ensemble: {results[pdb_id]['ensemble_gini_mean']:.4f} ± {results[pdb_id]['ensemble_gini_sd']:.4f})"
            print(f"║  {pdb_id} {info['label'][:40]:40s} Gini = {g:.4f}{ens_str}")
    
    print("╠"+"═"*68+"╣")
    
    # Verdict
    if "1CFD" in results and "1CLL" in results:
        g_apo = results["1CFD"]["gini"]
        g_holo = results["1CLL"]["gini"]
        dg = g_holo - g_apo
        
        if dg < 0:
            print("║  ✓ PREDICTION CONFIRMED: Apo Gini > Holo Gini (LOOSENING)      ║")
            print(f"║    ΔGini = {dg:+.4f}                                              ║")
            print("║    Consistent with NMR: S² drops from ~0.85 to ~0.72            ║")
            print("║    → Gini captures the same physics as order parameters          ║")
        elif dg > 0:
            print("║  ✗ PREDICTION FAILED: Holo Gini > Apo Gini (STIFFENING)         ║")
            print(f"║    ΔGini = {dg:+.4f}                                              ║")
            print("║    This contradicts the NMR data. Investigate why.               ║")
        else:
            print("║  ~ INCONCLUSIVE: ΔGini ≈ 0                                      ║")
    
    print("╚"+"═"*68+"╝")
    
    # Export per-residue data for further analysis
    export_path = DATA_DIR / "calmodulin_results.json"
    export = {}
    for pdb_id, result in results.items():
        export[pdb_id] = {
            "gini": result["gini"],
            "n_residues": result["n_residues"],
            "n_kappa": result["n_kappa"],
            "kappa_sq": {str(k): v for k, v in result["kappa_sq"].items()},
        }
        if "ensemble_ginis" in result:
            export[pdb_id]["ensemble_ginis"] = result["ensemble_ginis"]
    
    with open(export_path, "w") as f:
        json.dump(export, f, indent=2)
    print(f"\nPer-residue data exported to: {export_path}")
    print("Load with: json.load(open('calmodulin_pdb/calmodulin_results.json'))")


if __name__ == "__main__":
    main()

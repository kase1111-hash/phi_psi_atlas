#!/usr/bin/env python3
"""
Baseline Gradient Test: Four NMR-Characterized Proteins
========================================================
Tests the atlas prediction: ΔGini direction is determined by baseline Gini
relative to the ~0.81 attractor.

    Baseline < 0.81 → perturbation STIFFENS (ΔGini > 0)
    Baseline > 0.81 → perturbation LOOSENS  (ΔGini < 0)

Proteins (from Gemini's Section 18 NMR table):
    1. Calmodulin   — apo vs Ca²⁺-bound (already tested: ΔGini = +0.09)
    2. Ubiquitin     — apo vs partner-bound
    3. GB1 domain    — apo vs IgG Fc-bound
    4. Lysozyme HEWL — apo vs inhibitor-bound

For each protein we compute:
    - Baseline Gini (apo structure)
    - Perturbed Gini (bound structure, same chain)
    - ΔGini and direction
    - Whether direction matches gradient prediction

Usage:
    python gradient_test_4proteins.py

    Requires: numpy. Downloads PDB files from RCSB automatically.

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

ATTRACTOR = 0.81  # proteome-wide Gini attractor from atlas

PROTEINS = OrderedDict([
    # ── CALMODULIN ──
    ("calmodulin", {
        "apo": {
            "pdb": "1CFD", "chain": "A", "is_nmr": True,
            "label": "Apo CaM (NMR, Ca²⁺-free)",
        },
        "bound": {
            "pdb": "1CLL", "chain": "A", "is_nmr": False,
            "label": "Holo CaM (X-ray, Ca²⁺-saturated)",
            "perturbation": "Ca²⁺ binding (4 ions)",
        },
        "extra_bound": {
            "pdb": "1CDL", "chain": "A", "is_nmr": False,
            "label": "CaM-MLCK Complex",
            "perturbation": "Ca²⁺ + MLCK peptide",
        },
        "nmr_S2": 0.74,  # Gemini's value for Phe-92 in apo
        "nmr_ref": "Barbato et al. 1992",
    }),

    # ── UBIQUITIN ──
    # Apo: 1UBQ (X-ray, 1.8 Å) — the gold standard ubiquitin structure
    # Bound: 1NBF chain C/D (ubiquitin bound to HAUSP deubiquitinase)
    #        or 1S1Q chain B (ubiquitin bound to TSG101 UEV domain, 2.0 Å)
    # Alt bound: 1AAR (K48-linked diubiquitin, closed conformation)
    ("ubiquitin", {
        "apo": {
            "pdb": "1UBQ", "chain": "A", "is_nmr": False,
            "label": "Apo Ubiquitin (X-ray, 1.8 Å)",
        },
        "bound": {
            "pdb": "1S1Q", "chain": "B", "is_nmr": False,
            "label": "Ub-TSG101 UEV Complex (2.0 Å)",
            "perturbation": "TSG101 UEV domain binding (Ile44 patch)",
        },
        "extra_bound": {
            "pdb": "1NBF", "chain": "C", "is_nmr": False,
            "label": "Ub-HAUSP Complex",
            "perturbation": "HAUSP deubiquitinase binding",
        },
        "nmr_S2": 0.84,  # Gemini's value for Thr-9
        "nmr_ref": "Tjandra et al. 1995 / Schneider et al. 1992",
    }),

    # ── GB1 DOMAIN ──
    # Apo: 1PGA (X-ray, 2.07 Å) — B1 IgG-binding domain of Protein G
    # Bound: 1FCC chain C (Protein G C2 fragment bound to IgG Fc)
    #   Note: 1FCC has chains A,B = IgG Fc and C,D = Protein G C2 fragment.
    #   C2 is a different domain (Fc-binding) than B1 (Fab-binding).
    #   For B1 specifically: 1IGC has GB1 (domain III) bound to Fab.
    ("gb1", {
        "apo": {
            "pdb": "1PGA", "chain": "A", "is_nmr": False,
            "label": "Apo GB1 (X-ray, 2.07 Å)",
        },
        "bound": {
            "pdb": "1PGB", "chain": "A", "is_nmr": False,
            "label": "GB1 alternate crystal form (1.92 Å)",
            "perturbation": "Crystal packing (different space group)",
        },
        "extra_bound": {
            "pdb": "1IGC", "chain": "C", "is_nmr": False,
            "label": "GB1-Fab Complex",
            "perturbation": "IgG Fab fragment binding",
        },
        "nmr_S2": 0.88,  # Gemini's value for Tyr-3
        "nmr_ref": "Barchi et al. 1994",
    }),

    # ── LYSOZYME (HEWL) ──
    # Apo: 1AKI (X-ray, 1.5 Å) — hen egg-white lysozyme, apo form
    # Bound: 1LZB (X-ray) — HEWL with bound tri-NAG inhibitor
    #        or 1HEW (X-ray, 1.65 Å) — HEWL with tri-NAG at 1.65 Å
    #        or 4LZT (X-ray, 0.95 Å) — ultra-high resolution
    ("lysozyme", {
        "apo": {
            "pdb": "1AKI", "chain": "A", "is_nmr": False,
            "label": "Apo HEWL (X-ray, 1.5 Å)",
        },
        "bound": {
            "pdb": "1HEW", "chain": "A", "is_nmr": False,
            "label": "HEWL + tri-NAG (X-ray, 1.65 Å)",
            "perturbation": "tri-NAG inhibitor binding",
        },
        "extra_bound": {
            "pdb": "1LZB", "chain": "A", "is_nmr": False,
            "label": "HEWL + tri-NAG (alt structure)",
            "perturbation": "tri-NAG inhibitor binding (alt)",
        },
        "nmr_S2": 0.81,  # Gemini's value for Gly-102
        "nmr_ref": "Buck et al. 1995",
    }),
])

DATA_DIR = Path("gradient_test_pdb")


# ═══════════════════════════════════════════════════════════════════
# PDB DOWNLOAD & PARSING (same as calmodulin script)
# ═══════════════════════════════════════════════════════════════════

def download_pdb(pdb_id, data_dir=DATA_DIR):
    data_dir.mkdir(exist_ok=True)
    filepath = data_dir / f"{pdb_id.upper()}.pdb"
    if filepath.exists():
        return filepath
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"    Downloading {pdb_id} ... ", end="", flush=True)
    try:
        urlretrieve(url, filepath)
        print("ok")
    except Exception as e:
        print(f"FAIL ({e})")
        return None
    return filepath


def parse_backbone_atoms(filepath, chain_id="A", model_num=1):
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
                if in_target_model: break
                in_target_model = False
                continue
            if current_model == 0:
                in_target_model = True
            if not line.startswith("ATOM"):
                continue
            if not in_target_model:
                continue
            atom_name = line[12:16].strip()
            if atom_name not in ("N", "CA", "C"):
                continue
            chain = line[21].strip()
            if chain != chain_id:
                continue
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                continue
            icode = line[26].strip()
            if icode:
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            if resnum not in residues:
                residues[resnum] = {}
            residues[resnum][atom_name] = np.array([x, y, z])
    return residues


# ═══════════════════════════════════════════════════════════════════
# DIHEDRAL & CURVATURE COMPUTATION (identical to atlas methodology)
# ═══════════════════════════════════════════════════════════════════

def dihedral_angle(p1, p2, p3, p4):
    b1 = p2 - p1; b2 = p3 - p2; b3 = p4 - p3
    n1 = np.cross(b1, b2); n2 = np.cross(b2, b3)
    n1n = np.linalg.norm(n1); n2n = np.linalg.norm(n2)
    if n1n < 1e-10 or n2n < 1e-10: return np.nan
    n1 /= n1n; n2 /= n2n
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    return np.arctan2(np.dot(m1, n2), np.dot(n1, n2))


def compute_dihedrals(residues):
    sorted_rn = sorted(residues.keys())
    dihedrals = {}
    for idx, rn in enumerate(sorted_rn):
        r = residues[rn]
        if not all(k in r for k in ("N", "CA", "C")):
            continue
        phi = np.nan
        if idx > 0:
            prev_rn = sorted_rn[idx - 1]
            prev = residues.get(prev_rn, {})
            if "C" in prev and (rn - prev_rn) == 1:
                phi = dihedral_angle(prev["C"], r["N"], r["CA"], r["C"])
        psi = np.nan
        if idx < len(sorted_rn) - 1:
            next_rn = sorted_rn[idx + 1]
            nxt = residues.get(next_rn, {})
            if "N" in nxt and (next_rn - rn) == 1:
                psi = dihedral_angle(r["N"], r["CA"], r["C"], nxt["N"])
        dihedrals[rn] = (phi, psi)
    return dihedrals


def angular_diff(a, b):
    return np.arctan2(np.sin(a - b), np.cos(a - b))


def compute_kappa_squared(dihedrals):
    sorted_rn = sorted(dihedrals.keys())
    kappa_sq = {}
    for idx in range(1, len(sorted_rn) - 1):
        rn_p = sorted_rn[idx - 1]; rn = sorted_rn[idx]; rn_n = sorted_rn[idx + 1]
        if (rn - rn_p) != 1 or (rn_n - rn) != 1:
            continue
        pp, sp = dihedrals[rn_p]; pc, sc = dihedrals[rn]; pn, sn = dihedrals[rn_n]
        if any(np.isnan(x) for x in [pp, sp, pc, sc, pn, sn]):
            continue
        dp1 = angular_diff(pc, pp); ds1 = angular_diff(sc, sp)
        dp2 = angular_diff(pn, pc); ds2 = angular_diff(sn, sc)
        a1 = np.sqrt(dp1**2 + ds1**2); a2 = np.sqrt(dp2**2 + ds2**2)
        if a1 < 1e-10 or a2 < 1e-10:
            continue
        dTp = dp2/a2 - dp1/a1; dTs = ds2/a2 - ds1/a1
        ds_mid = 0.5 * (a1 + a2)
        kappa_sq[rn] = (dTp**2 + dTs**2) / (ds_mid**2)
    return kappa_sq


def gini_coefficient(values):
    v = np.array(values, dtype=float)
    v = v[~np.isnan(v)]
    if len(v) == 0: return np.nan
    v = np.sort(v); n = len(v); total = np.sum(v)
    if total < 1e-15: return 0.0
    return (2.0 * np.sum(np.arange(1, n+1) * v) / (n * total)) - (n+1) / n


def rama_class(phi, psi):
    if np.isnan(phi) or np.isnan(psi): return "unk"
    pd = np.degrees(phi); sd = np.degrees(psi)
    if -180 <= pd <= -30 and -80 <= sd <= -10: return "alpha"
    if -180 <= pd <= -30 and (80 <= sd <= 180 or -10 <= sd <= 80): return "beta"
    if 30 <= pd <= 120 and -60 <= sd <= 60: return "alphaL"
    return "coil"


# ═══════════════════════════════════════════════════════════════════
# ANALYSIS PIPELINE
# ═══════════════════════════════════════════════════════════════════

def analyze_one(filepath, chain_id, label, is_nmr=False, verbose=True):
    """Full pipeline: parse → dihedrals → κ² → Gini. Returns result dict."""
    
    models = [1]
    if is_nmr:
        mc = 0
        with open(filepath) as f:
            for line in f:
                if line.startswith("MODEL"): mc += 1
        if mc > 1:
            models = list(range(1, min(mc + 1, 21)))
    
    all_ginis = []
    primary = None
    
    for m in models:
        res = parse_backbone_atoms(filepath, chain_id, m)
        if len(res) < 5: continue
        dih = compute_dihedrals(res)
        ksq = compute_kappa_squared(dih)
        if len(ksq) < 3: continue
        g = gini_coefficient(list(ksq.values()))
        all_ginis.append(g)
        if m == 1:
            primary = {
                "gini": g, "kappa_sq": ksq, "dihedrals": dih,
                "n_res": len(res), "n_kappa": len(ksq),
            }
    
    if primary is None:
        if verbose: print(f"    [ERROR] {label}: no valid models")
        return None
    
    vals = list(primary["kappa_sq"].values())
    
    if verbose:
        print(f"    {label}")
        print(f"      Residues: {primary['n_res']}  κ² computed: {primary['n_kappa']}")
        print(f"      Σκ² = {np.sum(vals):.0f}  mean = {np.mean(vals):.1f}  max = {np.max(vals):.0f}")
        print(f"      *** Gini = {primary['gini']:.4f} ***")
    
    if len(all_ginis) > 1:
        primary["ensemble_mean"] = np.mean(all_ginis)
        primary["ensemble_sd"] = np.std(all_ginis)
        primary["ensemble_n"] = len(all_ginis)
        if verbose:
            print(f"      NMR ensemble: {np.mean(all_ginis):.4f} ± {np.std(all_ginis):.4f} (n={len(all_ginis)})")
    
    # Secondary structure breakdown
    ss_counts = {"alpha": 0, "beta": 0, "coil": 0}
    for rn in primary["kappa_sq"]:
        if rn in dih:
            ss = rama_class(*dih[rn])
            if ss in ss_counts: ss_counts[ss] += 1
    primary["ss_counts"] = ss_counts
    if verbose:
        total = sum(ss_counts.values())
        if total > 0:
            print(f"      SS: α={ss_counts['alpha']} ({100*ss_counts['alpha']/total:.0f}%)  "
                  f"β={ss_counts['beta']} ({100*ss_counts['beta']/total:.0f}%)  "
                  f"coil={ss_counts['coil']} ({100*ss_counts['coil']/total:.0f}%)")
    
    return primary


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  BASELINE GRADIENT TEST: Four NMR-Characterized Proteins       ║")
    print("║  Section 18: Kinetic Connection — Geometry vs. Dynamics        ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    print(f"\n  Attractor: G ≈ {ATTRACTOR}")
    print(f"  Prediction: baseline < {ATTRACTOR} → stiffens (ΔG > 0)")
    print(f"              baseline > {ATTRACTOR} → loosens  (ΔG < 0)")
    
    # Download all PDB files
    print("\n── Downloading PDB files ──")
    paths = {}
    for pname, pdata in PROTEINS.items():
        for state_key in ["apo", "bound", "extra_bound"]:
            if state_key not in pdata: continue
            state = pdata[state_key]
            pid = state["pdb"]
            if pid not in paths:
                p = download_pdb(pid)
                if p: paths[pid] = p
    
    # Analyze each protein
    results = {}
    
    for pname, pdata in PROTEINS.items():
        print(f"\n{'═'*70}")
        print(f"  {pname.upper()}")
        print(f"  NMR S² at hotspot: {pdata['nmr_S2']}  ({pdata['nmr_ref']})")
        print(f"{'═'*70}")
        
        results[pname] = {}
        
        for state_key in ["apo", "bound", "extra_bound"]:
            if state_key not in pdata: continue
            state = pdata[state_key]
            pid = state["pdb"]
            if pid not in paths:
                print(f"    [SKIP] {pid} not available")
                continue
            
            r = analyze_one(
                paths[pid], state["chain"], 
                f"{pid}: {state['label']}", 
                state.get("is_nmr", False)
            )
            if r:
                results[pname][state_key] = r
                results[pname][state_key]["pdb"] = pid
    
    # ═══════════════════════════════════════════════════════════════
    # GRADIENT TEST
    # ═══════════════════════════════════════════════════════════════
    
    print("\n\n" + "╔" + "═"*68 + "╗")
    print("║  GRADIENT TEST RESULTS" + " "*45 + "║")
    print("╠" + "═"*68 + "╣")
    
    test_results = []
    
    for pname, pdata in PROTEINS.items():
        if "apo" not in results[pname]:
            continue
        
        g_apo = results[pname]["apo"]["gini"]
        zone = "LOW (stiffens)" if g_apo < ATTRACTOR else "HIGH (loosens)"
        predicted_dir = "+" if g_apo < ATTRACTOR else "-"
        
        print(f"║")
        print(f"║  {pname.upper()}")
        print(f"║    Baseline Gini = {g_apo:.4f}  ({zone})")
        
        for bound_key in ["bound", "extra_bound"]:
            if bound_key not in results[pname]:
                continue
            
            g_bound = results[pname][bound_key]["gini"]
            dg = g_bound - g_apo
            actual_dir = "+" if dg > 0 else "-"
            match = actual_dir == predicted_dir
            perturbation = pdata[bound_key].get("perturbation", "unknown")
            
            status = "✓ MATCH" if match else "✗ MISS"
            dir_label = "STIFFENS" if dg > 0 else "LOOSENS"
            
            print(f"║    → {pdata[bound_key]['pdb']}: {perturbation}")
            print(f"║      Gini = {g_bound:.4f}  ΔGini = {dg:+.4f} ({dir_label})")
            print(f"║      Predicted: {predicted_dir}  Actual: {actual_dir}  {status}")
            
            test_results.append({
                "protein": pname,
                "comparison": f"{pdata['apo']['pdb']} → {pdata[bound_key]['pdb']}",
                "perturbation": perturbation,
                "g_apo": g_apo,
                "g_bound": g_bound,
                "delta_gini": dg,
                "predicted": predicted_dir,
                "actual": actual_dir,
                "match": match,
            })
    
    # Summary
    n_tests = len(test_results)
    n_match = sum(1 for t in test_results if t["match"])
    
    print(f"║")
    print(f"╠" + "═"*68 + "╣")
    print(f"║  SUMMARY: {n_match}/{n_tests} predictions correct ({100*n_match/n_tests:.0f}%)" + " "*20 + "║" if n_tests > 0 else "")
    print(f"╠" + "═"*68 + "╣")
    
    # Detailed table
    print(f"║")
    print(f"║  {'Protein':12s} {'Baseline':>8s} {'Zone':>6s} {'ΔGini':>8s} {'Pred':>5s} {'Act':>5s} {'':>6s}")
    print(f"║  {'─'*12} {'─'*8} {'─'*6} {'─'*8} {'─'*5} {'─'*5} {'─'*6}")
    
    for t in test_results:
        zone = "LOW" if t["g_apo"] < ATTRACTOR else "HIGH"
        status = "✓" if t["match"] else "✗"
        print(f"║  {t['protein']:12s} {t['g_apo']:8.4f} {zone:>6s} {t['delta_gini']:+8.4f} "
              f"{t['predicted']:>5s} {t['actual']:>5s} {status:>6s}")
    
    print(f"║")
    print(f"╚" + "═"*68 + "╝")
    
    # Cross-check: does fold class explain the baseline?
    print("\n── Fold Class Context ──")
    for pname in results:
        if "apo" not in results[pname]: continue
        ss = results[pname]["apo"].get("ss_counts", {})
        total = sum(ss.values())
        if total > 0:
            alpha_pct = 100 * ss.get("alpha", 0) / total
            beta_pct = 100 * ss.get("beta", 0) / total
            g = results[pname]["apo"]["gini"]
            print(f"  {pname:12s}: Gini={g:.4f}  α={alpha_pct:.0f}%  β={beta_pct:.0f}%  "
                  f"{'α-rich' if alpha_pct > 40 else 'β-rich' if beta_pct > 40 else 'mixed'}")
    
    # Export
    export_path = DATA_DIR / "gradient_test_results.json"
    export = {}
    for pname in results:
        export[pname] = {}
        for state_key in results[pname]:
            r = results[pname][state_key]
            export[pname][state_key] = {
                "pdb": r.get("pdb", "?"),
                "gini": r["gini"],
                "n_res": r["n_res"],
                "n_kappa": r["n_kappa"],
                "kappa_sq": {str(k): v for k, v in r["kappa_sq"].items()},
            }
            if "ensemble_mean" in r:
                export[pname][state_key]["ensemble_mean"] = r["ensemble_mean"]
                export[pname][state_key]["ensemble_sd"] = r["ensemble_sd"]
    
    with open(export_path, "w") as f:
        json.dump(export, f, indent=2)
    print(f"\nResults exported to: {export_path}")
    
    # Also export the test summary
    summary_path = DATA_DIR / "gradient_test_summary.json"
    with open(summary_path, "w") as f:
        json.dump({
            "attractor": ATTRACTOR,
            "n_tests": n_tests,
            "n_match": n_match,
            "accuracy": n_match / n_tests if n_tests > 0 else 0,
            "tests": test_results,
        }, f, indent=2)
    print(f"Summary exported to: {summary_path}")


if __name__ == "__main__":
    main()

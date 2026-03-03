#!/usr/bin/env python3
"""
Quadrant Fix — Corrected Ramachandran Basin Definitions
========================================================
Re-classifies all residues from atlas_v2_full.json using
proper α-helix and β-sheet basin boundaries.

The v2 script used overly narrow basins:
  α: (-160,-30) x (-80,-10)   → missed most α-helix residues
  β: (-180,-60) x (60,180)    → missed β at negative ψ wrap

Corrected basins (standard Ramachandran):
  α-helix:     φ ∈ [-150, -30],  ψ ∈ [-70, -15]
  β-sheet:     φ ∈ [-180, -60],  ψ ∈ [90, 180] OR ψ ∈ [-180, -120]
  PPII/ext:    φ ∈ [-100, -40],  ψ ∈ [100, 180]
  Left-α:      φ ∈ [30, 100],    ψ ∈ [15, 70]

"Structured" = α OR β OR left-α (regular secondary structure)
"Flexible"   = everything else (coil, turns, PPII borders)

Expected result: ~50-65% structured, matching known protein statistics.

Usage:
    python quadrant_fix.py [--input atlas_v2_full.json] [--pdb-dir gradient_v3/pdbs]
"""

import numpy as np
import json, sys
from pathlib import Path
from collections import defaultdict
import argparse


# ═══════════════════════════════════════════════════════════════════
# BACKBONE PARSING (same as atlas)
# ═══════════════════════════════════════════════════════════════════

def parse_backbone(filepath):
    all_chains = defaultdict(dict)
    with open(filepath) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            an = line[12:16].strip()
            if an not in ("N", "CA", "C"):
                continue
            ch = line[21]
            try:
                rn = int(line[22:26].strip())
            except:
                continue
            if line[26].strip():
                continue
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            if rn not in all_chains[ch]:
                all_chains[ch][rn] = {}
            all_chains[ch][rn][an] = np.array([x, y, z])
    if not all_chains:
        return {}, ""
    best = max(all_chains.keys(), key=lambda c: len(all_chains[c]))
    return all_chains[best], best


def dihedral(p1, p2, p3, p4):
    b1, b2, b3 = p2-p1, p3-p2, p4-p3
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    nn1, nn2 = np.linalg.norm(n1), np.linalg.norm(n2)
    if nn1 < 1e-10 or nn2 < 1e-10:
        return np.nan
    n1, n2 = n1/nn1, n2/nn2
    m1 = np.cross(n1, b2/np.linalg.norm(b2))
    return np.arctan2(np.dot(m1, n2), np.dot(n1, n2))


def get_dihedrals(residues):
    sr = sorted(residues.keys())
    out = {}
    for i, rn in enumerate(sr):
        r = residues[rn]
        if not all(k in r for k in ("N", "CA", "C")):
            continue
        phi = np.nan
        if i > 0:
            prn = sr[i-1]; pr = residues.get(prn, {})
            if "C" in pr and (rn - prn) == 1:
                phi = dihedral(pr["C"], r["N"], r["CA"], r["C"])
        psi = np.nan
        if i < len(sr) - 1:
            nrn = sr[i+1]; nr = residues.get(nrn, {})
            if "N" in nr and (nrn - rn) == 1:
                psi = dihedral(r["N"], r["CA"], r["C"], nr["N"])
        out[rn] = (phi, psi)
    return out


def adiff(a, b):
    return np.arctan2(np.sin(a-b), np.cos(a-b))


def get_kappa_sq(dih):
    sr = sorted(dih.keys())
    ksq = {}
    for i in range(1, len(sr)-1):
        rp, rc, rn = sr[i-1], sr[i], sr[i+1]
        if (rc-rp) != 1 or (rn-rc) != 1:
            continue
        pp, sp = dih[rp]; pc, sc = dih[rc]; pn, sn = dih[rn]
        if any(np.isnan(x) for x in [pp, sp, pc, sc, pn, sn]):
            continue
        dp1, ds1 = adiff(pc, pp), adiff(sc, sp)
        dp2, ds2 = adiff(pn, pc), adiff(sn, sc)
        a1 = np.sqrt(dp1**2 + ds1**2)
        a2 = np.sqrt(dp2**2 + ds2**2)
        if a1 < 1e-10 or a2 < 1e-10:
            continue
        dTp = dp2/a2 - dp1/a1
        dTs = ds2/a2 - ds1/a1
        dm = 0.5*(a1+a2)
        ksq[rc] = (dTp**2 + dTs**2) / (dm**2)
    return ksq


# ═══════════════════════════════════════════════════════════════════
# CORRECTED RAMACHANDRAN CLASSIFICATION
# ═══════════════════════════════════════════════════════════════════

def classify_rama(phi_deg, psi_deg):
    """
    Classify a residue's (φ,ψ) into Ramachandran basin.
    Returns: 'alpha', 'beta', 'left_alpha', 'ppii', or 'other'
    
    Basin definitions follow Lovell et al. (2003) and 
    standard structural biology conventions.
    """
    p, s = phi_deg, psi_deg
    
    # Right-handed α-helix: core basin
    # Generous bounds to capture 3_10 and π-helix variants
    if -150 <= p <= -30 and -70 <= s <= -15:
        return "alpha"
    
    # Extended α region (captures 3_10 helix at more negative phi)
    if -180 <= p <= -150 and -70 <= s <= -15:
        return "alpha"
    
    # β-sheet: two regions due to ψ wrapping
    # Main β region (parallel and antiparallel)
    if -180 <= p <= -60 and 90 <= s <= 180:
        return "beta"
    # Wrapped β region (ψ near ±180)
    if -180 <= p <= -60 and -180 <= s <= -120:
        return "beta"
    
    # Left-handed α-helix (glycine, Asn/Asp)
    if 30 <= p <= 120 and 15 <= s <= 90:
        return "left_alpha"
    
    # PPII (polyproline II) — between α and β
    if -100 <= p <= -40 and 110 <= s <= 180:
        return "ppii"
    
    return "other"


def is_structured(phi_deg, psi_deg):
    """Is this residue in a regular secondary structure basin?"""
    basin = classify_rama(phi_deg, psi_deg)
    return basin in ("alpha", "beta", "left_alpha")


# ═══════════════════════════════════════════════════════════════════
# QUADRANT CLASSIFICATION
# ═══════════════════════════════════════════════════════════════════

def classify_protein(ksq_dict, dih_dict):
    """
    Classify residues into quadrants using corrected basins.
    
    Beam:   low κ², structured → rigid backbone in regular SS
    Loop:   low κ², flexible   → smooth curves in coil/turn
    Pin:    high κ², structured → stress concentrator in SS
    Spring: high κ², flexible   → mechanical spring in coil
    """
    if not ksq_dict:
        return {}, {}, {}
    
    kvals = np.array(list(ksq_dict.values()))
    median_k = np.median(kvals)
    
    counts = {"Beam": 0, "Loop": 0, "Pin": 0, "Spring": 0}
    basin_counts = {"alpha": 0, "beta": 0, "left_alpha": 0, "ppii": 0, "other": 0}
    
    for rn, ksq in ksq_dict.items():
        if rn not in dih_dict:
            continue
        phi, psi = dih_dict[rn]
        if np.isnan(phi) or np.isnan(psi):
            continue
        
        # Negate both phi and psi: our dihedral convention produces
        # -standard_angle (cross product order). κ² unaffected (sign 
        # cancels in angular differences). Verified: without negation,
        # 39% falls in left-alpha instead of alpha.
        phi_d = -np.degrees(phi)
        psi_d = -np.degrees(psi)
        
        basin = classify_rama(phi_d, psi_d)
        basin_counts[basin] += 1
        
        structured = basin in ("alpha", "beta", "left_alpha")
        high_k = ksq > median_k
        
        if high_k and structured:
            counts["Pin"] += 1
        elif high_k and not structured:
            counts["Spring"] += 1
        elif not high_k and structured:
            counts["Beam"] += 1
        else:
            counts["Loop"] += 1
    
    total = sum(counts.values())
    fracs = {k: v/total if total > 0 else 0.0 for k, v in counts.items()}
    
    return counts, fracs, basin_counts


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-dir", default="gradient_v3/pdbs")
    parser.add_argument("--output-dir", default="crystal_atlas_v2")
    args = parser.parse_args()
    
    pdb_dir = Path(args.pdb_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(exist_ok=True)
    
    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    if not pdb_files:
        for alt in ["gradient_large_scale/pdbs", "gradient_temperature/pdbs", "pdbs"]:
            pdb_files = sorted(Path(alt).glob("*.pdb"))
            if pdb_files:
                pdb_dir = Path(alt)
                break
    
    if not pdb_files:
        print(f"No .pdb files in {pdb_dir}")
        sys.exit(1)
    
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  QUADRANT FIX — Corrected Basin Definitions                ║")
    print(f"║  {len(pdb_files):>5d} PDB files                                        ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    
    results = []
    failed = 0
    
    # Global accumulators
    g_counts = {"Beam": 0, "Loop": 0, "Pin": 0, "Spring": 0}
    g_basins = {"alpha": 0, "beta": 0, "left_alpha": 0, "ppii": 0, "other": 0}
    
    for idx, fp in enumerate(pdb_files):
        if idx % 500 == 0 and idx > 0:
            print(f"  {idx}/{len(pdb_files)} ({len(results)} valid)...")
        
        pdb_id = fp.stem
        residues, chain = parse_backbone(str(fp))
        if len(residues) < 15:
            failed += 1; continue
        
        dih = get_dihedrals(residues)
        ksq = get_kappa_sq(dih)
        if len(ksq) < 10:
            failed += 1; continue
        
        counts, fracs, basins = classify_protein(ksq, dih)
        
        # Accumulate globals
        for k in g_counts: g_counts[k] += counts[k]
        for k in g_basins: g_basins[k] += basins[k]
        
        # Gini (recompute — fast)
        kvals = list(ksq.values())
        v = np.sort(np.array(kvals))
        n = len(v); s = np.sum(v)
        gini_val = (2*np.sum(np.arange(1,n+1)*v) / (n*s)) - (n+1)/n if s > 0 else 0.0
        
        results.append({
            "pdb": pdb_id,
            "chain": chain,
            "n_res": len(residues),
            "n_kappa": len(ksq),
            "gini": gini_val,
            "quad_counts": counts,
            "quad_fracs": fracs,
            "basin_counts": basins,
            "pct_structured": (basins["alpha"]+basins["beta"]+basins["left_alpha"]) / 
                              max(sum(basins.values()), 1) * 100,
        })
    
    print(f"  {len(pdb_files)} files → {len(results)} valid ({failed} failed)\n")
    
    # ═══════════════════════════════════════════════════════════════
    # RESULTS
    # ═══════════════════════════════════════════════════════════════
    
    total_q = sum(g_counts.values())
    total_b = sum(g_basins.values())
    
    print("="*60)
    print(f"CORRECTED QUADRANT CLASSIFICATION (N={len(results)})")
    print("="*60)
    
    # Basin validation
    pct_struct = (g_basins["alpha"]+g_basins["beta"]+g_basins["left_alpha"])/total_b*100
    print(f"\n── Ramachandran Basin Populations ──")
    print(f"  Total residues: {total_b:,}")
    print(f"  α-helix:    {g_basins['alpha']:>8,} ({100*g_basins['alpha']/total_b:.1f}%)")
    print(f"  β-sheet:    {g_basins['beta']:>8,} ({100*g_basins['beta']/total_b:.1f}%)")
    print(f"  Left-α:     {g_basins['left_alpha']:>8,} ({100*g_basins['left_alpha']/total_b:.1f}%)")
    print(f"  PPII:       {g_basins['ppii']:>8,} ({100*g_basins['ppii']/total_b:.1f}%)")
    print(f"  Other:      {g_basins['other']:>8,} ({100*g_basins['other']/total_b:.1f}%)")
    print(f"  Structured: {pct_struct:.1f}%  (expected: 50-65%)")
    
    if pct_struct < 40:
        print(f"  ⚠ Still low — may need wider basins")
    elif pct_struct > 70:
        print(f"  ⚠ High — basins may be too generous")
    else:
        print(f"  ✓ Within expected range")
    
    # Quadrant populations
    print(f"\n── Quadrant Populations (global) ──")
    print(f"  Total classified: {total_q:,}")
    print(f"  Beam:   {g_counts['Beam']:>8,}  ({100*g_counts['Beam']/total_q:>5.1f}%)  low κ², structured")
    print(f"  Loop:   {g_counts['Loop']:>8,}  ({100*g_counts['Loop']/total_q:>5.1f}%)  low κ², flexible")
    print(f"  Pin:    {g_counts['Pin']:>8,}  ({100*g_counts['Pin']/total_q:>5.1f}%)  high κ², structured")
    print(f"  Spring: {g_counts['Spring']:>8,}  ({100*g_counts['Spring']/total_q:>5.1f}%)  high κ², flexible")
    
    # Per-protein stats
    beam_f = np.array([r["quad_fracs"]["Beam"] for r in results])
    loop_f = np.array([r["quad_fracs"]["Loop"] for r in results])
    pin_f = np.array([r["quad_fracs"]["Pin"] for r in results])
    spring_f = np.array([r["quad_fracs"]["Spring"] for r in results])
    ginis = np.array([r["gini"] for r in results])
    n_res = np.array([r["n_res"] for r in results])
    pct_s = np.array([r["pct_structured"] for r in results])
    
    print(f"\n── Per-Protein Mean Fractions ──")
    print(f"  Beam:   {np.mean(beam_f):.3f} ± {np.std(beam_f):.3f}")
    print(f"  Loop:   {np.mean(loop_f):.3f} ± {np.std(loop_f):.3f}")
    print(f"  Pin:    {np.mean(pin_f):.3f} ± {np.std(pin_f):.3f}")
    print(f"  Spring: {np.mean(spring_f):.3f} ± {np.std(spring_f):.3f}")
    
    # Correlations
    print(f"\n── Gini vs Quadrant Correlations ──")
    print(f"  r(Gini, %Beam):    {np.corrcoef(ginis, beam_f)[0,1]:+.4f}")
    print(f"  r(Gini, %Loop):    {np.corrcoef(ginis, loop_f)[0,1]:+.4f}")
    print(f"  r(Gini, %Pin):     {np.corrcoef(ginis, pin_f)[0,1]:+.4f}")
    print(f"  r(Gini, %Spring):  {np.corrcoef(ginis, spring_f)[0,1]:+.4f}")
    print(f"  r(Gini, %struct):  {np.corrcoef(ginis, pct_s)[0,1]:+.4f}")
    
    # Size bins
    print(f"\n── By Size ──")
    print(f"  {'Range':>12s} {'N':>5s} {'Gini':>7s} {'Beam':>6s} {'Loop':>6s} {'Pin':>6s} {'Spr':>6s} {'%SS':>5s}")
    print("  " + "─"*55)
    for lo, hi in [(50,100),(100,150),(150,200),(200,300),(300,500),(500,1000)]:
        m = (n_res >= lo) & (n_res < hi); n = np.sum(m)
        if n < 5: continue
        print(f"  [{lo:>3d},{hi:>3d}) {n:>5d} {np.mean(ginis[m]):>7.4f} "
              f"{np.mean(beam_f[m]):>6.3f} {np.mean(loop_f[m]):>6.3f} "
              f"{np.mean(pin_f[m]):>6.3f} {np.mean(spring_f[m]):>6.3f} "
              f"{np.mean(pct_s[m]):>5.1f}")
    
    # Export
    csv_path = out_dir / "atlas_v2_corrected.csv"
    with open(csv_path, "w") as f:
        f.write("pdb,chain,n_res,gini,beam_frac,loop_frac,pin_frac,spring_frac,"
                "pct_structured,pct_alpha,pct_beta\n")
        for r in sorted(results, key=lambda x: x["gini"]):
            tb = sum(r["basin_counts"].values())
            pa = 100*r["basin_counts"]["alpha"]/tb if tb > 0 else 0
            pb = 100*r["basin_counts"]["beta"]/tb if tb > 0 else 0
            f.write(f"{r['pdb']},{r['chain']},{r['n_res']},{r['gini']:.6f},"
                    f"{r['quad_fracs']['Beam']:.4f},{r['quad_fracs']['Loop']:.4f},"
                    f"{r['quad_fracs']['Pin']:.4f},{r['quad_fracs']['Spring']:.4f},"
                    f"{r['pct_structured']:.1f},{pa:.1f},{pb:.1f}\n")
    print(f"\n  CSV: {csv_path}")
    
    json_path = out_dir / "atlas_v2_corrected.json"
    summary = {
        "n_proteins": len(results),
        "basins_global": g_basins,
        "basins_global_pct": {k: round(100*v/total_b, 1) for k, v in g_basins.items()},
        "pct_structured_global": round(pct_struct, 1),
        "quadrants_global": g_counts,
        "quadrants_global_pct": {k: round(100*v/total_q, 1) for k, v in g_counts.items()},
        "quadrants_per_protein": {
            "Beam": {"mean": round(float(np.mean(beam_f)),4), "std": round(float(np.std(beam_f)),4)},
            "Loop": {"mean": round(float(np.mean(loop_f)),4), "std": round(float(np.std(loop_f)),4)},
            "Pin": {"mean": round(float(np.mean(pin_f)),4), "std": round(float(np.std(pin_f)),4)},
            "Spring": {"mean": round(float(np.mean(spring_f)),4), "std": round(float(np.std(spring_f)),4)},
        },
        "correlations": {
            "gini_beam": round(float(np.corrcoef(ginis, beam_f)[0,1]),4),
            "gini_loop": round(float(np.corrcoef(ginis, loop_f)[0,1]),4),
            "gini_pin": round(float(np.corrcoef(ginis, pin_f)[0,1]),4),
            "gini_spring": round(float(np.corrcoef(ginis, spring_f)[0,1]),4),
            "gini_pct_struct": round(float(np.corrcoef(ginis, pct_s)[0,1]),4),
        },
        "results": results,
    }
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  JSON: {json_path}")


if __name__ == "__main__":
    main()

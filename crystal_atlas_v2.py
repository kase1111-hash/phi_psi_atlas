#!/usr/bin/env python3
"""
Crystal Atlas v2 — Mechanical Taxonomy at Scale
=================================================
Computes the full mechanical taxonomy for every PDB file on disk:
  - Per-residue κ² (discrete backbone curvature)
  - Gini coefficient (curvature concentration)
  - Quadrant classification (Beam/Loop/Spring/Pin)
  - Population statistics across entire dataset

Replaces the original 997-structure atlas with 3000+ structures.

No network access needed. Runs on cached PDB files.

Usage:
    python crystal_atlas_v2.py --pdb-dir gradient_v3/pdbs

Author: Kase / True North Construction LLC
License: CC0 1.0 Universal (unpatentable prior art)
"""

import numpy as np
import os, sys, json
from pathlib import Path
from collections import defaultdict
import argparse

# ═══════════════════════════════════════════════════════════════════
# BACKBONE GEOMETRY (identical methodology to all prior scripts)
# ═══════════════════════════════════════════════════════════════════

def parse_backbone(filepath):
    """Parse backbone N, CA, C atoms for all chains. Return longest protein chain."""
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


def gini(values):
    v = np.sort(np.array(values, dtype=float))
    v = v[~np.isnan(v)]
    if len(v) == 0: return np.nan
    n = len(v); s = np.sum(v)
    if s < 1e-15: return 0.0
    return (2*np.sum(np.arange(1,n+1)*v) / (n*s)) - (n+1)/n


# ═══════════════════════════════════════════════════════════════════
# QUADRANT CLASSIFICATION
# ═══════════════════════════════════════════════════════════════════

def classify_residues(kappa_sq_dict, dih_dict):
    """
    Classify each residue into mechanical quadrant based on:
      - κ² value (high/low relative to median)
      - Ramachandran position (structured/flexible)
    
    Quadrants:
      Beam:   low κ², structured (α-helix, β-sheet)  → rigid rod
      Loop:   low κ², flexible (coil, turn)           → flexible hinge
      Pin:    high κ², structured                      → concentrated stress
      Spring: high κ², flexible                        → mechanical spring
    
    "Structured" = residue (φ,ψ) in α or β basin
    "Flexible" = everything else (coil, turn, PPII, left-handed)
    """
    if not kappa_sq_dict:
        return {}
    
    kvals = np.array(list(kappa_sq_dict.values()))
    median_k = np.median(kvals)
    
    classifications = {}
    for rn, ksq in kappa_sq_dict.items():
        if rn not in dih_dict:
            continue
        phi, psi = dih_dict[rn]
        if np.isnan(phi) or np.isnan(psi):
            continue
        
        # Convert to degrees for basin classification
        phi_d = np.degrees(phi)
        psi_d = np.degrees(psi)
        
        # Structured = in α-helix or β-sheet Ramachandran basin
        in_alpha = (-160 < phi_d < -30) and (-80 < psi_d < -10)
        in_beta = (-180 < phi_d < -60) and (60 < psi_d < 180)
        # Also include left-extended β
        in_beta2 = (-180 < phi_d < -60) and (-180 < psi_d < -120)
        structured = in_alpha or in_beta or in_beta2
        
        high_k = ksq > median_k
        
        if high_k and structured:
            q = "Pin"
        elif high_k and not structured:
            q = "Spring"
        elif not high_k and structured:
            q = "Beam"
        else:
            q = "Loop"
        
        classifications[rn] = {
            "quadrant": q,
            "kappa_sq": ksq,
            "phi": phi_d,
            "psi": psi_d,
            "structured": structured,
            "high_kappa": high_k,
        }
    
    return classifications


def quadrant_stats(classifications):
    """Compute quadrant population fractions."""
    counts = {"Beam": 0, "Loop": 0, "Pin": 0, "Spring": 0}
    for info in classifications.values():
        counts[info["quadrant"]] += 1
    total = sum(counts.values())
    if total == 0:
        return counts, {k: 0.0 for k in counts}
    fracs = {k: v/total for k, v in counts.items()}
    return counts, fracs


# ═══════════════════════════════════════════════════════════════════
# FULL PROTEIN ANALYSIS
# ═══════════════════════════════════════════════════════════════════

def analyze_protein(filepath):
    """Complete analysis: κ², Gini, quadrant classification."""
    residues, chain = parse_backbone(filepath)
    if len(residues) < 15:
        return None
    
    dih = get_dihedrals(residues)
    ksq = get_kappa_sq(dih)
    if len(ksq) < 10:
        return None
    
    kvals = list(ksq.values())
    g = gini(kvals)
    
    # Classify
    classifications = classify_residues(ksq, dih)
    counts, fracs = quadrant_stats(classifications)
    
    # κ² distribution stats
    karr = np.array(kvals)
    
    return {
        "chain": chain,
        "n_res": len(residues),
        "n_kappa": len(ksq),
        "gini": g,
        "kappa_mean": float(np.mean(karr)),
        "kappa_median": float(np.median(karr)),
        "kappa_std": float(np.std(karr)),
        "kappa_max": float(np.max(karr)),
        "kappa_p90": float(np.percentile(karr, 90)),
        "kappa_p95": float(np.percentile(karr, 95)),
        "kappa_p99": float(np.percentile(karr, 99)),
        "quad_counts": counts,
        "quad_fracs": fracs,
        "n_classified": sum(counts.values()),
    }


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Crystal Atlas v2")
    parser.add_argument("--pdb-dir", type=str, default="gradient_v3/pdbs")
    parser.add_argument("--output-dir", type=str, default="crystal_atlas_v2")
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
        print(f"No .pdb files in {pdb_dir}. Use --pdb-dir <path>")
        sys.exit(1)
    
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  CRYSTAL ATLAS v2 — Mechanical Taxonomy at Scale           ║")
    print(f"║  {len(pdb_files):>5d} PDB files                                        ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    
    # Process all files
    results = []
    failed = 0
    all_kappa_values = []  # pooled across all proteins
    
    for idx, fp in enumerate(pdb_files):
        if idx % 500 == 0 and idx > 0:
            print(f"  Processed {idx}/{len(pdb_files)} ({len(results)} valid)...")
        
        pdb_id = fp.stem
        r = analyze_protein(str(fp))
        if r is None:
            failed += 1
            continue
        
        r["pdb"] = pdb_id
        results.append(r)
    
    print(f"  Processed {len(pdb_files)} files → {len(results)} valid ({failed} failed)\n")
    
    if len(results) < 20:
        print("Too few results."); sys.exit(1)
    
    # ═══════════════════════════════════════════════════════════════
    # ANALYSIS
    # ═══════════════════════════════════════════════════════════════
    
    ginis = np.array([r["gini"] for r in results])
    n_res = np.array([r["n_res"] for r in results])
    
    print("="*65)
    print(f"CRYSTAL ATLAS v2: {len(results)} proteins")
    print(f"  (Original atlas: 997 proteins)")
    print("="*65)
    
    # ── 1. Gini Distribution ──
    print(f"\n{'─'*40}")
    print("1. GINI DISTRIBUTION")
    print(f"{'─'*40}")
    print(f"  N:       {len(ginis)}")
    print(f"  Mean:    {np.mean(ginis):.4f}")
    print(f"  Median:  {np.median(ginis):.4f}")
    print(f"  Std:     {np.std(ginis):.4f}")
    print(f"  IQR:     [{np.percentile(ginis,25):.4f}, {np.percentile(ginis,75):.4f}]")
    print(f"  95% CI:  [{np.percentile(ginis,2.5):.4f}, {np.percentile(ginis,97.5):.4f}]")
    print(f"  Range:   [{np.min(ginis):.4f}, {np.max(ginis):.4f}]")
    
    # Skewness/kurtosis
    m = np.mean(ginis); s = np.std(ginis)
    sk = float(np.mean(((ginis - m) / s) ** 3))
    ku = float(np.mean(((ginis - m) / s) ** 4) - 3)
    print(f"  Skew:    {sk:.4f}")
    print(f"  Kurt:    {ku:.4f}")
    
    # Histogram
    print(f"\n  Histogram:")
    bins = np.arange(0.55, 1.01, 0.025)
    counts, _ = np.histogram(ginis, bins=bins)
    mx = max(counts)
    for i in range(len(counts)):
        bar = "█" * int(40 * counts[i] / mx) if mx > 0 else ""
        print(f"    [{bins[i]:.3f},{bins[i+1]:.3f}) {counts[i]:>4d} {bar}")
    
    # Mode
    fine_bins = np.arange(0.60, 1.0, 0.005)
    fc, _ = np.histogram(ginis, bins=fine_bins)
    kernel = np.ones(5)/5
    sm = np.convolve(fc, kernel, mode='same')
    peak_idx = np.argmax(sm)
    mode_g = (fine_bins[peak_idx] + fine_bins[peak_idx+1]) / 2
    print(f"\n  Mode:    {mode_g:.3f}")
    
    # ── 2. Quadrant Populations ──
    print(f"\n{'─'*40}")
    print("2. QUADRANT POPULATIONS (global)")
    print(f"{'─'*40}")
    
    total_beam = sum(r["quad_counts"]["Beam"] for r in results)
    total_loop = sum(r["quad_counts"]["Loop"] for r in results)
    total_pin = sum(r["quad_counts"]["Pin"] for r in results)
    total_spring = sum(r["quad_counts"]["Spring"] for r in results)
    total_all = total_beam + total_loop + total_pin + total_spring
    
    print(f"  Total classified residues: {total_all:,}")
    print(f"  Beam:   {total_beam:>8,}  ({100*total_beam/total_all:>5.1f}%)  low κ², structured")
    print(f"  Loop:   {total_loop:>8,}  ({100*total_loop/total_all:>5.1f}%)  low κ², flexible")
    print(f"  Pin:    {total_pin:>8,}  ({100*total_pin/total_all:>5.1f}%)  high κ², structured")
    print(f"  Spring: {total_spring:>8,}  ({100*total_spring/total_all:>5.1f}%)  high κ², flexible")
    
    # Per-protein quadrant stats
    beam_fracs = np.array([r["quad_fracs"]["Beam"] for r in results])
    loop_fracs = np.array([r["quad_fracs"]["Loop"] for r in results])
    pin_fracs = np.array([r["quad_fracs"]["Pin"] for r in results])
    spring_fracs = np.array([r["quad_fracs"]["Spring"] for r in results])
    
    print(f"\n  Per-protein mean fractions:")
    print(f"    Beam:   {np.mean(beam_fracs):.3f} ± {np.std(beam_fracs):.3f}")
    print(f"    Loop:   {np.mean(loop_fracs):.3f} ± {np.std(loop_fracs):.3f}")
    print(f"    Pin:    {np.mean(pin_fracs):.3f} ± {np.std(pin_fracs):.3f}")
    print(f"    Spring: {np.mean(spring_fracs):.3f} ± {np.std(spring_fracs):.3f}")
    
    # ── 3. κ² Distribution ──
    print(f"\n{'─'*40}")
    print("3. κ² DISTRIBUTION (pooled)")
    print(f"{'─'*40}")
    
    all_means = np.array([r["kappa_mean"] for r in results])
    all_medians = np.array([r["kappa_median"] for r in results])
    all_p95 = np.array([r["kappa_p95"] for r in results])
    
    print(f"  Per-protein mean κ²:   {np.mean(all_means):.4f} ± {np.std(all_means):.4f}")
    print(f"  Per-protein median κ²: {np.mean(all_medians):.4f} ± {np.std(all_medians):.4f}")
    print(f"  Per-protein P95 κ²:    {np.mean(all_p95):.4f} ± {np.std(all_p95):.4f}")
    
    # ── 4. Size Dependence ──
    print(f"\n{'─'*40}")
    print("4. SIZE DEPENDENCE")
    print(f"{'─'*40}")
    
    r_size_gini = np.corrcoef(n_res, ginis)[0,1]
    r_size_beam = np.corrcoef(n_res, beam_fracs)[0,1]
    r_size_pin = np.corrcoef(n_res, pin_fracs)[0,1]
    print(f"  r(n_res, Gini):     {r_size_gini:.4f}")
    print(f"  r(n_res, %Beam):    {r_size_beam:.4f}")
    print(f"  r(n_res, %Pin):     {r_size_pin:.4f}")
    
    print(f"\n  By size:")
    print(f"  {'Range':>12s} {'N':>5s} {'Gini':>7s} {'Beam':>6s} {'Loop':>6s} {'Pin':>6s} {'Spr':>6s}")
    print("  " + "─"*50)
    for lo, hi in [(50,100),(100,150),(150,200),(200,300),(300,500)]:
        m = (n_res >= lo) & (n_res < hi); n = np.sum(m)
        if n < 5: continue
        print(f"  [{lo:>3d},{hi:>3d}) {n:>5d} {np.mean(ginis[m]):>7.4f} "
              f"{np.mean(beam_fracs[m]):>6.3f} {np.mean(loop_fracs[m]):>6.3f} "
              f"{np.mean(pin_fracs[m]):>6.3f} {np.mean(spring_fracs[m]):>6.3f}")
    
    # ── 5. Gini vs Quadrant Correlations ──
    print(f"\n{'─'*40}")
    print("5. GINI vs QUADRANT CORRELATIONS")
    print(f"{'─'*40}")
    
    r_g_beam = np.corrcoef(ginis, beam_fracs)[0,1]
    r_g_loop = np.corrcoef(ginis, loop_fracs)[0,1]
    r_g_pin = np.corrcoef(ginis, pin_fracs)[0,1]
    r_g_spring = np.corrcoef(ginis, spring_fracs)[0,1]
    print(f"  r(Gini, %Beam):    {r_g_beam:+.4f}")
    print(f"  r(Gini, %Loop):    {r_g_loop:+.4f}")
    print(f"  r(Gini, %Pin):     {r_g_pin:+.4f}")
    print(f"  r(Gini, %Spring):  {r_g_spring:+.4f}")
    print(f"\n  Interpretation:")
    if r_g_pin > 0.3:
        print(f"    High Gini ↔ more Pin residues (concentrated stress)")
    if r_g_beam > 0.3:
        print(f"    High Gini ↔ more Beam residues")
    if r_g_loop < -0.3:
        print(f"    High Gini ↔ fewer Loop residues")
    
    # ── 6. Comparison with Original Atlas ──
    print(f"\n{'─'*40}")
    print("6. COMPARISON WITH ORIGINAL ATLAS (997)")
    print(f"{'─'*40}")
    print(f"  {'Metric':>25s} {'Original':>10s} {'v2 (N={len(results)})':>15s}")
    print("  " + "─"*52)
    print(f"  {'Mean Gini':>25s} {'~0.835':>10s} {np.mean(ginis):>15.4f}")
    print(f"  {'Median Gini':>25s} {'~0.840':>10s} {np.median(ginis):>15.4f}")
    print(f"  {'Std Gini':>25s} {'~0.05':>10s} {np.std(ginis):>15.4f}")
    print(f"  {'Mode Gini':>25s} {'~0.86':>10s} {mode_g:>15.3f}")
    print(f"  {'%Beam (global)':>25s} {'~40%':>10s} {100*total_beam/total_all:>14.1f}%")
    print(f"  {'%Loop (global)':>25s} {'~12%':>10s} {100*total_loop/total_all:>14.1f}%")
    print(f"  {'%Pin (global)':>25s} {'~38%':>10s} {100*total_pin/total_all:>14.1f}%")
    print(f"  {'%Spring (global)':>25s} {'~10%':>10s} {100*total_spring/total_all:>14.1f}%")
    
    # ── Export ──
    print(f"\n{'─'*40}")
    print("EXPORT")
    print(f"{'─'*40}")
    
    # CSV
    csv_path = out_dir / "atlas_v2.csv"
    with open(csv_path, "w") as f:
        f.write("pdb,chain,n_res,gini,kappa_mean,kappa_median,kappa_p95,"
                "beam_frac,loop_frac,pin_frac,spring_frac,"
                "beam_n,loop_n,pin_n,spring_n\n")
        for r in sorted(results, key=lambda x: x["gini"]):
            f.write(f"{r['pdb']},{r['chain']},{r['n_res']},{r['gini']:.6f},"
                    f"{r['kappa_mean']:.6f},{r['kappa_median']:.6f},{r['kappa_p95']:.6f},"
                    f"{r['quad_fracs']['Beam']:.4f},{r['quad_fracs']['Loop']:.4f},"
                    f"{r['quad_fracs']['Pin']:.4f},{r['quad_fracs']['Spring']:.4f},"
                    f"{r['quad_counts']['Beam']},{r['quad_counts']['Loop']},"
                    f"{r['quad_counts']['Pin']},{r['quad_counts']['Spring']}\n")
    print(f"  CSV: {csv_path}")
    
    # JSON summary
    json_path = out_dir / "atlas_v2_summary.json"
    summary = {
        "n_proteins": len(results),
        "n_failed": failed,
        "gini": {
            "mean": float(np.mean(ginis)),
            "median": float(np.median(ginis)),
            "std": float(np.std(ginis)),
            "mode": float(mode_g),
            "skew": sk,
            "kurtosis": ku,
            "iqr": [float(np.percentile(ginis,25)), float(np.percentile(ginis,75))],
            "ci95": [float(np.percentile(ginis,2.5)), float(np.percentile(ginis,97.5))],
        },
        "quadrants_global": {
            "Beam": total_beam, "Loop": total_loop,
            "Pin": total_pin, "Spring": total_spring,
            "total": total_all,
        },
        "quadrants_global_pct": {
            "Beam": round(100*total_beam/total_all, 1),
            "Loop": round(100*total_loop/total_all, 1),
            "Pin": round(100*total_pin/total_all, 1),
            "Spring": round(100*total_spring/total_all, 1),
        },
        "quadrants_per_protein_mean": {
            "Beam": round(float(np.mean(beam_fracs)), 4),
            "Loop": round(float(np.mean(loop_fracs)), 4),
            "Pin": round(float(np.mean(pin_fracs)), 4),
            "Spring": round(float(np.mean(spring_fracs)), 4),
        },
        "correlations": {
            "size_vs_gini": round(r_size_gini, 4),
            "gini_vs_beam": round(r_g_beam, 4),
            "gini_vs_loop": round(r_g_loop, 4),
            "gini_vs_pin": round(r_g_pin, 4),
            "gini_vs_spring": round(r_g_spring, 4),
        },
        "kappa_stats": {
            "mean_of_means": round(float(np.mean(all_means)), 4),
            "mean_of_medians": round(float(np.mean(all_medians)), 4),
            "mean_of_p95": round(float(np.mean(all_p95)), 4),
        },
        "histogram_bins": bins.tolist(),
        "histogram_counts": counts.tolist(),
    }
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  JSON: {json_path}")
    
    # Full per-protein data
    full_path = out_dir / "atlas_v2_full.json"
    with open(full_path, "w") as f:
        json.dump(results, f)
    print(f"  Full: {full_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Gini Census — Population Distribution of Curvature Concentration
================================================================
Computes Gini(κ²) for every PDB file in the download cache.
No network access needed — works entirely on local files.

Outputs:
  - gini_census.csv:  PDB, chain, Gini, n_residues, n_kappa
  - gini_census.json: Full stats + histogram data
  - Console: Distribution analysis, normality tests, fold-size correlations

Usage:
    python gini_census.py [--pdb-dir gradient_v3/pdbs] [--chain A]

Author: Kase / True North Construction LLC
License: CC0 1.0 Universal
"""

import numpy as np
import os, sys, json, glob
from pathlib import Path
from collections import defaultdict

# ═══════════════════════════════════════════════════════════════════
# CURVATURE (identical to atlas methodology)
# ═══════════════════════════════════════════════════════════════════

def parse_backbone(filepath, chain_id=None):
    """Parse backbone N, CA, C atoms. If chain_id is None, use longest chain."""
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
            if line[26].strip():  # skip alt locations
                continue
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            if rn not in all_chains[ch]:
                all_chains[ch][rn] = {}
            all_chains[ch][rn][an] = np.array([x, y, z])
    
    if chain_id and chain_id in all_chains:
        return all_chains[chain_id], chain_id
    
    # Pick longest chain
    if not all_chains:
        return {}, ""
    best_ch = max(all_chains.keys(), key=lambda c: len(all_chains[c]))
    return all_chains[best_ch], best_ch


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

def compute_gini(filepath, chain_id=None):
    """Full pipeline: PDB → Gini. Returns dict or None."""
    residues, used_chain = parse_backbone(filepath, chain_id)
    if len(residues) < 10:
        return None
    dih = get_dihedrals(residues)
    ksq = get_kappa_sq(dih)
    if len(ksq) < 5:
        return None
    
    kvals = list(ksq.values())
    g = gini(kvals)
    
    return {
        "gini": g,
        "chain": used_chain,
        "n_res": len(residues),
        "n_kappa": len(ksq),
        "mean_kappa": float(np.mean(kvals)),
        "median_kappa": float(np.median(kvals)),
        "max_kappa": float(np.max(kvals)),
    }


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Gini census across PDB files")
    parser.add_argument("--pdb-dir", type=str, default="gradient_v3/pdbs",
                        help="Directory containing .pdb files")
    parser.add_argument("--output-dir", type=str, default="gini_census",
                        help="Output directory")
    args = parser.parse_args()
    
    pdb_dir = Path(args.pdb_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(exist_ok=True)
    
    # Find all PDB files
    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    if not pdb_files:
        # Try other common locations
        for alt in ["gradient_large_scale/pdbs", "gradient_temperature/pdbs", "pdbs", "."]:
            pdb_files = sorted(Path(alt).glob("*.pdb"))
            if pdb_files:
                print(f"  Found PDB files in {alt}/")
                break
    
    if not pdb_files:
        print(f"No .pdb files found in {pdb_dir}. Try --pdb-dir <path>")
        sys.exit(1)
    
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  GINI CENSUS — Population Distribution                     ║")
    print(f"║  {len(pdb_files)} PDB files                                          ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    
    # Compute Gini for each
    results = []
    failed = 0
    
    for idx, fp in enumerate(pdb_files):
        if idx % 200 == 0 and idx > 0:
            print(f"  Processed {idx}/{len(pdb_files)} ({len(results)} valid)...")
        
        pdb_id = fp.stem
        result = compute_gini(str(fp))
        if result is None:
            failed += 1
            continue
        
        results.append({
            "pdb": pdb_id,
            "chain": result["chain"],
            "gini": result["gini"],
            "n_res": result["n_res"],
            "n_kappa": result["n_kappa"],
            "mean_kappa": result["mean_kappa"],
            "median_kappa": result["median_kappa"],
            "max_kappa": result["max_kappa"],
        })
    
    print(f"  Processed {len(pdb_files)}/{len(pdb_files)} ({len(results)} valid, {failed} failed)")
    
    if len(results) < 10:
        print("Too few valid results.")
        sys.exit(1)
    
    # ═══════════════════════════════════════════════════════════════
    # ANALYSIS
    # ═══════════════════════════════════════════════════════════════
    
    ginis = np.array([r["gini"] for r in results])
    n_res = np.array([r["n_res"] for r in results])
    mean_k = np.array([r["mean_kappa"] for r in results])
    
    print(f"\n{'='*60}")
    print(f"GINI CENSUS: {len(results)} proteins")
    print(f"{'='*60}")
    
    print(f"\n── Distribution ──")
    print(f"  Mean:   {np.mean(ginis):.4f}")
    print(f"  Median: {np.median(ginis):.4f}")
    print(f"  Std:    {np.std(ginis):.4f}")
    print(f"  Min:    {np.min(ginis):.4f}")
    print(f"  Max:    {np.max(ginis):.4f}")
    print(f"  IQR:    [{np.percentile(ginis, 25):.4f}, {np.percentile(ginis, 75):.4f}]")
    print(f"  95%:    [{np.percentile(ginis, 2.5):.4f}, {np.percentile(ginis, 97.5):.4f}]")
    
    # Skewness and kurtosis
    from scipy.stats import skew, kurtosis, shapiro, normaltest
    try:
        sk = skew(ginis)
        ku = kurtosis(ginis)
        print(f"  Skewness: {sk:.4f}")
        print(f"  Kurtosis: {ku:.4f} (0 = normal)")
    except ImportError:
        # Manual skewness/kurtosis
        m = np.mean(ginis); s = np.std(ginis)
        sk = np.mean(((ginis - m) / s) ** 3)
        ku = np.mean(((ginis - m) / s) ** 4) - 3
        print(f"  Skewness: {sk:.4f}")
        print(f"  Excess kurtosis: {ku:.4f}")
    
    # Normality test
    try:
        if len(ginis) <= 5000:
            stat, p = shapiro(ginis)
            print(f"  Shapiro-Wilk: W={stat:.4f}, p={p:.4e}")
        stat2, p2 = normaltest(ginis)
        print(f"  D'Agostino:  stat={stat2:.4f}, p={p2:.4e}")
    except:
        pass
    
    # Histogram (text-based)
    print(f"\n── Histogram ──")
    bin_edges = np.arange(0.55, 1.01, 0.025)
    counts, _ = np.histogram(ginis, bins=bin_edges)
    max_count = max(counts)
    for i in range(len(counts)):
        lo, hi = bin_edges[i], bin_edges[i+1]
        bar = "█" * int(50 * counts[i] / max_count) if max_count > 0 else ""
        print(f"  [{lo:.3f},{hi:.3f}) {counts[i]:>4d} {bar}")
    
    # Size dependence
    print(f"\n── Size Dependence ──")
    corr_size = np.corrcoef(n_res, ginis)[0, 1]
    print(f"  Pearson r(n_res, Gini) = {corr_size:.4f}")
    
    # Bin by size
    size_bins = [(50, 100), (100, 150), (150, 200), (200, 300), (300, 500)]
    print(f"  {'Size range':>12s} {'N':>5s} {'mean G':>8s} {'std G':>8s}")
    print("  " + "─"*36)
    for lo, hi in size_bins:
        m = (n_res >= lo) & (n_res < hi)
        n = np.sum(m)
        if n < 5: continue
        print(f"  [{lo:>3d},{hi:>3d}) {n:>5d} {np.mean(ginis[m]):>8.4f} {np.std(ginis[m]):>8.4f}")
    
    # Mean curvature vs Gini
    print(f"\n── Mean κ² vs Gini ──")
    corr_k = np.corrcoef(mean_k, ginis)[0, 1]
    print(f"  Pearson r(mean_κ², Gini) = {corr_k:.4f}")
    
    # Where is the peak?
    print(f"\n── Mode Detection ──")
    # KDE-like approach: smoothed histogram
    fine_bins = np.arange(0.60, 1.0, 0.005)
    fine_counts, _ = np.histogram(ginis, bins=fine_bins)
    # Smooth with a 5-bin window
    kernel = np.ones(5) / 5
    smoothed = np.convolve(fine_counts, kernel, mode='same')
    peak_idx = np.argmax(smoothed)
    peak_gini = (fine_bins[peak_idx] + fine_bins[peak_idx + 1]) / 2
    print(f"  Distribution peak (mode): {peak_gini:.3f}")
    
    # Check for bimodality
    # Simple test: is there a valley between two peaks?
    peaks = []
    for i in range(1, len(smoothed) - 1):
        if smoothed[i] > smoothed[i-1] and smoothed[i] > smoothed[i+1]:
            peaks.append((fine_bins[i] + fine_bins[i+1]) / 2)
    print(f"  Detected peaks: {len(peaks)}")
    for p in peaks[:5]:
        print(f"    Peak at G ≈ {p:.3f}")
    
    # Attractor position comparison
    print(f"\n── Attractor Comparison ──")
    print(f"  Census mean:            {np.mean(ginis):.4f}")
    print(f"  Census median:          {np.median(ginis):.4f}")
    print(f"  Census mode:            {peak_gini:.3f}")
    print(f"  Gradient zero-crossing: 0.828 (from apo-holo regression)")
    print(f"  Gradient best threshold: 0.863 (from binary classification)")
    print(f"  AlphaFold mean:         0.835 (from prior analysis)")
    
    # Export CSV
    csv_path = out_dir / "gini_census.csv"
    with open(csv_path, "w") as f:
        f.write("pdb,chain,gini,n_res,n_kappa,mean_kappa,median_kappa,max_kappa\n")
        for r in sorted(results, key=lambda x: x["gini"]):
            f.write(f"{r['pdb']},{r['chain']},{r['gini']:.6f},{r['n_res']},"
                    f"{r['n_kappa']},{r['mean_kappa']:.6f},"
                    f"{r['median_kappa']:.6f},{r['max_kappa']:.6f}\n")
    print(f"\n  CSV: {csv_path}")
    
    # Export JSON
    json_path = out_dir / "gini_census.json"
    export = {
        "n_proteins": len(results),
        "mean": float(np.mean(ginis)),
        "median": float(np.median(ginis)),
        "std": float(np.std(ginis)),
        "mode": float(peak_gini),
        "iqr": [float(np.percentile(ginis, 25)), float(np.percentile(ginis, 75))],
        "ci95": [float(np.percentile(ginis, 2.5)), float(np.percentile(ginis, 97.5))],
        "skewness": float(sk),
        "kurtosis": float(ku),
        "corr_size": float(corr_size),
        "corr_mean_kappa": float(corr_k),
        "histogram": {
            "bin_edges": bin_edges.tolist(),
            "counts": counts.tolist(),
        },
        "results": results,
    }
    with open(json_path, "w") as f:
        json.dump(export, f, indent=2)
    print(f"  JSON: {json_path}")


if __name__ == "__main__":
    main()

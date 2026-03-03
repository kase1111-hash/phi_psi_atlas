#!/usr/bin/env python3
"""
Sensitivity Tests for Crystal Atlas v2
=======================================
Test 1: Quadrant threshold sensitivity
    - Re-classify residues using P30, P40, P50 (median), P60, P70
    - Check if Pin > Beam > Loop > Spring hierarchy survives
    - Check if size-independence holds

Test 2: Resolution stratification
    - Parse resolution from PDB REMARK 2 records
    - Compare Gini distributions: ≤1.5 Å vs 1.5–2.0 Å vs 2.0–2.5 Å

Usage:
    python sensitivity_tests.py --pdb-dir gradient_v3/pdbs

No network access needed. Runs on cached PDB files.
"""

import numpy as np
import json, sys, os
from pathlib import Path
from collections import defaultdict
import argparse


# ═══════════════════════════════════════════════════════════════════
# BACKBONE PARSING (identical to atlas)
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


def gini(values):
    v = np.sort(np.array(values))
    n = len(v); s = np.sum(v)
    if s < 1e-15 or n < 2:
        return 0.0
    return (2*np.sum(np.arange(1, n+1)*v) / (n*s)) - (n+1)/n


# ═══════════════════════════════════════════════════════════════════
# RAMACHANDRAN CLASSIFICATION (corrected — negate both phi and psi)
# ═══════════════════════════════════════════════════════════════════

def is_structured(phi_rad, psi_rad):
    """Classify residue as structured (alpha/beta/left-alpha) or flexible."""
    # Negate both: our dihedral convention is -standard
    p = -np.degrees(phi_rad)
    s = -np.degrees(psi_rad)
    # Alpha-helix
    if -180 <= p <= -30 and -70 <= s <= -15:
        return True
    if -180 <= p <= -150 and -70 <= s <= -15:
        return True
    # Beta-sheet
    if -180 <= p <= -60 and 90 <= s <= 180:
        return True
    if -180 <= p <= -60 and -180 <= s <= -120:
        return True
    # Left-handed alpha
    if 30 <= p <= 120 and 15 <= s <= 90:
        return True
    return False


# ═══════════════════════════════════════════════════════════════════
# RESOLUTION PARSER
# ═══════════════════════════════════════════════════════════════════

def parse_resolution(filepath):
    """Extract resolution from PDB REMARK 2 record."""
    with open(filepath) as f:
        for line in f:
            if line.startswith("REMARK   2 RESOLUTION."):
                # Format: REMARK   2 RESOLUTION.    1.80 ANGSTROMS.
                parts = line.split()
                for i, p in enumerate(parts):
                    if p == "RESOLUTION.":
                        # Next token should be the number
                        if i + 1 < len(parts):
                            try:
                                return float(parts[i+1])
                            except ValueError:
                                pass
                # Try alternate parsing
                try:
                    val = line[22:30].strip()
                    return float(val)
                except:
                    pass
    return None


# ═══════════════════════════════════════════════════════════════════
# TEST 1: THRESHOLD SENSITIVITY
# ═══════════════════════════════════════════════════════════════════

def classify_at_threshold(ksq_dict, dih_dict, percentile):
    """Classify residues using a given percentile as the high/low κ² boundary."""
    if not ksq_dict:
        return {}
    kvals = np.array(list(ksq_dict.values()))
    threshold = np.percentile(kvals, percentile)

    counts = {"Beam": 0, "Loop": 0, "Pin": 0, "Spring": 0}
    for rn, ksq in ksq_dict.items():
        if rn not in dih_dict:
            continue
        phi, psi = dih_dict[rn]
        if np.isnan(phi) or np.isnan(psi):
            continue
        structured = is_structured(phi, psi)
        high_k = ksq > threshold

        if high_k and structured:
            counts["Pin"] += 1
        elif high_k and not structured:
            counts["Spring"] += 1
        elif not high_k and structured:
            counts["Beam"] += 1
        else:
            counts["Loop"] += 1
    return counts


def run_threshold_test(proteins):
    """
    Test quadrant taxonomy at P30, P40, P50, P60, P70 thresholds.
    """
    print("=" * 70)
    print("  TEST 1: QUADRANT THRESHOLD SENSITIVITY")
    print("=" * 70)

    percentiles = [30, 40, 50, 60, 70]

    # Global accumulators per percentile
    global_counts = {p: {"Beam": 0, "Loop": 0, "Pin": 0, "Spring": 0} for p in percentiles}
    # Per-protein fracs for size-independence test
    per_protein = {p: [] for p in percentiles}

    for prot in proteins:
        ksq = prot["ksq"]
        dih = prot["dih"]
        n_res = prot["n_res"]
        if len(ksq) < 10:
            continue

        for pctl in percentiles:
            counts = classify_at_threshold(ksq, dih, pctl)
            total = sum(counts.values())
            if total < 5:
                continue
            for k in counts:
                global_counts[pctl][k] += counts[k]
            fracs = {k: v/total for k, v in counts.items()}
            per_protein[pctl].append({"n_res": n_res, "fracs": fracs})

    # Report global fractions
    print(f"\n── Global Quadrant Fractions ──")
    print(f"  {'Pctl':>5s}  {'Beam':>6s}  {'Loop':>6s}  {'Pin':>6s}  {'Spr':>6s}  {'Hierarchy':>30s}")
    print("  " + "─" * 65)

    for pctl in percentiles:
        gc = global_counts[pctl]
        total = sum(gc.values())
        if total == 0:
            continue
        fracs = {k: gc[k]/total for k in gc}
        ranked = sorted(fracs.keys(), key=lambda k: -fracs[k])
        hierarchy = " > ".join(f"{k}({100*fracs[k]:.0f}%)" for k in ranked)
        marker = " ✓" if ranked[0] == "Pin" and ranked[1] == "Beam" else " ✗"
        print(f"  P{pctl:>3d}  {100*fracs['Beam']:>5.1f}%  {100*fracs['Loop']:>5.1f}%  "
              f"{100*fracs['Pin']:>5.1f}%  {100*fracs['Spring']:>5.1f}%  {hierarchy}{marker}")

    # Size-independence at each threshold
    print(f"\n── Size Independence (Beam fraction across size bins) ──")
    print(f"  {'Pctl':>5s}  {'50-100':>7s}  {'100-200':>7s}  {'200-300':>7s}  {'300-500':>7s}  {'500+':>7s}  {'Range':>7s}")
    print("  " + "─" * 55)

    size_bins = [(50, 100), (100, 200), (200, 300), (300, 500), (500, 2000)]
    for pctl in percentiles:
        vals = []
        for lo, hi in size_bins:
            subset = [p["fracs"]["Beam"] for p in per_protein[pctl] if lo <= p["n_res"] < hi]
            if len(subset) >= 10:
                vals.append(np.mean(subset))
            else:
                vals.append(np.nan)
        valid = [v for v in vals if not np.isnan(v)]
        rng = max(valid) - min(valid) if len(valid) >= 2 else np.nan
        print(f"  P{pctl:>3d}  " + "  ".join(f"{v:>6.3f}" if not np.isnan(v) else "    --" for v in vals) +
              f"  {rng:>6.3f}")

    print(f"\n── Size Independence (Pin fraction across size bins) ──")
    print(f"  {'Pctl':>5s}  {'50-100':>7s}  {'100-200':>7s}  {'200-300':>7s}  {'300-500':>7s}  {'500+':>7s}  {'Range':>7s}")
    print("  " + "─" * 55)

    for pctl in percentiles:
        vals = []
        for lo, hi in size_bins:
            subset = [p["fracs"]["Pin"] for p in per_protein[pctl] if lo <= p["n_res"] < hi]
            if len(subset) >= 10:
                vals.append(np.mean(subset))
            else:
                vals.append(np.nan)
        valid = [v for v in vals if not np.isnan(v)]
        rng = max(valid) - min(valid) if len(valid) >= 2 else np.nan
        print(f"  P{pctl:>3d}  " + "  ".join(f"{v:>6.3f}" if not np.isnan(v) else "    --" for v in vals) +
              f"  {rng:>6.3f}")

    # Per-protein mean fracs
    print(f"\n── Per-Protein Mean Fractions ──")
    print(f"  {'Pctl':>5s}  {'Beam':>12s}  {'Loop':>12s}  {'Pin':>12s}  {'Spring':>12s}")
    print("  " + "─" * 55)

    for pctl in percentiles:
        pp = per_protein[pctl]
        if not pp:
            continue
        for q in ["Beam", "Loop", "Pin", "Spring"]:
            pass
        b = [p["fracs"]["Beam"] for p in pp]
        l = [p["fracs"]["Loop"] for p in pp]
        pi = [p["fracs"]["Pin"] for p in pp]
        s = [p["fracs"]["Spring"] for p in pp]
        print(f"  P{pctl:>3d}  {np.mean(b):>.3f}±{np.std(b):>.3f}  "
              f"{np.mean(l):>.3f}±{np.std(l):>.3f}  "
              f"{np.mean(pi):>.3f}±{np.std(pi):>.3f}  "
              f"{np.mean(s):>.3f}±{np.std(s):>.3f}")


# ═══════════════════════════════════════════════════════════════════
# TEST 2: RESOLUTION STRATIFICATION
# ═══════════════════════════════════════════════════════════════════

def run_resolution_test(proteins):
    """
    Stratify by resolution and compare Gini distributions.
    """
    print("\n" + "=" * 70)
    print("  TEST 2: RESOLUTION STRATIFICATION")
    print("=" * 70)

    # Bin by resolution
    bins = [
        ("≤1.5 Å", 0.0, 1.5),
        ("1.5–2.0 Å", 1.5, 2.0),
        ("2.0–2.5 Å", 2.0, 2.5),
        (">2.5 Å", 2.5, 10.0),
    ]

    has_res = [(p, p["resolution"]) for p in proteins if p["resolution"] is not None]
    no_res = [p for p in proteins if p["resolution"] is None]

    print(f"\n  Proteins with resolution: {len(has_res)}")
    print(f"  Proteins without resolution (NMR etc.): {len(no_res)}")

    if len(has_res) < 50:
        print("  [WARN] Too few proteins with resolution data")
        return

    print(f"\n  {'Bin':>12s}  {'N':>6s}  {'Mean Gini':>10s}  {'Median':>8s}  {'SD':>8s}  {'IQR':>16s}")
    print("  " + "─" * 65)

    all_ginis = {}
    for label, lo, hi in bins:
        ginis = np.array([p["gini"] for p, r in has_res if lo < r <= hi])
        if len(ginis) < 5:
            print(f"  {label:>12s}  {len(ginis):>6d}  {'--':>10s}")
            continue
        iqr = f"[{np.percentile(ginis,25):.3f}, {np.percentile(ginis,75):.3f}]"
        print(f"  {label:>12s}  {len(ginis):>6d}  {np.mean(ginis):>10.4f}  "
              f"{np.median(ginis):>8.4f}  {np.std(ginis):>8.4f}  {iqr:>16s}")
        all_ginis[label] = ginis

    # Statistical comparison between best two bins
    if "≤1.5 Å" in all_ginis and "1.5–2.0 Å" in all_ginis:
        g1 = all_ginis["≤1.5 Å"]
        g2 = all_ginis["1.5–2.0 Å"]
        diff = np.mean(g1) - np.mean(g2)
        # Welch's t-test (no scipy dependency)
        se = np.sqrt(np.var(g1)/len(g1) + np.var(g2)/len(g2))
        t_stat = diff / se if se > 0 else 0
        # Approximate p from t for large N
        from math import erfc, sqrt
        p_approx = erfc(abs(t_stat) / sqrt(2))
        print(f"\n  ≤1.5 Å vs 1.5–2.0 Å:")
        print(f"    ΔMean = {diff:+.4f}")
        print(f"    t = {t_stat:.2f}, p ≈ {p_approx:.4f}")
        if p_approx > 0.05:
            print(f"    → No significant difference (p > 0.05)")
        else:
            print(f"    → Significant difference (p < 0.05)")
            print(f"    → Effect size: {abs(diff):.4f} Gini units")

    # Quadrant fracs by resolution
    print(f"\n  Quadrant fractions by resolution:")
    print(f"  {'Bin':>12s}  {'N':>6s}  {'Beam':>6s}  {'Loop':>6s}  {'Pin':>6s}  {'Spr':>6s}")
    print("  " + "─" * 45)

    for label, lo, hi in bins:
        subset = [p for p, r in has_res if lo < r <= hi]
        if len(subset) < 10:
            continue
        beam = np.mean([p["quad_fracs"]["Beam"] for p in subset])
        loop = np.mean([p["quad_fracs"]["Loop"] for p in subset])
        pin = np.mean([p["quad_fracs"]["Pin"] for p in subset])
        spring = np.mean([p["quad_fracs"]["Spring"] for p in subset])
        print(f"  {label:>12s}  {len(subset):>6d}  {beam:>5.3f}  {loop:>5.3f}  "
              f"{pin:>5.3f}  {spring:>5.3f}")


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-dir", default="gradient_v3/pdbs")
    args = parser.parse_args()

    pdb_dir = Path(args.pdb_dir)
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
    print("║  SENSITIVITY TESTS — Threshold + Resolution                ║")
    print(f"║  {len(pdb_files):>5d} PDB files                                        ║")
    print("╚══════════════════════════════════════════════════════════════╝")

    # Parse all proteins
    proteins = []
    failed = 0

    for idx, fp in enumerate(pdb_files):
        if idx % 500 == 0 and idx > 0:
            print(f"  {idx}/{len(pdb_files)} ({len(proteins)} valid)...")

        pdb_id = fp.stem
        residues, chain = parse_backbone(str(fp))
        if len(residues) < 15:
            failed += 1; continue

        dih = get_dihedrals(residues)
        ksq = get_kappa_sq(dih)
        if len(ksq) < 10:
            failed += 1; continue

        # Gini
        kvals = list(ksq.values())
        g = gini(kvals)

        # Resolution
        res = parse_resolution(str(fp))

        # Quadrant classification at median (for resolution test)
        median_k = np.median(kvals)
        counts = {"Beam": 0, "Loop": 0, "Pin": 0, "Spring": 0}
        for rn, kv in ksq.items():
            if rn not in dih:
                continue
            phi, psi = dih[rn]
            if np.isnan(phi) or np.isnan(psi):
                continue
            structured = is_structured(phi, psi)
            high_k = kv > median_k
            if high_k and structured:
                counts["Pin"] += 1
            elif high_k and not structured:
                counts["Spring"] += 1
            elif not high_k and structured:
                counts["Beam"] += 1
            else:
                counts["Loop"] += 1
        total = sum(counts.values())
        fracs = {k: v/total if total > 0 else 0 for k, v in counts.items()}

        proteins.append({
            "pdb": pdb_id,
            "n_res": len(residues),
            "gini": g,
            "resolution": res,
            "ksq": ksq,
            "dih": dih,
            "quad_counts": counts,
            "quad_fracs": fracs,
        })

    print(f"\n  {len(pdb_files)} files → {len(proteins)} valid ({failed} failed)")
    print(f"  With resolution: {sum(1 for p in proteins if p['resolution'] is not None)}")
    print(f"  Without resolution: {sum(1 for p in proteins if p['resolution'] is None)}")

    # Run tests
    run_threshold_test(proteins)
    run_resolution_test(proteins)

    # Summary
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print("\n  Test 1 (Threshold): Does Pin > Beam > Loop > Spring hold at all")
    print("  percentiles? Are fractions size-independent?")
    print("\n  Test 2 (Resolution): Does Gini shift between resolution bins?")
    print("  Do quadrant fractions depend on resolution?")
    print("\n  If both pass → taxonomy is robust to methodological choices.")


if __name__ == "__main__":
    main()

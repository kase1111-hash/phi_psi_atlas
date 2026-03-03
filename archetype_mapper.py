#!/usr/bin/env python3
"""
Archetype Mapper — Cross-Species Mechanical Fingerprint Taxonomy
================================================================
Two modes:
  1. CENSUS: Single-structure κ² baseline (works with AlphaFold or crystal)
  2. COMPARE: Two-structure transition fingerprint (needs crystal pairs)

Plus:
  3. BATCH CENSUS: Scan all structures in a directory
  4. ARCHETYPE MAP: Compare fingerprints across species for same protein family

Uses strain_analyzer.py v15 as the math engine.
"""

import sys
import os
import json
import csv
import glob
import numpy as np
from pathlib import Path
from collections import defaultdict
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from strain_analyzer import (
    parse_mmcif, compute_torsion_angles, compute_kappa_squared,
    gini_coefficient, get_primary_ligand, rama_class
)

# ═══════════════════════════════════════════════════════════════
# MODE 1: CENSUS — Single structure baseline
# ═══════════════════════════════════════════════════════════════

def census(cif_path, chain="A", label=None):
    """Compute κ² baseline for a single structure.
    Returns dict with Gini, κ² stats, SS breakdown, top residues.
    Works with crystal structures AND AlphaFold models."""
    
    if label is None:
        label = Path(cif_path).stem
    
    try:
        bb, bf, lig = parse_mmcif(str(cif_path), chain)
    except Exception as e:
        return {"label": label, "chain": chain, "error": str(e)}
    
    n_res = len(set(a['resnum'] for a in bb))
    phi, psi, res_map = compute_torsion_angles(bb)
    ksq = compute_kappa_squared(phi, psi)
    
    # Filter zeros
    valid = {r: v for r, v in ksq.items() if v > 0}
    if len(valid) < 10:
        return {"label": label, "chain": chain, "error": "too few valid residues"}
    
    vals = list(valid.values())
    residues = sorted(valid.keys())
    
    # Gini
    g = gini_coefficient(vals)
    
    # SS composition
    ss_count = {"alpha": 0, "beta": 0, "loop": 0}
    for r in residues:
        ss = rama_class(phi.get(r, np.nan), psi.get(r, np.nan))
        if ss in ss_count:
            ss_count[ss] += 1
    
    # Top 10 highest curvature residues
    ranked = sorted(residues, key=lambda r: -valid[r])
    top10 = []
    for r in ranked[:10]:
        ss = rama_class(phi.get(r, np.nan), psi.get(r, np.nan))
        top10.append({
            "resnum": r, "kappa2": round(valid[r], 1), "ss": ss
        })
    
    # Ligands
    lig_coords, lig_name = get_primary_ligand(lig)
    
    # κ² by SS type
    kappa_by_ss = {"alpha": [], "beta": [], "loop": []}
    for r in residues:
        ss = rama_class(phi.get(r, np.nan), psi.get(r, np.nan))
        if ss in kappa_by_ss:
            kappa_by_ss[ss].append(valid[r])
    
    ss_means = {}
    for ss_type, vals_ss in kappa_by_ss.items():
        if vals_ss:
            ss_means[ss_type] = round(np.mean(vals_ss), 1)
    
    return {
        "label": label,
        "chain": chain,
        "pdb_id": Path(cif_path).stem.split("-")[0].upper(),
        "n_residues": n_res,
        "n_valid": len(valid),
        "gini": round(g, 4),
        "kappa2_mean": round(np.mean(vals), 1),
        "kappa2_median": round(np.median(vals), 1),
        "kappa2_max": round(max(vals), 1),
        "kappa2_std": round(np.std(vals), 1),
        "ss_composition": ss_count,
        "ss_fraction": {
            "alpha": round(ss_count["alpha"] / len(residues), 2),
            "beta": round(ss_count["beta"] / len(residues), 2),
            "loop": round(ss_count["loop"] / len(residues), 2),
        },
        "kappa2_by_ss": ss_means,
        "ligand": lig_name,
        "ligand_atoms": len(lig_coords) if lig_coords else 0,
        "top10_residues": top10,
        # Raw profile for cross-comparison
        "_profile": {r: valid[r] for r in residues},
    }


# ═══════════════════════════════════════════════════════════════
# MODE 2: COMPARE — Two-structure transition fingerprint
# ═══════════════════════════════════════════════════════════════

def compare(ref_path, pert_path, chain="A", label=None):
    """Full transition fingerprint between two structures.
    Returns the standard 8-step analysis as a dict."""
    
    if label is None:
        label = f"{Path(ref_path).stem}_vs_{Path(pert_path).stem}"
    
    try:
        bb_a, bf_a, lig_a = parse_mmcif(str(ref_path), chain)
        bb_d, bf_d, lig_d = parse_mmcif(str(pert_path), chain)
    except Exception as e:
        return {"label": label, "error": str(e)}
    
    n_a = len(set(a['resnum'] for a in bb_a))
    n_d = len(set(a['resnum'] for a in bb_d))
    
    phi_a, psi_a, res_a = compute_torsion_angles(bb_a)
    phi_d, psi_d, res_d = compute_torsion_angles(bb_d)
    ksq_a = compute_kappa_squared(phi_a, psi_a)
    ksq_d = compute_kappa_squared(phi_d, psi_d)
    
    common_all = sorted(set(ksq_a.keys()) & set(ksq_d.keys()))
    common = [r for r in common_all if ksq_a[r] > 0 or ksq_d[r] > 0]
    
    if len(common) < 10:
        return {"label": label, "error": "too few common residues"}
    
    k_a = [ksq_a[r] for r in common]
    k_d = [ksq_d[r] for r in common]
    g_a = gini_coefficient(k_a)
    g_d = gini_coefficient(k_d)
    dg = g_d - g_a
    
    lig_coords, lig_name = get_primary_ligand(lig_d)
    
    # Delta kappa squared
    dksq = {r: ksq_d[r] - ksq_a[r] for r in common}
    ranked = sorted(common, key=lambda r: -abs(dksq[r]))
    top20 = ranked[:20]
    
    # Pocket
    pocket = set()
    if lig_coords:
        for r in common:
            ca = res_d.get(r, {}).get('CA', res_a.get(r, {}).get('CA'))
            if ca is not None:
                if min(np.linalg.norm(ca - l) for l in lig_coords) < 10.0:
                    pocket.add(r)
    
    n_local = sum(1 for r in top20 if r in pocket)
    n_distal = 20 - n_local
    
    # Mode
    local_res = [r for r in common if r in pocket]
    distal_res = [r for r in common if r not in pocket]
    mode = "N/A"
    lam = 0.0
    dg_l = 0.0
    dg_dv = 0.0
    
    if len(local_res) >= 5 and len(distal_res) >= 5:
        g_la = gini_coefficient([ksq_a[r] for r in local_res])
        g_ld = gini_coefficient([ksq_d[r] for r in local_res])
        g_da = gini_coefficient([ksq_a[r] for r in distal_res])
        g_dd = gini_coefficient([ksq_d[r] for r in distal_res])
        dg_l = g_ld - g_la
        dg_dv = g_dd - g_da
        
        if abs(dg_l) < 0.01 and abs(dg_dv) < 0.01:
            mode = "MINIMAL"
        elif np.sign(dg_l) == np.sign(dg_dv):
            mode = "UNIFORM"
        elif dg_l > 0 and dg_dv < 0:
            mode = "CONCENTRATOR"
        elif dg_l < 0 and dg_dv > 0:
            mode = "EXPORTER"
        else:
            mode = "MIXED"
        lam = dg_l / dg if abs(dg) > 0.001 else 0
    
    # SS
    ss_count = {"alpha": 0, "beta": 0, "loop": 0}
    for r in common:
        ss = rama_class(phi_a.get(r, np.nan), psi_a.get(r, np.nan))
        if ss in ss_count:
            ss_count[ss] += 1
    
    # Hotspots
    hotspots = []
    for i, r in enumerate(top20):
        ss = rama_class(phi_a.get(r, np.nan), psi_a.get(r, np.nan))
        zone = "LOCAL" if r in pocket else "DISTAL"
        direction = "stiffen" if dksq[r] > 0 else "loosen"
        dist = None
        if lig_coords:
            ca = res_d.get(r, {}).get('CA', res_a.get(r, {}).get('CA'))
            if ca is not None:
                dist = round(min(np.linalg.norm(ca - l) for l in lig_coords), 1)
        hotspots.append({
            "rank": i + 1, "resnum": r, "ss": ss,
            "delta_kappa2": round(dksq[r], 1),
            "direction": direction, "zone": zone,
            "ligand_dist": dist
        })
    
    # Direction counts
    n_stiffen = sum(1 for h in hotspots if h["direction"] == "stiffen")
    n_loosen = 20 - n_stiffen
    
    return {
        "label": label,
        "chain": chain,
        "ref_pdb": Path(ref_path).stem.upper(),
        "pert_pdb": Path(pert_path).stem.upper(),
        "n_ref": n_a,
        "n_pert": n_d,
        "n_common": len(common),
        "gini_ref": round(g_a, 4),
        "gini_pert": round(g_d, 4),
        "delta_gini": round(dg, 4),
        "direction": "STIFFENED" if dg > 0 else "LOOSENED",
        "kappa2_mean_ref": round(np.mean(k_a), 1),
        "kappa2_mean_pert": round(np.mean(k_d), 1),
        "kappa2_ratio": round(np.mean(k_d) / np.mean(k_a), 2) if np.mean(k_a) > 0 else None,
        "ss_composition": ss_count,
        "ligand": lig_name,
        "ligand_atoms": len(lig_coords) if lig_coords else 0,
        "pocket_size": len(pocket),
        "pocket_fraction": round(len(pocket) / len(common), 3) if pocket else 0,
        "top20_local": n_local,
        "top20_distal": n_distal,
        "distal_pct": round(100 * n_distal / 20, 0),
        "mode": mode,
        "lambda": round(lam, 2),
        "delta_gini_local": round(dg_l, 4),
        "delta_gini_distal": round(dg_dv, 4),
        "hotspot_stiffen": n_stiffen,
        "hotspot_loosen": n_loosen,
        "hotspots": hotspots,
        # Raw profiles for correlation
        "_profile_ref": {r: ksq_a[r] for r in common},
        "_profile_pert": {r: ksq_d[r] for r in common},
        "_delta_profile": {r: dksq[r] for r in common},
    }


# ═══════════════════════════════════════════════════════════════
# MODE 3: BATCH CENSUS — Scan all structures in a directory
# ═══════════════════════════════════════════════════════════════

def batch_census(directory, chain="A", pattern="*.cif", workers=1,
                 output_csv=None, resume=True):
    """Run census on every CIF file in a directory.
    
    Args:
        directory: Path to directory of CIF files
        chain: Chain ID to analyze
        pattern: Glob pattern for files
        workers: Number of parallel threads (1=serial, 4-16 recommended for I/O bound)
        output_csv: If set, stream results to CSV incrementally (enables resume)
        resume: If True and output_csv exists, skip already-processed files
    
    Returns:
        List of census result dicts
    """
    import concurrent.futures
    import threading
    
    files = sorted(glob.glob(os.path.join(directory, pattern)))
    # Filter duplicates
    files = [f for f in files if "__1_" not in Path(f).stem]
    
    # ── Resume support: load already-done labels from existing CSV ──
    done_labels = set()
    existing_rows = []
    if resume and output_csv and os.path.exists(output_csv):
        try:
            with open(output_csv, 'r') as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    done_labels.add(row.get("label", ""))
                    existing_rows.append(row)
            print(f"  Resume: {len(done_labels):,} already done, skipping")
        except:
            pass
    
    todo = []
    for f in files:
        stem = Path(f).stem
        if stem in done_labels:
            continue
        todo.append(f)
    
    if not todo:
        print(f"  All {len(files)} files already processed")
        return existing_rows
    
    print(f"  Processing {len(todo):,} files ({len(done_labels):,} cached) with {workers} worker(s)")
    
    # ── Thread-safe state ──
    lock = threading.Lock()
    counter = {"done": 0, "ok": 0, "err": 0, "total": len(todo)}
    results = []
    
    # ── CSV writer setup for streaming ──
    csv_fh = None
    csv_writer = None
    csv_fields = [
        "label", "pdb_id", "chain", "n_residues", "n_valid", "gini",
        "kappa2_mean", "kappa2_median", "kappa2_max", "kappa2_std",
        "alpha", "beta", "loop", "frac_alpha", "frac_beta",
        "kappa2_alpha", "kappa2_beta", "kappa2_loop",
        "ligand", "ligand_atoms", "top1_res", "top1_kappa2"
    ]
    
    if output_csv:
        write_header = not os.path.exists(output_csv) or len(done_labels) == 0
        csv_fh = open(output_csv, 'a' if done_labels else 'w', newline='')
        csv_writer = csv.DictWriter(csv_fh, fieldnames=csv_fields)
        if write_header:
            csv_writer.writeheader()
            csv_fh.flush()
    
    def _process_one(filepath):
        """Worker function — runs census on one file."""
        stem = Path(filepath).stem
        try:
            r = census(filepath, chain=chain)
            is_ok = "error" not in r
            
            with lock:
                counter["done"] += 1
                if is_ok:
                    counter["ok"] += 1
                    results.append(r)
                    # Stream to CSV
                    if csv_writer:
                        csv_writer.writerow(census_to_row(r))
                        # Flush periodically
                        if counter["done"] % 50 == 0:
                            csv_fh.flush()
                else:
                    counter["err"] += 1
                    results.append(r)
                
                done = counter["done"]
                total = counter["total"]
                
                # Progress: print every N files based on total
                if total > 1000:
                    interval = 100
                elif total > 100:
                    interval = 10
                else:
                    interval = 1
                
                if done % interval == 0 or done == total:
                    pct = 100 * done / total
                    rate = ""
                    print(f"  [{done:>{len(str(total))}}/{total}] {pct:5.1f}%  ok={counter['ok']}  err={counter['err']}  {stem}")
            
            return r
        except Exception as e:
            with lock:
                counter["done"] += 1
                counter["err"] += 1
                err_result = {"label": stem, "chain": chain, "error": str(e)}
                results.append(err_result)
            return err_result
    
    # ── Execute ──
    import time
    t0 = time.time()
    
    if workers <= 1:
        # Serial mode
        for f in todo:
            _process_one(f)
    else:
        # Threaded mode — I/O bound so threads are fine (GIL released during file reads)
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_process_one, f): f for f in todo}
            # Just wait for all to complete — progress printed inside _process_one
            concurrent.futures.wait(futures)
    
    elapsed = time.time() - t0
    rate = len(todo) / elapsed if elapsed > 0 else 0
    
    # Flush and close CSV
    if csv_fh:
        csv_fh.flush()
        csv_fh.close()
    
    print(f"\n  Done: {counter['ok']:,} ok, {counter['err']:,} errors in {elapsed:.1f}s ({rate:.1f} files/sec)")
    if output_csv:
        total_rows = len(done_labels) + counter['ok']
        print(f"  CSV: {total_rows:,} total rows in {output_csv}")
    
    return results


def batch_compare(pairs, chain="A", workers=1, output_csv=None):
    """Run compare() on a list of (ref_path, pert_path, label) tuples.
    
    Args:
        pairs: List of (ref_path, pert_path, label) tuples
        chain: Chain ID
        workers: Number of parallel threads
        output_csv: If set, stream results to CSV
    
    Returns:
        List of comparison result dicts
    """
    import concurrent.futures
    import threading
    
    lock = threading.Lock()
    counter = {"done": 0, "ok": 0, "err": 0}
    results = []
    
    csv_fh = None
    csv_writer = None
    cmp_fields = [
        "label", "ref_pdb", "pert_pdb", "chain", "n_common",
        "gini_ref", "gini_pert", "delta_gini", "direction",
        "kappa2_mean_ref", "kappa2_mean_pert", "kappa2_ratio",
        "ligand", "pocket_size", "pocket_fraction",
        "distal_pct", "mode", "lambda", "dg_local", "dg_distal",
        "n_stiffen", "n_loosen", "archetype"
    ]
    
    if output_csv:
        csv_fh = open(output_csv, 'w', newline='')
        csv_writer = csv.DictWriter(csv_fh, fieldnames=cmp_fields)
        csv_writer.writeheader()
    
    def _process_pair(ref_path, pert_path, label):
        try:
            r = compare(ref_path, pert_path, chain=chain, label=label)
            is_ok = "error" not in r
            
            with lock:
                counter["done"] += 1
                if is_ok:
                    counter["ok"] += 1
                    results.append(r)
                    if csv_writer:
                        csv_writer.writerow(compare_to_row(r))
                else:
                    counter["err"] += 1
                    results.append(r)
                
                arch = classify_archetype(r) if is_ok else "ERROR"
                dg = r.get("delta_gini", 0)
                print(f"  [{counter['done']}/{len(pairs)}] {label:30s} ΔG={dg:+.4f} → {arch}")
            return r
        except Exception as e:
            with lock:
                counter["done"] += 1
                counter["err"] += 1
            return {"label": label, "error": str(e)}
    
    import time
    t0 = time.time()
    
    if workers <= 1:
        for ref, pert, label in pairs:
            _process_pair(ref, pert, label)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_process_pair, r, p, l): l for r, p, l in pairs}
            concurrent.futures.wait(futures)
    
    elapsed = time.time() - t0
    
    if csv_fh:
        csv_fh.flush()
        csv_fh.close()
    
    print(f"\n  Done: {counter['ok']} ok, {counter['err']} errors in {elapsed:.1f}s")
    return results


# ═══════════════════════════════════════════════════════════════
# MODE 4: PROFILE CORRELATION — Compare κ² shapes
# ═══════════════════════════════════════════════════════════════

def profile_correlation(census_a, census_b):
    """Compute correlation between two κ² profiles.
    Aligns by residue number (for homologs with same numbering)
    or by rank order (for structural alignment)."""
    
    prof_a = census_a.get("_profile", {})
    prof_b = census_b.get("_profile", {})
    
    # Method 1: Residue-number alignment
    common = sorted(set(prof_a.keys()) & set(prof_b.keys()))
    rho_resnum = None
    if len(common) >= 20:
        va = [prof_a[r] for r in common]
        vb = [prof_b[r] for r in common]
        if np.std(va) > 0 and np.std(vb) > 0:
            rho_resnum = round(np.corrcoef(va, vb)[0, 1], 3)
    
    # Method 2: Rank-order alignment (for different numbering)
    sorted_a = sorted(prof_a.values())
    sorted_b = sorted(prof_b.values())
    min_len = min(len(sorted_a), len(sorted_b))
    rho_rank = None
    if min_len >= 20:
        # Trim to same length, comparing from highest
        va = sorted_a[-min_len:]
        vb = sorted_b[-min_len:]
        if np.std(va) > 0 and np.std(vb) > 0:
            rho_rank = round(np.corrcoef(va, vb)[0, 1], 3)
    
    # Method 3: Gini + mean comparison (coarsest)
    gini_diff = abs(census_a.get("gini", 0) - census_b.get("gini", 0))
    mean_ratio = None
    ma = census_a.get("kappa2_mean", 0)
    mb = census_b.get("kappa2_mean", 0)
    if ma > 0 and mb > 0:
        mean_ratio = round(max(ma, mb) / min(ma, mb), 2)
    
    return {
        "label_a": census_a.get("label", "?"),
        "label_b": census_b.get("label", "?"),
        "common_residues": len(common),
        "rho_resnum": rho_resnum,
        "rho_rank": rho_rank,
        "gini_diff": round(gini_diff, 4),
        "mean_ratio": mean_ratio,
    }


# ═══════════════════════════════════════════════════════════════
# MODE 5: ARCHETYPE CLASSIFIER
# ═══════════════════════════════════════════════════════════════

def classify_archetype(comparison_result):
    """Classify a transition into mechanical archetypes based on fingerprint."""
    
    r = comparison_result
    if "error" in r:
        return "ERROR"
    
    dg = r.get("delta_gini", 0)
    mode = r.get("mode", "N/A")
    distal = r.get("distal_pct", 0)
    lam = r.get("lambda", 0)
    ratio = r.get("kappa2_ratio", 1)
    dg_l = r.get("delta_gini_local", 0)
    dg_d = r.get("delta_gini_distal", 0)
    
    # Noise floor
    if abs(dg) < 0.015:
        if any(abs(h["delta_kappa2"]) > 500 for h in r.get("hotspots", [])):
            return "REDISTRIBUTOR"  # Zero-sum strain shuffle (Myosin-like)
        return "MINIMAL"  # No significant change
    
    # Strong stiffening
    if dg > 0.05:
        if ratio > 2.5 and distal > 70:
            return "GLOBAL_STIFFENER"  # Calmodulin-like ion trigger
        if mode == "CONCENTRATOR":
            return "CONCENTRATOR"  # Strain focuses at pocket
        if mode == "EXPORTER":
            return "EXPORTER"  # Pocket clamps, periphery loosens
        if mode == "UNIFORM" and distal > 60:
            return "UNIFORM_STIFFENER"  # Everything tightens
        return "STIFFENER"
    
    # Strong loosening
    if dg < -0.05:
        if distal > 90 and ratio > 2:
            return "ELASTIC_SPRING"  # GCK-like global discharge
        if mode == "UNIFORM":
            return "UNIFORM_LOOSENER"  # Everything relaxes
        return "LOOSENER"
    
    # Moderate changes
    if mode == "MIXED" or (dg_l * dg_d < 0):
        return "ALLOSTERIC_RELAY"  # AdK-like, local/distal diverge
    
    if distal > 80:
        return "DISTAL_BROADCASTER"  # Signal propagates far
    
    return "MODERATE"  # Doesn't fit clean archetype


# ═══════════════════════════════════════════════════════════════
# OUTPUT FORMATTERS
# ═══════════════════════════════════════════════════════════════

def census_to_row(c):
    """Convert census dict to flat CSV row."""
    if "error" in c:
        return {"label": c.get("label"), "error": c.get("error")}
    return {
        "label": c["label"],
        "pdb_id": c.get("pdb_id", ""),
        "chain": c["chain"],
        "n_residues": c["n_residues"],
        "n_valid": c["n_valid"],
        "gini": c["gini"],
        "kappa2_mean": c["kappa2_mean"],
        "kappa2_median": c["kappa2_median"],
        "kappa2_max": c["kappa2_max"],
        "kappa2_std": c["kappa2_std"],
        "alpha": c["ss_composition"]["alpha"],
        "beta": c["ss_composition"]["beta"],
        "loop": c["ss_composition"]["loop"],
        "frac_alpha": c["ss_fraction"]["alpha"],
        "frac_beta": c["ss_fraction"]["beta"],
        "kappa2_alpha": c["kappa2_by_ss"].get("alpha", ""),
        "kappa2_beta": c["kappa2_by_ss"].get("beta", ""),
        "kappa2_loop": c["kappa2_by_ss"].get("loop", ""),
        "ligand": c["ligand"],
        "ligand_atoms": c["ligand_atoms"],
        "top1_res": c["top10_residues"][0]["resnum"] if c["top10_residues"] else "",
        "top1_kappa2": c["top10_residues"][0]["kappa2"] if c["top10_residues"] else "",
    }


def write_census_csv(results, outpath):
    """Write batch census results to CSV."""
    good = [census_to_row(r) for r in results if "error" not in r]
    if not good:
        print("No valid results to write")
        return
    
    fields = list(good[0].keys())
    with open(outpath, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in good:
            w.writerow(row)
    print(f"Wrote {len(good)} rows to {outpath}")


def compare_to_row(c):
    """Convert comparison dict to flat CSV row."""
    if "error" in c:
        return {"label": c.get("label"), "error": c.get("error")}
    return {
        "label": c["label"],
        "ref_pdb": c["ref_pdb"],
        "pert_pdb": c["pert_pdb"],
        "chain": c["chain"],
        "n_common": c["n_common"],
        "gini_ref": c["gini_ref"],
        "gini_pert": c["gini_pert"],
        "delta_gini": c["delta_gini"],
        "direction": c["direction"],
        "kappa2_mean_ref": c["kappa2_mean_ref"],
        "kappa2_mean_pert": c["kappa2_mean_pert"],
        "kappa2_ratio": c["kappa2_ratio"],
        "ligand": c["ligand"],
        "pocket_size": c["pocket_size"],
        "pocket_fraction": c["pocket_fraction"],
        "distal_pct": c["distal_pct"],
        "mode": c["mode"],
        "lambda": c["lambda"],
        "dg_local": c["delta_gini_local"],
        "dg_distal": c["delta_gini_distal"],
        "n_stiffen": c["hotspot_stiffen"],
        "n_loosen": c["hotspot_loosen"],
        "archetype": classify_archetype(c),
    }


# ═══════════════════════════════════════════════════════════════
# CROSS-SPECIES COMPARISON
# ═══════════════════════════════════════════════════════════════

def cross_species_report(comparisons, family_name="Unknown"):
    """Given a list of compare() results for the same protein family
    across species, generate a cross-species archetype report."""
    
    lines = []
    lines.append(f"\n{'='*72}")
    lines.append(f"  CROSS-SPECIES ARCHETYPE MAP: {family_name}")
    lines.append(f"{'='*72}\n")
    
    for c in comparisons:
        if "error" in c:
            lines.append(f"  {c['label']}: ERROR — {c['error']}")
            continue
        
        arch = classify_archetype(c)
        lines.append(f"  {c['label']}")
        lines.append(f"    Residues: {c['n_common']} | SS: α={c['ss_composition']['alpha']} β={c['ss_composition']['beta']} L={c['ss_composition']['loop']}")
        lines.append(f"    ΔGini: {c['delta_gini']:+.4f} ({c['direction']}) | κ² ratio: {c['kappa2_ratio']}×")
        lines.append(f"    Distal: {c['distal_pct']:.0f}% | Mode: {c['mode']} | λ={c['lambda']:.2f}")
        lines.append(f"    ΔGini local: {c['delta_gini_local']:+.4f} | ΔGini distal: {c['delta_gini_distal']:+.4f}")
        lines.append(f"    ─→ ARCHETYPE: {arch}")
        lines.append("")
    
    # Summary
    archetypes = [classify_archetype(c) for c in comparisons if "error" not in c]
    if archetypes:
        unique = set(archetypes)
        lines.append(f"  VERDICT: {len(unique)} unique archetype(s) across {len(archetypes)} species")
        if len(unique) == 1:
            lines.append(f"  ─→ UNIVERSAL TOOL: All species use '{archetypes[0]}'")
        elif len(unique) == len(archetypes):
            lines.append(f"  ─→ SPECIES-SPECIFIC: Every species uses a different mechanical solution")
        else:
            lines.append(f"  ─→ MIXED: Some convergence, some divergence")
            for a in unique:
                count = archetypes.count(a)
                lines.append(f"     {a}: {count}/{len(archetypes)}")
    
    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════
# MAIN — CLI interface
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Archetype Mapper — Mechanical Fingerprint Taxonomy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single structure census
  python archetype_mapper.py census structure.cif

  # Batch census with 8 threads, resume-capable
  python archetype_mapper.py batch /path/to/cifs -w 8 -o results.csv

  # Resume interrupted batch (skips already-done files)
  python archetype_mapper.py batch /path/to/cifs -w 8 -o results.csv

  # Two-structure comparison
  python archetype_mapper.py compare apo.cif drug.cif

  # Batch compare from a pairs file
  python archetype_mapper.py batch-compare pairs.txt -w 4 -o comparisons.csv
        """
    )
    sub = parser.add_subparsers(dest="command")
    
    # Census
    p_cen = sub.add_parser("census", help="Single structure baseline")
    p_cen.add_argument("cif", help="CIF file path")
    p_cen.add_argument("-c", "--chain", default="A")
    
    # Compare
    p_cmp = sub.add_parser("compare", help="Two-structure transition")
    p_cmp.add_argument("ref", help="Reference CIF")
    p_cmp.add_argument("pert", help="Perturbed CIF")
    p_cmp.add_argument("-c", "--chain", default="A")
    
    # Batch census
    p_bat = sub.add_parser("batch", help="Batch census of directory (threaded)")
    p_bat.add_argument("directory", help="Directory of CIF files")
    p_bat.add_argument("-c", "--chain", default="A")
    p_bat.add_argument("-o", "--output", default="census_results.csv",
                       help="Output CSV (supports resume)")
    p_bat.add_argument("-w", "--workers", type=int, default=1,
                       help="Number of parallel threads (default: 1, recommended: 4-16)")
    p_bat.add_argument("-p", "--pattern", default="*.cif",
                       help="File glob pattern (default: *.cif)")
    p_bat.add_argument("--no-resume", action="store_true",
                       help="Don't resume from existing CSV, start fresh")
    
    # Batch compare
    p_bcmp = sub.add_parser("batch-compare", help="Batch compare from pairs file")
    p_bcmp.add_argument("pairs_file", help="Tab-separated: ref_path\\tpert_path\\tlabel")
    p_bcmp.add_argument("-c", "--chain", default="A")
    p_bcmp.add_argument("-o", "--output", default="compare_results.csv")
    p_bcmp.add_argument("-w", "--workers", type=int, default=1)
    
    args = parser.parse_args()
    
    if args.command == "census":
        r = census(args.cif, chain=args.chain)
        print(json.dumps({k: v for k, v in r.items() if not k.startswith("_")}, indent=2))
    
    elif args.command == "compare":
        r = compare(args.ref, args.pert, chain=args.chain)
        arch = classify_archetype(r)
        print(json.dumps({k: v for k, v in r.items() if not k.startswith("_")}, indent=2))
        print(f"\n  ARCHETYPE: {arch}")
    
    elif args.command == "batch":
        results = batch_census(
            args.directory,
            chain=args.chain,
            pattern=args.pattern,
            workers=args.workers,
            output_csv=args.output,
            resume=not args.no_resume,
        )
    
    elif args.command == "batch-compare":
        pairs = []
        with open(args.pairs_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 3:
                    pairs.append((parts[0], parts[1], parts[2]))
                elif len(parts) == 2:
                    label = f"{Path(parts[0]).stem}_vs_{Path(parts[1]).stem}"
                    pairs.append((parts[0], parts[1], label))
        
        if pairs:
            batch_compare(pairs, chain=args.chain, workers=args.workers,
                         output_csv=args.output)
        else:
            print("  No valid pairs found in file")
    
    else:
        parser.print_help()

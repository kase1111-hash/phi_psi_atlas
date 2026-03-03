#!/usr/bin/env python3
"""
kinase_dfg_survey.py
====================
Downloads and analyzes kinase structures across DFG-in and DFG-out states.
For each kinase: apo + DFG-in inhibitor + DFG-out inhibitor.

Usage:
    python kinase_dfg_survey.py

Outputs:
    kinase_dfg_report.txt
    kinase_dfg_results.csv
"""

import os
import sys
import subprocess
import csv
import json
import statistics
import urllib.request

BASE_DIR = r"G:\Mechanical Taxonomy"
CACHE_DIR = os.path.join(BASE_DIR, "pdb_cache")
REPORT_FILE = os.path.join(BASE_DIR, "kinase_dfg_report.txt")
CSV_FILE = os.path.join(BASE_DIR, "kinase_dfg_results.csv")

# ─────────────────────────────────────────────────────────
# KINASE DATABASE
# Each entry: (PDB_ID, chain, condition, DFG_state)
# DFG states: "apo", "DFG-in", "DFG-out"
# ─────────────────────────────────────────────────────────

KINASES = {
    # Already tested (will re-run for consistency)
    "CDK2": [
        ("1HCL", "A", "Apo (inactive)", "apo"),
        ("1AQ1", "A", "Staurosporine (potent)", "DFG-in"),
        ("2VTP", "A", "LZ9 (inactive-stabilizer)", "DFG-out"),
    ],
    "Abl": [
        ("2G1T", "A", "Apo (DFG-out)", "apo"),
        ("2HYY", "A", "Dasatinib", "DFG-in"),
        ("1IEP", "A", "Imatinib", "DFG-out"),
    ],
    "p38 MAPK": [
        ("1P38", "A", "Apo", "apo"),
        ("1KV2", "A", "SB203580", "DFG-in"),
        ("3GI3", "A", "BIRB796", "DFG-out"),
    ],
    "EGFR": [
        ("2GS7", "A", "Apo (+ ANP analog)", "apo"),
        ("1M17", "A", "Erlotinib", "DFG-in"),
        ("2RF9", "A", "Lapatinib", "DFG-out"),
    ],
    # NEW kinases
    "B-Raf": [
        ("4E26", "A", "Apo-like", "apo"),
        ("3OG7", "A", "Vemurafenib", "DFG-in"),
        ("1UWH", "A", "Sorafenib (BAY43-9006)", "DFG-out"),
    ],
    "c-Met": [
        ("3DKC", "A", "Apo (unliganded)", "apo"),
        ("3DKG", "A", "Crizotinib-analog", "DFG-in"),
        ("3LQ8", "A", "ARQ197 (tivantinib)", "DFG-out"),
    ],
    "JAK2": [
        ("4FVQ", "A", "Apo-like", "apo"),
        ("3FUP", "A", "CMP6 (type I)", "DFG-in"),
        ("3E62", "A", "JAK inhibitor (type II)", "DFG-out"),
    ],
    "Aurora A": [
        ("4J8M", "A", "Apo", "apo"),
        ("3E5A", "A", "VX-680 (tozasertib)", "DFG-in"),
        ("3UOH", "A", "CD532 (DFG-out)", "DFG-out"),
    ],
    "PKA (cAMP-dep)": [
        ("1J3H", "A", "Apo (open)", "apo"),
        ("1ATP", "E", "PKI + ATP (active)", "DFG-in"),
        ("4WB5", "A", "Type II inhibitor", "DFG-out"),
    ],
    "c-Src": [
        ("1YOJ", "A", "Apo (inactive)", "apo"),
        ("3G5D", "A", "Dasatinib", "DFG-in"),
        ("2OIQ", "A", "Imatinib-like", "DFG-out"),
    ],
    "VEGFR2 (KDR)": [
        ("4AGD", "A", "Apo-like", "apo"),
        ("3VHE", "A", "Type I inhibitor", "DFG-in"),
        ("3VHK", "A", "Sorafenib-like (type II)", "DFG-out"),
    ],
    "FGFR1": [
        ("3GQI", "A", "Apo-like (unliganded)", "apo"),
        ("3C4F", "A", "SU5402 derivative", "DFG-in"),
        ("3RHX", "A", "Ponatinib", "DFG-out"),
    ],
    "IGF1R": [
        ("1P4O", "A", "Apo (inactive)", "apo"),
        ("1K3A", "A", "AMP-PNP (active)", "DFG-in"),
        ("3LVP", "A", "OSI-906 (linsitinib)", "DFG-out"),
    ],
    "PDK1": [
        ("3ORX", "A", "Apo-like", "apo"),
        ("3NAX", "A", "Compound 7 (type I)", "DFG-in"),
        ("3ORZ", "A", "RS1 (allosteric/type II)", "DFG-out"),
    ],
    "CDK6": [
        ("1BI7", "A", "Apo (with V-cyclin)", "apo"),
        ("2EUF", "A", "Fisetin (flavonoid)", "DFG-in"),
        ("5L2T", "A", "Palbociclib", "DFG-out"),
    ],
    "Nek2": [
        ("2JAV", "A", "Apo", "apo"),
        ("2W5A", "A", "Aminopyrazine inh.", "DFG-in"),
        ("5M53", "A", "JH295 (type II)", "DFG-out"),
    ],
}

# ─────────────────────────────────────────────────────────
# STEP 1: Download missing CIF files
# ─────────────────────────────────────────────────────────

def download_cif(pdb_id):
    """Download CIF file from RCSB if not already cached."""
    cif_path = os.path.join(CACHE_DIR, f"{pdb_id}.cif")
    if os.path.exists(cif_path):
        return True
    
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    print(f"  Downloading {pdb_id}...", end=" ", flush=True)
    try:
        urllib.request.urlretrieve(url, cif_path)
        print("OK")
        return True
    except Exception as e:
        print(f"FAILED ({e})")
        return False


# ─────────────────────────────────────────────────────────
# STEP 2: Run archetype_mapper on each structure
# ─────────────────────────────────────────────────────────

def run_mapper(pdb_id, chain):
    """Run archetype_mapper.py and return parsed JSON result."""
    cif_path = os.path.join(CACHE_DIR, f"{pdb_id}.cif")
    if not os.path.exists(cif_path):
        return None
    
    cmd = [
        sys.executable, 
        os.path.join(BASE_DIR, "archetype_mapper.py"),
        "census", cif_path, "--chain", chain
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            return None
        
        # Parse JSON from stdout
        output = result.stdout.strip()
        data = json.loads(output)
        
        if "error" in data:
            # Try alternative chains
            for alt_chain in ["B", "C", "H", "L"]:
                if alt_chain == chain:
                    continue
                cmd2 = [
                    sys.executable,
                    os.path.join(BASE_DIR, "archetype_mapper.py"),
                    "census", cif_path, "--chain", alt_chain
                ]
                r2 = subprocess.run(cmd2, capture_output=True, text=True, timeout=120)
                if r2.returncode == 0:
                    d2 = json.loads(r2.stdout.strip())
                    if "error" not in d2 and d2.get("n_valid", 0) > 100:
                        data = d2
                        data["_actual_chain"] = alt_chain
                        break
        
        return data
    except Exception as e:
        print(f"    Error running mapper on {pdb_id}: {e}")
        return None


# ─────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────

def main():
    report = []
    csv_rows = []
    
    def log(msg):
        print(msg)
        report.append(msg)
    
    log("=" * 70)
    log("  KINASE DFG-IN vs DFG-OUT GEOMETRIC SURVEY")
    log(f"  {len(KINASES)} kinases, {sum(len(v) for v in KINASES.values())} structures")
    log("=" * 70)
    
    # Step 1: Download all CIF files
    log("\nSTEP 1: Downloading structures...")
    all_pdbs = set()
    for kinase, entries in KINASES.items():
        for pdb_id, chain, cond, dfg in entries:
            all_pdbs.add(pdb_id)
    
    downloaded = 0
    failed = 0
    for pdb_id in sorted(all_pdbs):
        if download_cif(pdb_id):
            downloaded += 1
        else:
            failed += 1
    
    log(f"  {downloaded} available, {failed} failed")
    
    # Step 2: Run mapper on all structures
    log("\nSTEP 2: Computing Gini for all structures...")
    
    results = {}  # kinase -> [(pdb, chain, cond, dfg, gini, kmax, n_res)]
    
    for kinase, entries in KINASES.items():
        log(f"\n  {kinase}:")
        kinase_results = []
        
        for pdb_id, chain, cond, dfg in entries:
            data = run_mapper(pdb_id, chain)
            
            if data and "error" not in data and data.get("n_valid", 0) > 100:
                gini = data["gini"]
                kmax = data["kappa2_max"]
                n_res = data["n_valid"]
                actual_chain = data.get("_actual_chain", chain)
                log(f"    {pdb_id} chain {actual_chain}: Gini = {gini:.4f} ({cond}, {dfg})")
                kinase_results.append((pdb_id, actual_chain, cond, dfg, gini, kmax, n_res))
                
                csv_rows.append({
                    "kinase": kinase,
                    "pdb_id": pdb_id,
                    "chain": actual_chain,
                    "condition": cond,
                    "dfg_state": dfg,
                    "gini": gini,
                    "kappa2_max": kmax,
                    "n_residues": n_res,
                })
            else:
                err = data.get("error", "unknown") if data else "mapper failed"
                log(f"    {pdb_id} chain {chain}: FAILED ({err})")
        
        results[kinase] = kinase_results
    
    # Step 3: Compute delta-Gini relative to apo for each kinase
    log(f"\n\n{'='*70}")
    log(f"  RESULTS: DELTA-GINI BY DFG STATE")
    log(f"{'='*70}")
    
    all_dfg_in = []
    all_dfg_out = []
    
    log(f"\n  {'Kinase':>12s} {'PDB':>5s} {'Condition':>25s} {'DFG':>8s} {'Gini':>7s} {'dGini':>8s}")
    log(f"  {'─'*12} {'─'*5} {'─'*25} {'─'*8} {'─'*7} {'─'*8}")
    
    for kinase, entries in results.items():
        if not entries:
            continue
        
        # Find apo
        apo_entries = [e for e in entries if e[3] == "apo"]
        if not apo_entries:
            log(f"  {kinase:>12s}: No apo structure, skipping")
            continue
        
        apo_gini = apo_entries[0][4]
        
        for pdb_id, chain, cond, dfg, gini, kmax, n_res in entries:
            dg = gini - apo_gini
            marker = ""
            
            if dfg == "DFG-in":
                all_dfg_in.append(dg)
                marker = " ←in"
            elif dfg == "DFG-out":
                all_dfg_out.append(dg)
                marker = " ←out"
            
            log(f"  {kinase:>12s} {pdb_id:>5s} {cond:>25s} {dfg:>8s} {gini:>7.4f} {dg:>+8.4f}{marker}")
    
    # Step 4: Aggregate statistics
    log(f"\n\n{'='*70}")
    log(f"  AGGREGATE STATISTICS")
    log(f"{'='*70}")
    
    if all_dfg_in:
        in_pos = len([d for d in all_dfg_in if d > 0])
        log(f"\n  DFG-in (active-state) inhibitors:")
        log(f"    n = {len(all_dfg_in)}")
        log(f"    Mean delta-Gini: {statistics.mean(all_dfg_in):+.4f}")
        log(f"    Median: {statistics.median(all_dfg_in):+.4f}")
        if len(all_dfg_in) > 1:
            log(f"    SD: {statistics.stdev(all_dfg_in):.4f}")
        log(f"    Stiffen: {in_pos}/{len(all_dfg_in)} ({100*in_pos/len(all_dfg_in):.0f}%)")
    
    if all_dfg_out:
        out_neg = len([d for d in all_dfg_out if d < 0])
        log(f"\n  DFG-out (inactive-state) inhibitors:")
        log(f"    n = {len(all_dfg_out)}")
        log(f"    Mean delta-Gini: {statistics.mean(all_dfg_out):+.4f}")
        log(f"    Median: {statistics.median(all_dfg_out):+.4f}")
        if len(all_dfg_out) > 1:
            log(f"    SD: {statistics.stdev(all_dfg_out):.4f}")
        log(f"    Loosen: {out_neg}/{len(all_dfg_out)} ({100*out_neg/len(all_dfg_out):.0f}%)")
    
    if all_dfg_in and all_dfg_out:
        from scipy.stats import mannwhitneyu
        try:
            u, p = mannwhitneyu(all_dfg_in, all_dfg_out, alternative='greater')
            log(f"\n  Mann-Whitney (DFG-in > DFG-out): U={u:.0f}, p={p:.4f}")
        except:
            log(f"\n  Mann-Whitney: could not compute (install scipy)")
        
        # Effect size
        n1, n2 = len(all_dfg_in), len(all_dfg_out)
        mean_diff = statistics.mean(all_dfg_in) - statistics.mean(all_dfg_out)
        log(f"  Mean difference (DFG-in minus DFG-out): {mean_diff:+.4f}")
    
    # Write report
    with open(REPORT_FILE, 'w') as f:
        f.write('\n'.join(report))
    log(f"\n\nReport: {REPORT_FILE}")
    
    # Write CSV
    if csv_rows:
        with open(CSV_FILE, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                "kinase", "pdb_id", "chain", "condition", 
                "dfg_state", "gini", "kappa2_max", "n_residues"
            ])
            writer.writeheader()
            writer.writerows(csv_rows)
        log(f"CSV: {CSV_FILE}")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Download PDB files for TRUE POSITIVE validation pairs.

These confirm that annotated homologs at the SAME scores as false positives
DO have genuine structural similarity (TM > 0.5).

Usage:
    python download_tp_pairs.py
"""

import json
import shutil
import time
import urllib.request
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
PDB_DIR = SCRIPT_DIR / "dali_pdb"
PAIR_DIR = SCRIPT_DIR / "tp_validation_pairs"

# Strategic picks:
# 1. High confidence (>0.85): LDH-A → LDH-B (0.926) — Rossmann fold paralogs
# 2. Overlap zone: Hb-β → Hb-G1 (0.808) — globin, same score range as FPs  
# 3. Overlap zone: LDH-A → MDH2 (0.815) — Rossmann superfamily, different enzyme
# 4. Overlap zone: 14-3-3ε → YWHAB (0.799) — same fold family
# 5. Low TP: Hb-α → Myoglobin (0.724) — classic remote homologs
# 6. Overlap zone: CA2 → CA6 (0.798) — same fold, EXACT same score as highest FP

TP_PAIRS = [
    (11, "P00338", "P07195", "LDH-A",   "LDH-B",  0.926, "Rossmann paralogs — high confidence control"),
    (12, "P68871", "P69891", "Hb-beta",  "Hb-G1",  0.808, "Globins — overlap zone, same scores as FPs"),
    (13, "P00338", "P40926", "LDH-A",   "MDH2",   0.815, "Rossmann superfamily — different enzyme"),
    (14, "P62258", "P31946", "14-3-3e",  "YWHAB",  0.799, "14-3-3 family — overlap zone"),
    (15, "P69905", "P02144", "Hb-alpha", "Myoglob", 0.724, "Classic remote homologs — low score"),
    (16, "P00918", "P23280", "CA2",      "CA6",    0.798, "Carbonic anhydrases — matches FP top score"),
]


def download_pdb(uid, gene_name):
    """Download PDB file from AlphaFold."""
    outpath = PDB_DIR / f"{gene_name}_{uid}.pdb"

    if outpath.exists() and outpath.stat().st_size > 5000:
        with open(outpath) as f:
            first = f.read(6)
        if first.startswith(("HEADER", "ATOM  ", "REMARK", "TITLE ")):
            print(f"    {gene_name:>10s} ({uid}): cached ({outpath.stat().st_size:,} bytes)")
            return outpath

    headers = {'User-Agent': 'Mozilla/5.0 PhiPsiAtlas/1.0'}

    for version in ['v6', 'v4']:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_{version}.pdb"
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read()
                text = data[:100].decode('utf-8', errors='ignore')
                if any(text.startswith(tag) for tag in ['HEADER', 'ATOM', 'MODEL', 'REMARK', 'TITLE']):
                    with open(outpath, 'wb') as f:
                        f.write(data)
                    print(f"    {gene_name:>10s} ({uid}): downloaded {version} ({len(data):,} bytes)")
                    return outpath
        except:
            pass

    # API fallback
    api_url = f"https://alphafold.com/api/prediction/{uid}"
    try:
        req = urllib.request.Request(api_url, headers=headers)
        with urllib.request.urlopen(req, timeout=30) as resp:
            entries = json.loads(resp.read().decode('utf-8'))
            if isinstance(entries, list) and entries:
                pdb_url = entries[0].get('pdbUrl', '')
                if pdb_url:
                    req2 = urllib.request.Request(pdb_url, headers=headers)
                    with urllib.request.urlopen(req2, timeout=30) as resp2:
                        data = resp2.read()
                        text = data[:100].decode('utf-8', errors='ignore')
                        if any(text.startswith(tag) for tag in ['HEADER', 'ATOM', 'MODEL', 'REMARK', 'TITLE']):
                            with open(outpath, 'wb') as f:
                                f.write(data)
                            print(f"    {gene_name:>10s} ({uid}): via API ({len(data):,} bytes)")
                            return outpath
    except Exception as e:
        print(f"    {gene_name:>10s} ({uid}): FAILED: {e}")

    return None


def main():
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    PAIR_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\n  Downloading TRUE POSITIVE validation structures...\n")

    uids = {}
    for num, q_uid, h_uid, q_gene, h_gene, geom, note in TP_PAIRS:
        uids[q_uid] = q_gene
        uids[h_uid] = h_gene

    for uid, gene in sorted(uids.items(), key=lambda x: x[1]):
        download_pdb(uid, gene)
        time.sleep(0.3)

    print(f"\n  Organizing into pair folders...\n")

    for num, q_uid, h_uid, q_gene, h_gene, geom, note in TP_PAIRS:
        q_pdb = PDB_DIR / f"{q_gene}_{q_uid}.pdb"
        h_pdb = PDB_DIR / f"{h_gene}_{h_uid}.pdb"

        if not q_pdb.exists() or not h_pdb.exists():
            print(f"  Pair {num}: MISSING files")
            continue

        pair_folder = PAIR_DIR / f"pair{num:02d}_{q_gene}_vs_{h_gene}"
        pair_folder.mkdir(parents=True, exist_ok=True)

        dst_q = pair_folder / f"structure1_{q_gene}_{q_uid}.pdb"
        dst_h = pair_folder / f"structure2_{h_gene}_{h_uid}.pdb"

        shutil.copy2(q_pdb, dst_q)
        shutil.copy2(h_pdb, dst_h)

        status = "OK" if (dst_q.stat().st_size > 5000 and dst_h.stat().st_size > 5000) else "CHECK"
        print(f"  Pair {num:2d}: {q_gene:>10s} vs {h_gene:>10s}  [geom={geom:.3f}] {status}")
        print(f"          {note}")
        print(f"          → {pair_folder.name}/")
        print()

    print(f"  TP pair folders in: {PAIR_DIR}/")
    print()
    print(f"  SUBMIT TO: https://www.rcsb.org/alignment")
    print(f"  Left:  Upload structure1_*.pdb, Chain ID: A")
    print(f"  Right: Upload structure2_*.pdb, Chain ID: A")
    print(f"  Algorithm: TM-align")
    print()
    print(f"  EXPECTED RESULTS:")
    print(f"  Pair 11: LDH-A  vs LDH-B   (0.926) → expect TM > 0.8 (close paralogs)")
    print(f"  Pair 12: Hb-β   vs Hb-G1   (0.808) → expect TM > 0.7 (globin family)")
    print(f"  Pair 13: LDH-A  vs MDH2    (0.815) → expect TM > 0.5 (Rossmann superfam)")
    print(f"  Pair 14: 14-3-3 vs YWHAB   (0.799) → expect TM > 0.8 (same family)")
    print(f"  Pair 15: Hb-α   vs Myoglob (0.724) → expect TM > 0.5 (remote globins)")
    print(f"  Pair 16: CA2    vs CA6     (0.798) → expect TM > 0.7 (same fold)")
    print()
    print(f"  KEY COMPARISON: Pair 16 (CA2→CA6, geom=0.798, TP)")
    print(f"              vs  Pair 1  (CA2→KRTAP, geom=0.804, FP with TM=0.22)")
    print(f"  Same query, nearly identical geometric scores, one is real, one isn't.")


if __name__ == "__main__":
    main()

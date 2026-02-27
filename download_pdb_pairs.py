#!/usr/bin/env python3
"""
Download PDB files directly from AlphaFold and organize into pair folders.
No conversion needed — AlphaFold provides native PDB format.

Usage:
    python download_pdb_pairs.py
"""

import json
import shutil
import sys
import time
import urllib.request
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
PDB_DIR = SCRIPT_DIR / "dali_pdb"
PAIR_DIR = SCRIPT_DIR / "dali_pairs"

PAIRS = [
    (1,  "P00918", "P59991", "CA2",     "KRTAP12-2",  0.804),
    (2,  "P60174", "P35270", "TPI",     "SPR",        0.798),
    (3,  "P62258", "O43808", "14-3-3e", "SLC25A17",   0.783),
    (4,  "P11166", "P09914", "GLUT1",   "IFIT1",      0.783),
    (5,  "P02751", "P21333", "FN1",     "FLNA",       0.779),
    (6,  "P01857", "P32969", "IgG1",    "RPL9",       0.778),
    (7,  "P38646", "O15344", "GRP75",   "MID1",       0.767),
    (8,  "P07900", "O14829", "HSP90a",  "PPEF1",      0.765),
    (9,  "P04626", "O95302", "HER2",    "FKBP9",      0.757),
    (10, "P21860", "P28799", "ErbB3",   "GRN",        0.757),
]


def download_pdb(uid, gene_name):
    """Download PDB file from AlphaFold. Returns path or None."""
    outpath = PDB_DIR / f"{gene_name}_{uid}.pdb"
    
    if outpath.exists() and outpath.stat().st_size > 5000:
        # Verify it's actually a PDB file
        with open(outpath) as f:
            first = f.read(6)
        if first.startswith(("HEADER", "ATOM  ", "REMARK", "TITLE ")):
            print(f"    {gene_name:>10s} ({uid}): cached ({outpath.stat().st_size:,} bytes)")
            return outpath

    headers = {'User-Agent': 'Mozilla/5.0 PhiPsiAtlas/1.0'}

    # Try direct PDB URLs (v6 first, then v4)
    for version in ['v6', 'v4', 'v3']:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_{version}.pdb"
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read()
                # Verify it's PDB format (starts with HEADER, ATOM, MODEL, etc.)
                text = data[:100].decode('utf-8', errors='ignore')
                if any(text.startswith(tag) for tag in ['HEADER', 'ATOM', 'MODEL', 'REMARK', 'TITLE']):
                    with open(outpath, 'wb') as f:
                        f.write(data)
                    print(f"    {gene_name:>10s} ({uid}): downloaded {version} ({len(data):,} bytes)")
                    return outpath
        except urllib.error.HTTPError:
            pass
        except Exception as e:
            print(f"    {gene_name:>10s} ({uid}): {version} failed: {e}")

    # Fallback: query API for pdbUrl
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
                            print(f"    {gene_name:>10s} ({uid}): downloaded via API ({len(data):,} bytes)")
                            return outpath
                        else:
                            print(f"    {gene_name:>10s} ({uid}): API pdbUrl returned non-PDB data")
    except Exception as e:
        print(f"    {gene_name:>10s} ({uid}): API fallback failed: {e}")

    print(f"    {gene_name:>10s} ({uid}): FAILED all methods")
    return None


def main():
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    PAIR_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\n  Downloading PDB files from AlphaFold...\n")

    # Collect unique downloads needed
    uids = {}
    for num, q_uid, h_uid, q_gene, h_gene, geom in PAIRS:
        uids[q_uid] = q_gene
        uids[h_uid] = h_gene

    print(f"  {len(uids)} unique structures needed\n")

    # Download all
    success = 0
    for uid, gene in sorted(uids.items(), key=lambda x: x[1]):
        path = download_pdb(uid, gene)
        if path:
            success += 1
        time.sleep(0.3)

    print(f"\n  Downloaded: {success}/{len(uids)}")

    # Organize into pair folders
    print(f"\n  Organizing into pair folders...\n")

    for num, q_uid, h_uid, q_gene, h_gene, geom in PAIRS:
        q_pdb = PDB_DIR / f"{q_gene}_{q_uid}.pdb"
        h_pdb = PDB_DIR / f"{h_gene}_{h_uid}.pdb"

        if not q_pdb.exists():
            print(f"  Pair {num}: MISSING {q_gene} ({q_uid})")
            continue
        if not h_pdb.exists():
            print(f"  Pair {num}: MISSING {h_gene} ({h_uid})")
            continue

        pair_folder = PAIR_DIR / f"pair{num:02d}_{q_gene}_vs_{h_gene}"
        pair_folder.mkdir(parents=True, exist_ok=True)

        dst_q = pair_folder / f"structure1_{q_gene}_{q_uid}.pdb"
        dst_h = pair_folder / f"structure2_{h_gene}_{h_uid}.pdb"

        shutil.copy2(q_pdb, dst_q)
        shutil.copy2(h_pdb, dst_h)

        # Verify files
        q_ok = dst_q.stat().st_size > 5000
        h_ok = dst_h.stat().st_size > 5000

        status = "OK" if (q_ok and h_ok) else "PROBLEM"
        print(f"  Pair {num:2d}: {q_gene:>8s} vs {h_gene:>10s}  [{status}]"
              f"  {dst_q.stat().st_size:>8,} + {dst_h.stat().st_size:>8,} bytes")

    print(f"\n  Pair folders in: {PAIR_DIR}/")
    print(f"  Each folder has:")
    print(f"    structure1_<GENE>_<UID>.pdb  ← chain A")
    print(f"    structure2_<GENE>_<UID>.pdb  ← chain A")
    print()
    print(f"  SUBMIT TO: https://www.rcsb.org/alignment")
    print(f"  Left:  Upload structure1_*.pdb, Chain ID: A")
    print(f"  Right: Upload structure2_*.pdb, Chain ID: A")
    print(f"  Algorithm: TM-align")
    print(f"  Click Submit → record TM-score")
    print()
    print(f"  Then run: python structural_validate.py results")


if __name__ == "__main__":
    main()

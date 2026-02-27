#!/usr/bin/env python3
"""
Organize DALI pairs into clearly named folders and print a submission guide.

Usage:
    python dali_pair_guide.py           # Print guide + create pair folders
    python dali_pair_guide.py --copy    # Also copy CIF files into pair folders
"""

import shutil
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
CIF_DIR = SCRIPT_DIR / "dali_structures"
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


def find_cif(uid):
    """Find CIF file for a UniProt ID (checks multiple naming patterns)."""
    patterns = [
        CIF_DIR / f"AF-{uid}-F1.cif",
        CIF_DIR / f"AF-{uid}-F1-model_v6.cif",
        CIF_DIR / f"AF-{uid}-F1-model_v4.cif",
        CIF_DIR / f"{uid}.cif",
    ]
    for p in patterns:
        if p.exists() and p.stat().st_size > 1000:
            return p
    return None


def main():
    do_copy = "--copy" in sys.argv

    print()
    print("  ╔══════════════════════════════════════════════════════════════════╗")
    print("  ║              DALI PAIR SUBMISSION GUIDE                         ║")
    print("  ║  http://ekhidna2.biocenter.helsinki.fi/dali/                    ║")
    print("  ║  Select: 'Pairwise structure comparison'                        ║")
    print("  ╚══════════════════════════════════════════════════════════════════╝")
    print()

    if do_copy:
        PAIR_DIR.mkdir(parents=True, exist_ok=True)

    for num, q_uid, h_uid, q_gene, h_gene, geom in PAIRS:
        q_cif = find_cif(q_uid)
        h_cif = find_cif(h_uid)

        q_status = f"OK ({q_cif.name})" if q_cif else "MISSING"
        h_status = f"OK ({h_cif.name})" if h_cif else "MISSING"

        print(f"  ── Pair {num:2d} ──────────────────────────────────────────────")
        print(f"     Structure 1 (query):  {q_gene:>8s}  ({q_uid})  [{q_status}]")
        print(f"     Structure 2 (hit):    {h_gene:>8s}  ({h_uid})  [{h_status}]")
        print(f"     Geometric score: {geom:.3f}")
        print()

        if do_copy and q_cif and h_cif:
            pair_folder = PAIR_DIR / f"pair{num:02d}_{q_gene}_vs_{h_gene}"
            pair_folder.mkdir(parents=True, exist_ok=True)

            dst_q = pair_folder / f"structure1_{q_gene}_{q_uid}.cif"
            dst_h = pair_folder / f"structure2_{h_gene}_{h_uid}.cif"

            shutil.copy2(q_cif, dst_q)
            shutil.copy2(h_cif, dst_h)
            print(f"     Copied to: {pair_folder.name}/")
            print(f"       → structure1_{q_gene}_{q_uid}.cif")
            print(f"       → structure2_{h_gene}_{h_uid}.cif")
            print()

    if do_copy:
        print(f"  Pair folders created in: {PAIR_DIR}/")
        print(f"  Each folder has structure1_*.cif and structure2_*.cif")
        print()

    print(f"  WORKFLOW:")
    print(f"  1. For each pair, upload Structure 1 and Structure 2 to DALI")
    print(f"  2. Wait for results (30s - 3min per pair)")
    print(f"  3. Record the Z-score")
    print(f"  4. Run: python dali_batch.py manual")
    print()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Fetch AlphaFold structures for benchmark query panel proteins.

Downloads CIF files to alphafold/ cache, then runs analyze_barcodes.py
to generate barcodes for the new proteins.

Usage:
    python fetch_benchmark_proteins.py          # Fetch missing query panel proteins
    python fetch_benchmark_proteins.py --all    # Fetch all, even if already barcoded
    python fetch_benchmark_proteins.py --check  # Just check which are missing
"""

import json
import os
import sys
import urllib.request
import urllib.error
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent

# AlphaFold cache directories — try both
CACHE_DIRS = [
    SCRIPT_DIR / "alphafold_cache",
    SCRIPT_DIR / "alphafold",
]

BARCODE_DIR = SCRIPT_DIR / "results" / "barcodes"

# Benchmark query panel — 20 proteins spanning structural space
QUERY_PANEL = {
    # All-alpha
    "P69905": "Hemoglobin alpha (all-alpha, globin)",
    "P68871": "Hemoglobin beta (all-alpha, globin)",
    "P02144": "Myoglobin (all-alpha, globin)",
    "P00918": "Carbonic anhydrase 2 (alpha/beta)",
    # RTK
    "P00533": "EGFR (multi-domain, RTK)",
    "P04626": "ErbB2/HER2 (multi-domain, RTK)",
    "P21860": "ErbB3 (multi-domain, RTK)",
    # Beta-rich
    "P01857": "IgG1 Fc (all-beta, Ig fold)",
    "P01834": "Ig kappa chain C (all-beta, Ig fold)",
    "P02751": "Fibronectin (beta-rich, FN3 repeats)",
    # Alpha/beta
    "P00338": "L-lactate dehydrogenase A (Rossmann)",
    "P04406": "GAPDH (alpha/beta, Rossmann)",
    "P60174": "Triosephosphate isomerase (TIM barrel)",
    # Small/mixed
    "P01308": "Insulin (small, alpha+beta)",
    "P62258": "14-3-3 epsilon (all-alpha)",
    # Membrane
    "P11166": "GLUT1 (membrane transporter)",
    # Disordered
    "P06748": "Nucleophosmin (mixed)",
    "P10636": "Tau (largely disordered)",
    # Chaperones
    "P07900": "HSP90-alpha (alpha-rich)",
    "P38646": "Mortalin/GRP75 (HSP70 fold)",
}


def find_cache_dir():
    """Find the active AlphaFold cache directory."""
    for d in CACHE_DIRS:
        if d.exists() and any(d.iterdir()):
            return d
    # Default to first
    CACHE_DIRS[0].mkdir(parents=True, exist_ok=True)
    return CACHE_DIRS[0]


def has_barcode(uid):
    """Check if a barcode already exists."""
    bc_file = BARCODE_DIR / f"{uid}.json"
    return bc_file.exists()


def has_structure(uid, cache_dir):
    """Check if structure file exists in cache."""
    # Check all possible naming patterns
    patterns = [
        f"AF-{uid}-F1-model_v*.cif",
        f"AF-{uid}-F1-model_v*.cif.gz",
        f"AF-{uid}-F1*.cif",
        f"*{uid}*.cif",
        f"*{uid}*.pdb",
    ]
    for pattern in patterns:
        matches = list(cache_dir.glob(pattern))
        if matches:
            return True
    return False


def fetch_alphafold_structure(uid, cache_dir):
    """Download AlphaFold structure for a UniProt ID.
    
    Tries multiple URL formats since the API has changed over time:
    1. New alphafold.com API (v6/v4)
    2. Legacy alphafold.ebi.ac.uk (v4/v3/v2)
    3. Direct 3D-Beacons/EBI FTP fallback
    """
    # Check if already exists
    for version in [6, 4, 3, 2]:
        outfile = cache_dir / f"AF-{uid}-F1-model_v{version}.cif"
        if outfile.exists():
            return str(outfile)
    
    # URL patterns to try (newest first)
    url_patterns = [
        # New alphafold.com format
        "https://alphafold.com/files/AF-{uid}-F1-model_v4.cif",
        # EBI hosted (still works for many)
        "https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.cif",
        # Older versions
        "https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v3.cif",
        "https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v2.cif",
        # 3D-Beacons / ModelArchive fallback
        "https://www.ebi.ac.uk/pdbe/static/entry/AF-{uid}-F1-model_v4.cif",
    ]
    
    for url_template in url_patterns:
        url = url_template.format(uid=uid)
        # Extract version from URL for filename
        version = "4"
        for v in ["_v6", "_v4", "_v3", "_v2"]:
            if v in url:
                version = v[2:]
                break
        
        outfile = cache_dir / f"AF-{uid}-F1-model_v{version}.cif"
        
        try:
            req = urllib.request.Request(url, headers={
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) PhiPsiAtlas/1.0',
                'Accept': '*/*',
            })
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read()
                if len(data) > 100:  # Sanity check - not an error page
                    with open(outfile, 'wb') as f:
                        f.write(data)
                    return str(outfile)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                continue
            # Try next URL pattern
            continue
        except Exception:
            continue
    
    return None


def main():
    args = sys.argv[1:]
    check_only = '--check' in args
    fetch_all = '--all' in args
    
    cache_dir = find_cache_dir()
    print(f"  AlphaFold cache: {cache_dir}")
    print(f"  Barcode dir: {BARCODE_DIR}")
    print(f"  Query panel: {len(QUERY_PANEL)} proteins")
    print()
    
    missing_structure = []
    missing_barcode = []
    have_both = []
    
    for uid, desc in QUERY_PANEL.items():
        has_bc = has_barcode(uid)
        has_str = has_structure(uid, cache_dir)
        
        status = ""
        if has_bc and has_str:
            status = "✓ ready"
            have_both.append(uid)
        elif has_str and not has_bc:
            status = "⚠ structure exists, needs barcoding"
            missing_barcode.append(uid)
        else:
            status = "✗ needs fetch"
            missing_structure.append(uid)
        
        if fetch_all or not (has_bc and has_str):
            print(f"  {uid:>10s} {status:>35s}  {desc}")
    
    print(f"\n  Summary: {len(have_both)} ready, "
          f"{len(missing_barcode)} need barcoding, "
          f"{len(missing_structure)} need fetching")
    
    if check_only:
        return
    
    # Fetch missing structures
    if missing_structure:
        print(f"\n  Fetching {len(missing_structure)} structures from AlphaFold...")
        fetched = 0
        failed = []
        for uid in missing_structure:
            desc = QUERY_PANEL[uid]
            result = fetch_alphafold_structure(uid, cache_dir)
            if result:
                fetched += 1
                print(f"    ✓ {uid} ({desc[:30]})")
                missing_barcode.append(uid)  # Now needs barcoding
            else:
                failed.append(uid)
                print(f"    ✗ {uid} ({desc[:30]}) — NOT FOUND")
        
        print(f"\n  Fetched: {fetched}/{len(missing_structure)}")
        if failed:
            print(f"  Failed: {', '.join(failed)}")
    
    # Report what needs barcoding
    if missing_barcode:
        print(f"\n  {len(missing_barcode)} proteins need barcoding.")
        print(f"  Run: python analyze_barcodes.py")
        print(f"  This will process all structures in {cache_dir} that don't have barcodes yet.")
    else:
        print(f"\n  All query panel proteins are barcoded and ready!")
        print(f"  Run: python benchmark_homology.py run")


if __name__ == "__main__":
    main()

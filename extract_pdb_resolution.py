#!/usr/bin/env python3
"""
extract_pdb_resolution.py
Extracts resolution from PDB CIF/mmCIF files in pdb_cache,
merges with pdb_crystal_census.csv, and outputs an enriched CSV.

Usage:
    python extract_pdb_resolution.py

Expects:
    G:\Mechanical Taxonomy\pdb_cache\*.cif
    G:\Mechanical Taxonomy\pdb_crystal_census.csv

Outputs:
    G:\Mechanical Taxonomy\pdb_census_with_resolution.csv
"""

import os
import csv
import re
import sys

PDB_CACHE = r"G:\Mechanical Taxonomy\pdb_cache"
CENSUS_FILE = r"G:\Mechanical Taxonomy\pdb_crystal_census.csv"
OUTPUT_FILE = r"G:\Mechanical Taxonomy\pdb_census_with_resolution.csv"


def extract_resolution(cif_path):
    """Extract resolution from a CIF/mmCIF file.
    
    Tries multiple CIF tags used by different PDB depositions:
      _refine.ls_d_res_high          (most common)
      _reflns.d_resolution_high      (alternative)
      _em_3d_reconstruction.resolution (cryo-EM)
    """
    resolution = None
    try:
        with open(cif_path, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.strip()
                # X-ray resolution (most common)
                if line.startswith('_refine.ls_d_res_high'):
                    parts = line.split()
                    if len(parts) >= 2 and parts[1] not in ('?', '.'):
                        try:
                            resolution = float(parts[1])
                            return resolution
                        except ValueError:
                            pass
                # Alternative tag
                elif line.startswith('_reflns.d_resolution_high'):
                    parts = line.split()
                    if len(parts) >= 2 and parts[1] not in ('?', '.'):
                        try:
                            resolution = float(parts[1])
                            return resolution
                        except ValueError:
                            pass
                # Cryo-EM
                elif line.startswith('_em_3d_reconstruction.resolution'):
                    parts = line.split()
                    if len(parts) >= 2 and parts[1] not in ('?', '.'):
                        try:
                            resolution = float(parts[1])
                            return resolution
                        except ValueError:
                            pass
    except Exception as e:
        print(f"  Warning: could not read {cif_path}: {e}")
    return resolution


def main():
    # Step 1: Extract resolution from all CIF files
    print(f"Scanning {PDB_CACHE} for CIF files...")
    cif_files = [f for f in os.listdir(PDB_CACHE) if f.endswith('.cif')]
    print(f"Found {len(cif_files)} CIF files")
    
    # Map PDB ID -> resolution
    pdb_resolution = {}
    found = 0
    missing = 0
    
    for i, fname in enumerate(cif_files):
        # Extract PDB ID from filename (e.g., "1abc.cif" -> "1abc")
        pdb_id = fname.replace('.cif', '').lower()
        
        res = extract_resolution(os.path.join(PDB_CACHE, fname))
        if res is not None:
            pdb_resolution[pdb_id] = res
            found += 1
        else:
            missing += 1
        
        if (i + 1) % 100 == 0:
            print(f"  Processed {i+1}/{len(cif_files)} ({found} with resolution, {missing} without)")
    
    print(f"\nDone: {found} with resolution, {missing} without")
    
    # Step 2: Merge with census
    print(f"\nMerging with {CENSUS_FILE}...")
    
    rows_out = []
    matched = 0
    unmatched = 0
    
    with open(CENSUS_FILE, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames) + ['resolution']
        
        for row in reader:
            pdb_id = row.get('pdb_id', '').lower().strip()
            res = pdb_resolution.get(pdb_id, None)
            row['resolution'] = res if res is not None else ''
            rows_out.append(row)
            if res is not None:
                matched += 1
            else:
                unmatched += 1
    
    print(f"Matched: {matched}, Unmatched: {unmatched}")
    
    # Step 3: Write enriched CSV
    with open(OUTPUT_FILE, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)
    
    print(f"\nOutput written to: {OUTPUT_FILE}")
    
    # Step 4: Quick summary stats
    resolutions = [pdb_resolution[k] for k in pdb_resolution]
    if resolutions:
        resolutions.sort()
        import statistics
        print(f"\nResolution summary ({len(resolutions)} structures):")
        print(f"  Min:    {min(resolutions):.2f} A")
        print(f"  Max:    {max(resolutions):.2f} A")
        print(f"  Mean:   {statistics.mean(resolutions):.2f} A")
        print(f"  Median: {statistics.median(resolutions):.2f} A")
        print(f"  < 1.5A: {len([r for r in resolutions if r < 1.5])}")
        print(f"  < 2.0A: {len([r for r in resolutions if r < 2.0])}")
        print(f"  < 2.5A: {len([r for r in resolutions if r < 2.5])}")
        print(f"  > 3.0A: {len([r for r in resolutions if r > 3.0])}")


if __name__ == '__main__':
    main()

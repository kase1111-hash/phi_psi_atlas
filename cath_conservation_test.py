#!/usr/bin/env python3
"""
CATH Fold-Family Conservation Test
===================================
Downloads CATH-Gene3D domain annotations for UniProt IDs,
merges with curvature survey CSV, and computes within-family
conservation (CV of Σκ²/residue).

USAGE:
    python cath_conservation_test.py human_curvature_dssp.csv

    # Skip download if you already have the mapping file:
    python cath_conservation_test.py human_curvature_dssp.csv --mapping cath_mapping.csv

OUTPUT:
    cath_conservation_results.csv  — per-family conservation metrics
    cath_summary.txt               — summary statistics
    Prints results to stdout

REQUIRES: pandas, requests (for download)

LICENSE: CC0 1.0
"""

import os
import sys
import csv
import time
import argparse
import numpy as np
import pandas as pd


def download_cath_gene3d_mapping(uniprot_ids, output_path="cath_mapping.csv"):
    """Download CATH-Gene3D superfamily assignments for UniProt IDs.
    
    Uses the InterPro API to get Gene3D (CATH) annotations.
    Falls back to batch UniProt API if InterPro is slow.
    """
    import requests
    
    print(f"Downloading CATH-Gene3D annotations for {len(uniprot_ids)} proteins...")
    print(f"This may take 10-30 minutes depending on API speed.")
    
    t0 = time.time()
    results = []
    batch_size = 100
    ids = list(uniprot_ids)
    
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        
        for uid in batch:
            try:
                # UniProt API: get cross-references including Gene3D
                url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
                r = requests.get(url, timeout=10)
                if r.status_code != 200:
                    continue
                
                data = r.json()
                
                # Extract Gene3D (CATH) annotations
                xrefs = data.get('uniProtKBCrossReferences', [])
                for xref in xrefs:
                    if xref.get('database') == 'Gene3D':
                        cath_id = xref.get('id', '')
                        # Gene3D IDs look like: 1.10.10.10
                        # Format: Class.Architecture.Topology.Homology
                        parts = cath_id.split('.')
                        if len(parts) >= 3:
                            results.append({
                                'uniprot_id': uid,
                                'cath_superfamily': cath_id,
                                'cath_class': parts[0],
                                'cath_architecture': '.'.join(parts[:2]),
                                'cath_topology': '.'.join(parts[:3]),
                            })
                
                # Also try Pfam for broader family info
                for xref in xrefs:
                    if xref.get('database') == 'Pfam':
                        pfam_id = xref.get('id', '')
                        # Only add if we didn't get CATH
                        existing = [r for r in results if r['uniprot_id'] == uid]
                        if not existing:
                            results.append({
                                'uniprot_id': uid,
                                'cath_superfamily': f'Pfam:{pfam_id}',
                                'cath_class': 'Pfam',
                                'cath_architecture': 'Pfam',
                                'cath_topology': f'Pfam:{pfam_id}',
                            })
                
            except Exception as e:
                continue
        
        # Progress
        done = min(i + batch_size, len(ids))
        if done % 500 == 0 or done == len(ids):
            rate = done / max(time.time() - t0, 1)
            eta = (len(ids) - done) / max(rate, 0.1)
            print(f"  [{done}/{len(ids)}] {len(results)} annotations found "
                  f"({rate:.0f}/s, ETA {eta/60:.0f}min)")
    
    # Save
    df = pd.DataFrame(results)
    df.to_csv(output_path, index=False)
    print(f"Saved {len(df)} annotations to {output_path}")
    return df


def download_cath_gene3d_bulk(uniprot_ids, output_path="cath_mapping.csv"):
    """Alternative: download the CATH-Gene3D domain list directly.
    
    This is faster — downloads the pre-computed mapping file from Gene3D.
    """
    import requests
    
    print("Downloading Gene3D domain assignment file...")
    print("(This is a ~200MB file, may take a few minutes)")
    
    # Gene3D provides bulk downloads
    # URL for latest release
    urls_to_try = [
        "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry_protein_pairs.tsv.gz",
        "http://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-domain-boundaries-seqreschopping.txt",
    ]
    
    # Simpler approach: use InterPro's protein2ipr mapping
    url = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz"
    
    print(f"Attempting bulk download from InterPro...")
    print(f"URL: {url}")
    print(f"If this fails, use --manual mode (see below)")
    
    try:
        import gzip
        import io
        
        r = requests.get(url, stream=True, timeout=30)
        r.raise_for_status()
        
        # This file is large (~5GB uncompressed), so we filter on the fly
        target_ids = set(uniprot_ids)
        results = []
        
        print("Streaming and filtering...")
        with gzip.open(io.BytesIO(r.content), 'rt') as f:
            for line_num, line in enumerate(f):
                if line_num % 10_000_000 == 0:
                    print(f"  Processed {line_num/1e6:.0f}M lines, {len(results)} hits...")
                
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                
                uid = parts[0]
                if uid not in target_ids:
                    continue
                
                ipr_id = parts[1]
                ipr_name = parts[2]
                db_name = parts[3]
                
                # Filter for CATH-Gene3D or Pfam entries
                if 'Gene3D' in db_name or 'G3DSA' in db_name:
                    db_id = parts[4] if len(parts) > 4 else ''
                    cath_parts = db_id.replace('G3DSA:', '').split('.')
                    results.append({
                        'uniprot_id': uid,
                        'cath_superfamily': db_id,
                        'cath_class': cath_parts[0] if len(cath_parts) > 0 else '',
                        'cath_architecture': '.'.join(cath_parts[:2]) if len(cath_parts) >= 2 else '',
                        'cath_topology': '.'.join(cath_parts[:3]) if len(cath_parts) >= 3 else '',
                    })
        
        df = pd.DataFrame(results)
        df.to_csv(output_path, index=False)
        print(f"Saved {len(df)} annotations to {output_path}")
        return df
        
    except Exception as e:
        print(f"Bulk download failed: {e}")
        print(f"\nFalling back to per-protein API queries...")
        return download_cath_gene3d_mapping(uniprot_ids, output_path)


def analyze_conservation(curv_df, cath_df, min_family_size=20):
    """Compute within-family conservation metrics."""
    
    # Merge
    merged = curv_df.merge(cath_df, on='uniprot_id', how='inner')
    print(f"\nMerged: {len(merged)} protein-domain pairs "
          f"({merged['uniprot_id'].nunique()} unique proteins)")
    
    # Analyze at topology level (Class.Architecture.Topology)
    results = []
    
    for level_name, level_col in [('cath_topology', 'cath_topology'), 
                                    ('cath_architecture', 'cath_architecture'),
                                    ('cath_superfamily', 'cath_superfamily')]:
        for family_id, group in merged.groupby(level_col):
            # Use unique proteins only (some have multiple domains)
            unique = group.drop_duplicates(subset='uniprot_id')
            n = len(unique)
            if n < min_family_size:
                continue
            
            ksq = unique['kappa_sq_per_res']
            cv = ksq.std() / ksq.mean() * 100 if ksq.mean() > 0 else np.nan
            
            results.append({
                'level': level_name,
                'family_id': family_id,
                'n_proteins': n,
                'mean_kappa_sq_per_res': round(ksq.mean(), 1),
                'median_kappa_sq_per_res': round(ksq.median(), 1),
                'std_kappa_sq_per_res': round(ksq.std(), 1),
                'cv_percent': round(cv, 1),
                'mean_length': round(unique['length'].mean(), 0),
                'mean_helix': round(unique['ss_frac_H'].mean(), 3),
                'mean_sheet': round(unique['ss_frac_E'].mean(), 3),
                'mean_plddt': round(unique['mean_plddt'].mean(), 1),
            })
    
    return pd.DataFrame(results)


def print_results(results_df):
    """Print summary of conservation results."""
    
    for level in ['cath_topology', 'cath_architecture', 'cath_superfamily']:
        sub = results_df[results_df['level'] == level].sort_values('cv_percent')
        if len(sub) == 0:
            continue
        
        print(f"\n{'='*70}")
        print(f"  CONSERVATION AT {level.upper()} LEVEL ({len(sub)} families with n >= 20)")
        print(f"{'='*70}")
        
        print(f"\n  CV Distribution:")
        for pct in [10, 25, 50, 75, 90]:
            print(f"    P{pct:02d}: {sub['cv_percent'].quantile(pct/100):.1f}%")
        
        tight = sub[sub['cv_percent'] < 30]
        print(f"\n  Families with CV < 30%: {len(tight)} / {len(sub)} ({len(tight)/len(sub)*100:.1f}%)")
        
        print(f"\n  TOP 15 TIGHTEST FAMILIES:")
        print(f"  {'Family':>25s} {'n':>5s} {'mean':>7s} {'CV%':>6s} {'H%':>5s} {'E%':>5s} {'L':>5s} {'pLDDT':>6s}")
        for _, r in sub.head(15).iterrows():
            print(f"  {str(r['family_id']):>25s} {r['n_proteins']:>5.0f} {r['mean_kappa_sq_per_res']:>7.1f} "
                  f"{r['cv_percent']:>5.1f}% {r['mean_helix']*100:>4.0f}% {r['mean_sheet']*100:>4.0f}% "
                  f"{r['mean_length']:>5.0f} {r['mean_plddt']:>5.1f}")
        
        print(f"\n  BOTTOM 5 LOOSEST FAMILIES:")
        for _, r in sub.tail(5).iterrows():
            print(f"  {str(r['family_id']):>25s} {r['n_proteins']:>5.0f} {r['mean_kappa_sq_per_res']:>7.1f} "
                  f"{r['cv_percent']:>5.1f}% {r['mean_helix']*100:>4.0f}% {r['mean_sheet']*100:>4.0f}% "
                  f"{r['mean_length']:>5.0f} {r['mean_plddt']:>5.1f}")
    
    # Compare to hemoglobin
    print(f"\n{'='*70}")
    print(f"  REFERENCE: Hemoglobin crystal structures")
    print(f"  CV = 16%, n = 10, mean Σκ²/res = 107")
    print(f"{'='*70}")
    
    # Overall summary
    topo = results_df[results_df['level'] == 'cath_topology']
    if len(topo) > 0:
        print(f"\n  OVERALL SUMMARY (topology level):")
        print(f"  {len(topo)} fold families analyzed")
        print(f"  Median within-family CV: {topo['cv_percent'].median():.1f}%")
        tight = topo[topo['cv_percent'] < 20]
        print(f"  Families with CV < 20% (hemoglobin-like conservation): {len(tight)}")
        if len(tight) > 0:
            for _, r in tight.iterrows():
                print(f"    {r['family_id']}: CV={r['cv_percent']:.1f}%, n={r['n_proteins']:.0f}")


def main():
    parser = argparse.ArgumentParser(description='CATH fold-family conservation test')
    parser.add_argument('curvature_csv', help='Curvature survey CSV (from proteome_curvature_survey.py)')
    parser.add_argument('--mapping', help='Pre-existing CATH mapping CSV (skip download)')
    parser.add_argument('--min-family', type=int, default=20,
                        help='Minimum proteins per family (default: 20)')
    parser.add_argument('--manual', action='store_true',
                        help='Use per-protein API queries instead of bulk download')
    parser.add_argument('--output', default='cath_conservation_results.csv',
                        help='Output CSV for results')
    args = parser.parse_args()
    
    # Load curvature data
    print(f"Loading curvature data from {args.curvature_csv}...")
    curv = pd.read_csv(args.curvature_csv)
    
    # Quality filter
    hq = curv[(curv['mean_plddt'] >= 70) & (curv['length'] >= 50) & (curv['length'] <= 2000)].copy()
    print(f"  {len(hq)} high-quality proteins (pLDDT >= 70, length 50-2000)")
    
    # Get CATH annotations
    if args.mapping and os.path.exists(args.mapping):
        print(f"Loading pre-existing CATH mapping from {args.mapping}...")
        cath = pd.read_csv(args.mapping)
    else:
        t0 = time.time()
        uniprot_ids = hq['uniprot_id'].unique()
        
        if args.manual:
            cath = download_cath_gene3d_mapping(uniprot_ids)
        else:
            cath = download_cath_gene3d_bulk(uniprot_ids)
    
    # Analyze
    results = analyze_conservation(hq, cath, min_family_size=args.min_family)
    
    # Save
    results.to_csv(args.output, index=False)
    print(f"\nSaved results to {args.output}")
    
    # Print
    print_results(results)


if __name__ == '__main__':
    main()

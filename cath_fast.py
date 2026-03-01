#!/usr/bin/env python3
"""
Fast CATH Fold-Family Conservation Test
========================================
Uses UniProt batch API (500 IDs per request) to fetch CATH-Gene3D
annotations ~500x faster than per-protein queries.

USAGE:
    python cath_fast.py human_curvature_dssp.csv

    # If you already downloaded annotations:
    python cath_fast.py human_curvature_dssp.csv --mapping cath_mapping.csv

OUTPUT:
    cath_mapping.csv               — raw annotations
    cath_conservation_results.csv  — per-family conservation metrics

REQUIRES: pandas, requests
LICENSE: CC0 1.0
"""

import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
import requests


def fetch_cath_batch(uniprot_ids, output_path="cath_mapping.csv"):
    """Fetch CATH-Gene3D annotations using UniProt REST batch search.
    
    Queries ~500 IDs per request. ~13K proteins in ~3-5 minutes.
    """
    t0 = time.time()
    all_ids = list(uniprot_ids)
    total = len(all_ids)
    batch_size = 500
    results = []
    
    print(f"Fetching CATH-Gene3D annotations for {total} proteins...")
    print(f"Using UniProt batch API ({batch_size} IDs per request)")
    
    for i in range(0, total, batch_size):
        batch = all_ids[i:i+batch_size]
        query = ' OR '.join(f'accession:{uid}' for uid in batch)
        
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': query,
            'fields': 'accession,xref_gene3d,xref_pfam',
            'format': 'tsv',
            'size': batch_size,
        }
        
        try:
            r = requests.get(url, params=params, timeout=60)
            if r.status_code == 200:
                lines = r.text.strip().split('\n')
                if len(lines) > 1:
                    header = lines[0].split('\t')
                    for line in lines[1:]:
                        fields = line.split('\t')
                        if len(fields) < len(header):
                            continue
                        
                        row = dict(zip(header, fields))
                        uid = row.get('Entry', '')
                        
                        # Parse Gene3D (CATH) annotations
                        gene3d = row.get('Gene3D', '')
                        if gene3d and gene3d != '':
                            for entry in gene3d.split(';'):
                                entry = entry.strip()
                                if not entry:
                                    continue
                                # Format: "G3DSA:1.10.10.10" or just "1.10.10.10"
                                cath_id = entry.replace('G3DSA:', '').strip()
                                parts = cath_id.split('.')
                                if len(parts) >= 3:
                                    results.append({
                                        'uniprot_id': uid,
                                        'source': 'Gene3D',
                                        'cath_superfamily': cath_id,
                                        'cath_class': parts[0],
                                        'cath_architecture': '.'.join(parts[:2]),
                                        'cath_topology': '.'.join(parts[:3]),
                                    })
                        
                        # Parse Pfam as backup
                        pfam = row.get('Pfam', '')
                        if pfam and pfam != '':
                            for entry in pfam.split(';'):
                                entry = entry.strip()
                                if not entry:
                                    continue
                                results.append({
                                    'uniprot_id': uid,
                                    'source': 'Pfam',
                                    'cath_superfamily': f'Pfam:{entry}',
                                    'cath_class': 'Pfam',
                                    'cath_architecture': 'Pfam',
                                    'cath_topology': f'Pfam:{entry}',
                                })
            
            elif r.status_code == 400:
                # Query too long, split into smaller batches
                for uid in batch:
                    try:
                        r2 = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.json", timeout=10)
                        if r2.status_code == 200:
                            data = r2.json()
                            for xref in data.get('uniProtKBCrossReferences', []):
                                if xref.get('database') == 'Gene3D':
                                    cath_id = xref.get('id', '')
                                    parts = cath_id.split('.')
                                    if len(parts) >= 3:
                                        results.append({
                                            'uniprot_id': uid,
                                            'source': 'Gene3D',
                                            'cath_superfamily': cath_id,
                                            'cath_class': parts[0],
                                            'cath_architecture': '.'.join(parts[:2]),
                                            'cath_topology': '.'.join(parts[:3]),
                                        })
                    except:
                        pass
        
        except requests.exceptions.Timeout:
            print(f"  Timeout at batch {i//batch_size}, retrying with smaller batch...")
            # Retry one at a time
            for uid in batch:
                try:
                    r2 = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.json", timeout=10)
                    if r2.status_code == 200:
                        data = r2.json()
                        for xref in data.get('uniProtKBCrossReferences', []):
                            if xref.get('database') == 'Gene3D':
                                cath_id = xref.get('id', '')
                                parts = cath_id.split('.')
                                if len(parts) >= 3:
                                    results.append({
                                        'uniprot_id': uid,
                                        'source': 'Gene3D',
                                        'cath_superfamily': cath_id,
                                        'cath_class': parts[0],
                                        'cath_architecture': '.'.join(parts[:2]),
                                        'cath_topology': '.'.join(parts[:3]),
                                    })
                except:
                    pass
        
        except Exception as e:
            print(f"  Error at batch {i//batch_size}: {e}")
        
        # Progress
        done = min(i + batch_size, total)
        elapsed = time.time() - t0
        rate = done / max(elapsed, 1)
        eta = (total - done) / max(rate, 0.1)
        
        if done % 2000 == 0 or done >= total or done <= batch_size:
            n_unique = len(set(r['uniprot_id'] for r in results))
            print(f"  [{done}/{total}] {len(results)} annotations ({n_unique} proteins) "
                  f"({rate:.0f} IDs/s, ETA {eta/60:.1f}min)")
        
        # Small delay to be polite to the API
        time.sleep(0.5)
    
    df = pd.DataFrame(results)
    df.to_csv(output_path, index=False)
    
    elapsed = time.time() - t0
    n_unique = df['uniprot_id'].nunique() if len(df) > 0 else 0
    print(f"\nDone in {elapsed/60:.1f} minutes")
    print(f"  {len(df)} total annotations for {n_unique} unique proteins")
    print(f"  Saved to {output_path}")
    
    return df


def analyze_conservation(curv_df, cath_df, min_family_size=20):
    """Compute within-family conservation metrics."""
    
    # Keep only Gene3D (CATH) annotations for fold analysis
    cath_gene3d = cath_df[cath_df['source'] == 'Gene3D'].copy()
    print(f"\nGene3D annotations: {len(cath_gene3d)} ({cath_gene3d['uniprot_id'].nunique()} proteins)")
    
    # Merge
    merged = curv_df.merge(cath_gene3d, on='uniprot_id', how='inner')
    print(f"Merged: {len(merged)} protein-domain pairs "
          f"({merged['uniprot_id'].nunique()} unique proteins)")
    
    results = []
    
    for level_name in ['cath_topology', 'cath_architecture', 'cath_class']:
        for family_id, group in merged.groupby(level_name):
            unique = group.drop_duplicates(subset='uniprot_id')
            n = len(unique)
            if n < min_family_size:
                continue
            
            ksq = unique['kappa_sq_per_res']
            mean_val = ksq.mean()
            cv = ksq.std() / mean_val * 100 if mean_val > 0 else np.nan
            
            results.append({
                'level': level_name,
                'family_id': family_id,
                'n_proteins': n,
                'mean_kappa_sq_per_res': round(mean_val, 1),
                'median_kappa_sq_per_res': round(ksq.median(), 1),
                'std_kappa_sq_per_res': round(ksq.std(), 1),
                'cv_percent': round(cv, 1),
                'mean_length': round(unique['length'].mean(), 0),
                'std_length': round(unique['length'].std(), 0),
                'mean_helix': round(unique['ss_frac_H'].mean(), 3),
                'mean_sheet': round(unique['ss_frac_E'].mean(), 3),
                'mean_plddt': round(unique['mean_plddt'].mean(), 1),
            })
    
    return pd.DataFrame(results)


def print_results(results_df):
    """Print conservation analysis results."""
    
    for level in ['cath_topology', 'cath_architecture', 'cath_class']:
        sub = results_df[results_df['level'] == level].sort_values('cv_percent')
        if len(sub) == 0:
            continue
        
        print(f"\n{'='*75}")
        print(f"  CONSERVATION AT {level.upper()} LEVEL ({len(sub)} families with n >= 20)")
        print(f"{'='*75}")
        
        print(f"\n  CV Distribution:")
        for pct in [5, 10, 25, 50, 75, 90]:
            print(f"    P{pct:02d}: {sub['cv_percent'].quantile(pct/100):.1f}%")
        
        tight20 = sub[sub['cv_percent'] < 20]
        tight30 = sub[sub['cv_percent'] < 30]
        print(f"\n  Families with CV < 20%: {len(tight20)} / {len(sub)}")
        print(f"  Families with CV < 30%: {len(tight30)} / {len(sub)}")
        
        print(f"\n  TOP 20 TIGHTEST FAMILIES:")
        print(f"  {'Family':>25s} {'n':>5s} {'mean':>7s} {'med':>7s} {'CV%':>6s} "
              f"{'H%':>5s} {'E%':>5s} {'L':>5s} {'L_sd':>5s}")
        for _, r in sub.head(20).iterrows():
            print(f"  {str(r['family_id']):>25s} {r['n_proteins']:>5.0f} "
                  f"{r['mean_kappa_sq_per_res']:>7.1f} {r['median_kappa_sq_per_res']:>7.1f} "
                  f"{r['cv_percent']:>5.1f}% {r['mean_helix']*100:>4.0f}% "
                  f"{r['mean_sheet']*100:>4.0f}% {r['mean_length']:>5.0f} {r['std_length']:>5.0f}")
    
    # Hemoglobin reference
    print(f"\n{'='*75}")
    print(f"  REFERENCE: Hemoglobin \u03B2-globin crystal structures")
    print(f"  CV = 16%, n = 10, mean \u03A3\u03BA\u00B2/res = 107, length = 143")
    print(f"  CATH topology: 1.10.490 (Globin-like)")
    print(f"{'='*75}")
    
    # Look for globin family specifically
    topo = results_df[results_df['level'] == 'cath_topology']
    globin = topo[topo['family_id'].str.startswith('1.10.490')]
    if len(globin) > 0:
        print(f"\n  GLOBIN FAMILY IN DATA:")
        for _, r in globin.iterrows():
            print(f"    {r['family_id']}: n={r['n_proteins']:.0f}, CV={r['cv_percent']:.1f}%, "
                  f"mean={r['mean_kappa_sq_per_res']:.1f}")
    else:
        # Try broader search
        globin_like = topo[topo['family_id'].str.startswith('1.10')]
        if len(globin_like) > 0:
            print(f"\n  CATH 1.10.* families (Mainly Alpha, Orthogonal Bundle):")
            for _, r in globin_like.sort_values('cv_percent').head(10).iterrows():
                print(f"    {r['family_id']}: n={r['n_proteins']:.0f}, CV={r['cv_percent']:.1f}%, "
                      f"mean={r['mean_kappa_sq_per_res']:.1f}")


def main():
    parser = argparse.ArgumentParser(description='Fast CATH fold-family conservation test')
    parser.add_argument('curvature_csv', help='Curvature survey CSV')
    parser.add_argument('--mapping', help='Pre-existing CATH mapping CSV (skip download)')
    parser.add_argument('--min-family', type=int, default=20,
                        help='Minimum proteins per family (default: 20)')
    parser.add_argument('--output', default='cath_conservation_results.csv',
                        help='Output CSV for results')
    args = parser.parse_args()
    
    # Load curvature data
    print(f"Loading curvature data from {args.curvature_csv}...")
    curv = pd.read_csv(args.curvature_csv)
    hq = curv[(curv['mean_plddt'] >= 70) & (curv['length'] >= 50) & (curv['length'] <= 2000)].copy()
    print(f"  {len(hq)} high-quality proteins")
    
    # Get annotations
    if args.mapping and os.path.exists(args.mapping):
        print(f"Loading existing mapping from {args.mapping}...")
        cath = pd.read_csv(args.mapping)
    else:
        cath = fetch_cath_batch(hq['uniprot_id'].unique(), output_path='cath_mapping.csv')
    
    # Analyze
    results = analyze_conservation(hq, cath, min_family_size=args.min_family)
    results.to_csv(args.output, index=False)
    print(f"\nSaved results to {args.output}")
    
    # Print
    print_results(results)


if __name__ == '__main__':
    main()

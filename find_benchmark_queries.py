#!/usr/bin/env python3
"""
Find benchmark query candidates from the EXISTING barcode database.

Instead of fetching specific proteins, this scans all barcoded proteins,
fetches their classifications from UniProt, and picks diverse representatives
across structural classes.

Usage:
    python find_benchmark_queries.py              # Find and classify candidates
    python find_benchmark_queries.py --no-fetch   # Use cached classifications only
    python find_benchmark_queries.py --report     # Show cached results
"""

import json
import os
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path
from collections import defaultdict

SCRIPT_DIR = Path(__file__).parent
BARCODE_DIR = SCRIPT_DIR / "results" / "barcodes"
CACHE_DIR = SCRIPT_DIR / "data" / "classifications"
OUTPUT_FILE = SCRIPT_DIR / "data" / "benchmark_queries.json"

sys.path.insert(0, str(SCRIPT_DIR))
from geom_homology import load_all_barcodes


def fetch_uniprot_basic(uid, timeout=10):
    """Fetch basic protein info from UniProt API."""
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
    try:
        req = urllib.request.Request(url, headers={
            'User-Agent': 'Mozilla/5.0 PhiPsiAtlas/1.0',
            'Accept': 'application/json',
        })
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode('utf-8', errors='replace'))
        
        result = {
            'uniprot_id': uid,
            'protein_name': '',
            'gene_name': '',
            'scop2_superfamily': '',
            'cath_topology': '',
            'pfam_families': [],
            'keywords': [],
            'subcellular': [],
        }
        
        # Protein name
        rec = data.get('proteinDescription', {}).get('recommendedName', {})
        if rec:
            result['protein_name'] = rec.get('fullName', {}).get('value', '')
        
        # Gene
        genes = data.get('genes', [])
        if genes:
            result['gene_name'] = genes[0].get('geneName', {}).get('value', '')
        
        # Cross-refs
        for xref in data.get('uniProtKBCrossReferences', []):
            db = xref.get('database', '')
            xid = xref.get('id', '')
            props = {p.get('key',''): p.get('value','') for p in xref.get('properties',[])}
            
            if db == 'Pfam':
                result['pfam_families'].append({'id': xid, 'name': props.get('EntryName','')})
            elif db == 'SUPFAM':
                if not result['scop2_superfamily']:
                    result['scop2_superfamily'] = xid
            elif db == 'Gene3D':
                if not result['cath_topology']:
                    result['cath_topology'] = xid
        
        # Keywords (for structural class hints)
        for kw in data.get('keywords', []):
            result['keywords'].append(kw.get('name', ''))
        
        # Subcellular location
        for comment in data.get('comments', []):
            if comment.get('commentType') == 'SUBCELLULAR LOCATION':
                for sl in comment.get('subcellularLocations', []):
                    loc = sl.get('location', {}).get('value', '')
                    if loc:
                        result['subcellular'].append(loc)
        
        return result
    except Exception as e:
        return {'uniprot_id': uid, 'protein_name': '?', 'gene_name': '', 
                'error': str(e), 'pfam_families': [], 'keywords': []}


def classify_structural_class(bc, cls_data):
    """Assign a structural class based on barcode SS and keywords."""
    # From barcode
    ss = bc.get('ss_counts', {})
    total = sum(ss.values()) or 1
    h_frac = ss.get('H', 0) / total
    e_frac = ss.get('E', 0) / total
    
    keywords = set(cls_data.get('keywords', []))
    subcel = ' '.join(cls_data.get('subcellular', []))
    
    # Structural class assignment
    if 'Transmembrane' in keywords or 'transmembrane' in subcel.lower() or 'membrane' in subcel.lower():
        if h_frac > 0.4:
            return 'TM-alpha'
        elif e_frac > 0.2:
            return 'TM-beta'
        else:
            return 'TM-mixed'
    
    if h_frac > 0.5 and e_frac < 0.1:
        return 'all-alpha'
    elif e_frac > 0.3 and h_frac < 0.15:
        return 'all-beta'
    elif h_frac > 0.2 and e_frac > 0.15:
        return 'alpha-beta'
    elif h_frac < 0.15 and e_frac < 0.15:
        return 'coil-rich'
    else:
        return 'mixed'


def pick_diverse_queries(candidates, n_per_class=3, max_total=20):
    """Pick diverse representative queries across structural classes."""
    by_class = defaultdict(list)
    for uid, info in candidates.items():
        by_class[info['structural_class']].append((uid, info))
    
    # Sort each class by number of Pfam domains (prefer well-annotated)
    # and by protein length (prefer medium-sized, 200-800 aa)
    def score(item):
        uid, info = item
        n_pfam = len(info.get('pfam_families', []))
        length = info.get('length', 500)
        # Prefer 200-800 aa, well-annotated
        len_score = -abs(length - 500) / 500
        return n_pfam + len_score
    
    selected = {}
    for cls in sorted(by_class.keys()):
        entries = sorted(by_class[cls], key=score, reverse=True)
        for uid, info in entries[:n_per_class]:
            if len(selected) < max_total:
                selected[uid] = info
    
    return selected


def main():
    args = sys.argv[1:]
    no_fetch = '--no-fetch' in args
    report_only = '--report' in args
    
    if report_only:
        if OUTPUT_FILE.exists():
            with open(OUTPUT_FILE) as f:
                data = json.load(f)
            print(f"\n  Cached benchmark queries: {len(data)}")
            for uid, info in data.items():
                print(f"  {uid:>10s} {info.get('gene_name','?'):>8s} "
                      f"{info.get('structural_class','?'):>12s} "
                      f"L={info.get('length',0):5d} "
                      f"Pfam={len(info.get('pfam_families',[]))} "
                      f"{info.get('protein_name','?')[:40]}")
        else:
            print("  No cached queries. Run: python find_benchmark_queries.py")
        return
    
    # Load barcodes
    barcodes = load_all_barcodes()
    if not barcodes:
        print("  No barcodes found.")
        return
    print(f"  Loaded {len(barcodes)} barcodes")
    
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    
    # Classify a strategic sample — take every Nth protein to get broad coverage
    # Plus any that are already cached
    sample_size = min(300, len(barcodes))
    step = max(1, len(barcodes) // sample_size)
    sample_indices = list(range(0, len(barcodes), step))[:sample_size]
    
    # Also include proteins with interesting barcode properties
    # (high defect, high geodesic, extreme winding, various sizes)
    by_length = sorted(range(len(barcodes)), key=lambda i: barcodes[i].get('length', 0))
    # Add small, medium, large representatives
    for frac in [0.05, 0.15, 0.25, 0.5, 0.75, 0.85, 0.95]:
        idx = int(frac * len(by_length))
        if idx < len(by_length) and idx not in sample_indices:
            sample_indices.append(by_length[idx])
    
    sample_uids = []
    for i in sample_indices:
        uid = barcodes[i].get('uniprot_id', '')
        if uid:
            sample_uids.append((uid, barcodes[i]))
    
    print(f"  Classifying {len(sample_uids)} proteins...")
    
    candidates = {}
    n_fetched = 0
    n_cached = 0
    
    for uid, bc in sample_uids:
        # Check cache first
        cache_file = CACHE_DIR / f"{uid}.json"
        if cache_file.exists():
            try:
                with open(cache_file) as f:
                    cls_data = json.load(f)
                n_cached += 1
            except (json.JSONDecodeError, IOError):
                cls_data = None
        else:
            cls_data = None
        
        if cls_data is None and not no_fetch:
            cls_data = fetch_uniprot_basic(uid)
            if cls_data and cls_data.get('protein_name', '?') != '?':
                with open(cache_file, 'w') as f:
                    json.dump(cls_data, f, indent=1)
                n_fetched += 1
            time.sleep(0.3)  # Rate limit
        
        if cls_data is None:
            cls_data = {'uniprot_id': uid, 'protein_name': '?', 'pfam_families': [], 'keywords': []}
        
        struct_class = classify_structural_class(bc, cls_data)
        
        info = {
            'uniprot_id': uid,
            'protein_name': cls_data.get('protein_name', '?'),
            'gene_name': cls_data.get('gene_name', ''),
            'structural_class': struct_class,
            'length': bc.get('length', 0),
            'n_segments': bc.get('n_segments', len(bc.get('segments', []))),
            'scop2_superfamily': cls_data.get('scop2_superfamily', ''),
            'cath_topology': cls_data.get('cath_topology', ''),
            'pfam_families': cls_data.get('pfam_families', []),
            'geodesic_frac': bc.get('geodesic_fraction', 0),
            'defect_frac': bc.get('defect_fraction', 0),
            'Q_magnitude': bc.get('Q_magnitude', 0),
        }
        candidates[uid] = info
        
        if (n_fetched + n_cached) % 50 == 0 and (n_fetched + n_cached) > 0:
            print(f"    Processed {n_fetched + n_cached}/{len(sample_uids)} "
                  f"({n_cached} cached, {n_fetched} fetched)")
    
    print(f"  Done: {n_cached} cached, {n_fetched} fetched, "
          f"{len(sample_uids) - n_cached - n_fetched} unclassified")
    
    # Show class distribution
    from collections import Counter
    class_counts = Counter(info['structural_class'] for info in candidates.values()
                          if info.get('protein_name', '?') != '?')
    
    print(f"\n  STRUCTURAL CLASS DISTRIBUTION ({len(candidates)} classified):")
    for cls, count in sorted(class_counts.items(), key=lambda x: -x[1]):
        print(f"    {cls:>12s}: {count:4d}")
    
    # Pick diverse queries
    # Only from well-annotated proteins (have Pfam and SCOP)
    annotated = {uid: info for uid, info in candidates.items()
                 if info.get('protein_name', '?') != '?'
                 and len(info.get('pfam_families', [])) > 0}
    
    print(f"\n  Well-annotated candidates: {len(annotated)}")
    
    selected = pick_diverse_queries(annotated, n_per_class=3, max_total=20)
    
    print(f"\n  SELECTED BENCHMARK QUERIES ({len(selected)}):")
    print(f"  {'UniProt':>10s} {'Gene':>8s} {'Class':>12s} {'Length':>6s} "
          f"{'Pfam':>5s} {'SCOP':>12s} {'Name'}")
    print(f"  {'-'*10} {'-'*8} {'-'*12} {'-'*6} {'-'*5} {'-'*12} {'-'*40}")
    
    for uid, info in sorted(selected.items(), key=lambda x: x[1]['structural_class']):
        scop = info.get('scop2_superfamily', '-')[:12]
        print(f"  {uid:>10s} {info['gene_name']:>8s} {info['structural_class']:>12s} "
              f"{info['length']:6d} {len(info['pfam_families']):5d} "
              f"{scop:>12s} {info['protein_name'][:40]}")
    
    # Save
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(selected, f, indent=2)
    print(f"\n  Saved to {OUTPUT_FILE}")
    
    # Generate the query list for benchmark_homology.py
    query_str = ','.join(selected.keys())
    print(f"\n  To run benchmark:")
    print(f"  python benchmark_homology.py run --queries={query_str}")


if __name__ == "__main__":
    main()

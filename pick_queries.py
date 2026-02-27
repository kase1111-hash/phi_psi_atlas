#!/usr/bin/env python3
"""
Pick diverse benchmark queries using ONLY barcode properties (no API needed).

Selects proteins from the existing database spanning:
  - Structural classes (helix-rich, strand-rich, mixed, coil-rich)
  - Size bins (small, medium, large, giant)
  - Curvature extremes (high defect, high geodesic, high |Q|)

Then adds EGFR and insulin (already benchmarked) for continuity.

Usage:
    python pick_queries.py                 # Pick and list queries
    python pick_queries.py --run           # Pick queries and immediately run benchmark
"""

import json
import sys
from pathlib import Path
from collections import defaultdict

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from geom_homology import load_all_barcodes


def classify_by_ss(bc):
    """Classify structural class from barcode SS counts."""
    ss = bc.get('ss_counts', {})
    total = sum(ss.values()) or 1
    h_frac = ss.get('H', 0) / total
    e_frac = ss.get('E', 0) / total
    c_frac = ss.get('C', 0) / total
    
    if h_frac > 0.50 and e_frac < 0.10:
        return 'helix-rich'
    elif e_frac > 0.30 and h_frac < 0.15:
        return 'strand-rich'
    elif h_frac > 0.20 and e_frac > 0.15:
        return 'alpha-beta'
    elif c_frac > 0.65:
        return 'coil-rich'
    else:
        return 'mixed'


def size_bin(length):
    if length < 150:
        return 'small'
    elif length < 400:
        return 'medium'
    elif length < 800:
        return 'large'
    else:
        return 'giant'


def main():
    barcodes = load_all_barcodes()
    print(f"  Loaded {len(barcodes)} barcodes")
    
    # Classify all
    classified = []
    for bc in barcodes:
        uid = bc.get('uniprot_id', '')
        if not uid:
            continue
        length = bc.get('length', 0)
        n_segs = bc.get('n_segments', len(bc.get('segments', [])))
        
        # Skip tiny proteins with < 3 segments (not useful for alignment)
        if n_segs < 5:
            continue
        
        ss = bc.get('ss_counts', {})
        total = sum(ss.values()) or 1
        
        info = {
            'uid': uid,
            'length': length,
            'n_segments': n_segs,
            'h_frac': ss.get('H', 0) / total,
            'e_frac': ss.get('E', 0) / total,
            'c_frac': ss.get('C', 0) / total,
            'geo_frac': bc.get('geodesic_fraction', 0),
            'def_frac': bc.get('defect_fraction', 0),
            'Q_mag': bc.get('Q_magnitude', 0),
            'ss_class': classify_by_ss(bc),
            'size': size_bin(length),
        }
        classified.append(info)
    
    print(f"  Classified {len(classified)} proteins (≥5 segments)")
    
    # Show class distribution
    from collections import Counter
    class_counts = Counter(p['ss_class'] for p in classified)
    print(f"\n  SS-based structural classes:")
    for cls, n in sorted(class_counts.items(), key=lambda x: -x[1]):
        print(f"    {cls:>12s}: {n:4d}")
    
    # Strategy: pick 3-4 per structural class, spread across size bins
    # Prefer proteins with moderate segment counts (not too few, not too many)
    # to get meaningful alignments
    
    selected = {}
    
    # For each SS class, pick one per size bin where possible
    for ss_class in ['helix-rich', 'strand-rich', 'alpha-beta', 'coil-rich', 'mixed']:
        pool = [p for p in classified if p['ss_class'] == ss_class]
        
        # Group by size
        by_size = defaultdict(list)
        for p in pool:
            by_size[p['size']].append(p)
        
        # From each available size bin, pick the one with most segments
        # (well-structured, good for alignment testing)
        picked = 0
        for sz in ['medium', 'large', 'giant', 'small']:
            if picked >= 3:
                break
            candidates = by_size.get(sz, [])
            if not candidates:
                continue
            # Sort by segment count (prefer medium complexity)
            candidates.sort(key=lambda p: -p['n_segments'])
            # Pick from middle of the list (not extreme)
            idx = min(len(candidates) // 3, len(candidates) - 1)
            best = candidates[idx]
            if best['uid'] not in selected:
                selected[best['uid']] = best
                picked += 1
    
    # Also ensure EGFR (P00533) and Insulin (P01308) are included if available
    for uid in ['P00533', 'P01308']:
        matches = [p for p in classified if p['uid'] == uid]
        if matches and uid not in selected:
            selected[uid] = matches[0]
    
    # Add a couple of extremes for interest
    # Highest defect fraction (most perturbed)
    high_def = sorted(classified, key=lambda p: -p['def_frac'])
    for p in high_def[:2]:
        if p['uid'] not in selected and len(selected) < 22:
            p['note'] = 'high-defect'
            selected[p['uid']] = p
    
    # Highest geodesic fraction (smoothest)
    high_geo = sorted([p for p in classified if p['n_segments'] >= 10],
                      key=lambda p: -p['geo_frac'])
    for p in high_geo[:2]:
        if p['uid'] not in selected and len(selected) < 22:
            p['note'] = 'high-geodesic'
            selected[p['uid']] = p
    
    # Print results
    print(f"\n  SELECTED BENCHMARK QUERIES ({len(selected)}):")
    print(f"  {'UID':>12s} {'Class':>12s} {'Size':>7s} {'Length':>6s} {'Segs':>5s} "
          f"{'H%':>5s} {'E%':>5s} {'C%':>5s} {'Geo%':>5s} {'Def%':>5s} {'|Q|':>6s} {'Note'}")
    print(f"  {'-'*12} {'-'*12} {'-'*7} {'-'*6} {'-'*5} "
          f"{'-'*5} {'-'*5} {'-'*5} {'-'*5} {'-'*5} {'-'*6} {'-'*15}")
    
    for uid in sorted(selected, key=lambda u: (selected[u]['ss_class'], selected[u]['length'])):
        p = selected[uid]
        note = p.get('note', '')
        print(f"  {uid:>12s} {p['ss_class']:>12s} {p['size']:>7s} {p['length']:6d} "
              f"{p['n_segments']:5d} {p['h_frac']:5.1%} {p['e_frac']:5.1%} "
              f"{p['c_frac']:5.1%} {p['geo_frac']:5.1%} {p['def_frac']:5.1%} "
              f"{p['Q_mag']:6.1f} {note}")
    
    query_list = ','.join(selected.keys())
    
    print(f"\n  Query IDs: {query_list}")
    print(f"\n  Run benchmark:")
    print(f"  python benchmark_homology.py run --queries={query_list}")
    
    # Save
    outfile = SCRIPT_DIR / "data" / "benchmark_queries_auto.json"
    outfile.parent.mkdir(parents=True, exist_ok=True)
    with open(outfile, 'w') as f:
        json.dump(selected, f, indent=2)
    print(f"\n  Saved: {outfile}")
    
    # Optionally run immediately
    if '--run' in sys.argv:
        import subprocess
        cmd = f'python benchmark_homology.py run --queries={query_list}'
        print(f"\n  Running: {cmd}")
        subprocess.run(cmd, shell=True, cwd=str(SCRIPT_DIR))


if __name__ == "__main__":
    main()

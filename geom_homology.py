#!/usr/bin/env python3
"""
GEOMETRIC HOMOLOGY SEARCH — Find structural analogs through curvature barcodes.

Compares protein curvature barcodes on T² to find geometric analogs invisible
to sequence alignment. Two-tier approach:
  1. Winding filter: reject candidates with incompatible (p, q) topology
  2. Segment alignment: Smith-Waterman on curvature segment strings with
     a custom substitution matrix over (ss_type, model_class, parameters)

Usage:
    python geom_homology.py search P00533            # Search database for analogs
    python geom_homology.py search P00533 --top 20   # Top 20 hits
    python geom_homology.py compare P00533 P68871    # Pairwise comparison
    python geom_homology.py matrix                   # Show substitution matrix
    python geom_homology.py cluster                  # All-vs-all clustering
    python geom_homology.py cluster --dendrogram     # With ASCII dendrogram

Requires: numpy. Optional: scipy (for clustering).
"""

import json
import math
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np

# ─── Paths ───
SCRIPT_DIR = Path(__file__).parent
BARCODE_DIR = SCRIPT_DIR / "results" / "barcodes"
RESULTS_DIR = SCRIPT_DIR / "results"


# ═══════════════════════════════════════════════════════════════════════
#  SCORING SYSTEM
# ═══════════════════════════════════════════════════════════════════════

# Macro-categories for model classes
MACRO = {
    'geodesic': 'constant_k', 'circular_arc': 'constant_k',
    'clothoid': 'monotone', 'exponential': 'monotone', 'fermat': 'monotone',
    'linear': 'monotone',
    'quadratic': 'polynomial',
    'gauss_peak': 'localized',
    'sinusoidal': 'oscillatory', 'damped_oscillation': 'oscillatory',
    'damped_osc': 'oscillatory',
    'sigmoid': 'transition', 'step': 'transition',
}

# All possible macro categories
MACROS = ['constant_k', 'monotone', 'polynomial', 'localized', 'oscillatory', 'transition']

# SS types
SS_TYPES = ['H', 'E', 'C']

# ── Substitution matrix ──
# Score for matching two segments based on (ss_type, macro_class)
# Same SS + same macro = best; same SS + related macro = ok; different SS = penalty

# Macro-category similarity (symmetric)
_MACRO_SIM = {}
for m in MACROS:
    _MACRO_SIM[(m, m)] = 1.0

# Related macros get partial credit
_MACRO_SIM[('constant_k', 'monotone')] = 0.3   # arc → gentle trend
_MACRO_SIM[('monotone', 'polynomial')] = 0.5    # linear → quadratic
_MACRO_SIM[('oscillatory', 'localized')] = 0.3  # wave → single peak
_MACRO_SIM[('transition', 'localized')] = 0.3   # step → peak
_MACRO_SIM[('transition', 'monotone')] = 0.4    # step → trend
# Symmetrize
for (a, b), v in list(_MACRO_SIM.items()):
    _MACRO_SIM[(b, a)] = v

# Fill defaults
for m1 in MACROS:
    for m2 in MACROS:
        if (m1, m2) not in _MACRO_SIM:
            _MACRO_SIM[(m1, m2)] = 0.0


def segment_score(seg1, seg2, length_weight=0.3, param_weight=0.3):
    """Score the similarity of two barcode segments.
    
    Components:
      1. SS match (hard): same SS required for positive score
      2. Macro-class similarity (soft): from substitution matrix
      3. Length similarity (soft): penalize large length differences
      4. Parameter similarity (soft): for Gaussian peaks, compare (A, σ)
    
    Returns score in [-1, +1] range approximately.
    """
    ss1 = seg1.get('ss_type', '?')
    ss2 = seg2.get('ss_type', '?')
    
    mc1 = MACRO.get(seg1.get('model_class', ''), 'constant_k')
    mc2 = MACRO.get(seg2.get('model_class', ''), 'constant_k')
    
    # SS match — different SS is a strong negative
    if ss1 != ss2:
        return -0.5
    
    # Macro-class similarity
    macro_sim = _MACRO_SIM.get((mc1, mc2), 0.0)
    
    # Length similarity: exp(-|Δlen|/scale)
    len1 = seg1.get('length', 10)
    len2 = seg2.get('length', 10)
    length_sim = math.exp(-abs(len1 - len2) / max(max(len1, len2), 1))
    
    # Parameter similarity for Gaussian peaks
    param_sim = 0.5  # default neutral
    p1 = seg1.get('params', {})
    p2 = seg2.get('params', {})
    
    if mc1 == 'localized' and mc2 == 'localized':
        # Compare amplitude and width
        A1 = abs(p1.get('A', 0))
        A2 = abs(p2.get('A', 0))
        s1 = p1.get('sigma', 1)
        s2 = p2.get('sigma', 1)
        if A1 > 0 and A2 > 0:
            A_ratio = min(A1, A2) / max(A1, A2)
            s_ratio = min(s1, s2) / max(s1, s2) if max(s1, s2) > 0 else 0
            param_sim = 0.5 * (A_ratio + s_ratio)
    elif mc1 == mc2 == 'constant_k':
        # Compare mean curvature
        k1 = abs(seg1.get('mean_kappa', 0))
        k2 = abs(seg2.get('mean_kappa', 0))
        if max(k1, k2) > 0:
            param_sim = min(k1, k2) / max(k1, k2)
        else:
            param_sim = 1.0
    
    # Weighted combination
    base = macro_sim
    score = base + length_weight * length_sim + param_weight * param_sim
    
    # Normalize to roughly [-1, 1]
    return score


# ── Gap penalties ──
GAP_OPEN = -0.8
GAP_EXTEND = -0.2


def align_barcodes(segs1, segs2):
    """Smith-Waterman local alignment of two segment arrays.
    
    Returns:
        score: alignment score (higher = more similar)
        alignment: list of (i, j) matched pairs (-1 for gaps)
        n_matched: number of matched segment pairs
    """
    n = len(segs1)
    m = len(segs2)
    
    if n == 0 or m == 0:
        return 0.0, [], 0
    
    # Score matrix
    H = np.zeros((n + 1, m + 1))
    traceback = np.zeros((n + 1, m + 1), dtype=int)  # 0=stop, 1=diag, 2=up, 3=left
    
    max_score = 0.0
    max_pos = (0, 0)
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Match/mismatch
            match = H[i-1, j-1] + segment_score(segs1[i-1], segs2[j-1])
            
            # Gap in segs2 (deletion)
            gap_up = H[i-1, j] + (GAP_EXTEND if traceback[i-1, j] == 2 else GAP_OPEN)
            
            # Gap in segs1 (insertion)
            gap_left = H[i, j-1] + (GAP_EXTEND if traceback[i, j-1] == 3 else GAP_OPEN)
            
            best = max(0, match, gap_up, gap_left)
            H[i, j] = best
            
            if best == 0:
                traceback[i, j] = 0
            elif best == match:
                traceback[i, j] = 1
            elif best == gap_up:
                traceback[i, j] = 2
            else:
                traceback[i, j] = 3
            
            if best > max_score:
                max_score = best
                max_pos = (i, j)
    
    # Traceback
    alignment = []
    i, j = max_pos
    while i > 0 and j > 0 and H[i, j] > 0:
        if traceback[i, j] == 1:  # diagonal
            alignment.append((i-1, j-1))
            i -= 1
            j -= 1
        elif traceback[i, j] == 2:  # up (gap in segs2)
            alignment.append((i-1, -1))
            i -= 1
        elif traceback[i, j] == 3:  # left (gap in segs1)
            alignment.append((-1, j-1))
            j -= 1
        else:
            break
    
    alignment.reverse()
    n_matched = sum(1 for a, b in alignment if a >= 0 and b >= 0)
    
    return max_score, alignment, n_matched


# ═══════════════════════════════════════════════════════════════════════
#  WINDING VECTOR SIMILARITY
# ═══════════════════════════════════════════════════════════════════════

def winding_similarity(bc1, bc2):
    """Compare winding vectors (p, q) of two barcodes.
    
    Returns similarity in [0, 1]:
      1.0 = identical winding
      0.0 = orthogonal or very different
    
    Uses cosine similarity of normalized winding vectors plus
    magnitude similarity.
    """
    p1, q1 = bc1['p'], bc1['q']
    p2, q2 = bc2['p'], bc2['q']
    Q1 = bc1['Q_magnitude']
    Q2 = bc2['Q_magnitude']
    
    # Handle zero vectors
    if Q1 < 0.1 and Q2 < 0.1:
        return 1.0  # both near-zero winding
    if Q1 < 0.1 or Q2 < 0.1:
        return 0.3  # one is near-zero
    
    # Cosine similarity of (p, q) vectors
    dot = p1 * p2 + q1 * q2
    cos_sim = dot / (Q1 * Q2)
    cos_sim = max(-1, min(1, cos_sim))  # clamp
    
    # Map from [-1, 1] to [0, 1]
    cos_score = (cos_sim + 1) / 2
    
    # Magnitude similarity: exp(-|ΔQ|/scale)
    mag_sim = math.exp(-abs(Q1 - Q2) / max(Q1, Q2))
    
    # Combined: weight cosine more heavily
    return 0.7 * cos_score + 0.3 * mag_sim


def length_similarity(bc1, bc2):
    """Compare chain lengths. Returns similarity in [0, 1]."""
    l1 = bc1.get('length', 100)
    l2 = bc2.get('length', 100)
    return min(l1, l2) / max(l1, l2) if max(l1, l2) > 0 else 1.0


def ss_composition_similarity(bc1, bc2):
    """Compare secondary structure compositions. Returns similarity in [0, 1]."""
    ss1 = bc1.get('ss_composition', {})
    ss2 = bc2.get('ss_composition', {})
    
    t1 = sum(ss1.values()) or 1
    t2 = sum(ss2.values()) or 1
    
    # Euclidean distance of SS fractions
    dist = 0
    for ss in ['H', 'E', 'C']:
        f1 = ss1.get(ss, 0) / t1
        f2 = ss2.get(ss, 0) / t2
        dist += (f1 - f2) ** 2
    dist = math.sqrt(dist)
    
    # Max distance is sqrt(2) for completely different compositions
    return 1.0 - dist / math.sqrt(2)


# ═══════════════════════════════════════════════════════════════════════
#  COMPOSITE SIMILARITY
# ═══════════════════════════════════════════════════════════════════════

def geometric_similarity(bc1, bc2, align=True):
    """Compute composite geometric similarity between two barcodes.
    
    Components (weighted):
      - Winding similarity (20%): topological compatibility
      - SS composition (10%): overall structural class match
      - Length similarity (5%): size compatibility
      - Segment alignment (35%): detailed curvature profile match
      - Coverage (15%): fraction of both proteins spanned by alignment
      - Depth bonus (15%): rewards many matched segments (kills single-segment hits)
    
    Returns dict with total score and component scores.
    """
    w_wind = winding_similarity(bc1, bc2)
    w_ss = ss_composition_similarity(bc1, bc2)
    w_len = length_similarity(bc1, bc2)
    
    # Segment alignment (expensive — skip if winding is very different)
    if align and w_wind > 0.2:
        segs1 = bc1.get('segments', [])
        segs2 = bc2.get('segments', [])
        raw_score, alignment, n_matched = align_barcodes(segs1, segs2)
        
        # Normalize alignment score by max possible
        max_possible = min(len(segs1), len(segs2)) * 1.6  # rough max per-segment score
        w_align = raw_score / max_possible if max_possible > 0 else 0
        w_align = max(0, min(1, w_align))
        
        # Coverage: fraction of each protein's segments covered by the alignment
        n1 = len(segs1) or 1
        n2 = len(segs2) or 1
        cov1 = n_matched / n1  # fraction of query covered
        cov2 = n_matched / n2  # fraction of hit covered
        w_coverage = math.sqrt(cov1 * cov2)  # geometric mean
        
        # Depth bonus: sigmoid on n_matched, saturates around 20
        # 1 match → 0.05, 5 matches → 0.26, 10 → 0.50, 20 → 0.82, 50 → 0.98
        w_depth = 1.0 / (1.0 + math.exp(-0.2 * (n_matched - 10)))
    else:
        raw_score = 0
        alignment = []
        n_matched = 0
        w_align = 0
        w_coverage = 0
        w_depth = 0
    
    # Composite — depth and coverage prevent single-segment inflation
    total = (0.20 * w_wind + 0.10 * w_ss + 0.05 * w_len +
             0.35 * w_align + 0.15 * w_coverage + 0.15 * w_depth)
    
    return {
        'total': total,
        'winding': w_wind,
        'ss_comp': w_ss,
        'length': w_len,
        'alignment': w_align,
        'coverage': w_coverage,
        'depth': w_depth,
        'raw_align_score': raw_score,
        'n_matched': n_matched,
        'n_aligned': len(alignment),
    }


# ═══════════════════════════════════════════════════════════════════════
#  LOADING
# ═══════════════════════════════════════════════════════════════════════

def load_all_barcodes():
    """Load all barcode JSON files."""
    barcodes = []
    if not BARCODE_DIR.exists():
        return barcodes
    for f in sorted(BARCODE_DIR.glob("*_barcode.json")):
        try:
            with open(f) as fh:
                bc = json.load(fh)
                bc['_file'] = f.name
                barcodes.append(bc)
        except (json.JSONDecodeError, KeyError):
            pass
    return barcodes


def find_barcode(barcodes, query):
    """Find a barcode by UniProt ID prefix."""
    query = query.strip().upper()
    for bc in barcodes:
        if query in bc.get('uniprot_id', '').upper():
            return bc
    return None


# ═══════════════════════════════════════════════════════════════════════
#  COMMANDS
# ═══════════════════════════════════════════════════════════════════════

def cmd_search(barcodes, query_id, top_k=10):
    """Search the barcode database for geometric analogs of a query protein."""
    query = find_barcode(barcodes, query_id)
    if query is None:
        print(f"  Protein {query_id} not found in barcode database.")
        return
    
    quid = query.get('uniprot_id', query_id)
    print(f"\n  GEOMETRIC HOMOLOGY SEARCH")
    print(f"  Query: {quid} ({query['length']} aa, {query['n_segments']} segments, "
          f"|Q|={query['Q_magnitude']:.2f})")
    print(f"  SS: H={query['ss_composition']['H']}, E={query['ss_composition']['E']}, "
          f"C={query['ss_composition']['C']}")
    print(f"  {'='*72}")
    
    # Score against all others
    hits = []
    for bc in barcodes:
        uid = bc.get('uniprot_id', '')
        if uid == quid:
            continue
        
        # Quick winding filter
        w_sim = winding_similarity(query, bc)
        if w_sim < 0.15:
            continue
        
        sim = geometric_similarity(query, bc, align=True)
        hits.append((sim['total'], uid, bc, sim))
    
    hits.sort(key=lambda x: -x[0])
    
    print(f"\n  Searched {len(barcodes) - 1} proteins, {len(hits)} passed winding filter")
    print(f"\n  {'Rank':>4s} {'UniProt':>12s} {'Score':>7s} {'Wind':>6s} {'SS':>6s} "
          f"{'Align':>6s} {'Cov':>5s} {'Depth':>5s} {'Matched':>8s} {'Length':>6s}")
    print(f"  {'-'*4} {'-'*12} {'-'*7} {'-'*6} {'-'*6} {'-'*6} {'-'*5} {'-'*5} {'-'*8} {'-'*6}")
    
    for rank, (score, uid, bc, sim) in enumerate(hits[:top_k], 1):
        print(f"  {rank:4d} {uid:>12s} {score:7.3f} {sim['winding']:6.3f} "
              f"{sim['ss_comp']:6.3f} {sim['alignment']:6.3f} "
              f"{sim.get('coverage', 0):5.3f} {sim.get('depth', 0):5.2f} "
              f"{sim['n_matched']:4d}/{sim['n_aligned']:<3d} "
              f"{bc['length']:6d}")
    
    # Show detailed alignment for top hit
    if hits:
        best_score, best_uid, best_bc, best_sim = hits[0]
        print(f"\n  TOP HIT DETAIL: {quid} vs {best_uid}")
        print(f"  {'-'*60}")
        
        segs1 = query.get('segments', [])
        segs2 = best_bc.get('segments', [])
        _, alignment, _ = align_barcodes(segs1, segs2)
        
        print(f"  {'':>4s} {'Query':>20s} {'':>3s} {'Hit':>20s} {'Score':>7s}")
        print(f"  {'':>4s} {'-'*20} {'':>3s} {'-'*20} {'-'*7}")
        
        for idx, (i, j) in enumerate(alignment[:30]):
            if i >= 0 and j >= 0:
                s1 = segs1[i]
                s2 = segs2[j]
                sc = segment_score(s1, s2)
                q_str = f"{s1['ss_type']}/{s1['model_class'][:12]:12s}({s1['length']:2d})"
                h_str = f"{s2['ss_type']}/{s2['model_class'][:12]:12s}({s2['length']:2d})"
                marker = "✓" if sc > 0.5 else " "
                print(f"  {idx+1:4d} {q_str:>20s} <-> {h_str:>20s} {sc:+6.2f} {marker}")
            elif i >= 0:
                s1 = segs1[i]
                q_str = f"{s1['ss_type']}/{s1['model_class'][:12]:12s}({s1['length']:2d})"
                print(f"  {idx+1:4d} {q_str:>20s} <-> {'--- gap ---':>20s}")
            else:
                s2 = segs2[j]
                h_str = f"{s2['ss_type']}/{s2['model_class'][:12]:12s}({s2['length']:2d})"
                print(f"  {idx+1:4d} {'--- gap ---':>20s} <-> {h_str:>20s}")
        
        if len(alignment) > 30:
            print(f"  ... ({len(alignment) - 30} more aligned segments)")
    
    # Save
    out = {
        'query': quid,
        'n_searched': len(barcodes) - 1,
        'n_passed_filter': len(hits),
        'top_hits': [{'uid': uid, 'score': score, **sim}
                     for score, uid, bc, sim in hits[:top_k]]
    }
    out_path = RESULTS_DIR / f"homology_{quid}.json"
    try:
        with open(out_path, 'w') as f:
            json.dump(out, f, indent=2)
        print(f"\n  Results saved: {out_path}")
    except Exception:
        pass
    print()


def cmd_compare(barcodes, id1, id2):
    """Detailed pairwise comparison of two proteins."""
    bc1 = find_barcode(barcodes, id1)
    bc2 = find_barcode(barcodes, id2)
    
    if bc1 is None:
        print(f"  {id1} not found."); return
    if bc2 is None:
        print(f"  {id2} not found."); return
    
    uid1 = bc1.get('uniprot_id', id1)
    uid2 = bc2.get('uniprot_id', id2)
    
    print(f"\n  PAIRWISE GEOMETRIC COMPARISON")
    print(f"  {'='*60}")
    print(f"  Protein A: {uid1} ({bc1['length']} aa, {bc1['n_segments']} segs, |Q|={bc1['Q_magnitude']:.2f})")
    print(f"  Protein B: {uid2} ({bc2['length']} aa, {bc2['n_segments']} segs, |Q|={bc2['Q_magnitude']:.2f})")
    
    sim = geometric_similarity(bc1, bc2, align=True)
    
    print(f"\n  SIMILARITY SCORES:")
    print(f"    Total:     {sim['total']:.4f}")
    print(f"    Winding:   {sim['winding']:.4f} (topology)")
    print(f"    SS comp:   {sim['ss_comp']:.4f} (structural class)")
    print(f"    Length:     {sim['length']:.4f} (size)")
    print(f"    Alignment: {sim['alignment']:.4f} (curvature profile)")
    print(f"    Matched:   {sim['n_matched']}/{sim['n_aligned']} segments")
    
    # Winding detail
    print(f"\n  WINDING COMPARISON:")
    print(f"    {'':15s} {'Protein A':>12s} {'Protein B':>12s} {'Δ':>10s}")
    print(f"    {'p (φ-winding)':15s} {bc1['p']:+12.3f} {bc2['p']:+12.3f} {bc1['p']-bc2['p']:+10.3f}")
    print(f"    {'q (ψ-winding)':15s} {bc1['q']:+12.3f} {bc2['q']:+12.3f} {bc1['q']-bc2['q']:+10.3f}")
    print(f"    {'|Q|':15s} {bc1['Q_magnitude']:12.3f} {bc2['Q_magnitude']:12.3f} "
          f"{bc1['Q_magnitude']-bc2['Q_magnitude']:+10.3f}")
    
    # SS composition
    print(f"\n  SS COMPOSITION:")
    for ss in ['H', 'E', 'C']:
        v1 = bc1['ss_composition'].get(ss, 0)
        v2 = bc2['ss_composition'].get(ss, 0)
        t1 = sum(bc1['ss_composition'].values()) or 1
        t2 = sum(bc2['ss_composition'].values()) or 1
        print(f"    {ss}: {v1/t1:.1%} vs {v2/t2:.1%}")
    
    # Segment alignment
    segs1 = bc1.get('segments', [])
    segs2 = bc2.get('segments', [])
    _, alignment, _ = align_barcodes(segs1, segs2)
    
    if alignment:
        print(f"\n  SEGMENT ALIGNMENT ({len(alignment)} positions):")
        print(f"  {'':>4s} {'Protein A':>22s} {'':>3s} {'Protein B':>22s} {'Score':>7s}")
        
        for idx, (i, j) in enumerate(alignment):
            if i >= 0 and j >= 0:
                s1, s2 = segs1[i], segs2[j]
                sc = segment_score(s1, s2)
                a_str = f"{s1['ss_type']}/{s1['model_class'][:12]:12s}({s1['length']:2d})"
                b_str = f"{s2['ss_type']}/{s2['model_class'][:12]:12s}({s2['length']:2d})"
                print(f"  {idx+1:4d} {a_str:>22s} <-> {b_str:>22s} {sc:+6.2f}")
            elif i >= 0:
                s1 = segs1[i]
                a_str = f"{s1['ss_type']}/{s1['model_class'][:12]:12s}({s1['length']:2d})"
                print(f"  {idx+1:4d} {a_str:>22s} <-> {'--- gap ---':>22s}")
            else:
                s2 = segs2[j]
                b_str = f"{s2['ss_type']}/{s2['model_class'][:12]:12s}({s2['length']:2d})"
                print(f"  {idx+1:4d} {'--- gap ---':>22s} <-> {b_str:>22s}")
    
    print()


def cmd_cluster(barcodes, dendrogram=False):
    """All-vs-all geometric similarity clustering."""
    n = len(barcodes)
    print(f"\n  ALL-VS-ALL GEOMETRIC CLUSTERING")
    print(f"  {'='*60}")
    print(f"  Proteins: {n}")
    
    if n > 200:
        print(f"  Limiting to first 200 for performance.")
        barcodes = barcodes[:200]
        n = 200
    
    if n < 5:
        print("  Need at least 5 proteins.")
        return
    
    # Compute pairwise distance matrix (using 1 - similarity)
    # Use fast mode: no alignment, just winding + SS + length
    print(f"  Computing pairwise similarities (fast mode: no alignment)...")
    dist = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i + 1, n):
            w = winding_similarity(barcodes[i], barcodes[j])
            ss = ss_composition_similarity(barcodes[i], barcodes[j])
            le = length_similarity(barcodes[i], barcodes[j])
            sim = 0.5 * w + 0.3 * ss + 0.2 * le
            dist[i, j] = 1 - sim
            dist[j, i] = dist[i, j]
    
    # Find natural clusters using simple k-medoids or hierarchical
    # Start with finding most similar pairs
    print(f"\n  MOST SIMILAR PAIRS (geometric analogs):")
    print(f"  {'Protein A':>12s} {'Protein B':>12s} {'Similarity':>11s} {'Wind':>6s} {'SS':>6s}")
    print(f"  {'-'*12} {'-'*12} {'-'*11} {'-'*6} {'-'*6}")
    
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            sim = 1 - dist[i, j]
            pairs.append((sim, i, j))
    pairs.sort(key=lambda x: -x[0])
    
    for sim, i, j in pairs[:20]:
        uid_i = barcodes[i].get('uniprot_id', '?')
        uid_j = barcodes[j].get('uniprot_id', '?')
        w = winding_similarity(barcodes[i], barcodes[j])
        ss = ss_composition_similarity(barcodes[i], barcodes[j])
        print(f"  {uid_i:>12s} {uid_j:>12s} {sim:11.4f} {w:6.3f} {ss:6.3f}")
    
    # Most isolated proteins (least similar to anything)
    print(f"\n  MOST GEOMETRICALLY UNIQUE PROTEINS:")
    max_sim_per = []
    for i in range(n):
        best = max(1 - dist[i, j] for j in range(n) if j != i)
        max_sim_per.append((best, i))
    max_sim_per.sort()
    
    for best_sim, i in max_sim_per[:10]:
        uid = barcodes[i].get('uniprot_id', '?')
        print(f"    {uid:>12s}: max similarity to any other = {best_sim:.4f} "
              f"(len={barcodes[i]['length']}, |Q|={barcodes[i]['Q_magnitude']:.1f})")
    
    # Simple cluster assignment: group by SS composition class
    print(f"\n  STRUCTURAL CLASS GROUPS:")
    groups = defaultdict(list)
    for i, bc in enumerate(barcodes):
        ss = bc.get('ss_composition', {})
        total = sum(ss.values()) or 1
        h = ss.get('H', 0) / total
        e = ss.get('E', 0) / total
        if h > 0.5:
            groups['helix-rich'].append(i)
        elif e > 0.3:
            groups['strand-rich'].append(i)
        elif h + e < 0.3:
            groups['coil-rich'].append(i)
        else:
            groups['mixed'].append(i)
    
    for name, members in sorted(groups.items()):
        if not members:
            continue
        # Mean intra-group similarity
        if len(members) >= 2:
            intra_sims = []
            for a in range(len(members)):
                for b in range(a + 1, len(members)):
                    intra_sims.append(1 - dist[members[a], members[b]])
            mean_intra = sum(intra_sims) / len(intra_sims)
        else:
            mean_intra = 0
        print(f"    {name:15s}: n={len(members):4d}, mean intra-similarity={mean_intra:.4f}")
    
    # Save distance matrix
    out_path = RESULTS_DIR / "geometric_distances.json"
    try:
        uids = [bc.get('uniprot_id', '?') for bc in barcodes]
        with open(out_path, 'w') as f:
            json.dump({
                'n': n,
                'uids': uids,
                'top_pairs': [{'a': barcodes[i].get('uniprot_id'), 
                               'b': barcodes[j].get('uniprot_id'),
                               'similarity': round(sim, 4)}
                              for sim, i, j in pairs[:50]],
            }, f, indent=2)
        print(f"\n  Results saved: {out_path}")
    except Exception:
        pass
    print()


def cmd_matrix():
    """Display the substitution matrix."""
    print(f"\n  CURVATURE SEGMENT SUBSTITUTION MATRIX")
    print(f"  {'='*60}")
    print(f"  Macro-category similarity scores:")
    print()
    
    header = f"  {'':15s}" + "".join(f"{m[:8]:>9s}" for m in MACROS)
    print(header)
    print(f"  {'-'*15}" + "-" * 9 * len(MACROS))
    for m1 in MACROS:
        row = f"  {m1:15s}"
        for m2 in MACROS:
            v = _MACRO_SIM.get((m1, m2), 0)
            row += f"{v:9.2f}"
        print(row)
    
    print(f"\n  Additional scoring components:")
    print(f"    SS match required for positive score (different SS = -0.5)")
    print(f"    Length similarity: exp(-|ΔL|/max(L1,L2)), weight=0.3")
    print(f"    Parameter similarity (Gaussian peaks): min(A)/max(A) × min(σ)/max(σ), weight=0.3")
    print(f"    Gap open: {GAP_OPEN}, Gap extend: {GAP_EXTEND}")
    print()


# ─── Main ───
if __name__ == "__main__":
    args = sys.argv[1:]
    
    if not args:
        print(__doc__)
        sys.exit(0)
    
    cmd = args[0].lower()
    barcodes = load_all_barcodes()
    
    if not barcodes and cmd != 'matrix':
        print(f"  No barcodes found in {BARCODE_DIR}")
        print("  Run the main pipeline first to generate barcodes.")
        sys.exit(1)
    
    print(f"\n  Loaded {len(barcodes)} barcodes")
    
    if cmd == 'search':
        query = args[1] if len(args) > 1 else None
        if not query:
            print("  Usage: geom_homology.py search <UniProt_ID>")
            sys.exit(1)
        top_k = 10
        for i, a in enumerate(args):
            if a == '--top' and i + 1 < len(args):
                top_k = int(args[i + 1])
        cmd_search(barcodes, query, top_k=top_k)
    
    elif cmd == 'compare':
        if len(args) < 3:
            print("  Usage: geom_homology.py compare <ID1> <ID2>")
            sys.exit(1)
        cmd_compare(barcodes, args[1], args[2])
    
    elif cmd == 'matrix':
        cmd_matrix()
    
    elif cmd == 'cluster':
        dendrogram = '--dendrogram' in args
        cmd_cluster(barcodes, dendrogram=dendrogram)
    
    else:
        print(f"  Unknown command: {cmd}")
        print("  Commands: search, compare, matrix, cluster")

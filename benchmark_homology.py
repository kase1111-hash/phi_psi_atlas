#!/usr/bin/env python3
"""
GEOMETRIC HOMOLOGY BENCHMARK — Validate curvature-based search against SCOP/CATH/BLAST.

Runs a diverse panel of query proteins through the geometric homology search,
fetches structural classifications, and computes precision metrics.

Usage:
    python benchmark_homology.py run                # Full benchmark (fetches SCOP/CATH)
    python benchmark_homology.py run --no-fetch      # Offline (uses cached classifications)
    python benchmark_homology.py run --queries P00533,P69905,P01308  # Custom queries
    python benchmark_homology.py classify            # Just fetch/show classifications
    python benchmark_homology.py report              # Generate report from cached results

Requires: numpy. Network access for SCOP/CATH/BLAST fetching (or --no-fetch for offline).
"""

import json
import math
import os
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path
from collections import defaultdict, Counter

import numpy as np

# ─── Paths ───
SCRIPT_DIR = Path(__file__).parent
BARCODE_DIR = SCRIPT_DIR / "results" / "barcodes"
RESULTS_DIR = SCRIPT_DIR / "results"
CACHE_DIR = SCRIPT_DIR / "data" / "classifications"
BENCHMARK_OUT = RESULTS_DIR / "benchmark_results.json"

# Import from geom_homology
sys.path.insert(0, str(SCRIPT_DIR))
try:
    from geom_homology import (
        load_all_barcodes, find_barcode, geometric_similarity,
        winding_similarity, align_barcodes
    )
except ImportError:
    print("  [ERROR] geom_homology.py must be in the same directory.")
    sys.exit(1)


# ═══════════════════════════════════════════════════════════════════════
#  DEFAULT QUERY PANEL — diverse structural classes
# ═══════════════════════════════════════════════════════════════════════

# These are well-characterized human proteins spanning structural space.
# The benchmark will use whichever are present in the barcode database.
DEFAULT_QUERIES = {
    # All-alpha
    "P69905": "Hemoglobin alpha (all-alpha, globin fold)",
    "P68871": "Hemoglobin beta (all-alpha, globin fold)",
    "P02144": "Myoglobin (all-alpha, globin fold)",
    "P00918": "Carbonic anhydrase 2 (alpha/beta, CA fold)",
    # Receptor tyrosine kinases
    "P00533": "EGFR (multi-domain, RTK)",
    "P04626": "ErbB2/HER2 (multi-domain, RTK)",
    "P21860": "ErbB3 (multi-domain, RTK)",
    # Beta-rich
    "P01857": "IgG1 Fc (all-beta, Ig fold)",
    "P01834": "Ig kappa chain C (all-beta, Ig fold)",
    "P02751": "Fibronectin (beta-rich, FN3 repeats)",
    # Alpha/beta
    "P00338": "L-lactate dehydrogenase A (alpha/beta, Rossmann)",
    "P04406": "GAPDH (alpha/beta, Rossmann)",
    "P60174": "Triosephosphate isomerase (alpha/beta, TIM barrel)",
    # Small / mixed
    "P01308": "Insulin (small, alpha+beta)",
    "P62258": "14-3-3 epsilon (all-alpha, 14-3-3 fold)",
    # Membrane
    "P11166": "GLUT1 (membrane transporter)",
    # Disordered / coil-rich
    "P06748": "Nucleophosmin (mixed, partly disordered)",
    "P10636": "Tau (largely disordered)",
    # Enzymes
    "P07900": "HSP90-alpha (alpha-rich, chaperone)",
    "P38646": "Mortalin/GRP75 (alpha/beta, HSP70 fold)",
}


# ═══════════════════════════════════════════════════════════════════════
#  CLASSIFICATION FETCHING
# ═══════════════════════════════════════════════════════════════════════

def _fetch_url(url, timeout=10):
    """Fetch URL with timeout, return text or None."""
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'PhiPsiAtlas/1.0'})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read().decode('utf-8', errors='replace')
    except Exception:
        return None


def fetch_interpro_classification(uniprot_id):
    """Fetch structural classification from InterPro/UniProt API.
    
    Returns dict with: interpro_families, scop_superfamily, cath_topology,
    pfam_families, go_terms.
    """
    uid = uniprot_id.strip().upper()
    
    result = {
        'uniprot_id': uid,
        'protein_name': '',
        'scop2_superfamily': '',
        'cath_topology': '',
        'pfam_families': [],
        'interpro_entries': [],
        'gene_name': '',
    }
    
    # 1. UniProt basic info
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
    text = _fetch_url(url, timeout=15)
    if text:
        try:
            data = json.loads(text)
            # Protein name
            rec_name = data.get('proteinDescription', {}).get('recommendedName', {})
            if rec_name:
                full_name = rec_name.get('fullName', {}).get('value', '')
                result['protein_name'] = full_name
            # Gene name
            genes = data.get('genes', [])
            if genes:
                result['gene_name'] = genes[0].get('geneName', {}).get('value', '')
            
            # Cross-references for SCOP, CATH, Pfam
            xrefs = data.get('uniProtKBCrossReferences', [])
            for xref in xrefs:
                db = xref.get('database', '')
                xid = xref.get('id', '')
                props = {p.get('key', ''): p.get('value', '') for p in xref.get('properties', [])}
                
                if db == 'Pfam':
                    result['pfam_families'].append({
                        'id': xid,
                        'name': props.get('EntryName', ''),
                    })
                elif db == 'SUPFAM':
                    result['scop2_superfamily'] = xid
                elif db == 'Gene3D':
                    if not result['cath_topology']:
                        result['cath_topology'] = xid
                elif db == 'InterPro':
                    result['interpro_entries'].append({
                        'id': xid,
                        'name': props.get('EntryName', ''),
                    })
        except (json.JSONDecodeError, KeyError):
            pass
    
    return result


def fetch_blast_identity(uid1, uid2):
    """Estimate sequence identity between two proteins using UniProt alignment.
    
    This is a rough estimate — fetches sequences and does a simple identity count.
    For real benchmarking, use actual BLAST.
    """
    # Fetch both sequences
    seqs = {}
    for uid in [uid1, uid2]:
        url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
        text = _fetch_url(url, timeout=10)
        if text:
            lines = text.strip().split('\n')
            seqs[uid] = ''.join(l.strip() for l in lines[1:] if not l.startswith('>'))
    
    if len(seqs) < 2:
        return None
    
    # Simple: shared k-mer fraction as a proxy for sequence identity
    # (Real BLAST would be better, but this works offline)
    s1 = seqs[uid1]
    s2 = seqs[uid2]
    
    k = 3
    kmers1 = set(s1[i:i+k] for i in range(len(s1) - k + 1))
    kmers2 = set(s2[i:i+k] for i in range(len(s2) - k + 1))
    
    if not kmers1 or not kmers2:
        return None
    
    shared = len(kmers1 & kmers2)
    total = len(kmers1 | kmers2)
    jaccard = shared / total if total > 0 else 0
    
    return jaccard


def load_or_fetch_classifications(uids, no_fetch=False):
    """Load cached classifications or fetch from web."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    
    classifications = {}
    to_fetch = []
    
    for uid in uids:
        cache_file = CACHE_DIR / f"{uid}.json"
        if cache_file.exists():
            try:
                with open(cache_file) as f:
                    classifications[uid] = json.load(f)
                continue
            except (json.JSONDecodeError, IOError):
                pass
        to_fetch.append(uid)
    
    if to_fetch and not no_fetch:
        print(f"  Fetching classifications for {len(to_fetch)} proteins...")
        for i, uid in enumerate(to_fetch):
            try:
                cls = fetch_interpro_classification(uid)
                classifications[uid] = cls
                # Cache
                cache_file = CACHE_DIR / f"{uid}.json"
                with open(cache_file, 'w') as f:
                    json.dump(cls, f, indent=1)
            except Exception as e:
                print(f"    {uid}: fetch failed ({e})")
                classifications[uid] = {'uniprot_id': uid, 'protein_name': '?', 'pfam_families': []}
            
            if (i + 1) % 10 == 0:
                print(f"    {i+1}/{len(to_fetch)}")
            time.sleep(0.5)  # Rate limit
    elif to_fetch:
        for uid in to_fetch:
            classifications[uid] = {'uniprot_id': uid, 'protein_name': '?', 'pfam_families': []}
    
    return classifications


# ═══════════════════════════════════════════════════════════════════════
#  CLASSIFICATION COMPARISON
# ═══════════════════════════════════════════════════════════════════════

def pfam_overlap(cls1, cls2):
    """Compute Pfam domain overlap between two proteins.
    
    Returns: (n_shared, n_union, jaccard)
    """
    fams1 = set(f['id'] for f in cls1.get('pfam_families', []))
    fams2 = set(f['id'] for f in cls2.get('pfam_families', []))
    
    shared = len(fams1 & fams2)
    union = len(fams1 | fams2)
    jaccard = shared / union if union > 0 else 0
    
    return shared, union, jaccard


def same_scop_superfamily(cls1, cls2):
    """Check if two proteins share a SCOP superfamily."""
    s1 = cls1.get('scop2_superfamily', '')
    s2 = cls2.get('scop2_superfamily', '')
    if s1 and s2 and s1 == s2:
        return True
    return False


def same_cath_topology(cls1, cls2):
    """Check if two proteins share a CATH topology (first 3 levels)."""
    c1 = cls1.get('cath_topology', '')
    c2 = cls2.get('cath_topology', '')
    if not c1 or not c2:
        return None  # Unknown
    # CATH IDs like "1.10.510" — compare first 2 levels for topology
    parts1 = c1.split('.')[:3]
    parts2 = c2.split('.')[:3]
    return parts1 == parts2


# ═══════════════════════════════════════════════════════════════════════
#  BENCHMARK
# ═══════════════════════════════════════════════════════════════════════

def run_benchmark(barcodes, query_ids=None, top_k=20, no_fetch=False, fetch_blast=False):
    """Run full benchmark: search, classify, evaluate."""
    
    if query_ids is None:
        query_ids = list(DEFAULT_QUERIES.keys())
    
    # Filter to queries present in barcode database
    available = []
    for qid in query_ids:
        bc = find_barcode(barcodes, qid)
        if bc:
            available.append(qid)
    
    if not available:
        print("  No query proteins found in barcode database.")
        print(f"  Looked for: {query_ids[:5]}...")
        print(f"  Database has: {[bc.get('uniprot_id','?') for bc in barcodes[:5]]}...")
        return
    
    print(f"\n  GEOMETRIC HOMOLOGY BENCHMARK")
    print(f"  {'='*70}")
    print(f"  Query panel: {len(available)}/{len(query_ids)} available in database")
    print(f"  Database size: {len(barcodes)}")
    print(f"  Top-k: {top_k}")
    
    # ── Phase 1: Run searches ──
    print(f"\n  Phase 1: Running geometric searches...")
    search_results = {}
    
    for qid in available:
        query_bc = find_barcode(barcodes, qid)
        quid = query_bc.get('uniprot_id', qid)
        
        hits = []
        for bc in barcodes:
            uid = bc.get('uniprot_id', '')
            if uid == quid:
                continue
            w_sim = winding_similarity(query_bc, bc)
            if w_sim < 0.15:
                continue
            sim = geometric_similarity(query_bc, bc, align=True)
            hits.append((sim['total'], uid, sim))
        
        hits.sort(key=lambda x: -x[0])
        search_results[quid] = hits[:top_k]
        
        desc = DEFAULT_QUERIES.get(qid, '')
        n_segs = query_bc.get('n_segments', '?')
        print(f"    {quid}: {len(hits)} candidates, top score={hits[0][0]:.3f} "
              f"({n_segs} segs) — {desc}")
    
    # ── Phase 2: Collect all UIDs that need classification ──
    all_uids = set(available)
    for quid, hits in search_results.items():
        for score, uid, sim in hits:
            all_uids.add(uid)
    
    print(f"\n  Phase 2: Fetching classifications for {len(all_uids)} proteins...")
    classifications = load_or_fetch_classifications(list(all_uids), no_fetch=no_fetch)
    
    # ── Phase 3: Evaluate ──
    print(f"\n  Phase 3: Evaluating...")
    
    per_query = {}
    all_precisions = {1: [], 3: [], 5: [], 10: [], 20: []}
    all_pfam_overlaps = []
    geometry_only_hits = []  # High geom score, low Pfam overlap
    
    for quid in search_results:
        hits = search_results[quid]
        q_cls = classifications.get(quid, {})
        q_pfams = set(f['id'] for f in q_cls.get('pfam_families', []))
        
        hit_evaluations = []
        tp_at_k = {1: 0, 3: 0, 5: 0, 10: 0, 20: 0}
        scop_tp_at_k = {1: 0, 3: 0, 5: 0, 10: 0, 20: 0}
        
        for rank, (score, uid, sim) in enumerate(hits, 1):
            h_cls = classifications.get(uid, {})
            
            # Pfam overlap
            shared, union, jaccard = pfam_overlap(q_cls, h_cls)
            
            # SCOP/CATH match
            scop_match = same_scop_superfamily(q_cls, h_cls)
            cath_match = same_cath_topology(q_cls, h_cls)
            
            # Primary ground truth: SCOP superfamily match
            # Secondary: Pfam domain overlap
            is_scop_tp = scop_match is True
            is_pfam_tp = shared > 0
            is_tp = is_scop_tp or is_pfam_tp  # Either counts as TP
            
            eval_entry = {
                'rank': rank,
                'uid': uid,
                'geom_score': round(score, 4),
                'protein_name': h_cls.get('protein_name', '?')[:60],
                'gene_name': h_cls.get('gene_name', ''),
                'pfam_shared': shared,
                'pfam_jaccard': round(jaccard, 3),
                'scop_match': scop_match,
                'cath_match': cath_match,
                'is_tp': is_tp,
                'n_matched': sim.get('n_matched', 0),
                'coverage': round(sim.get('coverage', 0), 3),
            }
            hit_evaluations.append(eval_entry)
            
            all_pfam_overlaps.append(jaccard)
            
            # Track TP at various k
            for k in tp_at_k:
                if rank <= k and is_tp:
                    tp_at_k[k] += 1
                if rank <= k and is_scop_tp:
                    scop_tp_at_k[k] += 1
            
            # Geometry-only discovery: high score, zero Pfam overlap
            if score > 0.65 and shared == 0 and rank <= 10:
                geometry_only_hits.append({
                    'query': quid,
                    'query_name': q_cls.get('protein_name', '?')[:40],
                    'hit': uid,
                    'hit_name': h_cls.get('protein_name', '?')[:40],
                    'geom_score': round(score, 3),
                    'rank': rank,
                    'coverage': round(sim.get('coverage', 0), 3),
                    'n_matched': sim.get('n_matched', 0),
                })
        
        # Precision@k
        for k in all_precisions:
            p_at_k = tp_at_k[k] / min(k, len(hits)) if hits else 0
            all_precisions[k].append(p_at_k)
        
        per_query[quid] = {
            'description': DEFAULT_QUERIES.get(quid, ''),
            'n_hits': len(hits),
            'hits': hit_evaluations,
            'tp_at_k': tp_at_k,
            'scop_tp_at_k': scop_tp_at_k,
        }
    
    # ── Phase 4: Report ──
    print(f"\n  {'='*70}")
    print(f"  BENCHMARK RESULTS")
    print(f"  {'='*70}")
    
    # Precision@k summary
    print(f"\n  PRECISION@k (SCOP superfamily OR Pfam overlap = true positive):")
    print(f"  {'k':>5s} {'Mean P@k':>10s} {'Std':>8s} {'n queries':>10s}")
    for k in sorted(all_precisions.keys()):
        vals = all_precisions[k]
        if vals:
            mean_p = sum(vals) / len(vals)
            std_p = (sum((v - mean_p)**2 for v in vals) / len(vals)) ** 0.5
            print(f"  {k:5d} {mean_p:10.3f} {std_p:8.3f} {len(vals):10d}")
    
    # SCOP-only precision
    scop_precisions = {1: [], 3: [], 5: [], 10: [], 20: []}
    for quid in search_results:
        pq = per_query[quid]
        for k in scop_precisions:
            hits_k = min(k, pq['n_hits'])
            p_at_k = pq.get('scop_tp_at_k', {}).get(k, 0) / hits_k if hits_k else 0
            scop_precisions[k].append(p_at_k)
    
    print(f"\n  PRECISION@k (SCOP superfamily only):")
    print(f"  {'k':>5s} {'Mean P@k':>10s} {'n queries':>10s}")
    for k in sorted(scop_precisions.keys()):
        vals = scop_precisions[k]
        if vals:
            mean_p = sum(vals) / len(vals)
            print(f"  {k:5d} {mean_p:10.3f} {len(vals):10d}")
    
    # Per-query summary
    print(f"\n  PER-QUERY RESULTS:")
    print(f"  {'Query':>10s} {'Description':>40s} {'SCOP@5':>7s} {'TP@10':>6s} {'Top hit':>12s} {'SCOP?':>6s}")
    print(f"  {'-'*10} {'-'*40} {'-'*7} {'-'*6} {'-'*12} {'-'*6}")
    
    for quid in search_results:
        pq = per_query[quid]
        hits = pq['hits']
        top_uid = hits[0]['uid'] if hits else '?'
        top_tp = '✓' if hits and (hits[0].get('scop_match') or hits[0]['is_tp']) else '✗'
        desc = pq['description'][:40]
        print(f"  {quid:>10s} {desc:>40s} {pq['tp_at_k'][5]:5d} {pq['tp_at_k'][10]:6d} "
              f"{top_uid:>12s} {top_tp:>6s}")
    
    # Geometry-only discoveries
    if geometry_only_hits:
        print(f"\n  GEOMETRY-ONLY HITS (score > 0.65, zero Pfam overlap, rank ≤ 10):")
        print(f"  These are candidate novel geometric relationships.")
        print(f"  {'Query':>10s} {'Hit':>10s} {'Score':>7s} {'Cov':>5s} {'Match':>6s} {'Query name':>30s} {'Hit name':>30s}")
        print(f"  {'-'*10} {'-'*10} {'-'*7} {'-'*5} {'-'*6} {'-'*30} {'-'*30}")
        for gh in sorted(geometry_only_hits, key=lambda x: -x['geom_score']):
            print(f"  {gh['query']:>10s} {gh['hit']:>10s} {gh['geom_score']:7.3f} "
                  f"{gh['coverage']:5.3f} {gh['n_matched']:6d} "
                  f"{gh['query_name']:>30s} {gh['hit_name']:>30s}")
    else:
        print(f"\n  No geometry-only hits found (all top hits share Pfam domains).")
    
    # Top hit detail for each query
    print(f"\n  TOP 5 HITS PER QUERY:")
    for quid in search_results:
        pq = per_query[quid]
        q_cls = classifications.get(quid, {})
        print(f"\n  {quid} — {q_cls.get('protein_name', '?')[:50]}")
        pfams_q = [f['id'] for f in q_cls.get('pfam_families', [])]
        if pfams_q:
            print(f"    Query Pfam: {', '.join(pfams_q[:5])}")
        
        for h in pq['hits'][:5]:
            tp_mark = "✓" if h['is_tp'] else " "
            pfam_str = f"Pfam:{h['pfam_shared']}" if h['pfam_shared'] > 0 else "no Pfam"
            print(f"    {h['rank']:2d}. {h['uid']:>10s} {h['geom_score']:.3f} "
                  f"cov={h['coverage']:.2f} match={h['n_matched']:3d} "
                  f"{tp_mark} {pfam_str:>10s} {h['gene_name']:>8s} "
                  f"{h['protein_name'][:35]}")
    
    # Pfam overlap distribution
    print(f"\n  PFAM OVERLAP DISTRIBUTION (all query-hit pairs):")
    n_total = len(all_pfam_overlaps)
    n_zero = sum(1 for j in all_pfam_overlaps if j == 0)
    n_low = sum(1 for j in all_pfam_overlaps if 0 < j <= 0.3)
    n_med = sum(1 for j in all_pfam_overlaps if 0.3 < j <= 0.7)
    n_high = sum(1 for j in all_pfam_overlaps if j > 0.7)
    print(f"    Zero overlap:  {n_zero:4d} ({n_zero/n_total:.0%})")
    print(f"    Low (0-0.3):   {n_low:4d} ({n_low/n_total:.0%})")
    print(f"    Med (0.3-0.7): {n_med:4d} ({n_med/n_total:.0%})")
    print(f"    High (>0.7):   {n_high:4d} ({n_high/n_total:.0%})")
    
    # Mean reciprocal rank for first TP
    mrr_values = []
    for quid in search_results:
        pq = per_query[quid]
        for h in pq['hits']:
            if h['is_tp']:
                mrr_values.append(1.0 / h['rank'])
                break
        else:
            mrr_values.append(0.0)
    
    if mrr_values:
        mrr = sum(mrr_values) / len(mrr_values)
        print(f"\n  Mean Reciprocal Rank (first TP): {mrr:.3f}")
        n_found = sum(1 for v in mrr_values if v > 0)
        print(f"  Queries with at least 1 TP in top-{top_k}: {n_found}/{len(mrr_values)}")
    
    # ── Save ──
    output = {
        'n_queries': len(search_results),
        'n_database': len(barcodes),
        'top_k': top_k,
        'precision_at_k': {str(k): round(sum(v)/len(v), 4) if v else 0
                          for k, v in all_precisions.items()},
        'mean_reciprocal_rank': round(mrr, 4) if mrr_values else 0,
        'n_tp_found': n_found if mrr_values else 0,
        'geometry_only_hits': geometry_only_hits,
        'per_query': per_query,
        'pfam_overlap_distribution': {
            'zero': n_zero, 'low': n_low, 'med': n_med, 'high': n_high,
        },
    }
    
    try:
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        with open(BENCHMARK_OUT, 'w') as f:
            json.dump(output, f, indent=2)
        print(f"\n  Results saved: {BENCHMARK_OUT}")
    except Exception as e:
        print(f"  (Could not save: {e})")
    
    print()
    return output


def cmd_classify(barcodes, no_fetch=False):
    """Show classifications for all barcoded proteins."""
    uids = [bc.get('uniprot_id', '') for bc in barcodes if bc.get('uniprot_id')]
    
    # Start with query panel
    panel_uids = [uid for uid in DEFAULT_QUERIES if find_barcode(barcodes, uid)]
    print(f"\n  Available query proteins: {len(panel_uids)}/{len(DEFAULT_QUERIES)}")
    
    if not panel_uids:
        print("  None of the default query proteins found in database.")
        print("  Available UIDs:", [bc.get('uniprot_id','?') for bc in barcodes[:10]])
        return
    
    classifications = load_or_fetch_classifications(panel_uids, no_fetch=no_fetch)
    
    print(f"\n  {'UniProt':>10s} {'Gene':>8s} {'Pfam domains':>4s} {'SCOP SF':>12s} {'Name'}")
    print(f"  {'-'*10} {'-'*8} {'-'*4} {'-'*12} {'-'*40}")
    
    for uid in panel_uids:
        cls = classifications.get(uid, {})
        n_pfam = len(cls.get('pfam_families', []))
        scop = cls.get('scop2_superfamily', '-')[:12]
        gene = cls.get('gene_name', '-')[:8]
        name = cls.get('protein_name', '?')[:40]
        print(f"  {uid:>10s} {gene:>8s} {n_pfam:4d} {scop:>12s} {name}")
    print()


def cmd_report():
    """Load and display cached benchmark results."""
    if not BENCHMARK_OUT.exists():
        print(f"  No benchmark results found. Run: python benchmark_homology.py run")
        return
    
    with open(BENCHMARK_OUT) as f:
        results = json.load(f)
    
    print(f"\n  CACHED BENCHMARK RESULTS")
    print(f"  {'='*60}")
    print(f"  Queries: {results['n_queries']}")
    print(f"  Database: {results['n_database']}")
    print(f"  MRR: {results['mean_reciprocal_rank']}")
    
    print(f"\n  Precision@k:")
    for k, p in sorted(results['precision_at_k'].items(), key=lambda x: int(x[0])):
        print(f"    P@{k}: {p:.3f}")
    
    n_geom_only = len(results.get('geometry_only_hits', []))
    print(f"\n  Geometry-only discoveries: {n_geom_only}")
    for gh in results.get('geometry_only_hits', [])[:10]:
        print(f"    {gh['query']} → {gh['hit']}: score={gh['geom_score']:.3f}, "
              f"cov={gh['coverage']:.3f}")
    print()


# ─── Main ───
if __name__ == "__main__":
    args = sys.argv[1:]
    
    if not args:
        print(__doc__)
        sys.exit(0)
    
    cmd = args[0].lower()
    
    barcodes = load_all_barcodes()
    if not barcodes:
        print(f"  No barcodes found in {BARCODE_DIR}")
        sys.exit(1)
    print(f"  Loaded {len(barcodes)} barcodes")
    
    no_fetch = '--no-fetch' in args
    
    if cmd == 'run':
        # Custom queries?
        custom = None
        for a in args:
            if a.startswith('--queries='):
                custom = a.split('=')[1].split(',')
            elif a.startswith('--queries') and args.index(a) + 1 < len(args):
                next_a = args[args.index(a) + 1]
                if not next_a.startswith('--'):
                    custom = next_a.split(',')
        
        top_k = 20
        for i, a in enumerate(args):
            if a == '--top' and i + 1 < len(args):
                top_k = int(args[i + 1])
        
        run_benchmark(barcodes, query_ids=custom, top_k=top_k, no_fetch=no_fetch)
    
    elif cmd == 'classify':
        cmd_classify(barcodes, no_fetch=no_fetch)
    
    elif cmd == 'report':
        cmd_report()
    
    else:
        print(f"  Unknown command: {cmd}")
        print("  Commands: run, classify, report")

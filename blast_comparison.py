#!/usr/bin/env python3
"""
BLAST/sequence identity comparison for geometric homology hits.

For each benchmark query, fetches sequences for query + top N hits from UniProt,
computes pairwise sequence identity (global Needleman-Wunsch or local Smith-Waterman),
and plots geometric score vs sequence identity.

Key question: Do geometry-only hits (high geom score, no Pfam/SCOP) have LOW
sequence identity? If yes → method finds what BLAST can't.

Usage:
    python blast_comparison.py                   # Full analysis
    python blast_comparison.py --query P00533    # Single query
    python blast_comparison.py --no-fetch        # Use cached sequences
    python blast_comparison.py --report          # Show cached results

Requires: numpy. Optional: matplotlib (for plots), Biopython (for proper alignment).
Falls back to simple k-mer identity if Biopython unavailable.
"""

import json
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path
from collections import defaultdict

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "data"
SEQ_CACHE = DATA_DIR / "sequences"
FIG_DIR = SCRIPT_DIR / "figures"

sys.path.insert(0, str(SCRIPT_DIR))

import numpy as np

try:
    from Bio import pairwise2
    from Bio.Align import substitution_matrices
    HAS_BIO = True
except ImportError:
    HAS_BIO = False

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def fetch_fasta(uid, timeout=10):
    """Fetch FASTA sequence from UniProt."""
    cache_file = SEQ_CACHE / f"{uid}.fasta"
    if cache_file.exists():
        with open(cache_file) as f:
            return f.read()

    url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    try:
        req = urllib.request.Request(url, headers={
            'User-Agent': 'Mozilla/5.0 PhiPsiAtlas/1.0',
            'Accept': 'text/plain',
        })
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = resp.read().decode('utf-8', errors='replace')
            if data.startswith('>'):
                SEQ_CACHE.mkdir(parents=True, exist_ok=True)
                with open(cache_file, 'w') as f:
                    f.write(data)
                return data
    except Exception as e:
        pass
    return None


def parse_fasta(fasta_str):
    """Extract sequence from FASTA string."""
    if not fasta_str:
        return ''
    lines = fasta_str.strip().split('\n')
    seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return seq


def kmer_identity(seq1, seq2, k=3):
    """Quick approximate sequence identity using k-mer overlap.
    Not a proper alignment — just a fast proxy.
    Returns value between 0 and 1.
    """
    if not seq1 or not seq2:
        return 0.0

    kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))

    if not kmers1 or not kmers2:
        return 0.0

    overlap = len(kmers1 & kmers2)
    total = len(kmers1 | kmers2)
    return overlap / total if total > 0 else 0.0


def global_identity(seq1, seq2):
    """Compute global pairwise identity.
    Uses Biopython if available, otherwise falls back to k-mer.
    Returns fraction (0-1).
    """
    if not seq1 or not seq2:
        return 0.0

    # Truncate very long sequences to keep computation manageable
    max_len = 2000
    s1 = seq1[:max_len]
    s2 = seq2[:max_len]

    if HAS_BIO:
        try:
            # Use simple match scoring: +1 match, -1 mismatch, -2 gap open, -0.5 extend
            alignments = pairwise2.align.globalms(s1, s2, 1, -1, -2, -0.5,
                                                   one_alignment_only=True,
                                                   score_only=False)
            if alignments:
                a1, a2 = alignments[0].seqA, alignments[0].seqB
                matches = sum(1 for c1, c2 in zip(a1, a2) if c1 == c2 and c1 != '-')
                aligned_len = max(len(s1), len(s2))
                return matches / aligned_len if aligned_len > 0 else 0.0
        except Exception:
            pass

    # Fallback: k-mer identity
    return kmer_identity(s1, s2, k=3)


def analyze_query(query_uid, hits, no_fetch=False):
    """Analyze sequence identity for one query's hits.

    Returns list of dicts: {hit_uid, geom_score, seq_identity, is_tp, category}
    """
    results = []

    # Fetch query sequence
    if not no_fetch:
        q_fasta = fetch_fasta(query_uid)
        q_seq = parse_fasta(q_fasta)
        time.sleep(0.2)
    else:
        cache_file = SEQ_CACHE / f"{query_uid}.fasta"
        q_seq = parse_fasta(cache_file.read_text() if cache_file.exists() else '')

    if not q_seq:
        print(f"    Could not fetch sequence for {query_uid}")
        return results

    for hit in hits:
        hit_uid = hit.get('uid', '')
        geom_score = hit.get('geom_score', 0)
        is_tp = hit.get('is_tp', False) or hit.get('scop_match', False)
        pfam_shared = hit.get('pfam_shared', 0)

        # Categorize
        if is_tp:
            category = 'true_positive'
        elif pfam_shared > 0:
            category = 'pfam_match'
        else:
            category = 'geometry_only'

        # Fetch hit sequence
        if not no_fetch:
            h_fasta = fetch_fasta(hit_uid)
            h_seq = parse_fasta(h_fasta)
            time.sleep(0.2)
        else:
            cache_file = SEQ_CACHE / f"{hit_uid}.fasta"
            h_seq = parse_fasta(cache_file.read_text() if cache_file.exists() else '')

        if not h_seq:
            continue

        # Compute identity
        identity = global_identity(q_seq, h_seq)

        results.append({
            'query': query_uid,
            'hit': hit_uid,
            'gene_name': hit.get('gene_name', ''),
            'geom_score': geom_score,
            'seq_identity': identity,
            'is_tp': is_tp,
            'pfam_shared': pfam_shared,
            'category': category,
            'rank': hit.get('rank', 0),
        })

    return results


def plot_results(all_results, output_name='fig6_blast_comparison'):
    """Plot geometric score vs sequence identity."""
    if not HAS_MPL:
        print("  matplotlib not available for plotting.")
        return

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))

    colors = {
        'true_positive': '#27AE60',
        'pfam_match': '#3498DB',
        'geometry_only': '#E74C3C',
    }
    labels = {
        'true_positive': 'SCOP/Pfam TP',
        'pfam_match': 'Pfam match only',
        'geometry_only': 'Geometry-only hit',
    }

    plotted = set()
    for r in all_results:
        cat = r['category']
        label = labels[cat] if cat not in plotted else None
        plotted.add(cat)

        ax.scatter(r['seq_identity'] * 100, r['geom_score'],
                   c=colors[cat], s=40, alpha=0.7, edgecolors='black', lw=0.3,
                   label=label, zorder=5 if cat == 'geometry_only' else 3)

    ax.set_xlabel('Sequence identity (%)', fontsize=11)
    ax.set_ylabel('Geometric similarity score', fontsize=11)
    ax.set_title('Geometric homology vs sequence identity', fontsize=13)

    # Twilight zone
    ax.axvspan(0, 20, alpha=0.08, color='red', label='Twilight zone (<20%)')
    ax.axvline(20, color='red', ls='--', alpha=0.3)

    ax.legend(fontsize=9, loc='lower right')
    ax.set_xlim(-2, 100)
    ax.set_ylim(0.5, 1.0)

    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / f'{output_name}.svg', format='svg', bbox_inches='tight')
    fig.savefig(FIG_DIR / f'{output_name}.png', dpi=300, bbox_inches='tight')
    print(f"  Saved: {FIG_DIR / output_name}.svg/png")
    plt.close(fig)


def main():
    args = sys.argv[1:]
    no_fetch = '--no-fetch' in args
    report_only = '--report' in args
    single_query = None
    for i, a in enumerate(args):
        if a == '--query' and i + 1 < len(args):
            single_query = args[i + 1]

    results_file = DATA_DIR / "blast_comparison_results.json"

    if report_only:
        if results_file.exists():
            with open(results_file) as f:
                all_results = json.load(f)
            print(f"  Loaded {len(all_results)} cached results")
        else:
            print("  No cached results. Run: python blast_comparison.py")
            return
    else:
        # Load benchmark results
        bench_files = [
            DATA_DIR / "benchmark_results.json",
            SCRIPT_DIR / "benchmark_results.json",
            SCRIPT_DIR / "data" / "benchmark_results.json",
            Path("benchmark_results.json"),  # current working directory
            Path("data") / "benchmark_results.json",
        ]
        for p in Path("/mnt/user-data/uploads").glob("benchmark_results*.json"):
            bench_files.append(p)

        bench = None
        for bf in bench_files:
            if bf.exists():
                try:
                    with open(bf) as f:
                        bench = json.load(f)
                    break
                except (json.JSONDecodeError, IOError):
                    pass

        if bench is None:
            print("  No benchmark_results.json found. Run benchmark_homology.py first.")
            return

        all_results = []
        per_query = bench.get('per_query', {})

        for quid, pq in per_query.items():
            if single_query and quid != single_query:
                continue

            hits = pq.get('hits', [])
            print(f"  Analyzing {quid} ({len(hits)} hits)...")
            results = analyze_query(quid, hits, no_fetch=no_fetch)
            all_results.extend(results)
            print(f"    Got {len(results)} identity measurements")

        # Save
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        with open(results_file, 'w') as f:
            json.dump(all_results, f, indent=2)
        print(f"\n  Saved {len(all_results)} results to {results_file}")

    # Summary statistics
    if all_results:
        tp_ids = [r['seq_identity'] for r in all_results if r['category'] == 'true_positive']
        go_ids = [r['seq_identity'] for r in all_results if r['category'] == 'geometry_only']

        print(f"\n  SEQUENCE IDENTITY SUMMARY:")
        if tp_ids:
            print(f"    True positives (n={len(tp_ids)}): "
                  f"mean={np.mean(tp_ids)*100:.1f}%, "
                  f"median={np.median(tp_ids)*100:.1f}%, "
                  f"range=[{np.min(tp_ids)*100:.1f}%, {np.max(tp_ids)*100:.1f}%]")
        if go_ids:
            print(f"    Geometry-only (n={len(go_ids)}): "
                  f"mean={np.mean(go_ids)*100:.1f}%, "
                  f"median={np.median(go_ids)*100:.1f}%, "
                  f"range=[{np.min(go_ids)*100:.1f}%, {np.max(go_ids)*100:.1f}%]")
            twilight = sum(1 for x in go_ids if x < 0.20)
            print(f"    Geometry-only in twilight zone (<20%): "
                  f"{twilight}/{len(go_ids)} ({twilight/len(go_ids)*100:.0f}%)")

        # Plot
        plot_results(all_results)


if __name__ == "__main__":
    main()

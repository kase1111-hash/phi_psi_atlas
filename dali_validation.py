#!/usr/bin/env python3
"""
DALI structural superposition validation for geometry-only hits.

Selects the top geometry-only hits from the benchmark (high geometric score,
no SCOP/Pfam overlap) and prepares them for DALI server submission.

DALI Z-score interpretation:
  Z > 20  — very close structural match (same fold family)
  Z > 8   — significant structural similarity (same superfamily)
  Z > 2   — significant (worth investigating)
  Z < 2   — not significant

Usage:
    python dali_validation.py select                # Pick pairs for validation
    python dali_validation.py select --n 15         # Top 15 pairs
    python dali_validation.py parse results.txt     # Parse DALI output
    python dali_validation.py report                # Show cached results
    python dali_validation.py submit                # Show DALI submission URLs

Workflow:
  1. Run `python dali_validation.py select` to pick the best geometry-only pairs
  2. Submit pairs to DALI server (ekhidna2.biocenter.helsinki.fi/dali/)
  3. Download results
  4. Run `python dali_validation.py parse <results_file>` to analyze
"""

import json
import sys
from pathlib import Path
from collections import defaultdict

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "data"
FIG_DIR = SCRIPT_DIR / "figures"

sys.path.insert(0, str(SCRIPT_DIR))

try:
    import numpy as np
    HAS_NP = True
except ImportError:
    HAS_NP = False

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def load_benchmark():
    """Load benchmark results."""
    candidates = [
        DATA_DIR / "benchmark_results.json",
        SCRIPT_DIR / "benchmark_results.json",
    ]
    for p in Path("/mnt/user-data/uploads").glob("benchmark_results*.json"):
        candidates.append(p)

    for bf in candidates:
        if bf.exists():
            try:
                with open(bf) as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError):
                pass
    return None


def select_pairs(bench, n=10):
    """Select the best geometry-only pairs for DALI validation.

    Criteria:
    - No SCOP superfamily match
    - No Pfam overlap
    - Geometric score > 0.65
    - Rank <= 10

    Sorted by geometric score (highest first).
    Returns list of dicts: {query, hit, geom_score, rank, query_gene, hit_gene, ...}
    """
    pairs = []

    for quid, pq in bench.get('per_query', {}).items():
        for hit in pq.get('hits', []):
            is_tp = hit.get('is_tp', False) or hit.get('scop_match', False)
            pfam = hit.get('pfam_shared', 0)
            score = hit.get('geom_score', 0)
            rank = hit.get('rank', 99)

            # Geometry-only: no known homology
            if not is_tp and pfam == 0 and score > 0.65 and rank <= 10:
                pairs.append({
                    'query': quid,
                    'hit': hit.get('uid', ''),
                    'geom_score': score,
                    'rank': rank,
                    'n_matched': hit.get('n_matched', 0),
                    'coverage': hit.get('coverage', 0),
                    'query_gene': pq.get('query_gene', quid),
                    'hit_gene': hit.get('gene_name', ''),
                    'hit_name': hit.get('protein_name', '')[:50],
                })

    # Sort by geometric score
    pairs.sort(key=lambda x: -x['geom_score'])
    return pairs[:n]


def cmd_select(n=10):
    """Select and display pairs for DALI validation."""
    bench = load_benchmark()
    if bench is None:
        print("  No benchmark results found.")
        return

    pairs = select_pairs(bench, n=n)

    print(f"\n  TOP {len(pairs)} GEOMETRY-ONLY PAIRS FOR DALI VALIDATION")
    print(f"  (high geometric score, no SCOP/Pfam overlap)")
    print(f"  {'='*80}")
    print(f"  {'#':>3s} {'Query':>12s} {'Hit':>12s} {'Score':>6s} {'Cov':>5s} "
          f"{'Match':>5s} {'Q-gene':>8s} {'H-gene':>10s} {'Hit protein'}")
    print(f"  {'-'*3} {'-'*12} {'-'*12} {'-'*6} {'-'*5} {'-'*5} {'-'*8} {'-'*10} {'-'*30}")

    for i, p in enumerate(pairs):
        print(f"  {i+1:3d} {p['query']:>12s} {p['hit']:>12s} {p['geom_score']:6.3f} "
              f"{p['coverage']:5.2f} {p['n_matched']:5d} {p['query_gene']:>8s} "
              f"{p['hit_gene']:>10s} {p['hit_name'][:30]}")

    # Save
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    outfile = DATA_DIR / "dali_pairs.json"
    with open(outfile, 'w') as f:
        json.dump(pairs, f, indent=2)
    print(f"\n  Saved to {outfile}")

    # Generate AlphaFold structure IDs for DALI submission
    print(f"\n  DALI SUBMISSION INFO:")
    print(f"  Server: https://ekhidna2.biocenter.helsinki.fi/dali/")
    print(f"  Use 'pairwise' mode with AlphaFold structure IDs:")
    print()
    for i, p in enumerate(pairs):
        af_q = f"AF-{p['query']}-F1"
        af_h = f"AF-{p['hit']}-F1"
        print(f"  Pair {i+1}: {af_q}  vs  {af_h}  "
              f"(geom={p['geom_score']:.3f}, {p['query_gene']}→{p['hit_gene']})")

    return pairs


def cmd_submit():
    """Show URLs for manual DALI submission."""
    pairs_file = DATA_DIR / "dali_pairs.json"
    if not pairs_file.exists():
        print("  Run `python dali_validation.py select` first.")
        return

    with open(pairs_file) as f:
        pairs = json.load(f)

    print(f"\n  DALI SUBMISSION GUIDE ({len(pairs)} pairs)")
    print(f"  Server: https://ekhidna2.biocenter.helsinki.fi/dali/")
    print(f"  Mode: Pairwise structure comparison")
    print()
    print(f"  For each pair, upload both CIF files or use AlphaFold IDs.")
    print(f"  Record the Z-score from the DALI output.")
    print()

    for i, p in enumerate(pairs):
        print(f"  Pair {i+1}: {p['query']} ({p['query_gene']}) vs "
              f"{p['hit']} ({p['hit_gene']})")
        print(f"    Geometric score: {p['geom_score']:.3f}")
        print(f"    AlphaFold: AF-{p['query']}-F1 vs AF-{p['hit']}-F1")
        print()


def cmd_parse(results_file):
    """Parse DALI results file.

    Expects tab-separated or space-separated output with Z-scores.
    Can handle various DALI output formats.
    """
    results_path = Path(results_file)
    if not results_path.exists():
        print(f"  File not found: {results_file}")
        return

    # Load our pairs for reference
    pairs_file = DATA_DIR / "dali_pairs.json"
    pairs = []
    if pairs_file.exists():
        with open(pairs_file) as f:
            pairs = json.load(f)

    # Parse DALI output
    print(f"\n  Parsing DALI results from {results_file}...")

    with open(results_path) as f:
        content = f.read()

    # Try to extract Z-scores
    # DALI output varies, but Z-scores are typically labeled
    z_scores = []

    for line in content.split('\n'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        # Look for lines with Z-score values
        parts = line.split()

        # Common DALI format: "No:  Chain   Z    rmsd lali nres %id Description"
        for i, part in enumerate(parts):
            if part.upper() == 'Z' or part.upper() == 'Z-SCORE':
                # Next numeric value is likely the Z-score
                continue
            try:
                z = float(part)
                if 0 < z < 100 and i > 0:  # Plausible Z-score range
                    # Try to identify protein IDs in the same line
                    z_scores.append({'line': line, 'z_score': z, 'raw_parts': parts})
                    break
            except ValueError:
                continue

    if z_scores:
        print(f"  Found {len(z_scores)} potential Z-scores")
        for zs in z_scores:
            print(f"    Z = {zs['z_score']:.1f}  |  {zs['line'][:80]}")
    else:
        print(f"  Could not auto-parse Z-scores from file.")
        print(f"  Please enter Z-scores manually:")
        z_scores = manual_entry(pairs)

    return z_scores


def manual_entry(pairs):
    """Allow manual entry of DALI Z-scores."""
    results = []
    if not pairs:
        print("  No pairs loaded. Run `select` first.")
        return results

    print(f"\n  Enter Z-scores for each pair (or 'skip' / 'done'):")
    for i, p in enumerate(pairs):
        prompt = (f"  Pair {i+1}: {p['query_gene']} vs {p['hit_gene']} "
                  f"(geom={p['geom_score']:.3f}) Z-score = ")
        try:
            val = input(prompt).strip()
        except (EOFError, KeyboardInterrupt):
            break

        if val.lower() == 'done':
            break
        if val.lower() == 'skip' or not val:
            continue

        try:
            z = float(val)
            p['dali_z'] = z
            results.append(p)
        except ValueError:
            print(f"    Invalid value: {val}")

    if results:
        outfile = DATA_DIR / "dali_results.json"
        with open(outfile, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\n  Saved {len(results)} results to {outfile}")

    return results


def cmd_report():
    """Show cached DALI validation results."""
    results_file = DATA_DIR / "dali_results.json"
    if not results_file.exists():
        print("  No DALI results cached. Run validation workflow first.")
        return

    with open(results_file) as f:
        results = json.load(f)

    print(f"\n  DALI VALIDATION RESULTS ({len(results)} pairs)")
    print(f"  {'='*80}")

    significant = 0
    strong = 0

    for r in results:
        z = r.get('dali_z', 0)
        sig = '***' if z > 8 else '**' if z > 4 else '*' if z > 2 else 'ns'

        if z > 2:
            significant += 1
        if z > 8:
            strong += 1

        print(f"  {r.get('query_gene','?'):>8s} vs {r.get('hit_gene','?'):>10s}  "
              f"geom={r['geom_score']:.3f}  DALI Z={z:5.1f}  {sig}")

    print(f"\n  Summary:")
    print(f"    Significant (Z > 2): {significant}/{len(results)} "
          f"({significant/len(results)*100:.0f}%)" if results else "")
    print(f"    Strong (Z > 8): {strong}/{len(results)} "
          f"({strong/len(results)*100:.0f}%)" if results else "")

    # Correlation
    if HAS_NP and len(results) > 2:
        geom_scores = [r['geom_score'] for r in results]
        z_scores = [r.get('dali_z', 0) for r in results]
        corr = np.corrcoef(geom_scores, z_scores)[0, 1]
        print(f"    Correlation (geom vs Z): r = {corr:.3f}")

    # Plot if possible
    if HAS_MPL and results:
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))

        geom = [r['geom_score'] for r in results]
        z = [r.get('dali_z', 0) for r in results]
        labels_text = [f"{r.get('query_gene','')}-{r.get('hit_gene','')}" for r in results]

        colors = ['#27AE60' if zi > 2 else '#E74C3C' for zi in z]
        ax.scatter(geom, z, c=colors, s=60, edgecolors='black', lw=0.5, zorder=5)

        for i, label in enumerate(labels_text):
            ax.annotate(label, (geom[i], z[i]), textcoords='offset points',
                       xytext=(5, 5), fontsize=7)

        ax.axhline(2, color='red', ls='--', alpha=0.4, label='Z = 2 (significance)')
        ax.axhline(8, color='green', ls='--', alpha=0.4, label='Z = 8 (superfamily)')
        ax.set_xlabel('Geometric similarity score', fontsize=11)
        ax.set_ylabel('DALI Z-score', fontsize=11)
        ax.set_title('Geometric homology vs structural superposition', fontsize=12)
        ax.legend(fontsize=9)

        fig.tight_layout()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / 'fig7_dali_validation.svg', bbox_inches='tight')
        fig.savefig(FIG_DIR / 'fig7_dali_validation.png', dpi=300, bbox_inches='tight')
        print(f"  Saved: {FIG_DIR}/fig7_dali_validation.svg/png")
        plt.close(fig)


def main():
    args = sys.argv[1:] or ['select']
    cmd = args[0]

    if cmd == 'select':
        n = 10
        for i, a in enumerate(args):
            if a == '--n' and i + 1 < len(args):
                n = int(args[i + 1])
        cmd_select(n=n)

    elif cmd == 'submit':
        cmd_submit()

    elif cmd == 'parse':
        if len(args) < 2:
            print("  Usage: python dali_validation.py parse <results_file>")
            return
        cmd_parse(args[1])

    elif cmd == 'report':
        cmd_report()

    elif cmd == 'manual':
        # Manual Z-score entry
        pairs_file = DATA_DIR / "dali_pairs.json"
        if pairs_file.exists():
            with open(pairs_file) as f:
                pairs = json.load(f)
            manual_entry(pairs)
        else:
            print("  Run `select` first.")

    else:
        print(f"  Unknown command: {cmd}")
        print("  Commands: select, submit, parse, manual, report")


if __name__ == "__main__":
    main()

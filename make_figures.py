#!/usr/bin/env python3
"""
Generate publication-quality figures for the curvature modulation paper.

Figures:
  1. winding    — Winding vector scatter plot (p vs q), colored by SS class
  2. barcode    — Example curvature barcode (EGFR or any protein)
  3. null       — Basin-control null comparison (real vs null bar chart)
  4. precision  — Precision vs protein size from benchmark
  5. alignment  — Segment alignment heatmap between two proteins

Usage:
    python make_figures.py all                    # Generate all figures
    python make_figures.py winding                # Just the winding plot
    python make_figures.py barcode P00533         # Barcode for EGFR
    python make_figures.py barcode O75110         # Barcode for ATPase
    python make_figures.py null                   # Basin null comparison
    python make_figures.py precision              # Precision vs size
    python make_figures.py alignment O75110 O43861  # ATPase vs ATP9B alignment

Outputs: figures/ directory with SVG + PNG at 300 DPI
"""

import json
import math
import sys
from pathlib import Path
from collections import defaultdict

SCRIPT_DIR = Path(__file__).parent
BARCODE_DIR = SCRIPT_DIR / "results" / "barcodes"
DATA_DIR = SCRIPT_DIR / "data"
FIG_DIR = SCRIPT_DIR / "figures"

sys.path.insert(0, str(SCRIPT_DIR))

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("  matplotlib not found. Install: pip install matplotlib")
    print("  Trying: pip install matplotlib --break-system-packages")

# ── Style ──
COLORS = {
    'helix': '#E74C3C',      # red
    'strand': '#3498DB',     # blue
    'coil': '#95A5A6',       # gray
    'mixed': '#9B59B6',      # purple
    'helix_rich': '#E74C3C',
    'strand_rich': '#3498DB',
    'coil_rich': '#95A5A6',
    'alpha_beta': '#2ECC71',  # green
    'real': '#2C3E50',       # dark navy
    'null': '#BDC3C7',       # light gray
    'tp': '#27AE60',         # green
    'fp': '#E74C3C',         # red
    'geom_only': '#F39C12',  # orange
}

FONT = {'family': 'serif', 'size': 10}


def load_barcodes():
    """Load all barcodes."""
    from geom_homology import load_all_barcodes
    return load_all_barcodes()


def classify_ss(bc):
    """Classify by SS composition."""
    ss = bc.get('ss_counts', {})
    total = sum(ss.values()) or 1
    h = ss.get('H', 0) / total
    e = ss.get('E', 0) / total
    if h > 0.50 and e < 0.10:
        return 'helix_rich'
    elif e > 0.30 and h < 0.15:
        return 'strand_rich'
    elif h > 0.20 and e > 0.15:
        return 'alpha_beta'
    else:
        return 'mixed'


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 1: WINDING VECTOR SCATTER
# ═══════════════════════════════════════════════════════════════════════

def fig_winding(barcodes):
    """Winding vector scatter plot (p vs q), colored by structural class."""
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    groups = defaultdict(lambda: ([], []))
    for bc in barcodes:
        p = bc.get('p_winding', 0)
        q = bc.get('q_winding', 0)
        cls = classify_ss(bc)
        groups[cls][0].append(p)
        groups[cls][1].append(q)

    order = ['helix_rich', 'strand_rich', 'alpha_beta', 'mixed']
    labels = {'helix_rich': 'Helix-rich', 'strand_rich': 'Strand-rich',
              'alpha_beta': r'$\alpha/\beta$', 'mixed': 'Mixed'}

    for cls in order:
        ps, qs = groups[cls]
        ax.scatter(ps, qs, c=COLORS[cls], s=8, alpha=0.5, label=labels[cls],
                   edgecolors='none', rasterized=True)

    ax.set_xlabel(r'$p$-winding ($\phi$ net turns)', fontsize=11)
    ax.set_ylabel(r'$q$-winding ($\psi$ net turns)', fontsize=11)
    ax.set_title(r'Winding vectors on $T^2$', fontsize=13)
    ax.axhline(0, color='k', lw=0.5, ls='--', alpha=0.3)
    ax.axvline(0, color='k', lw=0.5, ls='--', alpha=0.3)
    ax.legend(fontsize=9, loc='upper left', framealpha=0.8)
    ax.set_xlim(-35, 35)
    ax.set_ylim(-35, 15)

    # Annotate quadrant percentages
    n = len(barcodes)
    quads = {'p+q-': 0, 'p-q-': 0, 'p+q+': 0, 'p-q+': 0}
    for bc in barcodes:
        p = bc.get('p_winding', 0)
        q = bc.get('q_winding', 0)
        key = f"p{'+' if p >= 0 else '-'}q{'+' if q >= 0 else '-'}"
        quads[key] += 1
    ax.text(20, -30, f"{quads['p+q-']/n:.0%}", fontsize=8, color='gray', ha='center')
    ax.text(-20, -30, f"{quads['p-q-']/n:.0%}", fontsize=8, color='gray', ha='center')
    ax.text(20, 10, f"{quads['p+q+']/n:.0%}", fontsize=8, color='gray', ha='center')
    ax.text(-20, 10, f"{quads['p-q+']/n:.0%}", fontsize=8, color='gray', ha='center')

    fig.tight_layout()
    save(fig, 'fig1_winding')
    return fig


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 2: EXAMPLE CURVATURE BARCODE
# ═══════════════════════════════════════════════════════════════════════

def fig_barcode(barcodes, uid='P00533'):
    """Render a single protein's curvature barcode as a figure."""
    from geom_homology import find_barcode, MACRO
    bc = find_barcode(barcodes, uid)
    if not bc:
        print(f"  Protein {uid} not found.")
        return None

    segments = bc.get('segments', [])
    if not segments:
        print(f"  No segments for {uid}.")
        return None

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 3.5), height_ratios=[3, 1],
                                     sharex=True, gridspec_kw={'hspace': 0.08})

    ss_colors = {'H': COLORS['helix'], 'E': COLORS['strand'], 'C': COLORS['coil']}

    # Macro-category alpha for visual distinction
    macro_alpha = {
        'constant_k': 0.4,
        'monotone': 0.6,
        'polynomial': 0.7,
        'localized': 0.9,
        'oscillatory': 0.8,
        'transition': 0.75,
    }

    # Compute cumulative residue positions from segment lengths
    # (start_residue is not stored in barcode JSON)
    positions = []
    pos = 0
    for seg in segments:
        seg_len = seg.get('length', 10)
        positions.append((pos, seg_len))
        pos += seg_len

    max_residue = pos
    kappa_values = []

    for i, seg in enumerate(segments):
        start, seg_len = positions[i]
        ss = seg.get('ss_type', 'C')
        model = seg.get('model_class', 'circular_arc')
        mean_k = abs(seg.get('mean_kappa', 0))
        kappa_values.append(mean_k)

        macro = MACRO.get(model, 'constant_k')
        color = ss_colors.get(ss, COLORS['coil'])
        alpha = macro_alpha.get(macro, 0.5)

        # Top panel: curvature bars at residue positions
        ax1.bar(start, mean_k, width=seg_len, align='edge',
                color=color, alpha=alpha, edgecolor='white', linewidth=0.3)

        # Bottom panel: SS track
        ax2.bar(start, 1, width=seg_len, align='edge',
                color=color, alpha=0.85, edgecolor='white', linewidth=0.3)

    # Mark Gaussian peaks with triangles
    for i, seg in enumerate(segments):
        if seg.get('model_class') == 'gauss_peak':
            start, seg_len = positions[i]
            mid = start + seg_len / 2
            mean_k = abs(seg.get('mean_kappa', 0))
            ax1.plot(mid, mean_k * 1.05, 'v', color='black', markersize=4, zorder=10)

    # Curvature panel
    ax1.set_ylabel(r'$|\kappa|$ (rad/residue)', fontsize=10)
    if kappa_values:
        ax1.set_ylim(0, max(kappa_values) * 1.15)

    name = bc.get('protein_name', uid)[:45]
    gene = bc.get('gene_name', '')
    length_aa = bc.get('length', max_residue)
    ax1.set_title(f'{uid} ({gene}) — {name}\n'
                  f'{length_aa} aa, {len(segments)} segments, |Q|={bc.get("Q_magnitude", 0):.1f}',
                  fontsize=10)

    # SS track
    ax2.set_ylabel('SS', fontsize=9)
    ax2.set_xlabel('Residue position', fontsize=10)
    ax2.set_ylim(0, 1)
    ax2.set_yticks([])

    ax1.set_xlim(0, max_residue)

    # Legend
    patches = [mpatches.Patch(color=COLORS['helix'], label='Helix', alpha=0.7),
               mpatches.Patch(color=COLORS['strand'], label='Strand', alpha=0.7),
               mpatches.Patch(color=COLORS['coil'], label='Coil', alpha=0.7)]
    ax1.legend(handles=patches, fontsize=8, loc='upper right', framealpha=0.8,
               ncol=3)

    fig.tight_layout()
    save(fig, f'fig2_barcode_{uid}')
    return fig


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 3: BASIN-CONTROL NULL COMPARISON
# ═══════════════════════════════════════════════════════════════════════

def fig_null():
    """Basin-control null comparison bar chart."""
    # Try loading from file first, then use hardcoded n=557 data
    null_file = DATA_DIR / "basin_null_results.json"
    if null_file.exists():
        with open(null_file) as f:
            data = json.load(f)
        results = data['results']
    else:
        # Hardcoded from n=557, Nrep=5 run
        results = {
            'H': {'real_nonconst_frac': 0.6509, 'null_nonconst_frac': 0.4486,
                   'p_binary': 1e-230, 'direction': 'enriched'},
            'E': {'real_nonconst_frac': 0.3171, 'null_nonconst_frac': 0.1141,
                   'p_binary': 1e-242, 'direction': 'enriched'},
            'C': {'real_nonconst_frac': 0.2839, 'null_nonconst_frac': 0.3307,
                   'p_binary': 1e-20, 'direction': 'suppressed'},
        }

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    ss_types = ['H', 'E', 'C']
    ss_labels = ['Helix', 'Strand', 'Coil']
    x = np.arange(len(ss_types))
    width = 0.35

    real_vals = [results[s]['real_nonconst_frac'] * 100 for s in ss_types]
    null_vals = [results[s]['null_nonconst_frac'] * 100 for s in ss_types]

    bars1 = ax.bar(x - width/2, real_vals, width, label='Real', color=COLORS['real'], alpha=0.85)
    bars2 = ax.bar(x + width/2, null_vals, width, label='Null (basin)', color=COLORS['null'], alpha=0.85)

    # Significance stars
    for i, s in enumerate(ss_types):
        p = results[s].get('p_binary', 0)
        direction = results[s].get('direction', '')
        if p < 1e-100:
            stars = '***'
        elif p < 1e-10:
            stars = '**'
        elif p < 0.001:
            stars = '*'
        else:
            stars = 'ns'

        y_max = max(real_vals[i], null_vals[i])
        arrow = '↑' if direction == 'enriched' else '↓'
        ax.text(i, y_max + 3, f'{stars} {arrow}', ha='center', fontsize=9, fontweight='bold')

    ax.set_ylabel('Non-constant segments (%)', fontsize=11)
    ax.set_xlabel('Secondary structure type', fontsize=11)
    ax.set_title('Curvature modulation vs basin-control null\n'
                 f'(n = 557 proteins, $N_{{rep}}$ = 5)', fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(ss_labels, fontsize=11)
    ax.set_ylim(0, 85)
    ax.legend(fontsize=10, loc='upper right')

    # Bootstrap CI for helix
    if 'bootstrap' in results.get('H', {}):
        ci_lo = results['H']['bootstrap']['ci_95_low'] * 100
        ci_hi = results['H']['bootstrap']['ci_95_high'] * 100
        ax.errorbar(0 - width/2, real_vals[0], 
                    yerr=[[real_vals[0] - ci_lo], [ci_hi - real_vals[0]]],
                    fmt='none', ecolor='black', capsize=3, lw=1.5)

    fig.tight_layout()
    save(fig, 'fig3_null_comparison')
    return fig


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 4: PRECISION VS PROTEIN SIZE
# ═══════════════════════════════════════════════════════════════════════

def fig_precision():
    """Precision vs protein size from benchmark results."""
    # Try loading benchmark results
    bench_files = [
        SCRIPT_DIR / "data" / "benchmark_results.json",
        DATA_DIR / "benchmark_results.json",
    ]
    # Also check uploads
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
        # Hardcoded from 20-query benchmark (5,563 proteins)
        queries = [
            {'uid': 'P69905', 'label': 'Hb-α', 'segs': 11, 'tp5': 2, 'tp20': 4, 'mrr': 1.0},
            {'uid': 'P68871', 'label': 'Hb-β', 'segs': 13, 'tp5': 4, 'tp20': 4, 'mrr': 1.0},
            {'uid': 'P02144', 'label': 'Myoglobin', 'segs': 12, 'tp5': 1, 'tp20': 2, 'mrr': 1.0},
            {'uid': 'P00918', 'label': 'CA2', 'segs': 30, 'tp5': 4, 'tp20': 6, 'mrr': 1.0},
            {'uid': 'P00533', 'label': 'EGFR', 'segs': 100, 'tp5': 4, 'tp20': 13, 'mrr': 1.0},
            {'uid': 'P04626', 'label': 'HER2', 'segs': 102, 'tp5': 5, 'tp20': 12, 'mrr': 1.0},
            {'uid': 'P21860', 'label': 'ErbB3', 'segs': 104, 'tp5': 5, 'tp20': 10, 'mrr': 1.0},
            {'uid': 'P01857', 'label': 'IgG1', 'segs': 38, 'tp5': 5, 'tp20': 12, 'mrr': 1.0},
            {'uid': 'P01834', 'label': 'Ig-κ', 'segs': 12, 'tp5': 5, 'tp20': 16, 'mrr': 1.0},
            {'uid': 'P02751', 'label': 'FN1', 'segs': 100, 'tp5': 3, 'tp20': 9, 'mrr': 1.0},
            {'uid': 'P00338', 'label': 'LDH-A', 'segs': 30, 'tp5': 5, 'tp20': 7, 'mrr': 1.0},
            {'uid': 'P04406', 'label': 'GAPDH', 'segs': 30, 'tp5': 1, 'tp20': 1, 'mrr': 1.0},
            {'uid': 'P62258', 'label': '14-3-3ε', 'segs': 22, 'tp5': 4, 'tp20': 5, 'mrr': 1.0},
            {'uid': 'P11166', 'label': 'GLUT1', 'segs': 45, 'tp5': 4, 'tp20': 6, 'mrr': 1.0},
            {'uid': 'P06748', 'label': 'NPM1', 'segs': 15, 'tp5': 1, 'tp20': 1, 'mrr': 1.0},
            {'uid': 'P07900', 'label': 'HSP90α', 'segs': 65, 'tp5': 2, 'tp20': 2, 'mrr': 1.0},
            {'uid': 'P38646', 'label': 'GRP75', 'segs': 55, 'tp5': 5, 'tp20': 9, 'mrr': 1.0},
            {'uid': 'P10636', 'label': 'Tau', 'segs': 30, 'tp5': 0, 'tp20': 1, 'mrr': 0.07},
            {'uid': 'P60174', 'label': 'TPI', 'segs': 20, 'tp5': 0, 'tp20': 0, 'mrr': 0.0},
            {'uid': 'P01308', 'label': 'Insulin', 'segs': 6, 'tp5': 0, 'tp20': 0, 'mrr': 0.0},
        ]
    else:
        queries = []
        barcodes = load_barcodes()
        bc_index = {bc.get('uniprot_id', ''): bc for bc in barcodes}
        for quid, pq in bench.get('per_query', {}).items():
            hits = pq.get('hits', [])
            tp5 = sum(1 for h in hits[:5] if h.get('is_tp') or h.get('scop_match'))
            tp20 = sum(1 for h in hits[:20] if h.get('is_tp') or h.get('scop_match'))
            # Get segment count
            bc = bc_index.get(quid, {})
            segs = bc.get('n_segments', len(bc.get('segments', [])))
            if segs == 0:
                segs = hits[0].get('n_matched', 10) if hits else 10
            # MRR
            mrr = 0
            for i, h in enumerate(hits):
                if h.get('is_tp') or h.get('scop_match'):
                    mrr = 1.0 / (i + 1)
                    break
            gene = bc.get('gene_name', quid[:6])
            queries.append({'uid': quid, 'label': gene, 'segs': segs,
                           'tp5': tp5, 'tp20': tp20, 'mrr': mrr})

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    segs = [q['segs'] for q in queries]
    tp5s = [q['tp5'] / 5.0 for q in queries]  # precision@5
    mrrs = [q['mrr'] for q in queries]
    labels = [q['label'] for q in queries]

    # Left: P@5 vs segments
    for i, q in enumerate(queries):
        color = COLORS['tp'] if q['mrr'] > 0 else COLORS['fp']
        ax1.scatter(q['segs'], q['tp5'] / 5.0, c=color, s=80, zorder=5, edgecolors='black', lw=0.5)
        ax1.annotate(q['label'], (q['segs'], q['tp5'] / 5.0),
                    textcoords='offset points', xytext=(5, 5), fontsize=8)

    ax1.set_xlabel('Query segments', fontsize=11)
    ax1.set_ylabel('Precision@5', fontsize=11)
    ax1.set_title('Retrieval precision vs query complexity', fontsize=12)
    ax1.set_ylim(-0.05, 1.1)
    ax1.axvline(20, color='red', ls='--', alpha=0.4, label='~20 seg threshold')
    ax1.legend(fontsize=9)

    # Right: MRR vs segments
    for i, q in enumerate(queries):
        color = COLORS['tp'] if q['mrr'] > 0 else COLORS['fp']
        ax2.scatter(q['segs'], q['mrr'], c=color, s=80, zorder=5, edgecolors='black', lw=0.5)
        ax2.annotate(q['label'], (q['segs'], q['mrr']),
                    textcoords='offset points', xytext=(5, 5), fontsize=8)

    ax2.set_xlabel('Query segments', fontsize=11)
    ax2.set_ylabel('Mean Reciprocal Rank', fontsize=11)
    ax2.set_title('First TP rank vs query complexity', fontsize=12)
    ax2.set_ylim(-0.05, 1.1)
    ax2.axvline(20, color='red', ls='--', alpha=0.4)

    fig.tight_layout()
    save(fig, 'fig4_precision_vs_size')
    return fig


# ═══════════════════════════════════════════════════════════════════════
#  FIGURE 5: ALIGNMENT HEATMAP
# ═══════════════════════════════════════════════════════════════════════

def fig_alignment(barcodes, uid1='O75110', uid2='O43861'):
    """Segment alignment between two proteins as a dot-plot style figure."""
    from geom_homology import find_barcode, align_barcodes, segment_score

    bc1 = find_barcode(barcodes, uid1)
    bc2 = find_barcode(barcodes, uid2)
    if not bc1 or not bc2:
        print(f"  Cannot find {uid1} or {uid2}")
        return None

    segs1 = bc1.get('segments', [])
    segs2 = bc2.get('segments', [])

    n1, n2 = len(segs1), len(segs2)
    if n1 == 0 or n2 == 0:
        print("  No segments to align.")
        return None

    # Compute full score matrix
    score_matrix = np.zeros((n1, n2))
    for i in range(n1):
        for j in range(n2):
            score_matrix[i, j] = segment_score(segs1[i], segs2[j])

    # Get the actual alignment path
    _, alignment, n_matched = align_barcodes(segs1, segs2)

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    # Heatmap of pairwise scores
    im = ax.imshow(score_matrix, cmap='YlOrRd', aspect='auto',
                   origin='lower', interpolation='nearest',
                   vmin=-0.5, vmax=1.5)

    # Overlay alignment path
    if alignment:
        ai = [a[0] for a in alignment]
        aj = [a[1] for a in alignment]
        ax.plot(aj, ai, 'k-', lw=1.5, alpha=0.6)
        ax.scatter(aj, ai, c='black', s=15, zorder=5)

    # Color bar for SS type along axes
    # Top axis: protein 2 SS
    for j, seg in enumerate(segs2):
        ss = seg.get('ss_type', 'C')
        color = {'H': COLORS['helix'], 'E': COLORS['strand'], 'C': COLORS['coil']}.get(ss, 'gray')
        ax.plot(j, n1 + 1, 's', color=color, markersize=3)

    # Left axis: protein 1 SS
    for i, seg in enumerate(segs1):
        ss = seg.get('ss_type', 'C')
        color = {'H': COLORS['helix'], 'E': COLORS['strand'], 'C': COLORS['coil']}.get(ss, 'gray')
        ax.plot(-2, i, 's', color=color, markersize=3)

    gene1 = bc1.get('gene_name', '') or uid1
    gene2 = bc2.get('gene_name', '') or uid2
    ax.set_xlabel(f'{uid2} ({gene2}) segments', fontsize=11)
    ax.set_ylabel(f'{uid1} ({gene1}) segments', fontsize=11)
    ax.set_title(f'Curvature segment alignment: {gene1} vs {gene2}\n'
                 f'{n_matched} matched segments', fontsize=12)

    plt.colorbar(im, ax=ax, label='Segment score', shrink=0.8)

    fig.tight_layout()
    save(fig, f'fig5_alignment_{uid1}_{uid2}')
    return fig


# ═══════════════════════════════════════════════════════════════════════
#  UTILITIES
# ═══════════════════════════════════════════════════════════════════════

def save(fig, name):
    """Save figure as SVG and PNG."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    svg_path = FIG_DIR / f'{name}.svg'
    png_path = FIG_DIR / f'{name}.png'
    fig.savefig(svg_path, format='svg', bbox_inches='tight')
    fig.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
    print(f"  Saved: {svg_path}")
    print(f"  Saved: {png_path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    if not HAS_MPL:
        sys.exit(1)

    matplotlib.rc('font', **FONT)

    args = sys.argv[1:] or ['all']
    cmd = args[0]

    # Load barcodes for most figures
    barcodes = None
    if cmd in ('all', 'winding', 'barcode', 'alignment'):
        barcodes = load_barcodes()
        print(f"  Loaded {len(barcodes)} barcodes")

    if cmd == 'all':
        fig_winding(barcodes)
        fig_barcode(barcodes, 'P00533')
        fig_barcode(barcodes, 'O75110')
        fig_null()
        fig_precision()
        fig_alignment(barcodes, 'O75110', 'O43861')
        print(f"\n  All figures saved to {FIG_DIR}/")

    elif cmd == 'winding':
        fig_winding(barcodes)

    elif cmd == 'barcode':
        uid = args[1] if len(args) > 1 else 'P00533'
        fig_barcode(barcodes, uid)

    elif cmd == 'null':
        fig_null()

    elif cmd == 'precision':
        fig_precision()

    elif cmd == 'alignment':
        uid1 = args[1] if len(args) > 1 else 'O75110'
        uid2 = args[2] if len(args) > 2 else 'O43861'
        if barcodes is None:
            barcodes = load_barcodes()
        fig_alignment(barcodes, uid1, uid2)

    else:
        print(f"  Unknown command: {cmd}")
        print("  Commands: all, winding, barcode, null, precision, alignment")


if __name__ == "__main__":
    main()

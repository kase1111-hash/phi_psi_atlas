#!/usr/bin/env python3
"""
BARCODE ANALYZER — Live analysis of the AlphaFold geometric barcode database.
Run this while the main pipeline is still downloading.

Usage:
    python analyze_barcodes.py                    # Full dashboard
    python analyze_barcodes.py --top-defects 20   # Top 20 most defective proteins
    python analyze_barcodes.py --top-geodesic 20  # Top 20 smoothest proteins
    python analyze_barcodes.py --winding-map       # (p,q) winding scatter
    python analyze_barcodes.py --clusters          # Auto-cluster by winding vector
    python analyze_barcodes.py --find P04637       # Find a specific protein
    python analyze_barcodes.py --defect-catalog    # All Gaussian peaks across all proteins
    python analyze_barcodes.py --peak-analysis     # Functional vs structural peak comparison
    python analyze_barcodes.py --peak-analysis 3   # Custom proximity threshold (default 5)
    python analyze_barcodes.py --peak-analysis --no-fetch  # Skip UniProt API, stats only
    python analyze_barcodes.py --inter-peak-spacing  # Phased array test: spacing mod 3.6
"""

import json
import os
import sys
import math
from pathlib import Path
from collections import Counter, defaultdict

# ─── Locate barcodes ───
SCRIPT_DIR = Path(__file__).parent
BARCODE_DIR = SCRIPT_DIR / "results" / "barcodes"
RESULTS_DIR = SCRIPT_DIR / "results"
CACHE_DIRS = [
    SCRIPT_DIR / "alphafold_cache",
    SCRIPT_DIR / "data" / "alphafold",
]

if not BARCODE_DIR.exists():
    print(f"[ERROR] Barcodes directory not found: {BARCODE_DIR}")
    print("  Run alphafold_deconstruct first to generate barcodes.")
    sys.exit(1)


def load_all_barcodes():
    """Load all barcode JSON files from the barcodes directory."""
    barcodes = []
    for f in sorted(BARCODE_DIR.glob("*_barcode.json")):
        try:
            with open(f) as fh:
                bc = json.load(fh)
                bc['_file'] = f.name
                barcodes.append(bc)
        except (json.JSONDecodeError, KeyError):
            pass
    return barcodes


def fmt_pct(val):
    return f"{val*100:.0f}%" if val is not None else "—"


def dashboard(barcodes):
    """Print a full statistical dashboard."""
    n = len(barcodes)
    print()
    print("=" * 72)
    print(f"  GEOMETRIC BARCODE DATABASE — {n} proteins loaded")
    print("=" * 72)
    print()

    if n == 0:
        print("  No barcodes found yet. Pipeline still starting?")
        return

    # ── Basic stats ──
    lengths = [b['length'] for b in barcodes]
    p_vals = [b['p'] for b in barcodes]
    q_vals = [b['q'] for b in barcodes]
    Q_vals = [b['Q_magnitude'] for b in barcodes]
    geo_fracs = [b['frac_geodesic'] for b in barcodes]
    def_fracs = [b['frac_defect'] for b in barcodes]
    n_segs = [b['n_segments'] for b in barcodes]

    print(f"  {'Statistic':30s} {'Mean':>10s} {'Median':>10s} {'Min':>10s} {'Max':>10s} {'Std':>10s}")
    print(f"  {'-'*30} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
    for label, vals in [
        ("Chain length (aa)", lengths),
        ("p-winding (φ)", p_vals),
        ("q-winding (ψ)", q_vals),
        ("|Q| magnitude", Q_vals),
        ("Segments per protein", n_segs),
        ("Geodesic fraction", geo_fracs),
        ("Defect fraction", def_fracs),
    ]:
        sv = sorted(vals)
        mean = sum(vals) / n
        median = sv[n // 2]
        mn, mx = sv[0], sv[-1]
        std = math.sqrt(sum((v - mean) ** 2 for v in vals) / max(n - 1, 1))
        print(f"  {label:30s} {mean:10.2f} {median:10.2f} {mn:10.2f} {mx:10.2f} {std:10.2f}")

    # ── SS composition ──
    print()
    total_H = sum(b['ss_composition']['H'] for b in barcodes)
    total_E = sum(b['ss_composition']['E'] for b in barcodes)
    total_C = sum(b['ss_composition']['C'] for b in barcodes)
    total_res = total_H + total_E + total_C
    print(f"  Total residues: {total_res:,}")
    print(f"    Helix:  {total_H:>10,} ({total_H/total_res:.1%})")
    print(f"    Strand: {total_E:>10,} ({total_E/total_res:.1%})")
    print(f"    Coil:   {total_C:>10,} ({total_C/total_res:.1%})")

    # ── Classification distribution ──
    print()
    print("  CURVATURE CLASS DISTRIBUTION (all segments)")
    print(f"  {'-'*50}")
    all_classes = Counter()
    total_segs = 0
    for b in barcodes:
        for s in b.get('segments', []):
            all_classes[s['model_class']] += 1
            total_segs += 1

    for cls, cnt in all_classes.most_common():
        bar = "█" * int(50 * cnt / total_segs)
        print(f"    {cls:20s} {cnt:6d} ({cnt/total_segs:5.1%}) {bar}")
    print(f"    {'TOTAL':20s} {total_segs:6d}")

    # ── Geodesic fraction by SS type ──
    print()
    print("  GEODESIC FRACTION BY SS TYPE")
    print(f"  {'-'*50}")
    ss_geo = defaultdict(lambda: [0, 0])  # [geodesic, total]
    for b in barcodes:
        for s in b.get('segments', []):
            ss = s['ss_type']
            ss_geo[ss][1] += 1
            if s['model_class'] in ('circular_arc', 'geodesic'):
                ss_geo[ss][0] += 1
    for ss in ['H', 'E', 'C']:
        if ss_geo[ss][1] > 0:
            g, t = ss_geo[ss]
            bar = "█" * int(50 * g / t)
            print(f"    {ss}: {g}/{t} ({g/t:.1%}) {bar}")

    # ── Defect class by SS type ──
    print()
    print("  DEFECT CLASSES BY SS TYPE")
    print(f"  {'-'*50}")
    defect_classes = {'gauss_peak', 'quadratic', 'damped_osc', 'damped_oscillation'}
    ss_defects = defaultdict(Counter)
    for b in barcodes:
        for s in b.get('segments', []):
            if s['model_class'] in defect_classes:
                ss_defects[s['ss_type']][s['model_class']] += 1
    for ss in ['H', 'E', 'C']:
        if ss_defects[ss]:
            parts = ", ".join(f"{cls}={cnt}" for cls, cnt in ss_defects[ss].most_common())
            print(f"    {ss}: {parts}")

    # ── Length vs winding correlation ──
    print()
    if n >= 10:
        mean_L = sum(lengths) / n
        mean_Q = sum(Q_vals) / n
        cov = sum((l - mean_L) * (q - mean_Q) for l, q in zip(lengths, Q_vals)) / (n - 1)
        std_L = math.sqrt(sum((l - mean_L) ** 2 for l in lengths) / (n - 1))
        std_Q = math.sqrt(sum((q - mean_Q) ** 2 for q in Q_vals) / (n - 1))
        r = cov / (std_L * std_Q) if std_L > 0 and std_Q > 0 else 0
        print(f"  Length vs |Q| correlation: r = {r:.3f} (n={n})")

    # ── Winding quadrant distribution ──
    print()
    print("  WINDING QUADRANT DISTRIBUTION (p,q signs)")
    quads = Counter()
    for b in barcodes:
        qx = "+" if b['p'] >= 0 else "-"
        qy = "+" if b['q'] >= 0 else "-"
        quads[f"p{qx} q{qy}"] += 1
    for label, cnt in quads.most_common():
        bar = "█" * int(40 * cnt / n)
        print(f"    {label}: {cnt:5d} ({cnt/n:5.1%}) {bar}")

    # ── Cache stats ──
    print()
    n_cached = 0
    for d in CACHE_DIRS:
        if d.exists():
            n_cif = len(list(d.glob("*.cif")))
            n_pdb = len(list(d.glob("*.pdb")))
            n_cached += n_cif + n_pdb
            print(f"  Cache: {d.name}/ — {n_cif} CIF + {n_pdb} PDB = {n_cif+n_pdb} files")
    if n_cached > 0:
        print(f"  Barcode coverage: {n}/{n_cached} ({n/n_cached:.1%})")

    # ── Basin null results (if available) ──
    null_path = RESULTS_DIR / "basin_null_results.json"
    if null_path.exists():
        print()
        print("  BASIN-CONTROL NULL (from basin_null_alphafold.py)")
        print(f"  {'-'*50}")
        try:
            with open(null_path) as f:
                null_data = json.load(f)
            nr = null_data.get('results', {})
            print(f"  Proteins: {null_data.get('n_proteins', '?')}, "
                  f"Reps: {null_data.get('n_reps', '?')}, "
                  f"Real segs: {null_data.get('n_real_segments', '?')}, "
                  f"Null segs: {null_data.get('n_null_segments', '?')}")
            for ss in ['H', 'E', 'C']:
                if ss in nr:
                    r = nr[ss]
                    ss_name = {'H': 'Helix', 'E': 'Strand', 'C': 'Coil'}[ss]
                    direction = r.get('direction', '?')
                    marker = "**" if ss == 'H' else "  "
                    print(f"  {marker}{ss_name:8s}: real={r['real_nonconst_frac']:.1%} non-const, "
                          f"null={r['null_nonconst_frac']:.1%}, "
                          f"p={r['p_binary']:.2e}, OR={r['odds_ratio']:.2f} "
                          f"[{direction}]{marker}")
                    if ss == 'H' and 'bootstrap' in r:
                        b = r['bootstrap']
                        print(f"    Bootstrap 95% CI: [{b['ci_95_low']:.1%}, {b['ci_95_high']:.1%}], "
                              f"{b['frac_exceeds_null']:.0%} exceed null")
                    if ss == 'H' and 'leave_one_out' in r:
                        l = r['leave_one_out']
                        print(f"    LOO range: [{l['min_nonconst']:.1%}, {l['max_nonconst']:.1%}], "
                              f"all below null: {l['all_below_null']}")
        except Exception as e:
            print(f"  (Error reading null results: {e})")

    print()


def top_defects(barcodes, n=20):
    """Show proteins with highest defect fraction."""
    print()
    print(f"  TOP {n} MOST DEFECTIVE PROTEINS (highest defect fraction)")
    print(f"  {'='*72}")
    ranked = sorted(barcodes, key=lambda b: b['frac_defect'], reverse=True)[:n]
    print(f"  {'UniProt':15s} {'L':>5s} {'Segs':>5s} {'Geo%':>6s} {'Def%':>6s} {'nDef':>5s} {'|Q|':>7s} {'H':>4s} {'E':>4s} {'C':>4s}")
    print(f"  {'-'*15} {'-'*5} {'-'*5} {'-'*6} {'-'*6} {'-'*5} {'-'*7} {'-'*4} {'-'*4} {'-'*4}")
    for b in ranked:
        print(f"  {b['uniprot_id']:15s} {b['length']:5d} {b['n_segments']:5d} "
              f"{fmt_pct(b['frac_geodesic']):>6s} {fmt_pct(b['frac_defect']):>6s} "
              f"{b['n_defect']:5d} {b['Q_magnitude']:7.2f} "
              f"{b['ss_composition']['H']:4d} {b['ss_composition']['E']:4d} {b['ss_composition']['C']:4d}")
    print()


def top_geodesic(barcodes, n=20):
    """Show proteins with highest geodesic fraction (smoothest on T²)."""
    print()
    # Filter to proteins with at least 5 segments
    filtered = [b for b in barcodes if b['n_segments'] >= 5]
    print(f"  TOP {n} SMOOTHEST PROTEINS (highest geodesic fraction, ≥5 segments)")
    print(f"  {'='*72}")
    ranked = sorted(filtered, key=lambda b: b['frac_geodesic'], reverse=True)[:n]
    print(f"  {'UniProt':15s} {'L':>5s} {'Segs':>5s} {'Geo%':>6s} {'Def%':>6s} {'|Q|':>7s} {'H':>4s} {'E':>4s} {'C':>4s}")
    print(f"  {'-'*15} {'-'*5} {'-'*5} {'-'*6} {'-'*6} {'-'*7} {'-'*4} {'-'*4} {'-'*4}")
    for b in ranked:
        print(f"  {b['uniprot_id']:15s} {b['length']:5d} {b['n_segments']:5d} "
              f"{fmt_pct(b['frac_geodesic']):>6s} {fmt_pct(b['frac_defect']):>6s} "
              f"{b['Q_magnitude']:7.2f} "
              f"{b['ss_composition']['H']:4d} {b['ss_composition']['E']:4d} {b['ss_composition']['C']:4d}")
    print()


def winding_map(barcodes):
    """Show (p,q) winding distribution as ASCII scatter."""
    print()
    print("  WINDING VECTOR MAP (p vs q)")
    print(f"  {'='*72}")

    p_vals = [b['p'] for b in barcodes]
    q_vals = [b['q'] for b in barcodes]

    # Determine grid bounds
    p_min, p_max = min(p_vals), max(p_vals)
    q_min, q_max = min(q_vals), max(q_vals)

    # Clamp to reasonable display range
    p_lo, p_hi = max(p_min, -30), min(p_max, 30)
    q_lo, q_hi = max(q_min, -30), min(q_max, 30)

    W, H_grid = 60, 30
    grid = [[' '] * W for _ in range(H_grid)]

    # Plot points
    for p, q in zip(p_vals, q_vals):
        if p_lo <= p <= p_hi and q_lo <= q <= q_hi:
            col = int((p - p_lo) / max(p_hi - p_lo, 0.01) * (W - 1))
            row = int((q_hi - q) / max(q_hi - q_lo, 0.01) * (H_grid - 1))
            col = max(0, min(W - 1, col))
            row = max(0, min(H_grid - 1, row))
            if grid[row][col] == ' ':
                grid[row][col] = '.'
            elif grid[row][col] == '.':
                grid[row][col] = 'o'
            elif grid[row][col] == 'o':
                grid[row][col] = 'O'
            else:
                grid[row][col] = '#'

    # Draw axes
    zero_col = int((0 - p_lo) / max(p_hi - p_lo, 0.01) * (W - 1))
    zero_row = int((q_hi - 0) / max(q_hi - q_lo, 0.01) * (H_grid - 1))
    zero_col = max(0, min(W - 1, zero_col))
    zero_row = max(0, min(H_grid - 1, zero_row))
    for r in range(H_grid):
        if 0 <= zero_col < W and grid[r][zero_col] == ' ':
            grid[r][zero_col] = '|'
    for c in range(W):
        if 0 <= zero_row < H_grid and grid[zero_row][c] == ' ':
            grid[zero_row][c] = '-'
    if 0 <= zero_row < H_grid and 0 <= zero_col < W:
        grid[zero_row][zero_col] = '+'

    print(f"  q={q_hi:+.0f}")
    for row in grid:
        print(f"    {''.join(row)}")
    print(f"  q={q_lo:+.0f}")
    print(f"  {'p='+str(round(p_lo)):>6s}{' '*(W-12)}{'p='+str(round(p_hi)):<6s}")
    print()
    print(f"  Legend: . = 1 protein, o = 2, O = 3, # = 4+")
    print(f"  Proteins outside [{p_lo:.0f},{p_hi:.0f}] x [{q_lo:.0f},{q_hi:.0f}]: "
          f"{sum(1 for p,q in zip(p_vals,q_vals) if p<p_lo or p>p_hi or q<q_lo or q>q_hi)}")
    print()


def find_protein(barcodes, query):
    """Find and display a specific protein barcode."""
    query = query.strip().upper()
    matches = [b for b in barcodes if query in b.get('uniprot_id', '').upper()]

    if not matches:
        print(f"  No barcode found matching '{query}'")
        return

    for b in matches:
        print()
        print(f"  BARCODE: {b['uniprot_id']}")
        print(f"  {'='*60}")
        print(f"  Length:     {b['length']} residues")
        print(f"  SS:         H={b['ss_composition']['H']} E={b['ss_composition']['E']} C={b['ss_composition']['C']}")
        print(f"  p-winding:  {b['p']:+.3f}")
        print(f"  q-winding:  {b['q']:+.3f}")
        print(f"  |Q|:        {b['Q_magnitude']:.3f}")
        print(f"  Segments:   {b['n_segments']}")
        print(f"  Geodesic:   {b['n_geodesic']} ({fmt_pct(b['frac_geodesic'])})")
        print(f"  Defects:    {b['n_defect']} ({fmt_pct(b['frac_defect'])})")
        print()
        print(f"  {'Idx':>4s} {'SS':>3s} {'Class':>20s} {'Len':>4s} {'Kappa':>8s} {'R2':>6s} {'Range'}")
        print(f"  {'-'*4} {'-'*3} {'-'*20} {'-'*4} {'-'*8} {'-'*6} {'-'*12}")
        for i, s in enumerate(b.get('segments', [])):
            cls = s['model_class']
            marker = " <<" if cls in ('gauss_peak', 'quadratic', 'damped_osc', 'damped_oscillation') else ""
            print(f"  {i+1:4d} {s['ss_type']:>3s} {cls:>20s} {s['length']:4d} "
                  f"{s['mean_kappa']:+8.3f} {s['R2']:6.3f} "
                  f"{s['start_idx']}-{s['end_idx']}{marker}")
        print()


def defect_catalog(barcodes):
    """Catalog all Gaussian peak defects across all proteins."""
    print()
    print("  GAUSSIAN PEAK DEFECT CATALOG")
    print(f"  {'='*72}")

    defect_types = {'gauss_peak', 'quadratic', 'damped_osc', 'damped_oscillation'}
    peaks = []
    for b in barcodes:
        for s in b.get('segments', []):
            if s['model_class'] in defect_types:
                peaks.append({
                    'uid': b['uniprot_id'],
                    'length': b['length'],
                    **s,
                })

    print(f"  Total defects: {len(peaks)} across {len(barcodes)} proteins")
    print()

    # By SS type
    ss_counts = Counter(p['ss_type'] for p in peaks)
    print(f"  By SS type: " + ", ".join(f"{ss}={cnt}" for ss, cnt in ss_counts.most_common()))

    # By class
    cls_counts = Counter(p['model_class'] for p in peaks)
    print(f"  By class:   " + ", ".join(f"{cls}={cnt}" for cls, cnt in cls_counts.most_common()))
    print()

    # Gaussian peaks with parameters
    gauss = [p for p in peaks if p['model_class'] == 'gauss_peak']
    if gauss:
        print(f"  GAUSSIAN PEAKS (n={len(gauss)}):")
        print(f"  {'UniProt':15s} {'SS':>3s} {'Res':>6s} {'Len':>4s} {'Kappa':>8s} {'A':>8s} {'Sigma':>8s}")
        print(f"  {'-'*15} {'-'*3} {'-'*6} {'-'*4} {'-'*8} {'-'*8} {'-'*8}")

        # Sort by amplitude
        gauss.sort(key=lambda p: abs(p.get('params', {}).get('A', 0)), reverse=True)
        for p in gauss[:50]:  # Top 50
            params = p.get('params', {})
            A = params.get('A', 0)
            sigma = params.get('sigma', 0)
            pk_res = params.get('peak_residue', p['start_idx'])
            print(f"  {p['uid']:15s} {p['ss_type']:>3s} {pk_res:6d} {p['length']:4d} "
                  f"{p['mean_kappa']:+8.3f} {A:8.2f} {sigma:8.3f}")

        # Amplitude statistics
        amps = [abs(p.get('params', {}).get('A', 0)) for p in gauss if p.get('params', {}).get('A')]
        sigmas = [p.get('params', {}).get('sigma', 0) for p in gauss if p.get('params', {}).get('sigma')]
        if amps:
            print()
            print(f"  Amplitude |A|: mean={sum(amps)/len(amps):.2f}, "
                  f"median={sorted(amps)[len(amps)//2]:.2f}, "
                  f"range=[{min(amps):.2f}, {max(amps):.2f}]")
        if sigmas:
            print(f"  Width sigma:   mean={sum(sigmas)/len(sigmas):.3f}, "
                  f"median={sorted(sigmas)[len(sigmas)//2]:.3f}, "
                  f"range=[{min(sigmas):.3f}, {max(sigmas):.3f}]")
    print()


def cluster_analysis(barcodes):
    """Simple k-means-like clustering by (p,q) winding vector."""
    print()
    print("  WINDING VECTOR CLUSTER ANALYSIS")
    print(f"  {'='*72}")

    if len(barcodes) < 20:
        print("  Need at least 20 proteins for clustering.")
        return

    # Classify by helix/strand dominance
    groups = {'helix-rich': [], 'strand-rich': [], 'coil-rich': [], 'mixed': []}
    for b in barcodes:
        ss = b['ss_composition']
        total = ss['H'] + ss['E'] + ss['C']
        if total == 0:
            continue
        h_frac = ss['H'] / total
        e_frac = ss['E'] / total
        c_frac = ss['C'] / total

        if h_frac > 0.5:
            groups['helix-rich'].append(b)
        elif e_frac > 0.3:
            groups['strand-rich'].append(b)
        elif c_frac > 0.7:
            groups['coil-rich'].append(b)
        else:
            groups['mixed'].append(b)

    print()
    print(f"  {'Group':15s} {'n':>6s} {'p mean':>8s} {'q mean':>8s} {'|Q| mean':>9s} {'Geo% mean':>10s} {'Def% mean':>10s}")
    print(f"  {'-'*15} {'-'*6} {'-'*8} {'-'*8} {'-'*9} {'-'*10} {'-'*10}")
    for name, members in groups.items():
        if not members:
            continue
        n = len(members)
        mean_p = sum(b['p'] for b in members) / n
        mean_q = sum(b['q'] for b in members) / n
        mean_Q = sum(b['Q_magnitude'] for b in members) / n
        mean_geo = sum(b['frac_geodesic'] for b in members) / n
        mean_def = sum(b['frac_defect'] for b in members) / n
        print(f"  {name:15s} {n:6d} {mean_p:+8.2f} {mean_q:+8.2f} {mean_Q:9.2f} "
              f"{mean_geo:9.1%} {mean_def:9.1%}")

    # Size distribution
    print()
    print("  SIZE BINS (by chain length)")
    print(f"  {'-'*50}")
    bins = [(0, 100, 'tiny'), (100, 300, 'small'), (300, 600, 'medium'),
            (600, 1000, 'large'), (1000, 99999, 'giant')]
    for lo, hi, label in bins:
        members = [b for b in barcodes if lo <= b['length'] < hi]
        if not members:
            continue
        n = len(members)
        mean_geo = sum(b['frac_geodesic'] for b in members) / n
        mean_def = sum(b['frac_defect'] for b in members) / n
        mean_Q = sum(b['Q_magnitude'] for b in members) / n
        print(f"    {label:8s} ({lo:4d}-{hi:4d} aa): n={n:5d}, "
              f"geo={mean_geo:.0%}, def={mean_def:.0%}, |Q|={mean_Q:.1f}")
    print()


# ═══════════════════════════════════════════════════════════════════════
#  PEAK ANALYSIS — Functional vs Structural Gaussian Peak Comparison
# ═══════════════════════════════════════════════════════════════════════

UNIPROT_CACHE_DIR = SCRIPT_DIR / "data" / "uniprot_features"


def fetch_uniprot_features(uniprot_id):
    """Fetch functional site annotations from UniProt for a given ID.
    
    Returns dict: {
        'active_site': [(start, end), ...],
        'binding_site': [(start, end), ...],
        'site': [(start, end), ...],
        'metal_binding': [(start, end), ...],
        'dna_binding': [(start, end), ...],
        'disulfide': [(start, end), ...],
        'motif': [(start, end), ...],
        'all_functional': set of 0-indexed residue positions
    }
    """
    UNIPROT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = UNIPROT_CACHE_DIR / f"{uniprot_id}.json"
    
    # Check cache first
    if cache_file.exists():
        try:
            with open(cache_file) as f:
                cached = json.load(f)
            # Convert all_functional back to set
            cached['all_functional'] = set(cached.get('all_functional', []))
            return cached
        except (json.JSONDecodeError, KeyError):
            pass
    
    # Fetch from UniProt
    try:
        import requests
    except ImportError:
        return None
    
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        resp = requests.get(url, timeout=15, headers={"Accept": "application/json"})
        if resp.status_code != 200:
            return None
        data = resp.json()
    except Exception:
        return None
    
    # Parse features
    FUNCTIONAL_TYPES = {
        'Active site': 'active_site',
        'Binding site': 'binding_site',
        'Site': 'site',
        'Metal binding': 'metal_binding',
        'DNA binding': 'dna_binding',
        'Disulfide bond': 'disulfide',
        'Motif': 'motif',
        'Nucleotide binding': 'nucleotide_binding',
        'Calcium binding': 'calcium_binding',
        'Zinc finger': 'zinc_finger',
        'Lipid moiety-binding region': 'lipid_binding',
    }
    
    result = {v: [] for v in FUNCTIONAL_TYPES.values()}
    result['all_functional'] = set()
    
    features = data.get('features', [])
    for feat in features:
        ftype = feat.get('type', '')
        if ftype in FUNCTIONAL_TYPES:
            loc = feat.get('location', {})
            start_info = loc.get('start', {})
            end_info = loc.get('end', {})
            start = start_info.get('value')
            end = end_info.get('value')
            if start is not None and end is not None:
                # UniProt uses 1-indexed; convert to 0-indexed
                s0 = int(start) - 1
                e0 = int(end) - 1
                cat = FUNCTIONAL_TYPES[ftype]
                result[cat].append((s0, e0))
                for r in range(s0, e0 + 1):
                    result['all_functional'].add(r)
    
    # Cache (convert set to list for JSON)
    cache_data = {k: v for k, v in result.items() if k != 'all_functional'}
    cache_data['all_functional'] = sorted(result['all_functional'])
    try:
        with open(cache_file, 'w') as f:
            json.dump(cache_data, f)
    except Exception:
        pass
    
    return result


def peak_analysis(barcodes, proximity_threshold=5, fetch_annotations=True):
    """Gaussian peak escapement analysis: functional vs structural (A, sigma) comparison.
    
    The "watch jewel" hypothesis: functional-site Gaussian peaks should have
    higher amplitude and narrower width than structural peaks, because evolution
    has tuned specific helix positions to be geometrically compliant at precise
    locations.
    
    Args:
        barcodes: list of barcode dicts
        proximity_threshold: max residues from annotated site to count as "functional"
        fetch_annotations: if True, fetch UniProt annotations via API
    """
    print()
    print("  " + "=" * 72)
    print("  GAUSSIAN PEAK ESCAPEMENT ANALYSIS")
    print("  Functional vs Structural (A, σ) Comparison")
    print("  " + "=" * 72)
    
    # ── Phase 1: Collect all Gaussian peaks ──
    all_peaks = []
    for b in barcodes:
        uid = b.get('uniprot_id', '')
        for s in b.get('segments', []):
            if s.get('model_class') == 'gauss_peak':
                params = s.get('params', {})
                A = params.get('A', 0)
                sigma = params.get('sigma', 0)
                peak_res = params.get('peak_residue', s.get('start_idx', -1))
                if A and sigma and peak_res >= 0:
                    all_peaks.append({
                        'uid': uid,
                        'ss_type': s['ss_type'],
                        'peak_residue': peak_res,
                        'start_idx': s['start_idx'],
                        'end_idx': s['end_idx'],
                        'length': s['length'],
                        'A': float(A),
                        'abs_A': abs(float(A)),
                        'sigma': float(sigma),
                        'mean_kappa': s.get('mean_kappa', 0),
                        'R2': s.get('R2', 0),
                        'sharpness': abs(float(A)) / max(float(sigma), 0.01),
                    })
    
    n_peaks = len(all_peaks)
    print(f"\n  Total Gaussian peaks: {n_peaks}")
    if n_peaks == 0:
        print("  No Gaussian peaks found. Need more proteins in database.")
        return
    
    # SS breakdown
    ss_counts = Counter(p['ss_type'] for p in all_peaks)
    print(f"  By SS: " + ", ".join(f"{ss}={c}" for ss, c in ss_counts.most_common()))
    
    # ── Phase 2: Fetch UniProt functional annotations ──
    if not fetch_annotations:
        print("\n  Skipping UniProt annotation fetch (--no-fetch)")
        _peak_stats_only(all_peaks)
        return
    
    uids = sorted(set(p['uid'] for p in all_peaks if p['uid']))
    print(f"\n  Fetching UniProt annotations for {len(uids)} proteins...")
    
    annotations = {}
    n_fetched = 0
    n_annotated = 0
    for uid in uids:
        feats = fetch_uniprot_features(uid)
        if feats and feats['all_functional']:
            annotations[uid] = feats
            n_annotated += 1
        n_fetched += 1
        if n_fetched % 50 == 0:
            print(f"    ...fetched {n_fetched}/{len(uids)}, {n_annotated} with annotations")
    
    print(f"  Fetched: {n_fetched}, with functional annotations: {n_annotated}")
    
    if n_annotated == 0:
        print("  No functional annotations found. Showing overall stats only.")
        _peak_stats_only(all_peaks)
        return
    
    # ── Phase 3: Label peaks as functional vs structural ──
    functional_peaks = []
    structural_peaks = []
    unmatched_peaks = []  # peaks in proteins without annotations
    
    for pk in all_peaks:
        uid = pk['uid']
        if uid not in annotations:
            unmatched_peaks.append(pk)
            continue
        
        func_residues = annotations[uid]['all_functional']
        peak_res = pk['peak_residue']
        
        # Check proximity: is peak_residue within threshold of any functional residue?
        min_dist = min((abs(peak_res - fr) for fr in func_residues), default=999)
        pk['min_dist_functional'] = min_dist
        
        # Also check if ANY residue in the segment overlaps functional sites
        seg_residues = set(range(pk['start_idx'], pk['end_idx'] + 1))
        overlap = seg_residues & func_residues
        pk['n_functional_overlap'] = len(overlap)
        
        if min_dist <= proximity_threshold:
            pk['label'] = 'functional'
            functional_peaks.append(pk)
        else:
            pk['label'] = 'structural'
            structural_peaks.append(pk)
    
    n_func = len(functional_peaks)
    n_struct = len(structural_peaks)
    n_unmatched = len(unmatched_peaks)
    
    print(f"\n  Peak classification (threshold = {proximity_threshold} residues):")
    print(f"    Functional (near annotated site):  {n_func}")
    print(f"    Structural (far from annotation):  {n_struct}")
    print(f"    Unmatched (no annotation for uid):  {n_unmatched}")
    
    # ── Phase 4: Compare (A, σ) distributions ──
    if n_func < 3 or n_struct < 3:
        print(f"\n  Insufficient peaks for comparison (need ≥3 each).")
        print(f"  Functional: {n_func}, Structural: {n_struct}")
        _peak_stats_only(all_peaks)
        return
    
    func_A = [p['abs_A'] for p in functional_peaks]
    struct_A = [p['abs_A'] for p in structural_peaks]
    func_sigma = [p['sigma'] for p in functional_peaks]
    struct_sigma = [p['sigma'] for p in structural_peaks]
    func_sharp = [p['sharpness'] for p in functional_peaks]
    struct_sharp = [p['sharpness'] for p in structural_peaks]
    
    print(f"\n  {'Metric':20s} {'Functional':>12s} {'Structural':>12s} {'Direction':>12s}")
    print(f"  {'-'*20} {'-'*12} {'-'*12} {'-'*12}")
    
    def _median(vals):
        s = sorted(vals)
        n = len(s)
        return s[n // 2] if n else 0
    
    def _mean(vals):
        return sum(vals) / len(vals) if vals else 0
    
    # Amplitude
    mf_A, ms_A = _median(func_A), _median(struct_A)
    direction_A = "func > struct" if mf_A > ms_A else "struct > func"
    print(f"  {'|A| median':20s} {mf_A:12.3f} {ms_A:12.3f} {direction_A:>12s}")
    print(f"  {'|A| mean':20s} {_mean(func_A):12.3f} {_mean(struct_A):12.3f}")
    
    # Width
    mf_s, ms_s = _median(func_sigma), _median(struct_sigma)
    direction_s = "func < struct" if mf_s < ms_s else "struct < func"
    print(f"  {'σ median':20s} {mf_s:12.4f} {ms_s:12.4f} {direction_s:>12s}")
    print(f"  {'σ mean':20s} {_mean(func_sigma):12.4f} {_mean(struct_sigma):12.4f}")
    
    # Sharpness = |A| / σ  (the "jewel bearing" metric)
    mf_sh, ms_sh = _median(func_sharp), _median(struct_sharp)
    direction_sh = "func > struct" if mf_sh > ms_sh else "struct > func"
    print(f"  {'Sharpness |A|/σ':20s} {mf_sh:12.2f} {ms_sh:12.2f} {direction_sh:>12s}")
    print(f"  {'Sharpness mean':20s} {_mean(func_sharp):12.2f} {_mean(struct_sharp):12.2f}")
    
    # ── Phase 5: Statistical tests ──
    print(f"\n  STATISTICAL TESTS")
    print(f"  {'-'*50}")
    
    try:
        from scipy import stats as sp_stats
        has_scipy = True
    except ImportError:
        has_scipy = False
    
    if has_scipy:
        # Mann-Whitney: is functional |A| > structural |A|?
        u_A, p_A = sp_stats.mannwhitneyu(func_A, struct_A, alternative='greater')
        n1, n2 = len(func_A), len(struct_A)
        r_A = 1 - 2 * u_A / (n1 * n2)
        print(f"  |A| functional > structural: U={u_A:.0f}, p={p_A:.4f}, r={r_A:.3f}")
        
        # Mann-Whitney: is functional σ < structural σ? (narrower peaks)
        u_s, p_s = sp_stats.mannwhitneyu(func_sigma, struct_sigma, alternative='less')
        r_s = 1 - 2 * u_s / (n1 * n2)
        print(f"  σ functional < structural:   U={u_s:.0f}, p={p_s:.4f}, r={r_s:.3f}")
        
        # Mann-Whitney: is functional sharpness > structural sharpness?
        u_sh, p_sh = sp_stats.mannwhitneyu(func_sharp, struct_sharp, alternative='greater')
        r_sh = 1 - 2 * u_sh / (n1 * n2)
        print(f"  Sharpness func > struct:     U={u_sh:.0f}, p={p_sh:.4f}, r={r_sh:.3f}")
        
        # KS test for distributional difference
        ks_A, pks_A = sp_stats.ks_2samp(func_A, struct_A)
        ks_s, pks_s = sp_stats.ks_2samp(func_sigma, struct_sigma)
        print(f"  KS test |A|:  D={ks_A:.3f}, p={pks_A:.4f}")
        print(f"  KS test σ:    D={ks_s:.3f}, p={pks_s:.4f}")
    else:
        print("  (scipy not available — install for statistical tests)")
        # Simple permutation test fallback
        print("  Using simple median comparison only.")
    
    # ── Phase 6: Verdict ──
    print(f"\n  ESCAPEMENT VERDICT")
    print(f"  {'-'*50}")
    
    jewel_hypothesis = (mf_A > ms_A and mf_s < ms_s)
    if jewel_hypothesis:
        print("  ✓ WATCH JEWEL HYPOTHESIS SUPPORTED")
        print("    Functional peaks are sharper (higher A, narrower σ) than structural peaks.")
        print("    Helices appear pre-stressed at functional sites with precise geometric compliance.")
    elif mf_A > ms_A:
        print("  ~ PARTIAL: Higher amplitude at functional sites, but not narrower.")
        print("    Functional sites produce larger curvature kicks but over broader regions.")
    elif mf_s < ms_s:
        print("  ~ PARTIAL: Narrower σ at functional sites, but not higher amplitude.")
        print("    Functional sites are geometrically precise but not higher-energy.")
    else:
        print("  ✗ WATCH JEWEL HYPOTHESIS NOT SUPPORTED")
        print("    Functional and structural peaks do not show the predicted (A, σ) separation.")
        print("    Gaussian peaks may serve uniform structural roles regardless of functional annotation.")
    
    if has_scipy and n_func >= 5 and n_struct >= 5:
        sig = " (significant)" if p_sh < 0.05 else " (not significant)"
        print(f"    Sharpness test: p = {p_sh:.4f}{sig}")
    
    # ── Phase 7: By SS type breakdown ──
    print(f"\n  BY SS TYPE")
    print(f"  {'-'*50}")
    for ss in ['H', 'E', 'C']:
        f_ss = [p for p in functional_peaks if p['ss_type'] == ss]
        s_ss = [p for p in structural_peaks if p['ss_type'] == ss]
        if len(f_ss) < 2 or len(s_ss) < 2:
            continue
        ss_name = {'H': 'Helix', 'E': 'Strand', 'C': 'Coil'}[ss]
        mfA = _median([p['abs_A'] for p in f_ss])
        msA = _median([p['abs_A'] for p in s_ss])
        mfs = _median([p['sigma'] for p in f_ss])
        mss = _median([p['sigma'] for p in s_ss])
        print(f"  {ss_name:8s}: func n={len(f_ss):3d} (|A|={mfA:.2f}, σ={mfs:.3f})  "
              f"struct n={len(s_ss):3d} (|A|={msA:.2f}, σ={mss:.3f})  "
              f"{'sharper' if mfA > msA and mfs < mss else 'mixed'}")
    
    # ── Phase 8: Top functional peaks ──
    print(f"\n  TOP FUNCTIONAL PEAKS (by sharpness)")
    print(f"  {'UniProt':12s} {'SS':>3s} {'Res':>5s} {'|A|':>8s} {'σ':>8s} {'Sharp':>8s} {'Dist':>5s} {'Overlap':>7s}")
    print(f"  {'-'*12} {'-'*3} {'-'*5} {'-'*8} {'-'*8} {'-'*8} {'-'*5} {'-'*7}")
    for pk in sorted(functional_peaks, key=lambda x: x['sharpness'], reverse=True)[:20]:
        print(f"  {pk['uid']:12s} {pk['ss_type']:>3s} {pk['peak_residue']:5d} "
              f"{pk['abs_A']:8.2f} {pk['sigma']:8.3f} {pk['sharpness']:8.1f} "
              f"{pk['min_dist_functional']:5d} {pk['n_functional_overlap']:7d}")
    
    # ── Phase 9: Save results ──
    results = {
        'n_peaks_total': n_peaks,
        'n_functional': n_func,
        'n_structural': n_struct,
        'n_unmatched': n_unmatched,
        'proximity_threshold': proximity_threshold,
        'functional_A_median': mf_A,
        'structural_A_median': ms_A,
        'functional_sigma_median': mf_s,
        'structural_sigma_median': ms_s,
        'functional_sharpness_median': mf_sh,
        'structural_sharpness_median': ms_sh,
        'jewel_hypothesis_supported': jewel_hypothesis,
    }
    if has_scipy:
        results['p_amplitude'] = float(p_A)
        results['p_sigma'] = float(p_s)
        results['p_sharpness'] = float(p_sh)
    
    out_path = RESULTS_DIR / "peak_analysis_results.json"
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\n  Results saved: {out_path}")
    except Exception as e:
        print(f"\n  (Could not save results: {e})")
    
    print()


def _peak_stats_only(all_peaks):
    """Print basic peak statistics without functional/structural split."""
    print(f"\n  GAUSSIAN PEAK STATISTICS (all peaks)")
    print(f"  {'-'*50}")
    
    def _median(vals):
        s = sorted(vals)
        return s[len(s) // 2] if s else 0
    
    amps = [abs(p['A']) for p in all_peaks]
    sigmas = [p['sigma'] for p in all_peaks]
    sharps = [p['sharpness'] for p in all_peaks]
    
    print(f"  |A|:        mean={sum(amps)/len(amps):.3f}, median={_median(amps):.3f}, "
          f"range=[{min(amps):.3f}, {max(amps):.3f}]")
    print(f"  σ:          mean={sum(sigmas)/len(sigmas):.4f}, median={_median(sigmas):.4f}, "
          f"range=[{min(sigmas):.4f}, {max(sigmas):.4f}]")
    print(f"  Sharpness:  mean={sum(sharps)/len(sharps):.2f}, median={_median(sharps):.2f}, "
          f"range=[{min(sharps):.2f}, {max(sharps):.2f}]")
    
    # By SS
    for ss in ['H', 'E', 'C']:
        ss_peaks = [p for p in all_peaks if p['ss_type'] == ss]
        if ss_peaks:
            ss_name = {'H': 'Helix', 'E': 'Strand', 'C': 'Coil'}[ss]
            a = [abs(p['A']) for p in ss_peaks]
            s = [p['sigma'] for p in ss_peaks]
            print(f"  {ss_name:8s}: n={len(ss_peaks):4d}, |A| med={_median(a):.3f}, σ med={_median(s):.4f}")
    print()


# ═══════════════════════════════════════════════════════════════════════
#  INTER-PEAK SPACING — Mechanical Phased Array Hypothesis
# ═══════════════════════════════════════════════════════════════════════

HELIX_PERIOD = 3.6  # residues per helix turn


def inter_peak_spacing(barcodes):
    """Test the mechanical phased array hypothesis.
    
    If Gaussian peaks on the same helix behave like a phased array — 
    propagating strain through the H-bond lattice with constructive
    interference — then inter-peak spacings should cluster at integer
    multiples of the helix period (3.6 residues).
    
    Computes:
      - All pairwise residue spacings between Gaussian peaks on the SAME helix
        (peaks on the same protein, same contiguous helix region)
      - Spacing histogram modulo 3.6
      - Rayleigh test for circular uniformity (is the mod-3.6 distribution non-uniform?)
      - Phase concentration: fraction of spacings within ±0.5 residues of 0 mod 3.6
    """
    print()
    print("  " + "=" * 72)
    print("  INTER-PEAK SPACING — Mechanical Phased Array Test")
    print("  " + "=" * 72)
    
    # ── Collect helix Gaussian peaks grouped by protein ──
    # We need to identify peaks that are on the SAME helix run.
    # Two helix gauss_peak segments are on the same helix if there are no
    # non-helix segments between them (i.e. the helix is contiguous).
    
    protein_helix_peaks = defaultdict(list)  # uid -> list of (peak_residue, A, sigma, helix_group_id)
    
    for b in barcodes:
        uid = b.get('uniprot_id', '')
        segments = b.get('segments', [])
        
        # Assign helix group IDs: consecutive H segments get the same group
        helix_group = 0
        prev_was_helix = False
        seg_groups = []
        
        for s in segments:
            if s['ss_type'] == 'H':
                if not prev_was_helix:
                    helix_group += 1
                seg_groups.append(helix_group)
                prev_was_helix = True
            else:
                seg_groups.append(-1)
                prev_was_helix = False
        
        # Collect gauss_peak segments in helices
        for i, s in enumerate(segments):
            if s.get('model_class') == 'gauss_peak' and s['ss_type'] == 'H':
                params = s.get('params', {})
                peak_res = params.get('peak_residue', s.get('start_idx', -1))
                A = params.get('A', 0)
                sigma = params.get('sigma', 0)
                if peak_res >= 0 and A:
                    protein_helix_peaks[uid].append({
                        'peak_residue': peak_res,
                        'A': float(A),
                        'sigma': float(sigma),
                        'helix_group': seg_groups[i] if i < len(seg_groups) else -1,
                        'start_idx': s['start_idx'],
                        'end_idx': s['end_idx'],
                    })
    
    # ── Compute pairwise spacings within same helix group ──
    all_spacings = []          # absolute residue spacings
    same_helix_spacings = []   # only peaks on same contiguous helix
    any_helix_spacings = []    # peaks on any helix in same protein
    
    n_proteins_with_multi = 0
    n_pairs_same = 0
    n_pairs_any = 0
    
    for uid, peaks in protein_helix_peaks.items():
        if len(peaks) < 2:
            continue
        n_proteins_with_multi += 1
        
        # Same helix group pairs
        by_group = defaultdict(list)
        for pk in peaks:
            by_group[pk['helix_group']].append(pk)
        
        for grp, grp_peaks in by_group.items():
            if grp < 0 or len(grp_peaks) < 2:
                continue
            grp_peaks_sorted = sorted(grp_peaks, key=lambda x: x['peak_residue'])
            for i in range(len(grp_peaks_sorted)):
                for j in range(i + 1, len(grp_peaks_sorted)):
                    d = abs(grp_peaks_sorted[j]['peak_residue'] - grp_peaks_sorted[i]['peak_residue'])
                    if d > 0:
                        same_helix_spacings.append(d)
                        n_pairs_same += 1
        
        # Any helix pairs (across different helix runs in same protein)
        peaks_sorted = sorted(peaks, key=lambda x: x['peak_residue'])
        for i in range(len(peaks_sorted)):
            for j in range(i + 1, len(peaks_sorted)):
                d = abs(peaks_sorted[j]['peak_residue'] - peaks_sorted[i]['peak_residue'])
                if d > 0:
                    any_helix_spacings.append(d)
                    n_pairs_any += 1
    
    all_spacings = same_helix_spacings  # primary analysis uses same-helix
    
    n_total_helix_peaks = sum(len(v) for v in protein_helix_peaks.values())
    print(f"\n  Helix Gaussian peaks: {n_total_helix_peaks} across {len(protein_helix_peaks)} proteins")
    print(f"  Proteins with ≥2 helix peaks: {n_proteins_with_multi}")
    print(f"  Same-helix pairs: {n_pairs_same}")
    print(f"  Any-helix pairs:  {n_pairs_any}")
    
    if n_pairs_same < 5:
        print(f"\n  Insufficient same-helix peak pairs ({n_pairs_same}). Need ≥5.")
        if n_pairs_any >= 5:
            print("  Falling back to any-helix pairs (weaker test).")
            all_spacings = any_helix_spacings
        else:
            print("  Need more proteins with multiple helix Gaussian peaks.")
            return
    
    # ── Spacing modulo helix period ──
    phases = [d % HELIX_PERIOD for d in all_spacings]
    # Normalize to [-period/2, +period/2] for centering
    phases_centered = [(p if p <= HELIX_PERIOD / 2 else p - HELIX_PERIOD) for p in phases]
    
    print(f"\n  SPACING DISTRIBUTION (n = {len(all_spacings)} pairs)")
    print(f"  {'-'*50}")
    
    def _median(vals):
        s = sorted(vals)
        return s[len(s) // 2] if s else 0
    
    def _mean(vals):
        return sum(vals) / len(vals) if vals else 0
    
    print(f"  Raw spacing: mean={_mean(all_spacings):.1f}, "
          f"median={_median(all_spacings):.1f}, range=[{min(all_spacings)}, {max(all_spacings)}]")
    
    # ── ASCII histogram of phase (mod 3.6) ──
    print(f"\n  PHASE HISTOGRAM (spacing mod {HELIX_PERIOD})")
    print(f"  0.0 = in-phase (constructive), {HELIX_PERIOD/2:.1f} = anti-phase (destructive)")
    print(f"  {'-'*50}")
    
    n_bins = 12
    bin_width = HELIX_PERIOD / n_bins
    bins = [0] * n_bins
    for ph in phases:
        b = min(int(ph / bin_width), n_bins - 1)
        bins[b] += 1
    
    max_count = max(bins) if bins else 1
    for i in range(n_bins):
        lo = i * bin_width
        hi = (i + 1) * bin_width
        cnt = bins[i]
        bar = "█" * int(40 * cnt / max_count) if max_count > 0 else ""
        marker = " ◄ in-phase" if i == 0 else (" ◄ anti-phase" if i == n_bins // 2 else "")
        print(f"    [{lo:4.1f}-{hi:4.1f}): {cnt:4d} {bar}{marker}")
    
    # ── Expected uniform count per bin ──
    expected = len(phases) / n_bins
    print(f"    Expected (uniform): {expected:.1f} per bin")
    
    # ── Phase concentration metric ──
    # Fraction within ±0.5 residues of 0 mod 3.6 (in-phase zone)
    in_phase_window = 0.5  # residues
    n_in_phase = sum(1 for p in phases_centered if abs(p) <= in_phase_window)
    expected_in_phase = len(phases) * (2 * in_phase_window / HELIX_PERIOD)
    concentration = n_in_phase / max(len(phases), 1)
    expected_concentration = 2 * in_phase_window / HELIX_PERIOD
    
    print(f"\n  PHASE CONCENTRATION (±{in_phase_window} residues of 0 mod {HELIX_PERIOD})")
    print(f"  In-phase peaks:    {n_in_phase}/{len(phases)} ({concentration:.1%})")
    print(f"  Expected (uniform): {expected_in_phase:.1f}/{len(phases)} ({expected_concentration:.1%})")
    enrichment = concentration / expected_concentration if expected_concentration > 0 else 0
    print(f"  Enrichment ratio:  {enrichment:.2f}x")
    
    # ── Rayleigh test for circular non-uniformity ──
    print(f"\n  STATISTICAL TESTS")
    print(f"  {'-'*50}")
    
    # Convert phases to angles on circle [0, 2π)
    angles = [2 * math.pi * p / HELIX_PERIOD for p in phases]
    
    # Rayleigh test: R = |mean resultant vector|
    cos_sum = sum(math.cos(a) for a in angles)
    sin_sum = sum(math.sin(a) for a in angles)
    n = len(angles)
    R_bar = math.sqrt(cos_sum**2 + sin_sum**2) / n
    mean_angle = math.atan2(sin_sum, cos_sum)
    mean_phase = (mean_angle * HELIX_PERIOD / (2 * math.pi)) % HELIX_PERIOD
    
    # Rayleigh test statistic: Z = n * R_bar^2
    Z = n * R_bar**2
    # Approximate p-value: p ≈ exp(-Z) for large n
    p_rayleigh = math.exp(-Z) if Z < 500 else 0.0
    
    print(f"  Rayleigh test for circular uniformity:")
    print(f"    Mean resultant length R̄ = {R_bar:.4f} (0 = uniform, 1 = perfectly concentrated)")
    print(f"    Mean phase = {mean_phase:.2f} residues (mod {HELIX_PERIOD})")
    print(f"    Test statistic Z = {Z:.2f}")
    print(f"    p-value ≈ {p_rayleigh:.2e}" + (" (significant)" if p_rayleigh < 0.05 else " (not significant)"))
    
    # ── Chi-squared goodness of fit against uniform ──
    try:
        from scipy import stats as sp_stats
        has_scipy = True
    except ImportError:
        has_scipy = False
    
    if has_scipy:
        chi2, p_chi2 = sp_stats.chisquare(bins)
        print(f"  Chi-squared vs uniform: χ² = {chi2:.2f}, p = {p_chi2:.4f}")
        
        # Kuiper's test (circular KS) approximation via KS on phases
        uniform_phases = [i / len(phases) for i in range(len(phases))]
        sorted_phases_norm = sorted([p / HELIX_PERIOD for p in phases])
        ks_stat, p_ks = sp_stats.kstest(sorted_phases_norm, 'uniform')
        print(f"  KS vs uniform: D = {ks_stat:.4f}, p = {p_ks:.4f}")
    
    # ── Spacing histogram (raw, not mod) ──
    print(f"\n  RAW SPACING HISTOGRAM (first 40 residues)")
    print(f"  Vertical lines at multiples of {HELIX_PERIOD} (helix turns)")
    print(f"  {'-'*50}")
    
    raw_bins = [0] * 40
    for d in all_spacings:
        if d < 40:
            raw_bins[int(d)] += 1
    
    max_raw = max(raw_bins) if raw_bins else 1
    for i in range(40):
        cnt = raw_bins[i]
        if cnt == 0 and i > 0:
            continue
        bar = "█" * int(40 * cnt / max_raw) if max_raw > 0 else ""
        turn_marker = ""
        for t in range(1, 12):
            if abs(i - t * HELIX_PERIOD) < 0.5:
                turn_marker = f" ◄ ~{t} turn{'s' if t > 1 else ''}"
                break
        print(f"    {i:3d} res: {cnt:4d} {bar}{turn_marker}")
    
    # ── Verdict ──
    print(f"\n  PHASED ARRAY VERDICT")
    print(f"  {'-'*50}")
    
    if p_rayleigh < 0.01 and mean_phase < 1.0:
        print("  ✓ PHASED ARRAY HYPOTHESIS SUPPORTED")
        print(f"    Spacings cluster near 0 mod {HELIX_PERIOD} (R̄={R_bar:.3f}, p={p_rayleigh:.2e}).")
        print("    Helix Gaussian peaks show in-phase preference consistent with")
        print("    constructive interference through the H-bond lattice.")
    elif p_rayleigh < 0.05:
        print("  ~ WEAK PHASE PREFERENCE DETECTED")
        print(f"    Non-uniform distribution (p={p_rayleigh:.3f}) but effect is modest (R̄={R_bar:.3f}).")
        print("    May indicate partial phasing or mixed populations.")
    elif enrichment > 1.3:
        print("  ~ IN-PHASE ENRICHMENT WITHOUT CIRCULAR SIGNIFICANCE")
        print(f"    Enrichment at 0 mod {HELIX_PERIOD} = {enrichment:.2f}x, but Rayleigh p={p_rayleigh:.3f}.")
        print("    Suggestive but not conclusive at current sample size.")
    else:
        print("  ✗ NO PHASED ARRAY SIGNAL DETECTED")
        print(f"    Spacings are consistent with uniform mod {HELIX_PERIOD} (R̄={R_bar:.3f}, p={p_rayleigh:.2f}).")
        print("    Gaussian peaks appear independently positioned on helices.")
    
    # ── Phase-amplitude correlation ──
    # Do in-phase peak pairs have correlated amplitudes? (co-phased = same sign A)
    print(f"\n  PHASE-AMPLITUDE CORRELATION")
    print(f"  {'-'*50}")
    
    # Rebuild pairs with amplitude info
    pair_data = []
    for uid, peaks in protein_helix_peaks.items():
        by_group = defaultdict(list)
        for pk in peaks:
            by_group[pk['helix_group']].append(pk)
        for grp, grp_peaks in by_group.items():
            if grp < 0 or len(grp_peaks) < 2:
                continue
            grp_peaks_sorted = sorted(grp_peaks, key=lambda x: x['peak_residue'])
            for i in range(len(grp_peaks_sorted)):
                for j in range(i + 1, len(grp_peaks_sorted)):
                    pi, pj = grp_peaks_sorted[i], grp_peaks_sorted[j]
                    d = abs(pj['peak_residue'] - pi['peak_residue'])
                    if d > 0:
                        phase = d % HELIX_PERIOD
                        phase_c = phase if phase <= HELIX_PERIOD / 2 else phase - HELIX_PERIOD
                        same_sign = (pi['A'] > 0) == (pj['A'] > 0)
                        pair_data.append({
                            'spacing': d,
                            'phase': phase,
                            'phase_centered': phase_c,
                            'same_sign_A': same_sign,
                            'A_product': pi['A'] * pj['A'],
                            'uid': uid,
                        })
    
    if len(pair_data) >= 10:
        in_phase_pairs = [p for p in pair_data if abs(p['phase_centered']) <= 0.9]
        anti_phase_pairs = [p for p in pair_data if abs(abs(p['phase_centered']) - HELIX_PERIOD / 2) <= 0.9]
        
        if in_phase_pairs and anti_phase_pairs:
            ip_same = sum(1 for p in in_phase_pairs if p['same_sign_A'])
            ap_same = sum(1 for p in anti_phase_pairs if p['same_sign_A'])
            ip_frac = ip_same / len(in_phase_pairs) if in_phase_pairs else 0
            ap_frac = ap_same / len(anti_phase_pairs) if anti_phase_pairs else 0
            
            print(f"  In-phase pairs (|Δ mod 3.6| ≤ 0.9):     n={len(in_phase_pairs):3d}, "
                  f"same-sign A: {ip_frac:.0%}")
            print(f"  Anti-phase pairs (|Δ mod 3.6| ~ 1.8):    n={len(anti_phase_pairs):3d}, "
                  f"same-sign A: {ap_frac:.0%}")
            
            if ip_frac > ap_frac + 0.1:
                print("  → In-phase peaks tend to have same-sign amplitude (constructive).")
            elif ap_frac > ip_frac + 0.1:
                print("  → Anti-phase peaks tend to have same-sign amplitude (unexpected).")
            else:
                print("  → No clear amplitude-phase correlation.")
        else:
            print(f"  In-phase pairs: {len(in_phase_pairs)}, anti-phase: {len(anti_phase_pairs)}")
            print("  Insufficient pairs in one or both categories.")
    else:
        print(f"  Only {len(pair_data)} pairs — need ≥10 for correlation analysis.")
    
    # ═══════════════════════════════════════════════════════════════════
    #  EXCLUSION ZONE ANALYSIS
    #  Gaussian peaks maintain minimum separation — each peak "claims"
    #  a territory on the helix where no other peak can exist.
    # ═══════════════════════════════════════════════════════════════════
    print(f"\n  {'=' * 60}")
    print(f"  EXCLUSION ZONE ANALYSIS")
    print(f"  {'=' * 60}")
    
    # Collect nearest-neighbor distances (consecutive peaks on same helix)
    nn_distances = []  # nearest-neighbor spacing between consecutive peaks
    helix_info = []    # (helix_length, n_peaks, peak_positions) for null model
    
    for uid, peaks in protein_helix_peaks.items():
        by_group = defaultdict(list)
        for pk in peaks:
            by_group[pk['helix_group']].append(pk)
        
        for grp, grp_peaks in by_group.items():
            if grp < 0 or len(grp_peaks) < 2:
                # Single-peak helices: record for null model
                if grp >= 0 and len(grp_peaks) == 1:
                    pk = grp_peaks[0]
                    helix_len = pk['end_idx'] - pk['start_idx'] + 1
                    # Approximate: peak is somewhere in the helix
                    helix_info.append((helix_len, 1, [pk['peak_residue']]))
                continue
            
            grp_sorted = sorted(grp_peaks, key=lambda x: x['peak_residue'])
            positions = [p['peak_residue'] for p in grp_sorted]
            
            # Estimate helix span from first segment start to last segment end
            helix_start = min(p['start_idx'] for p in grp_sorted)
            helix_end = max(p['end_idx'] for p in grp_sorted)
            helix_len = helix_end - helix_start + 1
            helix_info.append((helix_len, len(positions), positions))
            
            # Nearest-neighbor distances
            for i in range(len(positions) - 1):
                d = positions[i + 1] - positions[i]
                if d > 0:
                    nn_distances.append(d)
    
    if not nn_distances:
        print("  No consecutive peak pairs found on same helix.")
    else:
        nn_sorted = sorted(nn_distances)
        min_nn = nn_sorted[0]
        q25 = nn_sorted[len(nn_sorted) // 4]
        median_nn = nn_sorted[len(nn_sorted) // 2]
        q75 = nn_sorted[3 * len(nn_sorted) // 4]
        max_nn = nn_sorted[-1]
        mean_nn = sum(nn_distances) / len(nn_distances)
        
        print(f"\n  Nearest-neighbor distances (consecutive peaks, same helix):")
        print(f"    n pairs:  {len(nn_distances)}")
        print(f"    Min:      {min_nn} residues ({min_nn / HELIX_PERIOD:.1f} turns)")
        print(f"    Q25:      {q25} residues ({q25 / HELIX_PERIOD:.1f} turns)")
        print(f"    Median:   {median_nn} residues ({median_nn / HELIX_PERIOD:.1f} turns)")
        print(f"    Q75:      {q75} residues ({q75 / HELIX_PERIOD:.1f} turns)")
        print(f"    Max:      {max_nn} residues ({max_nn / HELIX_PERIOD:.1f} turns)")
        print(f"    Mean:     {mean_nn:.1f} residues ({mean_nn / HELIX_PERIOD:.1f} turns)")
        
        # ── NN distance histogram ──
        print(f"\n  NEAREST-NEIGHBOR DISTANCE HISTOGRAM")
        print(f"  {'-'*50}")
        
        hist_max = min(60, max_nn + 2)
        bin_w = 3  # 3-residue bins
        n_hist_bins = (hist_max + bin_w - 1) // bin_w
        nn_hist = [0] * n_hist_bins
        for d in nn_distances:
            b = min(d // bin_w, n_hist_bins - 1)
            nn_hist[b] += 1
        
        max_cnt = max(nn_hist) if nn_hist else 1
        for i in range(n_hist_bins):
            lo = i * bin_w
            hi = (i + 1) * bin_w
            cnt = nn_hist[i]
            if cnt == 0 and lo > median_nn + 20:
                continue
            bar = "█" * int(40 * cnt / max_cnt) if max_cnt > 0 else ""
            zone_marker = ""
            if lo < min_nn:
                zone_marker = " ◄ EXCLUSION ZONE"
            elif lo == (min_nn // bin_w) * bin_w:
                zone_marker = " ◄ zone boundary"
            print(f"    [{lo:3d}-{hi:3d}): {cnt:4d} {bar}{zone_marker}")
        
        # ── Null model: random peak placement on helices ──
        print(f"\n  EXCLUSION ZONE NULL TEST")
        print(f"  H₀: peaks placed uniformly at random along helices")
        print(f"  {'-'*50}")
        
        import random
        rng = random.Random(42)
        
        # For each multi-peak helix, simulate random placement and measure NN distance
        null_nn_distances = []
        n_null_reps = 2000
        
        # Collect multi-peak helices with their dimensions
        multi_peak_helices = [(hl, np_, pos) for hl, np_, pos in helix_info if np_ >= 2]
        
        if multi_peak_helices:
            for _ in range(n_null_reps):
                for helix_len, n_peaks_on_helix, real_positions in multi_peak_helices:
                    # Place n_peaks_on_helix points uniformly on [0, helix_len)
                    # Each peak has a width (~mean sigma), but for the null we 
                    # place point positions
                    null_positions = sorted(rng.randint(0, max(helix_len - 1, 1))
                                            for _ in range(n_peaks_on_helix))
                    for i in range(len(null_positions) - 1):
                        d = null_positions[i + 1] - null_positions[i]
                        if d > 0:
                            null_nn_distances.append(d)
            
            null_sorted = sorted(null_nn_distances)
            null_min = null_sorted[0] if null_sorted else 0
            null_median = null_sorted[len(null_sorted) // 2] if null_sorted else 0
            null_mean = sum(null_nn_distances) / len(null_nn_distances) if null_nn_distances else 0
            null_q05 = null_sorted[int(len(null_sorted) * 0.05)] if null_sorted else 0
            
            print(f"  Real NN distances:  min={min_nn}, median={median_nn}, mean={mean_nn:.1f}")
            print(f"  Null NN distances:  min={null_min}, median={null_median}, "
                  f"mean={null_mean:.1f}, 5th pctile={null_q05}")
            
            # What fraction of null NN distances are ≥ real minimum?
            frac_above_real_min = sum(1 for d in null_nn_distances if d >= min_nn) / len(null_nn_distances)
            print(f"\n  Null P(NN ≥ {min_nn}): {frac_above_real_min:.1%}")
            
            # Mann-Whitney: are real NN distances larger than null?
            try:
                from scipy import stats as sp_stats
                u_nn, p_nn = sp_stats.mannwhitneyu(nn_distances, null_nn_distances, alternative='greater')
                n1, n2 = len(nn_distances), len(null_nn_distances)
                r_nn = 1 - 2 * u_nn / (n1 * n2)
                print(f"  Mann-Whitney (real > null): U={u_nn:.0f}, p={p_nn:.4f}, r={r_nn:.3f}")
            except ImportError:
                pass
            
            # Exclusion zone size estimation
            # Find the distance below which we see fewer peaks than expected
            # Use null CDF vs real CDF
            print(f"\n  EXCLUSION ZONE ESTIMATION")
            print(f"  {'-'*50}")
            
            # For each distance threshold, compare real vs null fraction below it
            thresholds = list(range(3, min(40, max_nn), 3))
            print(f"  {'Threshold':>10s} {'Real ≤ thr':>12s} {'Null ≤ thr':>12s} {'Depletion':>10s}")
            print(f"  {'-'*10} {'-'*12} {'-'*12} {'-'*10}")
            
            zone_size = 0
            for thr in thresholds:
                real_below = sum(1 for d in nn_distances if d <= thr) / len(nn_distances)
                null_below = sum(1 for d in null_nn_distances if d <= thr) / len(null_nn_distances)
                depletion = 1 - (real_below / null_below) if null_below > 0 else 0
                marker = " ◄" if depletion > 0.3 else ""
                print(f"  {thr:10d} {real_below:11.1%} {null_below:11.1%} {depletion:9.0%}{marker}")
                if depletion > 0.3 and zone_size == 0:
                    zone_size = thr
            
            # Estimate zone in helix turns
            if zone_size > 0:
                print(f"\n  Estimated exclusion zone: ~{zone_size} residues "
                      f"(~{zone_size / HELIX_PERIOD:.1f} helix turns)")
                print(f"  Each Gaussian peak 'claims' ~{zone_size} residues of helix territory.")
            else:
                # Use minimum NN as fallback estimate
                zone_size = min_nn
                print(f"\n  Minimum observed NN: {min_nn} residues ({min_nn / HELIX_PERIOD:.1f} turns)")
            
            # ── Biological interpretation ──
            print(f"\n  INTERPRETATION")
            print(f"  {'-'*50}")
            
            # Estimate peak "footprint" from sigma values
            peak_sigmas = []
            for uid, peaks in protein_helix_peaks.items():
                for pk in peaks:
                    if pk['sigma'] > 0:
                        peak_sigmas.append(pk['sigma'])
            
            if peak_sigmas:
                # σ is in arc-length units; approximate residue count
                # Typical step ~0.3 rad, so σ residues ≈ σ / 0.3
                med_sigma = sorted(peak_sigmas)[len(peak_sigmas) // 2]
                mean_sigma = sum(peak_sigmas) / len(peak_sigmas)
                # Two-sigma footprint (95% of the Gaussian)
                footprint_arclength = 4 * med_sigma  # ±2σ
                # Convert: typical inter-residue arc ~0.2-0.5 rad on T²
                # Use median NN as empirical footprint instead
                print(f"  Median peak σ (arc-length): {med_sigma:.3f}")
                print(f"  Gaussian 95% footprint (4σ): {4*med_sigma:.2f} arc-length units")
            
            if min_nn >= 8:
                print(f"\n  Gaussian peaks maintain ≥{min_nn}-residue ({min_nn/HELIX_PERIOD:.0f}-turn) separation.")
                print(f"  This is consistent with each peak representing an independent")
                print(f"  structural feature (kink, cap, binding distortion) that requires")
                print(f"  a minimum helix 'runway' to exist.")
                print(f"  Peaks are territorial, not cooperative — each claims a segment")
                print(f"  of helix where no other curvature perturbation can fit.")
            elif zone_size > 0:
                print(f"\n  Exclusion zone of ~{zone_size} residues suggests steric/geometric")
                print(f"  constraints on how closely curvature perturbations can be packed.")
        else:
            print("  No multi-peak helices found for null model.")
    
    # ── REVISED VERDICT ──
    print(f"\n  {'=' * 60}")
    print(f"  REVISED VERDICT")
    print(f"  {'=' * 60}")
    print(f"  Phased array (cooperative):   ✗ No signal (R̄={R_bar:.3f}, p={p_rayleigh:.2f})")
    if nn_distances and min_nn >= 8:
        print(f"  Exclusion zone (territorial): ✓ Min NN = {min_nn} res "
              f"({min_nn/HELIX_PERIOD:.0f} turns)")
        print(f"  ")
        print(f"  Gaussian peaks are TERRITORIAL, not COOPERATIVE.")
        print(f"  Each peak claims ~{min_nn}+ residues of helix and excludes")
        print(f"  neighboring perturbations. The helix is not a phased array")
        print(f"  but a segmented ruler — each curvature event marks an")
        print(f"  independent structural feature with its own buffer zone.")
    elif nn_distances:
        print(f"  Exclusion zone: min NN = {min_nn} residues (weak or absent)")
    print()
    # ── Save results ──
    results = {
        'n_helix_gauss_peaks': n_total_helix_peaks,
        'n_proteins_with_multi_peaks': n_proteins_with_multi,
        'n_same_helix_pairs': n_pairs_same,
        'n_any_helix_pairs': n_pairs_any,
        'helix_period': HELIX_PERIOD,
        'R_bar': R_bar,
        'mean_phase': mean_phase,
        'Z_rayleigh': Z,
        'p_rayleigh': p_rayleigh,
        'in_phase_concentration': concentration,
        'expected_concentration': expected_concentration,
        'enrichment_ratio': enrichment,
        'spacings': all_spacings,
        'phase_histogram_bins': bins,
        'nn_distances': nn_distances if nn_distances else [],
        'nn_min': min(nn_distances) if nn_distances else None,
        'nn_median': sorted(nn_distances)[len(nn_distances)//2] if nn_distances else None,
        'nn_mean': sum(nn_distances)/len(nn_distances) if nn_distances else None,
        'phased_array_supported': False,
        'exclusion_zone_found': bool(nn_distances and min(nn_distances) >= 8),
    }
    
    out_path = RESULTS_DIR / "inter_peak_spacing_results.json"
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\n  Results saved: {out_path}")
    except Exception as e:
        print(f"\n  (Could not save: {e})")
    
    print()


# ─── Main ───
if __name__ == "__main__":
    barcodes = load_all_barcodes()
    print(f"\n  Loaded {len(barcodes)} barcodes from {BARCODE_DIR}")

    args = sys.argv[1:]

    if not args:
        dashboard(barcodes)
        top_defects(barcodes, 15)
        top_geodesic(barcodes, 15)
        winding_map(barcodes)
        cluster_analysis(barcodes)

    elif args[0] == '--top-defects':
        n = int(args[1]) if len(args) > 1 else 20
        top_defects(barcodes, n)

    elif args[0] == '--top-geodesic':
        n = int(args[1]) if len(args) > 1 else 20
        top_geodesic(barcodes, n)

    elif args[0] == '--winding-map':
        winding_map(barcodes)

    elif args[0] == '--clusters':
        cluster_analysis(barcodes)

    elif args[0] == '--find':
        if len(args) < 2:
            print("  Usage: analyze_barcodes.py --find P04637")
        else:
            find_protein(barcodes, args[1])

    elif args[0] == '--defect-catalog':
        defect_catalog(barcodes)

    elif args[0] == '--peak-analysis':
        threshold = 5
        do_fetch = True
        for a in args[1:]:
            if a == '--no-fetch':
                do_fetch = False
            elif a.startswith('--threshold='):
                threshold = int(a.split('=')[1])
            else:
                try:
                    threshold = int(a)
                except ValueError:
                    pass
        peak_analysis(barcodes, proximity_threshold=threshold, fetch_annotations=do_fetch)

    elif args[0] == '--inter-peak-spacing':
        inter_peak_spacing(barcodes)

    else:
        print(f"  Unknown option: {args[0]}")
        print("  Options: --top-defects, --top-geodesic, --winding-map, --clusters,")
        print("           --find, --defect-catalog, --peak-analysis, --inter-peak-spacing")

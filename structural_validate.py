#!/usr/bin/env python3
"""
Structural validation via RCSB Pairwise Alignment Tool.

Opens browser tabs for each pair at rcsb.org/alignment where you can
upload PDB files and get TM-scores instantly.

Also provides direct Foldseek search links.

Usage:
    python structural_validate.py           # Print URLs + open guide
    python structural_validate.py results   # Enter TM-scores manually
    python structural_validate.py report    # Show report + generate fig
"""

import json
import sys
import webbrowser
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
PDB_DIR = SCRIPT_DIR / "dali_pdb"
PAIR_DIR = SCRIPT_DIR / "dali_pairs"
RESULTS_DIR = SCRIPT_DIR / "dali_results"

PAIRS = [
    {"id": 1,  "query": "P00918", "hit": "P59991", "q_name": "CA2",     "h_name": "KRTAP12-2",  "geom": 0.804},
    {"id": 2,  "query": "P60174", "hit": "P35270", "q_name": "TPI",     "h_name": "SPR",        "geom": 0.798},
    {"id": 3,  "query": "P62258", "hit": "O43808", "q_name": "14-3-3e", "h_name": "SLC25A17",   "geom": 0.783},
    {"id": 4,  "query": "P11166", "hit": "P09914", "q_name": "GLUT1",   "h_name": "IFIT1",      "geom": 0.783},
    {"id": 5,  "query": "P02751", "hit": "P21333", "q_name": "FN1",     "h_name": "FLNA",       "geom": 0.779},
    {"id": 6,  "query": "P01857", "hit": "P32969", "q_name": "IgG1",    "h_name": "RPL9",       "geom": 0.778},
    {"id": 7,  "query": "P38646", "hit": "O15344", "q_name": "GRP75",   "h_name": "MID1",       "geom": 0.767},
    {"id": 8,  "query": "P07900", "hit": "O14829", "q_name": "HSP90a",  "h_name": "PPEF1",      "geom": 0.765},
    {"id": 9,  "query": "P04626", "hit": "O95302", "q_name": "HER2",    "h_name": "FKBP9",      "geom": 0.757},
    {"id": 10, "query": "P21860", "hit": "P28799", "q_name": "ErbB3",   "h_name": "GRN",        "geom": 0.757},
]


def cmd_guide():
    """Print step-by-step guide for structural validation."""
    print()
    print("  ╔══════════════════════════════════════════════════════════════════╗")
    print("  ║        STRUCTURAL VALIDATION — STEP BY STEP                     ║")
    print("  ╚══════════════════════════════════════════════════════════════════╝")
    print()
    print("  METHOD: RCSB PDB Comparison Tool (uses TM-align, CE, FATCAT)")
    print("  URL:    https://www.rcsb.org/alignment")
    print()
    print("  FOR EACH PAIR:")
    print("  1. Open https://www.rcsb.org/alignment")
    print("  2. Left side  → 'Upload file' → pick structure1_*.pdb")
    print("  3. Right side → 'Upload file' → pick structure2_*.pdb")
    print("  4. Select algorithm: 'TM-align' (best for fold comparison)")
    print("  5. Click 'Submit'")
    print("  6. Record the TM-score (0-1 scale)")
    print("     TM > 0.5 = same fold   |   TM > 0.7 = very similar")
    print("  7. Also note RMSD and alignment length if shown")
    print()
    print("  YOUR 10 PAIRS:")
    print("  Each pair folder has structure1_*.pdb and structure2_*.pdb")
    print()

    for p in PAIRS:
        pair_folder = None
        for d in PAIR_DIR.iterdir() if PAIR_DIR.exists() else []:
            if d.is_dir() and f"pair{p['id']:02d}" in d.name:
                pair_folder = d
                break

        files = ""
        if pair_folder:
            pdbs = sorted(pair_folder.glob("*.pdb"))
            if len(pdbs) >= 2:
                files = f"\n         File 1: {pdbs[0].name}\n         File 2: {pdbs[1].name}"

        print(f"  Pair {p['id']:2d}: {p['q_name']:>8s}  vs  {p['h_name']:>10s}  "
              f"(geom={p['geom']:.3f}){files}")
        print()

    print("  ALTERNATIVE: Foldseek web server")
    print("  URL: https://search.foldseek.com/search")
    print("  Upload structure1_*.pdb, search against AlphaFold Swiss-Prot")
    print("  Check if structure2's protein appears in results")
    print()
    print("  After recording TM-scores:")
    print("  Run: python structural_validate.py results")


def cmd_results():
    """Enter TM-scores manually."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    results_file = RESULTS_DIR / "structural_validation.json"

    if results_file.exists():
        with open(results_file) as f:
            existing = json.load(f)
    else:
        existing = {}

    print()
    print("  Enter TM-scores for each pair (from RCSB or Foldseek).")
    print("  TM > 0.5 = same fold  |  TM > 0.7 = very similar  |  TM < 0.3 = unrelated")
    print("  Press Enter to skip, 'q' to quit.")
    print("  You can also enter DALI Z-scores (values > 2 detected automatically)")
    print()

    for p in PAIRS:
        pid = str(p["id"])
        existing_tm = existing.get(pid, {}).get("tm_score")
        existing_z = existing.get(pid, {}).get("dali_z")
        note = ""
        if existing_tm is not None:
            note += f"  [TM={existing_tm}]"
        if existing_z is not None:
            note += f"  [Z={existing_z}]"

        prompt = (f"  Pair {p['id']:2d}: {p['q_name']:>8s} → {p['h_name']:>10s} "
                  f"(geom={p['geom']:.3f}){note}\n"
                  f"         TM-score (or Z-score) = ")
        try:
            val = input(prompt).strip()
        except (EOFError, KeyboardInterrupt):
            print("\n  Stopped.")
            break

        if val.lower() == 'q':
            break
        if not val:
            continue

        try:
            score = float(val)

            # Auto-detect if it's a TM-score (0-1) or DALI Z-score (>1)
            if score <= 1.0:
                tm = score
                z = None
                sig = "same fold" if tm > 0.5 else "unrelated"
            else:
                z = score
                tm = None
                sig = "significant" if z > 2 else "not significant"

            existing[pid] = {
                "pair_id": p["id"],
                "query": p["query"],
                "hit": p["hit"],
                "q_name": p["q_name"],
                "h_name": p["h_name"],
                "geom_score": p["geom"],
                "tm_score": tm,
                "dali_z": z,
                "same_fold": (tm is not None and tm > 0.5) or (z is not None and z > 2),
            }
            print(f"         → {sig}")
        except ValueError:
            # Maybe they entered RMSD or other text
            print(f"         → Invalid number, skipped")

    with open(results_file, 'w') as f:
        json.dump(existing, f, indent=2)
    print(f"\n  Saved to {results_file}")
    cmd_report()


def cmd_report():
    """Show results and generate figure."""
    results_file = RESULTS_DIR / "structural_validation.json"

    if not results_file.exists():
        print("  No results yet. Run: python structural_validate.py results")
        return

    with open(results_file) as f:
        results = json.load(f)

    if not results:
        print("  No scores recorded yet.")
        return

    print()
    print("  ╔══════════════════════════════════════════════════════════════════╗")
    print("  ║        STRUCTURAL VALIDATION RESULTS                            ║")
    print("  ╚══════════════════════════════════════════════════════════════════╝")
    print()
    print(f"  {'#':>3s} {'Query':>8s} → {'Hit':>10s}  {'Geom':>5s}  {'TM':>5s}  {'Z':>5s}  {'Fold?':>5s}")
    print(f"  {'─'*55}")

    n_same_fold = 0
    n_total = 0
    geom_scores = []
    validation_scores = []  # TM or normalized Z

    for p in PAIRS:
        pid = str(p["id"])
        r = results.get(pid, {})
        tm = r.get("tm_score")
        z = r.get("dali_z")
        same = r.get("same_fold", False)

        tm_str = f"{tm:.3f}" if tm is not None else "  -  "
        z_str = f"{z:.1f}" if z is not None else "  -  "
        fold_str = " YES " if same else " no  " if (tm is not None or z is not None) else "  ?  "

        print(f"  {p['id']:3d} {p['q_name']:>8s} → {p['h_name']:>10s}  "
              f"{p['geom']:5.3f}  {tm_str}  {z_str}  {fold_str}")

        if tm is not None or z is not None:
            n_total += 1
            if same:
                n_same_fold += 1
            geom_scores.append(p["geom"])
            # Use TM if available, else normalize Z to 0-1 range
            if tm is not None:
                validation_scores.append(tm)
            elif z is not None:
                validation_scores.append(min(z / 20.0, 1.0))  # Rough normalization

    if n_total > 0:
        print()
        print(f"  SUMMARY:")
        print(f"    Pairs validated:        {n_total}")
        print(f"    Same fold (TM>0.5):     {n_same_fold}/{n_total} ({n_same_fold/n_total*100:.0f}%)")
        print(f"    Different fold:         {n_total - n_same_fold}/{n_total}")

        if n_total >= 3:
            import numpy as np
            corr = np.corrcoef(geom_scores, validation_scores)[0, 1]
            print(f"    Correlation (geom↔TM):  r = {corr:.3f}")

        # Generate figure
        try:
            _plot_validation(geom_scores, validation_scores, results)
        except Exception as e:
            print(f"    (Figure failed: {e})")

    # Save for paper
    paper_results = []
    for p in PAIRS:
        pid = str(p["id"])
        r = results.get(pid, {})
        if r.get("tm_score") is not None or r.get("dali_z") is not None:
            paper_results.append({
                "query": p["query"],
                "hit": p["hit"],
                "query_gene": p["q_name"],
                "hit_gene": p["h_name"],
                "geom_score": p["geom"],
                "tm_score": r.get("tm_score"),
                "dali_z": r.get("dali_z"),
                "same_fold": r.get("same_fold", False),
            })

    from pathlib import Path
    data_dir = SCRIPT_DIR / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    with open(data_dir / "structural_validation_results.json", 'w') as f:
        json.dump(paper_results, f, indent=2)


def _plot_validation(geom_scores, val_scores, results):
    """Generate validation scatter plot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))

    colors = ['#27AE60' if v > 0.5 else '#E74C3C' for v in val_scores]
    ax.scatter(geom_scores, val_scores, c=colors, s=80, edgecolors='black',
               lw=0.5, zorder=5)

    # Labels
    for p in PAIRS:
        pid = str(p["id"])
        r = results.get(pid, {})
        tm = r.get("tm_score")
        z = r.get("dali_z")
        val = tm if tm is not None else (min(z/20, 1.0) if z is not None else None)
        if val is not None:
            label = f"{p['q_name']}-{p['h_name']}"
            ax.annotate(label, (p["geom"], val), textcoords='offset points',
                       xytext=(5, 5), fontsize=7)

    ax.axhline(0.5, color='#E74C3C', ls='--', alpha=0.4, label='TM = 0.5 (same fold)')
    ax.axhline(0.7, color='#27AE60', ls='--', alpha=0.4, label='TM = 0.7 (very similar)')

    ax.set_xlabel('Geometric similarity score', fontsize=11)
    ax.set_ylabel('TM-score (structural validation)', fontsize=11)
    ax.set_title('Geometry-only hits: 3D structural validation', fontsize=12)
    ax.legend(fontsize=9, loc='upper left')
    ax.set_xlim(0.74, 0.82)
    ax.set_ylim(0, 1.05)

    if len(geom_scores) >= 3:
        corr = np.corrcoef(geom_scores, val_scores)[0, 1]
        ax.text(0.95, 0.05, f'r = {corr:.2f}', transform=ax.transAxes,
                fontsize=10, ha='right', va='bottom',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    fig.tight_layout()

    fig_dir = SCRIPT_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / 'fig7_structural_validation.svg', bbox_inches='tight')
    fig.savefig(fig_dir / 'fig7_structural_validation.png', dpi=300, bbox_inches='tight')
    print(f"    Figure saved: figures/fig7_structural_validation.svg/png")
    plt.close(fig)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "guide"

    if cmd in ("guide", "help"):
        cmd_guide()
    elif cmd in ("results", "enter", "manual"):
        cmd_results()
    elif cmd in ("report",):
        cmd_report()
    else:
        print(f"  Unknown: {cmd}. Try: guide, results, report")

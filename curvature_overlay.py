#!/usr/bin/env python3
"""
Curvature Profile Overlay — Three Validation Pairs
====================================================

Extracts per-residue κ² profiles from AlphaFold CIF files and generates
overlay plots showing backbone curvature alignment between:

  1. Human GBA1 (P04062) ↔ Soybean Q39817  [Gaucher disease]
  2. Human TM2D3 (Q9BRN9) ↔ E. coli P39310 [Alzheimer's disease]
  3. Human COX11 (Q9Y6N1) ↔ Pseudomonas O86062 [Mitochondrial disease]

Each plot shows the κ² trace for both proteins on the same axes,
with the x-axis normalized to relative position (0-1) so different-length
proteins can be directly compared.

Usage:
  python curvature_overlay.py --cache alphafold_cache --output overlay_plots

Requires: matplotlib, numpy
  pip install matplotlib numpy
"""

import os
import sys
import math
import argparse
import glob
import numpy as np

# ═══════════════════════════════════════════════════════════════
# CIF PARSER — Extract Cα coordinates
# ═══════════════════════════════════════════════════════════════

def parse_cif_ca(path):
    """Extract Cα coordinates from mmCIF file. Returns list of (x,y,z) and pLDDT."""
    coords = []
    plddts = []
    
    with open(path) as f:
        in_atom = False
        for line in f:
            if line.startswith("_atom_site."):
                in_atom = True
                continue
            if in_atom and not line.startswith("_") and not line.startswith("#") and line.strip():
                parts = line.split()
                if len(parts) < 10:
                    if line.strip() == "#":
                        in_atom = False
                    continue
                
                # mmCIF atom_site columns (AlphaFold format):
                # 0:group 1:id 2:type_symbol 3:label_atom_id 4:label_alt_id
                # 5:label_comp_id 6:label_asym_id 7:label_entity_id 8:label_seq_id
                # 9:Cartn_x 10:Cartn_y 11:Cartn_z 12:occupancy 13:B_iso
                # 14:auth_seq_id 15:auth_asym_id
                
                try:
                    atom_name = parts[3].strip('"')
                    if atom_name != "CA":
                        continue
                    
                    x = float(parts[9])
                    y = float(parts[10])
                    z = float(parts[11])
                    b_factor = float(parts[13])  # pLDDT in AlphaFold
                    
                    coords.append((x, y, z))
                    plddts.append(b_factor)
                except (ValueError, IndexError):
                    continue
            elif in_atom and (line.startswith("#") or line.startswith("_") or line.startswith("loop_")):
                in_atom = False
    
    return np.array(coords), np.array(plddts)


# ═══════════════════════════════════════════════════════════════
# CURVATURE COMPUTATION — Frenet-Serret κ²
# ═══════════════════════════════════════════════════════════════

def compute_kappa2(coords):
    """Compute discrete backbone curvature κ² at each residue using
    Frenet-Serret formalism on Cα coordinates.
    
    κ² = |T'(s)|² where T is the unit tangent vector.
    Discretized as: κ²_i = |T_{i+1} - T_{i-1}|² / |Δs|²
    """
    n = len(coords)
    if n < 3:
        return np.zeros(n)
    
    kappa2 = np.zeros(n)
    
    # Compute tangent vectors
    tangents = np.zeros((n, 3))
    for i in range(n - 1):
        diff = coords[i + 1] - coords[i]
        norm = np.linalg.norm(diff)
        if norm > 0:
            tangents[i] = diff / norm
    tangents[-1] = tangents[-2]
    
    # Compute curvature from tangent changes
    for i in range(1, n - 1):
        dt = tangents[i] - tangents[i - 1]
        ds = np.linalg.norm(coords[i + 1] - coords[i - 1]) / 2.0
        if ds > 0:
            kappa2[i] = np.dot(dt, dt) / (ds * ds)
    
    return kappa2


# ═══════════════════════════════════════════════════════════════
# SECONDARY STRUCTURE ASSIGNMENT (simplified)
# ═══════════════════════════════════════════════════════════════

def assign_ss(coords):
    """Simple SS assignment based on Cα distances.
    Helix: i to i+3 distance ~5.0-5.5 Å
    Sheet: i to i+1 distance ~3.3-3.5 Å and extended
    """
    n = len(coords)
    ss = ['L'] * n  # L=loop, H=helix, E=strand
    
    for i in range(n - 3):
        d = np.linalg.norm(coords[i + 3] - coords[i])
        if 4.8 < d < 5.8:
            for j in range(i, min(i + 4, n)):
                if ss[j] == 'L':
                    ss[j] = 'H'
    
    # Simple strand detection: extended backbone
    for i in range(1, n - 1):
        d_prev = np.linalg.norm(coords[i] - coords[i - 1])
        d_next = np.linalg.norm(coords[i + 1] - coords[i])
        if 3.6 < d_prev < 4.0 and 3.6 < d_next < 4.0 and ss[i] == 'L':
            ss[i] = 'E'
    
    return ss


# ═══════════════════════════════════════════════════════════════
# PROFILE STATISTICS
# ═══════════════════════════════════════════════════════════════

def profile_stats(kappa2, label=""):
    """Compute summary statistics for a κ² profile."""
    valid = kappa2[kappa2 > 0]
    if len(valid) == 0:
        return {}
    
    sorted_k = np.sort(valid)
    n = len(sorted_k)
    cumsum = np.cumsum(sorted_k)
    total = cumsum[-1]
    
    # Gini coefficient
    gini = 1.0 - 2.0 * np.sum(cumsum) / (n * total) + 1.0 / n if total > 0 else 0
    
    # Peak
    peak_idx = np.argmax(kappa2)
    peak_val = kappa2[peak_idx]
    mean_val = np.mean(valid)
    concentration = peak_val / mean_val if mean_val > 0 else 0
    
    return {
        "label": label,
        "n_residues": len(kappa2),
        "gini": gini,
        "mean": mean_val,
        "median": np.median(valid),
        "max": peak_val,
        "peak_residue": peak_idx,
        "concentration": concentration,
    }


# ═══════════════════════════════════════════════════════════════
# CORRELATION ANALYSIS
# ═══════════════════════════════════════════════════════════════

def profile_correlation(k1, k2):
    """Compute correlation between two κ² profiles after length normalization."""
    # Interpolate both to 1000 points
    n_points = 1000
    x1 = np.linspace(0, 1, len(k1))
    x2 = np.linspace(0, 1, len(k2))
    x_common = np.linspace(0, 1, n_points)
    
    k1_interp = np.interp(x_common, x1, k1)
    k2_interp = np.interp(x_common, x2, k2)
    
    # Pearson correlation
    if np.std(k1_interp) > 0 and np.std(k2_interp) > 0:
        pearson = np.corrcoef(k1_interp, k2_interp)[0, 1]
    else:
        pearson = 0
    
    # Log-space correlation (better for high-dynamic-range data)
    lk1 = np.log1p(k1_interp)
    lk2 = np.log1p(k2_interp)
    if np.std(lk1) > 0 and np.std(lk2) > 0:
        log_pearson = np.corrcoef(lk1, lk2)[0, 1]
    else:
        log_pearson = 0
    
    # Rank correlation (Spearman)
    from scipy.stats import spearmanr
    try:
        spearman, sp_pval = spearmanr(k1_interp, k2_interp)
    except:
        spearman, sp_pval = 0, 1
    
    return {
        "pearson": pearson,
        "log_pearson": log_pearson,
        "spearman": spearman,
        "spearman_p": sp_pval,
    }


# ═══════════════════════════════════════════════════════════════
# PLOTTING
# ═══════════════════════════════════════════════════════════════

def plot_overlay(k1, k2, plddt1, plddt2, ss1, ss2, stats1, stats2, 
                 corr, title, label1, label2, outpath):
    """Generate publication-quality overlay plot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.gridspec import GridSpec
    
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(4, 2, height_ratios=[3, 1, 1, 0.5], hspace=0.3, wspace=0.3)
    
    # Color scheme
    c1 = '#2166AC'  # Blue for human
    c2 = '#B2182B'  # Red for match
    c1_light = '#92C5DE'
    c2_light = '#F4A582'
    
    # ─── Main overlay plot (normalized position) ───
    ax1 = fig.add_subplot(gs[0, :])
    
    x1 = np.linspace(0, 1, len(k1))
    x2 = np.linspace(0, 1, len(k2))
    
    # Log scale for better visibility
    lk1 = np.log10(k1 + 1)
    lk2 = np.log10(k2 + 1)
    
    ax1.fill_between(x1, 0, lk1, alpha=0.3, color=c1_light)
    ax1.fill_between(x2, 0, lk2, alpha=0.3, color=c2_light)
    ax1.plot(x1, lk1, color=c1, linewidth=1.0, label=label1, alpha=0.8)
    ax1.plot(x2, lk2, color=c2, linewidth=1.0, label=label2, alpha=0.8)
    
    # Mark peaks
    peak1 = np.argmax(k1)
    peak2 = np.argmax(k2)
    ax1.axvline(x=peak1/len(k1), color=c1, linestyle='--', alpha=0.5, linewidth=0.8)
    ax1.axvline(x=peak2/len(k2), color=c2, linestyle='--', alpha=0.5, linewidth=0.8)
    
    ax1.set_xlabel('Relative position (0 = N-terminus, 1 = C-terminus)', fontsize=11)
    ax1.set_ylabel('log₁₀(κ² + 1)', fontsize=11)
    ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10, loc='upper right')
    ax1.set_xlim(0, 1)
    
    # Stats box
    stats_text = (
        f"Pearson r = {corr['pearson']:.3f}   "
        f"Spearman ρ = {corr['spearman']:.3f}   "
        f"Log-Pearson r = {corr['log_pearson']:.3f}\n"
        f"{label1}: Gini={stats1['gini']:.4f}  peak={stats1['concentration']:.0f}×  n={stats1['n_residues']}\n"
        f"{label2}: Gini={stats2['gini']:.4f}  peak={stats2['concentration']:.0f}×  n={stats2['n_residues']}"
    )
    ax1.text(0.02, 0.97, stats_text, transform=ax1.transAxes, fontsize=9,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # ─── Secondary structure tracks ───
    ax_ss1 = fig.add_subplot(gs[1, :])
    
    ss_colors = {'H': '#FF6B6B', 'E': '#4ECDC4', 'L': '#F7F7F7'}
    
    for i, s in enumerate(ss1):
        ax_ss1.axvspan(i/len(ss1), (i+1)/len(ss1), color=ss_colors[s], alpha=0.7)
    ax_ss1.set_xlim(0, 1)
    ax_ss1.set_yticks([])
    ax_ss1.set_ylabel(label1[:12], fontsize=9, rotation=0, labelpad=60, va='center')
    ax_ss1.set_title('Secondary Structure', fontsize=10)
    
    ax_ss2 = fig.add_subplot(gs[2, :])
    for i, s in enumerate(ss2):
        ax_ss2.axvspan(i/len(ss2), (i+1)/len(ss2), color=ss_colors[s], alpha=0.7)
    ax_ss2.set_xlim(0, 1)
    ax_ss2.set_yticks([])
    ax_ss2.set_ylabel(label2[:12], fontsize=9, rotation=0, labelpad=60, va='center')
    ax_ss2.set_xlabel('Relative position', fontsize=10)
    
    # SS legend
    legend_elements = [
        Patch(facecolor='#FF6B6B', alpha=0.7, label='α-helix'),
        Patch(facecolor='#4ECDC4', alpha=0.7, label='β-sheet'),
        Patch(facecolor='#F7F7F7', edgecolor='gray', alpha=0.7, label='Loop/coil'),
    ]
    ax_ss2.legend(handles=legend_elements, loc='lower right', fontsize=8, ncol=3)
    
    # ─── pLDDT confidence ───
    ax_plddt = fig.add_subplot(gs[3, :])
    if len(plddt1) > 0:
        px1 = np.linspace(0, 1, len(plddt1))
        ax_plddt.fill_between(px1, 0, plddt1, alpha=0.4, color=c1_light)
    if len(plddt2) > 0:
        px2 = np.linspace(0, 1, len(plddt2))
        ax_plddt.fill_between(px2, 0, plddt2, alpha=0.4, color=c2_light)
    ax_plddt.axhline(y=70, color='orange', linestyle=':', linewidth=0.8, label='pLDDT=70 threshold')
    ax_plddt.set_xlim(0, 1)
    ax_plddt.set_ylim(0, 100)
    ax_plddt.set_ylabel('pLDDT', fontsize=9)
    ax_plddt.set_xlabel('Relative position', fontsize=10)
    ax_plddt.legend(fontsize=8)
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

PAIRS = [
    {
        "title": "PAIR 1: Gaucher Disease — Human GBA1 ↔ Soybean",
        "human_uid": "P04062",
        "human_label": "Human GBA1 (Gaucher)",
        "match_uid": "Q39817",
        "match_label": "Soybean Q39817",
        "disease": "Gaucher disease",
        "clinical": "ERT validated (imiglucerase). Soybean protein is a mechanical twin.",
    },
    {
        "title": "PAIR 2: Alzheimer's Disease — Human TM2D3 ↔ E. coli",
        "human_uid": "Q9BRN9",
        "human_label": "Human TM2D3 (Alzheimer's)",
        "match_uid": "P39310",
        "match_label": "E. coli P39310",
        "disease": "Alzheimer disease",
        "clinical": "TM2D3 involved in Aβ processing. E. coli mimic in every human gut.",
    },
    {
        "title": "PAIR 3: Mitochondrial Disease — Human COX11 ↔ Pseudomonas",
        "human_uid": "Q9Y6N1",
        "human_label": "Human COX11 (Mito disease)",
        "match_uid": "O86062",
        "match_label": "Pseudomonas O86062",
        "disease": "Mitochondrial complex IV deficiency",
        "clinical": "COX11 essential for cytochrome c oxidase. Perfect cross-kingdom match.",
    },
]


def find_cif(cache_dir, uid):
    """Find CIF file for a UniProt ID in the alphafold cache."""
    pattern = os.path.join(cache_dir, "**", f"AF-{uid}-F1-model_v*.cif")
    matches = glob.glob(pattern, recursive=True)
    if matches:
        return matches[0]
    
    # Try without subdirectories
    pattern = os.path.join(cache_dir, f"AF-{uid}-F1-model_v*.cif")
    matches = glob.glob(pattern)
    if matches:
        return matches[0]
    
    # Try just searching
    for root, dirs, files in os.walk(cache_dir):
        for f in files:
            if uid in f and f.endswith('.cif'):
                return os.path.join(root, f)
    
    return None


def main():
    parser = argparse.ArgumentParser(description="Generate curvature overlay plots for validation pairs")
    parser.add_argument("--cache", required=True, help="Path to alphafold_cache directory")
    parser.add_argument("--output", default="overlay_plots", help="Output directory for plots")
    
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    print("=" * 72)
    print("  CURVATURE PROFILE OVERLAY — THREE VALIDATION PAIRS")
    print("=" * 72)
    
    results = []
    
    for pair in PAIRS:
        print(f"\n  {'═' * 60}")
        print(f"  {pair['title']}")
        print(f"  {pair['clinical']}")
        print(f"  {'═' * 60}")
        
        # Find CIF files
        h_cif = find_cif(args.cache, pair["human_uid"])
        m_cif = find_cif(args.cache, pair["match_uid"])
        
        if not h_cif:
            print(f"  ❌ Cannot find CIF for {pair['human_uid']}")
            continue
        if not m_cif:
            print(f"  ❌ Cannot find CIF for {pair['match_uid']}")
            continue
        
        print(f"  Human:  {h_cif}")
        print(f"  Match:  {m_cif}")
        
        # Parse structures
        h_coords, h_plddt = parse_cif_ca(h_cif)
        m_coords, m_plddt = parse_cif_ca(m_cif)
        
        print(f"  Human Cα: {len(h_coords)} residues")
        print(f"  Match Cα: {len(m_coords)} residues")
        
        if len(h_coords) < 10 or len(m_coords) < 10:
            print(f"  ❌ Too few residues, skipping")
            continue
        
        # Compute curvature
        h_kappa2 = compute_kappa2(h_coords)
        m_kappa2 = compute_kappa2(m_coords)
        
        # Stats
        h_stats = profile_stats(h_kappa2, pair["human_label"])
        m_stats = profile_stats(m_kappa2, pair["match_label"])
        
        print(f"\n  {pair['human_label']}:")
        print(f"    Gini = {h_stats['gini']:.4f}, peak = {h_stats['concentration']:.0f}× at residue {h_stats['peak_residue']}")
        print(f"  {pair['match_label']}:")
        print(f"    Gini = {m_stats['gini']:.4f}, peak = {m_stats['concentration']:.0f}× at residue {m_stats['peak_residue']}")
        
        # Correlation
        try:
            corr = profile_correlation(h_kappa2, m_kappa2)
            print(f"\n  Correlation:")
            print(f"    Pearson r = {corr['pearson']:.4f}")
            print(f"    Log-Pearson r = {corr['log_pearson']:.4f}")
            print(f"    Spearman ρ = {corr['spearman']:.4f} (p = {corr['spearman_p']:.2e})")
        except ImportError:
            print(f"\n  [scipy not available — skipping Spearman, using Pearson only]")
            corr = {"pearson": 0, "log_pearson": 0, "spearman": 0, "spearman_p": 1}
            # Manual Pearson
            n_pts = 1000
            x_c = np.linspace(0, 1, n_pts)
            k1i = np.interp(x_c, np.linspace(0,1,len(h_kappa2)), h_kappa2)
            k2i = np.interp(x_c, np.linspace(0,1,len(m_kappa2)), m_kappa2)
            corr["pearson"] = np.corrcoef(k1i, k2i)[0,1]
            corr["log_pearson"] = np.corrcoef(np.log1p(k1i), np.log1p(k2i))[0,1]
            print(f"    Pearson r = {corr['pearson']:.4f}")
            print(f"    Log-Pearson r = {corr['log_pearson']:.4f}")
        
        # SS assignment
        h_ss = assign_ss(h_coords)
        m_ss = assign_ss(m_coords)
        
        # Plot
        safe_name = pair["human_uid"] + "_vs_" + pair["match_uid"]
        outpath = os.path.join(args.output, f"overlay_{safe_name}.png")
        
        try:
            plot_overlay(h_kappa2, m_kappa2, h_plddt, m_plddt,
                        h_ss, m_ss, h_stats, m_stats, corr,
                        pair["title"], pair["human_label"], pair["match_label"],
                        outpath)
        except Exception as e:
            print(f"  ❌ Plot failed: {e}")
            # Save raw data instead
            np.savetxt(os.path.join(args.output, f"kappa2_{pair['human_uid']}.csv"),
                      h_kappa2, delimiter=',', header='kappa2')
            np.savetxt(os.path.join(args.output, f"kappa2_{pair['match_uid']}.csv"),
                      m_kappa2, delimiter=',', header='kappa2')
            print(f"  Saved raw κ² data as CSV instead")
        
        results.append({
            "pair": pair["title"],
            "h_gini": h_stats["gini"],
            "m_gini": m_stats["gini"],
            "h_conc": h_stats["concentration"],
            "m_conc": m_stats["concentration"],
            "pearson": corr["pearson"],
            "log_pearson": corr["log_pearson"],
        })
    
    # Summary
    print(f"\n{'=' * 72}")
    print(f"  SUMMARY OF THREE VALIDATION PAIRS")
    print(f"{'=' * 72}")
    print(f"  {'Pair':50s} {'H_Gini':>7s} {'M_Gini':>7s} {'Pearson':>8s} {'LogR':>8s}")
    print(f"  {'─' * 82}")
    for r in results:
        print(f"  {r['pair']:50s} {r['h_gini']:>7.4f} {r['m_gini']:>7.4f} {r['pearson']:>8.4f} {r['log_pearson']:>8.4f}")
    
    print(f"\n  Plots saved to {args.output}/")
    print(f"  These figures show per-residue backbone curvature alignment")
    print(f"  between human disease proteins and their cross-kingdom mechanical twins.")


if __name__ == "__main__":
    main()

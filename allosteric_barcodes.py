#!/usr/bin/env python3
"""
ALLOSTERIC BARCODES — Curvature Modulation as Allosteric Coordinate
====================================================================

Compares geometric barcodes between two conformational states of the
same protein to detect allosteric signal propagation through backbone
curvature changes.

Test system: Human hemoglobin T-state (deoxy) vs R-state (oxy)
  - 2DN2: deoxy T-state, 1.25 Å
  - 2DN1: oxy R-state, 1.25 Å

The hypothesis: allosteric transitions redistribute backbone curvature
at mechanically coupled sites. If true, κ(s) differences between T and
R states should localize to known allosteric hot spots.

Usage:
    python allosteric_barcodes.py                # Full analysis
    python allosteric_barcodes.py download       # Download structures only
    python allosteric_barcodes.py analyze        # Run comparison
    python allosteric_barcodes.py figure         # Generate figures only

Author: Kase Knochenhauer / True North Construction LLC
"""

import json
import os
import sys
import math
import urllib.request
import urllib.parse
from pathlib import Path
from collections import defaultdict

import numpy as np
from scipy.ndimage import gaussian_filter1d

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "data" / "allosteric"
RESULTS_DIR = SCRIPT_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"

for d in [DATA_DIR, RESULTS_DIR, FIGURES_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Import the core pipeline
sys.path.insert(0, str(SCRIPT_DIR))
from cornu_spirals_on_T2 import (
    extract_dihedrals_manual, extract_dihedrals_mmcif,
    dssp_assign_pure_numpy, assign_ss_from_dihedrals,
    torus_curvature, classify_spiral,
)


# ═══════════════════════════════════════════════════════════════════════
#  CONFIGURATION — Allosteric test systems
# ═══════════════════════════════════════════════════════════════════════

SYSTEMS = {
    'hemoglobin': {
        'name': 'Human Hemoglobin',
        'state1': {'pdb': '2DN2', 'label': 'T-state (deoxy)', 'color': '#2166ac'},
        'state2': {'pdb': '2DN1', 'label': 'R-state (oxy)', 'color': '#b2182b'},
        'chains': {
            'A': {'name': 'α1 subunit', 'length': 141},
            'B': {'name': 'β1 subunit', 'length': 146},
        },
        # Known allosteric sites (1-indexed residue numbers)
        'known_sites': {
            'A': {
                'FG_corner': list(range(91, 100)),     # α FG corner (91-99)
                'F_helix_prox': list(range(82, 92)),   # α F helix proximal
                'C_helix': list(range(40, 50)),        # α C helix
                'heme_contact': [58, 62, 63, 87, 92],  # α heme pocket
                'alpha1_beta2': [38, 40, 41, 44, 92, 94, 95, 96, 97, 140, 141],  # interface
            },
            'B': {
                'FG_corner': list(range(92, 101)),     # β FG corner
                'F_helix_prox': list(range(83, 93)),   # β F helix proximal  
                'C_helix': list(range(41, 51)),        # β C helix
                'heme_contact': [63, 67, 68, 92, 95],  # β heme pocket
                'alpha1_beta2': [37, 40, 41, 97, 99, 100, 101, 102, 145, 146],  # interface
                'DPG_site': [1, 2, 82, 143, 144, 146],  # 2,3-DPG binding
                'salt_bridges': [94, 146],             # His146-Asp94 T-state salt bridge
            },
        },
    },
}


# ═══════════════════════════════════════════════════════════════════════
#  STEP 1: DOWNLOAD
# ═══════════════════════════════════════════════════════════════════════

def download_pdb(pdb_id, output_dir=None):
    """Download a PDB file from RCSB."""
    if output_dir is None:
        output_dir = DATA_DIR
    outpath = output_dir / f"{pdb_id.lower()}.pdb"
    if outpath.exists() and outpath.stat().st_size > 5000:
        return outpath

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    headers = {'User-Agent': 'Mozilla/5.0 PhiPsiAtlas/1.0'}
    try:
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = resp.read()
            with open(outpath, 'wb') as f:
                f.write(data)
            return outpath
    except Exception as e:
        print(f"  Download failed for {pdb_id}: {e}")
        return None


def cmd_download():
    """Download all structure pairs."""
    print("\n  DOWNLOADING ALLOSTERIC STRUCTURE PAIRS\n")
    for sys_name, sys_info in SYSTEMS.items():
        for state in ['state1', 'state2']:
            pdb_id = sys_info[state]['pdb']
            path = download_pdb(pdb_id)
            if path:
                print(f"  ✓ {pdb_id} ({sys_info[state]['label']}): {path}")
            else:
                print(f"  ✗ {pdb_id}: FAILED")


# ═══════════════════════════════════════════════════════════════════════
#  STEP 2: EXTRACT PER-RESIDUE CURVATURE
# ═══════════════════════════════════════════════════════════════════════

def extract_chain_dihedrals(pdb_path, chain_id):
    """Extract dihedrals for a specific chain from a multi-chain PDB."""
    pdb_path = str(pdb_path)
    
    # Parse PDB manually to extract specific chain
    atoms = defaultdict(dict)  # {(chain, resnum): {'N': coord, 'CA': coord, ...}}
    
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                chain = line[21]
                if chain != chain_id:
                    continue
                
                atom_name = line[12:16].strip()
                if atom_name not in ('N', 'CA', 'C', 'O'):
                    continue
                
                # Skip alternate conformations (take A or first)
                altloc = line[16]
                if altloc not in (' ', '', 'A'):
                    continue
                
                resname = line[17:20].strip()
                # Skip non-standard residues
                if resname in ('HOH', 'HEM', 'SO4', 'PO4', 'GOL'):
                    continue
                
                resnum = int(line[22:26].strip())
                icode = line[26].strip()
                key = (resnum, icode)
                
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                atoms[key][atom_name] = np.array([x, y, z])
    
    if not atoms:
        return None, None, None, None, None
    
    # Sort by residue number
    sorted_keys = sorted(atoms.keys())
    n = len(sorted_keys)
    
    # Build arrays
    phi = np.full(n, np.nan)
    psi = np.full(n, np.nan)
    res_ids = []
    backbone_coords = np.zeros((n, 4, 3))
    has_O = np.zeros(n, dtype=bool)
    
    for i, key in enumerate(sorted_keys):
        res_atoms = atoms[key]
        res_ids.append(key[0])
        
        if 'N' in res_atoms:
            backbone_coords[i, 0] = res_atoms['N']
        if 'CA' in res_atoms:
            backbone_coords[i, 1] = res_atoms['CA']
        if 'C' in res_atoms:
            backbone_coords[i, 2] = res_atoms['C']
        if 'O' in res_atoms:
            backbone_coords[i, 3] = res_atoms['O']
            has_O[i] = True
    
    # Compute dihedrals
    for i in range(n):
        a_i = atoms[sorted_keys[i]]
        if i > 0:
            a_prev = atoms[sorted_keys[i-1]]
            if 'C' in a_prev and 'N' in a_i and 'CA' in a_i and 'C' in a_i:
                phi[i] = _dihedral(a_prev['C'], a_i['N'], a_i['CA'], a_i['C'])
        if i < n - 1:
            a_next = atoms[sorted_keys[i+1]]
            if 'N' in a_i and 'CA' in a_i and 'C' in a_i and 'N' in a_next:
                psi[i] = _dihedral(a_i['N'], a_i['CA'], a_i['C'], a_next['N'])
    
    return phi, psi, np.array(res_ids), backbone_coords, has_O


def _dihedral(p1, p2, p3, p4):
    """Compute dihedral angle between 4 points (radians)."""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    if n1_norm < 1e-10 or n2_norm < 1e-10:
        return 0.0
    n1 = n1 / n1_norm
    n2 = n2 / n2_norm
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return np.arctan2(y, x)


def per_residue_curvature(phi, psi):
    """Compute per-residue curvature κ(i) on the flat torus.
    
    Returns array of length n-2 (curvature defined at interior points).
    """
    kappa_result = torus_curvature(phi, psi)
    if isinstance(kappa_result, tuple):
        kappa = kappa_result[0]
    else:
        kappa = kappa_result
    return np.asarray(kappa, dtype=float)


def extract_curvature_profile(pdb_path, chain_id):
    """Extract full curvature profile for a chain.
    
    Returns: (residue_numbers, kappa, ss, phi, psi) or None on failure.
    """
    phi, psi, res_ids, bb_coords, has_O = extract_chain_dihedrals(pdb_path, chain_id)
    
    if phi is None or len(phi) < 10:
        return None
    
    # Filter valid residues
    valid = ~(np.isnan(phi) | np.isnan(psi))
    phi_v = phi[valid]
    psi_v = psi[valid]
    res_v = res_ids[valid]
    bb_v = bb_coords[valid]
    ho_v = has_O[valid]
    
    # SS assignment via DSSP
    ss = None
    if np.sum(ho_v) > 0.8 * len(ho_v):
        try:
            ss = dssp_assign_pure_numpy(bb_v)
        except:
            pass
    if ss is None:
        ss = assign_ss_from_dihedrals(phi_v, psi_v)
    
    # Per-residue curvature
    kappa = per_residue_curvature(phi_v, psi_v)
    
    if len(kappa) < 5:
        return None
    
    # Align residue numbers with kappa (kappa is n-2 for n residues)
    # kappa[i] corresponds to residue i+1 (0-indexed)
    res_kappa = res_v[1:-1] if len(kappa) == len(res_v) - 2 else res_v[:len(kappa)]
    ss_kappa = ss[1:-1] if len(kappa) == len(ss) - 2 else ss[:len(kappa)]
    
    return {
        'res_ids': res_kappa,
        'kappa': kappa,
        'ss': ss_kappa,
        'phi': phi_v,
        'psi': psi_v,
    }


# ═══════════════════════════════════════════════════════════════════════
#  STEP 3: COMPARE STATES
# ═══════════════════════════════════════════════════════════════════════

def compare_states(profile1, profile2, smooth_sigma=1.5):
    """Compare curvature profiles between two conformational states.
    
    Returns dict with per-residue delta-kappa and statistical measures.
    """
    res1 = profile1['res_ids']
    res2 = profile2['res_ids']
    k1 = profile1['kappa']
    k2 = profile2['kappa']
    
    # Find common residues
    set1 = set(res1)
    set2 = set(res2)
    common = sorted(set1 & set2)
    
    if len(common) < 10:
        return None
    
    # Build aligned arrays
    idx1 = {r: i for i, r in enumerate(res1)}
    idx2 = {r: i for i, r in enumerate(res2)}
    
    aligned_res = []
    aligned_k1 = []
    aligned_k2 = []
    aligned_ss1 = []
    aligned_ss2 = []
    
    for r in common:
        if r in idx1 and r in idx2:
            aligned_res.append(r)
            aligned_k1.append(k1[idx1[r]])
            aligned_k2.append(k2[idx2[r]])
            aligned_ss1.append(profile1['ss'][idx1[r]])
            aligned_ss2.append(profile2['ss'][idx2[r]])
    
    aligned_res = np.array(aligned_res)
    aligned_k1 = np.array(aligned_k1)
    aligned_k2 = np.array(aligned_k2)
    
    # Delta curvature
    delta_k = aligned_k2 - aligned_k1
    abs_delta = np.abs(delta_k)
    
    # Smooth for visualization
    if smooth_sigma > 0:
        delta_smooth = gaussian_filter1d(delta_k, sigma=smooth_sigma)
        k1_smooth = gaussian_filter1d(aligned_k1, sigma=smooth_sigma)
        k2_smooth = gaussian_filter1d(aligned_k2, sigma=smooth_sigma)
    else:
        delta_smooth = delta_k
        k1_smooth = aligned_k1
        k2_smooth = aligned_k2
    
    # Identify significant changes (> 1.5 σ from mean)
    mean_delta = np.mean(abs_delta)
    std_delta = np.std(abs_delta)
    threshold = mean_delta + 1.5 * std_delta
    
    significant = abs_delta > threshold
    sig_residues = aligned_res[significant]
    
    # SS transition sites
    ss_changes = np.array([1 if s1 != s2 else 0 
                           for s1, s2 in zip(aligned_ss1, aligned_ss2)])
    
    return {
        'residues': aligned_res,
        'kappa1': aligned_k1,
        'kappa2': aligned_k2,
        'delta_kappa': delta_k,
        'delta_smooth': delta_smooth,
        'kappa1_smooth': k1_smooth,
        'kappa2_smooth': k2_smooth,
        'abs_delta': abs_delta,
        'threshold': threshold,
        'significant_residues': sig_residues,
        'ss1': aligned_ss1,
        'ss2': aligned_ss2,
        'ss_changes': ss_changes,
        'mean_abs_delta': mean_delta,
        'std_abs_delta': std_delta,
        'n_significant': int(np.sum(significant)),
        'n_residues': len(aligned_res),
    }


def score_allosteric_overlap(comparison, known_sites):
    """Score overlap between significant curvature changes and known sites."""
    sig_set = set(comparison['significant_residues'])
    
    results = {}
    total_known = set()
    
    for site_name, site_residues in known_sites.items():
        site_set = set(site_residues)
        total_known |= site_set
        overlap = sig_set & site_set
        n_site = len(site_set)
        n_overlap = len(overlap)
        
        # Enrichment vs random expectation
        n_total = comparison['n_residues']
        n_sig = comparison['n_significant']
        expected = n_site * n_sig / n_total if n_total > 0 else 0
        enrichment = n_overlap / expected if expected > 0 else float('inf')
        
        results[site_name] = {
            'site_residues': sorted(site_set),
            'overlap_residues': sorted(overlap),
            'n_site': n_site,
            'n_overlap': n_overlap,
            'expected': expected,
            'enrichment': enrichment,
        }
    
    # Global enrichment
    global_overlap = sig_set & total_known
    n_total = comparison['n_residues']
    n_sig = comparison['n_significant']
    n_known = len(total_known)
    expected_global = n_known * n_sig / n_total if n_total > 0 else 0
    
    results['_global'] = {
        'n_known': n_known,
        'n_overlap': len(global_overlap),
        'expected': expected_global,
        'enrichment': len(global_overlap) / expected_global if expected_global > 0 else 0,
    }
    
    return results


# ═══════════════════════════════════════════════════════════════════════
#  STEP 4: FIGURES
# ═══════════════════════════════════════════════════════════════════════

def make_allosteric_figure(sys_name, sys_info, comparisons, overlaps):
    """Generate multi-panel figure showing allosteric barcode comparison."""
    chains = list(comparisons.keys())
    n_chains = len(chains)
    
    fig = plt.figure(figsize=(14, 5 * n_chains + 2))
    gs = GridSpec(n_chains * 2, 1, hspace=0.35, height_ratios=[3, 1.5] * n_chains)
    
    state1_label = sys_info['state1']['label']
    state2_label = sys_info['state2']['label']
    c1 = sys_info['state1']['color']
    c2 = sys_info['state2']['color']
    
    fig.suptitle(f"Allosteric Curvature Modulation: {sys_info['name']}\n"
                 f"{state1_label} ({sys_info['state1']['pdb']}) vs "
                 f"{state2_label} ({sys_info['state2']['pdb']})",
                 fontsize=13, fontweight='bold', y=0.98)
    
    for ci, chain_id in enumerate(chains):
        comp = comparisons[chain_id]
        chain_name = sys_info['chains'][chain_id]['name']
        
        # Panel A: Overlaid curvature profiles
        ax1 = fig.add_subplot(gs[ci * 2])
        
        res = comp['residues']
        ax1.fill_between(res, comp['kappa1_smooth'], alpha=0.3, color=c1, label=state1_label)
        ax1.fill_between(res, comp['kappa2_smooth'], alpha=0.3, color=c2, label=state2_label)
        ax1.plot(res, comp['kappa1_smooth'], color=c1, linewidth=1.2)
        ax1.plot(res, comp['kappa2_smooth'], color=c2, linewidth=1.2)
        
        # Mark known allosteric sites
        if chain_id in sys_info.get('known_sites', {}):
            known = sys_info['known_sites'][chain_id]
            for site_name, site_res in known.items():
                for r in site_res:
                    if r in set(res):
                        ax1.axvline(r, color='gold', alpha=0.15, linewidth=2)
        
        ax1.set_ylabel('κ (rad/residue)')
        ax1.set_title(f'Chain {chain_id}: {chain_name}', fontsize=11, loc='left')
        ax1.legend(loc='upper right', fontsize=8)
        ax1.set_xlim(res[0], res[-1])
        
        # Add SS track at bottom of panel
        ss_colors = {'H': '#ff6b6b', 'E': '#4ecdc4', 'C': '#95a5a6'}
        for i, (r, s) in enumerate(zip(res, comp['ss1'])):
            ax1.axvspan(r - 0.5, r + 0.5, ymin=0, ymax=0.03,
                       color=ss_colors.get(s, '#95a5a6'), alpha=0.8)
        
        # Panel B: Delta-curvature (signed difference)
        ax2 = fig.add_subplot(gs[ci * 2 + 1])
        
        # Color by direction
        pos = comp['delta_smooth'] > 0
        neg = ~pos
        ax2.bar(res[pos], comp['delta_smooth'][pos], width=1, color=c2, alpha=0.6, label='R > T')
        ax2.bar(res[neg], comp['delta_smooth'][neg], width=1, color=c1, alpha=0.6, label='T > R')
        
        # Threshold lines
        ax2.axhline(comp['threshold'], color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax2.axhline(-comp['threshold'], color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax2.axhline(0, color='black', linewidth=0.5)
        
        # Mark significant residues
        sig = comp['significant_residues']
        for r in sig:
            idx = np.where(comp['residues'] == r)[0]
            if len(idx) > 0:
                ax2.plot(r, comp['delta_smooth'][idx[0]], 'v', color='red',
                        markersize=4, zorder=5)
        
        # Mark known sites
        if chain_id in overlaps:
            ol = overlaps[chain_id]
            for site_name, site_data in ol.items():
                if site_name.startswith('_'):
                    continue
                for r in site_data.get('overlap_residues', []):
                    idx = np.where(comp['residues'] == r)[0]
                    if len(idx) > 0:
                        ax2.plot(r, comp['delta_smooth'][idx[0]], '*',
                                color='gold', markersize=10, zorder=6,
                                markeredgecolor='black', markeredgewidth=0.5)
        
        ax2.set_ylabel('Δκ (R − T)')
        ax2.set_xlabel('Residue number')
        ax2.legend(loc='upper right', fontsize=7)
        ax2.set_xlim(res[0], res[-1])
    
    # Save
    outpath = FIGURES_DIR / f"fig_allosteric_{sys_name}.png"
    fig.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Figure saved: {outpath}")
    return outpath


# ═══════════════════════════════════════════════════════════════════════
#  MAIN ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

def cmd_analyze():
    """Run full allosteric barcode comparison."""
    print("\n  ALLOSTERIC BARCODE ANALYSIS")
    print("  " + "═" * 60)
    
    all_results = {}
    
    for sys_name, sys_info in SYSTEMS.items():
        print(f"\n  System: {sys_info['name']}")
        print(f"  {sys_info['state1']['label']} ({sys_info['state1']['pdb']}) vs "
              f"{sys_info['state2']['label']} ({sys_info['state2']['pdb']})")
        
        pdb1 = DATA_DIR / f"{sys_info['state1']['pdb'].lower()}.pdb"
        pdb2 = DATA_DIR / f"{sys_info['state2']['pdb'].lower()}.pdb"
        
        if not pdb1.exists() or not pdb2.exists():
            print(f"  Missing structure files. Run 'download' first.")
            continue
        
        comparisons = {}
        overlaps = {}
        
        for chain_id, chain_info in sys_info['chains'].items():
            print(f"\n  Chain {chain_id}: {chain_info['name']}")
            
            # Extract curvature profiles
            prof1 = extract_curvature_profile(pdb1, chain_id)
            prof2 = extract_curvature_profile(pdb2, chain_id)
            
            if prof1 is None or prof2 is None:
                print(f"    FAILED to extract profiles")
                continue
            
            print(f"    State 1: {len(prof1['kappa'])} residues with curvature")
            print(f"    State 2: {len(prof2['kappa'])} residues with curvature")
            
            # Compare
            comp = compare_states(prof1, prof2)
            if comp is None:
                print(f"    FAILED to compare (insufficient common residues)")
                continue
            
            comparisons[chain_id] = comp
            
            print(f"    Common residues: {comp['n_residues']}")
            print(f"    Mean |Δκ|: {comp['mean_abs_delta']:.4f} rad/residue")
            print(f"    Significant changes (>{comp['threshold']:.3f}): "
                  f"{comp['n_significant']}/{comp['n_residues']} "
                  f"({comp['n_significant']/comp['n_residues']*100:.1f}%)")
            print(f"    Significant residues: {list(comp['significant_residues'])}")
            
            # Score against known sites
            known = sys_info.get('known_sites', {}).get(chain_id, {})
            if known:
                ol = score_allosteric_overlap(comp, known)
                overlaps[chain_id] = ol
                
                print(f"\n    Known allosteric site overlap:")
                for site_name, site_data in ol.items():
                    if site_name.startswith('_'):
                        continue
                    e = site_data['enrichment']
                    n_ol = site_data['n_overlap']
                    n_site = site_data['n_site']
                    print(f"      {site_name:20s}: {n_ol}/{n_site} residues "
                          f"(enrichment={e:.1f}×)")
                
                gl = ol['_global']
                print(f"      {'GLOBAL':20s}: {gl['n_overlap']}/{gl['n_known']} "
                      f"(enrichment={gl['enrichment']:.1f}×)")
        
        # Generate figure
        if comparisons:
            fig_path = make_allosteric_figure(sys_name, sys_info, comparisons, overlaps)
        
        # Save results
        result = {
            'system': sys_name,
            'state1': sys_info['state1']['pdb'],
            'state2': sys_info['state2']['pdb'],
            'chains': {},
        }
        for chain_id, comp in comparisons.items():
            result['chains'][chain_id] = {
                'n_residues': comp['n_residues'],
                'n_significant': comp['n_significant'],
                'mean_abs_delta': float(comp['mean_abs_delta']),
                'std_abs_delta': float(comp['std_abs_delta']),
                'threshold': float(comp['threshold']),
                'significant_residues': [int(r) for r in comp['significant_residues']],
                'overlap': {k: {kk: (vv if not isinstance(vv, (np.integer, np.floating)) 
                                     else float(vv))
                                for kk, vv in v.items() 
                                if not isinstance(vv, (list, np.ndarray))}
                           for k, v in overlaps.get(chain_id, {}).items()},
            }
        
        all_results[sys_name] = result
    
    # Save
    with open(RESULTS_DIR / "allosteric_results.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\n  Results saved to {RESULTS_DIR}/allosteric_results.json")


# ═══════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    args = sys.argv[1:] or ['full']
    cmd = args[0]
    
    if cmd == 'download':
        cmd_download()
    elif cmd == 'analyze':
        cmd_analyze()
    elif cmd == 'figure':
        cmd_analyze()  # Figures generated as part of analysis
    elif cmd == 'full':
        cmd_download()
        cmd_analyze()
    else:
        print(f"  Unknown: {cmd}")
        print(f"  Commands: download, analyze, figure, full")


if __name__ == "__main__":
    main()

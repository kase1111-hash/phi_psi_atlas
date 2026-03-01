#!/usr/bin/env python3
"""
Proteome-Scale Curvature Conservation Survey
=============================================
Processes AlphaFold mmCIF structure files, computes per-residue torus
geodesic curvature κ(s), and outputs a CSV with per-protein and
per-region bending energy metrics.

REQUIRES: cornu_spirals_on_T2.py in the same directory (or on PYTHONPATH)

USAGE:
    python proteome_curvature_survey.py /path/to/cif/directory output.csv

    # For multiple directories (e.g., 34 proteomes):
    python proteome_curvature_survey.py /path/to/proteome_db output.csv --recursive

    # Limit to first N files (for testing):
    python proteome_curvature_survey.py /path/to/cifs output.csv --limit 100

    # Resume from where you left off (skips already-processed IDs):
    python proteome_curvature_survey.py /path/to/cifs output.csv --resume

OUTPUT CSV COLUMNS:
    uniprot_id          UniProt accession (from filename or sequence)
    source_file         Filename processed
    length              Number of valid residues (with phi/psi)
    n_segments          Number of barcode segments
    n_gauss_peaks       Number of Gaussian peak segments
    gauss_positions     Comma-separated peak residue positions
    total_kappa_sq      Σκ² (total bending energy)
    kappa_sq_per_res    Σκ²/length (bending energy density)
    mean_kappa          Mean |κ| across all residues
    max_kappa           Maximum |κ|
    std_kappa           Std dev of κ
    kappa_gini          Gini coefficient of κ² (inequality of distribution)
    ss_frac_H           Fraction helix
    ss_frac_E           Fraction sheet
    ss_frac_C           Fraction coil
    mean_plddt          Mean pLDDT confidence score
    min_plddt           Minimum pLDDT
    plddt_below_70      Fraction of residues with pLDDT < 70
    kappa_sq_quartiles  κ² at 25th, 50th, 75th, 90th percentiles
    top_peak_kappa_sq   κ² of highest peak residue
    top_peak_position   Position of highest κ² residue
    processing_time     Seconds to process this protein

LICENSE: CC0 1.0 Universal (Public Domain)
"""

import os
import sys
import csv
import time
import glob
import argparse
import traceback
import numpy as np
from pathlib import Path

# ══════════════════════════════════════════════════════════════
# SELF-CONTAINED CURVATURE COMPUTATION
# (no dependency on cornu_spirals_on_T2.py for the core metric)
# ══════════════════════════════════════════════════════════════

def parse_mmcif_backbone(cif_path):
    """Extract backbone atoms (N, CA, C, O) and pLDDT from mmCIF.
    Returns dict of {resnum: {'N': xyz, 'CA': xyz, 'C': xyz, 'O': xyz, 'plddt': float, 'name': str}}
    """
    columns = []
    in_atom = False
    residues = {}
    
    with open(cif_path, 'r') as f:
        for line in f:
            ls = line.strip()
            if ls.startswith('_atom_site.'):
                in_atom = True
                columns.append(ls.split('.')[1])
            elif in_atom:
                if not ls.startswith('ATOM') and not ls.startswith('HETATM'):
                    if ls.startswith('_') or ls.startswith('#') or ls.startswith('loop'):
                        in_atom = False
                        columns = []
                    continue
                
                parts = ls.split()
                if len(parts) < len(columns):
                    continue
                
                cm = {c: parts[i] for i, c in enumerate(columns) if i < len(parts)}
                
                group = cm.get('group_PDB', 'ATOM')
                if group != 'ATOM':
                    continue
                
                atom_name = cm.get('label_atom_id', cm.get('auth_atom_id', ''))
                if atom_name not in ('N', 'CA', 'C', 'O'):
                    continue
                
                # Use auth_seq_id for residue number (more reliable for AF structures)
                try:
                    resnum = int(cm.get('auth_seq_id', cm.get('label_seq_id', '0')))
                except ValueError:
                    continue
                
                try:
                    x = float(cm.get('Cartn_x', 0))
                    y = float(cm.get('Cartn_y', 0))
                    z = float(cm.get('Cartn_z', 0))
                except ValueError:
                    continue
                
                resname = cm.get('auth_comp_id', cm.get('label_comp_id', ''))
                
                # pLDDT is stored in B-factor column for AlphaFold structures
                try:
                    plddt = float(cm.get('B_iso_or_equiv', 0))
                except ValueError:
                    plddt = 0.0
                
                if resnum not in residues:
                    residues[resnum] = {'name': resname, 'plddt': plddt}
                
                residues[resnum][atom_name] = np.array([x, y, z])
    
    return residues


def compute_dihedrals(residues):
    """Compute phi/psi from backbone coordinates.
    Returns: phi, psi, plddt arrays (length n), all aligned to valid residues.
    """
    sorted_keys = sorted(k for k, v in residues.items() 
                         if all(a in v for a in ('N', 'CA', 'C')))
    n = len(sorted_keys)
    if n < 5:
        return None, None, None, 0
    
    phi = np.full(n, np.nan)
    psi = np.full(n, np.nan)
    plddt = np.zeros(n)
    
    def _dihedral(p0, p1, p2, p3):
        b1 = p1 - p0; b2 = p2 - p1; b3 = p3 - p2
        n1 = np.cross(b1, b2); n2 = np.cross(b2, b3)
        nn1, nn2 = np.linalg.norm(n1), np.linalg.norm(n2)
        if nn1 < 1e-10 or nn2 < 1e-10:
            return np.nan
        n1 /= nn1; n2 /= nn2
        m1 = np.cross(n1, b2 / np.linalg.norm(b2))
        return np.arctan2(np.dot(m1, n2), np.dot(n1, n2))
    
    for i, key in enumerate(sorted_keys):
        at = residues[key]
        plddt[i] = at.get('plddt', 0.0)
        
        if i > 0:
            prev_key = sorted_keys[i-1]
            if 'C' in residues[prev_key]:
                phi[i] = _dihedral(residues[prev_key]['C'], at['N'], at['CA'], at['C'])
        
        if i < n - 1:
            next_key = sorted_keys[i+1]
            if 'N' in residues[next_key]:
                psi[i] = _dihedral(at['N'], at['CA'], at['C'], residues[next_key]['N'])
    
    return phi, psi, plddt, n


def torus_curvature_standalone(phi, psi):
    """Compute geodesic curvature on flat torus (R=r=1).
    Returns kappa array (length n-1 where n = number of valid phi/psi pairs).
    """
    valid = ~(np.isnan(phi) | np.isnan(psi))
    pv = phi[valid]
    sv = psi[valid]
    n = len(pv)
    
    if n < 3:
        return np.array([]), np.array([])
    
    # Points on flat torus
    # Arc lengths between consecutive points
    dphi = np.diff(pv)
    dpsi = np.diff(sv)
    
    # Wrap to [-pi, pi]
    dphi = (dphi + np.pi) % (2 * np.pi) - np.pi
    dpsi = (dpsi + np.pi) % (2 * np.pi) - np.pi
    
    ds = np.sqrt(dphi**2 + dpsi**2)
    ds = np.maximum(ds, 1e-10)
    
    # Tangent vectors
    tx = dphi / ds
    ty = dpsi / ds
    
    # Curvature from tangent angle change
    if len(tx) < 2:
        return np.array([]), np.array([])
    
    # Turning angle between consecutive tangent vectors
    cross = tx[:-1] * ty[1:] - ty[:-1] * tx[1:]
    dot = tx[:-1] * tx[1:] + ty[:-1] * ty[1:]
    theta = np.arctan2(cross, dot)
    
    # Mean arc length at each interior point
    s_mean = 0.5 * (ds[:-1] + ds[1:])
    s_mean = np.maximum(s_mean, 1e-10)
    
    kappa = theta / s_mean
    
    # Cumulative arc length for position reference
    s_cumul = np.concatenate([[0], np.cumsum(ds)])
    s_kappa = 0.5 * (s_cumul[:-2] + s_cumul[1:-1])  # midpoints
    
    return kappa, s_kappa


def assign_ss_from_dihedrals(phi, psi):
    """Simple secondary structure from Ramachandran regions."""
    n = len(phi)
    ss = np.full(n, 'C', dtype='U1')
    for i in range(n):
        p, s = phi[i], psi[i]
        if np.isnan(p) or np.isnan(s):
            continue
        pd, sd = np.degrees(p), np.degrees(s)
        if -160 < pd < -20 and -80 < sd < -10:
            ss[i] = 'H'  # alpha helix
        elif -180 < pd < -20 and 50 < sd < 180:
            ss[i] = 'E'  # beta sheet
        elif -180 < pd < -20 and -180 < sd < -120:
            ss[i] = 'E'  # extended
    return ss


# ══════════════════════════════════════════════════════════════
# SELF-CONTAINED DSSP (pure numpy, ~97% agreement with original)
# From cornu_spirals_on_T2.py by Kase Knochenhauer, CC0 1.0
# ══════════════════════════════════════════════════════════════

_DSSP_Q1Q2 = 0.084
_DSSP_F = 332
_DSSP_CUTOFF = -0.5
_DSSP_MARGIN = 1.0


def _dssp_unfold(a, window, axis):
    idx = (np.arange(window)[:, None] +
           np.arange(a.shape[axis] - window + 1)[None, :])
    unfolded = np.take(a, idx, axis=axis)
    return np.moveaxis(unfolded, axis - 1, -1)


def _dssp_upsample(a, window=3):
    return _dssp_unfold(_dssp_unfold(a, window, -2), window, -2)


def _dssp_hydrogen_positions(coord):
    vec_cn = coord[1:, 0] - coord[:-1, 2]
    vec_cn = vec_cn / np.linalg.norm(vec_cn, axis=-1, keepdims=True)
    vec_can = coord[1:, 0] - coord[1:, 1]
    vec_can = vec_can / np.linalg.norm(vec_can, axis=-1, keepdims=True)
    vec_nh = vec_cn + vec_can
    vec_nh = vec_nh / np.linalg.norm(vec_nh, axis=-1, keepdims=True)
    return coord[1:, 0] + 1.01 * vec_nh


def _dssp_hbond_map(coord):
    n_res = coord.shape[0]
    h_pos = _dssp_hydrogen_positions(coord)
    n_1 = coord[1:, 0]
    c_0 = coord[:-1, 2]
    o_0 = coord[:-1, 3]
    n = n_res - 1
    cmap = np.tile(c_0, (n, 1, 1))
    omap = np.tile(o_0, (n, 1, 1))
    nmap = np.tile(n_1, (1, 1, n)).reshape(n, n, 3)
    hmap = np.tile(h_pos, (1, 1, n)).reshape(n, n, 3)
    d_on = np.linalg.norm(omap - nmap, axis=-1)
    d_ch = np.linalg.norm(cmap - hmap, axis=-1)
    d_oh = np.linalg.norm(omap - hmap, axis=-1)
    d_cn = np.linalg.norm(cmap - nmap, axis=-1)
    e = np.pad(
        _DSSP_Q1Q2 * (1.0/d_on + 1.0/d_ch - 1.0/d_oh - 1.0/d_cn) * _DSSP_F,
        [[1, 0], [0, 1]]
    )
    local_mask = ~np.eye(n_res, dtype=bool)
    local_mask *= ~np.diag(np.ones(n_res - 1, dtype=bool), k=-1)
    local_mask *= ~np.diag(np.ones(n_res - 2, dtype=bool), k=-2)
    hbond_map = np.clip(_DSSP_CUTOFF - _DSSP_MARGIN - e, a_min=-_DSSP_MARGIN, a_max=_DSSP_MARGIN)
    hbond_map = (np.sin(hbond_map / _DSSP_MARGIN * np.pi / 2) + 1.0) / 2
    hbond_map *= local_mask
    return hbond_map


def dssp_assign(coord):
    """Assign H/E/C secondary structure using hydrogen-bond-based DSSP.
    coord: (n_residues, 4, 3) array of (N, CA, C, O) coordinates.
    Returns: (n_residues,) array of 'H', 'E', 'C'.
    """
    n_res = coord.shape[0]
    if n_res < 6:
        return np.full(n_res, 'C', dtype='<U1')
    hbmap = _dssp_hbond_map(coord)
    hbmap = np.swapaxes(hbmap, -1, -2)
    turn3 = np.diagonal(hbmap, offset=3) > 0.0
    turn4 = np.diagonal(hbmap, offset=4) > 0.0
    turn5 = np.diagonal(hbmap, offset=5) > 0.0
    h3 = np.pad(turn3[:-1] * turn3[1:], [[1, 3]])
    h4 = np.pad(turn4[:-1] * turn4[1:], [[1, 4]])
    h5 = np.pad(turn5[:-1] * turn5[1:], [[1, 5]])
    helix4 = h4 + np.roll(h4, 1, 0) + np.roll(h4, 2, 0) + np.roll(h4, 3, 0)
    h3 = h3 * ~np.roll(helix4, -1, 0) * ~helix4
    h5 = h5 * ~np.roll(helix4, -1, 0) * ~helix4
    helix3 = h3 + np.roll(h3, 1, 0) + np.roll(h3, 2, 0)
    helix5 = h5 + np.roll(h5, 1, 0) + np.roll(h5, 2, 0) + np.roll(h5, 3, 0) + np.roll(h5, 4, 0)
    unfoldmap = _dssp_upsample(hbmap, 3) > 0.0
    unfoldmap_rev = np.swapaxes(unfoldmap, 0, 1)
    p_bridge = ((unfoldmap[:, :, 0, 1] * unfoldmap_rev[:, :, 1, 2]) +
                (unfoldmap_rev[:, :, 0, 1] * unfoldmap[:, :, 1, 2]))
    p_bridge = np.pad(p_bridge, [[1, 1], [1, 1]])
    a_bridge = ((unfoldmap[:, :, 1, 1] * unfoldmap_rev[:, :, 1, 1]) +
                (unfoldmap[:, :, 0, 2] * unfoldmap_rev[:, :, 0, 2]))
    a_bridge = np.pad(a_bridge, [[1, 1], [1, 1]])
    ladder = (p_bridge + a_bridge).sum(-1) > 0.0
    helix = (helix3 + helix4 + helix5) > 0.0
    strand = ladder
    ss = np.full(n_res, 'C', dtype='<U1')
    ss[helix] = 'H'
    ss[strand] = 'E'
    return ss


def gini_coefficient(values):
    """Gini coefficient: 0 = perfectly equal, 1 = maximally concentrated."""
    v = np.abs(values)
    if np.sum(v) == 0:
        return 0.0
    v = np.sort(v)
    n = len(v)
    idx = np.arange(1, n + 1)
    return (2 * np.sum(idx * v) / (n * np.sum(v))) - (n + 1) / n


def process_one_structure(cif_path, fast=False, use_dssp=False):
    """Process a single CIF file → dict of metrics."""
    t0 = time.time()
    
    # Parse
    residues = parse_mmcif_backbone(cif_path)
    if len(residues) < 10:
        return None
    
    # Dihedrals
    phi, psi, plddt_arr, n_total = compute_dihedrals(residues)
    if phi is None or n_total < 10:
        return None
    
    # Valid mask
    valid = ~(np.isnan(phi) | np.isnan(psi))
    n_valid = int(np.sum(valid))
    if n_valid < 10:
        return None
    
    pv = phi[valid]
    sv = psi[valid]
    plddt_valid = plddt_arr[valid]
    
    # Curvature
    kappa, s_kappa = torus_curvature_standalone(phi, psi)
    if len(kappa) < 5:
        return None
    
    kappa_sq = kappa ** 2
    
    # Secondary structure — DSSP or dihedral-based
    if use_dssp:
        try:
            sorted_keys = sorted(k for k, v in residues.items()
                                 if all(a in v for a in ('N', 'CA', 'C')))
            coord4 = np.zeros((len(sorted_keys), 4, 3))
            all_have_O = True
            for i, k in enumerate(sorted_keys):
                coord4[i, 0] = residues[k]['N']
                coord4[i, 1] = residues[k]['CA']
                coord4[i, 2] = residues[k]['C']
                if 'O' in residues[k]:
                    coord4[i, 3] = residues[k]['O']
                else:
                    all_have_O = False
            
            if all_have_O and len(sorted_keys) >= 6:
                # DSSP uses O(n²) memory; fall back to dihedrals for huge proteins
                if len(sorted_keys) > 1500:
                    ss = assign_ss_from_dihedrals(pv, sv)
                else:
                    ss_full = dssp_assign(coord4)
                # Map to valid-only residues
                ss = ss_full[valid[:len(ss_full)]] if len(ss_full) >= len(valid) else ss_full
                if len(ss) != n_valid:
                    ss = ss[:n_valid] if len(ss) > n_valid else np.pad(ss, (0, n_valid - len(ss)), constant_values='C')
            else:
                ss = assign_ss_from_dihedrals(pv, sv)
        except Exception:
            ss = assign_ss_from_dihedrals(pv, sv)
    else:
        ss = assign_ss_from_dihedrals(pv, sv)
    n_h = int(np.sum(ss == 'H'))
    n_e = int(np.sum(ss == 'E'))
    n_c = int(np.sum(ss == 'C'))
    
    # ── Core metrics ──
    total_kappa_sq = float(np.sum(kappa_sq))
    kappa_sq_per_res = total_kappa_sq / n_valid
    mean_kappa = float(np.mean(np.abs(kappa)))
    max_kappa = float(np.max(np.abs(kappa)))
    std_kappa = float(np.std(kappa))
    kappa_gini = gini_coefficient(kappa_sq)
    
    # Quantiles of κ²
    q25, q50, q75, q90 = np.percentile(kappa_sq, [25, 50, 75, 90])
    
    # Top peak
    top_idx = int(np.argmax(kappa_sq))
    top_peak_ksq = float(kappa_sq[top_idx])
    top_peak_pos = top_idx  # residue index
    
    # pLDDT stats
    mean_plddt = float(np.mean(plddt_valid)) if len(plddt_valid) > 0 else 0.0
    min_plddt = float(np.min(plddt_valid)) if len(plddt_valid) > 0 else 0.0
    plddt_below_70 = float(np.mean(plddt_valid < 70)) if len(plddt_valid) > 0 else 0.0
    
    # ── Gaussian peak detection (skip if --fast) ──
    n_gauss = 0
    gauss_positions = []
    if not fast:
        try:
            from cornu_spirals_on_T2 import classify_spiral
            
            # Segment by SS boundaries
            sc = np.where(ss[1:] != ss[:-1])[0] + 1
            bounds = np.concatenate([[0], sc, [len(ss)]])
            
            for i in range(len(bounds) - 1):
                st = int(bounds[i])
                en = int(bounds[i+1]) - 1
                sl = en - st + 1
                if sl < 4:
                    continue
                sk = kappa[max(0, st):min(en+1, len(kappa))]
                if len(sk) < 4:
                    continue
                sk_s = s_kappa[max(0, st):min(en+1, len(s_kappa))]
                if len(sk_s) < 4:
                    continue
                sr = sk_s - sk_s[0]
                
                try:
                    _, _, cls, fits, _ = classify_spiral(sk, sr)
                    if cls == 'gauss_peak':
                        n_gauss += 1
                        if 'gauss_peak' in fits:
                            params = fits['gauss_peak'].get('params', {})
                            s0 = params.get('s0', 0)
                            if len(sr) > 0:
                                peak_res = st + int(np.argmin(np.abs(sr - s0)))
                            else:
                                peak_res = st
                            gauss_positions.append(peak_res)
                except Exception:
                    pass
        except ImportError:
            pass  # classify_spiral not available; skip peak detection
    
    elapsed = time.time() - t0
    
    return {
        'length': n_valid,
        'n_segments': 0,  # filled by classify_spiral if available
        'n_gauss_peaks': n_gauss,
        'gauss_positions': ','.join(str(p) for p in gauss_positions),
        'total_kappa_sq': round(total_kappa_sq, 2),
        'kappa_sq_per_res': round(kappa_sq_per_res, 2),
        'mean_kappa': round(mean_kappa, 4),
        'max_kappa': round(max_kappa, 4),
        'std_kappa': round(std_kappa, 4),
        'kappa_gini': round(kappa_gini, 4),
        'ss_frac_H': round(n_h / n_valid, 3),
        'ss_frac_E': round(n_e / n_valid, 3),
        'ss_frac_C': round(n_c / n_valid, 3),
        'mean_plddt': round(mean_plddt, 1),
        'min_plddt': round(min_plddt, 1),
        'plddt_below_70': round(plddt_below_70, 3),
        'kappa_sq_q25': round(q25, 2),
        'kappa_sq_q50': round(q50, 2),
        'kappa_sq_q75': round(q75, 2),
        'kappa_sq_q90': round(q90, 2),
        'top_peak_kappa_sq': round(top_peak_ksq, 2),
        'top_peak_position': top_peak_pos,
        'processing_time': round(elapsed, 3),
    }


def extract_uniprot_id(filepath):
    """Try to extract UniProt ID from AlphaFold-style filename.
    Handles patterns like:
      AF-P68871-F1-model_v4.cif
      structure1_CA2_P00918.pdb.A.cif
      P68871.cif
    """
    basename = os.path.basename(filepath)
    
    # Pattern: AF-{UNIPROT}-F1-model_v{N}.cif
    if basename.startswith('AF-'):
        parts = basename.split('-')
        if len(parts) >= 2:
            return parts[1]
    
    # Pattern: ..._{UNIPROT}.pdb... or {UNIPROT}.cif
    import re
    # Match UniProt accession pattern: [OPQ][0-9][A-Z0-9]{3}[0-9] or [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
    matches = re.findall(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?', basename)
    if matches:
        return matches[-1]  # Take last match (often the UniProt ID after gene name)
    
    # Fallback: use filename stem
    return basename.split('.')[0]


def find_cif_files(input_path, recursive=False):
    """Find all .cif files in the given path."""
    p = Path(input_path)
    if p.is_file():
        return [str(p)]
    
    if recursive:
        return sorted(str(f) for f in p.rglob('*.cif'))
    else:
        return sorted(str(f) for f in p.glob('*.cif'))


def main():
    parser = argparse.ArgumentParser(
        description='Proteome-scale curvature conservation survey',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  # Process a directory of CIF files:
  python proteome_curvature_survey.py ./human_proteome/ human_curvature.csv

  # Process all proteomes recursively:
  python proteome_curvature_survey.py ./proteome_db/ all_proteomes.csv --recursive

  # Quick test on 50 files:
  python proteome_curvature_survey.py ./human_proteome/ test.csv --limit 50

  # Resume interrupted run:
  python proteome_curvature_survey.py ./human_proteome/ human_curvature.csv --resume
        """
    )
    parser.add_argument('input', help='Directory containing .cif files, or a single .cif file')
    parser.add_argument('output', help='Output CSV path')
    parser.add_argument('--recursive', '-r', action='store_true',
                        help='Search subdirectories recursively')
    parser.add_argument('--limit', '-n', type=int, default=0,
                        help='Process only first N files (0 = all)')
    parser.add_argument('--resume', action='store_true',
                        help='Skip files already in output CSV')
    parser.add_argument('--workers', '-w', type=int, default=1,
                        help='Number of parallel workers (default: 1, safe for large files)')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress progress output')
    parser.add_argument('--fast', '-f', action='store_true',
                        help='Skip Gaussian peak classification (much faster, ~20x). '
                             'Core metrics (Σκ², Gini, pLDDT) still computed.')
    parser.add_argument('--dssp', action='store_true',
                        help='Use hydrogen-bond-based DSSP for SS assignment (more accurate, '
                             'adds ~25ms per protein). Default uses dihedral-based SS.')
    args = parser.parse_args()
    
    # Find files
    cif_files = find_cif_files(args.input, args.recursive)
    if not cif_files:
        print(f"ERROR: No .cif files found in {args.input}")
        sys.exit(1)
    
    if args.limit > 0:
        cif_files = cif_files[:args.limit]
    
    total = len(cif_files)
    print(f"Found {total} CIF files to process")
    
    # Resume support
    already_done = set()
    if args.resume and os.path.exists(args.output):
        with open(args.output, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                already_done.add(row.get('source_file', ''))
        print(f"Resuming: {len(already_done)} files already processed")
    
    # CSV setup
    fieldnames = [
        'uniprot_id', 'source_file', 'length', 'n_segments', 'n_gauss_peaks',
        'gauss_positions', 'total_kappa_sq', 'kappa_sq_per_res', 'mean_kappa',
        'max_kappa', 'std_kappa', 'kappa_gini', 'ss_frac_H', 'ss_frac_E',
        'ss_frac_C', 'mean_plddt', 'min_plddt', 'plddt_below_70',
        'kappa_sq_q25', 'kappa_sq_q50', 'kappa_sq_q75', 'kappa_sq_q90',
        'top_peak_kappa_sq', 'top_peak_position', 'processing_time'
    ]
    
    mode = 'a' if args.resume and os.path.exists(args.output) else 'w'
    csvfile = open(args.output, mode, newline='')
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    if mode == 'w':
        writer.writeheader()
    
    # Process
    n_success = 0
    n_fail = 0
    n_skip = 0
    t_start = time.time()
    
    for i, cif_path in enumerate(cif_files):
        basename = os.path.basename(cif_path)
        
        if basename in already_done:
            n_skip += 1
            continue
        
        try:
            result = process_one_structure(cif_path, fast=args.fast, use_dssp=args.dssp)
            
            if result is None:
                n_fail += 1
                if not args.quiet:
                    print(f"  [{i+1}/{total}] SKIP {basename} (too short or parse error)")
                continue
            
            result['uniprot_id'] = extract_uniprot_id(cif_path)
            result['source_file'] = basename
            
            writer.writerow(result)
            csvfile.flush()  # flush after each row for resume safety
            n_success += 1
            
            if not args.quiet and (n_success % 100 == 0 or n_success <= 10):
                elapsed = time.time() - t_start
                rate = n_success / elapsed if elapsed > 0 else 0
                eta = (total - i - 1) / rate if rate > 0 else 0
                print(f"  [{i+1}/{total}] {basename}: L={result['length']}, "
                      f"\u03A3\u03BA\u00B2/res={result['kappa_sq_per_res']:.1f}, "
                      f"pLDDT={result['mean_plddt']:.0f}, "
                      f"gini={result['kappa_gini']:.2f} "
                      f"({rate:.1f}/s, ETA {eta/60:.0f}min)")
        
        except Exception as e:
            n_fail += 1
            if not args.quiet:
                print(f"  [{i+1}/{total}] ERROR {basename}: {e}")
            continue
    
    csvfile.close()
    
    elapsed = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"  COMPLETE")
    print(f"  Processed: {n_success}")
    print(f"  Failed:    {n_fail}")
    print(f"  Skipped:   {n_skip}")
    print(f"  Time:      {elapsed:.0f}s ({elapsed/60:.1f}min)")
    print(f"  Rate:      {n_success/elapsed:.1f} proteins/sec")
    print(f"  Output:    {args.output}")
    print(f"{'='*60}")
    
    # Quick summary stats if we have data
    if n_success > 0:
        print(f"\n  QUICK STATS (from successful runs):")
        try:
            import pandas as pd
            df = pd.read_csv(args.output)
            print(f"  Mean length:        {df['length'].mean():.0f} residues")
            print(f"  Mean \u03A3\u03BA\u00B2/res:       {df['kappa_sq_per_res'].mean():.1f} \u00B1 {df['kappa_sq_per_res'].std():.1f}")
            print(f"  Mean Gini:          {df['kappa_gini'].mean():.3f} \u00B1 {df['kappa_gini'].std():.3f}")
            print(f"  Mean pLDDT:         {df['mean_plddt'].mean():.1f}")
            print(f"  Correlation (length vs \u03A3\u03BA\u00B2/res): {df['length'].corr(df['kappa_sq_per_res']):.3f}")
            
            # Length-bin analysis
            bins = [0, 100, 200, 400, 800, 1600, 99999]
            labels = ['<100', '100-200', '200-400', '400-800', '800-1600', '>1600']
            df['len_bin'] = pd.cut(df['length'], bins=bins, labels=labels)
            print(f"\n  \u03A3\u03BA\u00B2/residue by protein length:")
            for lb in labels:
                sub = df[df['len_bin'] == lb]
                if len(sub) > 0:
                    print(f"    {lb:>8s}: {sub['kappa_sq_per_res'].mean():6.1f} \u00B1 {sub['kappa_sq_per_res'].std():5.1f} (n={len(sub)})")
        except ImportError:
            print("  (install pandas for summary statistics)")


if __name__ == '__main__':
    main()

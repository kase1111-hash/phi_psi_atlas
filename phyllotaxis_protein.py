#!/usr/bin/env python3
"""
PHYLLOTAXIS ↔ PROTEIN BACKBONE GEOMETRY ANALYZER
=================================================
Kase / True North Construction LLC — 2026-02-26

Analyzes protein backbone trajectories on the Ramachandran torus T²
for phyllotactic structure: winding numbers, golden-angle proximity,
parastichy-like basin-visit patterns, and quasi-periodic angular advances.

Uses cached AlphaFold CIF/PDB files. No network required if structures
are already downloaded.

Usage:
    python phyllotaxis_protein.py --cache ./alphafold_cache
    python phyllotaxis_protein.py --cache ./alphafold_cache --uniprot P0CG47
    python phyllotaxis_protein.py --cache ./alphafold_cache --organism human --max 50
    python phyllotaxis_protein.py --demo  # synthetic demo with no files needed

Output:
    phyllotaxis_report.png   — Multi-panel analysis figure
    phyllotaxis_results.json — Machine-readable results
"""

import argparse
import json
import gzip
import math
import re
import sys
from pathlib import Path
from collections import Counter, defaultdict
from dataclasses import dataclass, field, asdict
from typing import Optional

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import Normalize
from scipy import stats, signal

# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

GOLDEN_ANGLE_DEG = 137.507764  # 360 * (1 - 1/φ)
GOLDEN_RATIO = (1 + math.sqrt(5)) / 2
RAD = np.pi / 180

# Ramachandran basin centers (degrees)
BASINS = {
    'α':   (-63, -43),     # α-helix
    'β':   (-120, 130),    # β-strand
    'PPII': (-75, 145),    # polyproline II
    'αL':  (57, 47),       # left-handed helix
}

# α-helix angular advance per residue
ALPHA_STEP = (-63 * RAD, -43 * RAD)  # ~ (-1.1, -0.75) rad
# → per-turn advance: 3.6 residues × step = (-226.8°, -154.8°)
# Effective angular advance on T²: ~100° per residue (|Δ| on torus)

# Noble numbers (continued fraction [1,1,1,...] = φ, [2,1,1,...], etc.)
NOBLE_NUMBERS = {
    'φ (golden)': GOLDEN_RATIO,
    'silver': 1 + math.sqrt(2),
    'bronze': (3 + math.sqrt(13)) / 2,
    'φ²': GOLDEN_RATIO**2,
}

# ═══════════════════════════════════════════════════════════════════════
# CIF / PDB PARSER (no BioPython dependency)
# ═══════════════════════════════════════════════════════════════════════

def open_structure_file(path):
    """Open CIF/PDB file, handling .gz transparently."""
    p = Path(path)
    if p.suffix == '.gz':
        return gzip.open(p, 'rt')
    return open(p, 'r')


def extract_backbone_atoms(filepath):
    """
    Extract backbone atom coordinates from CIF or PDB file.
    Returns dict: residue_index → {'N': xyz, 'CA': xyz, 'C': xyz}
    """
    p = Path(filepath)
    ext = p.suffixes[0] if p.suffixes else p.suffix  # handle .cif.gz

    if ext == '.cif':
        return _parse_cif_backbone(filepath)
    elif ext == '.pdb':
        return _parse_pdb_backbone(filepath)
    else:
        # Try CIF first, fall back to PDB
        try:
            return _parse_cif_backbone(filepath)
        except Exception:
            return _parse_pdb_backbone(filepath)


def _parse_cif_backbone(filepath):
    """Parse mmCIF format for backbone N, CA, C atoms."""
    atoms = {}
    in_atom_site = False
    col_map = {}
    col_count = 0

    with open_structure_file(filepath) as f:
        for line in f:
            line = line.strip()

            if line.startswith('_atom_site.'):
                in_atom_site = True
                col_name = line.split('.')[1].split()[0]
                col_map[col_name] = col_count
                col_count += 1
                continue

            if in_atom_site and not line.startswith('_') and not line.startswith('#') and line:
                if line.startswith('loop_') or line.startswith('data_'):
                    in_atom_site = False
                    continue

                parts = line.split()
                if len(parts) < col_count:
                    in_atom_site = False
                    continue

                try:
                    group = parts[col_map.get('group_PDB', 0)]
                    if group not in ('ATOM', 'HETATM'):
                        continue

                    atom_name = parts[col_map.get('label_atom_id', 0)].strip('"')
                    if atom_name not in ('N', 'CA', 'C'):
                        continue

                    # Use label_seq_id for residue number
                    seq_id_key = 'label_seq_id'
                    if seq_id_key not in col_map:
                        seq_id_key = 'auth_seq_id'
                    res_seq = int(parts[col_map[seq_id_key]])

                    # Only chain A / first model
                    model_key = 'pdbx_PDB_model_num'
                    if model_key in col_map:
                        model = parts[col_map[model_key]]
                        if model != '1':
                            continue

                    x = float(parts[col_map['Cartn_x']])
                    y = float(parts[col_map['Cartn_y']])
                    z = float(parts[col_map['Cartn_z']])

                    if res_seq not in atoms:
                        atoms[res_seq] = {}
                    atoms[res_seq][atom_name] = np.array([x, y, z])

                except (KeyError, ValueError, IndexError):
                    continue

    return atoms


def _parse_pdb_backbone(filepath):
    """Parse PDB format for backbone N, CA, C atoms."""
    atoms = {}
    with open_structure_file(filepath) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                if line.startswith('ENDMDL'):
                    break  # first model only
                continue

            atom_name = line[12:16].strip()
            if atom_name not in ('N', 'CA', 'C'):
                continue

            try:
                res_seq = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                if res_seq not in atoms:
                    atoms[res_seq] = {}
                atoms[res_seq][atom_name] = np.array([x, y, z])
            except (ValueError, IndexError):
                continue

    return atoms


# ═══════════════════════════════════════════════════════════════════════
# DIHEDRAL COMPUTATION
# ═══════════════════════════════════════════════════════════════════════

def dihedral_angle(p1, p2, p3, p4):
    """Compute dihedral angle (radians) from four atom positions."""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    norm1 = np.linalg.norm(n1)
    norm2 = np.linalg.norm(n2)

    if norm1 < 1e-10 or norm2 < 1e-10:
        return None

    n1 /= norm1
    n2 /= norm2

    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    return np.arctan2(y, x)


def compute_dihedrals(backbone_atoms):
    """
    Compute (φ, ψ) dihedral angles for each residue.
    φ_i = dihedral(C_{i-1}, N_i, CA_i, C_i)
    ψ_i = dihedral(N_i, CA_i, C_i, N_{i+1})
    Returns: arrays of phi, psi in radians (length L-2 for interior residues)
    """
    sorted_res = sorted(backbone_atoms.keys())
    phi_list = []
    psi_list = []
    res_indices = []

    for idx in range(1, len(sorted_res) - 1):
        r_prev = sorted_res[idx - 1]
        r_curr = sorted_res[idx]
        r_next = sorted_res[idx + 1]

        prev = backbone_atoms.get(r_prev, {})
        curr = backbone_atoms.get(r_curr, {})
        nxt = backbone_atoms.get(r_next, {})

        # Need: C_{i-1}, N_i, CA_i, C_i, N_{i+1}
        if not all(k in prev for k in ['C']):
            continue
        if not all(k in curr for k in ['N', 'CA', 'C']):
            continue
        if not all(k in nxt for k in ['N']):
            continue

        phi = dihedral_angle(prev['C'], curr['N'], curr['CA'], curr['C'])
        psi = dihedral_angle(curr['N'], curr['CA'], curr['C'], nxt['N'])

        if phi is not None and psi is not None:
            phi_list.append(phi)
            psi_list.append(psi)
            res_indices.append(r_curr)

    return np.array(phi_list), np.array(psi_list), res_indices


# ═══════════════════════════════════════════════════════════════════════
# PHYLLOTAXIS ANALYSIS CORE
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class PhyllotaxisResult:
    """Complete phyllotaxis analysis for one protein."""
    protein_id: str
    length: int
    # Winding numbers
    winding_phi: float = 0.0        # total phi winding / 2π
    winding_psi: float = 0.0        # total psi winding / 2π
    winding_ratio: float = 0.0      # |w_phi / w_psi|
    nearest_noble: str = ""         # closest noble number name
    noble_distance: float = 0.0     # distance to nearest noble number
    # Angular advance statistics
    mean_advance_phi: float = 0.0   # mean Δφ per residue (degrees)
    mean_advance_psi: float = 0.0   # mean Δψ per residue (degrees)
    advance_std_phi: float = 0.0
    advance_std_psi: float = 0.0
    effective_angle: float = 0.0    # |Δ| on torus per step (degrees)
    golden_proximity: float = 0.0   # how close effective angle is to 137.5°
    # Quasi-periodicity
    dominant_frequency: float = 0.0  # strongest spectral peak (cycles/residue)
    spectral_entropy: float = 0.0    # periodicity measure (lower = more periodic)
    autocorr_peak: float = 0.0       # first autocorrelation peak (residues)
    # Basin-visit pattern
    basin_sequence: str = ""
    basin_counts: dict = field(default_factory=dict)
    parastichy_counts: list = field(default_factory=list)
    return_lengths: list = field(default_factory=list)  # steps to return to same basin
    fibonacci_proximity: float = 0.0  # how close return lengths are to Fibonacci numbers
    # BPS (from TorusFold framework)
    bps_energy: float = 0.0
    bps_per_residue: float = 0.0


def classify_basin(phi_deg, psi_deg):
    """Assign (φ,ψ) point to nearest Ramachandran basin."""
    def torus_dist(a1, a2):
        """Angular distance on torus (degrees)."""
        d1 = abs(a1[0] - a2[0])
        d2 = abs(a1[1] - a2[1])
        d1 = min(d1, 360 - d1)
        d2 = min(d2, 360 - d2)
        return math.sqrt(d1**2 + d2**2)

    point = (phi_deg, psi_deg)
    best_basin = 'L'  # loop/coil default
    best_dist = 60     # threshold: must be within 60° of a basin center

    for name, center in BASINS.items():
        d = torus_dist(point, center)
        if d < best_dist:
            best_dist = d
            best_basin = name

    return best_basin


def angular_advance(angles):
    """Compute per-step angular advance, wrapping correctly on S¹."""
    diff = np.diff(angles)
    # Wrap to [-π, π]
    diff = (diff + np.pi) % (2 * np.pi) - np.pi
    return diff


def winding_number(angles):
    """
    Compute winding number: total unwrapped advance / 2π.
    This counts how many times the trajectory wraps around S¹.
    """
    advances = angular_advance(angles)
    total = np.sum(advances)
    return total / (2 * np.pi)


def spectral_analysis(advances):
    """Analyze frequency content of angular advances."""
    if len(advances) < 8:
        return 0.0, 1.0

    # Remove mean
    centered = advances - np.mean(advances)

    # Power spectrum via FFT
    fft = np.fft.rfft(centered)
    power = np.abs(fft)**2
    freqs = np.fft.rfftfreq(len(centered))

    # Skip DC component
    power = power[1:]
    freqs = freqs[1:]

    if len(power) == 0 or np.sum(power) == 0:
        return 0.0, 1.0

    # Dominant frequency
    dominant_idx = np.argmax(power)
    dominant_freq = freqs[dominant_idx]

    # Spectral entropy (normalized)
    p_norm = power / np.sum(power)
    p_norm = p_norm[p_norm > 0]
    entropy = -np.sum(p_norm * np.log(p_norm)) / np.log(len(p_norm)) if len(p_norm) > 1 else 1.0

    return dominant_freq, entropy


def autocorrelation_peak(advances):
    """Find first significant autocorrelation peak (periodicity estimate)."""
    if len(advances) < 10:
        return 0.0

    centered = advances - np.mean(advances)
    acf = np.correlate(centered, centered, mode='full')
    acf = acf[len(acf)//2:]  # take positive lags only

    if acf[0] == 0:
        return 0.0

    acf = acf / acf[0]  # normalize

    # Find first peak after initial decay
    for i in range(2, len(acf) - 1):
        if acf[i] > acf[i-1] and acf[i] > acf[i+1] and acf[i] > 0.1:
            return float(i)

    return 0.0


def basin_return_analysis(basin_seq):
    """
    Analyze return times: how many steps until each basin type recurs.
    In phyllotaxis, parastichy counts are consecutive Fibonacci numbers.
    """
    if not basin_seq:
        return [], [], 0.0

    # Track last occurrence of each basin
    last_seen = {}
    return_lengths = []

    for i, b in enumerate(basin_seq):
        if b in last_seen:
            return_lengths.append(i - last_seen[b])
        last_seen[b] = i

    if not return_lengths:
        return [], [], 0.0

    # Count frequency of each return length
    counts = Counter(return_lengths)
    parastichy = sorted(counts.items(), key=lambda x: -x[1])

    # Fibonacci sequence for comparison
    fibs = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89]

    # How close are the top return lengths to Fibonacci numbers?
    top_returns = [p[0] for p in parastichy[:5]]
    fib_distances = []
    for r in top_returns:
        min_dist = min(abs(r - f) for f in fibs)
        fib_distances.append(min_dist / max(r, 1))

    fib_proximity = 1.0 - np.mean(fib_distances) if fib_distances else 0.0

    return return_lengths, parastichy, fib_proximity


def compute_bps_energy(phi, psi):
    """
    Compute BPS energy from superpotential W = -√p(φ,ψ).
    Uses von Mises mixture for Ramachandran density (from TorusFold).
    """
    # Von Mises mixture components (simplified 4-basin model)
    components = [
        {'mu_phi': -63*RAD, 'mu_psi': -43*RAD, 'kappa': 8.0, 'weight': 0.33},   # α
        {'mu_phi': -120*RAD, 'mu_psi': 130*RAD, 'kappa': 6.0, 'weight': 0.28},   # β
        {'mu_phi': -75*RAD, 'mu_psi': 145*RAD, 'kappa': 5.0, 'weight': 0.12},    # PPII
        {'mu_phi': 57*RAD, 'mu_psi': 47*RAD, 'kappa': 4.0, 'weight': 0.05},      # αL
        {'mu_phi': 0, 'mu_psi': 0, 'kappa': 0.5, 'weight': 0.22},                 # background
    ]

    def rama_density(ph, ps):
        p = 0.0
        for c in components:
            von_mises = np.exp(
                c['kappa'] * (np.cos(ph - c['mu_phi']) + np.cos(ps - c['mu_psi']))
            )
            p += c['weight'] * von_mises
        return p / (2 * np.pi)  # normalize roughly

    # W = -√p at each residue
    W = np.array([-np.sqrt(rama_density(ph, ps)) for ph, ps in zip(phi, psi)])

    # BPS energy = total variation of W
    bps = np.sum(np.abs(np.diff(W)))

    return bps, W


def analyze_protein(phi, psi, protein_id="unknown"):
    """
    Full phyllotaxis analysis on a single protein's (φ,ψ) trajectory.
    phi, psi: arrays in radians
    Returns: PhyllotaxisResult
    """
    L = len(phi)
    result = PhyllotaxisResult(protein_id=protein_id, length=L)

    if L < 5:
        return result

    # ─── Winding numbers ───
    result.winding_phi = winding_number(phi)
    result.winding_psi = winding_number(psi)

    if abs(result.winding_psi) > 0.01:
        result.winding_ratio = abs(result.winding_phi / result.winding_psi)
    else:
        result.winding_ratio = float('inf')

    # Nearest noble number
    if result.winding_ratio != float('inf'):
        best_name, best_dist = "", float('inf')
        for name, val in NOBLE_NUMBERS.items():
            d = abs(result.winding_ratio - val)
            if d < best_dist:
                best_dist = d
                best_name = name
        result.nearest_noble = best_name
        result.noble_distance = best_dist

    # ─── Angular advances ───
    dphi = angular_advance(phi)
    dpsi = angular_advance(psi)

    result.mean_advance_phi = np.degrees(np.mean(dphi))
    result.mean_advance_psi = np.degrees(np.mean(dpsi))
    result.advance_std_phi = np.degrees(np.std(dphi))
    result.advance_std_psi = np.degrees(np.std(dpsi))

    # Effective angular advance on T²
    eff_angles = np.degrees(np.sqrt(dphi**2 + dpsi**2))
    result.effective_angle = np.mean(eff_angles)

    # Golden angle proximity
    result.golden_proximity = 1.0 - abs(result.effective_angle - GOLDEN_ANGLE_DEG) / 180.0

    # ─── Spectral analysis ───
    # Analyze the magnitude of angular advances
    advance_mag = np.sqrt(dphi**2 + dpsi**2)
    result.dominant_frequency, result.spectral_entropy = spectral_analysis(advance_mag)
    result.autocorr_peak = autocorrelation_peak(advance_mag)

    # ─── Basin-visit sequence ───
    basin_seq = [classify_basin(np.degrees(ph), np.degrees(ps)) for ph, ps in zip(phi, psi)]
    result.basin_sequence = ''.join(b[0] for b in basin_seq)  # compact encoding
    result.basin_counts = dict(Counter(basin_seq))

    # Return-length analysis (parastichy analog)
    result.return_lengths, parastichy, result.fibonacci_proximity = basin_return_analysis(basin_seq)
    result.parastichy_counts = parastichy[:10]  # top 10

    # ─── BPS energy ───
    bps, W = compute_bps_energy(phi, psi)
    result.bps_energy = bps
    result.bps_per_residue = bps / L if L > 0 else 0.0

    return result


# ═══════════════════════════════════════════════════════════════════════
# MARKOV NULL MODEL
# ═══════════════════════════════════════════════════════════════════════

def markov_null_model(phi, psi, n_shuffles=100):
    """
    Generate Markov-shuffled controls preserving transition statistics.
    Tests whether quasi-periodic structure exceeds what nearest-neighbor
    transitions alone produce.
    """
    L = len(phi)
    if L < 10:
        return {}

    dphi = angular_advance(phi)
    dpsi = angular_advance(psi)

    null_golden_prox = []
    null_spectral_ent = []
    null_fib_prox = []

    rng = np.random.default_rng(42)

    for _ in range(n_shuffles):
        # Shuffle the advances (breaks long-range correlations, preserves marginals)
        perm = rng.permutation(len(dphi))
        shuf_dphi = dphi[perm]
        shuf_dpsi = dpsi[perm]

        # Reconstruct angles from shuffled advances
        shuf_phi = np.cumsum(np.concatenate([[phi[0]], shuf_dphi]))
        shuf_psi = np.cumsum(np.concatenate([[psi[0]], shuf_dpsi]))

        # Wrap to [-π, π]
        shuf_phi = (shuf_phi + np.pi) % (2 * np.pi) - np.pi
        shuf_psi = (shuf_psi + np.pi) % (2 * np.pi) - np.pi

        r = analyze_protein(shuf_phi, shuf_psi, "null")
        null_golden_prox.append(r.golden_proximity)
        null_spectral_ent.append(r.spectral_entropy)
        null_fib_prox.append(r.fibonacci_proximity)

    return {
        'golden_proximity': {
            'null_mean': np.mean(null_golden_prox),
            'null_std': np.std(null_golden_prox),
        },
        'spectral_entropy': {
            'null_mean': np.mean(null_spectral_ent),
            'null_std': np.std(null_spectral_ent),
        },
        'fibonacci_proximity': {
            'null_mean': np.mean(null_fib_prox),
            'null_std': np.std(null_fib_prox),
        },
    }


# ═══════════════════════════════════════════════════════════════════════
# VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════

def plot_phyllotaxis_report(results, null_models, output_path):
    """Generate multi-panel analysis figure."""
    fig = plt.figure(figsize=(20, 24))
    fig.suptitle('Phyllotaxis ↔ Protein Backbone Geometry',
                 fontsize=18, fontweight='bold', y=0.98)
    fig.text(0.5, 0.965,
             'Quasi-periodic structure on the Ramachandran torus T²',
             ha='center', fontsize=12, color='#555')

    gs = GridSpec(4, 3, hspace=0.35, wspace=0.3,
                  left=0.07, right=0.95, top=0.94, bottom=0.04)

    # Color scheme
    C_GOLD = '#DAA520'
    C_ALPHA = '#2196F3'
    C_BETA = '#E53935'
    C_PPII = '#4CAF50'
    C_LOOP = '#9E9E9E'
    C_NULL = '#BDBDBD'

    # ─── Panel 1: Winding number ratio distribution ───
    ax1 = fig.add_subplot(gs[0, 0])
    ratios = [r.winding_ratio for r in results if 0 < r.winding_ratio < 10]
    if ratios:
        ax1.hist(ratios, bins=40, color=C_ALPHA, alpha=0.7, edgecolor='white')
        for name, val in NOBLE_NUMBERS.items():
            if val < 10:
                ax1.axvline(val, color=C_GOLD, ls='--', lw=1.5, alpha=0.8)
                ax1.text(val, ax1.get_ylim()[1]*0.9, name, rotation=90,
                        va='top', ha='right', fontsize=7, color=C_GOLD)
    ax1.set_xlabel('Winding ratio |w_φ / w_ψ|')
    ax1.set_ylabel('Count')
    ax1.set_title('Winding Number Ratios on T²')

    # ─── Panel 2: Effective angle vs golden angle ───
    ax2 = fig.add_subplot(gs[0, 1])
    eff_angles = [r.effective_angle for r in results]
    if eff_angles:
        ax2.hist(eff_angles, bins=40, color=C_ALPHA, alpha=0.7, edgecolor='white')
        ax2.axvline(GOLDEN_ANGLE_DEG, color=C_GOLD, lw=2.5, ls='--',
                    label=f'Golden angle ({GOLDEN_ANGLE_DEG:.1f}°)')
        # α-helix effective advance
        alpha_eff = np.degrees(np.sqrt((-63*RAD)**2 + (-43*RAD)**2))
        ax2.axvline(alpha_eff, color=C_ALPHA, lw=2, ls=':',
                    label=f'α-helix advance ({alpha_eff:.1f}°)')
    ax2.set_xlabel('Effective angular advance (°/residue)')
    ax2.set_ylabel('Count')
    ax2.set_title('Angular Advance vs Golden Angle')
    ax2.legend(fontsize=8)

    # ─── Panel 3: Spectral entropy (real vs null) ───
    ax3 = fig.add_subplot(gs[0, 2])
    real_ent = [r.spectral_entropy for r in results]
    null_ent = [nm.get('spectral_entropy', {}).get('null_mean', 1.0) for nm in null_models]
    if real_ent and null_ent:
        ax3.scatter(null_ent, real_ent, c=C_ALPHA, alpha=0.6, s=40, zorder=3)
        lim = [0, 1.05]
        ax3.plot(lim, lim, 'k--', alpha=0.3, label='Real = Null')
        ax3.set_xlim(lim)
        ax3.set_ylim(lim)
    ax3.set_xlabel('Null model spectral entropy')
    ax3.set_ylabel('Real spectral entropy')
    ax3.set_title('Quasi-Periodicity: Real vs Markov Null')
    ax3.legend(fontsize=8)

    # ─── Panel 4: Ramachandran trajectory example (first protein) ───
    ax4 = fig.add_subplot(gs[1, 0])
    if results:
        r0 = results[0]
        # Re-extract angles for plotting if we have them stored
        ax4.text(0.5, 0.5, f'Protein: {r0.protein_id}\nL = {r0.length}\n'
                 f'w_φ = {r0.winding_phi:.2f}\nw_ψ = {r0.winding_psi:.2f}\n'
                 f'ratio = {r0.winding_ratio:.3f}\n'
                 f'nearest noble: {r0.nearest_noble}',
                 ha='center', va='center', transform=ax4.transAxes,
                 fontsize=11, family='monospace',
                 bbox=dict(boxstyle='round', facecolor='lightyellow'))
    ax4.set_title('Example Protein Summary')
    ax4.axis('off')

    # ─── Panel 5: Basin-visit sequence visualization ───
    ax5 = fig.add_subplot(gs[1, 1:])
    if results:
        basin_colors = {'α': C_ALPHA, 'β': C_BETA, 'P': C_PPII,
                       'a': '#9C27B0', 'L': C_LOOP}
        r0 = results[0]
        seq = r0.basin_sequence[:200]  # first 200 residues
        for i, b in enumerate(seq):
            color = basin_colors.get(b, C_LOOP)
            ax5.barh(0, 1, left=i, color=color, edgecolor='none', height=0.6)
        ax5.set_xlim(0, len(seq))
        ax5.set_yticks([])
        ax5.set_xlabel('Residue index')

        # Legend
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=C_ALPHA, label='α-helix'),
                          Patch(facecolor=C_BETA, label='β-strand'),
                          Patch(facecolor=C_PPII, label='PPII'),
                          Patch(facecolor=C_LOOP, label='Loop/coil')]
        ax5.legend(handles=legend_elements, loc='upper right', fontsize=8, ncol=4)
    ax5.set_title('Basin-Visit Sequence (Parastichy Pattern)')

    # ─── Panel 6: Return-length histogram (parastichy analog) ───
    ax6 = fig.add_subplot(gs[2, 0])
    all_returns = []
    for r in results:
        all_returns.extend(r.return_lengths)
    if all_returns:
        max_r = min(max(all_returns), 55)
        bins = np.arange(0.5, max_r + 1.5, 1)
        ax6.hist(all_returns, bins=bins, color=C_ALPHA, alpha=0.7, edgecolor='white')
        # Mark Fibonacci numbers
        fibs = [1, 2, 3, 5, 8, 13, 21, 34]
        for f in fibs:
            if f <= max_r:
                ax6.axvline(f, color=C_GOLD, ls='--', lw=1.2, alpha=0.7)
        ax6.text(0.95, 0.95, 'Gold lines = Fibonacci',
                transform=ax6.transAxes, ha='right', va='top',
                fontsize=8, color=C_GOLD)
    ax6.set_xlabel('Return length (residues)')
    ax6.set_ylabel('Count')
    ax6.set_title('Basin Return Times (cf. Parastichies)')

    # ─── Panel 7: Fibonacci proximity real vs null ───
    ax7 = fig.add_subplot(gs[2, 1])
    real_fib = [r.fibonacci_proximity for r in results]
    null_fib = [nm.get('fibonacci_proximity', {}).get('null_mean', 0.5) for nm in null_models]
    if real_fib and null_fib:
        ax7.scatter(null_fib, real_fib, c=C_GOLD, alpha=0.6, s=40, zorder=3)
        lim = [0, 1.05]
        ax7.plot(lim, lim, 'k--', alpha=0.3)
    ax7.set_xlabel('Null model Fibonacci proximity')
    ax7.set_ylabel('Real Fibonacci proximity')
    ax7.set_title('Fibonacci Structure in Basin Returns')

    # ─── Panel 8: BPS/L distribution ───
    ax8 = fig.add_subplot(gs[2, 2])
    bps_l = [r.bps_per_residue for r in results if r.bps_per_residue > 0]
    if bps_l:
        ax8.hist(bps_l, bins=30, color=C_ALPHA, alpha=0.7, edgecolor='white')
        mean_bps = np.mean(bps_l)
        ax8.axvline(mean_bps, color='red', lw=2, ls='--',
                    label=f'Mean = {mean_bps:.4f}')
        ax8.axvline(0.20, color=C_GOLD, lw=2, ls=':',
                    label='Universal ≈ 0.20')
    ax8.set_xlabel('BPS / L')
    ax8.set_ylabel('Count')
    ax8.set_title('BPS Energy Per Residue')
    ax8.legend(fontsize=8)

    # ─── Panel 9: Golden proximity vs protein length ───
    ax9 = fig.add_subplot(gs[3, 0])
    lengths = [r.length for r in results]
    gold_prox = [r.golden_proximity for r in results]
    if lengths and gold_prox:
        ax9.scatter(lengths, gold_prox, c=C_ALPHA, alpha=0.5, s=30)
        ax9.axhline(0.5, color=C_GOLD, ls=':', alpha=0.5)
    ax9.set_xlabel('Protein length (residues)')
    ax9.set_ylabel('Golden angle proximity')
    ax9.set_title('Golden Proximity vs Length')

    # ─── Panel 10: Autocorrelation peak distribution ───
    ax10 = fig.add_subplot(gs[3, 1])
    acf_peaks = [r.autocorr_peak for r in results if r.autocorr_peak > 0]
    if acf_peaks:
        ax10.hist(acf_peaks, bins=30, color=C_ALPHA, alpha=0.7, edgecolor='white')
        ax10.axvline(3.6, color=C_GOLD, lw=2, ls='--',
                    label='α-helix period (3.6)')
        ax10.axvline(2.0, color=C_BETA, lw=2, ls=':',
                    label='β-strand period (2.0)')
    ax10.set_xlabel('Autocorrelation peak (residues)')
    ax10.set_ylabel('Count')
    ax10.set_title('Periodicity of Angular Advances')
    ax10.legend(fontsize=8)

    # ─── Panel 11: Conceptual diagram ───
    ax11 = fig.add_subplot(gs[3, 2])
    ax11.text(0.5, 0.85, 'THE ANALOGY', ha='center', va='top',
             fontsize=13, fontweight='bold', transform=ax11.transAxes)
    analogy_text = (
        "PHYLLOTAXIS\n"
        "• Leaves at ~137.5° (golden angle)\n"
        "• On cylindrical stem (S¹ × ℝ)\n"
        "• Local auxin transport physics\n"
        "• Parastichy counts → Fibonacci\n"
        "\n"
        "PROTEIN BACKBONE\n"
        "• Residues at (Δφ, Δψ) offset\n"
        "• On Ramachandran torus (S¹ × S¹)\n"
        "• Local steric/H-bond physics\n"
        "• Basin returns → Fibonacci?"
    )
    ax11.text(0.5, 0.45, analogy_text, ha='center', va='center',
             fontsize=9, family='monospace', transform=ax11.transAxes,
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    ax11.axis('off')

    plt.savefig(output_path, dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"  Report saved: {output_path}")


def plot_torus_trajectory(phi, psi, protein_id, output_path):
    """
    Plot protein trajectory on the Ramachandran torus as a
    colored path, with basin regions shaded.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle(f'Torus Trajectory: {protein_id}', fontsize=14, fontweight='bold')

    phi_deg = np.degrees(phi)
    psi_deg = np.degrees(psi)

    # Left panel: Ramachandran plot with trajectory
    colors_by_pos = np.arange(len(phi_deg))

    # Basin shading
    for name, (cp, cs) in BASINS.items():
        circle = plt.Circle((cp, cs), 40, alpha=0.08,
                           color={'α': '#2196F3', 'β': '#E53935',
                                  'PPII': '#4CAF50', 'αL': '#9C27B0'}[name])
        ax1.add_patch(circle)
        ax1.text(cp, cs, name, ha='center', va='center', fontsize=10,
                fontweight='bold', alpha=0.4)

    sc = ax1.scatter(phi_deg, psi_deg, c=colors_by_pos, cmap='viridis',
                     s=8, alpha=0.7, zorder=3)
    ax1.plot(phi_deg, psi_deg, 'k-', alpha=0.15, lw=0.5)

    ax1.set_xlim(-180, 180)
    ax1.set_ylim(-180, 180)
    ax1.set_xlabel('φ (degrees)')
    ax1.set_ylabel('ψ (degrees)')
    ax1.set_title('Ramachandran Trajectory')
    ax1.set_aspect('equal')
    cb = plt.colorbar(sc, ax=ax1, label='Residue index')

    # Right panel: Unwrapped winding (cumulative angle)
    dphi = angular_advance(phi)
    dpsi = angular_advance(psi)
    cum_phi = np.cumsum(dphi) / (2 * np.pi)
    cum_psi = np.cumsum(dpsi) / (2 * np.pi)

    ax2.plot(cum_phi, cum_psi, '-', lw=1.0, alpha=0.7, color='#2196F3')
    ax2.scatter(cum_phi[::10], cum_psi[::10], c=np.arange(len(cum_phi))[::10],
               cmap='viridis', s=15, zorder=3)
    ax2.set_xlabel('Cumulative φ winding (turns)')
    ax2.set_ylabel('Cumulative ψ winding (turns)')
    ax2.set_title('Unwrapped Trajectory on T²')
    ax2.set_aspect('equal')

    # Annotate winding ratio
    if abs(cum_psi[-1]) > 0.01:
        ratio = abs(cum_phi[-1] / cum_psi[-1])
        ax2.text(0.05, 0.95, f'w_φ/w_ψ = {ratio:.3f}',
                transform=ax2.transAxes, fontsize=11, va='top',
                bbox=dict(boxstyle='round', facecolor='lightyellow'))

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"  Trajectory plot saved: {output_path}")


# ═══════════════════════════════════════════════════════════════════════
# FILE DISCOVERY & BATCH PROCESSING
# ═══════════════════════════════════════════════════════════════════════

def find_cached_structures(cache_dir, recursive=True):
    """
    Find all AlphaFold CIF/PDB files in cache directory.
    Searches recursively to handle any directory layout.
    Also finds plain .cif/.pdb files not from AlphaFold.
    """
    cache = Path(cache_dir)
    if not cache.exists():
        # Try common locations
        common_paths = [
            Path('./alphafold_cache'),
            Path('../alphafold_cache'),
            Path.home() / 'alphafold_cache',
            Path.home() / 'Documents' / 'alphafold_cache',
        ]
        for cp in common_paths:
            if cp.exists():
                print(f"  Auto-discovered cache at: {cp}")
                cache = cp
                break
        else:
            return []

    files = []
    glob_fn = cache.rglob if recursive else cache.glob
    for ext in ('*.cif', '*.cif.gz', '*.pdb', '*.pdb.gz'):
        files.extend(glob_fn(ext))

    # Extract UniProt IDs from filenames
    structures = []
    seen_ids = set()
    for f in sorted(files):
        match = re.match(r'AF-([A-Z0-9]+)-F\d+-model', f.name)
        if match:
            uid = match.group(1)
            if uid not in seen_ids:
                structures.append((uid, f))
                seen_ids.add(uid)
        else:
            # Non-AlphaFold file — use filename stem as ID
            stem = f.stem.replace('.cif', '').replace('.pdb', '')
            if stem not in seen_ids:
                structures.append((stem, f))
                seen_ids.add(stem)

    return structures


# ═══════════════════════════════════════════════════════════════════════
# AUTO-DOWNLOAD (when --download flag used)
# ═══════════════════════════════════════════════════════════════════════

# Representative proteins spanning fold classes
DEFAULT_PROTEINS = {
    # Plaxco 23 two-state folders
    'P01088': 'CI2 (α+β)',
    'P0CG47': 'Ubiquitin (α/β)',
    'P03034': 'λ-repressor (all-α)',
    'P02640': 'Villin HP (all-α)',
    'P02836': 'Engrailed HD (all-α)',
    'P06654': 'Protein G (α/β)',
    'P62942': 'FKBP12 (all-β)',
    'P00648': 'Barnase (α+β)',
    'P11540': 'Barstar (α+β)',
    'P0AE67': 'CheY (α/β)',
    'P02185': 'Myoglobin (all-α)',
    'P00004': 'Cytochrome c (all-α)',
    'P61823': 'RNase A (α+β)',
    'P00720': 'T4 lysozyme (α+β)',
    'P08032': 'SH3 domain (all-β)',
    'P24821': 'Tenascin (all-β)',
    # Additional fold diversity
    'P68871': 'Hemoglobin β (all-α)',
    'P01308': 'Insulin (sm. α)',
    'P10636': 'Tau (IDP)',
    'P0A855': 'Trp cage (mini-α)',
    'P69905': 'Hemoglobin α (all-α)',
    'P02768': 'Albumin (all-α, large)',
    'P04637': 'p53 (multi-domain)',
    'P0ABB4': 'ATP synthase α (α/β, large)',
}


def download_structures(cache_dir, uniprot_ids=None, max_proteins=None):
    """Download AlphaFold structures. Requires network access."""
    import urllib.request
    import urllib.error

    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)

    if uniprot_ids is None:
        uniprot_ids = list(DEFAULT_PROTEINS.keys())

    if max_proteins:
        uniprot_ids = uniprot_ids[:max_proteins]

    downloaded = 0
    skipped = 0
    failed = 0

    for uid in uniprot_ids:
        name = DEFAULT_PROTEINS.get(uid, uid)

        # Check if already cached
        found = False
        for v in (4, 3, 2, 1):
            p = cache / f'AF-{uid}-F1-model_v{v}.cif'
            if p.exists() and p.stat().st_size > 100:
                print(f"  CACHED  {uid} ({name})")
                skipped += 1
                found = True
                break
        if found:
            continue

        # Try downloading
        for v in (4, 3, 2):
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v{v}.cif"
            cif_path = cache / f"AF-{uid}-F1-model_v{v}.cif"
            try:
                print(f"  GET     {uid} ({name}) v{v}...", end=' ', flush=True)
                urllib.request.urlretrieve(url, str(cif_path))
                if cif_path.stat().st_size > 100:
                    print(f"OK ({cif_path.stat().st_size:,} bytes)")
                    downloaded += 1
                    break
                else:
                    cif_path.unlink(missing_ok=True)
                    print("empty")
            except urllib.error.HTTPError as e:
                cif_path.unlink(missing_ok=True)
                if e.code == 404:
                    print(f"404")
                    continue
                elif e.code == 429:
                    print(f"rate-limited, waiting 10s...")
                    import time
                    time.sleep(10)
                    break
                else:
                    print(f"HTTP {e.code}")
                    break
            except Exception as e:
                cif_path.unlink(missing_ok=True)
                print(f"error: {e}")
                break
        else:
            failed += 1
            print(f"  FAIL    {uid} ({name})")

    print(f"\n  Download summary: {downloaded} new, {skipped} cached, {failed} failed")
    return downloaded + skipped


def process_structure(filepath, protein_id, run_null=True, null_shuffles=50):
    """Process a single structure file → PhyllotaxisResult + null model."""
    try:
        atoms = extract_backbone_atoms(filepath)
        if len(atoms) < 10:
            return None, None, None, None

        phi, psi, res_idx = compute_dihedrals(atoms)
        if len(phi) < 5:
            return None, None, None, None

        result = analyze_protein(phi, psi, protein_id)

        null = {}
        if run_null:
            null = markov_null_model(phi, psi, n_shuffles=null_shuffles)

        return result, null, phi, psi

    except Exception as e:
        print(f"  ERROR processing {protein_id}: {e}")
        return None, None, None, None


# ═══════════════════════════════════════════════════════════════════════
# SYNTHETIC DEMO (no files needed)
# ═══════════════════════════════════════════════════════════════════════

def generate_synthetic_proteins(n=20):
    """
    Generate synthetic protein-like (φ,ψ) trajectories for demo.
    Mix of α-helical, β-strand, and mixed segments.
    """
    rng = np.random.default_rng(2026)
    proteins = []

    templates = {
        'pure_alpha': {'segments': [('α', 1.0)], 'lengths': (80, 200)},
        'pure_beta':  {'segments': [('β', 1.0)], 'lengths': (60, 150)},
        'alpha_beta': {'segments': [('α', 0.4), ('β', 0.3), ('L', 0.3)], 'lengths': (100, 300)},
        'mixed':      {'segments': [('α', 0.25), ('β', 0.25), ('PPII', 0.15), ('L', 0.35)], 'lengths': (80, 250)},
    }

    basin_params = {
        'α':    (-63*RAD, -43*RAD, 0.15),   # mu_phi, mu_psi, noise_scale
        'β':    (-120*RAD, 130*RAD, 0.20),
        'PPII': (-75*RAD, 145*RAD, 0.18),
        'αL':   (57*RAD, 47*RAD, 0.15),
        'L':    (0, 0, 0.8),  # loop — wide distribution
    }

    for i in range(n):
        # Pick a random template
        tname = rng.choice(list(templates.keys()))
        tmpl = templates[tname]
        L = rng.integers(tmpl['lengths'][0], tmpl['lengths'][1])

        phi = np.zeros(L)
        psi = np.zeros(L)

        pos = 0
        while pos < L:
            # Pick segment type according to weights
            seg_types = [s[0] for s in tmpl['segments']]
            seg_weights = [s[1] for s in tmpl['segments']]
            seg_type = rng.choice(seg_types, p=seg_weights)

            seg_len = rng.integers(5, 30)
            seg_len = min(seg_len, L - pos)

            mu_phi, mu_psi, noise = basin_params[seg_type]
            phi[pos:pos+seg_len] = mu_phi + rng.normal(0, noise, seg_len)
            psi[pos:pos+seg_len] = mu_psi + rng.normal(0, noise, seg_len)

            pos += seg_len

        # Wrap to [-π, π]
        phi = (phi + np.pi) % (2 * np.pi) - np.pi
        psi = (psi + np.pi) % (2 * np.pi) - np.pi

        proteins.append((f'SYNTH_{i:03d}_{tname}', phi, psi))

    return proteins


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def _print_no_structures_help(cache_dir):
    """Print detailed diagnostic when no structures found."""
    cache = Path(cache_dir)
    print(f"\n  ╔══════════════════════════════════════════════════════════╗")
    print(f"  ║  NO STRUCTURE FILES FOUND                                ║")
    print(f"  ╚══════════════════════════════════════════════════════════╝")
    print(f"\n  Searched: {cache.resolve()}")

    if not cache.exists():
        print(f"  ⚠ Directory does not exist!")
        # Check what IS in the current directory
        cwd = Path('.')
        print(f"\n  Current directory: {cwd.resolve()}")
        contents = list(cwd.iterdir())
        if contents:
            print(f"  Contents ({len(contents)} items):")
            for item in sorted(contents)[:20]:
                marker = '📁' if item.is_dir() else '📄'
                print(f"    {marker} {item.name}")
        # Look for any .cif files anywhere nearby
        cifs_nearby = list(cwd.rglob('*.cif'))[:5]
        if cifs_nearby:
            print(f"\n  Found .cif files nearby:")
            for f in cifs_nearby:
                print(f"    → {f}")
            print(f"\n  Try: python phyllotaxis_protein.py --cache \"{cifs_nearby[0].parent}\"")
    else:
        # Directory exists but no matching files
        all_files = list(cache.rglob('*'))
        file_types = Counter(f.suffix for f in all_files if f.is_file())
        print(f"\n  Directory exists with {len(all_files)} items:")
        for ext, count in file_types.most_common(10):
            print(f"    {ext or '(no ext)':>10}: {count}")

        # Check for CIF/PDB in subdirectories
        cifs = list(cache.rglob('*.cif')) + list(cache.rglob('*.cif.gz'))
        pdbs = list(cache.rglob('*.pdb')) + list(cache.rglob('*.pdb.gz'))
        if cifs or pdbs:
            print(f"\n  Found structure files in subdirectories!")
            for f in (cifs + pdbs)[:5]:
                print(f"    → {f}")

    print(f"\n  OPTIONS:")
    print(f"    1. Download proteins:  python phyllotaxis_protein.py --download")
    print(f"    2. Point to cache:     python phyllotaxis_protein.py --cache /path/to/cif/files")
    print(f"    3. Run synthetic demo: python phyllotaxis_protein.py --demo")
    print(f"\n  Expected filenames: AF-<UniProtID>-F1-model_v4.cif")
    print(f"  (Also accepts plain .cif or .pdb files)")


def main():
    parser = argparse.ArgumentParser(
        description='Phyllotaxis ↔ Protein Backbone Geometry Analyzer',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --demo                                     # Synthetic demo
  %(prog)s --cache ./alphafold_cache                  # All cached structures
  %(prog)s --cache ./alphafold_cache --uniprot P0CG47 # Specific protein
  %(prog)s --cache ./alphafold_cache --max 50         # First 50 structures
        """)

    parser.add_argument('--cache', type=str, default='./alphafold_cache',
                       help='Directory containing AlphaFold CIF/PDB files')
    parser.add_argument('--uniprot', type=str, default=None,
                       help='Specific UniProt ID to analyze')
    parser.add_argument('--max', type=int, default=None,
                       help='Max number of proteins to process')
    parser.add_argument('--demo', action='store_true',
                       help='Run with synthetic proteins (no files needed)')
    parser.add_argument('--download', action='store_true',
                       help='Download default protein set from AlphaFold (needs internet)')
    parser.add_argument('--null-shuffles', type=int, default=50,
                       help='Number of Markov null shuffles per protein')
    parser.add_argument('--output', type=str, default='.',
                       help='Output directory for results')
    parser.add_argument('--no-null', action='store_true',
                       help='Skip null model computation (faster)')
    parser.add_argument('--trajectory-plots', action='store_true',
                       help='Generate individual trajectory plots')
    parser.add_argument('--scan', action='store_true',
                       help='Scan cache dir and report what was found (no analysis)')

    args = parser.parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("PHYLLOTAXIS ↔ PROTEIN BACKBONE GEOMETRY ANALYZER")
    print("=" * 70)

    results = []
    null_models = []
    trajectories = {}  # protein_id → (phi, psi)

    # ─── Scan mode ───
    if args.scan:
        print(f"\n  Scanning for structure files...")
        structures = find_cached_structures(args.cache)
        if structures:
            print(f"\n  Found {len(structures)} structures in/under {args.cache}:")
            for uid, fpath in structures[:50]:
                size = fpath.stat().st_size
                print(f"    {uid:>12}  {fpath.name}  ({size:,} bytes)")
            if len(structures) > 50:
                print(f"    ... and {len(structures)-50} more")
        else:
            _print_no_structures_help(args.cache)
        sys.exit(0)

    # ─── Download mode ───
    if args.download:
        print(f"\n  Downloading AlphaFold structures to {args.cache}...")
        n = download_structures(args.cache, max_proteins=args.max)
        if n == 0:
            print("  No structures downloaded. Check your internet connection.")
            sys.exit(1)
        print()  # fall through to analysis

    if args.demo:
        # ─── Synthetic demo ───
        print("\n[DEMO MODE] Generating synthetic proteins...")
        synthetics = generate_synthetic_proteins(20)

        for pid, phi, psi in synthetics:
            print(f"  Analyzing {pid} (L={len(phi)})...")
            r = analyze_protein(phi, psi, pid)
            results.append(r)

            nm = {}
            if not args.no_null:
                nm = markov_null_model(phi, psi, n_shuffles=args.null_shuffles)
            null_models.append(nm)
            trajectories[pid] = (phi, psi)

    else:
        # ─── Real structures from cache ───
        if args.uniprot:
            # Look for specific protein
            cache = Path(args.cache)
            found = None
            for ext in ('.cif', '.cif.gz', '.pdb', '.pdb.gz'):
                for v in (4, 3, 2, 1):
                    p = cache / f'AF-{args.uniprot}-F1-model_v{v}{ext}'
                    if p.exists():
                        found = [(args.uniprot, p)]
                        break
                if found:
                    break
            if not found:
                print(f"  ERROR: No cached structure for {args.uniprot} in {cache}")
                print(f"  Available files: {list(cache.glob('AF-*'))[:5]}")
                sys.exit(1)
            structures = found
        else:
            structures = find_cached_structures(args.cache)

        if not structures:
            _print_no_structures_help(args.cache)
            sys.exit(1)

        if args.max:
            structures = structures[:args.max]

        print(f"\n  Found {len(structures)} structures in {args.cache}")

        for i, (uid, fpath) in enumerate(structures):
            print(f"  [{i+1}/{len(structures)}] {uid}...", end=' ', flush=True)
            r, nm, phi, psi = process_structure(
                fpath, uid,
                run_null=not args.no_null,
                null_shuffles=args.null_shuffles
            )
            if r:
                results.append(r)
                null_models.append(nm or {})
                trajectories[uid] = (phi, psi)
                print(f"L={r.length}, w_φ/w_ψ={r.winding_ratio:.3f}, "
                      f"eff_angle={r.effective_angle:.1f}°, "
                      f"BPS/L={r.bps_per_residue:.4f}")
            else:
                print("SKIPPED")

    if not results:
        print("\nNo proteins successfully processed.")
        sys.exit(1)

    # ─── Summary statistics ───
    print(f"\n{'='*70}")
    print(f"SUMMARY ({len(results)} proteins)")
    print(f"{'='*70}")

    ratios = [r.winding_ratio for r in results if 0 < r.winding_ratio < 100]
    eff_angles = [r.effective_angle for r in results]
    bps_l = [r.bps_per_residue for r in results]
    gold_prox = [r.golden_proximity for r in results]
    fib_prox = [r.fibonacci_proximity for r in results]

    print(f"\n  Winding ratio |w_φ/w_ψ|:")
    print(f"    Mean: {np.mean(ratios):.3f} ± {np.std(ratios):.3f}")
    print(f"    Median: {np.median(ratios):.3f}")
    for name, val in NOBLE_NUMBERS.items():
        near = sum(1 for r in ratios if abs(r - val) < 0.2)
        print(f"    Near {name} ({val:.3f}): {near}/{len(ratios)}")

    print(f"\n  Effective angular advance:")
    print(f"    Mean: {np.mean(eff_angles):.1f}° ± {np.std(eff_angles):.1f}°")
    print(f"    Golden angle = {GOLDEN_ANGLE_DEG:.1f}°")
    print(f"    Mean golden proximity: {np.mean(gold_prox):.3f}")

    print(f"\n  BPS / L:")
    print(f"    Mean: {np.mean(bps_l):.4f} ± {np.std(bps_l):.4f}")

    print(f"\n  Fibonacci proximity of basin returns:")
    print(f"    Mean: {np.mean(fib_prox):.3f} ± {np.std(fib_prox):.3f}")

    if null_models and any(null_models):
        valid_nulls = [nm for nm in null_models if nm]
        if valid_nulls:
            null_fib_means = [nm.get('fibonacci_proximity', {}).get('null_mean', 0.5)
                            for nm in valid_nulls]
            print(f"    Null model mean: {np.mean(null_fib_means):.3f}")
            real_vs_null = np.mean(fib_prox[:len(valid_nulls)]) - np.mean(null_fib_means)
            print(f"    Real - Null: {real_vs_null:+.3f} "
                  f"({'ABOVE' if real_vs_null > 0 else 'BELOW'} null)")

    # ─── Phyllotaxis classification ───
    print(f"\n  PHYLLOTACTIC CLASSIFICATION:")
    for r in results[:5]:  # show first 5
        print(f"    {r.protein_id}: "
              f"ratio={r.winding_ratio:.3f} "
              f"({r.nearest_noble}, d={r.noble_distance:.3f}), "
              f"eff={r.effective_angle:.1f}°, "
              f"S_ent={r.spectral_entropy:.3f}")

    # ─── Generate plots ───
    print(f"\n  Generating report...")
    report_path = out_dir / 'phyllotaxis_report.png'
    plot_phyllotaxis_report(results, null_models, report_path)

    # Individual trajectory plots
    if args.trajectory_plots:
        traj_dir = out_dir / 'trajectories'
        traj_dir.mkdir(exist_ok=True)
        for pid, (phi, psi) in list(trajectories.items())[:10]:
            plot_torus_trajectory(phi, psi, pid,
                                traj_dir / f'{pid}_trajectory.png')

    # ─── Save JSON results ───
    json_path = out_dir / 'phyllotaxis_results.json'
    json_data = {
        'metadata': {
            'date': '2026-02-26',
            'n_proteins': len(results),
            'mode': 'demo' if args.demo else 'alphafold_cache',
            'golden_angle_deg': GOLDEN_ANGLE_DEG,
        },
        'summary': {
            'winding_ratio_mean': float(np.mean(ratios)) if ratios else 0,
            'winding_ratio_std': float(np.std(ratios)) if ratios else 0,
            'effective_angle_mean': float(np.mean(eff_angles)),
            'effective_angle_std': float(np.std(eff_angles)),
            'bps_per_residue_mean': float(np.mean(bps_l)),
            'bps_per_residue_std': float(np.std(bps_l)),
            'golden_proximity_mean': float(np.mean(gold_prox)),
            'fibonacci_proximity_mean': float(np.mean(fib_prox)),
        },
        'proteins': [],
    }

    for r in results:
        d = asdict(r)
        # Clean up non-serializable items
        d['parastichy_counts'] = [(k, v) for k, v in d.get('parastichy_counts', [])]
        d.pop('return_lengths', None)  # too verbose for JSON
        json_data['proteins'].append(d)

    with open(json_path, 'w') as f:
        json.dump(json_data, f, indent=2, default=str)
    print(f"  Results saved: {json_path}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()

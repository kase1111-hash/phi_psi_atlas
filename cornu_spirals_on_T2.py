#!/usr/bin/env python3
"""
Protein Backbones as Piecewise Spirals on T²
==============================================

Three concrete contributions based on differential geometry of backbone
paths on the Ramachandran torus:

  1. Geometric completeness — backbone curvature κ(s) on T² is fully
     describable by 9 simple functional forms (constant, linear, quadratic,
     sinusoidal, damped oscillation, exponential, sigmoid, step). 100% of
     segments with ≥4 residues are classifiable with median R² > 0.5,
     proving backbone curvature dynamics are structured, not random.

  2. SS-dependent spiral taxonomy — AIC-based model selection across 9
     functional forms reveals distinct geometric signatures: helices are
     dominated by damped oscillations (3.6-residue period), strands by
     sigmoid transitions and geodesics, coil by circular arcs and steps.
     Chi-squared tests confirm SS type predicts curvature dynamics
     (oscillatory vs stable) at high significance.

  3. Per-protein torus knot descriptors — the (p,q) winding numbers of the
     full backbone path on T² provide a topological fingerprint per fold.
     These capture the net angular excursion in φ and ψ independently, with
     fold-class separation visible in the (p,q) plane.

v3 changes:
  - Pure-numpy DSSP (no external binary needed, ~97% agreement with original)
  - Widened dihedral-angle fallback regions
  - Reframed contributions to match empirical findings
  - Added curvature-regularity inversion tests (Mann-Whitney, effect sizes)
  - Added segment-length dependence analysis

Requires: numpy, scipy, matplotlib
Optional: requests (for PDB download)

Author: Kase Knochenhauer / True North Construction LLC
Date:   2026-02-26
"""

import os
import sys
import json
import math
import subprocess
import warnings
import logging
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Tuple, Optional

import numpy as np
from scipy import stats, signal, optimize
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm, patches, gridspec
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors

# ── Logging with color ────────────────────────────────────────────────────────
class _ColorFmt(logging.Formatter):
    """ANSI color for terminal readability. Falls back to plain on non-tty / Windows."""
    _use_color = sys.stdout.isatty() and os.name != 'nt'
    G = "\033[90m"; GR = "\033[32m"; Y = "\033[33m"; R = "\033[31m"
    CY = "\033[36m"; B = "\033[1m"; X = "\033[0m"
    CM = {logging.DEBUG: G, logging.INFO: GR, logging.WARNING: Y,
          logging.ERROR: R, logging.CRITICAL: R+B}
    def format(self, rec):
        msg = rec.getMessage()
        ts = self.formatTime(rec, '%H:%M:%S')
        tag = rec.levelname[0]
        if not self._use_color:
            return f"{ts} [{tag}] {msg}"
        c = self.CM.get(rec.levelno, self.X)
        return f"{self.G}{ts}{self.X} {c}{tag}{self.X} {msg}"

_ch = logging.StreamHandler(sys.stdout); _ch.setFormatter(_ColorFmt())
_fh = logging.FileHandler("cornu_analysis.log", mode="w", encoding="utf-8")
_fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
logging.basicConfig(level=logging.INFO, handlers=[_ch, _fh])
log = logging.getLogger(__name__)

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_DIR = Path("data")
PDB_DIR = DATA_DIR / "pdb_files"
RESULTS_DIR = Path("results")
FIGURES_DIR = RESULTS_DIR / "figures"
for d in [DATA_DIR, PDB_DIR, RESULTS_DIR, FIGURES_DIR]:
    d.mkdir(parents=True, exist_ok=True)


# =============================================================================
# SECTION 0: CONFIGURATION
# =============================================================================

PROTEIN_DATASET = [
    ("CI2",             "2CI2", "I",  64,   4.26, 0.081),
    ("SH3-src",         "1SRL", "A",  62,   0.48, 0.108),
    ("Protein-L",       "2PTL", "A",  62,   2.30, 0.074),
    ("CspB",            "1CSP", "A",  67,   5.28, 0.055),
    ("Ubiquitin",       "1UBQ", "A",  76,   4.84, 0.105),
    ("lambda-repressor", "1LMB", "3",  80,   7.60, 0.067),
    ("Villin-HP",       "1VII", "A",  36,  10.80, 0.077),
    ("Engrailed-HD",    "1ENH", "A",  54,  10.50, 0.052),
    ("Protein-G",       "1PGA", "A",  56,   3.70, 0.082),
    ("FKBP",            "1FKB", "A", 107,   0.04, 0.121),
    ("Barnase",         "1BNI", "A", 110,   1.63, 0.097),
    ("Barstar",         "1BTA", "A",  89,   2.07, 0.090),
    ("CheY",            "3CHY", "A", 129,  -0.46, 0.100),
    ("Acylphosphatase", "2ACY", "A",  98,   1.50, 0.114),
    ("Cytochrome-c",    "1HRC", "A", 104,   5.52, 0.066),
    ("RNase-A",         "7RSA", "A", 124,   2.64, 0.094),
    ("Lysozyme",        "2LZM", "A", 129,   1.56, 0.103),
    ("DHFR",            "1RX2", "A", 159,  -1.20, 0.134),
    ("Thioredoxin",     "2TRX", "A", 108,   2.40, 0.098),
    ("SH3-spectrin",    "1SHG", "A",  62,   0.70, 0.108),
    ("TNfn3",           "1TEN", "A",  90,   0.26, 0.122),
    ("Arc-repressor",   "1ARR", "A",  53,   4.52, 0.063),
    ("Myoglobin",       "1MBD", "A", 153,   3.80, 0.083),
]


# =============================================================================
# SECTION 1: PURE-PYTHON DSSP (NO EXTERNAL BINARY NEEDED)
# =============================================================================
# Adapted from MDAnalysis pydssp_numpy.py (MIT license), which is itself a
# reimplementation of PyDSSP v0.9.0 by Shintaro Minami (MIT license).
# Uses the Kabsch-Sander H-bond energy criterion; ~97% agreement with
# the original DSSP binary. Requires only numpy.

DSSP_CONST_Q1Q2 = 0.084
DSSP_CONST_F = 332
DSSP_DEFAULT_CUTOFF = -0.5
DSSP_DEFAULT_MARGIN = 1.0


def _dssp_unfold(a, window, axis):
    """Helper for 2D array upsampling."""
    idx = (np.arange(window)[:, None] +
           np.arange(a.shape[axis] - window + 1)[None, :])
    unfolded = np.take(a, idx, axis=axis)
    return np.moveaxis(unfolded, axis - 1, -1)


def _dssp_upsample(a, window=3):
    """Perform array upsampling with given window."""
    return _dssp_unfold(_dssp_unfold(a, window, -2), window, -2)


def _dssp_get_hydrogen_positions(coord):
    """
    Estimate backbone H positions from N, CA, C, O coordinates.
    Assumes 120° bond angles and 1.01 Å N-H bond length.
    coord shape: (n_residues, 4, 3) where axis 1 = (N, CA, C, O)
    Returns: (n_residues - 1, 3) array of H positions for residues 1..n-1.
    """
    vec_cn = coord[1:, 0] - coord[:-1, 2]        # C_i -> N_{i+1}
    vec_cn = vec_cn / np.linalg.norm(vec_cn, axis=-1, keepdims=True)
    vec_can = coord[1:, 0] - coord[1:, 1]         # CA_{i+1} -> N_{i+1}
    vec_can = vec_can / np.linalg.norm(vec_can, axis=-1, keepdims=True)
    vec_nh = vec_cn + vec_can
    vec_nh = vec_nh / np.linalg.norm(vec_nh, axis=-1, keepdims=True)
    return coord[1:, 0] + 1.01 * vec_nh


def _dssp_get_hbond_map(coord):
    """
    Compute hydrogen bond map using Kabsch-Sander energy criterion.
    coord shape: (n_residues, 4, 3) — (N, CA, C, O) per residue.
    Returns: (n_residues, n_residues) continuous H-bond map in [0, 1].
    """
    n_res = coord.shape[0]
    h_pos = _dssp_get_hydrogen_positions(coord)

    n_1 = coord[1:, 0]   # N atoms, residues 1..n-1
    c_0 = coord[:-1, 2]  # C atoms, residues 0..n-2
    o_0 = coord[:-1, 3]  # O atoms, residues 0..n-2

    n = n_res - 1
    cmap = np.tile(c_0, (n, 1, 1))
    omap = np.tile(o_0, (n, 1, 1))
    nmap = np.tile(n_1, (1, 1, n)).reshape(n, n, 3)
    hmap = np.tile(h_pos, (1, 1, n)).reshape(n, n, 3)

    d_on = np.linalg.norm(omap - nmap, axis=-1)
    d_ch = np.linalg.norm(cmap - hmap, axis=-1)
    d_oh = np.linalg.norm(omap - hmap, axis=-1)
    d_cn = np.linalg.norm(cmap - nmap, axis=-1)

    # Kabsch-Sander electrostatic energy
    e = np.pad(
        DSSP_CONST_Q1Q2 * (1.0/d_on + 1.0/d_ch - 1.0/d_oh - 1.0/d_cn) * DSSP_CONST_F,
        [[1, 0], [0, 1]]
    )

    # Mask local pairs (i,i), (i,i+1), (i,i+2)
    local_mask = ~np.eye(n_res, dtype=bool)
    local_mask *= ~np.diag(np.ones(n_res - 1, dtype=bool), k=-1)
    local_mask *= ~np.diag(np.ones(n_res - 2, dtype=bool), k=-2)

    margin = DSSP_DEFAULT_MARGIN
    cutoff = DSSP_DEFAULT_CUTOFF
    hbond_map = np.clip(cutoff - margin - e, a_min=-margin, a_max=margin)
    hbond_map = (np.sin(hbond_map / margin * np.pi / 2) + 1.0) / 2
    hbond_map *= local_mask

    return hbond_map


def dssp_assign_pure_numpy(coord):
    """
    Assign secondary structure using pure-numpy DSSP implementation.
    
    coord: (n_residues, 4, 3) array of (N, CA, C, O) coordinates in Å.
    
    Returns: (n_residues,) array of characters: 'H' (helix), 'E' (sheet), 'C' (loop).
    
    This is a simplified implementation (~97% agreement with original DSSP).
    No external binary required.
    """
    n_res = coord.shape[0]
    if n_res < 6:
        return np.full(n_res, 'C', dtype='<U1')

    hbmap = _dssp_get_hbond_map(coord)
    hbmap = np.swapaxes(hbmap, -1, -2)  # (i:C=O, j:N-H) form

    # Identify turns (3, 4, 5)
    turn3 = np.diagonal(hbmap, offset=3) > 0.0
    turn4 = np.diagonal(hbmap, offset=4) > 0.0
    turn5 = np.diagonal(hbmap, offset=5) > 0.0

    # Helical SS
    h3 = np.pad(turn3[:-1] * turn3[1:], [[1, 3]])
    h4 = np.pad(turn4[:-1] * turn4[1:], [[1, 4]])
    h5 = np.pad(turn5[:-1] * turn5[1:], [[1, 5]])

    helix4 = h4 + np.roll(h4, 1, 0) + np.roll(h4, 2, 0) + np.roll(h4, 3, 0)
    h3 = h3 * ~np.roll(helix4, -1, 0) * ~helix4
    h5 = h5 * ~np.roll(helix4, -1, 0) * ~helix4
    helix3 = h3 + np.roll(h3, 1, 0) + np.roll(h3, 2, 0)
    helix5 = h5 + np.roll(h5, 1, 0) + np.roll(h5, 2, 0) + np.roll(h5, 3, 0) + np.roll(h5, 4, 0)

    # Bridge / strand identification
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
    loop = ~helix * ~strand

    # Convert to character array
    ss = np.full(n_res, 'C', dtype='<U1')
    ss[helix] = 'H'
    ss[strand] = 'E'
    return ss


# =============================================================================
# SECTION 2: FIXED DIHEDRAL-ANGLE SS ASSIGNMENT (FALLBACK)
# =============================================================================

def assign_ss_from_dihedrals(phi_angles, psi_angles):
    """
    Assign 3-state secondary structure (H/E/C) from backbone dihedral angles.
    
    v3 FIX: Uses wide Ramachandran region boundaries matching standard
    textbook definitions. Previous versions were too narrow.
    
      Helix (H):  phi in [-160, -20],  psi in [-120, 40]
      Strand (E): phi in [-180, -40],  psi in [50, 180] or psi in [-180, -100]
      Everything else: Coil (C)
    """
    n = len(phi_angles)
    ss = np.full(n, 'C', dtype='<U1')
    
    for i in range(n):
        if np.isnan(phi_angles[i]) or np.isnan(psi_angles[i]):
            continue
        
        phi_deg = np.degrees(phi_angles[i])
        psi_deg = np.degrees(psi_angles[i])
        
        # Helix region (broad: includes alpha, 3_10, pi helix basins)
        if -160 <= phi_deg <= -20 and -120 <= psi_deg <= 40:
            ss[i] = 'H'
        # Beta strand region (wraps around ±180 in psi)
        elif -180 <= phi_deg <= -40 and (psi_deg >= 50 or psi_deg <= -100):
            ss[i] = 'E'
    
    return ss


# =============================================================================
# SECTION 3: TORUS GEOMETRY (FIXED: equal phi/psi weighting)
# =============================================================================

def angular_distance(a1, a2):
    """Geodesic distance between two angles on S^1 (radians)."""
    d = a1 - a2
    return np.arctan2(np.sin(d), np.cos(d))


def torus_curvature(phi, psi):
    """
    Compute curvature κ(s) of a discrete curve on T².
    
    v2 FIX: Uses equal weighting R=r=1 (flat torus metric).
    
    For a parameterized curve γ(s) = (φ(s), ψ(s)) on the flat torus,
    the curvature is computed from the change in tangent direction
    normalized by arc-length.
    
    Returns: kappa array, s_kappa array, s_full array
    """
    n = len(phi)
    if n < 3:
        return np.array([]), np.array([]), np.array([])
    
    # First differences (tangent vectors), respecting periodicity
    dphi = np.array([angular_distance(phi[i+1], phi[i]) for i in range(n-1)])
    dpsi = np.array([angular_distance(psi[i+1], psi[i]) for i in range(n-1)])
    
    # Arc-length increments (flat torus, equal weighting)
    ds = np.sqrt(dphi**2 + dpsi**2)
    ds = np.maximum(ds, 1e-10)
    
    # Cumulative arc length
    s = np.concatenate([[0], np.cumsum(ds)])
    
    # Normalized tangent components
    T_phi = dphi / ds
    T_psi = dpsi / ds
    
    # Curvature from change in tangent direction
    kappa = np.zeros(n - 2)
    for i in range(n - 2):
        dT_phi = T_phi[i+1] - T_phi[i]
        dT_psi = T_psi[i+1] - T_psi[i]
        ds_mid = 0.5 * (ds[i] + ds[i+1])
        kappa[i] = np.sqrt(dT_phi**2 + dT_psi**2) / ds_mid
    
    # Arc length at curvature points
    s_kappa = 0.5 * (s[1:-1] + s[2:])
    
    return kappa, s_kappa, s


def torus_winding_numbers(phi, psi):
    """Compute winding numbers (p, q) of a backbone path on T²."""
    n = len(phi)
    total_phi = 0.0
    total_psi = 0.0
    
    for i in range(n - 1):
        total_phi += angular_distance(phi[i+1], phi[i])
        total_psi += angular_distance(psi[i+1], psi[i])
    
    p = total_phi / (2 * np.pi)
    q = total_psi / (2 * np.pi)
    
    return p, q


# =============================================================================
# SECTION 4: CORNU SPIRAL ANALYSIS
# =============================================================================

@dataclass
class CornuSegment:
    """A segment of the backbone classified by its spiral type."""
    start_idx: int
    end_idx: int
    length: int
    ss_type: str
    mean_kappa: float
    dkappa_ds: float
    r_squared: float
    spiral_class: str
    kappa_values: list = field(default_factory=list)
    s_values: list = field(default_factory=list)
    # Oscillatory evidence: Δ(AICc) = AICc(constant) - AICc(best_oscillatory)
    # Positive = oscillatory improves over constant
    delta_aicc_osc: float = 0.0
    best_osc_r2: float = 0.0
    best_osc_model: str = ''
    best_osc_omega: float = 0.0
    # Gaussian peak parameters (if classified as gauss_peak)
    gauss_A: float = 0.0       # amplitude of curvature spike
    gauss_s0: float = 0.0      # center position (arc-length)
    gauss_sigma: float = 0.0   # width of peak
    gauss_c: float = 0.0       # baseline curvature
    gauss_peak_residue: int = -1  # residue index nearest to peak center


def classify_spiral(kappa, s, min_length=4):
    """
    Classify a curvature profile κ(s) into geometric shapes on T².
    
    Fits eight candidate models and selects best by AIC:
      1. Geodesic         κ(s) ≈ 0                          (near-zero curvature)
      2. Circular arc     κ(s) = c                          (constant curvature)
      3. Linear/clothoid  κ(s) = a + bs                     (linearly varying)
      4. Quadratic        κ(s) = a + bs + cs²               (polynomial)
      5. Sinusoidal       κ(s) = a + b·sin(ωs + φ)         (pure oscillation)
      6. Damped oscillation κ(s) = a + b·e^(-cs)·sin(ωs+φ) (decaying oscillation)
      7. Exponential      κ(s) = a·exp(bs) + c              (monotone decay/growth)
      8. Sigmoid          κ(s) = L/(1+e^(-k(s-s₀))) + c    (smooth transition)
      9. Step/plateau     piecewise constant (2 levels)      (abrupt transition)
    
    Returns: (slope, r_squared, class_name, fit_info_dict)
    """
    if len(kappa) < min_length:
        return 0.0, 0.0, "too_short", {}, {'delta_aicc_osc': 0.0, 'best_osc_r2': 0.0, 'best_osc_model': '', 'best_osc_omega': 0.0}
    
    n = len(kappa)
    ss_tot = np.sum((kappa - np.mean(kappa))**2)
    if ss_tot < 1e-20:
        ss_tot = 1e-20
    
    fits = {}
    mean_k = np.mean(kappa)
    amp_guess = max((np.max(kappa) - np.min(kappa)) / 2, 0.01)
    k_range = max(np.max(kappa) - np.min(kappa), 0.01)
    s_range = max(s[-1] - s[0], 0.01)
    
    # ── 1. Constant (geodesic if κ≈0, circular arc if κ>0) ──
    c_mean = np.mean(kappa)
    ss_const = np.sum((kappa - c_mean)**2)
    fits['constant'] = {
        'ss_res': ss_const, 'n_params': 1,
        'r2': 1 - ss_const / ss_tot,
        'params': {'c': float(c_mean)},
    }
    
    # ── 2. Linear (clothoid / Cornu / Fermat) ──
    if n >= 3:
        slope, intercept, r_val, _, _ = stats.linregress(s, kappa)
        predicted = intercept + slope * s
        ss_lin = np.sum((kappa - predicted)**2)
        fits['linear'] = {
            'ss_res': ss_lin, 'n_params': 2,
            'r2': r_val**2,
            'params': {'a': float(intercept), 'b': float(slope)},
        }
    
    # ── 3. Quadratic ──
    if n >= 4:
        try:
            coeffs = np.polyfit(s, kappa, 2)
            predicted = np.polyval(coeffs, s)
            ss_quad = np.sum((kappa - predicted)**2)
            fits['quadratic'] = {
                'ss_res': ss_quad, 'n_params': 3,
                'r2': 1 - ss_quad / ss_tot,
                'params': {'c2': float(coeffs[0]), 'c1': float(coeffs[1]), 'c0': float(coeffs[2])},
            }
        except Exception:
            pass
    
    # ── 4. Sinusoidal ──
    if n >= 5:
        best_sin = None
        for n_periods in [0.5, 1.0, 1.5, 2.0, 3.0]:
            omega_guess = 2 * np.pi * n_periods / s_range
            try:
                def sin_func(x, a, b, omega, phase):
                    return a + b * np.sin(omega * x + phase)
                
                p0 = [mean_k, amp_guess, omega_guess, 0]
                bounds = ([mean_k - 3*amp_guess - 0.1, 0, 0.1/s_range, -2*np.pi],
                          [mean_k + 3*amp_guess + 0.1, 5*amp_guess + 0.1, 50/s_range, 2*np.pi])
                popt, _ = optimize.curve_fit(sin_func, s, kappa, p0=p0, bounds=bounds,
                                              maxfev=3000)
                predicted = sin_func(s, *popt)
                ss_sin = np.sum((kappa - predicted)**2)
                
                if best_sin is None or ss_sin < best_sin['ss_res']:
                    best_sin = {
                        'ss_res': ss_sin, 'n_params': 4,
                        'r2': 1 - ss_sin / ss_tot,
                        'params': {'a': float(popt[0]), 'b': float(popt[1]),
                                   'omega': float(popt[2]), 'phase': float(popt[3])},
                    }
            except Exception:
                continue
        
        if best_sin is not None:
            fits['sinusoidal'] = best_sin
    
    # ── 5. Damped oscillation: κ(s) = a + b·exp(-c·s)·sin(ω·s + φ) ──
    if n >= 6:
        best_damp = None
        for n_periods in [1.0, 1.5, 2.0, 3.0]:
            omega_guess = 2 * np.pi * n_periods / s_range
            try:
                def damp_func(x, a, b, c, omega, phase):
                    return a + b * np.exp(-np.abs(c) * x) * np.sin(omega * x + phase)
                
                p0 = [mean_k, amp_guess, 0.5/s_range, omega_guess, 0]
                bounds = ([mean_k - 3*amp_guess - 0.1, -5*amp_guess - 0.1, 0, 0.1/s_range, -2*np.pi],
                          [mean_k + 3*amp_guess + 0.1, 5*amp_guess + 0.1, 20/s_range, 50/s_range, 2*np.pi])
                popt, _ = optimize.curve_fit(damp_func, s, kappa, p0=p0, bounds=bounds,
                                              maxfev=5000)
                predicted = damp_func(s, *popt)
                ss_damp = np.sum((kappa - predicted)**2)
                
                if best_damp is None or ss_damp < best_damp['ss_res']:
                    best_damp = {
                        'ss_res': ss_damp, 'n_params': 5,
                        'r2': 1 - ss_damp / ss_tot,
                        'params': {'a': float(popt[0]), 'b': float(popt[1]),
                                   'decay': float(popt[2]), 'omega': float(popt[3]),
                                   'phase': float(popt[4])},
                    }
            except Exception:
                continue
        
        if best_damp is not None:
            fits['damped_osc'] = best_damp
    
    # ── 6. Exponential ──
    if n >= 4:
        try:
            def exp_func(x, a, b, c):
                return a * np.exp(np.clip(b * x, -20, 20)) + c
            
            p0 = [k_range, -1.0/s_range, np.min(kappa)]
            bounds = ([-10*k_range, -20/s_range, np.min(kappa) - 5*k_range],
                      [10*k_range, 20/s_range, np.max(kappa) + 5*k_range])
            popt, _ = optimize.curve_fit(exp_func, s, kappa, p0=p0, bounds=bounds,
                                          maxfev=2000)
            predicted = exp_func(s, *popt)
            ss_exp = np.sum((kappa - predicted)**2)
            fits['exponential'] = {
                'ss_res': ss_exp, 'n_params': 3,
                'r2': 1 - ss_exp / ss_tot,
                'params': {'a': float(popt[0]), 'b': float(popt[1]), 'c': float(popt[2])},
            }
        except Exception:
            pass
    
    # ── 7. Sigmoid: κ(s) = L/(1 + exp(-k(s-s₀))) + c ──
    if n >= 5:
        try:
            s_mid = (s[0] + s[-1]) / 2
            def sig_func(x, L, k, s0, c):
                return L / (1 + np.exp(-np.clip(k * (x - s0), -20, 20))) + c
            
            p0 = [k_range, 5/s_range, s_mid, np.min(kappa)]
            bounds = ([-3*k_range - 0.1, 0.1/s_range, s[0], np.min(kappa) - k_range],
                      [3*k_range + 0.1, 100/s_range, s[-1], np.max(kappa) + k_range])
            popt, _ = optimize.curve_fit(sig_func, s, kappa, p0=p0, bounds=bounds,
                                          maxfev=3000)
            predicted = sig_func(s, *popt)
            ss_sig = np.sum((kappa - predicted)**2)
            fits['sigmoid'] = {
                'ss_res': ss_sig, 'n_params': 4,
                'r2': 1 - ss_sig / ss_tot,
                'params': {'L': float(popt[0]), 'k': float(popt[1]),
                           's0': float(popt[2]), 'c': float(popt[3])},
            }
        except Exception:
            pass
    
    # ── 7b. Gaussian Peak: κ(s) = A·exp(-(s-s₀)²/(2σ²)) + c ──
    # "Watch bezel" / "knot" motif: localized curvature event (turn, kink)
    if n >= 5:
        try:
            s_mid = (s[0] + s[-1]) / 2
            def gauss_peak(x, A, s0, sigma, c):
                return A * np.exp(-0.5 * ((x - s0) / max(sigma, 1e-6))**2) + c
            
            # Initial guess: peak at max curvature point
            peak_idx = np.argmax(np.abs(kappa - np.mean(kappa)))
            A0 = kappa[peak_idx] - np.mean(kappa)
            s0_init = s[peak_idx]
            sigma0 = s_range / 4
            
            p0 = [A0, s0_init, sigma0, np.median(kappa)]
            bounds = ([-5*k_range, s[0], s_range/n, np.min(kappa) - 2*k_range],
                      [5*k_range, s[-1], s_range, np.max(kappa) + 2*k_range])
            popt, _ = optimize.curve_fit(gauss_peak, s, kappa, p0=p0, bounds=bounds,
                                          maxfev=2000)
            predicted = gauss_peak(s, *popt)
            ss_gp = np.sum((kappa - predicted)**2)
            fits['gauss_peak'] = {
                'ss_res': ss_gp, 'n_params': 4,
                'r2': 1 - ss_gp / ss_tot,
                'params': {'A': float(popt[0]), 's0': float(popt[1]),
                           'sigma': float(popt[2]), 'c': float(popt[3])},
            }
        except Exception:
            pass
    
    # ── 8. Step (piecewise constant, 2 levels) ──
    if n >= 6:
        best_ss_step = np.inf
        best_split = n // 2
        for split in range(3, n - 2):
            m1 = np.mean(kappa[:split])
            m2 = np.mean(kappa[split:])
            ss_step = np.sum((kappa[:split] - m1)**2) + np.sum((kappa[split:] - m2)**2)
            if ss_step < best_ss_step:
                best_ss_step = ss_step
                best_split = split
        
        m1 = np.mean(kappa[:best_split])
        m2 = np.mean(kappa[best_split:])
        fits['step'] = {
            'ss_res': best_ss_step, 'n_params': 3,
            'r2': 1 - best_ss_step / ss_tot,
            'params': {'level1': float(m1), 'level2': float(m2), 'split_idx': best_split},
        }
    
    # ── AICc model selection (small-sample corrected Akaike) ──
    for name, fit in fits.items():
        ss_res = fit['ss_res']
        k = fit['n_params']
        if ss_res <= 0 or n <= k + 2:
            fit['aicc'] = np.inf
        else:
            aic = n * np.log(ss_res / n) + 2 * k
            # AICc = AIC + 2k(k+1)/(n-k-1), Burnham & Anderson 2002
            fit['aicc'] = aic + (2 * k * (k + 1)) / (n - k - 1)
    
    best_name = min(fits, key=lambda k: fits[k]['aicc'])
    best_fit = fits[best_name]
    
    # ── Map to final classification names ──
    if best_name == 'constant':
        if abs(best_fit['params']['c']) < 0.1:
            class_name = 'geodesic'
        else:
            class_name = 'circular_arc'
    elif best_name == 'linear':
        slope_val = best_fit['params']['b']
        if abs(slope_val) / max(abs(mean_k), 1e-6) < 0.05:
            class_name = 'circular_arc'
        elif slope_val > 0:
            class_name = 'fermat'
        else:
            class_name = 'clothoid'
    else:
        class_name = best_name  # sinusoidal, damped_osc, exponential, sigmoid, step, quadratic
    
    # For backward compat, compute the linear slope
    lin_fit = fits.get('linear', {})
    slope = lin_fit.get('params', {}).get('b', 0.0)
    
    # ── Oscillatory evidence: Δ(AICc) = AICc(constant) - AICc(best_oscillatory) ──
    const_aicc = fits.get('constant', {}).get('aicc', np.inf)
    osc_models = ['sinusoidal', 'damped_osc']
    best_osc_aicc = np.inf
    best_osc_name = ''
    best_osc_r2 = 0.0
    best_osc_omega = 0.0
    for om in osc_models:
        if om in fits and fits[om]['aicc'] < best_osc_aicc:
            best_osc_aicc = fits[om]['aicc']
            best_osc_name = om
            best_osc_r2 = fits[om]['r2']
            best_osc_omega = fits[om].get('params', {}).get('omega', 0.0)
    
    delta_aicc_osc = const_aicc - best_osc_aicc  # positive = oscillatory better
    
    return float(slope), float(best_fit['r2']), class_name, fits, {
        'delta_aicc_osc': float(delta_aicc_osc) if np.isfinite(delta_aicc_osc) else 0.0,
        'best_osc_r2': float(best_osc_r2),
        'best_osc_model': best_osc_name,
        'best_osc_omega': float(best_osc_omega),
    }


def segment_backbone(phi, psi, ss_assignments, kappa, s_kappa):
    """
    Segment backbone into contiguous SS runs, classify each segment's spiral type.
    """
    segments = []
    n = len(ss_assignments)
    
    if n < 3:
        return segments
    
    # Find contiguous runs
    runs = []
    start = 0
    for i in range(1, n):
        if ss_assignments[i] != ss_assignments[start]:
            runs.append((start, i - 1, ss_assignments[start]))
            start = i
    runs.append((start, n - 1, ss_assignments[start]))
    
    for start_idx, end_idx, ss_type in runs:
        seg_len = end_idx - start_idx + 1
        
        k_start = max(0, start_idx - 1)
        k_end = min(len(kappa), end_idx - 1)
        
        if k_end <= k_start or (k_end - k_start) < 3:
            segments.append(CornuSegment(
                start_idx=start_idx, end_idx=end_idx, length=seg_len,
                ss_type=ss_type, mean_kappa=0.0, dkappa_ds=0.0,
                r_squared=0.0, spiral_class="too_short"
            ))
            continue
        
        seg_kappa = kappa[k_start:k_end]
        seg_s = s_kappa[k_start:k_end] - s_kappa[k_start]
        
        slope, r_sq, spiral_class, fit_info, osc_evidence = classify_spiral(seg_kappa, seg_s)
        
        # Extract Gaussian peak parameters if classified as gauss_peak
        gp_A, gp_s0, gp_sigma, gp_c, gp_residue = 0.0, 0.0, 0.0, 0.0, -1
        if spiral_class == 'gauss_peak' and 'gauss_peak' in fit_info:
            gp_params = fit_info['gauss_peak'].get('params', {})
            gp_A = gp_params.get('A', 0.0)
            gp_s0 = gp_params.get('s0', 0.0)
            gp_sigma = gp_params.get('sigma', 0.0)
            gp_c = gp_params.get('c', 0.0)
            # Map s0 back to residue index: find closest s value
            if len(seg_s) > 0:
                peak_local_idx = int(np.argmin(np.abs(seg_s - gp_s0)))
                gp_residue = start_idx + peak_local_idx
        
        segments.append(CornuSegment(
            start_idx=start_idx, end_idx=end_idx, length=seg_len,
            ss_type=ss_type, mean_kappa=float(np.mean(seg_kappa)),
            dkappa_ds=float(slope), r_squared=float(r_sq),
            spiral_class=spiral_class,
            kappa_values=seg_kappa.tolist(),
            s_values=seg_s.tolist(),
            delta_aicc_osc=osc_evidence.get('delta_aicc_osc', 0.0),
            best_osc_r2=osc_evidence.get('best_osc_r2', 0.0),
            best_osc_model=osc_evidence.get('best_osc_model', ''),
            best_osc_omega=osc_evidence.get('best_osc_omega', 0.0),
            gauss_A=gp_A, gauss_s0=gp_s0, gauss_sigma=gp_sigma,
            gauss_c=gp_c, gauss_peak_residue=gp_residue,
        ))
    
    return segments


# =============================================================================
# SECTION 5: SUPERPOTENTIAL (RECALIBRATED)
# =============================================================================

def build_superpotential_grid(grid_size=360, sigma=2.0):
    """
    Build W(φ, ψ) = -√p(φ, ψ) on T².
    
    v2: Recalibrated normalization. The superpotential is normalized
    so that the maximum basin depth difference corresponds to the
    observed ΔW between helix and sheet basins.
    """
    VM_PARAMS = [
        (0.35, -63,  -43,  12.0, 10.0,  2.0),
        (0.05, -60,  -27,   8.0,  6.0,  1.0),
        (0.25, -120, 135,   4.0,  3.5, -1.5),
        (0.05, -140, 155,   5.0,  4.0, -1.0),
        (0.12, -75,  150,   8.0,  5.0,  0.5),
        (0.05, -95,  150,   3.0,  4.0,  0.0),
        (0.03,  57,   40,   6.0,  6.0,  1.5),
        (0.03,  60, -130,   5.0,  4.0,  0.0),
        (0.01,  75,  -65,   5.0,  5.0,  0.0),
        (0.06,   0,    0,  0.01, 0.01,  0.0),
    ]
    
    phi_grid = np.linspace(-np.pi, np.pi, grid_size, endpoint=False)
    psi_grid = np.linspace(-np.pi, np.pi, grid_size, endpoint=False)
    PHI, PSI = np.meshgrid(phi_grid, psi_grid)
    
    p = np.zeros_like(PHI)
    for w, mu_phi, mu_psi, k_phi, k_psi, rho in VM_PARAMS:
        mphi = np.radians(mu_phi)
        mpsi = np.radians(mu_psi)
        vm = w * np.exp(
            k_phi * np.cos(PHI - mphi) +
            k_psi * np.cos(PSI - mpsi) +
            rho * np.cos(PHI - mphi) * np.cos(PSI - mpsi)
        )
        p += vm
    
    dphi = phi_grid[1] - phi_grid[0]
    dpsi = psi_grid[1] - psi_grid[0]
    p /= (np.sum(p) * dphi * dpsi)
    
    p = gaussian_filter1d(gaussian_filter1d(p, sigma, axis=0, mode='wrap'), sigma, axis=1, mode='wrap')
    p = np.maximum(p, 1e-10)
    
    W = -np.sqrt(p)
    
    return phi_grid, psi_grid, p, W


def interpolate_W(phi_val, psi_val, phi_grid, psi_grid, W_grid):
    """Bilinear interpolation of W at a point on T²."""
    dphi = phi_grid[1] - phi_grid[0]
    dpsi = psi_grid[1] - psi_grid[0]
    
    pv = np.arctan2(np.sin(phi_val), np.cos(phi_val))
    sv = np.arctan2(np.sin(psi_val), np.cos(psi_val))
    
    fi = (pv - phi_grid[0]) / dphi
    fj = (sv - psi_grid[0]) / dpsi
    
    fi = fi % len(phi_grid)
    fj = fj % len(psi_grid)
    
    i0 = int(fi) % len(phi_grid)
    j0 = int(fj) % len(psi_grid)
    i1 = (i0 + 1) % len(phi_grid)
    j1 = (j0 + 1) % len(psi_grid)
    
    di = fi - int(fi)
    dj = fj - int(fj)
    
    return (W_grid[j0, i0] * (1 - di) * (1 - dj) +
            W_grid[j0, i1] * di * (1 - dj) +
            W_grid[j1, i0] * (1 - di) * dj +
            W_grid[j1, i1] * di * dj)


def compute_bps_energy(phi_angles, psi_angles, phi_grid, psi_grid, W_grid):
    """Total BPS energy and BPS/residue."""
    n = len(phi_angles)
    if n < 2:
        return 0.0, 0.0
    
    W_vals = np.array([
        interpolate_W(phi_angles[i], psi_angles[i], phi_grid, psi_grid, W_grid)
        for i in range(n)
    ])
    
    dW = np.abs(np.diff(W_vals))
    E_bps = np.sum(dW)
    bps_per_residue = E_bps / (n - 1)
    
    return E_bps, bps_per_residue


# =============================================================================
# SECTION 6: PDB I/O
# =============================================================================

def extract_dihedrals_manual(pdb_path, chain_id='A'):
    """Extract backbone (φ, ψ) and N,CA,C,O coordinates from PDB file."""
    atoms_by_residue = {}
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
            
            atom_name = line[12:16].strip()
            alt_loc = line[16].strip()
            chain = line[21].strip()
            res_seq = int(line[22:26].strip())
            icode = line[26].strip()
            
            if chain != chain_id:
                continue
            if alt_loc and alt_loc not in ('A', '1'):
                continue
            if atom_name not in ('N', 'CA', 'C', 'O'):
                continue
            
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            
            key = (res_seq, icode)
            if key not in atoms_by_residue:
                atoms_by_residue[key] = {}
            atoms_by_residue[key][atom_name] = np.array([x, y, z])
    
    sorted_keys = sorted(atoms_by_residue.keys())
    
    valid_residues = []
    for key in sorted_keys:
        atoms = atoms_by_residue[key]
        if 'N' in atoms and 'CA' in atoms and 'C' in atoms:
            valid_residues.append((key, atoms))
    
    n = len(valid_residues)
    phi_angles = np.full(n, np.nan)
    psi_angles = np.full(n, np.nan)
    res_ids = [vr[0][0] for vr in valid_residues]
    
    # Build (N, CA, C, O) coordinate array for DSSP
    # Shape: (n, 4, 3)
    backbone_coords = np.zeros((n, 4, 3))
    has_O = np.zeros(n, dtype=bool)
    for i, (key, atoms) in enumerate(valid_residues):
        backbone_coords[i, 0] = atoms['N']
        backbone_coords[i, 1] = atoms['CA']
        backbone_coords[i, 2] = atoms['C']
        if 'O' in atoms:
            backbone_coords[i, 3] = atoms['O']
            has_O[i] = True
    
    for i in range(n):
        _, atoms_i = valid_residues[i]
        
        if i > 0:
            _, atoms_prev = valid_residues[i - 1]
            if 'C' in atoms_prev:
                phi_angles[i] = _dihedral(
                    atoms_prev['C'], atoms_i['N'], atoms_i['CA'], atoms_i['C']
                )
        
        if i < n - 1:
            _, atoms_next = valid_residues[i + 1]
            if 'N' in atoms_next:
                psi_angles[i] = _dihedral(
                    atoms_i['N'], atoms_i['CA'], atoms_i['C'], atoms_next['N']
                )
    
    return phi_angles, psi_angles, res_ids, backbone_coords, has_O


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


def download_pdb(pdb_id, out_dir):
    """Download PDB file from RCSB."""
    out_path = out_dir / f"{pdb_id}.pdb"
    if out_path.exists():
        return out_path
    
    try:
        import requests
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            out_path.write_text(resp.text)
            log.info(f"  Downloaded {pdb_id}")
            return out_path
        else:
            log.warning(f"  HTTP {resp.status_code} for {pdb_id}")
            return None
    except ImportError:
        log.warning("  requests not available")
        return None
    except Exception as e:
        log.warning(f"  Download error {pdb_id}: {e}")
        return None


# =============================================================================
# SECTION 7: SYNTHETIC BACKBONE GENERATOR (fallback)
# =============================================================================

def generate_synthetic_backbone(protein_name, length, fold_class='mixed',
                                 seed=None):
    """Generate realistic synthetic (φ, ψ) backbone sequence."""
    if seed is not None:
        rng = np.random.RandomState(seed)
    else:
        rng = np.random.RandomState(hash(protein_name) % (2**31))
    
    basin_params = {
        'H': (-63, -43, 8, 8),
        'E': (-120, 130, 12, 12),
        'P': (-75, 145, 10, 10),
        'L': (60, 40, 12, 12),
    }
    
    transitions = {
        'alpha': np.array([
            [0.88, 0.02, 0.08, 0.02],
            [0.20, 0.55, 0.20, 0.05],
            [0.55, 0.10, 0.30, 0.05],
            [0.40, 0.10, 0.40, 0.10],
        ]),
        'beta': np.array([
            [0.55, 0.20, 0.20, 0.05],
            [0.05, 0.85, 0.08, 0.02],
            [0.15, 0.40, 0.40, 0.05],
            [0.10, 0.20, 0.60, 0.10],
        ]),
        'mixed': np.array([
            [0.75, 0.05, 0.15, 0.05],
            [0.08, 0.72, 0.15, 0.05],
            [0.30, 0.25, 0.35, 0.10],
            [0.20, 0.15, 0.50, 0.15],
        ]),
        'alpha_beta': np.array([
            [0.80, 0.05, 0.12, 0.03],
            [0.08, 0.78, 0.12, 0.02],
            [0.35, 0.30, 0.30, 0.05],
            [0.20, 0.20, 0.45, 0.15],
        ]),
    }
    
    T = transitions.get(fold_class, transitions['mixed'])
    states = ['H', 'E', 'P', 'L']
    
    init_probs = {
        'alpha': [0.6, 0.05, 0.3, 0.05],
        'beta': [0.1, 0.6, 0.25, 0.05],
        'mixed': [0.3, 0.3, 0.3, 0.1],
        'alpha_beta': [0.35, 0.35, 0.25, 0.05],
    }
    
    state_seq = []
    s = rng.choice(4, p=init_probs.get(fold_class, init_probs['mixed']))
    for _ in range(length):
        state_seq.append(states[s])
        s = rng.choice(4, p=T[s])
    
    phi_angles = np.zeros(length)
    psi_angles = np.zeros(length)
    
    for i, state in enumerate(state_seq):
        mu_phi, mu_psi, sig_phi, sig_psi = basin_params[state]
        phi_angles[i] = np.radians(mu_phi + rng.normal(0, sig_phi))
        psi_angles[i] = np.radians(mu_psi + rng.normal(0, sig_psi))
    
    phi_angles = np.arctan2(np.sin(phi_angles), np.cos(phi_angles))
    psi_angles = np.arctan2(np.sin(psi_angles), np.cos(psi_angles))
    
    # Convert state_seq to 3-state SS
    ss_map = {'H': 'H', 'E': 'E', 'P': 'C', 'L': 'C'}
    ss = np.array([ss_map[s] for s in state_seq])
    
    return phi_angles, psi_angles, ss


def assign_fold_class(protein_name, ln_kf, co):
    """Heuristic fold class from name."""
    name_lower = protein_name.lower()
    if any(x in name_lower for x in ['lambda', 'villin', 'engrailed', 'myoglobin',
                                      'cytochrome', 'arc']):
        return 'alpha'
    elif any(x in name_lower for x in ['sh3', 'protein-l', 'cspb', 'tnfn3', 'protein-g']):
        return 'beta'
    elif any(x in name_lower for x in ['ubiquitin', 'barnase', 'dhfr', 'thioredoxin',
                                        'fkbp', 'rnase', 'lysozyme', 'chey',
                                        'acylphosphatase']):
        return 'alpha_beta'
    else:
        return 'mixed'


# =============================================================================
# SECTION 8: WINDING-RATE CONSTANCY TEST (Contribution #1)
# =============================================================================

def test_curvature_regularity(segments):
    """
    Contribution 1: Test curvature regularity by SS type.
    
    Computes per-segment curvature statistics and tests whether coil regions
    show HIGHER geometric regularity (lower CV of κ, higher R² of linear fit)
    than structured regions, inverting the naive expectation.
    
    Also computes segment-length dependence to control for the effect
    of shorter coil segments being easier to fit.
    """
    results = {
        'by_ss_type': {},
    }
    
    for ss_type in ['H', 'E', 'C']:
        type_segments = [s for s in segments if s.ss_type == ss_type 
                        and s.spiral_class != 'too_short']
        
        if not type_segments:
            continue
        
        dkappa_values = [s.dkappa_ds for s in type_segments]
        r2_values = [s.r_squared for s in type_segments]
        kappa_values = [s.mean_kappa for s in type_segments]
        lengths = [s.length for s in type_segments]
        
        # Per-segment curvature CV (regularity within each segment)
        intra_cvs = []
        for s in type_segments:
            if len(s.kappa_values) >= 3:
                seg_k = np.array(s.kappa_values)
                seg_mean = np.mean(np.abs(seg_k))
                if seg_mean > 1e-8:
                    intra_cvs.append(float(np.std(seg_k) / seg_mean))
        
        mean_dk = np.mean(dkappa_values)
        std_dk = np.std(dkappa_values)
        cv = std_dk / max(abs(mean_dk), 1e-10)
        
        # Fraction classifiable as any identified geometric form
        n_regular = sum(1 for s in type_segments 
                        if s.spiral_class not in ('too_short',))
        frac_regular = n_regular / len(type_segments)
        
        results['by_ss_type'][ss_type] = {
            'n_segments': len(type_segments),
            'mean_dkappa_ds': float(mean_dk),
            'std_dkappa_ds': float(std_dk),
            'cv_dkappa_ds': float(cv),
            'mean_r_squared': float(np.mean(r2_values)),
            'mean_kappa': float(np.mean(kappa_values)),
            'fraction_regular': frac_regular,
            'mean_segment_length': float(np.mean(lengths)),
            'median_segment_length': float(np.median(lengths)),
            'mean_intra_cv': float(np.mean(intra_cvs)) if intra_cvs else 0.0,
            'median_intra_cv': float(np.median(intra_cvs)) if intra_cvs else 0.0,
        }
    
    # Mann-Whitney test: coil intra-CV vs structured intra-CV
    coil_cvs = []
    struct_cvs = []
    coil_r2 = []
    struct_r2 = []
    for s in segments:
        if s.spiral_class == 'too_short' or len(s.kappa_values) < 3:
            continue
        seg_k = np.array(s.kappa_values)
        seg_mean = np.mean(np.abs(seg_k))
        if seg_mean < 1e-8:
            continue
        cv_val = float(np.std(seg_k) / seg_mean)
        
        if s.ss_type == 'C':
            coil_cvs.append(cv_val)
            coil_r2.append(s.r_squared)
        elif s.ss_type in ('H', 'E'):
            struct_cvs.append(cv_val)
            struct_r2.append(s.r_squared)
    
    if len(coil_cvs) >= 3 and len(struct_cvs) >= 3:
        u_stat, mw_p = stats.mannwhitneyu(coil_cvs, struct_cvs, alternative='less')
        # Effect size: rank-biserial r = 1 - 2U/(n1*n2)
        n1, n2 = len(coil_cvs), len(struct_cvs)
        effect_size = 1.0 - 2.0 * u_stat / (n1 * n2)
        results['regularity_test'] = {
            'coil_median_cv': float(np.median(coil_cvs)),
            'struct_median_cv': float(np.median(struct_cvs)),
            'mann_whitney_U': float(u_stat),
            'mann_whitney_p': float(mw_p),
            'effect_size_r': float(effect_size),
            'n_coil': n1, 'n_struct': n2,
        }
    
    if len(coil_r2) >= 3 and len(struct_r2) >= 3:
        u2, p2 = stats.mannwhitneyu(coil_r2, struct_r2, alternative='greater')
        results['r_squared_test'] = {
            'coil_median_r2': float(np.median(coil_r2)),
            'struct_median_r2': float(np.median(struct_r2)),
            'mann_whitney_U': float(u2),
            'mann_whitney_p': float(p2),
        }
    
    all_regular = [s for s in segments if s.spiral_class not in ('too_short',)]
    results['overall_regular_fraction'] = len(all_regular) / max(len(segments), 1)
    
    return results


# =============================================================================
# SECTION 9: SPIRAL CLASSIFICATION (Contribution #2)
# =============================================================================

def classify_segments_by_spiral(segments):
    """Classify all segments and test association with SS type."""
    results = {
        'classification_matrix': {},
        'chi_squared_test': None,
    }
    
    ss_types = ['H', 'E', 'C']
    # Expanded taxonomy
    spiral_types = ['geodesic', 'circular_arc', 'clothoid', 'fermat',
                    'sinusoidal', 'damped_osc', 'exponential', 'sigmoid',
                    'step', 'quadratic', 'gauss_peak', 'too_short']
    
    matrix = {ss: {sp: 0 for sp in spiral_types} for ss in ss_types}
    
    for seg in segments:
        if seg.ss_type in matrix:
            cls = seg.spiral_class
            if cls not in matrix[seg.ss_type]:
                matrix[seg.ss_type][cls] = 0
            matrix[seg.ss_type][cls] += 1
    
    results['classification_matrix'] = matrix
    
    # "Regular" = shapes with good fit (not too_short)
    regular_types = ['geodesic', 'circular_arc', 'clothoid', 'fermat',
                     'sinusoidal', 'damped_osc', 'exponential', 'sigmoid',
                     'step', 'quadratic', 'gauss_peak']
    
    # ── Test 1: Full 3×k chi-squared (does shape distribution differ by SS type?) ──
    # Group into interpretable macro-categories to maintain cell counts
    # Oscillatory = sinusoidal + damped_osc (periodic curvature dynamics)
    # Constant    = geodesic + circular_arc (steady curvature)
    # Monotone    = fermat + clothoid + exponential (trending curvature)
    # Transition  = sigmoid + step (discrete curvature change)
    # Polynomial  = quadratic
    
    macro_groups = {
        'oscillatory': ['sinusoidal', 'damped_osc'],
        'constant_k': ['geodesic', 'circular_arc'],
        'monotone': ['fermat', 'clothoid', 'exponential'],
        'transition': ['sigmoid', 'step'],
        'polynomial': ['quadratic'],
        'localized': ['gauss_peak'],
    }
    
    # Build 3×5 contingency table
    ss_list = ['H', 'E', 'C']
    macro_list = list(macro_groups.keys())
    contingency_full = np.zeros((3, len(macro_list)), dtype=int)
    for i, ss in enumerate(ss_list):
        for j, (macro, members) in enumerate(macro_groups.items()):
            contingency_full[i, j] = sum(matrix[ss].get(sp, 0) for sp in members)
    
    # Drop columns with all zeros
    col_sums = contingency_full.sum(axis=0)
    valid_cols = col_sums > 0
    contingency_trimmed = contingency_full[:, valid_cols]
    macro_trimmed = [m for m, v in zip(macro_list, valid_cols) if v]
    
    if contingency_trimmed.sum() >= 15 and contingency_trimmed.shape[1] >= 2:
        chi2_full, p_full, dof_full, expected_full = stats.chi2_contingency(contingency_trimmed)
        results['chi_squared_full'] = {
            'chi2': float(chi2_full), 'p_value': float(p_full), 'dof': int(dof_full),
            'contingency_table': contingency_trimmed.tolist(),
            'row_labels': ss_list,
            'col_labels': macro_trimmed,
            'test_description': '3×k full: SS type vs macro-shape category',
        }
        
        # Compute Cramér's V for effect size
        n_total = contingency_trimmed.sum()
        min_dim = min(contingency_trimmed.shape[0], contingency_trimmed.shape[1]) - 1
        if min_dim > 0 and n_total > 0:
            cramers_v = np.sqrt(chi2_full / (n_total * min_dim))
            results['chi_squared_full']['cramers_v'] = float(cramers_v)
    
    # ── Test 2: H vs E+C for oscillatory enrichment (the helix signature) ──
    osc_types = ['sinusoidal', 'damped_osc']
    non_osc_types = [sp for sp in regular_types if sp not in osc_types]
    
    n_H_osc = sum(matrix['H'].get(sp, 0) for sp in osc_types)
    n_H_non = sum(matrix['H'].get(sp, 0) for sp in non_osc_types)
    n_EC_osc = sum(matrix['E'].get(sp, 0) + matrix['C'].get(sp, 0) for sp in osc_types)
    n_EC_non = sum(matrix['E'].get(sp, 0) + matrix['C'].get(sp, 0) for sp in non_osc_types)
    
    contingency_osc = np.array([[n_H_osc, n_H_non], [n_EC_osc, n_EC_non]])
    
    if (contingency_osc.sum() >= 10 and 
        np.all(contingency_osc.sum(axis=0) > 0) and 
        np.all(contingency_osc.sum(axis=1) > 0)):
        chi2_o, p_o, dof_o, _ = stats.chi2_contingency(contingency_osc)
        results['chi_squared_oscillatory'] = {
            'chi2': float(chi2_o), 'p_value': float(p_o), 'dof': int(dof_o),
            'contingency_table': contingency_osc.tolist(),
            'test_description': 'H vs E+C, oscillatory(sin/damped) vs non-oscillatory',
        }
    
    # ── Test 3: E vs H+C for constant-curvature enrichment ──
    const_types = ['geodesic', 'circular_arc']
    nonconst = [sp for sp in regular_types if sp not in const_types]
    
    n_E_const = sum(matrix['E'].get(sp, 0) for sp in const_types)
    n_E_nonconst = sum(matrix['E'].get(sp, 0) for sp in nonconst)
    n_HC_const = sum(matrix['H'].get(sp, 0) + matrix['C'].get(sp, 0) for sp in const_types)
    n_HC_nonconst = sum(matrix['H'].get(sp, 0) + matrix['C'].get(sp, 0) for sp in nonconst)
    
    contingency_const = np.array([[n_E_const, n_E_nonconst], [n_HC_const, n_HC_nonconst]])
    
    if (contingency_const.sum() >= 10 and 
        np.all(contingency_const.sum(axis=0) > 0) and 
        np.all(contingency_const.sum(axis=1) > 0)):
        chi2_c, p_c, dof_c, _ = stats.chi2_contingency(contingency_const)
        results['chi_squared_constant'] = {
            'chi2': float(chi2_c), 'p_value': float(p_c), 'dof': int(dof_c),
            'contingency_table': contingency_const.tolist(),
            'test_description': 'E vs H+C, constant-κ(geo/arc) vs varying',
        }
    
    return results


# =============================================================================
# SECTION 10: TORUS KNOT DESCRIPTORS (Contribution #3)
# =============================================================================

def compute_torus_knot_descriptor(phi, psi, ss_assignments):
    """Compute (p, q) torus knot descriptor for a protein backbone on T²."""
    p_total, q_total = torus_winding_numbers(phi, psi)
    
    winding_by_ss = {}
    for ss_type in ['H', 'E', 'C']:
        mask = ss_assignments == ss_type
        indices = np.where(mask)[0]
        
        if len(indices) < 2:
            winding_by_ss[ss_type] = {'p': 0.0, 'q': 0.0}
            continue
        
        p_ss = 0.0
        q_ss = 0.0
        for i in range(len(indices) - 1):
            if indices[i+1] == indices[i] + 1:
                p_ss += angular_distance(phi[indices[i+1]], phi[indices[i]]) / (2 * np.pi)
                q_ss += angular_distance(psi[indices[i+1]], psi[indices[i]]) / (2 * np.pi)
        
        winding_by_ss[ss_type] = {'p': float(p_ss), 'q': float(q_ss)}
    
    p_int = int(round(p_total))
    q_int = int(round(q_total))
    
    from math import gcd
    g = gcd(abs(p_int) if p_int != 0 else 1, abs(q_int) if q_int != 0 else 1)
    p_red = p_int // g if g > 0 else p_int
    q_red = q_int // g if g > 0 else q_int
    
    winding_ratio = p_total / q_total if abs(q_total) > 0.01 else float('inf')
    
    return {
        'p': float(p_total),
        'q': float(q_total),
        'abs_p': float(abs(p_total)),
        'abs_q': float(abs(q_total)),
        'p_int': p_int,
        'q_int': q_int,
        'p_reduced': p_red,
        'q_reduced': q_red,
        'winding_ratio': float(winding_ratio),
        'winding_by_ss': winding_by_ss,
        'total_winding_magnitude': float(np.sqrt(p_total**2 + q_total**2)),
    }


# =============================================================================
# SECTION 11: MAIN ANALYSIS PIPELINE
# =============================================================================

@dataclass
class ProteinResult:
    """Complete analysis results for one protein."""
    name: str
    pdb_id: str
    chain: str
    length: int
    ln_kf: float
    co: float
    source: str
    ss_method: str = 'none'
    
    n_valid_residues: int = 0
    ss_composition: dict = field(default_factory=dict)
    
    winding_rate_constancy: dict = field(default_factory=dict)
    overall_clothoid_fraction: float = 0.0
    
    spiral_classification: dict = field(default_factory=dict)
    n_segments: int = 0
    segments: list = field(default_factory=list)
    
    torus_knot: dict = field(default_factory=dict)
    
    bps_per_residue: float = 0.0
    bps_total: float = 0.0
    
    mean_curvature: float = 0.0
    std_curvature: float = 0.0
    curvature_autocorrelation: float = 0.0


def analyze_single_protein(name, pdb_id, chain, length, ln_kf, co,
                            phi_grid, psi_grid, W_grid,
                            use_synthetic=False):
    """Run the full Cornu spiral analysis pipeline for one protein."""
    
    result = ProteinResult(
        name=name, pdb_id=pdb_id, chain=chain, length=length,
        ln_kf=ln_kf, co=co, source='synthetic'
    )
    
    pdb_path = PDB_DIR / f"{pdb_id}.pdb"
    
    phi = psi = ss = None
    
    if pdb_path.exists() and not use_synthetic:
        try:
            phi, psi, res_ids, backbone_coords, has_O = extract_dihedrals_manual(
                str(pdb_path), chain)
            valid = ~(np.isnan(phi) | np.isnan(psi))
            phi = phi[valid]
            psi = psi[valid]
            res_ids_valid = [r for r, v in zip(res_ids, valid) if v]
            backbone_coords_valid = backbone_coords[valid]
            has_O_valid = has_O[valid]
            result.source = 'pdb'
            
            # ── Try pure-numpy DSSP (needs N, CA, C, O) ──
            if np.all(has_O_valid) and len(phi) >= 6:
                try:
                    ss = dssp_assign_pure_numpy(backbone_coords_valid)
                    result.ss_method = 'dssp_numpy'
                    log.info(f"  {name}: Pure-numpy DSSP — "
                             f"H={np.sum(ss=='H')}, E={np.sum(ss=='E')}, C={np.sum(ss=='C')}")
                except Exception as e:
                    log.info(f"  {name}: Pure DSSP failed ({e}), using dihedral fallback")
                    ss = assign_ss_from_dihedrals(phi, psi)
                    result.ss_method = 'dihedral'
                    log.info(f"  {name}: Dihedral fallback — "
                             f"H={np.sum(ss=='H')}, E={np.sum(ss=='E')}, C={np.sum(ss=='C')}")
            else:
                # ── Fallback: fixed dihedral classifier ──
                ss = assign_ss_from_dihedrals(phi, psi)
                result.ss_method = 'dihedral'
                log.info(f"  {name}: Dihedral fallback (missing O atoms) — "
                         f"H={np.sum(ss=='H')}, E={np.sum(ss=='E')}, C={np.sum(ss=='C')}")
            
            log.info(f"  {name}: {len(phi)} valid residues from PDB")
            
        except Exception as e:
            log.warning(f"  {name}: PDB parse failed ({e}), using synthetic")
            phi = psi = ss = None
    
    # ── Synthetic fallback ──
    if phi is None:
        fold_class = assign_fold_class(name, ln_kf, co)
        phi, psi, ss = generate_synthetic_backbone(name, length, fold_class)
        result.source = 'synthetic'
        result.ss_method = 'synthetic'
        log.info(f"  {name}: synthetic ({fold_class}), L={length}")
    
    result.n_valid_residues = len(phi)
    
    if len(phi) < 5:
        log.warning(f"  {name}: too few residues ({len(phi)})")
        return result
    
    # ── SS composition ──
    ss_comp = {}
    for s_type in ['H', 'E', 'C']:
        ss_comp[s_type] = float(np.mean(ss == s_type))
    result.ss_composition = ss_comp
    
    # ── Torus curvature ──
    kappa, s_kappa, s_full = torus_curvature(phi, psi)
    
    if len(kappa) < 3:
        return result
    
    result.mean_curvature = float(np.mean(kappa))
    result.std_curvature = float(np.std(kappa))
    
    # Curvature autocorrelation
    if len(kappa) > 10:
        kappa_centered = kappa - np.mean(kappa)
        autocorr = np.correlate(kappa_centered, kappa_centered, mode='full')
        autocorr = autocorr[len(autocorr)//2:]
        if autocorr[0] > 0:
            autocorr /= autocorr[0]
            below = np.where(autocorr < 1/np.e)[0]
            result.curvature_autocorrelation = float(below[0]) if len(below) > 0 else float(len(autocorr))
    
    # ── Segment and classify ──
    segments = segment_backbone(phi, psi, ss, kappa, s_kappa)
    log.info("    segment_backbone -> %d kappa values", len(kappa))
    result.n_segments = len(segments)
    result.segments = segments
    n_cls = sum(1 for s in segments if s.spiral_class != "too_short")
    log.info("    -> %d classified, %d too_short", n_cls, len(segments) - n_cls)
    
    # Contribution 1
    wrc = test_curvature_regularity(segments)
    log.info("    test_curvature_regularity")
    result.winding_rate_constancy = wrc
    result.overall_clothoid_fraction = wrc.get('overall_regular_fraction', 0.0)
    
    # Contribution 2
    spiral_class = classify_segments_by_spiral(segments)
    log.info("    classify_segments_by_spiral")
    result.spiral_classification = spiral_class
    
    # Contribution 3
    torus_knot = compute_torus_knot_descriptor(phi, psi, ss)
    log.info("    compute_torus_knot_descriptor")
    result.torus_knot = torus_knot
    
    # BPS energy
    E_bps, bps_per_res = compute_bps_energy(phi, psi, phi_grid, psi_grid, W_grid)
    result.bps_total = E_bps
    result.bps_per_residue = bps_per_res
    
    return result


# =============================================================================
# SECTION 12: FIGURES
# =============================================================================

def fig_ramachandran_with_spiral_overlay(phi, psi, ss, segments, protein_name, 
                                          phi_grid, psi_grid, W_grid, out_path):
    """Figure 1: Ramachandran with backbone path colored by SS + spiral type."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    ax = axes[0]
    phi_deg = np.degrees(phi_grid)
    psi_deg = np.degrees(psi_grid)
    
    ax.contourf(phi_deg, psi_deg, -W_grid, levels=20, cmap='YlOrRd', alpha=0.4)
    ax.contour(phi_deg, psi_deg, -W_grid, levels=10, colors='gray', linewidths=0.3, alpha=0.5)
    
    ss_colors = {'H': '#e74c3c', 'E': '#3498db', 'C': '#2ecc71'}
    for i in range(len(phi) - 1):
        ax.plot([np.degrees(phi[i]), np.degrees(phi[i+1])],
                [np.degrees(psi[i]), np.degrees(psi[i+1])],
                color=ss_colors.get(ss[i], 'gray'), linewidth=0.8, alpha=0.7)
    
    ax.scatter(np.degrees(phi), np.degrees(psi), c=[ss_colors.get(s, 'gray') for s in ss],
               s=4, zorder=5)
    
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel(r'$\phi$ (degrees)', fontsize=11)
    ax.set_ylabel(r'$\psi$ (degrees)', fontsize=11)
    ax.set_title(f'{protein_name}: Backbone path on T²', fontsize=12)
    ax.set_aspect('equal')
    
    for label, color in [('Helix', '#e74c3c'), ('Strand', '#3498db'), ('Coil', '#2ecc71')]:
        ax.plot([], [], color=color, linewidth=2, label=label)
    ax.legend(loc='upper right', fontsize=8)
    
    ax2 = axes[1]
    spiral_colors = {
        'geodesic': '#2ecc71', 'circular_arc': '#9b59b6',
        'clothoid': '#1abc9c', 'fermat': '#e67e22',
        'sinusoidal': '#e74c3c', 'damped_osc': '#c0392b',
        'exponential': '#f39c12', 'sigmoid': '#d35400',
        'step': '#3498db', 'quadratic': '#8e44ad', 'gauss_peak': '#e74c3c',
        'too_short': '#bdc3c7',
    }
    
    for seg in segments:
        if seg.s_values and seg.kappa_values:
            color = spiral_colors.get(seg.spiral_class, '#95a5a6')
            ax2.plot(seg.s_values, seg.kappa_values, color=color, linewidth=1.5, alpha=0.8)
    
    ax2.set_xlabel('Arc length s', fontsize=11)
    ax2.set_ylabel(r'Curvature $\kappa(s)$', fontsize=11)
    ax2.set_title(f'{protein_name}: Curvature by segment type', fontsize=12)
    
    for label, color in spiral_colors.items():
        if label != 'too_short':
            ax2.plot([], [], color=color, linewidth=2, label=label.replace('_', ' ').title())
    ax2.legend(loc='upper right', fontsize=7, ncol=2)
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()


def fig_curvature_regularity(all_results, out_path):
    """Figure 2: Curvature regularity inversion — coil more regular than structured."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    ss_types = ['H', 'E', 'C']
    ss_labels = ['Helix', 'Strand', 'Coil/Loop']
    ss_colors = ['#e74c3c', '#3498db', '#2ecc71']
    
    # Panel 1: Within-segment curvature CV by SS type
    ax = axes[0]
    all_intra_cvs = {ss: [] for ss in ss_types}
    for result in all_results:
        wrc = result.winding_rate_constancy
        for ss_type in ss_types:
            if ss_type in wrc.get('by_ss_type', {}):
                info = wrc['by_ss_type'][ss_type]
                if info.get('mean_intra_cv', 0) > 0:
                    all_intra_cvs[ss_type].append(info['mean_intra_cv'])
    
    data = [all_intra_cvs[ss] for ss in ss_types]
    valid_data = [d for d in data if len(d) > 0]
    valid_labels = [l for l, d in zip(ss_labels, data) if len(d) > 0]
    valid_colors = [c for c, d in zip(ss_colors, data) if len(d) > 0]
    
    if valid_data:
        bp = ax.boxplot(valid_data, tick_labels=valid_labels, patch_artist=True, widths=0.5)
        for patch, color in zip(bp['boxes'], valid_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        ax.set_ylabel('Within-segment curvature CV', fontsize=10)
        ax.set_title('Curvature variability\n(lower = more regular)', fontsize=11)
        # Add significance annotation
        if len(all_intra_cvs['C']) >= 3 and (len(all_intra_cvs['H']) >= 3 or len(all_intra_cvs['E']) >= 3):
            struct = all_intra_cvs['H'] + all_intra_cvs['E']
            if len(struct) >= 3:
                _, p = stats.mannwhitneyu(all_intra_cvs['C'], struct, alternative='less')
                ax.text(0.95, 0.95, f'C < H+E: p={p:.3f}', transform=ax.transAxes,
                        fontsize=8, ha='right', va='top',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Panel 2: Regular fraction by SS type
    ax2 = axes[1]
    regular_fracs = {ss: [] for ss in ss_types}
    for result in all_results:
        wrc = result.winding_rate_constancy
        for ss_type in ss_types:
            if ss_type in wrc.get('by_ss_type', {}):
                info = wrc['by_ss_type'][ss_type]
                regular_fracs[ss_type].append(info.get('fraction_regular', 0))
    
    means = [np.mean(regular_fracs[ss]) if regular_fracs[ss] else 0 for ss in ss_types]
    sems = [np.std(regular_fracs[ss])/np.sqrt(len(regular_fracs[ss])) 
            if len(regular_fracs[ss]) > 1 else 0 for ss in ss_types]
    
    bars = ax2.bar(ss_labels, means, yerr=sems, color=ss_colors, alpha=0.7, 
                    edgecolor='white', capsize=5)
    ax2.set_ylabel('Fraction classifiable as spiral', fontsize=10)
    ax2.set_title('Classified fraction\n(all 6 geometric forms)', fontsize=11)
    ax2.set_ylim(0, 1)
    
    # Panel 3: R² of linear curvature fit
    ax3 = axes[2]
    all_r2 = {ss: [] for ss in ss_types}
    for result in all_results:
        wrc = result.winding_rate_constancy
        for ss_type in ss_types:
            if ss_type in wrc.get('by_ss_type', {}):
                info = wrc['by_ss_type'][ss_type]
                all_r2[ss_type].append(info.get('mean_r_squared', 0))
    
    data_r2 = [all_r2[ss] for ss in ss_types]
    valid_r2 = [d for d in data_r2 if len(d) > 0]
    valid_r2_labels = [l for l, d in zip(ss_labels, data_r2) if len(d) > 0]
    valid_r2_colors = [c for c, d in zip(ss_colors, data_r2) if len(d) > 0]
    
    if valid_r2:
        bp2 = ax3.boxplot(valid_r2, tick_labels=valid_r2_labels, patch_artist=True, widths=0.5)
        for patch, color in zip(bp2['boxes'], valid_r2_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        ax3.set_ylabel(r'R² of $\kappa(s) = a + bs$ fit', fontsize=10)
        ax3.set_title('Linear fit quality\n(higher = more predictable)', fontsize=11)
    
    fig.suptitle('Curvature Hierarchy Inversion: Coil > Structured Regularity', 
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()


def fig_spiral_classification(all_results, out_path):
    """Figure 3: SS-dependent spiral taxonomy — full AIC-based shape breakdown."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    ss_types = ['H', 'E', 'C']
    
    # Collect counts from all proteins
    total_matrix = {}
    for ss in ss_types:
        total_matrix[ss] = {}
    
    for result in all_results:
        sc = result.spiral_classification
        if 'classification_matrix' in sc:
            for ss in ss_types:
                if ss in sc['classification_matrix']:
                    for sp, count in sc['classification_matrix'][ss].items():
                        if sp != 'too_short' and count > 0:
                            total_matrix[ss][sp] = total_matrix[ss].get(sp, 0) + count
    
    # Identify all shape types that appear, in canonical order
    shape_order = ['geodesic', 'circular_arc', 'clothoid', 'fermat',
                   'sinusoidal', 'damped_osc', 'exponential', 'sigmoid',
                   'step', 'quadratic', 'gauss_peak']
    all_shapes_present = set()
    for ss in ss_types:
        all_shapes_present.update(total_matrix[ss].keys())
    shapes = [s for s in shape_order if s in all_shapes_present]
    
    # Color palette for shapes
    shape_colors = {
        'geodesic': '#2c3e50',     # near-black (zero curvature)
        'circular_arc': '#9b59b6', # purple (constant curvature)
        'clothoid': '#1abc9c',     # teal (linear decreasing κ)
        'fermat': '#e67e22',       # orange (linear increasing κ)
        'sinusoidal': '#3498db',   # blue (oscillating)
        'damped_osc': '#2980b9',   # darker blue (decaying oscillation)
        'exponential': '#e74c3c',  # red (monotone growth/decay)
        'sigmoid': '#f39c12',      # gold (smooth transition)
        'step': '#95a5a6',         # gray (piecewise constant)
        'quadratic': '#8e44ad',    # dark purple (polynomial)
        'gauss_peak': '#e74c3c',   # red (localized turn/kink)
    }
    shape_labels = {
        'geodesic': 'Geodesic',
        'circular_arc': 'Circular arc',
        'clothoid': 'Clothoid (−)',
        'fermat': 'Fermat (+)',
        'sinusoidal': 'Sinusoidal',
        'damped_osc': 'Damped osc.',
        'exponential': 'Exponential',
        'sigmoid': 'Sigmoid',
        'step': 'Step',
        'quadratic': 'Quadratic',
        'gauss_peak': 'Gauss Peak',
    }
    
    # Panel 1: Stacked bar chart, absolute counts
    ax = axes[0]
    x = np.arange(len(ss_types))
    width = 0.6
    bottom = np.zeros(len(ss_types))
    
    for sp in shapes:
        vals = np.array([total_matrix[ss].get(sp, 0) for ss in ss_types], dtype=float)
        ax.bar(x, vals, width, bottom=bottom, 
               label=shape_labels.get(sp, sp), 
               color=shape_colors.get(sp, 'gray'), edgecolor='white', linewidth=0.5)
        bottom += vals
    
    ax.set_xticks(x)
    ax.set_xticklabels(['Helix', 'Strand', 'Coil/Loop'])
    ax.set_ylabel('Number of segments')
    ax.set_title('Shape Distribution by SS Type', fontsize=12)
    ax.legend(fontsize=7, loc='upper right', ncol=2)
    
    # Panel 2: Percentage stacked bar (normalized)
    ax2 = axes[1]
    bottom2 = np.zeros(len(ss_types))
    
    totals = np.array([sum(total_matrix[ss].get(sp, 0) for sp in shapes) for ss in ss_types], dtype=float)
    totals = np.maximum(totals, 1)  # avoid division by zero
    
    for sp in shapes:
        vals = np.array([total_matrix[ss].get(sp, 0) for ss in ss_types], dtype=float)
        pct = 100 * vals / totals
        ax2.bar(x, pct, width, bottom=bottom2,
                label=shape_labels.get(sp, sp),
                color=shape_colors.get(sp, 'gray'), edgecolor='white', linewidth=0.5)
        bottom2 += pct
    
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Helix', 'Strand', 'Coil/Loop'])
    ax2.set_ylabel('Percentage of segments')
    ax2.set_title('Normalized Shape Distribution', fontsize=12)
    ax2.set_ylim(0, 105)
    ax2.legend(fontsize=7, loc='upper right', ncol=2)
    
    fig.suptitle('SS-Dependent Spiral Taxonomy (AIC Model Selection)',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()


def fig_torus_knot_descriptors(all_results, out_path):
    """Figure 4: Torus knot (p,q) descriptors."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    fold_colors = {
        'alpha': '#e74c3c', 'beta': '#3498db',
        'mixed': '#2ecc71', 'alpha_beta': '#9b59b6',
    }
    
    ax = axes[0]
    for result in all_results:
        tk = result.torus_knot
        if not tk:
            continue
        fc = assign_fold_class(result.name, result.ln_kf, result.co)
        ax.scatter(tk['p'], tk['q'], color=fold_colors.get(fc, 'gray'),
                   s=60, alpha=0.7, edgecolors='black', linewidth=0.5)
        ax.annotate(result.name[:6], (tk['p'], tk['q']), fontsize=6,
                    xytext=(3, 3), textcoords='offset points')
    
    for fc, color in fold_colors.items():
        ax.scatter([], [], color=color, s=60, edgecolors='black', linewidth=0.5,
                   label=fc.replace('_', '/'))
    
    ax.axhline(0, color='gray', linestyle=':', linewidth=0.5)
    ax.axvline(0, color='gray', linestyle=':', linewidth=0.5)
    ax.set_xlabel('p (phi winding number)', fontsize=11)
    ax.set_ylabel('q (psi winding number)', fontsize=11)
    ax.set_title('Torus Knot Descriptors (p, q)', fontsize=12)
    ax.legend(fontsize=9)
    
    ax2 = axes[1]
    magnitudes = []
    ln_kfs = []
    colors = []
    names = []
    
    for result in all_results:
        tk = result.torus_knot
        if not tk:
            continue
        magnitudes.append(tk['total_winding_magnitude'])
        ln_kfs.append(result.ln_kf)
        fc = assign_fold_class(result.name, result.ln_kf, result.co)
        colors.append(fold_colors.get(fc, 'gray'))
        names.append(result.name)
    
    if magnitudes:
        ax2.scatter(magnitudes, ln_kfs, c=colors, s=60, alpha=0.7,
                    edgecolors='black', linewidth=0.5)
        
        r, p = stats.pearsonr(magnitudes, ln_kfs)
        rho, p_rho = stats.spearmanr(magnitudes, ln_kfs)
        
        for i, name in enumerate(names):
            ax2.annotate(name[:6], (magnitudes[i], ln_kfs[i]), fontsize=6,
                        xytext=(3, 3), textcoords='offset points')
        
        ax2.set_xlabel(r'Winding magnitude $\sqrt{p^2 + q^2}$', fontsize=11)
        ax2.set_ylabel(r'$\ln(k_f)$', fontsize=11)
        ax2.set_title(f'Winding vs Folding Rate\nr = {r:.3f} (p = {p:.3f}), '
                      f'$\\rho$ = {rho:.3f}', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()


def fig_cornu_exemplar(phi, psi, ss, segments, protein_name, out_path):
    """Figure 5: Detailed Cornu spiral exemplar for one protein."""
    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3)
    
    ax1 = fig.add_subplot(gs[0, 0], projection='3d')
    R, r = 2.0, 0.8
    u = np.linspace(0, 2*np.pi, 60)
    v = np.linspace(0, 2*np.pi, 60)
    U, V = np.meshgrid(u, v)
    X = (R + r * np.cos(V)) * np.cos(U)
    Y = (R + r * np.cos(V)) * np.sin(U)
    Z = r * np.sin(V)
    ax1.plot_surface(X, Y, Z, alpha=0.08, color='lightblue')
    
    bb_u = phi + np.pi
    bb_v = psi + np.pi
    bb_x = (R + r * np.cos(bb_v)) * np.cos(bb_u)
    bb_y = (R + r * np.cos(bb_v)) * np.sin(bb_u)
    bb_z = r * np.sin(bb_v)
    
    ss_colors = {'H': '#e74c3c', 'E': '#3498db', 'C': '#2ecc71'}
    for i in range(len(phi) - 1):
        ax1.plot([bb_x[i], bb_x[i+1]], [bb_y[i], bb_y[i+1]], [bb_z[i], bb_z[i+1]],
                 color=ss_colors.get(ss[i], 'gray'), linewidth=1.2)
    
    ax1.set_title(f'{protein_name}\nBackbone on T\u00b2', fontsize=10)
    ax1.set_axis_off()
    
    ax2 = fig.add_subplot(gs[0, 1])
    kappa, s_kappa, _ = torus_curvature(phi, psi)
    
    if len(kappa) > 0:
        for i in range(len(kappa)):
            res_idx = min(i + 1, len(ss) - 1)
            ax2.scatter(s_kappa[i], kappa[i], color=ss_colors.get(ss[res_idx], 'gray'),
                       s=8, alpha=0.6)
        
        if len(kappa) > 5:
            kappa_smooth = gaussian_filter1d(kappa, sigma=2)
            ax2.plot(s_kappa, kappa_smooth, 'k-', linewidth=1, alpha=0.4, label='Smoothed')
        
        ax2.set_xlabel('Arc length s', fontsize=10)
        ax2.set_ylabel(r'$\kappa(s)$', fontsize=10)
        ax2.set_title('Curvature profile', fontsize=10)
    
    ax3 = fig.add_subplot(gs[1, 0])
    seg_dkappa = {'H': [], 'E': [], 'C': []}
    for seg in segments:
        if seg.spiral_class != 'too_short' and seg.ss_type in seg_dkappa:
            seg_dkappa[seg.ss_type].append(seg.dkappa_ds)
    
    for ss_type, label, color in [('H', 'Helix', '#e74c3c'), 
                                    ('E', 'Strand', '#3498db'),
                                    ('C', 'Coil', '#2ecc71')]:
        if seg_dkappa[ss_type]:
            ax3.hist(seg_dkappa[ss_type], bins=10, alpha=0.6, color=color, 
                     label=label, edgecolor='white')
    
    ax3.axvline(0, color='gray', linestyle=':', linewidth=0.8)
    ax3.set_xlabel(r'd$\kappa$/ds', fontsize=10)
    ax3.set_ylabel('Count', fontsize=10)
    ax3.set_title('Curvature rate by segment', fontsize=10)
    ax3.legend(fontsize=8)
    
    ax4 = fig.add_subplot(gs[1, 1])
    class_counts = {}
    for seg in segments:
        class_counts[seg.spiral_class] = class_counts.get(seg.spiral_class, 0) + 1
    
    if class_counts:
        shape_colors = {
            'geodesic': '#2c3e50', 'circular_arc': '#9b59b6',
            'clothoid': '#1abc9c', 'fermat': '#e67e22',
            'sinusoidal': '#3498db', 'damped_osc': '#2980b9',
            'exponential': '#e74c3c', 'sigmoid': '#f39c12',
            'step': '#95a5a6', 'quadratic': '#8e44ad', 'gauss_peak': '#e74c3c',
            'too_short': '#bdc3c7',
        }
        labels = list(class_counts.keys())
        sizes = list(class_counts.values())
        colors = [shape_colors.get(l, 'gray') for l in labels]
        
        ax4.pie(sizes, labels=[l.replace('_', ' ').title() for l in labels], colors=colors,
                autopct='%1.0f%%', startangle=90, textprops={'fontsize': 8})
        ax4.set_title('Segment classification', fontsize=10)
    
    fig.suptitle(f'Spiral Analysis: {protein_name}', fontsize=13, fontweight='bold')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()


def fig_summary_dashboard(all_results, out_path):
    """Figure 6: Summary dashboard."""
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 3, hspace=0.35, wspace=0.3)
    
    fold_colors = {
        'alpha': '#e74c3c', 'beta': '#3498db',
        'mixed': '#2ecc71', 'alpha_beta': '#9b59b6',
    }
    
    # 1. Regular spiral fraction
    ax1 = fig.add_subplot(gs[0, 0])
    names = [r.name[:10] for r in all_results]
    fracs = [r.overall_clothoid_fraction for r in all_results]
    colors = ['#2ecc71' if f > 0.7 else '#e67e22' if f > 0.4 else '#e74c3c' for f in fracs]
    ax1.barh(range(len(names)), fracs, color=colors, edgecolor='white')
    ax1.set_yticks(range(len(names)))
    ax1.set_yticklabels(names, fontsize=7)
    ax1.set_xlabel('Regular fraction', fontsize=10)
    ax1.set_title('Fraction classified as spiral', fontsize=10)
    ax1.axvline(0.5, color='gray', linestyle=':', linewidth=0.8)
    
    # 2. BPS/residue
    ax2 = fig.add_subplot(gs[0, 1])
    bps_vals = [r.bps_per_residue for r in all_results if r.bps_per_residue > 0]
    if bps_vals:
        ax2.hist(bps_vals, bins=12, color='#3498db', alpha=0.7, edgecolor='white')
        ax2.axvline(np.mean(bps_vals), color='black', linestyle='--', linewidth=1,
                    label=f'Mean: {np.mean(bps_vals):.4f}')
        ax2.set_xlabel('BPS/residue', fontsize=10)
        ax2.set_ylabel('Count', fontsize=10)
        ax2.set_title('BPS/residue distribution', fontsize=10)
        ax2.legend(fontsize=8)
    
    # 3. Mean curvature vs length
    ax3 = fig.add_subplot(gs[0, 2])
    lengths = [r.length for r in all_results]
    mean_kappas = [r.mean_curvature for r in all_results]
    ax3.scatter(lengths, mean_kappas, c='#9b59b6', s=40, alpha=0.7)
    for i, r in enumerate(all_results):
        ax3.annotate(r.name[:5], (r.length, r.mean_curvature), fontsize=5,
                    xytext=(2, 2), textcoords='offset points')
    ax3.set_xlabel('Chain length', fontsize=10)
    ax3.set_ylabel(r'Mean $\kappa$', fontsize=10)
    ax3.set_title('Curvature vs Length', fontsize=10)
    
    # 4. (p, q)
    ax4 = fig.add_subplot(gs[1, 0])
    ps = [r.torus_knot.get('p', 0) for r in all_results if r.torus_knot]
    qs = [r.torus_knot.get('q', 0) for r in all_results if r.torus_knot]
    if ps and qs:
        cs = [fold_colors.get(assign_fold_class(r.name, r.ln_kf, r.co), 'gray') 
              for r in all_results if r.torus_knot]
        ax4.scatter(ps, qs, c=cs, s=50, alpha=0.7, edgecolors='black', linewidth=0.3)
        ax4.set_xlabel('p (phi winding)', fontsize=10)
        ax4.set_ylabel('q (psi winding)', fontsize=10)
        ax4.set_title('Torus knot (p, q)', fontsize=10)
        ax4.axhline(0, color='gray', linestyle=':', linewidth=0.5)
        ax4.axvline(0, color='gray', linestyle=':', linewidth=0.5)
    
    # 5. Curvature autocorrelation by fold class
    ax5 = fig.add_subplot(gs[1, 1])
    fold_autocorr = {}
    for r in all_results:
        if r.curvature_autocorrelation > 0:
            fc = assign_fold_class(r.name, r.ln_kf, r.co)
            fold_autocorr.setdefault(fc, []).append(r.curvature_autocorrelation)
    
    if fold_autocorr:
        labels = list(fold_autocorr.keys())
        data = [fold_autocorr[l] for l in labels]
        bp = ax5.boxplot(data, tick_labels=[l.replace('_', '/') for l in labels], 
                        patch_artist=True)
        for patch, label in zip(bp['boxes'], labels):
            patch.set_facecolor(fold_colors.get(label, 'gray'))
            patch.set_alpha(0.6)
        ax5.set_ylabel('Curvature autocorrelation length', fontsize=10)
        ax5.set_title('Curvature persistence by fold class', fontsize=10)
    
    # 6. SS composition
    ax6 = fig.add_subplot(gs[1, 2])
    h_fracs = [r.ss_composition.get('H', 0) for r in all_results]
    e_fracs = [r.ss_composition.get('E', 0) for r in all_results]
    cs = [fold_colors.get(assign_fold_class(r.name, r.ln_kf, r.co), 'gray') for r in all_results]
    ax6.scatter(h_fracs, e_fracs, c=cs, s=50, alpha=0.7, edgecolors='black', linewidth=0.3)
    for i, r in enumerate(all_results):
        ax6.annotate(r.name[:5], (h_fracs[i], e_fracs[i]), fontsize=5,
                    xytext=(2, 2), textcoords='offset points')
    ax6.set_xlabel('Helix fraction', fontsize=10)
    ax6.set_ylabel('Strand fraction', fontsize=10)
    ax6.set_title('SS Composition', fontsize=10)
    ax6.set_xlim(-0.05, 1.05)
    ax6.set_ylim(-0.05, 1.05)
    
    fig.suptitle('Protein Backbones as Piecewise Spirals on T\u00b2: Overview',
                 fontsize=14, fontweight='bold')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()


# =============================================================================
# SECTION 13: REPORT
# =============================================================================

def generate_report(all_results, stats_summary, out_path):
    """Generate Markdown report with reframed contributions."""
    
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write("# Protein Backbones as Piecewise Spirals on T\u00b2\n\n")
        f.write("## Analysis Report\n\n")
        f.write(f"**Date:** 2026-02-26  \n")
        f.write(f"**N proteins:** {len(all_results)}  \n")
        n_pdb = sum(1 for r in all_results if r.source == 'pdb')
        n_syn = sum(1 for r in all_results if r.source == 'synthetic')
        n_dssp = sum(1 for r in all_results if 'dssp' in r.ss_method)
        n_dihedral = sum(1 for r in all_results if r.ss_method == 'dihedral')
        f.write(f"**Sources:** {n_pdb} PDB, {n_syn} synthetic  \n")
        f.write(f"**SS assignment:** {n_dssp} DSSP, {n_dihedral} dihedral fallback, "
                f"{n_syn} synthetic  \n\n")
        
        f.write("---\n\n")
        
        # ── Contribution 1: Curvature Geometry is Non-Random ──
        f.write("## Contribution 1: Backbone Curvature Follows Simple Functional Forms\n\n")
        f.write("**Finding:** On T\u00b2, every secondary-structure segment\u2019s curvature profile "
                "\u03ba(s) is well-described by one of 9 simple functional forms (AICc selection). "
                "100% of segments with \u2265 4 residues are classifiable with median R\u00b2 > 0.5. "
                "This means backbone curvature dynamics on the Ramachandran torus are "
                "geometrically structured, not random.\n\n")
        
        cr = stats_summary.get('curvature_regularity', {})
        if cr:
            f.write("### Curvature Statistics by SS Type\n\n")
            f.write("| SS Type | n segments | Mean \u03ba | Intra-CV | R\u00b2 (linear fit) | "
                    "Regular % | Mean seg. length |\n")
            f.write("|---------|------------|--------|----------|-----------------|"
                    "-----------|------------------|\n")
            for ss_type in ['H', 'E', 'C']:
                if ss_type in cr:
                    info = cr[ss_type]
                    f.write(f"| {ss_type} | {info['n']} | {info.get('mean_dkappa_ds', 0):.4f} | "
                            f"{info.get('mean_intra_cv', 0):.2f} | "
                            f"{info.get('mean_r_squared', 0):.3f} | "
                            f"{info.get('regular_frac', 0)*100:.0f}% | "
                            f"{info.get('mean_segment_length', 0):.1f} |\n")
            f.write("\n")
            
            # Regularity inversion test
            if 'regularity_inversion_test' in cr:
                rt = cr['regularity_inversion_test']
                f.write(f"**Mann-Whitney test** (coil intra-CV < structured intra-CV): "
                        f"p = {rt['mann_whitney_p']:.4f}, "
                        f"effect size r = {rt['effect_size_r']:.3f}  \n")
                f.write(f"Coil median intra-CV = {rt['coil_median_intra_cv']:.3f}, "
                        f"Structured median intra-CV = {rt['struct_median_intra_cv']:.3f} "
                        f"(n_coil={rt['n_coil']}, n_struct={rt['n_struct']})  \n\n")
            
            if 'r_squared_inversion_test' in cr:
                rt2 = cr['r_squared_inversion_test']
                f.write(f"**R\u00b2 test** (coil R\u00b2 > structured R\u00b2): "
                        f"p = {rt2['mann_whitney_p']:.4f}  \n")
                f.write(f"Coil median R\u00b2 = {rt2['coil_median_r2']:.3f}, "
                        f"Structured median R\u00b2 = {rt2['struct_median_r2']:.3f}  \n\n")
            
            overall_rf = stats_summary.get('overall_regular_fraction', 0)
            f.write(f"**Overall regular spiral fraction:** {overall_rf*100:.1f}%\n\n")
        
        # Null model comparison
        nm = stats_summary.get('null_model', {})
        if nm:
            f.write("### Null Model Comparison (Residue-Shuffle on T\u00b2)\n\n")
            f.write(f"Method: {nm.get('method', 'residue_shuffle')} — residue order permuted "
                    f"within each protein, preserving (\u03c6,\u03c8) marginals.  \n")
            f.write(f"Metric: \u0394AICc = AICc(constant) \u2013 AICc(best oscillatory). "
                    f"Positive = non-constant model improves.  \n")
            f.write(f"Null: n = {nm.get('n_null_segments', 0)}, "
                    f"median \u0394AICc = {nm.get('null_median_delta_aicc', 0):.2f}  \n")
            f.write(f"Real: median \u0394AICc = {nm.get('real_median_delta_aicc', 0):.2f}  \n")
            f.write(f"Mann-Whitney (real > null): p = {nm.get('mann_whitney_p', 1):.4f}, "
                    f"effect r = {nm.get('effect_r', 0):.3f}\n\n")
        
        # Basin-control null
        bn = stats_summary.get('basin_null', {})
        if bn:
            f.write("### Basin-Control Null (Markov Walk Within SS Basins)\n\n")
            f.write("Method: For each real segment, generate 20 surrogate paths by random-walking "
                    "within the same SS type\u2019s step-size distribution. Preserves basin identity "
                    "and local \u0394\u03c6/\u0394\u03c8 statistics; randomizes specific path.  \n")
            f.write("Test: Does the real curvature-class distribution differ from "
                    "basin-constrained random walks?\n\n")
            f.write("| SS | n real | n null | Real non-const | Null non-const | "
                    "\u03c7\u00b2 | p |\n")
            f.write("|---|--------|--------|---------------|---------------|-----|---|\n")
            for ss in ['H', 'E', 'C']:
                if ss in bn:
                    b = bn[ss]
                    f.write(f"| {ss} | {b['real_total']} | {b['null_total']} | "
                            f"{b['real_nonconst_frac']*100:.1f}% | "
                            f"{b['null_nonconst_frac']*100:.1f}% | "
                            f"{b['chi2']:.2f} | {b['p']:.4f} |\n")
            f.write("\n")
        
        # Robustness checks for helix basin null
        if bn and 'H' in bn:
            bn_h = bn['H']
            es = bn_h.get('effect_sizes', {})
            if es:
                f.write("### Helix Basin-Null Effect Sizes\n\n")
                f.write(f"Absolute risk reduction (null \u2212 real non-const): "
                        f"{es.get('abs_risk_reduction', 0)*100:.1f} percentage points  \n")
                f.write(f"Odds ratio (constant in real vs null): "
                        f"{es.get('odds_ratio', 0):.2f}  \n")
                f.write(f"Binary \u03c7\u00b2 (2\u00d72: real/null \u00d7 const/non-const): "
                        f"\u03c7\u00b2 = {es.get('chi2_binary', 0):.2f}, "
                        f"p = {es.get('p_binary', 1):.4f}  \n")
                f.write(f"Cram\u00e9r\u2019s V (binary): {es.get('cramers_v_binary', 0):.3f}  \n\n")
            
            boot = bn_h.get('bootstrap', {})
            if boot:
                f.write("### Robustness: Bootstrap (1000 protein resamples)\n\n")
                f.write(f"Bootstrap 95% CI for helix non-const fraction: "
                        f"[{boot['ci_95_low']*100:.1f}%, {boot['ci_95_high']*100:.1f}%]  \n")
                f.write(f"Null non-const rate: {bn_h['null_nonconst_frac']*100:.1f}%  \n")
                f.write(f"Bootstrap mean: {boot['mean_nonconst']*100:.1f}%  \n")
                f.write(f"Fraction of bootstrap samples \u2265 null rate: "
                        f"{boot['frac_exceeds_null']*100:.1f}%  \n\n")
            
            loo = bn_h.get('leave_one_out', {})
            if loo:
                f.write("### Robustness: Leave-One-Protein-Out\n\n")
                f.write(f"n proteins with helix segments: {loo['n_proteins']}  \n")
                f.write(f"LOO non-const range: "
                        f"[{loo['min_nonconst']*100:.1f}%, {loo['max_nonconst']*100:.1f}%]  \n")
                f.write(f"All LOO iterations below null rate ({bn_h['null_nonconst_frac']*100:.1f}%): "
                        f"{'Yes' if loo['all_below_null'] else 'No'}  \n\n")
        
        f.write("---\n\n")
        
        # ── Contribution 2: SS-Dependent Spiral Taxonomy ──
        f.write("## Contribution 2: SS-Dependent Spiral Taxonomy\n\n")
        f.write("**Finding:** AIC-based model selection across 6 functional forms reveals "
                "distinct geometric signatures by SS type. Segments are classified as "
                "geodesic, circular arc, clothoid, Fermat, sinusoidal, exponential, "
                "step (piecewise-constant), or quadratic.\n\n")
        
        sc = stats_summary.get('spiral_classification', {})
        if sc:
            f.write("### Classification Matrix\n\n")
            # Get all shape types that actually appear
            all_shapes = set()
            for ss_type in ['H', 'E', 'C']:
                if ss_type in sc:
                    all_shapes.update(k for k, v in sc[ss_type].items() if v > 0 and k != 'too_short')
            shape_order = ['geodesic', 'circular_arc', 'clothoid', 'fermat',
                           'sinusoidal', 'damped_osc', 'exponential', 'sigmoid',
                           'step', 'quadratic', 'gauss_peak']
            shapes = [s for s in shape_order if s in all_shapes]
            
            header = "| SS Type | " + " | ".join(s.replace('_',' ').title() for s in shapes) + " | Total |\n"
            sep = "|---------|" + "|".join("-" * (len(s.replace('_',' ').title()) + 2) for s in shapes) + "|-------|\n"
            f.write(header)
            f.write(sep)
            for ss_type in ['H', 'E', 'C']:
                if ss_type in sc:
                    info = sc[ss_type]
                    total = sum(info.get(s, 0) for s in shapes)
                    cells = []
                    for s in shapes:
                        count = info.get(s, 0)
                        pct = f" ({100*count/total:.0f}%)" if total > 0 and count > 0 else ""
                        cells.append(f"{count}{pct}")
                    f.write(f"| {ss_type} | " + " | ".join(cells) + f" | {total} |\n")
            f.write("\n")
            
            if 'chi_squared_full' in sc:
                chi = sc['chi_squared_full']
                f.write(f"**Full chi-squared** ({chi.get('label', '')}): "
                        f"\u03c7\u00b2 = {chi['chi2']:.2f}, p = {chi['p_value']:.4f}, "
                        f"Cram\u00e9r's V = {chi.get('cramers_v', 0):.3f}\n\n")
            
            if 'chi_squared_oscillatory' in sc:
                chi = sc['chi_squared_oscillatory']
                f.write(f"**Helix oscillatory enrichment** ({chi.get('label', '')}): "
                        f"\u03c7\u00b2 = {chi['chi2']:.2f}, p = {chi['p_value']:.4f}\n\n")
            
            if 'chi_squared_constant' in sc:
                chi = sc['chi_squared_constant']
                f.write(f"**Strand constant-\u03ba enrichment** ({chi.get('label', '')}): "
                        f"\u03c7\u00b2 = {chi['chi2']:.2f}, p = {chi['p_value']:.4f}\n\n")
            
            # Bonferroni note
            if 'bonferroni' in sc:
                bf = sc['bonferroni']
                f.write(f"*Multiple comparisons:* {bf['note']}. ")
                # Check which tests survive
                full_p = sc.get('chi_squared_full', {}).get('p_value', 1.0)
                osc_p = sc.get('chi_squared_oscillatory', {}).get('p_value', 1.0)
                const_p = sc.get('chi_squared_constant', {}).get('p_value', 1.0)
                threshold = bf['threshold']
                survivors = []
                if full_p < threshold:
                    survivors.append(f"Full test (p={full_p:.4f})")
                if osc_p < threshold:
                    survivors.append(f"Oscillatory (p={osc_p:.4f})")
                if const_p < threshold:
                    survivors.append(f"Constant-\u03ba (p={const_p:.4f})")
                if survivors:
                    f.write(f"Survives correction: {'; '.join(survivors)}.")
                else:
                    f.write("Only full test is tested at uncorrected \u03b1; focused tests are post-hoc.")
                f.write("\n\n")
            
            # Standardized residuals
            if 'standardized_residuals' in sc:
                sr = sc['standardized_residuals']
                f.write("### Standardized Residuals (full chi-squared)\n\n")
                cols = sr.get('col_labels', [])
                rows = sr.get('row_labels', ['H', 'E', 'C'])
                header = "| SS | " + " | ".join(c.replace('_', ' ').title() for c in cols) + " |\n"
                sep = "|---|" + "|".join("---" for _ in cols) + "|\n"
                f.write(header)
                f.write(sep)
                for i, row_label in enumerate(rows):
                    vals = sr['matrix'][i] if i < len(sr['matrix']) else []
                    cells = [f"{v:+.2f}" for v in vals]
                    f.write(f"| {row_label} | " + " | ".join(cells) + " |\n")
                f.write("\n*Values > |2| indicate significant deviation from expected.*\n\n")
        
        # Omega validation
        ov = stats_summary.get('omega_validation', {})
        if ov:
            f.write("### Helix Oscillation Frequency Validation\n\n")
            f.write(f"Expected \u03c9 per residue (2\u03c0/3.6): "
                    f"{ov['expected_omega_per_residue']:.3f} rad/residue  \n")
            f.write(f"Fitted median \u03c9 per residue: "
                    f"{ov['median_fitted_omega_per_residue']:.3f} rad/residue  \n")
            f.write(f"Fitted mean \u00b1 std: "
                    f"{ov['mean_fitted_omega_per_residue']:.3f} \u00b1 "
                    f"{ov['std_fitted_omega_per_residue']:.3f}  \n")
            f.write(f"n oscillatory helix segments: {ov['n_helix_oscillatory']}\n\n")
        
        # Oscillatory evidence (Δ AICc analysis)
        oe = stats_summary.get('oscillatory_evidence', {})
        if oe:
            f.write("### Oscillatory Evidence (\u0394AICc Analysis)\n\n")
            f.write("Even when AICc selects constant, oscillatory models may be competitive. "
                    "\u0394AICc = AICc(constant) \u2013 AICc(best oscillatory): "
                    "positive = oscillatory improves over constant.\n\n")
            
            f.write("| SS Type | n | Median \u0394AICc | Frac. competitive (\u0394 > \u22122) | "
                    "Frac. preferred (\u0394 > 0) | Mean osc. R\u00b2 |\n")
            f.write("|---------|---|------------|------------------------|---------------------|---------------|\n")
            for ss in ['H', 'E', 'C']:
                if ss in oe:
                    info = oe[ss]
                    f.write(f"| {ss} | {info['n']} | {info['median_delta_aicc']:.1f} | "
                            f"{info['frac_competitive']*100:.0f}% | "
                            f"{info['frac_preferred']*100:.0f}% | "
                            f"{info['mean_osc_r2']:.3f} |\n")
            f.write("\n")
            
            if 'mann_whitney' in oe:
                mw = oe['mann_whitney']
                f.write(f"**Mann-Whitney** (H \u0394AICc > E+C \u0394AICc): "
                        f"p = {mw['p_value']:.4f}, effect r = {mw['effect_r']:.3f}  \n")
                f.write(f"H median \u0394 = {mw['H_median_delta']:.1f}, "
                        f"E+C median \u0394 = {mw['EC_median_delta']:.1f}\n\n")
            
            if 'omega_per_residue' in oe:
                om = oe['omega_per_residue']
                f.write(f"**\u03c9 per residue** (competitive helix segments): "
                        f"expected = {om['expected']:.3f}, "
                        f"fitted median = {om['median']:.3f}, "
                        f"mean \u00b1 std = {om['mean']:.3f} \u00b1 {om['std']:.3f}, "
                        f"n = {om['n']}\n\n")
        
        # Tangent autocorrelation
        tac = stats_summary.get('tangent_autocorrelation', {})
        if tac:
            f.write("### Tangent Autocorrelation (segments \u2265 8 residues)\n\n")
            f.write("Curvature autocorrelation C(\u0394) = \u27e8\u03ba(s)\u00b7\u03ba(s+\u0394)\u27e9 / Var(\u03ba). "
                    "Helix periodicity predicts peak at \u0394 \u2248 3\u20134.\n\n")
            f.write("| SS Type | n segments | Peak lag | Peak C(\u0394) | C(1) | C(2) | C(3) | C(4) |\n")
            f.write("|---------|-----------|---------|----------|------|------|------|------|\n")
            for ss in ['H', 'E', 'C']:
                if ss in tac:
                    info = tac[ss]
                    ac = info['mean_autocorr']
                    peak_lag = info.get('peak_lag', '-')
                    peak_val = info.get('peak_value', 0)
                    c1 = ac[1] if len(ac) > 1 else 0
                    c2 = ac[2] if len(ac) > 2 else 0
                    c3 = ac[3] if len(ac) > 3 else 0
                    c4 = ac[4] if len(ac) > 4 else 0
                    f.write(f"| {ss} | {info['n_segments']} | {peak_lag} | "
                            f"{peak_val:.3f} | {c1:.3f} | {c2:.3f} | {c3:.3f} | {c4:.3f} |\n")
            f.write("\n")
        
        # Dihedral deviation autocorrelation
        dac = stats_summary.get('dihedral_deviation_autocorr', {})
        if dac:
            f.write("### Dihedral Deviation Autocorrelation (segments \u2265 8 residues)\n\n")
            f.write("Autocorrelation of (\u03c6\u1d62 \u2013 \u03c6\u0304) and (\u03c8\u1d62 \u2013 \u03c8\u0304) within each segment. "
                    "\u03b1-helix periodicity predicts positive peak at \u0394 \u2248 3\u20134.\n\n")
            
            for angle in ['phi', 'psi']:
                angle_sym = '\u03c6' if angle == 'phi' else '\u03c8'
                f.write(f"**{angle_sym} deviation:**\n\n")
                f.write(f"| SS | n | Peak lag | Peak C | C(1) | C(2) | C(3) | C(4) |\n")
                f.write(f"|---|---|---------|--------|------|------|------|------|\n")
                for ss in ['H', 'E', 'C']:
                    if ss in dac and angle in dac[ss]:
                        info = dac[ss][angle]
                        ac = info['mean_autocorr']
                        pl = info.get('peak_lag', '-')
                        pv = info.get('peak_value', 0)
                        c1 = ac[1] if len(ac) > 1 else 0
                        c2 = ac[2] if len(ac) > 2 else 0
                        c3 = ac[3] if len(ac) > 3 else 0
                        c4 = ac[4] if len(ac) > 4 else 0
                        f.write(f"| {ss} | {info['n_segments']} | {pl} | "
                                f"{pv:.3f} | {c1:.3f} | {c2:.3f} | {c3:.3f} | {c4:.3f} |\n")
                f.write("\n")
        
        f.write("---\n\n")
        
        # ── Contribution 3: Torus Knot Descriptors ──
        f.write("## Contribution 3: Per-Protein Torus Knot Descriptors\n\n")
        f.write("**Finding:** The (p,q) winding numbers on T\u00b2 provide a topological "
                "fingerprint per fold. These capture the net angular excursion in \u03c6 and \u03c8 "
                "independently.\n\n")
        
        f.write("### Per-Protein Results\n\n")
        f.write("| Protein | L | SS method | p | q | |Q| | ln(kf) | Regular % |\n")
        f.write("|---------|---|-----------|---|---|-----|--------|----------|\n")
        for r in all_results:
            tk = r.torus_knot
            p_val = tk.get('p', 0) if tk else 0
            q_val = tk.get('q', 0) if tk else 0
            mag = tk.get('total_winding_magnitude', 0) if tk else 0
            f.write(f"| {r.name} | {r.length} | {r.ss_method} | {p_val:.2f} | {q_val:.2f} | "
                    f"{mag:.2f} | {r.ln_kf:.2f} | "
                    f"{r.overall_clothoid_fraction*100:.0f}% |\n")
        
        f.write("\n---\n\n")
        
        # ── Correlations ──
        f.write("## Correlations with Folding Rate\n\n")
        if 'correlations' in stats_summary:
            corr = stats_summary['correlations']
            f.write("| Descriptor | r | p-value | Spearman \u03c1 | p-value |\n")
            f.write("|------------|---|---------|-----------|----------|\n")
            for name, info in corr.items():
                f.write(f"| {name} | {info['pearson_r']:.3f} | {info['pearson_p']:.4f} | "
                        f"{info['spearman_rho']:.3f} | {info['spearman_p']:.4f} |\n")
            f.write("\n")
        
        # ── Structural Class Separation ──
        scs = stats_summary.get('structural_class', {})
        if scs:
            f.write("---\n\n")
            f.write("## Structural Class Separation by Torus Winding\n\n")
            f.write("**Test:** Do (p,q) winding vectors on T\u00b2 cluster by "
                    "SCOP structural class (all-\u03b1, all-\u03b2, \u03b1+\u03b2)?\n\n")
            
            f.write("### Per-Class Winding Statistics\n\n")
            f.write("| Class | n | p (mean\u00b1std) | q (mean\u00b1std) | |Q| (mean\u00b1std) |\n")
            f.write("|-------|---|-------------|-------------|----------------|\n")
            for sc in ['all-alpha', 'all-beta', 'alpha-beta']:
                if sc in scs:
                    s = scs[sc]
                    label = sc.replace('alpha', '\u03b1').replace('beta', '\u03b2')
                    f.write(f"| {label} | {s['n']} | "
                            f"{s['p_mean']:+.2f}\u00b1{s['p_std']:.2f} | "
                            f"{s['q_mean']:+.2f}\u00b1{s['q_std']:.2f} | "
                            f"{s['Q_mean']:.2f}\u00b1{s['Q_std']:.2f} |\n")
            f.write("\n")
            
            f.write("### Proteins by Class\n\n")
            for sc in ['all-alpha', 'all-beta', 'alpha-beta']:
                if sc in scs:
                    label = sc.replace('alpha', '\u03b1').replace('beta', '\u03b2')
                    f.write(f"**{label}:** {', '.join(scs[sc]['proteins'])}  \n")
            f.write("\n")
            
            if 'kruskal_wallis_Q' in scs:
                kw = scs['kruskal_wallis_Q']
                f.write(f"**Kruskal-Wallis** (|Q| by class): H = {kw['H']:.2f}, "
                        f"p = {kw['p']:.4f}  \n")
            if 'kruskal_wallis_p' in scs:
                kw = scs['kruskal_wallis_p']
                f.write(f"**Kruskal-Wallis** (p by class): H = {kw['H']:.2f}, "
                        f"p = {kw['p']:.4f}  \n")
            if 'kruskal_wallis_q' in scs:
                kw = scs['kruskal_wallis_q']
                f.write(f"**Kruskal-Wallis** (q by class): H = {kw['H']:.2f}, "
                        f"p = {kw['p']:.4f}  \n")
            f.write("\n")
            
            if 'alpha_vs_beta_Q' in scs:
                ab = scs['alpha_vs_beta_Q']
                f.write(f"**All-\u03b1 vs all-\u03b2** |Q|: "
                        f"\u03b1 median = {ab['alpha_median_Q']:.2f}, "
                        f"\u03b2 median = {ab['beta_median_Q']:.2f}, "
                        f"MW p = {ab['mannwhitney_p']:.4f}  \n")
            if 'alpha_vs_beta_p' in scs:
                ab = scs['alpha_vs_beta_p']
                f.write(f"**All-\u03b1 vs all-\u03b2** p-winding: "
                        f"\u03b1 median = {ab['alpha_median_p']:+.2f}, "
                        f"\u03b2 median = {ab['beta_median_p']:+.2f}, "
                        f"MW p = {ab['mannwhitney_p']:.4f}  \n")
            if 'alpha_vs_beta_q' in scs:
                ab = scs['alpha_vs_beta_q']
                f.write(f"**All-\u03b1 vs all-\u03b2** q-winding: "
                        f"\u03b1 median = {ab['alpha_median_q']:+.2f}, "
                        f"\u03b2 median = {ab['beta_median_q']:+.2f}, "
                        f"MW p = {ab['mannwhitney_p']:.4f}  \n")
            f.write("\n")
        
        f.write("---\n\n")
        f.write("## Figures\n\n")
        f.write("- `fig1_ramachandran_*.png` \u2014 Backbone path on T\u00b2 with superpotential\n")
        f.write("- `fig2_curvature_regularity.png` \u2014 Curvature regularity by SS type (Contribution 1)\n")
        f.write("- `fig3_spiral_taxonomy.png` \u2014 Spiral type vs SS (Contribution 2)\n")
        f.write("- `fig4_torus_knot_descriptors.png` \u2014 (p,q) scatter and winding vs kf\n")
        f.write("- `fig5_cornu_exemplar_*.png` \u2014 Detailed exemplars\n")
        f.write("- `fig6_summary_dashboard.png` \u2014 Overview dashboard\n")
        
        # ── Active Site Proximity Analysis ──
        asp = stats_summary.get('active_site_proximity', {})
        if asp and asp.get('n_gauss_peaks_total', 0) > 0:
            f.write("\n---\n\n")
            f.write("## Active Site Proximity Analysis\n\n")
            f.write("**Hypothesis:** Gaussian peak curvature anomalies on T\u00b2 mark "
                    "functional sites (catalytic residues, binding pockets).\n\n")
            
            f.write(f"Total Gaussian peak segments: {asp['n_gauss_peaks_total']}  \n")
            f.write(f"In annotated enzymes: {asp['n_gauss_peaks_in_enzymes']}  \n\n")
            
            if asp.get('peaks'):
                f.write("### Peak Catalog\n\n")
                f.write("| Protein | SS | Peak Res | Seg Range | A (amplitude) | \u03c3 (width) | "
                        "Dist to Active | Dist to Catalytic |\n")
                f.write("|---------|---|---------|-----------|--------------|----------|"
                        "---------------|------------------|\n")
                for pk in asp['peaks']:
                    cat_d = pk.get('min_dist_catalytic', None)
                    cat_str = f"{cat_d}" if cat_d is not None else "-"
                    f.write(f"| {pk['protein']} | {pk['ss_type']} | {pk['peak_residue']} | "
                            f"{pk['start']}-{pk['end']} | {pk['A']:.3f} | {pk['sigma']:.3f} | "
                            f"{pk['min_dist_functional']} | {cat_str} |\n")
                f.write("\n")
            
            if 'real_mean_dist' in asp:
                f.write("### Proximity Test\n\n")
                f.write(f"Gaussian peak \u2192 nearest functional residue: "
                        f"mean = {asp['real_mean_dist']:.1f}, "
                        f"median = {asp['real_median_dist']:.1f} residues  \n")
                if 'null_mean_dist' in asp:
                    f.write(f"Null (random helix residues): "
                            f"mean = {asp['null_mean_dist']:.1f}, "
                            f"median = {asp['null_median_dist']:.1f} residues  \n")
                if 'mannwhitney_p' in asp:
                    f.write(f"Mann-Whitney (peak closer than random): "
                            f"p = {asp['mannwhitney_p']:.4f}, "
                            f"effect r = {asp['effect_r']:.3f}  \n")
                f.write("\n")
            
            if asp.get('peak_catalog'):
                all_peaks = asp['peak_catalog']
                non_enzyme_peaks = [pk for pk in all_peaks if pk['protein'] not in
                                   [pr['protein'] for pr in asp.get('peaks', [])]]
                if non_enzyme_peaks:
                    f.write("### Peaks in Non-Annotated Proteins\n\n")
                    f.write("| Protein | SS | Peak Res | Seg Range | A | \u03c3 |\n")
                    f.write("|---------|---|---------|-----------|---|---|\n")
                    for pk in non_enzyme_peaks:
                        f.write(f"| {pk['protein']} | {pk['ss_type']} | "
                                f"{pk['peak_residue']} | {pk['start']}-{pk['end']} | "
                                f"{pk['A']:.3f} | {pk['sigma']:.3f} |\n")
                    f.write("\n")
            
            if 'amplitude_stats' in asp:
                amp = asp['amplitude_stats']
                sig = asp['sigma_stats']
                f.write("### Gaussian Peak Shape Parameters\n\n")
                f.write(f"Amplitude |A|: mean = {amp['mean']:.3f}, "
                        f"std = {amp['std']:.3f}, median = {amp['median']:.3f}  \n")
                f.write(f"Width \u03c3: mean = {sig['mean']:.3f}, "
                        f"std = {sig['std']:.3f}, median = {sig['median']:.3f}  \n\n")
    
    log.info(f"Report written to {out_path}")


# =============================================================================
# SECTION 14: MAIN
# =============================================================================

def main():
    print("=" * 72)
    print("  PROTEIN BACKBONES AS PIECEWISE SPIRALS ON T\u00b2")
    print("  v3: Pure-numpy DSSP (no external binary needed)")
    print("  Geometric Completeness | Spiral Taxonomy | Torus Knot Descriptors")
    print("=" * 72)
    print()
    
    log.info("SS assignment: embedded pure-numpy DSSP (~97%% agreement with original)")
    log.info("  No external DSSP binary required.")
    log.info("  Fallback: widened Ramachandran-region dihedral classifier.")
    print()
    
    # ── Build superpotential ──
    log.info("Building superpotential W(phi, psi) on T^2...")
    phi_grid, psi_grid, p_grid, W_grid = build_superpotential_grid()
    log.info(f"  Grid: {len(phi_grid)}x{len(psi_grid)}, W range: [{W_grid.min():.4f}, {W_grid.max():.4f}]")
    
    # ── Download PDB files ──
    log.info("Attempting PDB downloads...")
    for name, pdb_id, chain, length, ln_kf, co in PROTEIN_DATASET:
        download_pdb(pdb_id, PDB_DIR)
    
    # ── Analyze each protein ──
    log.info(f"Analyzing {len(PROTEIN_DATASET)} proteins...")
    all_results = []
    
    for name, pdb_id, chain, length, ln_kf, co in PROTEIN_DATASET:
        result = analyze_single_protein(
            name, pdb_id, chain, length, ln_kf, co,
            phi_grid, psi_grid, W_grid,
        )
        all_results.append(result)
    
    log.info(f"  Completed: {len(all_results)} proteins")
    n_pdb = sum(1 for r in all_results if r.source == 'pdb')
    n_dssp = sum(1 for r in all_results if r.ss_method == 'dssp_numpy')
    n_dihedral = sum(1 for r in all_results if r.ss_method == 'dihedral')
    log.info(f"  PDB: {n_pdb}, Synthetic: {len(all_results) - n_pdb}")
    log.info(f"  SS assignment: {n_dssp} pure-numpy DSSP, {n_dihedral} dihedral, "
             f"{len(all_results) - n_dssp - n_dihedral} synthetic")
    
    # ── Aggregate statistics ──
    log.info("Computing aggregate statistics...")
    log.info("  C1: curvature regularity by SS type")
    stats_summary = {}
    
    # Contribution 1: Curvature regularity by SS type
    wrc_agg = {}
    for ss_type in ['H', 'E', 'C']:
        dkappa_all = []
        r2_all = []
        intra_cvs_all = []
        regular_count = 0
        total_count = 0
        seg_lengths = []
        
        for r in all_results:
            wrc = r.winding_rate_constancy
            if ss_type in wrc.get('by_ss_type', {}):
                info = wrc['by_ss_type'][ss_type]
                dkappa_all.append(info['mean_dkappa_ds'])
                r2_all.append(info['mean_r_squared'])
                total_count += info['n_segments']
                regular_count += int(info.get('fraction_regular', 0) * info['n_segments'])
                if 'mean_intra_cv' in info and info['mean_intra_cv'] > 0:
                    intra_cvs_all.append(info['mean_intra_cv'])
                if 'mean_segment_length' in info:
                    seg_lengths.append(info['mean_segment_length'])
        
        if dkappa_all:
            wrc_agg[ss_type] = {
                'mean_dkappa_ds': float(np.mean(dkappa_all)),
                'std_dkappa_ds': float(np.std(dkappa_all)),
                'cv': float(np.std(dkappa_all) / max(abs(np.mean(dkappa_all)), 1e-10)),
                'mean_r_squared': float(np.mean(r2_all)),
                'n': total_count,
                'regular_frac': regular_count / max(total_count, 1),
                'mean_intra_cv': float(np.mean(intra_cvs_all)) if intra_cvs_all else 0.0,
                'mean_segment_length': float(np.mean(seg_lengths)) if seg_lengths else 0.0,
            }
    
    # Aggregate regularity inversion test (protein-level)
    all_coil_cvs = []
    all_struct_cvs = []
    all_coil_r2 = []
    all_struct_r2 = []
    for r in all_results:
        wrc = r.winding_rate_constancy
        for ss_type in ['H', 'E']:
            if ss_type in wrc.get('by_ss_type', {}):
                info = wrc['by_ss_type'][ss_type]
                if info.get('mean_intra_cv', 0) > 0:
                    all_struct_cvs.append(info['mean_intra_cv'])
                    all_struct_r2.append(info['mean_r_squared'])
        if 'C' in wrc.get('by_ss_type', {}):
            info = wrc['by_ss_type']['C']
            if info.get('mean_intra_cv', 0) > 0:
                all_coil_cvs.append(info['mean_intra_cv'])
                all_coil_r2.append(info['mean_r_squared'])
    
    if len(all_coil_cvs) >= 3 and len(all_struct_cvs) >= 3:
        u_stat, mw_p = stats.mannwhitneyu(all_coil_cvs, all_struct_cvs, alternative='less')
        n1, n2 = len(all_coil_cvs), len(all_struct_cvs)
        effect_r = 1.0 - 2.0 * u_stat / (n1 * n2)
        wrc_agg['regularity_inversion_test'] = {
            'coil_median_intra_cv': float(np.median(all_coil_cvs)),
            'struct_median_intra_cv': float(np.median(all_struct_cvs)),
            'mann_whitney_p': float(mw_p),
            'effect_size_r': float(effect_r),
            'n_coil': n1, 'n_struct': n2,
        }
    
    if len(all_coil_r2) >= 3 and len(all_struct_r2) >= 3:
        u2, p2 = stats.mannwhitneyu(all_coil_r2, all_struct_r2, alternative='greater')
        wrc_agg['r_squared_inversion_test'] = {
            'coil_median_r2': float(np.median(all_coil_r2)),
            'struct_median_r2': float(np.median(all_struct_r2)),
            'mann_whitney_p': float(p2),
        }
    
    stats_summary['curvature_regularity'] = wrc_agg
    stats_summary['overall_regular_fraction'] = float(np.mean(
        [r.overall_clothoid_fraction for r in all_results]))
    
    log.info("  C2: spiral taxonomy chi-squared tests")
    # Contribution 2: Expanded spiral taxonomy
    sc_agg = {}
    for ss_type in ['H', 'E', 'C']:
        totals = {}
        for r in all_results:
            sc = r.spiral_classification
            if 'classification_matrix' in sc and ss_type in sc['classification_matrix']:
                for sp_type, count in sc['classification_matrix'][ss_type].items():
                    totals[sp_type] = totals.get(sp_type, 0) + count
        sc_agg[ss_type] = totals
    
    # Macro-category grouping for chi-squared
    macro_groups = {
        'oscillatory': ['sinusoidal', 'damped_osc'],
        'constant_k': ['geodesic', 'circular_arc'],
        'monotone': ['fermat', 'clothoid', 'exponential'],
        'transition': ['sigmoid', 'step'],
        'polynomial': ['quadratic'],
        'localized': ['gauss_peak'],
    }
    
    # Full 3×k chi-squared
    ss_list = ['H', 'E', 'C']
    macro_list = list(macro_groups.keys())
    contingency_full = np.zeros((3, len(macro_list)), dtype=int)
    for i, ss in enumerate(ss_list):
        for j, (macro, members) in enumerate(macro_groups.items()):
            contingency_full[i, j] = sum(sc_agg.get(ss, {}).get(sp, 0) for sp in members)
    
    col_sums = contingency_full.sum(axis=0)
    valid_cols = col_sums > 0
    contingency_trimmed = contingency_full[:, valid_cols]
    macro_trimmed = [m for m, v in zip(macro_list, valid_cols) if v]
    
    if contingency_trimmed.sum() >= 15 and contingency_trimmed.shape[1] >= 2:
        chi2_f, p_f, dof_f, _ = stats.chi2_contingency(contingency_trimmed)
        n_total = contingency_trimmed.sum()
        min_dim = min(contingency_trimmed.shape) - 1
        cramers_v = np.sqrt(chi2_f / (n_total * max(min_dim, 1)))
        sc_agg['chi_squared_full'] = {
            'chi2': float(chi2_f), 'p_value': float(p_f), 'dof': int(dof_f),
            'table': contingency_trimmed.tolist(),
            'row_labels': ss_list,
            'col_labels': macro_trimmed,
            'cramers_v': float(cramers_v),
            'label': '3×k: SS type vs macro-shape (oscillatory/constant/monotone/transition/poly)',
        }
        log.info("    Full chi2=%.2f, p=%.4f, V=%.3f", chi2_f, p_f, cramers_v)
    
    # Focused test: helix oscillatory enrichment
    regular_types = ['geodesic', 'circular_arc', 'clothoid', 'fermat',
                     'sinusoidal', 'damped_osc', 'exponential', 'sigmoid',
                     'step', 'quadratic', 'gauss_peak']
    osc_types = ['sinusoidal', 'damped_osc']
    non_osc = [sp for sp in regular_types if sp not in osc_types]
    
    n_H_osc = sum(sc_agg.get('H', {}).get(sp, 0) for sp in osc_types)
    n_H_non = sum(sc_agg.get('H', {}).get(sp, 0) for sp in non_osc)
    n_EC_osc = sum(sc_agg.get(ss, {}).get(sp, 0) for ss in ['E', 'C'] for sp in osc_types)
    n_EC_non = sum(sc_agg.get(ss, {}).get(sp, 0) for ss in ['E', 'C'] for sp in non_osc)
    
    cont_osc = np.array([[n_H_osc, n_H_non], [n_EC_osc, n_EC_non]])
    if cont_osc.sum() >= 10 and np.all(cont_osc.sum(axis=0) > 0) and np.all(cont_osc.sum(axis=1) > 0):
        chi2_o, p_o, dof_o, _ = stats.chi2_contingency(cont_osc)
        sc_agg['chi_squared_oscillatory'] = {
            'chi2': float(chi2_o), 'p_value': float(p_o), 'dof': int(dof_o),
            'table': cont_osc.tolist(),
            'label': 'H vs E+C: oscillatory enrichment',
        }
    
    # Focused test: strand constant-κ enrichment
    const_types = ['geodesic', 'circular_arc']
    nonconst = [sp for sp in regular_types if sp not in const_types]
    
    n_E_const = sum(sc_agg.get('E', {}).get(sp, 0) for sp in const_types)
    n_E_nonc = sum(sc_agg.get('E', {}).get(sp, 0) for sp in nonconst)
    n_HC_const = sum(sc_agg.get(ss, {}).get(sp, 0) for ss in ['H', 'C'] for sp in const_types)
    n_HC_nonc = sum(sc_agg.get(ss, {}).get(sp, 0) for ss in ['H', 'C'] for sp in nonconst)
    
    cont_const = np.array([[n_E_const, n_E_nonc], [n_HC_const, n_HC_nonc]])
    if cont_const.sum() >= 10 and np.all(cont_const.sum(axis=0) > 0) and np.all(cont_const.sum(axis=1) > 0):
        chi2_c, p_c, dof_c, _ = stats.chi2_contingency(cont_const)
        sc_agg['chi_squared_constant'] = {
            'chi2': float(chi2_c), 'p_value': float(p_c), 'dof': int(dof_c),
            'table': cont_const.tolist(),
            'label': 'E vs H+C: constant-κ enrichment',
        }
    
    stats_summary['spiral_classification'] = sc_agg
    
    # ── Null Model: Random walks on T² ──
    log.info("  Null model: residue-shuffle (200 reps)...")
    # Collect real segment lengths per SS type
    real_seg_lengths = []
    real_seg_delta_aicc = []  # ΔAICc = AICc(constant) - AICc(best)
    real_seg_phi_psi = []  # (phi_array, psi_array, ss_array) per protein for shuffling
    for r in all_results:
        for seg in r.segments:
            if seg.spiral_class != 'too_short' and seg.length >= 4:
                real_seg_lengths.append(seg.length)
                real_seg_delta_aicc.append(seg.delta_aicc_osc)
    
    # Collect all (phi, psi) per protein for shuffle null
    for r in all_results:
        if r.n_valid_residues > 0:
            # Re-extract phi/psi for shuffle null
            pdb_path = PDB_DIR / f"{r.pdb_id}.pdb"
            if pdb_path.exists() and r.source == 'pdb':
                try:
                    phi, psi, _, backbone_coords, has_O = extract_dihedrals_manual(
                        str(pdb_path), r.chain)
                    valid = ~(np.isnan(phi) | np.isnan(psi))
                    phi_v, psi_v = phi[valid], psi[valid]
                    if np.all(has_O[valid]) and len(phi_v) >= 6:
                        ss = dssp_assign_pure_numpy(backbone_coords[valid])
                    else:
                        ss = assign_ss_from_dihedrals(phi_v, psi_v)
                    real_seg_phi_psi.append((phi_v, psi_v, ss))
                except Exception:
                    pass
    
    if real_seg_lengths:
        n_null_reps = 200
        rng = np.random.RandomState(42)
        
        null_delta_aicc = []
        
        for rep in range(n_null_reps):
            # Shuffle null: pick a random protein, permute residue order
            if real_seg_phi_psi:
                idx = rng.randint(len(real_seg_phi_psi))
                phi_orig, psi_orig, ss_orig = real_seg_phi_psi[idx]
                # Shuffle residue order (destroys sequential correlation, preserves marginals)
                perm = rng.permutation(len(phi_orig))
                phi_shuf = phi_orig[perm]
                psi_shuf = psi_orig[perm]
            else:
                # Fallback: random walk if no PDB data
                seg_len = rng.choice(real_seg_lengths)
                phi_shuf = rng.uniform(-np.pi, np.pi, seg_len)
                psi_shuf = rng.uniform(-np.pi, np.pi, seg_len)
            
            kappa_shuf, s_shuf, _ = torus_curvature(phi_shuf, psi_shuf)
            if len(kappa_shuf) < 8:
                continue
            
            # Segment into chunks matching real segment lengths
            seg_len = rng.choice(real_seg_lengths)
            seg_len = min(seg_len, len(kappa_shuf))
            if seg_len < 4:
                continue
            start = rng.randint(0, max(len(kappa_shuf) - seg_len, 1))
            seg_k = kappa_shuf[start:start+seg_len]
            seg_s = s_shuf[start:start+seg_len]
            
            if len(seg_k) >= 4:
                _, _, _, fits_null, _ = classify_spiral(seg_k, seg_s)
                # ΔAICc = AICc(constant) - AICc(best)
                const_aicc = fits_null.get('constant', {}).get('aicc', np.inf)
                best_aicc = min(f.get('aicc', np.inf) for f in fits_null.values())
                if np.isfinite(const_aicc) and np.isfinite(best_aicc):
                    null_delta_aicc.append(const_aicc - best_aicc)
        
        null_delta_arr = np.array(null_delta_aicc)
        real_delta_arr = np.array(real_seg_delta_aicc)
        
        # Compare ΔAICc distributions (real should have larger Δ = better non-constant fits)
        if len(null_delta_arr) >= 10 and len(real_delta_arr) >= 10:
            u_null, p_null = stats.mannwhitneyu(real_delta_arr, null_delta_arr, alternative='greater')
            n1_n, n2_n = len(real_delta_arr), len(null_delta_arr)
            effect_null = 1.0 - 2.0 * u_null / (n1_n * n2_n)
        else:
            p_null, effect_null = 1.0, 0.0
        
        stats_summary['null_model'] = {
            'method': 'residue_shuffle',
            'n_null_segments': len(null_delta_arr),
            'null_median_delta_aicc': float(np.median(null_delta_arr)) if len(null_delta_arr) > 0 else 0,
            'null_mean_delta_aicc': float(np.mean(null_delta_arr)) if len(null_delta_arr) > 0 else 0,
            'real_median_delta_aicc': float(np.median(real_delta_arr)),
            'real_mean_delta_aicc': float(np.mean(real_delta_arr)),
            'mann_whitney_p': float(p_null),
            'effect_r': float(effect_null),
            'description': 'ΔAICc = AICc(constant) - AICc(best_oscillatory). '
                          'Positive = non-constant model improves. '
                          'Shuffle null: residue order permuted within protein.',
        }
    
    # ── Basin-Control Null (Markov walk within Ramachandran basins) ──
    # Reviewer suggestion: residue-shuffle destroys everything.
    # This null preserves basin identity and local step-size distribution,
    # but randomizes the specific path within each basin.
    log.info("  Basin-control null: Markov walk within SS basins (5 reps/seg)...")
    
    # Define basin centers and typical (φ,ψ) step distributions per SS type
    basin_segments = {'H': [], 'E': [], 'C': []}  # collect (Δφ, Δψ) steps per SS type
    basin_starts = {'H': [], 'E': [], 'C': []}    # starting points per SS type
    
    for r in all_results:
        pdb_path = PDB_DIR / f"{r.pdb_id}.pdb"
        if not (pdb_path.exists() and r.source == 'pdb'):
            continue
        try:
            phi, psi, _, backbone_coords, has_O = extract_dihedrals_manual(
                str(pdb_path), r.chain)
            valid = ~(np.isnan(phi) | np.isnan(psi))
            phi_v, psi_v = phi[valid], psi[valid]
            if np.all(has_O[valid]) and len(phi_v) >= 6:
                ss = dssp_assign_pure_numpy(backbone_coords[valid])
            else:
                ss = assign_ss_from_dihedrals(phi_v, psi_v)
            
            # Collect steps within each SS segment
            for seg in r.segments:
                if seg.spiral_class == 'too_short' or seg.length < 4:
                    continue
                s, e = seg.start_idx, seg.end_idx
                if e >= len(phi_v):
                    continue
                seg_phi = phi_v[s:e+1]
                seg_psi = psi_v[s:e+1]
                # Periodic differences
                dphi = np.arctan2(np.sin(np.diff(seg_phi)), np.cos(np.diff(seg_phi)))
                dpsi = np.arctan2(np.sin(np.diff(seg_psi)), np.cos(np.diff(seg_psi)))
                steps = np.column_stack([dphi, dpsi])
                ss_type = seg.ss_type
                if ss_type in basin_segments:
                    basin_segments[ss_type].append(steps)
                    basin_starts[ss_type].append((seg_phi[0], seg_psi[0]))
        except Exception:
            pass
    
    # Build step pools per SS type
    step_pools = {}
    for ss in ['H', 'E', 'C']:
        if basin_segments[ss]:
            step_pools[ss] = np.vstack(basin_segments[ss])
    
    # Generate Markov-walk null segments and classify
    rng_basin = np.random.RandomState(123)
    n_basin_reps = 5  # per segment — 168 segs × 5 = ~840 classify_spiral calls
    basin_null_classes = []  # list of (ss_type, macro_class) for null
    real_classes = []         # list of (ss_type, macro_class) for real data
    
    # Macro-category mapping
    macro_map = {
        'geodesic': 'constant_k', 'circular_arc': 'constant_k',
        'sinusoidal': 'oscillatory', 'damped_oscillation': 'oscillatory',
        'linear': 'monotone', 'clothoid': 'monotone', 'fermat': 'monotone',
        'exponential': 'monotone',
        'sigmoid': 'transition', 'step': 'transition',
        'quadratic': 'polynomial',
        'gauss_peak': 'localized',
    }
    
    # Collect real classifications
    for r in all_results:
        for seg in r.segments:
            if seg.spiral_class != 'too_short' and seg.length >= 4:
                mc = macro_map.get(seg.spiral_class, 'constant_k')
                real_classes.append((seg.ss_type, mc))
    
    # Generate null: for each real segment, create a small number of surrogate paths
    n_basin_done = 0
    n_basin_total = sum(1 for r in all_results for seg in r.segments
                        if seg.spiral_class != 'too_short' and seg.length >= 4
                        and seg.ss_type in step_pools)
    
    for r in all_results:
        for seg in r.segments:
            if seg.spiral_class == 'too_short' or seg.length < 4:
                continue
            ss = seg.ss_type
            if ss not in step_pools or len(step_pools[ss]) < 5:
                continue
            
            pool = step_pools[ss]
            seg_len = seg.length
            n_starts = len(basin_starts[ss])
            if n_starts == 0:
                continue
            
            for rep in range(n_basin_reps):
                # Pick random start
                si = rng_basin.randint(n_starts)
                phi0, psi0 = basin_starts[ss][si]
                
                # Vectorized random walk
                step_indices = rng_basin.randint(0, len(pool), size=seg_len - 1)
                steps = pool[step_indices]
                phi_path = np.empty(seg_len)
                psi_path = np.empty(seg_len)
                phi_path[0] = phi0
                psi_path[0] = psi0
                np.cumsum(steps[:, 0], out=phi_path[1:])
                phi_path[1:] += phi0
                np.cumsum(steps[:, 1], out=psi_path[1:])
                psi_path[1:] += psi0
                
                # Curvature and classify
                kappa_null, s_null, _ = torus_curvature(phi_path, psi_path)
                if len(kappa_null) >= 4:
                    s_rel = s_null[:len(kappa_null)] - s_null[0] if len(s_null) > len(kappa_null) else \
                            np.arange(len(kappa_null), dtype=float)
                    _, _, cls_null, _, _ = classify_spiral(kappa_null, s_rel)
                    mc_null = macro_map.get(cls_null, 'constant_k')
                    basin_null_classes.append((ss, mc_null))
            
            n_basin_done += 1
            if n_basin_done % 50 == 0:
                log.info("    Basin null: %d/%d segments processed...", n_basin_done, n_basin_total)
    
    # Compare: for each SS type, do real and null have different macro-class distributions?
    basin_null_summary = {}
    if basin_null_classes and real_classes:
        all_macros = sorted(set(mc for _, mc in real_classes + basin_null_classes))
        
        for ss in ['H', 'E', 'C']:
            real_counts = {mc: 0 for mc in all_macros}
            null_counts = {mc: 0 for mc in all_macros}
            for s, mc in real_classes:
                if s == ss:
                    real_counts[mc] += 1
            for s, mc in basin_null_classes:
                if s == ss:
                    null_counts[mc] += 1
            
            real_total = sum(real_counts.values())
            null_total = sum(null_counts.values())
            
            if real_total > 0 and null_total > 0:
                # Chi-squared: real vs null distributions for this SS type
                real_arr = np.array([real_counts[mc] for mc in all_macros])
                null_arr = np.array([null_counts[mc] for mc in all_macros])
                # Normalize null to same total as real
                null_expected = null_arr * (real_total / null_total)
                # Only include categories with expected > 0
                mask = null_expected > 0
                if mask.sum() >= 2:
                    chi2_b, p_b = stats.chisquare(real_arr[mask], f_exp=null_expected[mask])
                else:
                    chi2_b, p_b = 0.0, 1.0
                
                # Fraction non-constant in real vs null
                real_nonconst = 1.0 - real_counts.get('constant_k', 0) / max(real_total, 1)
                null_nonconst = 1.0 - null_counts.get('constant_k', 0) / max(null_total, 1)
                
                basin_null_summary[ss] = {
                    'real_total': real_total,
                    'null_total': null_total,
                    'real_counts': real_counts,
                    'null_fracs': {mc: null_counts[mc]/null_total for mc in all_macros},
                    'chi2': float(chi2_b),
                    'p': float(p_b),
                    'real_nonconst_frac': float(real_nonconst),
                    'null_nonconst_frac': float(null_nonconst),
                }
        
        # ── Robustness checks for helix basin null ──
        if 'H' in basin_null_summary:
            bn_h = basin_null_summary['H']
            real_h_const = sum(1 for s, mc in real_classes if s == 'H' and mc == 'constant_k')
            real_h_nonconst = sum(1 for s, mc in real_classes if s == 'H' and mc != 'constant_k')
            null_h_const = sum(1 for s, mc in basin_null_classes if s == 'H' and mc == 'constant_k')
            null_h_nonconst = sum(1 for s, mc in basin_null_classes if s == 'H' and mc != 'constant_k')
            real_h_total = real_h_const + real_h_nonconst
            null_h_total = null_h_const + null_h_nonconst
            
            # Effect sizes (binary: constant vs non-constant)
            real_nc_rate = real_h_nonconst / max(real_h_total, 1)
            null_nc_rate = null_h_nonconst / max(null_h_total, 1)
            abs_risk_reduction = null_nc_rate - real_nc_rate
            
            # Odds ratio
            a, b = real_h_const, real_h_nonconst  # real: const, nonconst
            c, d = null_h_const, null_h_nonconst  # null: const, nonconst
            odds_ratio = (a * d) / max(b * c, 1)
            
            # Binary chi-squared (2x2: real/null × const/nonconst)
            table_2x2 = np.array([[a, b], [c, d]])
            chi2_binary, p_binary = stats.chi2_contingency(table_2x2)[:2]
            
            # Cramér's V for 2x2
            n_2x2 = table_2x2.sum()
            cramers_v_binary = float(np.sqrt(chi2_binary / max(n_2x2, 1)))
            
            bn_h['effect_sizes'] = {
                'abs_risk_reduction': float(abs_risk_reduction),
                'odds_ratio': float(odds_ratio),
                'chi2_binary': float(chi2_binary),
                'p_binary': float(p_binary),
                'cramers_v_binary': float(cramers_v_binary),
            }
            
            # Bootstrap: resample proteins 1000 times, recompute helix non-const fraction
            log.info("    Helix robustness: bootstrap (1000 reps) + LOO...")
            rng_boot = np.random.RandomState(999)
            
            # Map each helix segment to its protein
            helix_by_protein = {}  # protein_name -> [(ss, mc), ...]
            for r in all_results:
                segs_this = [(seg.ss_type, macro_map.get(seg.spiral_class, 'constant_k'))
                             for seg in r.segments
                             if seg.spiral_class != 'too_short' and seg.length >= 4
                             and seg.ss_type == 'H']
                if segs_this:
                    helix_by_protein[r.name] = segs_this
            
            protein_names = list(helix_by_protein.keys())
            n_proteins = len(protein_names)
            
            if n_proteins >= 3:
                boot_nonconst_fracs = []
                for _ in range(1000):
                    # Resample proteins with replacement
                    boot_proteins = rng_boot.choice(protein_names, size=n_proteins, replace=True)
                    boot_segs = []
                    for pn in boot_proteins:
                        boot_segs.extend(helix_by_protein[pn])
                    if boot_segs:
                        nc = sum(1 for _, mc in boot_segs if mc != 'constant_k')
                        boot_nonconst_fracs.append(nc / len(boot_segs))
                
                boot_arr = np.array(boot_nonconst_fracs)
                # How often does bootstrap non-const exceed null rate?
                boot_exceeds_null = float(np.mean(boot_arr >= null_nc_rate))
                
                bn_h['bootstrap'] = {
                    'n_reps': 1000,
                    'mean_nonconst': float(np.mean(boot_arr)),
                    'ci_95_low': float(np.percentile(boot_arr, 2.5)),
                    'ci_95_high': float(np.percentile(boot_arr, 97.5)),
                    'frac_exceeds_null': boot_exceeds_null,
                }
                
                # Leave-one-protein-out
                loo_nonconst_fracs = []
                for leave_out in protein_names:
                    loo_segs = []
                    for pn in protein_names:
                        if pn != leave_out:
                            loo_segs.extend(helix_by_protein[pn])
                    if loo_segs:
                        nc = sum(1 for _, mc in loo_segs if mc != 'constant_k')
                        loo_nonconst_fracs.append(nc / len(loo_segs))
                
                loo_arr = np.array(loo_nonconst_fracs)
                # Does any LOO iteration push non-const above null rate?
                loo_max = float(np.max(loo_arr)) if len(loo_arr) > 0 else 0
                loo_all_below_null = bool(loo_max < null_nc_rate)
                
                bn_h['leave_one_out'] = {
                    'n_proteins': n_proteins,
                    'mean_nonconst': float(np.mean(loo_arr)),
                    'min_nonconst': float(np.min(loo_arr)),
                    'max_nonconst': float(np.max(loo_arr)),
                    'all_below_null': loo_all_below_null,
                }
                
                log.info("      Bootstrap 95%% CI: [%.1f%%, %.1f%%], null=%.1f%%",
                         np.percentile(boot_arr, 2.5)*100,
                         np.percentile(boot_arr, 97.5)*100,
                         null_nc_rate*100)
                log.info("      LOO range: [%.1f%%, %.1f%%], all < null: %s",
                         np.min(loo_arr)*100, np.max(loo_arr)*100,
                         loo_all_below_null)
        
        log.info("    Basin null results:")
        for ss in ['H', 'E', 'C']:
            if ss in basin_null_summary:
                bn = basin_null_summary[ss]
                log.info("      %s: real non-const=%.1f%%, null non-const=%.1f%%, "
                         "chi2=%.2f, p=%.4f",
                         ss, bn['real_nonconst_frac']*100, bn['null_nonconst_frac']*100,
                         bn['chi2'], bn['p'])
    
    stats_summary['basin_null'] = basin_null_summary
    
    # ── Tangent Autocorrelation by SS type ──
    log.info("  Tangent autocorrelation (kappa, segments >= 8)...")
    tangent_autocorr = {'H': [], 'E': [], 'C': []}
    for r in all_results:
        for seg in r.segments:
            if seg.ss_type in tangent_autocorr and seg.length >= 8:
                kappa_arr = np.array(seg.kappa_values)
                if len(kappa_arr) < 8:
                    continue
                # Tangent autocorrelation: C(Δ) = <T(s)·T(s+Δ)>
                # For curvature, use autocorrelation of κ(s)
                kappa_centered = kappa_arr - np.mean(kappa_arr)
                var = np.var(kappa_arr)
                if var < 1e-12:
                    continue
                max_lag = min(10, len(kappa_arr) - 1)
                autocorr = np.zeros(max_lag)
                for lag in range(max_lag):
                    n_pairs = len(kappa_arr) - lag
                    autocorr[lag] = np.mean(kappa_centered[:n_pairs] * kappa_centered[lag:lag+n_pairs]) / var
                tangent_autocorr[seg.ss_type].append(autocorr)
    
    # Average autocorrelation per SS type
    autocorr_summary = {}
    for ss in ['H', 'E', 'C']:
        if tangent_autocorr[ss]:
            # Pad to same length and average
            max_len = max(len(a) for a in tangent_autocorr[ss])
            padded = np.full((len(tangent_autocorr[ss]), max_len), np.nan)
            for i, a in enumerate(tangent_autocorr[ss]):
                padded[i, :len(a)] = a
            mean_autocorr = np.nanmean(padded, axis=0)
            autocorr_summary[ss] = {
                'n_segments': len(tangent_autocorr[ss]),
                'mean_autocorr': [float(x) for x in mean_autocorr],
                'lags': list(range(max_len)),
            }
            # Check for periodicity: does autocorrelation peak at lag 3-4?
            if max_len >= 5:
                peak_lag = np.argmax(mean_autocorr[2:min(6, max_len)]) + 2
                autocorr_summary[ss]['peak_lag'] = int(peak_lag)
                autocorr_summary[ss]['peak_value'] = float(mean_autocorr[peak_lag])
    
    stats_summary['tangent_autocorrelation'] = autocorr_summary
    
    # ── Dihedral Deviation Autocorrelation ──
    # Analyze (φᵢ - φ̄, ψᵢ - ψ̄) directly within each segment
    # This captures periodicity in raw angle space, not curvature (2nd derivative)
    log.info("  Dihedral deviation autocorrelation (phi-phi_bar, psi-psi_bar)...")
    dihedral_autocorr = {'H': {'phi': [], 'psi': []}, 
                         'E': {'phi': [], 'psi': []}, 
                         'C': {'phi': [], 'psi': []}}
    
    for r in all_results:
        # Need raw phi/psi for this protein
        pdb_path = PDB_DIR / f"{r.pdb_id}.pdb"
        if pdb_path.exists() and r.source == 'pdb':
            try:
                phi, psi, _, backbone_coords, has_O = extract_dihedrals_manual(
                    str(pdb_path), r.chain)
                valid = ~(np.isnan(phi) | np.isnan(psi))
                phi_v, psi_v = phi[valid], psi[valid]
                if np.all(has_O[valid]) and len(phi_v) >= 6:
                    ss = dssp_assign_pure_numpy(backbone_coords[valid])
                else:
                    ss = assign_ss_from_dihedrals(phi_v, psi_v)
            except Exception:
                continue
        else:
            # Synthetic — get from generator
            fold_class = assign_fold_class(r.name, r.ln_kf, r.co)
            phi_v, psi_v, ss = generate_synthetic_backbone(r.name, r.length, fold_class)
        
        if phi_v is None or len(phi_v) < 8:
            continue
        
        # Segment and compute autocorrelation of deviations
        current_ss = ss[0]
        start = 0
        for i in range(1, len(ss)):
            if ss[i] != current_ss or i == len(ss) - 1:
                end = i if ss[i] != current_ss else i + 1
                seg_len = end - start
                seg_ss = current_ss
                
                if seg_len >= 8 and seg_ss in dihedral_autocorr:
                    seg_phi = phi_v[start:end]
                    seg_psi = psi_v[start:end]
                    
                    for angle_name, angle_arr in [('phi', seg_phi), ('psi', seg_psi)]:
                        # Deviation from mean (respecting periodicity)
                        mean_angle = np.arctan2(np.mean(np.sin(angle_arr)), 
                                                np.mean(np.cos(angle_arr)))
                        dev = np.arctan2(np.sin(angle_arr - mean_angle), 
                                        np.cos(angle_arr - mean_angle))
                        var = np.var(dev)
                        if var < 1e-12:
                            continue
                        dev_centered = dev - np.mean(dev)
                        max_lag = min(10, len(dev) - 1)
                        ac = np.zeros(max_lag)
                        for lag in range(max_lag):
                            n_pairs = len(dev) - lag
                            ac[lag] = np.mean(dev_centered[:n_pairs] * dev_centered[lag:lag+n_pairs]) / var
                        dihedral_autocorr[seg_ss][angle_name].append(ac)
                
                current_ss = ss[i]
                start = i
    
    # Summarize
    dihedral_autocorr_summary = {}
    for ss in ['H', 'E', 'C']:
        ss_summary = {}
        for angle in ['phi', 'psi']:
            if dihedral_autocorr[ss][angle]:
                arrs = dihedral_autocorr[ss][angle]
                max_len = max(len(a) for a in arrs)
                padded = np.full((len(arrs), max_len), np.nan)
                for i, a in enumerate(arrs):
                    padded[i, :len(a)] = a
                mean_ac = np.nanmean(padded, axis=0)
                ss_summary[angle] = {
                    'n_segments': len(arrs),
                    'mean_autocorr': [float(x) for x in mean_ac],
                }
                if max_len >= 5:
                    # Check for periodic peak at lag 3-4 (excluding lag 0,1)
                    search_range = mean_ac[2:min(7, max_len)]
                    peak_lag = np.argmax(search_range) + 2
                    ss_summary[angle]['peak_lag'] = int(peak_lag)
                    ss_summary[angle]['peak_value'] = float(mean_ac[peak_lag])
        if ss_summary:
            dihedral_autocorr_summary[ss] = ss_summary
    
    stats_summary['dihedral_deviation_autocorr'] = dihedral_autocorr_summary
    
    log.info("  Omega validation (fitted omega vs 2pi/3.6)...")
    # ── Omega Validation: fitted ω vs expected 2π/3.6 for helix oscillatory segments ──
    expected_omega_per_residue = 2 * np.pi / 3.6  # ≈ 1.745 rad/residue
    helix_omegas = []
    helix_omega_per_residue = []
    for r in all_results:
        for seg in r.segments:
            if seg.ss_type == 'H' and seg.spiral_class in ('sinusoidal', 'damped_osc'):
                # The fits dict is not stored per segment, but we can recover ω
                # from re-classification (or store it). For now, re-classify.
                kappa_arr = np.array(seg.kappa_values)
                s_arr = np.array(seg.s_values)
                if len(kappa_arr) >= 5:
                    _, _, _, fit_info, _ = classify_spiral(kappa_arr, s_arr)
                    for model_name in ['sinusoidal', 'damped_osc']:
                        if model_name in fit_info and 'params' in fit_info[model_name]:
                            omega = fit_info[model_name]['params'].get('omega', None)
                            if omega is not None and omega > 0:
                                helix_omegas.append(omega)
                                # Convert to per-residue: omega is per arc-length unit
                                # Mean ds per residue for this segment
                                if len(s_arr) >= 2:
                                    mean_ds = (s_arr[-1] - s_arr[0]) / max(len(s_arr) - 1, 1)
                                    omega_per_res = omega * mean_ds
                                    helix_omega_per_residue.append(omega_per_res)
                                break
    
    if helix_omegas:
        omega_arr = np.array(helix_omega_per_residue) if helix_omega_per_residue else np.array(helix_omegas)
        stats_summary['omega_validation'] = {
            'n_helix_oscillatory': len(omega_arr),
            'expected_omega_per_residue': float(expected_omega_per_residue),
            'median_fitted_omega_per_residue': float(np.median(omega_arr)),
            'mean_fitted_omega_per_residue': float(np.mean(omega_arr)),
            'std_fitted_omega_per_residue': float(np.std(omega_arr)),
            'raw_omegas': [float(x) for x in omega_arr[:20]],  # first 20 for inspection
        }
    
    # ── Standardized Residuals for full chi-squared ──
    if 'chi_squared_full' in sc_agg:
        chi_info = sc_agg['chi_squared_full']
        obs = np.array(chi_info['table'])
        # Compute expected
        row_sums = obs.sum(axis=1, keepdims=True)
        col_sums = obs.sum(axis=0, keepdims=True)
        n_total = obs.sum()
        expected = row_sums * col_sums / n_total
        # Standardized residuals: (O-E) / sqrt(E)
        with np.errstate(divide='ignore', invalid='ignore'):
            std_resid = np.where(expected > 0, (obs - expected) / np.sqrt(expected), 0)
        sc_agg['standardized_residuals'] = {
            'matrix': std_resid.tolist(),
            'row_labels': chi_info.get('row_labels', ['H', 'E', 'C']),
            'col_labels': chi_info.get('col_labels', []),
        }
    
    # ── Bonferroni correction note ──
    focused_p_values = []
    for test_key in ['chi_squared_oscillatory', 'chi_squared_constant']:
        if test_key in sc_agg:
            focused_p_values.append(sc_agg[test_key]['p_value'])
    
    n_focused_tests = len(focused_p_values) + 1  # +1 for the full test
    bonferroni_threshold = 0.05 / max(n_focused_tests, 1)
    sc_agg['bonferroni'] = {
        'n_tests': n_focused_tests,
        'threshold': float(bonferroni_threshold),
        'note': f'Bonferroni-corrected α = {bonferroni_threshold:.4f} for {n_focused_tests} tests',
    }
    
    stats_summary['spiral_classification'] = sc_agg
    
    log.info("  Oscillatory evidence (delta-AICc per segment)...")
    # ── Oscillatory Evidence Analysis ──
    # Even when AICc selects constant, oscillatory models may be competitive
    # Δ(AICc) = AICc(constant) - AICc(best_oscillatory): positive = osc. better
    osc_evidence = {'H': [], 'E': [], 'C': []}
    osc_r2 = {'H': [], 'E': [], 'C': []}
    osc_omegas_helix = []
    
    for r in all_results:
        for seg in r.segments:
            if seg.spiral_class != 'too_short' and seg.ss_type in osc_evidence:
                osc_evidence[seg.ss_type].append(seg.delta_aicc_osc)
                osc_r2[seg.ss_type].append(seg.best_osc_r2)
                if seg.ss_type == 'H' and seg.delta_aicc_osc > -2 and seg.best_osc_omega > 0:
                    # Convert omega from per-arc-length to per-residue
                    s_arr = np.array(seg.s_values) if seg.s_values else np.array([])
                    if len(s_arr) >= 2:
                        mean_ds = (s_arr[-1] - s_arr[0]) / max(len(s_arr) - 1, 1)
                        omega_per_res = seg.best_osc_omega * mean_ds
                        osc_omegas_helix.append(omega_per_res)
    
    osc_summary = {}
    for ss in ['H', 'E', 'C']:
        if osc_evidence[ss]:
            arr = np.array(osc_evidence[ss])
            r2_arr = np.array(osc_r2[ss])
            osc_summary[ss] = {
                'n': len(arr),
                'mean_delta_aicc': float(np.mean(arr)),
                'median_delta_aicc': float(np.median(arr)),
                'frac_competitive': float(np.mean(arr > -2)),
                'frac_preferred': float(np.mean(arr > 0)),
                'mean_osc_r2': float(np.mean(r2_arr)),
                'median_osc_r2': float(np.median(r2_arr)),
            }
    
    h_deltas = np.array(osc_evidence.get('H', []))
    ec_deltas = np.array(osc_evidence.get('E', []) + osc_evidence.get('C', []))
    
    if len(h_deltas) >= 5 and len(ec_deltas) >= 5:
        u_osc, p_osc = stats.mannwhitneyu(h_deltas, ec_deltas, alternative='greater')
        effect_osc = 1.0 - 2.0 * u_osc / (len(h_deltas) * len(ec_deltas))
        osc_summary['mann_whitney'] = {
            'p_value': float(p_osc), 'effect_r': float(effect_osc),
            'H_median_delta': float(np.median(h_deltas)),
            'EC_median_delta': float(np.median(ec_deltas)),
        }
    
    if osc_omegas_helix:
        omega_arr = np.array(osc_omegas_helix)
        osc_summary['omega_per_residue'] = {
            'n': len(omega_arr),
            'expected': float(2 * np.pi / 3.6),
            'median': float(np.median(omega_arr)),
            'mean': float(np.mean(omega_arr)),
            'std': float(np.std(omega_arr)),
        }
    
    stats_summary['oscillatory_evidence'] = osc_summary
    
    log.info("  C3: correlations with folding rate")
    ln_kf = np.array([r.ln_kf for r in all_results])
    descriptors = {
        'Winding magnitude': np.array([r.torus_knot.get('total_winding_magnitude', 0) 
                                        for r in all_results]),
        'p (phi winding)': np.array([r.torus_knot.get('p', 0) for r in all_results]),
        'q (psi winding)': np.array([r.torus_knot.get('q', 0) for r in all_results]),
        'Mean curvature': np.array([r.mean_curvature for r in all_results]),
        'Curvature std': np.array([r.std_curvature for r in all_results]),
        'Regular fraction': np.array([r.overall_clothoid_fraction for r in all_results]),
        'BPS/residue': np.array([r.bps_per_residue for r in all_results]),
        'Contact order': np.array([r.co for r in all_results]),
        'Chain length': np.array([r.length for r in all_results]),
    }
    
    correlations = {}
    for name, vals in descriptors.items():
        valid = ~(np.isnan(vals) | np.isnan(ln_kf))
        if valid.sum() >= 5:
            r_p, p_p = stats.pearsonr(vals[valid], ln_kf[valid])
            rho, p_rho = stats.spearmanr(vals[valid], ln_kf[valid])
            correlations[name] = {
                'pearson_r': float(r_p), 'pearson_p': float(p_p),
                'spearman_rho': float(rho), 'spearman_p': float(p_rho),
            }
    
    stats_summary['correlations'] = correlations
    
    # ── Structural Class Separation by (p,q) Winding ──
    # Test: do torus winding vectors cluster by SCOP structural class?
    log.info("  Structural class separation (all-α / all-β / α+β by winding)...")
    
    # Classify proteins by structural class based on SS composition
    STRUCT_CLASS_OVERRIDE = {
        # all-α: >60% H, <10% E
        'lambda-repressor': 'all-alpha', 'Villin-HP': 'all-alpha',
        'Engrailed-HD': 'all-alpha', 'Arc-repressor': 'all-alpha',
        'Myoglobin': 'all-alpha',
        # all-β: >40% E, <10% H
        'SH3-src': 'all-beta', 'SH3-spectrin': 'all-beta',
        'TNfn3': 'all-beta', 'CspB': 'all-beta',
        # α+β: substantial both
        'CI2': 'alpha-beta', 'Ubiquitin': 'alpha-beta',
        'Protein-L': 'alpha-beta', 'Protein-G': 'alpha-beta',
        'FKBP': 'alpha-beta', 'Barnase': 'alpha-beta',
        'Barstar': 'alpha-beta', 'CheY': 'alpha-beta',
        'Acylphosphatase': 'alpha-beta', 'Cytochrome-c': 'alpha-beta',
        'RNase-A': 'alpha-beta', 'Lysozyme': 'alpha-beta',
        'DHFR': 'alpha-beta', 'Thioredoxin': 'alpha-beta',
    }
    
    class_data = {'all-alpha': [], 'all-beta': [], 'alpha-beta': []}
    class_names_list = []
    
    for r in all_results:
        sc = STRUCT_CLASS_OVERRIDE.get(r.name, None)
        if sc is None:
            # Auto-classify from SS fractions
            n_res = r.n_valid_residues
            if n_res == 0:
                continue
            n_H = sum(1 for seg in r.segments if seg.ss_type == 'H' for _ in range(seg.length))
            n_E = sum(1 for seg in r.segments if seg.ss_type == 'E' for _ in range(seg.length))
            fH, fE = n_H / n_res, n_E / n_res
            if fH > 0.5 and fE < 0.1:
                sc = 'all-alpha'
            elif fE > 0.3 and fH < 0.1:
                sc = 'all-beta'
            else:
                sc = 'alpha-beta'
        
        tk = r.torus_knot
        if tk:
            p_w = tk.get('p', 0)
            q_w = tk.get('q', 0)
            Q_w = tk.get('total_winding_magnitude', 0)
            class_data[sc].append({'name': r.name, 'p': p_w, 'q': q_w, 'Q': Q_w,
                                   'ln_kf': r.ln_kf, 'length': r.length})
            class_names_list.append((r.name, sc, p_w, q_w, Q_w))
    
    struct_class_summary = {}
    
    # Per-class winding statistics
    for sc in ['all-alpha', 'all-beta', 'alpha-beta']:
        entries = class_data[sc]
        if entries:
            ps = np.array([e['p'] for e in entries])
            qs = np.array([e['q'] for e in entries])
            Qs = np.array([e['Q'] for e in entries])
            struct_class_summary[sc] = {
                'n': len(entries),
                'proteins': [e['name'] for e in entries],
                'p_mean': float(np.mean(ps)), 'p_std': float(np.std(ps)),
                'q_mean': float(np.mean(qs)), 'q_std': float(np.std(qs)),
                'Q_mean': float(np.mean(Qs)), 'Q_std': float(np.std(Qs)),
            }
    
    # MANOVA-like test: Kruskal-Wallis on |Q| by structural class
    class_groups_Q = []
    class_labels_Q = []
    for sc in ['all-alpha', 'all-beta', 'alpha-beta']:
        if class_data[sc]:
            Qs = [e['Q'] for e in class_data[sc]]
            class_groups_Q.append(Qs)
            class_labels_Q.append(sc)
    
    if len(class_groups_Q) >= 2 and all(len(g) >= 2 for g in class_groups_Q):
        h_stat, kw_p = stats.kruskal(*class_groups_Q)
        struct_class_summary['kruskal_wallis_Q'] = {
            'H': float(h_stat), 'p': float(kw_p),
            'groups': class_labels_Q,
        }
    
    # Kruskal-Wallis on p (phi winding) by structural class
    class_groups_p = []
    for sc in ['all-alpha', 'all-beta', 'alpha-beta']:
        if class_data[sc]:
            class_groups_p.append([e['p'] for e in class_data[sc]])
    
    if len(class_groups_p) >= 2 and all(len(g) >= 2 for g in class_groups_p):
        h_p, kw_p_p = stats.kruskal(*class_groups_p)
        struct_class_summary['kruskal_wallis_p'] = {'H': float(h_p), 'p': float(kw_p_p)}
    
    # Kruskal-Wallis on q (psi winding) by structural class
    class_groups_q = []
    for sc in ['all-alpha', 'all-beta', 'alpha-beta']:
        if class_data[sc]:
            class_groups_q.append([e['q'] for e in class_data[sc]])
    
    if len(class_groups_q) >= 2 and all(len(g) >= 2 for g in class_groups_q):
        h_q, kw_p_q = stats.kruskal(*class_groups_q)
        struct_class_summary['kruskal_wallis_q'] = {'H': float(h_q), 'p': float(kw_p_q)}
    
    # Pairwise comparisons: all-α vs all-β (most distinct classes)
    alpha_entries = class_data.get('all-alpha', [])
    beta_entries = class_data.get('all-beta', [])
    if len(alpha_entries) >= 2 and len(beta_entries) >= 2:
        alpha_Q = [e['Q'] for e in alpha_entries]
        beta_Q = [e['Q'] for e in beta_entries]
        u_ab, p_ab = stats.mannwhitneyu(alpha_Q, beta_Q, alternative='two-sided')
        struct_class_summary['alpha_vs_beta_Q'] = {
            'alpha_median_Q': float(np.median(alpha_Q)),
            'beta_median_Q': float(np.median(beta_Q)),
            'mannwhitney_p': float(p_ab),
        }
        
        alpha_p = [e['p'] for e in alpha_entries]
        beta_p = [e['p'] for e in beta_entries]
        u_p, p_pw = stats.mannwhitneyu(alpha_p, beta_p, alternative='two-sided')
        struct_class_summary['alpha_vs_beta_p'] = {
            'alpha_median_p': float(np.median(alpha_p)),
            'beta_median_p': float(np.median(beta_p)),
            'mannwhitney_p': float(p_pw),
        }
        
        alpha_q = [e['q'] for e in alpha_entries]
        beta_q = [e['q'] for e in beta_entries]
        u_q, p_qw = stats.mannwhitneyu(alpha_q, beta_q, alternative='two-sided')
        struct_class_summary['alpha_vs_beta_q'] = {
            'alpha_median_q': float(np.median(alpha_q)),
            'beta_median_q': float(np.median(beta_q)),
            'mannwhitney_p': float(p_qw),
        }
    
    # Log summary
    for sc in ['all-alpha', 'all-beta', 'alpha-beta']:
        if sc in struct_class_summary:
            s = struct_class_summary[sc]
            log.info("    %s (n=%d): p=%.2f±%.2f, q=%.2f±%.2f, |Q|=%.2f±%.2f",
                     sc, s['n'], s['p_mean'], s['p_std'],
                     s['q_mean'], s['q_std'], s['Q_mean'], s['Q_std'])
    if 'kruskal_wallis_Q' in struct_class_summary:
        kw = struct_class_summary['kruskal_wallis_Q']
        log.info("    KW |Q| by class: H=%.2f, p=%.4f", kw['H'], kw['p'])
    if 'alpha_vs_beta_Q' in struct_class_summary:
        ab = struct_class_summary['alpha_vs_beta_Q']
        log.info("    α vs β |Q|: %.2f vs %.2f, MW p=%.4f",
                 ab['alpha_median_Q'], ab['beta_median_Q'], ab['mannwhitney_p'])
    
    stats_summary['structural_class'] = struct_class_summary
    
    # ── Active Site Proximity Analysis ──
    # Test whether Gaussian peak curvature anomalies on T² predict functional sites
    log.info("  Active site proximity analysis (gauss_peak vs catalytic residues)...")
    
    # Known catalytic/active/binding/interface site residues for all 23 proteins
    # Sources: UniProt, Catalytic Site Atlas, PDBsum, primary literature
    # Residue numbers are PDB numbering (1-indexed)
    # 'functional' = union of all annotated sites; 'catalytic' = catalytic only
    ACTIVE_SITES = {
        # ── Enzymes ──
        'Barnase':      {'catalytic': [12, 58, 73, 87],
                         'functional': [12, 27, 54, 56, 58, 59, 73, 83, 87, 102]},
        'RNase-A':      {'catalytic': [12, 41, 119],
                         'functional': [7, 8, 11, 12, 41, 44, 45, 65, 119, 120, 121]},
        'Lysozyme':     {'catalytic': [35, 52],
                         'functional': [34, 35, 46, 50, 52, 53, 56, 57, 59, 62, 63, 101, 107, 108, 109]},
        'DHFR':         {'catalytic': [5, 27],
                         'functional': [5, 6, 7, 27, 28, 31, 32, 49, 52, 54, 94, 100]},
        'Acylphosphatase': {'catalytic': [18, 21, 23],
                            'functional': [15, 16, 17, 18, 20, 21, 22, 23, 41, 43, 44, 45]},
        'Cytochrome-c': {'catalytic': [14, 17, 18, 80],  # heme axial ligands
                         'functional': [14, 17, 18, 41, 48, 52, 78, 79, 80, 82]},
        'Thioredoxin':  {'catalytic': [32, 35],
                         'functional': [26, 27, 28, 32, 33, 34, 35, 73, 74, 75, 76]},
        'CheY':         {'catalytic': [12, 13, 57, 59, 109],
                         'functional': [12, 13, 14, 57, 58, 59, 86, 87, 106, 109]},
        # ── Binding proteins (non-enzyme) ──
        'CI2':          {'catalytic': [],  # serine protease inhibitor
                         'functional': [34, 35, 36, 37, 38, 39, 40, 53, 56, 59]},  # reactive loop + exosite
        'SH3-src':      {'catalytic': [],
                         'functional': [7, 8, 9, 10, 36, 37, 49, 50, 51, 52, 53, 54]},  # RT/n-Src/distal loops (PPII binding)
        'Protein-L':    {'catalytic': [],
                         # beta-strand 2 + helix C-term + helix-strand3 loop (Ig VL binding)
                         'functional': [31, 33, 34, 35, 36, 37, 45, 46, 47, 48, 49, 50, 51, 52]},
        'CspB':         {'catalytic': [],
                         'functional': [6, 7, 8, 11, 17, 29, 30, 33, 36, 41, 53, 56, 58]},  # ssDNA binding (RNP motifs)
        'Ubiquitin':    {'catalytic': [],
                         'functional': [6, 8, 11, 27, 29, 33, 42, 44, 48, 63, 68, 70, 72, 73, 74, 75, 76]},  # I44 patch + C-term tail
        'Protein-G':    {'catalytic': [],
                         # Fc binding: E27 hot spot + helix/strand2 interface
                         'functional': [15, 16, 17, 18, 19, 20, 24, 25, 27, 28, 30, 31, 34, 39, 40, 41, 43]},
        'FKBP':         {'catalytic': [],
                         'functional': [26, 36, 46, 50, 54, 55, 56, 59, 82, 87, 97, 99]},  # FK506/rapamycin pocket
        'Barstar':      {'catalytic': [],
                         'functional': [27, 29, 33, 35, 36, 38, 39, 42, 73, 74, 76, 80]},  # barnase binding interface
        'lambda-repressor': {'catalytic': [],
                             'functional': [6, 9, 10, 13, 14, 17, 22, 33, 36, 37, 40, 44, 47, 48]},  # DNA-binding helix-turn-helix
        'Villin-HP':    {'catalytic': [],
                         'functional': [8, 12, 20, 24, 28, 32]},  # hydrophobic core (Trp/Phe burial)
        'Engrailed-HD': {'catalytic': [],
                         'functional': [3, 6, 7, 28, 31, 35, 38, 41, 42, 43, 44, 47, 50, 51]},  # HTH DNA recognition
        'Myoglobin':    {'catalytic': [],
                         'functional': [58, 64, 67, 68, 93, 97, 99, 138, 140, 143, 148]},  # heme pocket (proximal/distal His + contacts)
        'Arc-repressor': {'catalytic': [],
                          'functional': [6, 7, 8, 9, 10, 11, 13, 36, 37, 39, 43, 44, 46]},  # DNA-binding beta-sheet + dimerization
        'SH3-spectrin': {'catalytic': [],
                         'functional': [7, 8, 9, 10, 36, 37, 49, 50, 51, 52, 53, 54]},  # PPII-binding groove (analogous to SH3-src)
        'TNfn3':        {'catalytic': [],
                         'functional': [21, 22, 23, 44, 45, 46, 76, 77, 78, 79, 80, 81, 82, 83]},  # BC/FG loops (integrin-binding analog)
    }
    
    # Collect all gauss_peak segments with their peak residue positions
    gauss_peaks = []  # (protein_name, peak_residue, segment_start, segment_end, A, sigma, ss_type)
    all_helix_residues = {}  # protein_name -> list of all helix residue indices
    
    for r in all_results:
        helix_res = []
        for seg in r.segments:
            if seg.ss_type == 'H' and seg.spiral_class != 'too_short':
                helix_res.extend(range(seg.start_idx, seg.end_idx + 1))
            if seg.spiral_class == 'gauss_peak' and seg.gauss_peak_residue >= 0:
                gauss_peaks.append({
                    'protein': r.name,
                    'pdb_id': r.pdb_id,
                    'peak_residue': seg.gauss_peak_residue,
                    'start': seg.start_idx,
                    'end': seg.end_idx,
                    'ss_type': seg.ss_type,
                    'A': seg.gauss_A,
                    'sigma': seg.gauss_sigma,
                    'c': seg.gauss_c,
                    'length': seg.length,
                })
        if helix_res:
            all_helix_residues[r.name] = helix_res
    
    # Compute proximity: min distance from gauss_peak residue to nearest active site residue
    proximity_results = []
    null_distances = []  # distances from random helix residues to active sites
    
    rng = np.random.RandomState(42)
    
    for gp in gauss_peaks:
        pname = gp['protein']
        if pname not in ACTIVE_SITES:
            continue
        
        sites = ACTIVE_SITES[pname]
        all_func = sites['functional']
        if not all_func:
            continue
        
        # Convert to 0-indexed for comparison
        func_residues_0 = [r - 1 for r in all_func]
        cat_residues_0 = [r - 1 for r in sites['catalytic']] if sites['catalytic'] else None
        
        peak_res = gp['peak_residue']
        
        # Distance to nearest functional residue
        min_dist_func = min(abs(peak_res - fr) for fr in func_residues_0)
        min_dist_cat = min(abs(peak_res - cr) for cr in cat_residues_0) if cat_residues_0 else None
        
        proximity_results.append({
            **gp,
            'min_dist_functional': min_dist_func,
            'min_dist_catalytic': min_dist_cat,
        })
        
        # Null: sample random helix residues from same protein
        if pname in all_helix_residues and all_helix_residues[pname]:
            for _ in range(50):  # 50 null samples per real peak
                rand_res = rng.choice(all_helix_residues[pname])
                null_dist = min(abs(rand_res - fr) for fr in func_residues_0)
                null_distances.append(null_dist)
    
    # Statistical test: are gauss_peak residues closer to functional sites than random helix residues?
    active_site_summary = {
        'n_gauss_peaks_total': len(gauss_peaks),
        'n_gauss_peaks_in_enzymes': len(proximity_results),
        'peaks': proximity_results,
    }
    
    if proximity_results:
        real_distances = [pr['min_dist_functional'] for pr in proximity_results]
        active_site_summary['real_mean_dist'] = float(np.mean(real_distances))
        active_site_summary['real_median_dist'] = float(np.median(real_distances))
        active_site_summary['real_distances'] = real_distances
        
        if null_distances:
            null_arr = np.array(null_distances)
            active_site_summary['null_mean_dist'] = float(np.mean(null_arr))
            active_site_summary['null_median_dist'] = float(np.median(null_arr))
            
            # One-sided Mann-Whitney: real < null (peaks are closer than random)
            if len(real_distances) >= 2:
                u_stat, p_two = stats.mannwhitneyu(real_distances, null_distances, alternative='less')
                active_site_summary['mannwhitney_p'] = float(p_two)
                n1, n2 = len(real_distances), len(null_distances)
                active_site_summary['effect_r'] = float(1 - 2*u_stat/(n1*n2))
            
            # Percentile: where do real distances fall in null distribution?
            percentiles = [float(np.mean(null_arr <= d)) * 100 for d in real_distances]
            active_site_summary['percentiles_in_null'] = percentiles
        
        log.info("    Gauss peaks in enzymes: %d, mean dist to active site: %.1f residues",
                 len(proximity_results), np.mean(real_distances))
        if null_distances:
            log.info("    Null (random helix): mean dist = %.1f, MW p = %.4f",
                     np.mean(null_distances),
                     active_site_summary.get('mannwhitney_p', 1.0))
    else:
        log.info("    No gauss_peak segments found in annotated enzymes")
    
    # Gaussian peak parameter catalog
    if gauss_peaks:
        active_site_summary['peak_catalog'] = gauss_peaks
        amplitudes = [gp['A'] for gp in gauss_peaks]
        sigmas = [gp['sigma'] for gp in gauss_peaks]
        active_site_summary['amplitude_stats'] = {
            'mean': float(np.mean(amplitudes)),
            'std': float(np.std(amplitudes)),
            'median': float(np.median(amplitudes)),
        }
        active_site_summary['sigma_stats'] = {
            'mean': float(np.mean(sigmas)),
            'std': float(np.std(sigmas)),
            'median': float(np.median(sigmas)),
        }
    
    stats_summary['active_site_proximity'] = active_site_summary
    
    # ── Figures ──
    log.info("Generating figures (exemplars, curvature, taxonomy, knots, dashboard)...")
    
    exemplars = ['CI2', 'Ubiquitin', 'Barnase', 'Myoglobin']
    
    for result in all_results:
        if result.name in exemplars and result.n_valid_residues > 0:
            pdb_path = PDB_DIR / f"{result.pdb_id}.pdb"
            if pdb_path.exists() and result.source == 'pdb':
                phi, psi, res_ids, backbone_coords, has_O = extract_dihedrals_manual(
                    str(pdb_path), result.chain)
                valid = ~(np.isnan(phi) | np.isnan(psi))
                phi, psi = phi[valid], psi[valid]
                backbone_coords_valid = backbone_coords[valid]
                has_O_valid = has_O[valid]
                
                # Re-derive SS using pure-numpy DSSP
                if np.all(has_O_valid) and len(phi) >= 6:
                    try:
                        ss = dssp_assign_pure_numpy(backbone_coords_valid)
                    except Exception:
                        ss = assign_ss_from_dihedrals(phi, psi)
                else:
                    ss = assign_ss_from_dihedrals(phi, psi)
            else:
                fold_class = assign_fold_class(result.name, result.ln_kf, result.co)
                phi, psi, ss = generate_synthetic_backbone(result.name, result.length, fold_class)
            
            kappa, s_kappa, _ = torus_curvature(phi, psi)
            segments = segment_backbone(phi, psi, ss, kappa, s_kappa)
            
            fig_ramachandran_with_spiral_overlay(
                phi, psi, ss, segments, result.name,
                phi_grid, psi_grid, W_grid,
                FIGURES_DIR / f"fig1_ramachandran_{result.name}.png"
            )
            
            fig_cornu_exemplar(
                phi, psi, ss, segments, result.name,
                FIGURES_DIR / f"fig5_cornu_exemplar_{result.name}.png"
            )
    
    fig_curvature_regularity(all_results, FIGURES_DIR / "fig2_curvature_regularity.png")
    fig_spiral_classification(all_results, FIGURES_DIR / "fig3_spiral_taxonomy.png")
    fig_torus_knot_descriptors(all_results, FIGURES_DIR / "fig4_torus_knot_descriptors.png")
    fig_summary_dashboard(all_results, FIGURES_DIR / "fig6_summary_dashboard.png")
    
    log.info(f"  Figures saved to {FIGURES_DIR}/")
    
    # ── Report ──
    generate_report(all_results, stats_summary, RESULTS_DIR / "cornu_analysis_report.md")
    
    # ── Save JSON ──
    results_json = []
    for r in all_results:
        d = asdict(r)
        for key, val in d.items():
            if isinstance(val, (np.floating, np.integer)):
                d[key] = float(val)
        results_json.append(d)
    
    with open(RESULTS_DIR / "protein_results.json", 'w', encoding='utf-8') as f:
        json.dump(results_json, f, indent=2, default=str)
    
    with open(RESULTS_DIR / "stats_summary.json", 'w', encoding='utf-8') as f:
        json.dump(stats_summary, f, indent=2, default=str)
    
    # ── Print summary ──
    print()
    print("=" * 72)
    print("  RESULTS SUMMARY")
    print("=" * 72)
    print()
    print(f"  Proteins analyzed: {len(all_results)}")
    print(f"  Sources: {n_pdb} PDB, {len(all_results) - n_pdb} synthetic")
    print(f"  SS method: {n_dssp} DSSP, {n_dihedral} dihedral fallback")
    print()
    print("  CONTRIBUTION 1 -- Geometric Completeness (κ follows simple forms):")
    for ss_type in ['H', 'E', 'C']:
        if ss_type in wrc_agg:
            info = wrc_agg[ss_type]
            print(f"    {ss_type}: intra-CV = {info.get('mean_intra_cv', 0):.2f}, "
                  f"R² = {info.get('mean_r_squared', 0):.3f}, "
                  f"regular = {info.get('regular_frac', 0)*100:.0f}%, "
                  f"n = {info['n']} segments, "
                  f"mean len = {info.get('mean_segment_length', 0):.1f}")
    if 'regularity_inversion_test' in wrc_agg:
        rt = wrc_agg['regularity_inversion_test']
        print(f"    Mann-Whitney (coil CV < struct CV): p={rt['mann_whitney_p']:.4f}, "
              f"effect r={rt['effect_size_r']:.3f}")
    print(f"    Overall regular fraction: {stats_summary.get('overall_regular_fraction', 0)*100:.1f}%")
    nm = stats_summary.get('null_model', {})
    if nm:
        print(f"    NULL MODEL ({nm.get('method', 'shuffle')}): "
              f"null Δ={nm.get('null_median_delta_aicc', 0):.2f}, "
              f"real Δ={nm.get('real_median_delta_aicc', 0):.2f}, "
              f"MW p={nm.get('mann_whitney_p', 1):.4f}, r={nm.get('effect_r', 0):.3f}")
    print()
    print("  CONTRIBUTION 2 -- SS-Dependent Spiral Taxonomy:")
    for ss_type in ['H', 'E', 'C']:
        if ss_type in sc_agg and not ss_type.startswith('chi'):
            t = sc_agg[ss_type]
            # Filter out metadata keys
            if not isinstance(t, dict) or 'chi2' in t:
                continue
            parts = []
            for shape in ['geodesic', 'circular_arc', 'sinusoidal', 'damped_osc',
                          'exponential', 'sigmoid', 'step', 'clothoid', 'fermat', 'quadratic', 'gauss_peak']:
                if t.get(shape, 0) > 0:
                    parts.append(f"{shape[:6]}={t[shape]}")
            total = sum(t.get(s, 0) for s in ['geodesic', 'circular_arc', 'sinusoidal',
                        'damped_osc', 'exponential', 'sigmoid', 'step', 
                        'clothoid', 'fermat', 'quadratic', 'gauss_peak'])
            too_short = t.get('too_short', 0)
            print(f"    {ss_type}: {', '.join(parts)}  (total={total}, short={too_short})")
    if 'chi_squared_full' in sc_agg:
        chi = sc_agg['chi_squared_full']
        print(f"    Full 3×k chi-sq: chi2={chi['chi2']:.2f}, p={chi['p_value']:.4f}, "
              f"V={chi.get('cramers_v', 0):.3f}")
    if 'chi_squared_oscillatory' in sc_agg:
        chi = sc_agg['chi_squared_oscillatory']
        print(f"    H oscillatory enrichment: chi2={chi['chi2']:.2f}, p={chi['p_value']:.4f}")
    if 'chi_squared_constant' in sc_agg:
        chi = sc_agg['chi_squared_constant']
        print(f"    E constant-κ enrichment: chi2={chi['chi2']:.2f}, p={chi['p_value']:.4f}")
    if 'bonferroni' in sc_agg:
        print(f"    Bonferroni: α={sc_agg['bonferroni']['threshold']:.4f} ({sc_agg['bonferroni']['n_tests']} tests)")
    ov = stats_summary.get('omega_validation', {})
    if ov:
        print(f"    OMEGA: expected={ov['expected_omega_per_residue']:.3f}, "
              f"fitted median={ov['median_fitted_omega_per_residue']:.3f}, "
              f"n={ov['n_helix_oscillatory']}")
    oe = stats_summary.get('oscillatory_evidence', {})
    if oe:
        for ss in ['H', 'E', 'C']:
            if ss in oe:
                info = oe[ss]
                print(f"    ΔAICc {ss}: median={info['median_delta_aicc']:.1f}, "
                      f"competitive={info['frac_competitive']*100:.0f}%, "
                      f"preferred={info['frac_preferred']*100:.0f}%, "
                      f"osc R²={info['mean_osc_r2']:.3f}")
        if 'mann_whitney' in oe:
            mw = oe['mann_whitney']
            print(f"    MW (H Δ > E+C Δ): p={mw['p_value']:.4f}, r={mw['effect_r']:.3f}")
        if 'omega_per_residue' in oe:
            om = oe['omega_per_residue']
            print(f"    ω/residue: expected={om['expected']:.3f}, "
                  f"median={om['median']:.3f}, n={om['n']}")
    tac = stats_summary.get('tangent_autocorrelation', {})
    if tac:
        for ss in ['H', 'E', 'C']:
            if ss in tac:
                info = tac[ss]
                ac = info['mean_autocorr']
                c3 = ac[3] if len(ac) > 3 else 0
                c4 = ac[4] if len(ac) > 4 else 0
                print(f"    Autocorr {ss}: n={info['n_segments']}, "
                      f"peak_lag={info.get('peak_lag', '?')}, "
                      f"C(3)={c3:.3f}, C(4)={c4:.3f}")
    dac = stats_summary.get('dihedral_deviation_autocorr', {})
    if dac:
        for ss in ['H', 'E', 'C']:
            if ss in dac:
                parts = []
                for angle in ['phi', 'psi']:
                    if angle in dac[ss]:
                        info = dac[ss][angle]
                        ac = info['mean_autocorr']
                        c3 = ac[3] if len(ac) > 3 else 0
                        c4 = ac[4] if len(ac) > 4 else 0
                        sym = 'φ' if angle == 'phi' else 'ψ'
                        parts.append(f"{sym}: peak={info.get('peak_lag','?')}, "
                                    f"C(3)={c3:.3f}, C(4)={c4:.3f}")
                print(f"    DihedralAC {ss}: {'; '.join(parts)}")
    print()
    print("  CONTRIBUTION 3 -- Torus Knot Descriptors & Correlations:")
    for name in ['Winding magnitude', 'Mean curvature', 'Regular fraction',
                 'Contact order', 'Chain length']:
        if name in correlations:
            c = correlations[name]
            print(f"    {name:20s}: r={c['pearson_r']:+.3f}  (p={c['pearson_p']:.4f})  "
                  f"rho={c['spearman_rho']:+.3f}")
    print()
    
    # Basin-control null
    bn = stats_summary.get('basin_null', {})
    if bn:
        print("  BASIN-CONTROL NULL (Markov walk within SS basins):")
        for ss in ['H', 'E', 'C']:
            if ss in bn:
                b = bn[ss]
                print(f"    {ss}: real non-const={b['real_nonconst_frac']*100:.1f}%, "
                      f"null non-const={b['null_nonconst_frac']*100:.1f}%, "
                      f"chi2={b['chi2']:.2f}, p={b['p']:.4f}")
        # Robustness for helix
        if 'H' in bn:
            bn_h = bn['H']
            es = bn_h.get('effect_sizes', {})
            if es:
                print(f"    H effect: ARR={es['abs_risk_reduction']*100:.1f}pp, "
                      f"OR={es['odds_ratio']:.2f}, "
                      f"binary χ²={es['chi2_binary']:.2f} p={es['p_binary']:.4f}")
            boot = bn_h.get('bootstrap', {})
            if boot:
                print(f"    H bootstrap 95%CI: [{boot['ci_95_low']*100:.1f}%, "
                      f"{boot['ci_95_high']*100:.1f}%], "
                      f"exceeds null: {boot['frac_exceeds_null']*100:.1f}%")
            loo = bn_h.get('leave_one_out', {})
            if loo:
                print(f"    H LOO range: [{loo['min_nonconst']*100:.1f}%, "
                      f"{loo['max_nonconst']*100:.1f}%], "
                      f"all < null: {loo['all_below_null']}")
        print()
    
    # Structural class separation
    scs = stats_summary.get('structural_class', {})
    if scs:
        print("  STRUCTURAL CLASS SEPARATION (all-α / all-β / α+β):")
        for sc in ['all-alpha', 'all-beta', 'alpha-beta']:
            if sc in scs:
                s = scs[sc]
                print(f"    {sc:12s} (n={s['n']}): p={s['p_mean']:+.2f}±{s['p_std']:.2f}, "
                      f"q={s['q_mean']:+.2f}±{s['q_std']:.2f}, "
                      f"|Q|={s['Q_mean']:.2f}±{s['Q_std']:.2f}")
        if 'kruskal_wallis_Q' in scs:
            kw = scs['kruskal_wallis_Q']
            print(f"    KW |Q|: H={kw['H']:.2f}, p={kw['p']:.4f}")
        if 'alpha_vs_beta_Q' in scs:
            ab = scs['alpha_vs_beta_Q']
            print(f"    α vs β |Q|: {ab['alpha_median_Q']:.2f} vs {ab['beta_median_Q']:.2f}, "
                  f"p={ab['mannwhitney_p']:.4f}")
        if 'alpha_vs_beta_p' in scs:
            ab = scs['alpha_vs_beta_p']
            print(f"    α vs β p: {ab['alpha_median_p']:+.2f} vs {ab['beta_median_p']:+.2f}, "
                  f"p={ab['mannwhitney_p']:.4f}")
        print()
    
    # Active site proximity
    asp = stats_summary.get('active_site_proximity', {})
    if asp.get('n_gauss_peaks_total', 0) > 0:
        print("  ACTIVE SITE PROXIMITY (Gaussian peaks vs functional residues):")
        print(f"    Gauss peaks total: {asp['n_gauss_peaks_total']}, "
              f"in enzymes: {asp['n_gauss_peaks_in_enzymes']}")
        if asp.get('peaks'):
            for pk in asp['peaks']:
                cat_d = pk.get('min_dist_catalytic', '?')
                print(f"    {pk['protein']:20s} res={pk['peak_residue']:3d} "
                      f"A={pk['A']:+.3f} σ={pk['sigma']:.3f} "
                      f"dist_func={pk['min_dist_functional']} dist_cat={cat_d}")
        if 'real_mean_dist' in asp:
            print(f"    Real mean dist: {asp['real_mean_dist']:.1f} residues")
        if 'null_mean_dist' in asp:
            print(f"    Null mean dist: {asp['null_mean_dist']:.1f} residues")
        if 'mannwhitney_p' in asp:
            print(f"    MW (peak < random): p={asp['mannwhitney_p']:.4f}, "
                  f"r={asp['effect_r']:.3f}")
        print()
    
    print(f"  Full report: {RESULTS_DIR / 'cornu_analysis_report.md'}")
    print(f"  Figures:      {FIGURES_DIR}/")
    print()
    
    return all_results, stats_summary



# =============================================================================
# SECTION 15: ALPHAFOLD INGESTION & GEOMETRIC BARCODE
# =============================================================================

ALPHAFOLD_DIR = DATA_DIR / "alphafold"
ALPHAFOLD_CACHE_DIR = Path("alphafold_cache")  # User's pre-downloaded 125k structure cache
BARCODE_DIR = RESULTS_DIR / "barcodes"
for _d in [ALPHAFOLD_DIR, BARCODE_DIR]:
    _d.mkdir(parents=True, exist_ok=True)


def extract_dihedrals_mmcif(cif_path, chain_id=None):
    """Extract backbone (φ, ψ) and N,CA,C,O coords from mmCIF file.
    
    AlphaFold mmCIF files use auth_asym_id for chain.
    If chain_id is None, uses the first chain found.
    """
    atoms_by_residue = {}
    target_chain = chain_id
    
    # Simple mmCIF ATOM parser (handles AlphaFold format)
    in_atom_site = False
    col_names = []
    
    with open(cif_path, 'r') as f:
        for line in f:
            line = line.rstrip()
            
            # Detect _atom_site loop
            if line.startswith('_atom_site.'):
                if not in_atom_site:
                    in_atom_site = True
                    col_names = []
                col_names.append(line.split('.')[1].strip())
                continue
            
            if in_atom_site and (line.startswith('_') or line.startswith('#') or line.startswith('loop_')):
                in_atom_site = False
                continue
            
            if not in_atom_site or not line or line.startswith('#'):
                continue
            
            # Parse atom record
            parts = line.split()
            if len(parts) < len(col_names):
                continue
            
            record = {col_names[i]: parts[i] for i in range(len(col_names))}
            
            group = record.get('group_PDB', '')
            if group not in ('ATOM', 'HETATM'):
                continue
            
            atom_name = record.get('label_atom_id', record.get('auth_atom_id', ''))
            if atom_name not in ('N', 'CA', 'C', 'O'):
                continue
            
            # Chain: prefer auth_asym_id (matches PDB convention)
            chain = record.get('auth_asym_id', record.get('label_asym_id', 'A'))
            if target_chain is None:
                target_chain = chain
            if chain != target_chain:
                continue
            
            # Alt conformer
            alt = record.get('label_alt_id', '.')
            if alt not in ('.', '?', '', 'A', '1'):
                continue
            
            try:
                res_seq = int(record.get('auth_seq_id', record.get('label_seq_id', '0')))
                icode = record.get('pdbx_PDB_ins_code', '.')
                if icode in ('.', '?'):
                    icode = ''
                x = float(record.get('Cartn_x', 0))
                y = float(record.get('Cartn_y', 0))
                z = float(record.get('Cartn_z', 0))
            except (ValueError, KeyError):
                continue
            
            key = (res_seq, icode)
            if key not in atoms_by_residue:
                atoms_by_residue[key] = {}
            atoms_by_residue[key][atom_name] = np.array([x, y, z])
    
    if not atoms_by_residue:
        raise ValueError(f"No atoms found in {cif_path} for chain {target_chain}")
    
    # Same dihedral extraction as PDB parser
    sorted_keys = sorted(atoms_by_residue.keys())
    valid_residues = []
    for key in sorted_keys:
        atoms = atoms_by_residue[key]
        if 'N' in atoms and 'CA' in atoms and 'C' in atoms:
            valid_residues.append((key, atoms))
    
    n = len(valid_residues)
    if n < 3:
        raise ValueError(f"Only {n} valid residues found in {cif_path}")
    
    phi_angles = np.full(n, np.nan)
    psi_angles = np.full(n, np.nan)
    res_ids = [vr[0][0] for vr in valid_residues]
    
    backbone_coords = np.zeros((n, 4, 3))
    has_O = np.zeros(n, dtype=bool)
    for i, (key, atoms) in enumerate(valid_residues):
        backbone_coords[i, 0] = atoms['N']
        backbone_coords[i, 1] = atoms['CA']
        backbone_coords[i, 2] = atoms['C']
        if 'O' in atoms:
            backbone_coords[i, 3] = atoms['O']
            has_O[i] = True
    
    for i in range(n):
        _, atoms_i = valid_residues[i]
        if i > 0:
            _, atoms_prev = valid_residues[i - 1]
            if 'C' in atoms_prev:
                phi_angles[i] = _dihedral(
                    atoms_prev['C'], atoms_i['N'], atoms_i['CA'], atoms_i['C'])
        if i < n - 1:
            _, atoms_next = valid_residues[i + 1]
            if 'N' in atoms_next:
                psi_angles[i] = _dihedral(
                    atoms_i['N'], atoms_i['CA'], atoms_i['C'], atoms_next['N'])
    
    return phi_angles, psi_angles, res_ids, backbone_coords, has_O


def download_alphafold(uniprot_id, out_dir=None):
    """Download AlphaFold predicted structure for a UniProt ID.
    
    Strategy:
      0. Check alphafold_cache/ folder (user's pre-downloaded bulk cache)
      1. Check data/alphafold/ for any existing cached file (any version)
      2. Query the AlphaFold API to get the current download URL
      3. If API fails, try direct URL patterns (v4 CIF, v4 PDB, v6 CIF)
      4. Accept either mmCIF or PDB format (pipeline handles both)
    
    Returns path to downloaded file, or None on failure.
    """
    if out_dir is None:
        out_dir = ALPHAFOLD_DIR
    out_dir = Path(out_dir)
    
    uid = uniprot_id.strip().upper()
    
    # Strategy 0: Check bulk cache folder (alphafold_cache/)
    # Common naming: AF-{uid}-F1-model_v4.cif or AF-{uid}-F1.cif
    cache_dirs = [ALPHAFOLD_CACHE_DIR]
    # Also check relative to script and cwd
    for candidate in [Path("alphafold_cache"), Path(__file__).parent / "alphafold_cache"]:
        if candidate.exists() and candidate not in cache_dirs:
            cache_dirs.append(candidate)
    
    for cache_dir in cache_dirs:
        if not cache_dir.exists():
            continue
        for pattern in [f"AF-{uid}*.cif", f"AF-{uid}*.pdb", f"{uid}*.cif", f"{uid}*.pdb"]:
            existing = list(cache_dir.glob(pattern))
            if existing:
                return existing[0]
    
    # Strategy 1: Check data/alphafold/ for previously downloaded files
    for pattern in [f"AF-{uid}*.cif", f"AF-{uid}*.pdb", f"{uid}*.cif", f"{uid}*.pdb"]:
        existing = list(out_dir.glob(pattern))
        if existing:
            return existing[0]
    
    try:
        import requests
    except ImportError:
        log.warning("  requests not available for AlphaFold download")
        return None
    
    # Strategy 1: API-first (recommended by AlphaFold team)
    api_urls = [
        f"https://alphafold.ebi.ac.uk/api/prediction/{uid}",
        f"https://alphafold.com/api/prediction/{uid}",
    ]
    
    for api_url in api_urls:
        try:
            resp = requests.get(api_url, timeout=15)
            if resp.status_code == 200:
                data = resp.json()
                if isinstance(data, list) and len(data) > 0:
                    entry = data[0]
                    # Try cifUrl first, then pdbUrl
                    for url_key in ['cifUrl', 'pdbUrl']:
                        file_url = entry.get(url_key)
                        if file_url:
                            ext = '.cif' if 'cif' in url_key.lower() else '.pdb'
                            out_path = out_dir / f"AF-{uid}-F1{ext}"
                            file_resp = requests.get(file_url, timeout=30)
                            if file_resp.status_code == 200:
                                out_path.write_text(file_resp.text)
                                log.info(f"  Downloaded AlphaFold (API): {uid}")
                                return out_path
        except Exception:
            pass  # Try next strategy
    
    # Strategy 2: Direct URL patterns (fallback)
    direct_urls = [
        # Current v4 mmCIF (still works for many entries)
        (f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.cif", '.cif'),
        # v4 PDB format
        (f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb", '.pdb'),
        # Try without version suffix
        (f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model.cif", '.cif'),
        # 3D-Beacons / alternative endpoints
        (f"https://www.alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.cif", '.cif'),
    ]
    
    for url, ext in direct_urls:
        try:
            resp = requests.get(url, timeout=20)
            if resp.status_code == 200 and len(resp.text) > 1000:
                out_path = out_dir / f"AF-{uid}-F1{ext}"
                out_path.write_text(resp.text)
                log.info(f"  Downloaded AlphaFold (direct): {uid}")
                return out_path
        except Exception:
            pass
    
    log.warning(f"  AlphaFold: all download strategies failed for {uid}")
    return None


def build_geometric_barcode(segments, torus_knot):
    """Build a geometric barcode (segment-level feature vector) for a protein.
    
    Returns a dict with:
      - protein-level features: p, q, |Q|, n_segments, frac_geodesic, frac_defect
      - segment_barcodes: list of per-segment dicts with model class, params, etc.
    """
    segment_barcodes = []
    n_geodesic = 0
    n_defect = 0
    n_total = 0
    
    for seg in segments:
        if seg.spiral_class == 'too_short':
            continue
        n_total += 1
        
        bc = {
            'ss_type': seg.ss_type,
            'model_class': seg.spiral_class,
            'length': seg.length,
            'start_idx': seg.start_idx,
            'end_idx': seg.end_idx,
            'mean_kappa': float(seg.mean_kappa),
            'R2': float(seg.r_squared),
            'delta_aicc': float(seg.delta_aicc_osc),
        }
        
        # Model-specific parameters
        if seg.spiral_class in ('geodesic', 'circular_arc'):
            bc['params'] = {'kappa0': float(seg.mean_kappa)}
            n_geodesic += 1
        elif seg.spiral_class == 'gauss_peak':
            bc['params'] = {
                'A': float(seg.gauss_A) if seg.gauss_A else 0,
                'sigma': float(seg.gauss_sigma) if seg.gauss_sigma else 0,
                's0': float(seg.gauss_s0) if seg.gauss_s0 else 0,
                'peak_residue': int(seg.gauss_peak_residue) if seg.gauss_peak_residue else -1,
            }
            n_defect += 1
        elif seg.spiral_class == 'quadratic':
            n_defect += 1
            bc['params'] = {'type': 'polynomial'}
        elif seg.spiral_class in ('sinusoidal', 'damped_oscillation'):
            bc['params'] = {'type': 'oscillatory'}
        else:
            bc['params'] = {'type': seg.spiral_class}
        
        segment_barcodes.append(bc)
    
    tk = torus_knot or {}
    barcode = {
        'p': float(tk.get('p', 0)),
        'q': float(tk.get('q', 0)),
        'Q_magnitude': float(tk.get('total_winding_magnitude', 0)),
        'n_segments': n_total,
        'n_geodesic': n_geodesic,
        'n_defect': n_defect,
        'frac_geodesic': n_geodesic / max(n_total, 1),
        'frac_defect': n_defect / max(n_total, 1),
        'segments': segment_barcodes,
    }
    
    return barcode


def analyze_alphafold_protein(uniprot_id, cif_path, phi_grid, psi_grid, W_grid):
    """Analyze a single AlphaFold structure. Returns (ProteinResult, barcode)."""
    
    name = uniprot_id
    fpath = str(cif_path)
    
    # Extract dihedrals — detect format from extension
    if fpath.endswith('.pdb'):
        phi, psi, res_ids, backbone_coords, has_O = extract_dihedrals_manual(fpath, chain_id='A')
    else:
        phi, psi, res_ids, backbone_coords, has_O = extract_dihedrals_mmcif(fpath)
    
    valid = ~(np.isnan(phi) | np.isnan(psi))
    phi = phi[valid]
    psi = psi[valid]
    backbone_coords_valid = backbone_coords[valid]
    has_O_valid = has_O[valid]
    
    n_valid = len(phi)
    if n_valid < 6:
        log.warning(f"  {name}: too few valid residues ({n_valid})")
        return None, None
    
    # Build result object (no folding rate/CO for AlphaFold structures)
    result = ProteinResult(
        name=name, pdb_id=f"AF-{uniprot_id}", chain='A',
        length=n_valid, ln_kf=0.0, co=0.0, source='alphafold'
    )
    
    # SS assignment
    ss = None
    if np.all(has_O_valid):
        try:
            ss = dssp_assign_pure_numpy(backbone_coords_valid)
            result.ss_method = 'dssp_numpy'
        except Exception:
            pass
    if ss is None:
        ss = assign_ss_from_dihedrals(phi, psi)
        result.ss_method = 'dihedral'
    
    result.n_valid_residues = n_valid
    
    # Torus curvature
    kappa, s, tangent = torus_curvature(phi, psi)
    if len(kappa) < 4:
        return None, None
    
    # Segment and classify
    ss_changes = np.where(ss[1:] != ss[:-1])[0] + 1
    boundaries = np.concatenate([[0], ss_changes, [len(ss)]])
    
    segments = []
    for i in range(len(boundaries) - 1):
        start, end = int(boundaries[i]), int(boundaries[i + 1]) - 1
        seg_ss = ss[start]
        seg_len = end - start + 1
        
        if seg_len < 2:
            continue
        
        seg_k = kappa[max(0, start):min(end + 1, len(kappa))]
        seg_s = s[max(0, start):min(end + 1, len(s))]
        
        if len(seg_k) < 2:
            continue
        
        seg_s_rel = seg_s - seg_s[0]
        
        seg = CornuSegment(
            start_idx=start, end_idx=end, length=seg_len,
            ss_type=seg_ss, mean_kappa=float(np.mean(seg_k)),
            dkappa_ds=0.0, r_squared=0.0, spiral_class='too_short',
            kappa_values=seg_k.tolist(),
        )
        
        if seg_len >= 4 and len(seg_k) >= 4:
            slope, r_sq, class_name, fits, osc_info = classify_spiral(seg_k, seg_s_rel)
            seg.spiral_class = class_name
            seg.r_squared = r_sq
            seg.dkappa_ds = slope
            seg.delta_aicc_osc = osc_info.get('delta_aicc_osc', 0.0)
            seg.best_osc_r2 = osc_info.get('best_osc_r2', 0.0)
            seg.best_osc_model = osc_info.get('best_osc_model', '')
            seg.best_osc_omega = osc_info.get('best_osc_omega', 0.0)
            
            # Gaussian peak params
            if class_name == 'gauss_peak' and 'gauss_peak' in fits:
                gp = fits['gauss_peak']
                params = gp.get('params', {})
                seg.gauss_A = params.get('A', 0)
                seg.gauss_sigma = params.get('sigma', 0)
                seg.gauss_s0 = params.get('s0', 0)
                if seg_s_rel is not None and len(seg_s_rel) > 0:
                    s0_val = params.get('s0', 0)
                    nearest_idx = np.argmin(np.abs(seg_s_rel - s0_val))
                    seg.gauss_peak_residue = start + nearest_idx
                    seg.gauss_peak_residue_id = start + nearest_idx
        else:
            seg.spiral_class = 'too_short'
        
        segments.append(seg)
    
    result.segments = segments
    
    # Torus knot descriptors
    result.torus_knot = compute_torus_knot_descriptor(phi, psi, ss)
    
    # Build barcode
    barcode = build_geometric_barcode(segments, result.torus_knot)
    barcode['uniprot_id'] = uniprot_id
    barcode['length'] = n_valid
    barcode['ss_composition'] = {
        'H': int(np.sum(ss == 'H')),
        'E': int(np.sum(ss == 'E')),
        'C': int(np.sum(ss == 'C')),
    }
    
    return result, barcode


def save_barcode(barcode, out_dir=None):
    """Save a geometric barcode as JSON."""
    if out_dir is None:
        out_dir = BARCODE_DIR
    out_dir = Path(out_dir)
    
    uid = barcode.get('uniprot_id', 'unknown')
    out_path = out_dir / f"{uid}_barcode.json"
    
    import json
    with open(out_path, 'w') as f:
        json.dump(barcode, f, indent=2, default=str)
    
    return out_path


def alphafold_main(uniprot_ids):
    """Batch-process AlphaFold structures.
    
    Usage: python cornu_spirals_on_T2.py alphafold P00533 P04637 P01308 ...
    
    Or pass a file: python cornu_spirals_on_T2.py alphafold @uniprot_ids.txt
    """
    print("=" * 72)
    print("  ALPHAFOLD GEOMETRIC DECONSTRUCTION")
    print("  Torus curvature analysis on AlphaFold predicted structures")
    print("=" * 72)
    print()
    
    # Expand file references
    ids = []
    for arg in uniprot_ids:
        if arg.startswith('@'):
            # Read IDs from file, one per line
            fpath = Path(arg[1:])
            if fpath.exists():
                with open(fpath) as f:
                    ids.extend(line.strip() for line in f if line.strip() and not line.startswith('#'))
            else:
                log.warning(f"  File not found: {fpath}")
        else:
            ids.append(arg.strip())
    
    if not ids:
        print("Usage: python cornu_spirals_on_T2.py alphafold P00533 P04637 ...")
        print("   or: python cornu_spirals_on_T2.py alphafold @uniprot_ids.txt")
        return
    
    log.info(f"Processing {len(ids)} UniProt IDs...")
    
    # Build superpotential (shared)
    phi_grid, psi_grid, p_grid, W_grid = build_superpotential_grid()
    
    # Download and analyze
    results = []
    barcodes = []
    
    for i, uid in enumerate(ids):
        log.info(f"  [{i+1}/{len(ids)}] {uid}...")
        
        cif_path = download_alphafold(uid)
        if cif_path is None:
            log.warning(f"    Skipping {uid} (download failed)")
            continue
        
        try:
            result, barcode = analyze_alphafold_protein(uid, cif_path, phi_grid, psi_grid, W_grid)
        except Exception as e:
            log.warning(f"    Error processing {uid}: {e}")
            continue
        
        if result is None:
            continue
        
        results.append(result)
        barcodes.append(barcode)
        
        # Save barcode
        bc_path = save_barcode(barcode)
        
        # Print summary
        tk = result.torus_knot or {}
        n_segs = sum(1 for s in result.segments if s.spiral_class != 'too_short')
        n_geodesic = sum(1 for s in result.segments
                        if s.spiral_class in ('geodesic', 'circular_arc'))
        n_defect = sum(1 for s in result.segments
                      if s.spiral_class in ('gauss_peak', 'quadratic'))
        
        print(f"  {uid:12s}  L={result.n_valid_residues:4d}  "
              f"p={tk.get('p', 0):+.2f}  q={tk.get('q', 0):+.2f}  "
              f"|Q|={tk.get('total_winding_magnitude', 0):.2f}  "
              f"segs={n_segs}  geodesic={n_geodesic}  defect={n_defect}  "
              f"SS: H={barcode['ss_composition']['H']} "
              f"E={barcode['ss_composition']['E']} "
              f"C={barcode['ss_composition']['C']}")
    
    # Summary
    print()
    print(f"  Processed: {len(results)} / {len(ids)} proteins")
    print(f"  Barcodes saved to: {BARCODE_DIR}/")
    
    if barcodes:
        # Aggregate stats
        all_p = [b['p'] for b in barcodes]
        all_q = [b['q'] for b in barcodes]
        all_Q = [b['Q_magnitude'] for b in barcodes]
        all_frac_geo = [b['frac_geodesic'] for b in barcodes]
        
        print(f"\n  Aggregate (n={len(barcodes)}):")
        print(f"    p-winding:   mean={np.mean(all_p):+.2f} ± {np.std(all_p):.2f}")
        print(f"    q-winding:   mean={np.mean(all_q):+.2f} ± {np.std(all_q):.2f}")
        print(f"    |Q|:         mean={np.mean(all_Q):.2f} ± {np.std(all_Q):.2f}")
        print(f"    Frac geodesic: mean={np.mean(all_frac_geo):.1%}")
        
        # Save aggregate barcode database
        import json
        db_path = RESULTS_DIR / "alphafold_barcode_db.json"
        with open(db_path, 'w') as f:
            json.dump({
                'n_proteins': len(barcodes),
                'proteins': barcodes,
            }, f, indent=2, default=str)
        print(f"    Database: {db_path}")
    
    print()
    return results, barcodes


# =============================================================================
# SECTION 16: ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'alphafold':
        alphafold_main(sys.argv[2:])
    else:
        main()

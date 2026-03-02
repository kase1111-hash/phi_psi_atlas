#!/usr/bin/env python3
"""
BACKBONE STRAIN ANALYZER
========================
Interactive walkthrough of torsional curvature analysis.

Usage: python strain_analyzer.py <apo.cif> <drug.cif> <chain> [--pocket-ref <ref.cif>]

Walks through 8 steps with explanations. Each step prints results
and offers [A/B/C/D] options for deeper exploration.

Requires: numpy, scipy (pip install numpy scipy)
PDB files: download mmCIF from https://www.rcsb.org/

Example:
  python strain_analyzer.py 1HHP.cif 1HSG.cif A
  python strain_analyzer.py 2CBA.cif 1YDA.cif A
"""

import sys
import os
import numpy as np
from collections import defaultdict

# ═══════════════════════════════════════════════════════
# CORE COMPUTATION ENGINE
# ═══════════════════════════════════════════════════════

def parse_mmcif(filepath, chain_id):
    """Parse mmCIF file. Returns backbone atoms, B-factors, and ligand coords."""
    backbone = []
    bfacs = {}
    ligands_by_comp = defaultdict(list)
    
    in_atom_site = False
    col_names = []
    
    excl = {'HOH','DOD','WAT','SO4','PO4','GOL','EDO','ACT','DMS',
            'CL','NA','ZN','CA','MG','MN','FE','CU','CO','NI','CD',
            'BME','MPD','PEG','EPE','TRS','MES','CIT','SCN','IPA',
            'FMT','IOD','BR','NH4'}
    
    with open(filepath) as f:
        for line in f:
            if line.startswith('_atom_site.'):
                col_names.append(line.strip().split('.')[1])
                in_atom_site = True
                continue
            if in_atom_site and not line.startswith('_') and not line.startswith('#') and not line.startswith('loop_'):
                if line.strip() == '' or line.startswith('#'):
                    in_atom_site = False
                    continue
                parts = line.split()
                if len(parts) < len(col_names):
                    continue
                rec = {col_names[i]: parts[i] for i in range(len(col_names))}
                if rec.get('pdbx_PDB_model_num', '1') != '1':
                    continue
                try:
                    x = float(rec['Cartn_x'])
                    y = float(rec['Cartn_y'])
                    z = float(rec['Cartn_z'])
                except (KeyError, ValueError):
                    continue
                
                group = rec.get('group_PDB', '')
                auth_chain = rec.get('auth_asym_id', rec.get('label_asym_id', ''))
                
                # Filter alternate conformations: only use first (A) or none (.)
                alt_id = rec.get('label_alt_id', '.')
                if alt_id not in ('.', '?', '', 'A'):
                    continue
                
                if group == 'HETATM':
                    comp = rec.get('auth_comp_id', rec.get('label_comp_id', ''))
                    if comp not in excl:
                        ligands_by_comp[comp].append(np.array([x, y, z]))
                
                if group == 'ATOM' and auth_chain == chain_id:
                    atom = rec.get('auth_atom_id', rec.get('label_atom_id', ''))
                    rn = int(rec.get('auth_seq_id', 0))
                    if atom == 'CA':
                        try:
                            bfacs[rn] = float(rec.get('B_iso_or_equiv', '0'))
                        except ValueError:
                            pass
                    if atom in ('N', 'CA', 'C'):
                        backbone.append({
                            'atom': atom, 'resnum': rn,
                            'x': x, 'y': y, 'z': z
                        })
    
    return backbone, bfacs, ligands_by_comp


def dihedral(p0, p1, p2, p3):
    """Compute dihedral angle between four points."""
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    n1 = np.linalg.norm(b1)
    if n1 < 1e-10:
        return np.nan
    b1 = b1 / n1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    if np.linalg.norm(v) < 1e-10 or np.linalg.norm(w) < 1e-10:
        return np.nan
    return np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))


def compute_torsion_angles(backbone):
    """Compute phi/psi for all residues."""
    res = {}
    for a in backbone:
        rn = a['resnum']
        if rn not in res:
            res[rn] = {}
        res[rn][a['atom']] = np.array([a['x'], a['y'], a['z']])
    
    sr = sorted(res.keys())
    phi = {}
    psi = {}
    
    for i, r in enumerate(sr):
        d = res[r]
        if 'N' not in d or 'CA' not in d or 'C' not in d:
            continue
        if i > 0 and 'C' in res.get(sr[i-1], {}):
            phi[r] = dihedral(res[sr[i-1]]['C'], d['N'], d['CA'], d['C'])
        else:
            phi[r] = np.nan
        if i < len(sr) - 1 and 'N' in res.get(sr[i+1], {}):
            psi[r] = dihedral(d['N'], d['CA'], d['C'], res[sr[i+1]]['N'])
        else:
            psi[r] = np.nan
    
    return phi, psi, res


def compute_kappa_squared(phi, psi):
    """Compute discrete torsional curvature for all residues."""
    residues = sorted(set(phi.keys()) & set(psi.keys()))
    ksq = {}
    
    for idx, r in enumerate(residues):
        if idx < 1 or idx >= len(residues) - 1:
            ksq[r] = 0.0
            continue
        
        prev = residues[idx - 1]
        nxt = residues[idx + 1]
        
        vals = [
            phi.get(prev, np.nan), phi.get(r, np.nan), phi.get(nxt, np.nan),
            psi.get(prev, np.nan), psi.get(r, np.nan), psi.get(nxt, np.nan)
        ]
        
        if any(np.isnan(x) for x in vals):
            ksq[r] = 0.0
            continue
        
        # Tangent changes in phi and psi
        dp1 = np.arctan2(np.sin(vals[1] - vals[0]), np.cos(vals[1] - vals[0]))
        ds1 = np.arctan2(np.sin(vals[4] - vals[3]), np.cos(vals[4] - vals[3]))
        dp2 = np.arctan2(np.sin(vals[2] - vals[1]), np.cos(vals[2] - vals[1]))
        ds2 = np.arctan2(np.sin(vals[5] - vals[4]), np.cos(vals[5] - vals[4]))
        
        # Arc lengths
        a1 = max(np.sqrt(dp1**2 + ds1**2), 1e-10)
        a2 = max(np.sqrt(dp2**2 + ds2**2), 1e-10)
        
        # Curvature squared (Frenet-Serret, arc-length normalized)
        ksq[r] = ((dp2/a2 - dp1/a1)**2 + (ds2/a2 - ds1/a1)**2) / (0.5*(a1+a2))**2
    
    return ksq


def gini_coefficient(values):
    """Compute Gini coefficient of a distribution."""
    x = np.sort(np.abs(np.array(values)))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    return float((2 * np.sum(np.arange(1, n+1) * x) - (n+1) * np.sum(x)) / (n * np.sum(x)))


def get_primary_ligand(ligands_by_comp):
    """Return the largest non-cofactor ligand."""
    excl_cof = {'HEM', 'NAG', 'BOG', 'MAN', 'CEL'}
    cands = {c: v for c, v in ligands_by_comp.items() if c not in excl_cof}
    if not cands:
        cands = dict(ligands_by_comp)
    if not cands:
        return [], "NONE"
    best = max(cands.items(), key=lambda x: len(x[1]))
    return best[1], best[0]


def rama_class(phi_val, psi_val):
    """Classify residue by Ramachandran region."""
    if np.isnan(phi_val) or np.isnan(psi_val):
        return 'other'
    pd = np.degrees(phi_val)
    sd = np.degrees(psi_val)
    if -120 < pd < -20 and -80 < sd < -10:
        return 'alpha'
    if -180 < pd < -50 and (sd > 50 or sd < -120):
        return 'beta'
    return 'loop'


# ═══════════════════════════════════════════════════════
# DISPLAY ENGINE
# ═══════════════════════════════════════════════════════

def banner(text):
    w = 72
    print()
    print("═" * w)
    print(f"  {text}")
    print("═" * w)


def step_header(num, title):
    print()
    print(f"  ┌{'─'*68}┐")
    print(f"  │  STEP {num}: {title:<58s}│")
    print(f"  └{'─'*68}┘")


def explain(text):
    """Print wrapped explanation text."""
    words = text.split()
    line = "  "
    for w in words:
        if len(line) + len(w) + 1 > 72:
            print(line)
            line = "  " + w
        else:
            line += " " + w if line.strip() else "  " + w
    if line.strip():
        print(line)


def key_finding(text):
    print()
    print(f"  ▸ {text}")


def offer_options(options):
    """Show A/B/C/D options and handle user choice."""
    print()
    for key, desc in options.items():
        print(f"    [{key}] {desc}")
    print(f"    [Enter] Continue to next step")
    print()
    
    while True:
        try:
            choice = input("  Choice: ").strip().upper()
        except (EOFError, KeyboardInterrupt):
            choice = ""
        if choice == "" or choice not in options:
            return choice
        return choice


def show_histogram(values, label, bins=20):
    """Simple ASCII histogram."""
    vals = np.array(values)
    counts, edges = np.histogram(vals[vals > 0], bins=bins)
    max_c = max(counts) if max(counts) > 0 else 1
    print(f"\n  Distribution of {label}:")
    for i in range(len(counts)):
        bar = "█" * int(40 * counts[i] / max_c)
        print(f"  {edges[i]:>8.1f} │{bar}")


# ═══════════════════════════════════════════════════════
# MAIN WALKTHROUGH
# ═══════════════════════════════════════════════════════

def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    
    apo_file = sys.argv[1]
    drug_file = sys.argv[2]
    chain = sys.argv[3]
    
    pocket_ref_file = None
    if '--pocket-ref' in sys.argv:
        idx = sys.argv.index('--pocket-ref')
        if idx + 1 < len(sys.argv):
            pocket_ref_file = sys.argv[idx + 1]
    
    banner("BACKBONE STRAIN ANALYZER")
    print(f"  Apo structure:   {apo_file}")
    print(f"  Drug structure:  {drug_file}")
    print(f"  Chain:           {chain}")
    if pocket_ref_file:
        print(f"  Pocket ref:      {pocket_ref_file}")
    
    # ═══ STEP 1: LOAD ═══
    step_header(1, "Load Structures")
    explain("Parsing mmCIF files to extract backbone atoms (N, CA, C), "
            "B-factors, and ligand coordinates. Each residue needs all "
            "three backbone atoms for torsion angle computation.")
    
    bb_apo, bf_apo, lig_apo = parse_mmcif(apo_file, chain)
    bb_drug, bf_drug, lig_drug = parse_mmcif(drug_file, chain)
    
    n_res_apo = len(set(a['resnum'] for a in bb_apo))
    n_res_drug = len(set(a['resnum'] for a in bb_drug))
    lig_coords, lig_name = get_primary_ligand(lig_drug)
    
    print(f"\n  Apo:  {n_res_apo} residues, {len(bb_apo)} backbone atoms")
    print(f"  Drug: {n_res_drug} residues, {len(bb_drug)} backbone atoms")
    print(f"  Primary ligand: {lig_name} ({len(lig_coords)} atoms)")
    
    choice = offer_options({
        'A': 'Show all ligand compounds found',
        'B': 'Show residue number range',
    })
    if choice == 'A':
        for comp, coords in sorted(lig_drug.items()):
            print(f"    {comp}: {len(coords)} atoms")
    elif choice == 'B':
        apo_res = sorted(set(a['resnum'] for a in bb_apo))
        drug_res = sorted(set(a['resnum'] for a in bb_drug))
        print(f"    Apo:  {apo_res[0]} to {apo_res[-1]}")
        print(f"    Drug: {drug_res[0]} to {drug_res[-1]}")
    
    # ═══ STEP 2: TORSION ANGLES ═══
    step_header(2, "Compute Torsion Angles (phi/psi)")
    explain("For each residue i, phi is the dihedral C(i-1)-N(i)-CA(i)-C(i), "
            "and psi is N(i)-CA(i)-C(i)-N(i+1). These define the backbone "
            "trajectory in Ramachandran conformational space.")
    
    phi_apo, psi_apo, res_apo = compute_torsion_angles(bb_apo)
    phi_drug, psi_drug, res_drug = compute_torsion_angles(bb_drug)
    
    valid_apo = sum(1 for r in phi_apo if not np.isnan(phi_apo[r]) and not np.isnan(psi_apo.get(r, np.nan)))
    valid_drug = sum(1 for r in phi_drug if not np.isnan(phi_drug[r]) and not np.isnan(psi_drug.get(r, np.nan)))
    
    print(f"\n  Apo:  {valid_apo} residues with valid phi/psi")
    print(f"  Drug: {valid_drug} residues with valid phi/psi")
    
    # SS composition
    ss_apo = defaultdict(int)
    for r in phi_apo:
        ss_apo[rama_class(phi_apo.get(r, np.nan), psi_apo.get(r, np.nan))] += 1
    
    print(f"\n  Secondary structure (apo, by Ramachandran):")
    for ss in ['alpha', 'beta', 'loop']:
        pct = 100 * ss_apo[ss] / sum(ss_apo.values()) if ss_apo else 0
        print(f"    {ss:>5s}: {ss_apo[ss]:>4d} ({pct:.0f}%)")
    
    choice = offer_options({
        'A': 'Show Ramachandran statistics',
        'B': 'Show largest phi/psi changes between structures',
    })
    if choice == 'A':
        phi_vals = [np.degrees(phi_apo[r]) for r in phi_apo if not np.isnan(phi_apo[r])]
        psi_vals = [np.degrees(psi_apo[r]) for r in psi_apo if not np.isnan(psi_apo[r])]
        print(f"    Phi: mean={np.mean(phi_vals):.1f}, SD={np.std(phi_vals):.1f}")
        print(f"    Psi: mean={np.mean(psi_vals):.1f}, SD={np.std(psi_vals):.1f}")
    elif choice == 'B':
        common = sorted(set(phi_apo.keys()) & set(phi_drug.keys()))
        changes = []
        for r in common:
            pa, pd = phi_apo.get(r, np.nan), phi_drug.get(r, np.nan)
            sa, sd = psi_apo.get(r, np.nan), psi_drug.get(r, np.nan)
            if any(np.isnan(x) for x in [pa, pd, sa, sd]):
                continue
            dp = abs(np.arctan2(np.sin(pd-pa), np.cos(pd-pa)))
            ds = abs(np.arctan2(np.sin(sd-sa), np.cos(sd-sa)))
            changes.append((r, np.degrees(dp), np.degrees(ds)))
        changes.sort(key=lambda x: -max(x[1], x[2]))
        print(f"    Top 10 largest phi/psi changes (circular-wrapped):")
        print(f"    {'Res':>4s}  {'|Δphi|':>8s}  {'|Δpsi|':>8s}")
        for r, dp, ds in changes[:10]:
            print(f"    {r:>4d}  {dp:>7.1f}°  {ds:>7.1f}°")
    
    # ═══ STEP 3: KAPPA-SQUARED ═══
    step_header(3, "Compute Kappa-Squared (Torsional Curvature)")
    explain("Kappa-squared measures how sharply the backbone trajectory "
            "bends in conformational space at each residue, normalized by "
            "arc-length. It approximates the geometric component of bending "
            "energy density in the Kirchhoff elastic rod framework.")
    
    ksq_apo = compute_kappa_squared(phi_apo, psi_apo)
    ksq_drug = compute_kappa_squared(phi_drug, psi_drug)
    
    common_all = sorted(set(ksq_apo.keys()) & set(ksq_drug.keys()))
    # Exclude residues where κ²=0 in either structure (terminal/NaN residues)
    common = [r for r in common_all if ksq_apo[r] > 0 or ksq_drug[r] > 0]
    
    k_apo_vals = [ksq_apo[r] for r in common]
    k_drug_vals = [ksq_drug[r] for r in common]
    
    n_excluded = len(common_all) - len(common)
    print(f"\n  Common residues: {len(common)} ({n_excluded} terminal/invalid excluded)")
    print(f"\n  Apo kappa-sq:  mean={np.mean(k_apo_vals):.1f}  "
          f"median={np.median(k_apo_vals):.1f}  max={max(k_apo_vals):.1f}")
    print(f"  Drug kappa-sq: mean={np.mean(k_drug_vals):.1f}  "
          f"median={np.median(k_drug_vals):.1f}  max={max(k_drug_vals):.1f}")
    
    # Top kappa-sq residues
    top_apo = sorted(common, key=lambda r: -ksq_apo[r])[:5]
    print(f"\n  Top 5 kappa-sq residues (apo):")
    for r in top_apo:
        ss = rama_class(phi_apo.get(r, np.nan), psi_apo.get(r, np.nan))
        print(f"    Residue {r:>4d} ({ss:>5s}): kappa-sq = {ksq_apo[r]:.1f}")
    
    choice = offer_options({
        'A': 'Show kappa-sq distribution histogram',
        'B': 'Show kappa-sq by secondary structure',
        'C': 'Explain the Frenet-Serret normalization',
    })
    if choice == 'A':
        show_histogram(k_apo_vals, "kappa-squared (apo)")
    elif choice == 'B':
        for ss in ['alpha', 'beta', 'loop']:
            ss_vals = [ksq_apo[r] for r in common 
                       if rama_class(phi_apo.get(r,np.nan), psi_apo.get(r,np.nan)) == ss]
            if ss_vals:
                print(f"    {ss:>5s}: mean={np.mean(ss_vals):.1f}  n={len(ss_vals)}")
    elif choice == 'C':
        explain("The Frenet-Serret framework defines curvature as |dT/ds| where "
                "T is the unit tangent and s is arc-length. We compute this discretely: "
                "tangent = direction of (delta-phi, delta-psi), arc-length = magnitude of "
                "(delta-phi, delta-psi). Dividing by arc-length squared amplifies deviations "
                "in regions with small step sizes (helices: ~38 deg) versus large step sizes "
                "(sheets: ~84 deg). This is physically meaningful: a 1-degree deviation in "
                "a helix represents more mechanical strain than in a sheet.")
    
    # ═══ STEP 4: GINI COEFFICIENT ═══
    step_header(4, "Compute Gini Coefficient")
    explain("The Gini coefficient measures inequality. Gini = 0 means all residues "
            "carry equal kappa-squared. Gini = 1 means one residue carries all of it. "
            "In the Kirchhoff framework, this measures how unevenly bending energy "
            "is distributed along the backbone.")
    
    g_apo = gini_coefficient(k_apo_vals)
    g_drug = gini_coefficient(k_drug_vals)
    
    print(f"\n  Gini (apo):  {g_apo:.4f}")
    print(f"  Gini (drug): {g_drug:.4f}")
    
    key_finding(f"Population baseline: 0.805 +/- 0.065 (107 crystal structures)")
    
    choice = offer_options({
        'A': 'Explain what high/low Gini means physically',
        'B': 'Show how Gini relates to rod mechanics',
    })
    if choice == 'A':
        explain("High Gini (>0.85): Bending energy concentrated in few hot spots. "
                "The backbone has sharp kinks separated by smooth regions. Brittle rod. "
                "Low Gini (<0.75): Bending energy distributed broadly. The backbone "
                "curves gently everywhere. Resilient rod.")
    elif choice == 'B':
        explain("In a clamped elastic beam, bending moment is maximum at the center "
                "and zero at the free end. This creates unequal curvature distribution "
                "(high Gini). A uniformly loaded beam has more even curvature (lower Gini). "
                "The Gini of kappa-squared is the structural analog: it tells you whether "
                "the protein behaves more like a clamped beam or a uniformly loaded one.")
    
    # ═══ STEP 5: DELTA-GINI ═══
    step_header(5, "Delta-Gini: Direction and Magnitude")
    
    dg = g_drug - g_apo
    direction = "STIFFENED" if dg > 0 else "LOOSENED" if dg < 0 else "UNCHANGED"
    
    print(f"\n  Delta-Gini = {dg:+.4f}")
    print(f"  Direction:   {direction}")
    print(f"  Magnitude:   |{dg:.4f}|")
    
    if abs(dg) < 0.012:
        key_finding("Below noise floor (0.012). Not significant.")
    elif abs(dg) < 0.03:
        key_finding("Small effect. May not survive bootstrap.")
    else:
        key_finding(f"Substantial effect ({abs(dg):.3f}).")
    
    explain(f"The drug {direction.lower()} the backbone: "
            f"{'concentrated' if dg > 0 else 'distributed'} bending energy "
            f"{'into fewer' if dg > 0 else 'across more'} positions.")
    
    choice = offer_options({
        'A': 'Run bootstrap confidence interval (slow)',
        'B': 'Compare to known target ranges',
    })
    if choice == 'A':
        np.random.seed(42)
        boot_dg = []
        for _ in range(2000):
            idx = np.random.choice(len(common), len(common), replace=True)
            boot_k_a = [k_apo_vals[i] for i in idx]
            boot_k_d = [k_drug_vals[i] for i in idx]
            boot_dg.append(gini_coefficient(boot_k_d) - gini_coefficient(boot_k_a))
        lo, hi = np.percentile(boot_dg, [2.5, 97.5])
        print(f"    95% CI: [{lo:+.4f}, {hi:+.4f}]")
        if lo > 0 or hi < 0:
            print(f"    *** SIGNIFICANT (CI excludes zero)")
        else:
            print(f"    Not significant (CI includes zero)")
    elif choice == 'B':
        print(f"    CA-II panel (n=16):  +0.004 to +0.165 (all stiffen)")
        print(f"    HIV-PR panel (n=27): -0.036 to +0.214 (92% stiffen)")
        print(f"    CDK2 panel (n=19):   -0.020 to +0.005 (71% loosen)")
        print(f"    Your result:         {dg:+.4f}")
    
    # ═══ STEP 6: HOTSPOT RESIDUES ═══
    step_header(6, "Identify Hotspot Residues")
    explain("Hotspots are residues where |delta-kappa-squared| is largest. "
            "These are positions where the drug most changed the backbone's "
            "bending energy density.")
    
    dksq = {r: ksq_drug[r] - ksq_apo[r] for r in common}
    ranked = sorted(common, key=lambda r: -abs(dksq[r]))
    
    print(f"\n  Top 20 most-perturbed residues:")
    print(f"  {'Rank':>4s}  {'Res':>4s}  {'SS':>5s}  {'delta-kappa-sq':>14s}  {'Direction':>10s}")
    print(f"  {'─'*4}  {'─'*4}  {'─'*5}  {'─'*14}  {'─'*10}")
    
    for i in range(min(20, len(ranked))):
        r = ranked[i]
        ss = rama_class(phi_apo.get(r, np.nan), psi_apo.get(r, np.nan))
        d = "stiffen" if dksq[r] > 0 else "loosen"
        print(f"  {i+1:>4d}  {r:>4d}  {ss:>5s}  {dksq[r]:>+14.1f}  {d:>10s}")
    
    choice = offer_options({
        'A': 'Show secondary structure composition of hotspots',
        'B': 'Show direction consistency of top 20',
    })
    if choice == 'A':
        ss_counts = defaultdict(int)
        for r in ranked[:20]:
            ss_counts[rama_class(phi_apo.get(r,np.nan), psi_apo.get(r,np.nan))] += 1
        for ss, ct in sorted(ss_counts.items(), key=lambda x: -x[1]):
            print(f"    {ss:>5s}: {ct:>2d}/20 ({100*ct/20:.0f}%)")
    elif choice == 'B':
        n_pos = sum(1 for r in ranked[:20] if dksq[r] > 0)
        n_neg = 20 - n_pos
        print(f"    Stiffen: {n_pos}/20 ({100*n_pos/20:.0f}%)")
        print(f"    Loosen:  {n_neg}/20 ({100*n_neg/20:.0f}%)")
    
    # ═══ STEP 7: LOCAL vs DISTAL ═══
    step_header(7, "Classify Local vs Distal")
    explain("Local = within 10 Angstroms of the primary ligand. "
            "Distal = everything else. This tells you whether the drug's "
            "mechanical effect stays at the binding site or propagates.")
    
    if not lig_coords:
        print("\n  No ligand atoms found. Skipping local/distal analysis.")
    else:
        # Define pocket
        pocket = set()
        for r in common:
            ca = res_drug.get(r, {}).get('CA', res_apo.get(r, {}).get('CA'))
            if ca is not None:
                min_d = min(np.linalg.norm(ca - l) for l in lig_coords)
                if min_d < 10:
                    pocket.add(r)
        
        sphere_pct = 100 * len(pocket) / len(common)
        
        n_local = sum(1 for r in ranked[:20] if r in pocket)
        n_distal = 20 - n_local
        
        print(f"\n  Pocket size: {len(pocket)} residues ({sphere_pct:.0f}% of protein)")
        print(f"  Top 20 hotspots: {n_local} local, {n_distal} distal "
              f"({100*n_distal/20:.0f}% distal)")
        
        if sphere_pct >= 30:
            key_finding(f"WARNING: Pocket covers {sphere_pct:.0f}% of protein. "
                        f"Local/distal distinction has reduced power.")
        
        choice = offer_options({
            'A': 'Show distance to ligand for each hotspot',
            'B': 'Explain the 30% sphere rule',
            'C': 'Compare to known baselines',
        })
        if choice == 'A':
            for i in range(min(20, len(ranked))):
                r = ranked[i]
                ca = res_drug.get(r, {}).get('CA', res_apo.get(r, {}).get('CA'))
                if ca is not None:
                    min_d = min(np.linalg.norm(ca - l) for l in lig_coords)
                    zone = "LOCAL" if r in pocket else "DISTAL"
                    print(f"    Res {r:>4d}: {min_d:>6.1f} Angstroms  {zone}")
        elif choice == 'B':
            explain("When >30% of residues are within 10 Angstroms of the ligand, "
                    "nearly half the protein is 'local'. Random hotspots would be "
                    "~50% distal just by chance. The test loses statistical power. "
                    "This happens with large cofactors (heme) or small proteins (<100 res).")
        elif choice == 'C':
            print(f"    Clean drug comparisons (sphere <30%): 86% distal")
            print(f"    Apo-apo controls: 99% distal")
            print(f"    Your result: {100*n_distal/20:.0f}% distal (sphere {sphere_pct:.0f}%)")
    
    # ═══ STEP 8: MECHANICAL MODE ═══
    step_header(8, "Mechanical Mode Classification")
    explain("The local/distal decomposition reveals how the drug redistributes "
            "strain between compartments. This is the most informative output: "
            "it tells you whether the drug clamps the pocket, exports strain, "
            "or shifts the whole protein uniformly.")
    
    if lig_coords and pocket:
        local_res = [r for r in common if r in pocket]
        distal_res = [r for r in common if r not in pocket]
        
        if len(local_res) >= 5 and len(distal_res) >= 5:
            g_local_apo = gini_coefficient([ksq_apo[r] for r in local_res])
            g_local_drug = gini_coefficient([ksq_drug[r] for r in local_res])
            g_distal_apo = gini_coefficient([ksq_apo[r] for r in distal_res])
            g_distal_drug = gini_coefficient([ksq_drug[r] for r in distal_res])
            
            dg_local = g_local_drug - g_local_apo
            dg_distal = g_distal_drug - g_distal_apo
            
            print(f"\n  Delta-Gini (local):   {dg_local:+.4f}  "
                  f"({'stiffen' if dg_local > 0 else 'loosen'})")
            print(f"  Delta-Gini (distal):  {dg_distal:+.4f}  "
                  f"({'stiffen' if dg_distal > 0 else 'loosen'})")
            print(f"  Delta-Gini (global):  {dg:+.4f}")
            
            # Classify mode
            if abs(dg_local) < 0.01 and abs(dg_distal) < 0.01:
                mode = "MINIMAL EFFECT"
            elif np.sign(dg_local) == np.sign(dg_distal):
                mode = "UNIFORM"
            elif dg_local > 0 and dg_distal < 0:
                mode = "CONCENTRATOR"
            elif dg_local < 0 and dg_distal > 0:
                mode = "EXPORTER"
            else:
                mode = "MIXED"
            
            lam = dg_local / dg if abs(dg) > 0.001 else 0
            
            print(f"\n  Mechanical mode: {mode}")
            print(f"  Lambda (local/global): {lam:.2f}")
            
            key_finding(f"This drug is a {mode.lower()}.")
            
            if mode == "CONCENTRATOR":
                explain("The drug stiffens its binding environment while loosening "
                        "the rest of the protein. Like tafamidis on TTR: "
                        "kinetic stabilization by local clamping.")
            elif mode == "EXPORTER":
                explain("The drug loosens or relaxes the binding pocket while "
                        "stiffening remote regions. Like darunavir on HIV-PR.")
            elif mode == "UNIFORM":
                explain("The drug shifts the whole protein's strain pattern "
                        "in the same direction. Most common mode (43%).")
            
            choice = offer_options({
                'A': 'Show Lambda test details',
                'B': 'Compare to known drug modes',
            })
            if choice == 'A':
                print(f"    Local residues:  {len(local_res)}")
                print(f"    Distal residues: {len(distal_res)}")
                print(f"    Local Gini apo:  {g_local_apo:.4f}")
                print(f"    Local Gini drug: {g_local_drug:.4f}")
                print(f"    Distal Gini apo: {g_distal_apo:.4f}")
                print(f"    Distal Gini drug:{g_distal_drug:.4f}")
            elif choice == 'B':
                print(f"    Tafamidis/TTR:     CONCENTRATOR (local +0.155, distal -0.011)")
                print(f"    Darunavir/HIV-PR:  EXPORTER (local -0.156, distal +0.103)")
                print(f"    KNI-272/HIV-PR:    CONCENTRATOR (Lambda = 22.4)")
                print(f"    Your drug:         {mode} (Lambda = {lam:.2f})")
        else:
            print(f"\n  Not enough residues in local ({len(local_res)}) "
                  f"or distal ({len(distal_res)}) compartment.")
    else:
        print(f"\n  No ligand found. Cannot classify mechanical mode.")
    
    # ═══ SUMMARY ═══
    banner("ANALYSIS COMPLETE")
    print(f"""
  Structure pair: {os.path.basename(apo_file)} vs {os.path.basename(drug_file)}
  Chain {chain}, {len(common)} common residues
  
  Gini (apo):    {g_apo:.4f}
  Gini (drug):   {g_drug:.4f}
  Delta-Gini:    {dg:+.4f} ({direction})
  
  Ligand:        {lig_name}""")
    
    if lig_coords and pocket:
        print(f"  Pocket:        {len(pocket)} residues ({sphere_pct:.0f}%)")
        print(f"  Top 20:        {n_distal}/20 distal ({100*n_distal/20:.0f}%)")
        if len(local_res) >= 5 and len(distal_res) >= 5:
            print(f"  Mode:          {mode}")
            print(f"  Lambda:        {lam:.2f}")
    
    print(f"""
  ─────────────────────────────────────────
  For more context, see: The Backbone Strain Atlas v15
  Method: Branham 2026
""")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Backbone Curvature Atlas — Extended Validation (Round 2)
=========================================================
Adds:
  - 8 more temperature cryo/RT pairs (→ 20 total)
  - 6 more ATP synthase structures (→ 9 total, 27 β-chains)
  - α-subunit negative control for ATP synthase
  - Cross-species ATP synthase (yeast)

Usage:
  python extended_validation.py                # Run all
  python extended_validation.py --resume       # Resume from checkpoint
  python extended_validation.py --test temp    # Temperature only
  python extended_validation.py --test atp     # ATP synthase only

Requires: numpy, requests
"""

import numpy as np
import os, sys, json, time, argparse
from pathlib import Path

CACHE_DIR = Path("./pdb_cache")
CHECKPOINT_FILE = Path("./extended_checkpoint.json")
RESULTS_FILE = Path("./extended_results.json")

# ══════════════════════════════════════════════════════════════
# CORE PIPELINE (identical to atlas — single source of truth)
# ══════════════════════════════════════════════════════════════

def parse_cif_chain(filepath, chain_id):
    atoms = []; in_atom_site = False; col_names = []
    with open(filepath) as f:
        for line in f:
            if line.startswith('_atom_site.'):
                col_names.append(line.strip().split('.')[1]); in_atom_site = True; continue
            if in_atom_site and not line.startswith('_') and not line.startswith('#') and not line.startswith('loop_'):
                if line.strip() == '' or line.startswith('#'): in_atom_site = False; continue
                parts = line.split()
                if len(parts) >= len(col_names):
                    rec = {col_names[i]: parts[i] for i in range(len(col_names))}
                    if rec.get('group_PDB') != 'ATOM': continue
                    if rec.get('auth_asym_id', rec.get('label_asym_id', '')) != chain_id: continue
                    if rec.get('pdbx_PDB_model_num', '1') != '1': continue
                    atom_name = rec.get('auth_atom_id', rec.get('label_atom_id', ''))
                    if atom_name not in ('N', 'CA', 'C'): continue
                    try:
                        x, y, z = float(rec['Cartn_x']), float(rec['Cartn_y']), float(rec['Cartn_z'])
                        resnum = int(rec.get('auth_seq_id', 0))
                    except (ValueError, KeyError): continue
                    atoms.append({'atom': atom_name, 'resnum': resnum, 'x': x, 'y': y, 'z': z})
    return atoms


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    n1 = np.linalg.norm(b1)
    if n1 < 1e-10: return np.nan
    b1 = b1 / n1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    if np.linalg.norm(v) < 1e-10 or np.linalg.norm(w) < 1e-10: return np.nan
    return np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))


def compute_profile(atoms):
    res = {}
    for a in atoms:
        rn = a['resnum']
        if rn not in res: res[rn] = {}
        res[rn][a['atom']] = np.array([a['x'], a['y'], a['z']])
    sr = sorted(res.keys())
    phi_arr, psi_arr, rn_arr = [], [], []
    for i, r in enumerate(sr):
        d = res[r]
        if 'N' not in d or 'CA' not in d or 'C' not in d: continue
        pp = np.nan
        if i > 0 and 'C' in res.get(sr[i-1], {}):
            pp = dihedral(res[sr[i-1]]['C'], d['N'], d['CA'], d['C'])
        ps = np.nan
        if i < len(sr)-1 and 'N' in res.get(sr[i+1], {}):
            ps = dihedral(d['N'], d['CA'], d['C'], res[sr[i+1]]['N'])
        phi_arr.append(pp); psi_arr.append(ps); rn_arr.append(r)
    phi, psi = np.array(phi_arr), np.array(psi_arr)
    n = len(phi)
    ksq = np.zeros(n)
    for i in range(1, n-1):
        vals = [phi[i-1], phi[i], phi[i+1], psi[i-1], psi[i], psi[i+1]]
        if any(np.isnan(x) for x in vals): continue
        dp1 = np.arctan2(np.sin(phi[i]-phi[i-1]), np.cos(phi[i]-phi[i-1]))
        ds1 = np.arctan2(np.sin(psi[i]-psi[i-1]), np.cos(psi[i]-psi[i-1]))
        dp2 = np.arctan2(np.sin(phi[i+1]-phi[i]), np.cos(phi[i+1]-phi[i]))
        ds2 = np.arctan2(np.sin(psi[i+1]-psi[i]), np.cos(psi[i+1]-psi[i]))
        a1 = max(np.sqrt(dp1**2 + ds1**2), 1e-10)
        a2 = max(np.sqrt(dp2**2 + ds2**2), 1e-10)
        ksq[i] = ((dp2/a2 - dp1/a1)**2 + (ds2/a2 - ds1/a1)**2) / (0.5*(a1+a2))**2
    h_count = sum(1 for p, s in zip(phi, psi) if not np.isnan(p) and not np.isnan(s) and -160 < np.degrees(p) < -20 and -80 < np.degrees(s) < -10)
    return {'ksq': ksq, 'rn': rn_arr, 'n': n, 'h_frac': h_count/n if n > 0 else 0}


def gini(x):
    x = np.sort(np.abs(x)); n = len(x)
    if n == 0 or x.sum() == 0: return 0.0
    return float((2 * np.sum(np.arange(1, n+1) * x) - (n+1) * np.sum(x)) / (n * np.sum(x)))


def compare_pair(f1, c1, f2, c2, n_boot=2000):
    a1, a2 = parse_cif_chain(f1, c1), parse_cif_chain(f2, c2)
    nca1 = len([a for a in a1 if a['atom'] == 'CA'])
    nca2 = len([a for a in a2 if a['atom'] == 'CA'])
    if nca1 < 30 or nca2 < 30: return None
    p1, p2 = compute_profile(a1), compute_profile(a2)
    i1 = {r: i for i, r in enumerate(p1['rn'])}
    i2 = {r: i for i, r in enumerate(p2['rn'])}
    common = sorted(set(p1['rn']) & set(p2['rn']))
    if len(common) < 30: return None
    dk = np.array([p2['ksq'][i2[r]] - p1['ksq'][i1[r]] for r in common])
    k1c = np.array([p1['ksq'][i1[r]] for r in common])
    k2c = np.array([p2['ksq'][i2[r]] for r in common])
    mk = float(np.mean(np.abs(dk)))
    g1, g2 = gini(k1c), gini(k2c)
    dg = g2 - g1
    ksq1 = np.sum(p1['ksq']) / p1['n']
    ksq2 = np.sum(p2['ksq']) / p2['n']
    cv = float(np.std([ksq1, ksq2]) / np.mean([ksq1, ksq2]) * 100) if np.mean([ksq1, ksq2]) > 0 else 0
    np.random.seed(42)
    n = len(common)
    dg_boot = np.zeros(n_boot)
    for b in range(n_boot):
        idx = np.random.randint(0, n, n)
        dg_boot[b] = gini(k2c[idx]) - gini(k1c[idx])
    dg_ci = [float(np.percentile(dg_boot, 2.5)), float(np.percentile(dg_boot, 97.5))]
    dg_sig = (dg_ci[0] > 0 and dg_ci[1] > 0) or (dg_ci[0] < 0 and dg_ci[1] < 0)
    return {'mk': mk, 'dg': dg, 'cv': cv, 'g1': g1, 'g2': g2,
            'dg_ci': dg_ci, 'dg_sig': dg_sig, 'n_common': len(common), 'nca1': nca1, 'nca2': nca2}


# ══════════════════════════════════════════════════════════════
# DOWNLOAD
# ══════════════════════════════════════════════════════════════

def download_cif(pdb_id, cache_dir):
    import requests
    pdb_id = pdb_id.upper()
    filepath = cache_dir / f"{pdb_id}.cif"
    if filepath.exists() and filepath.stat().st_size > 1000: return filepath
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200 and len(r.content) > 500:
            filepath.write_bytes(r.content)
            return filepath
    except: pass
    return None


def batch_download(pdb_ids, cache_dir, label=""):
    cache_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for i, pid in enumerate(pdb_ids):
        pid = pid.upper()
        p = download_cif(pid, cache_dir)
        if p: paths[pid] = p
        if (i+1) % 10 == 0:
            print(f"  [{label}] {i+1}/{len(pdb_ids)} downloaded", flush=True)
            time.sleep(0.3)
    return paths


# ══════════════════════════════════════════════════════════════
# CHECKPOINT
# ══════════════════════════════════════════════════════════════

def load_checkpoint():
    if CHECKPOINT_FILE.exists():
        with open(CHECKPOINT_FILE) as f: return json.load(f)
    return {"temp": {}, "atp": {}}

def save_checkpoint(cp):
    with open(CHECKPOINT_FILE, 'w') as f: json.dump(cp, f, indent=2)


# ══════════════════════════════════════════════════════════════
# EXTENDED TEMPERATURE PAIRS
# ══════════════════════════════════════════════════════════════

TEMP_PAIRS_NEW = [
    # (label, pdb_cryo, chain, pdb_rt, chain, fold_class)
    ("Ribonuclease T1", "9RNT", "A", "1RGE", "A", "α+β"),
    ("Cytochrome c", "1CRC", "A", "1AKK", "A", "all-α"),
    ("BPTI", "4PTI", "A", "5PTI", "A", "small/disulfide"),
    ("Thermolysin", "8TLN", "A", "4TMN", "A", "metalloprotease"),
    ("Phospholipase A2", "1POA", "A", "1PSJ", "A", "α+β/Ca²⁺"),
    ("Chymotrypsinogen", "1CHG", "A", "2CGA", "A", "serine protease"),
    ("Aldose reductase", "1US0", "A", "2ACR", "A", "TIM barrel"),
    ("Cutinase", "1CEX", "A", "1CUS", "A", "α/β hydrolase"),
]


def run_temp(skip_download=False):
    print("\n" + "="*80)
    print("  EXTENDED TEMPERATURE PAIRS (8 new → 20 total)")
    print("="*80)

    cp = load_checkpoint()
    all_pdbs = set()
    for label, cryo, cc, rt, rc, fold in TEMP_PAIRS_NEW:
        all_pdbs.add(cryo.upper()); all_pdbs.add(rt.upper())

    if not skip_download:
        paths = batch_download(list(all_pdbs), CACHE_DIR, "Temp")
    else:
        paths = {p: CACHE_DIR / f"{p}.cif" for p in all_pdbs if (CACHE_DIR / f"{p}.cif").exists()}

    results = []
    for label, cryo_id, cryo_ch, rt_id, rt_ch, fold in TEMP_PAIRS_NEW:
        key = f"{cryo_id}_{rt_id}"
        if key in cp.get("temp", {}):
            r = cp["temp"][key]
            results.append(r)
            sig = "†" if r['dg_sig'] else " "
            d = "NEG" if r['dg'] < 0 else "POS"
            print(f"  [cached] {label:>22s} ({fold:>16s}): ΔGini={r['dg']:+.3f}{sig} {d}")
            continue

        f1 = CACHE_DIR / f"{cryo_id.upper()}.cif"
        f2 = CACHE_DIR / f"{rt_id.upper()}.cif"
        if not f1.exists() or not f2.exists():
            print(f"  {label:>22s}: MISSING ({cryo_id} or {rt_id})")
            continue

        r = compare_pair(str(f1), cryo_ch, str(f2), rt_ch)
        if r is None:
            # Try other chains
            for alt_ch in ['B', 'C', 'D']:
                r = compare_pair(str(f1), alt_ch, str(f2), alt_ch)
                if r is not None: break
        if r is None:
            print(f"  {label:>22s}: FAILED (insufficient residues)")
            continue

        r['label'] = label; r['cryo'] = cryo_id; r['rt'] = rt_id; r['fold'] = fold
        results.append(r)
        cp.setdefault("temp", {})[key] = r
        save_checkpoint(cp)

        sig = "†" if r['dg_sig'] else " "
        d = "NEG ✓" if r['dg'] < 0 else "POS ✗"
        print(f"  {label:>22s} ({fold:>16s}): ΔGini={r['dg']:+.3f}{sig} {d} |Δκ²|/r={r['mk']:.1f} n={r['n_common']}")

    # Combined summary with original 12
    print(f"\n  NEW PAIRS: {len(results)} computed")
    if results:
        n_neg = sum(1 for r in results if r['dg'] < 0)
        print(f"  Negative: {n_neg}/{len(results)}")
        serine = [r for r in results if 'serine' in r.get('fold','').lower() or r.get('label','') in ('Chymotrypsinogen',)]
        non_serine = [r for r in results if r not in serine]
        if serine:
            print(f"  Serine proteases: {sum(1 for r in serine if r['dg']<0)}/{len(serine)} negative")
        if non_serine:
            print(f"  Non-serine: {sum(1 for r in non_serine if r['dg']<0)}/{len(non_serine)} negative")

    return results


# ══════════════════════════════════════════════════════════════
# EXTENDED ATP SYNTHASE
# ══════════════════════════════════════════════════════════════

ATP_STRUCTURES = [
    # (label, pdb, source, description)
    # Original 3
    ("Ground state", "1BMF", "Abrahams 1994", "ADP + Pi + AMPPNP"),
    ("Transition (AlF4)", "1E1R", "Braig 2000", "ADP·AlF4 analog"),
    ("AMPPNP bound", "1H8E", "Menz 2001", "AMPPNP all sites"),
    # New
    ("ADP crystal form 2", "1E79", "Braig 2000", "ADP different form"),
    ("TNP-AMP bound", "1NBM", "Menz 2001", "Fluorescent analog"),
    ("Mixed occupancy", "1OHH", "Kagawa 2003", "AMPPNP + ADP-AlF4"),
    ("Pre-hydrolysis", "2JDI", "Bowler 2006", "ADP·BeF3"),
    ("Post-hydrolysis", "2JJ2", "Bowler 2007", "ADP + Pi"),
    ("Yeast F1", "2V7Q", "Kabaleeswaran 2006", "S. cerevisiae, cross-species"),
]

# β-subunit chain IDs (bovine F1: D=βDP, E=βTP, F=βE)
# α-subunit chain IDs (bovine F1: A, B, C) — negative control
BETA_CHAINS = ['D', 'E', 'F']
ALPHA_CHAINS = ['A', 'B', 'C']


def run_atp(skip_download=False):
    print(f"\n\n{'='*80}")
    print(f"  EXTENDED ATP SYNTHASE — {len(ATP_STRUCTURES)} structures")
    print(f"  β-subunits (catalytic) + α-subunits (negative control)")
    print(f"{'='*80}")

    cp = load_checkpoint()
    all_pdbs = [s[1].upper() for s in ATP_STRUCTURES]

    if not skip_download:
        paths = batch_download(all_pdbs, CACHE_DIR, "ATP")
    else:
        paths = {p: CACHE_DIR / f"{p}.cif" for p in all_pdbs if (CACHE_DIR / f"{p}.cif").exists()}

    beta_results = {}
    alpha_results = {}

    for label, pdb_id, source, desc in ATP_STRUCTURES:
        fp = CACHE_DIR / f"{pdb_id.upper()}.cif"
        if not fp.exists():
            print(f"  {pdb_id}: MISSING")
            continue

        state_map = {"D": "βDP(ADP)", "E": "βTP(ATP)", "F": "βE(empty)"}

        # β-subunits
        print(f"\n  {pdb_id} — {label} ({desc})")
        print(f"    β-subunits:")
        for ch in BETA_CHAINS:
            key = f"{pdb_id}_{ch}"
            if key in cp.get("atp", {}):
                r = cp["atp"][key]
                beta_results[key] = r
                st = state_map.get(ch, ch)
                print(f"      [cached] ch{ch} [{st:>12s}]: Gini={r['gini']:.3f}")
                continue

            atoms = parse_cif_chain(str(fp), ch)
            nca = len([a for a in atoms if a['atom'] == 'CA'])
            if nca < 100:
                print(f"      ch{ch}: too short ({nca} CA)")
                continue

            prof = compute_profile(atoms)
            g = gini(prof['ksq'])
            entry = {
                'pdb': pdb_id, 'chain': ch, 'label': label,
                'source': source, 'desc': desc,
                'n': prof['n'], 'gini': g,
                'mean_ksq': float(np.mean(prof['ksq'])),
                'h_frac': prof['h_frac'],
                'subunit': 'beta',
            }
            beta_results[key] = entry
            cp.setdefault("atp", {})[key] = entry
            save_checkpoint(cp)

            st = state_map.get(ch, ch)
            print(f"      ch{ch} [{st:>12s}]: Gini={g:.3f}, n={prof['n']}, mean_κ²={np.mean(prof['ksq']):.1f}")

        # α-subunits (negative control)
        print(f"    α-subunits (negative control):")
        for ch in ALPHA_CHAINS:
            key = f"{pdb_id}_alpha_{ch}"
            if key in cp.get("atp", {}):
                r = cp["atp"][key]
                alpha_results[key] = r
                print(f"      [cached] ch{ch}: Gini={r['gini']:.3f}")
                continue

            atoms = parse_cif_chain(str(fp), ch)
            nca = len([a for a in atoms if a['atom'] == 'CA'])
            if nca < 100:
                print(f"      ch{ch}: too short ({nca} CA)")
                continue

            prof = compute_profile(atoms)
            g = gini(prof['ksq'])
            entry = {
                'pdb': pdb_id, 'chain': ch, 'label': label,
                'n': prof['n'], 'gini': g,
                'mean_ksq': float(np.mean(prof['ksq'])),
                'h_frac': prof['h_frac'],
                'subunit': 'alpha',
            }
            alpha_results[key] = entry
            cp.setdefault("atp", {})[key] = entry
            save_checkpoint(cp)
            print(f"      ch{ch}: Gini={g:.3f}, n={prof['n']}")

    # ── ANALYSIS ──
    print(f"\n\n  {'═'*60}")
    print(f"  ATP SYNTHASE ANALYSIS")
    print(f"  {'═'*60}")

    # Group β by structure
    by_struct_beta = {}
    for key, r in beta_results.items():
        pid = r['pdb']
        if pid not in by_struct_beta: by_struct_beta[pid] = {}
        by_struct_beta[pid][r['chain']] = r

    print(f"\n  β-SUBUNIT GINI BY STRUCTURE AND STATE:")
    print(f"  {'PDB':>6s} {'Label':>25s} {'βDP(ADP)':>10s} {'βTP(ATP)':>10s} {'βE(empty)':>10s} {'Range':>8s} {'Gradient':>10s}")
    print(f"  {'─'*6} {'─'*25} {'─'*10} {'─'*10} {'─'*10} {'─'*8} {'─'*10}")

    for pid in sorted(by_struct_beta.keys()):
        chains = by_struct_beta[pid]
        label = chains[list(chains.keys())[0]].get('label', '')
        gD = chains.get('D', {}).get('gini', np.nan)
        gE = chains.get('E', {}).get('gini', np.nan)
        gF = chains.get('F', {}).get('gini', np.nan)
        vals = [v for v in [gD, gE, gF] if not np.isnan(v)]
        rng = max(vals) - min(vals) if len(vals) >= 2 else 0
        # Check if gradient is D < E < F (ADP < ATP < empty)
        gradient = ""
        if not np.isnan(gD) and not np.isnan(gE) and not np.isnan(gF):
            if gD < gE < gF: gradient = "D<E<F ✓"
            elif gD < gE: gradient = "D<E"
            elif gE < gF: gradient = "E<F"
            else: gradient = "—"
        gD_s = f"{gD:.3f}" if not np.isnan(gD) else "—"
        gE_s = f"{gE:.3f}" if not np.isnan(gE) else "—"
        gF_s = f"{gF:.3f}" if not np.isnan(gF) else "—"
        print(f"  {pid:>6s} {label:>25s} {gD_s:>10s} {gE_s:>10s} {gF_s:>10s} {rng:>8.3f} {gradient:>10s}")

    # α-subunit control
    by_struct_alpha = {}
    for key, r in alpha_results.items():
        pid = r['pdb']
        if pid not in by_struct_alpha: by_struct_alpha[pid] = {}
        by_struct_alpha[pid][r['chain']] = r

    print(f"\n  α-SUBUNIT NEGATIVE CONTROL (should show NO state discrimination):")
    print(f"  {'PDB':>6s} {'αA':>8s} {'αB':>8s} {'αC':>8s} {'Range':>8s}")
    print(f"  {'─'*6} {'─'*8} {'─'*8} {'─'*8} {'─'*8}")

    for pid in sorted(by_struct_alpha.keys()):
        chains = by_struct_alpha[pid]
        gA = chains.get('A', {}).get('gini', np.nan)
        gB = chains.get('B', {}).get('gini', np.nan)
        gC = chains.get('C', {}).get('gini', np.nan)
        vals = [v for v in [gA, gB, gC] if not np.isnan(v)]
        rng = max(vals) - min(vals) if len(vals) >= 2 else 0
        gA_s = f"{gA:.3f}" if not np.isnan(gA) else "—"
        gB_s = f"{gB:.3f}" if not np.isnan(gB) else "—"
        gC_s = f"{gC:.3f}" if not np.isnan(gC) else "—"
        print(f"  {pid:>6s} {gA_s:>8s} {gB_s:>8s} {gC_s:>8s} {rng:>8.3f}")

    # Compare β vs α ranges
    beta_ranges = []
    alpha_ranges = []
    for pid in sorted(set(by_struct_beta.keys()) & set(by_struct_alpha.keys())):
        bvals = [by_struct_beta[pid][ch]['gini'] for ch in by_struct_beta[pid]]
        avals = [by_struct_alpha[pid][ch]['gini'] for ch in by_struct_alpha[pid]]
        if len(bvals) >= 2: beta_ranges.append(max(bvals) - min(bvals))
        if len(avals) >= 2: alpha_ranges.append(max(avals) - min(avals))

    if beta_ranges and alpha_ranges:
        print(f"\n  β-subunit mean Gini range: {np.mean(beta_ranges):.3f} ± {np.std(beta_ranges):.3f}")
        print(f"  α-subunit mean Gini range: {np.mean(alpha_ranges):.3f} ± {np.std(alpha_ranges):.3f}")
        ratio = np.mean(beta_ranges) / max(np.mean(alpha_ranges), 0.001)
        print(f"  Ratio (β/α): {ratio:.1f}×")
        if ratio > 2:
            print(f"  → β-subunits show MORE state discrimination than α ✓")
        else:
            print(f"  → β and α show similar variation (negative control fails)")

    # Same state across structures
    print(f"\n  SAME STATE CONSISTENCY ACROSS STRUCTURES:")
    for ch in BETA_CHAINS:
        state = {"D": "βDP(ADP)", "E": "βTP(ATP)", "F": "βE(empty)"}[ch]
        vals = [(pid, by_struct_beta[pid][ch]['gini'])
                for pid in sorted(by_struct_beta.keys())
                if ch in by_struct_beta[pid]]
        if len(vals) >= 2:
            gs = [v[1] for v in vals]
            print(f"  {state:>12s}: n={len(vals)}, mean={np.mean(gs):.3f} ± {np.std(gs):.3f}")
            for pid, g in vals:
                print(f"    {pid}: {g:.3f}")

    return beta_results, alpha_results


# ══════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--test', choices=['temp', 'atp'])
    parser.add_argument('--skip-download', action='store_true')
    args = parser.parse_args()

    print("╔════════════════════════════════════════════════════════════╗")
    print("║  EXTENDED VALIDATION — Temperature + ATP Synthase        ║")
    print("╚════════════════════════════════════════════════════════════╝")

    results = {}

    if args.test is None or args.test == 'temp':
        results['temp'] = run_temp(skip_download=args.skip_download)

    if args.test is None or args.test == 'atp':
        beta, alpha = run_atp(skip_download=args.skip_download)
        results['atp_beta'] = {k: {kk: vv for kk, vv in v.items()} for k, v in beta.items()}
        results['atp_alpha'] = {k: {kk: vv for kk, vv in v.items()} for k, v in alpha.items()}

    # Save
    try:
        with open(RESULTS_FILE, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"\n  Results saved to {RESULTS_FILE}")
    except Exception as e:
        print(f"\n  WARNING: Could not save: {e}")
        print(f"  All results printed above.")


if __name__ == '__main__':
    main()

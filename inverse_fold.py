#!/usr/bin/env python3
"""
INVERSE FOLD — Map geometric barcodes to amino acid probability distributions.

Given a protein's (φ, ψ) trajectory on T², predict the most likely amino acid
at each position using a reference frequency table built from AlphaFold structures.

Usage:
    python inverse_fold.py build                  # Build reference table from cache
    python inverse_fold.py predict P00533         # Predict sequence for a barcoded protein
    python inverse_fold.py predict P00533 --top 3 # Show top-3 AAs per position
    python inverse_fold.py validate               # Validate predictions against real sequences
    python inverse_fold.py peak-decode            # Decode Gaussian peaks to residue identity
    python inverse_fold.py stats                  # Summary statistics of the reference table

The reference table maps 10°×10° Ramachandran bins to amino acid frequency
distributions. Each bin stores P(AA | φ_bin, ψ_bin) computed from all cached
AlphaFold structures.

Requires: numpy. Optional: scipy (for entropy calculations).
"""

import json
import math
import os
import sys
import time
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np

# ─── Paths ───
SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "data"
CACHE_DIRS = [
    SCRIPT_DIR / "alphafold_cache",
    DATA_DIR / "alphafold",
]
BARCODE_DIR = SCRIPT_DIR / "results" / "barcodes"
TABLE_PATH = DATA_DIR / "rama_aa_frequency_table.json"

# 3-letter → 1-letter mapping
AA3_TO_1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    # Non-standard → closest standard
    'MSE': 'M', 'SEC': 'C', 'SEP': 'S', 'TPO': 'T', 'PTR': 'Y',
    'HYP': 'P', 'CSO': 'C', 'MLY': 'K', 'M3L': 'K',
}

STANDARD_AA = 'ACDEFGHIKLMNPQRSTVWY'
BIN_SIZE = 10  # degrees
N_BINS = 360 // BIN_SIZE  # 36 bins per axis


def _dihedral(p1, p2, p3, p4):
    """Compute dihedral angle in degrees from four points."""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    if n1_norm < 1e-10 or n2_norm < 1e-10:
        return np.nan
    n1 /= n1_norm
    n2 /= n2_norm
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return np.degrees(np.arctan2(y, x))


def extract_dihedrals_and_residues(cif_path):
    """Extract (φ, ψ, residue_name) from mmCIF or PDB file.
    
    Returns:
        phi: array of phi angles (degrees)
        psi: array of psi angles (degrees)  
        res_names: list of 3-letter residue names
        res_ids: list of residue sequence numbers
    """
    path = Path(cif_path)
    is_cif = path.suffix.lower() == '.cif'
    
    atoms_by_residue = {}  # key -> {'N': xyz, 'CA': xyz, 'C': xyz}
    res_name_map = {}      # key -> 3-letter residue name
    target_chain = None
    
    if is_cif:
        # mmCIF parser
        in_atom_site = False
        col_names = []
        
        with open(path, 'r') as f:
            for line in f:
                line = line.rstrip()
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
                
                parts = line.split()
                if len(parts) < len(col_names):
                    continue
                
                record = {col_names[i]: parts[i] for i in range(len(col_names))}
                
                if record.get('group_PDB', '') not in ('ATOM', 'HETATM'):
                    continue
                
                atom_name = record.get('label_atom_id', record.get('auth_atom_id', ''))
                if atom_name not in ('N', 'CA', 'C'):
                    continue
                
                chain = record.get('auth_asym_id', record.get('label_asym_id', 'A'))
                if target_chain is None:
                    target_chain = chain
                if chain != target_chain:
                    continue
                
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
                    res_name = record.get('label_comp_id', record.get('auth_comp_id', 'UNK'))
                except (ValueError, KeyError):
                    continue
                
                key = (res_seq, icode)
                if key not in atoms_by_residue:
                    atoms_by_residue[key] = {}
                atoms_by_residue[key][atom_name] = np.array([x, y, z])
                res_name_map[key] = res_name.upper()
    else:
        # PDB parser
        with open(path, 'r') as f:
            for line in f:
                if not (line.startswith('ATOM') or line.startswith('HETATM')):
                    continue
                atom_name = line[12:16].strip()
                if atom_name not in ('N', 'CA', 'C'):
                    continue
                
                chain = line[21]
                if target_chain is None:
                    target_chain = chain
                if chain != target_chain:
                    continue
                
                alt = line[16].strip()
                if alt and alt not in ('A', '1'):
                    continue
                
                try:
                    res_seq = int(line[22:26])
                    icode = line[26].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    res_name = line[17:20].strip().upper()
                except (ValueError, IndexError):
                    continue
                
                key = (res_seq, icode)
                if key not in atoms_by_residue:
                    atoms_by_residue[key] = {}
                atoms_by_residue[key][atom_name] = np.array([x, y, z])
                res_name_map[key] = res_name
    
    if not atoms_by_residue:
        return None, None, None, None
    
    # Compute dihedrals
    sorted_keys = sorted(atoms_by_residue.keys())
    valid = []
    for key in sorted_keys:
        atoms = atoms_by_residue[key]
        if 'N' in atoms and 'CA' in atoms and 'C' in atoms:
            valid.append((key, atoms))
    
    n = len(valid)
    if n < 3:
        return None, None, None, None
    
    phi = np.full(n, np.nan)
    psi = np.full(n, np.nan)
    res_names = []
    res_ids = []
    
    for i, (key, atoms) in enumerate(valid):
        res_names.append(res_name_map.get(key, 'UNK'))
        res_ids.append(key[0])
        
        if i > 0:
            _, prev_atoms = valid[i - 1]
            if 'C' in prev_atoms:
                phi[i] = _dihedral(prev_atoms['C'], atoms['N'], atoms['CA'], atoms['C'])
        if i < n - 1:
            _, next_atoms = valid[i + 1]
            if 'N' in next_atoms:
                psi[i] = _dihedral(atoms['N'], atoms['CA'], atoms['C'], next_atoms['N'])
    
    return phi, psi, res_names, res_ids


def phi_psi_to_bin(phi_deg, psi_deg):
    """Convert (φ, ψ) in degrees to bin indices."""
    # Map from [-180, 180) to [0, 360)
    phi_shifted = (phi_deg + 180) % 360
    psi_shifted = (psi_deg + 180) % 360
    phi_bin = int(phi_shifted // BIN_SIZE)
    psi_bin = int(psi_shifted // BIN_SIZE)
    return min(phi_bin, N_BINS - 1), min(psi_bin, N_BINS - 1)


def bin_to_phi_psi(phi_bin, psi_bin):
    """Convert bin indices back to center (φ, ψ) in degrees."""
    phi_center = phi_bin * BIN_SIZE + BIN_SIZE / 2 - 180
    psi_center = psi_bin * BIN_SIZE + BIN_SIZE / 2 - 180
    return phi_center, psi_center


# ═══════════════════════════════════════════════════════════════════════
#  BUILD — Construct the (φ, ψ) → AA frequency table from cache
# ═══════════════════════════════════════════════════════════════════════

def build_table(max_structures=None, verbose=True):
    """Build the Ramachandran amino acid frequency table from cached structures.
    
    For each 10°×10° bin in Ramachandran space, compute P(AA | φ_bin, ψ_bin).
    """
    # Find all structures
    structure_files = []
    for d in CACHE_DIRS:
        if d.exists():
            structure_files.extend(d.glob("*.cif"))
            structure_files.extend(d.glob("*.pdb"))
    
    # Deduplicate
    seen = set()
    unique = []
    for f in sorted(structure_files):
        if f.stem not in seen:
            seen.add(f.stem)
            unique.append(f)
    structure_files = unique
    
    if max_structures:
        structure_files = structure_files[:max_structures]
    
    if verbose:
        print(f"  Building frequency table from {len(structure_files)} structures...")
    
    # Count matrix: [phi_bin][psi_bin][aa] = count
    counts = [[Counter() for _ in range(N_BINS)] for _ in range(N_BINS)]
    total_residues = 0
    n_processed = 0
    n_errors = 0
    aa_total = Counter()
    
    t0 = time.time()
    for i, sf in enumerate(structure_files):
        try:
            phi, psi, res_names, res_ids = extract_dihedrals_and_residues(sf)
            if phi is None:
                n_errors += 1
                continue
            
            for j in range(len(phi)):
                if np.isnan(phi[j]) or np.isnan(psi[j]):
                    continue
                
                rname = res_names[j]
                aa1 = AA3_TO_1.get(rname)
                if aa1 is None or aa1 not in STANDARD_AA:
                    continue
                
                pb, qb = phi_psi_to_bin(phi[j], psi[j])
                counts[pb][qb][aa1] += 1
                total_residues += 1
                aa_total[aa1] += 1
            
            n_processed += 1
        except Exception:
            n_errors += 1
        
        if verbose and (i + 1) % 500 == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            print(f"    {i+1}/{len(structure_files)} ({rate:.0f}/sec), "
                  f"{total_residues:,} residues, {n_errors} errors")
    
    elapsed = time.time() - t0
    if verbose:
        print(f"  Done: {n_processed} structures, {total_residues:,} residues, "
              f"{n_errors} errors, {elapsed:.1f}s")
    
    # Convert to frequency table
    table = {
        'bin_size': BIN_SIZE,
        'n_bins': N_BINS,
        'n_structures': n_processed,
        'n_residues': total_residues,
        'aa_totals': dict(aa_total),
        'bins': {}
    }
    
    n_occupied = 0
    for pb in range(N_BINS):
        for qb in range(N_BINS):
            if not counts[pb][qb]:
                continue
            n_occupied += 1
            total = sum(counts[pb][qb].values())
            freqs = {aa: cnt / total for aa, cnt in counts[pb][qb].items()}
            # Also store count for confidence
            key = f"{pb},{qb}"
            table['bins'][key] = {
                'total': total,
                'freqs': freqs,
            }
    
    table['n_occupied_bins'] = n_occupied
    
    # Compute per-bin entropy
    entropies = []
    for key, bdata in table['bins'].items():
        freqs = bdata['freqs']
        H = -sum(f * math.log2(f) for f in freqs.values() if f > 0)
        bdata['entropy'] = round(H, 3)
        entropies.append(H)
    
    if entropies:
        table['mean_entropy'] = round(sum(entropies) / len(entropies), 3)
        table['min_entropy'] = round(min(entropies), 3)
        table['max_entropy'] = round(max(entropies), 3)
    
    # Save
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(TABLE_PATH, 'w') as f:
        json.dump(table, f, indent=1)
    
    if verbose:
        print(f"  Saved: {TABLE_PATH}")
        print(f"  Occupied bins: {n_occupied}/{N_BINS*N_BINS} ({n_occupied/(N_BINS*N_BINS):.1%})")
        if entropies:
            print(f"  Entropy (bits): mean={table['mean_entropy']:.2f}, "
                  f"min={table['min_entropy']:.2f}, max={table['max_entropy']:.2f}")
            # Max possible entropy = log2(20) = 4.32 bits
            print(f"  Max possible entropy: {math.log2(20):.2f} bits (20 AAs)")
            low_entropy = sum(1 for e in entropies if e < 2.0)
            print(f"  Low-entropy bins (<2.0 bits): {low_entropy} "
                  f"({low_entropy/len(entropies):.1%}) — strong AA constraint")
    
    return table


def load_table():
    """Load the frequency table from disk."""
    if not TABLE_PATH.exists():
        print(f"  [ERROR] Frequency table not found: {TABLE_PATH}")
        print(f"  Run: python inverse_fold.py build")
        sys.exit(1)
    with open(TABLE_PATH) as f:
        return json.load(f)


# ═══════════════════════════════════════════════════════════════════════
#  PREDICT — Map a barcode to amino acid probabilities
# ═══════════════════════════════════════════════════════════════════════

def predict_sequence(uniprot_id, table=None, top_k=3, verbose=True):
    """Predict amino acid probabilities for each position of a protein.
    
    Loads the structure file, computes (φ, ψ), and maps each position
    to P(AA | φ_bin, ψ_bin) from the reference table.
    
    Returns list of dicts per residue:
        {res_id, res_name_real, phi, psi, top_predictions: [(aa, prob), ...],
         entropy, is_constrained, predicted_1}
    """
    if table is None:
        table = load_table()
    
    uid = uniprot_id.strip().upper()
    
    # Find structure file
    structure_file = None
    for d in CACHE_DIRS:
        if not d.exists():
            continue
        for pattern in [f"AF-{uid}*.cif", f"AF-{uid}*.pdb"]:
            matches = list(d.glob(pattern))
            if matches:
                structure_file = matches[0]
                break
        if structure_file:
            break
    
    if not structure_file:
        if verbose:
            print(f"  Structure not found for {uid}")
        return None
    
    # Extract dihedrals and residue names
    phi, psi, res_names, res_ids = extract_dihedrals_and_residues(structure_file)
    if phi is None:
        if verbose:
            print(f"  Failed to extract dihedrals from {structure_file}")
        return None
    
    # Map each position
    results = []
    n_correct_top1 = 0
    n_correct_top3 = 0
    n_predicted = 0
    
    for i in range(len(phi)):
        rname3 = res_names[i]
        rname1 = AA3_TO_1.get(rname3, '?')
        
        entry = {
            'res_id': res_ids[i],
            'res_name_real': rname1,
            'res_name_3': rname3,
            'phi': float(phi[i]) if not np.isnan(phi[i]) else None,
            'psi': float(psi[i]) if not np.isnan(psi[i]) else None,
            'top_predictions': [],
            'entropy': None,
            'is_constrained': False,
            'predicted_1': '?',
        }
        
        if np.isnan(phi[i]) or np.isnan(psi[i]):
            results.append(entry)
            continue
        
        pb, qb = phi_psi_to_bin(phi[i], psi[i])
        key = f"{pb},{qb}"
        
        bin_data = table['bins'].get(key)
        if bin_data is None:
            results.append(entry)
            continue
        
        freqs = bin_data['freqs']
        sorted_freqs = sorted(freqs.items(), key=lambda x: -x[1])
        
        entry['top_predictions'] = [(aa, round(p, 4)) for aa, p in sorted_freqs[:top_k]]
        entry['entropy'] = bin_data.get('entropy')
        entry['predicted_1'] = sorted_freqs[0][0] if sorted_freqs else '?'
        entry['is_constrained'] = bin_data.get('entropy', 5) < 2.5
        entry['bin_count'] = bin_data.get('total', 0)
        
        # Check accuracy
        n_predicted += 1
        if sorted_freqs and sorted_freqs[0][0] == rname1:
            n_correct_top1 += 1
        top_k_aas = [aa for aa, _ in sorted_freqs[:top_k]]
        if rname1 in top_k_aas:
            n_correct_top3 += 1
        
        results.append(entry)
    
    if verbose and n_predicted > 0:
        print(f"\n  INVERSE FOLD: {uid}")
        print(f"  {'='*60}")
        print(f"  Residues: {len(results)}, predicted: {n_predicted}")
        print(f"  Top-1 accuracy: {n_correct_top1}/{n_predicted} ({n_correct_top1/n_predicted:.1%})")
        print(f"  Top-{top_k} accuracy: {n_correct_top3}/{n_predicted} ({n_correct_top3/n_predicted:.1%})")
        
        # Entropy distribution
        entropies = [r['entropy'] for r in results if r['entropy'] is not None]
        if entropies:
            constrained = sum(1 for e in entropies if e < 2.5)
            print(f"  Constrained positions (H < 2.5 bits): {constrained}/{len(entropies)} "
                  f"({constrained/len(entropies):.1%})")
        
        # Show highly constrained positions
        print(f"\n  MOST CONSTRAINED POSITIONS (lowest entropy):")
        print(f"  {'Res':>5s} {'Real':>5s} {'Pred':>5s} {'Prob':>6s} {'H(bits)':>8s} "
              f"{'φ':>7s} {'ψ':>7s}  Top-{top_k}")
        print(f"  {'-'*5} {'-'*5} {'-'*5} {'-'*6} {'-'*8} {'-'*7} {'-'*7}  {'-'*20}")
        
        constrained_results = sorted(
            [r for r in results if r['entropy'] is not None],
            key=lambda x: x['entropy']
        )
        for r in constrained_results[:30]:
            preds = r['top_predictions']
            pred1 = preds[0][0] if preds else '?'
            prob1 = preds[0][1] if preds else 0
            match = "✓" if pred1 == r['res_name_real'] else " "
            top_str = ", ".join(f"{aa}:{p:.0%}" for aa, p in preds)
            print(f"  {r['res_id']:5d} {r['res_name_real']:>5s} {pred1:>5s} {prob1:5.0%} "
                  f"{r['entropy']:8.2f} {r['phi']:7.1f} {r['psi']:7.1f} {match} {top_str}")
        
        # Show Gly/Pro predictions
        print(f"\n  PRO/GLY DETECTION:")
        gly_pro = [r for r in results if r['res_name_real'] in ('G', 'P') and r['entropy'] is not None]
        if gly_pro:
            for r in gly_pro[:20]:
                preds = r['top_predictions']
                pred1 = preds[0][0] if preds else '?'
                top_str = ", ".join(f"{aa}:{p:.0%}" for aa, p in preds[:3])
                match = "✓" if r['res_name_real'] in [aa for aa, _ in preds[:3]] else "✗"
                print(f"    Res {r['res_id']:4d} = {r['res_name_real']} → {match} [{top_str}] H={r['entropy']:.2f}")
    
    return results


# ═══════════════════════════════════════════════════════════════════════
#  PEAK-DECODE — Map Gaussian peaks to residue identity
# ═══════════════════════════════════════════════════════════════════════

def peak_decode(table=None, max_proteins=200):
    """For each Gaussian peak in the barcode database, predict the residue identity.
    
    Tests whether peaks are Pro/Gly-driven (intrinsic kinks) vs other residues
    (environment-driven distortions).
    """
    if table is None:
        table = load_table()
    
    if not BARCODE_DIR.exists():
        print("  No barcode directory found.")
        return
    
    # Load barcodes
    barcodes = []
    for f in sorted(BARCODE_DIR.glob("*_barcode.json")):
        try:
            with open(f) as fh:
                barcodes.append(json.load(fh))
        except (json.JSONDecodeError, KeyError):
            pass
    
    if not barcodes:
        print("  No barcodes found.")
        return
    
    print(f"\n  GAUSSIAN PEAK RESIDUE DECODER")
    print(f"  {'='*60}")
    
    # Collect all Gaussian peaks with their peak_residue positions
    peaks = []
    for b in barcodes[:max_proteins]:
        uid = b.get('uniprot_id', '')
        for s in b.get('segments', []):
            if s.get('model_class') == 'gauss_peak':
                params = s.get('params', {})
                peak_res = params.get('peak_residue', -1)
                if peak_res >= 0:
                    peaks.append({
                        'uid': uid,
                        'ss_type': s['ss_type'],
                        'peak_residue': peak_res,
                        'A': params.get('A', 0),
                        'sigma': params.get('sigma', 0),
                        'start_idx': s['start_idx'],
                        'end_idx': s['end_idx'],
                    })
    
    print(f"  Gaussian peaks to decode: {len(peaks)}")
    if not peaks:
        return
    
    # Group by UniProt ID and batch-process
    by_uid = defaultdict(list)
    for pk in peaks:
        by_uid[pk['uid']].append(pk)
    
    # For each protein, load structure and map peak residues to AA predictions
    peak_residue_types = Counter()  # actual residue at peak
    peak_predicted_types = Counter()  # predicted residue at peak
    peak_results = []
    
    n_processed = 0
    for uid, uid_peaks in by_uid.items():
        # Find structure
        structure_file = None
        for d in CACHE_DIRS:
            if not d.exists():
                continue
            for pattern in [f"AF-{uid}*.cif", f"AF-{uid}*.pdb"]:
                matches = list(d.glob(pattern))
                if matches:
                    structure_file = matches[0]
                    break
            if structure_file:
                break
        
        if not structure_file:
            continue
        
        try:
            phi, psi, res_names, res_ids = extract_dihedrals_and_residues(structure_file)
            if phi is None:
                continue
        except Exception:
            continue
        
        # Map residue IDs to index
        resid_to_idx = {rid: i for i, rid in enumerate(res_ids)}
        
        for pk in uid_peaks:
            # peak_residue is 0-indexed position in the chain
            idx = pk['peak_residue']
            if idx >= len(res_names):
                continue
            
            rname3 = res_names[idx]
            rname1 = AA3_TO_1.get(rname3, '?')
            
            pk['real_aa'] = rname1
            pk['real_aa3'] = rname3
            peak_residue_types[rname1] += 1
            
            # Get prediction from table
            if idx < len(phi) and not np.isnan(phi[idx]) and not np.isnan(psi[idx]):
                pb, qb = phi_psi_to_bin(phi[idx], psi[idx])
                key = f"{pb},{qb}"
                bin_data = table['bins'].get(key)
                if bin_data:
                    freqs = bin_data['freqs']
                    sorted_freqs = sorted(freqs.items(), key=lambda x: -x[1])
                    pk['predicted_aa'] = sorted_freqs[0][0] if sorted_freqs else '?'
                    pk['top_3'] = sorted_freqs[:3]
                    pk['entropy'] = bin_data.get('entropy')
                    peak_predicted_types[pk['predicted_aa']] += 1
            
            peak_results.append(pk)
        
        n_processed += 1
        if n_processed % 100 == 0:
            print(f"    ...{n_processed}/{len(by_uid)} proteins")
    
    # ── Analysis ──
    print(f"\n  Decoded {len(peak_results)} peaks from {n_processed} proteins")
    
    # Real AA distribution at peak positions
    print(f"\n  REAL RESIDUE AT GAUSSIAN PEAK POSITIONS:")
    total_peaks = sum(peak_residue_types.values())
    for aa, cnt in peak_residue_types.most_common():
        bar = "█" * int(30 * cnt / total_peaks)
        print(f"    {aa}: {cnt:5d} ({cnt/total_peaks:5.1%}) {bar}")
    
    # Pro + Gly fraction
    n_pro = peak_residue_types.get('P', 0)
    n_gly = peak_residue_types.get('G', 0)
    n_pro_gly = n_pro + n_gly
    print(f"\n  Pro + Gly at peaks: {n_pro_gly}/{total_peaks} ({n_pro_gly/total_peaks:.1%})")
    print(f"    Pro: {n_pro} ({n_pro/total_peaks:.1%})")
    print(f"    Gly: {n_gly} ({n_gly/total_peaks:.1%})")
    print(f"  Expected from proteome: Pro ~5%, Gly ~7% → ~12%")
    
    if n_pro_gly / total_peaks > 0.15:
        print(f"  → Pro/Gly ENRICHED at Gaussian peaks ({n_pro_gly/total_peaks:.1%} vs ~12% baseline)")
        print(f"    Many peaks are intrinsic kinks (Pro) or flexibility points (Gly)")
    else:
        print(f"  → Pro/Gly NOT enriched — peaks primarily driven by environment, not sequence")
    
    # By SS type
    print(f"\n  PEAK RESIDUE BY SS TYPE:")
    for ss in ['H', 'E', 'C']:
        ss_peaks = [p for p in peak_results if p['ss_type'] == ss and 'real_aa' in p]
        if not ss_peaks:
            continue
        ss_aa = Counter(p['real_aa'] for p in ss_peaks)
        ss_name = {'H': 'Helix', 'E': 'Strand', 'C': 'Coil'}[ss]
        top5 = ss_aa.most_common(5)
        top_str = ", ".join(f"{aa}:{c}" for aa, c in top5)
        pg = ss_aa.get('P', 0) + ss_aa.get('G', 0)
        print(f"    {ss_name:8s} (n={len(ss_peaks):4d}): Pro+Gly={pg} ({pg/len(ss_peaks):.0%})  Top: {top_str}")
    
    # Save
    out_path = SCRIPT_DIR / "results" / "peak_decode_results.json"
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        summary = {
            'n_peaks_decoded': len(peak_results),
            'n_proteins': n_processed,
            'residue_distribution': dict(peak_residue_types),
            'pro_gly_fraction': n_pro_gly / total_peaks if total_peaks else 0,
            'pro_fraction': n_pro / total_peaks if total_peaks else 0,
            'gly_fraction': n_gly / total_peaks if total_peaks else 0,
        }
        with open(out_path, 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"\n  Saved: {out_path}")
    except Exception as e:
        print(f"  (Could not save: {e})")
    
    print()
    return peak_results


# ═══════════════════════════════════════════════════════════════════════
#  VALIDATE — Batch accuracy test across all barcoded proteins
# ═══════════════════════════════════════════════════════════════════════

def validate(table=None, max_proteins=100):
    """Validate inverse fold predictions against real sequences."""
    if table is None:
        table = load_table()
    
    if not BARCODE_DIR.exists():
        print("  No barcode directory found.")
        return
    
    barcodes = []
    for f in sorted(BARCODE_DIR.glob("*_barcode.json")):
        try:
            with open(f) as fh:
                barcodes.append(json.load(fh))
        except (json.JSONDecodeError, KeyError):
            pass
    
    print(f"\n  INVERSE FOLD VALIDATION")
    print(f"  {'='*60}")
    
    total_1 = 0
    total_3 = 0
    total_5 = 0
    total_n = 0
    per_protein = []
    
    for b in barcodes[:max_proteins]:
        uid = b.get('uniprot_id', '')
        if not uid:
            continue
        results = predict_sequence(uid, table=table, top_k=5, verbose=False)
        if not results:
            continue
        
        n = sum(1 for r in results if r['top_predictions'])
        c1 = sum(1 for r in results if r['top_predictions'] and r['top_predictions'][0][0] == r['res_name_real'])
        c3 = sum(1 for r in results if r['top_predictions'] and r['res_name_real'] in [a for a, _ in r['top_predictions'][:3]])
        c5 = sum(1 for r in results if r['top_predictions'] and r['res_name_real'] in [a for a, _ in r['top_predictions'][:5]])
        
        if n > 0:
            per_protein.append({'uid': uid, 'n': n, 'top1': c1/n, 'top3': c3/n, 'top5': c5/n})
            total_1 += c1
            total_3 += c3
            total_5 += c5
            total_n += n
    
    if total_n == 0:
        print("  No predictions made.")
        return
    
    print(f"  Proteins validated: {len(per_protein)}")
    print(f"  Total residues: {total_n:,}")
    print(f"\n  ACCURACY:")
    print(f"    Top-1: {total_1/total_n:.1%} ({total_1:,}/{total_n:,})")
    print(f"    Top-3: {total_3/total_n:.1%} ({total_3:,}/{total_n:,})")
    print(f"    Top-5: {total_5/total_n:.1%} ({total_5:,}/{total_n:,})")
    print(f"    Random baseline (top-1): {1/20:.1%}")
    print(f"    Random baseline (top-3): {3/20:.1%}")
    
    # Best and worst proteins
    per_protein.sort(key=lambda x: -x['top1'])
    print(f"\n  BEST PROTEINS (top-1 accuracy):")
    for p in per_protein[:5]:
        print(f"    {p['uid']}: {p['top1']:.1%} ({p['n']} residues)")
    print(f"  WORST PROTEINS:")
    for p in per_protein[-5:]:
        print(f"    {p['uid']}: {p['top1']:.1%} ({p['n']} residues)")
    
    print()


# ═══════════════════════════════════════════════════════════════════════
#  STATS — Summary of the frequency table
# ═══════════════════════════════════════════════════════════════════════

def table_stats(table=None):
    """Print summary statistics of the frequency table."""
    if table is None:
        table = load_table()
    
    print(f"\n  FREQUENCY TABLE STATISTICS")
    print(f"  {'='*60}")
    print(f"  Built from: {table['n_structures']} structures, {table['n_residues']:,} residues")
    print(f"  Bin size: {table['bin_size']}° × {table['bin_size']}°")
    print(f"  Occupied bins: {table['n_occupied_bins']}/{N_BINS*N_BINS} "
          f"({table['n_occupied_bins']/(N_BINS*N_BINS):.1%})")
    print(f"  Entropy: mean={table.get('mean_entropy', '?')}, "
          f"min={table.get('min_entropy', '?')}, max={table.get('max_entropy', '?')}")
    print(f"  Max possible: {math.log2(20):.2f} bits")
    
    # AA totals
    print(f"\n  AMINO ACID TOTALS:")
    aa_totals = table.get('aa_totals', {})
    total = sum(aa_totals.values())
    for aa in sorted(aa_totals, key=lambda x: -aa_totals[x]):
        cnt = aa_totals[aa]
        bar = "█" * int(30 * cnt / total)
        print(f"    {aa}: {cnt:>10,} ({cnt/total:5.1%}) {bar}")
    
    # Lowest entropy bins (most constrained regions of Ramachandran space)
    print(f"\n  MOST CONSTRAINED BINS (lowest entropy):")
    print(f"  {'φ center':>9s} {'ψ center':>9s} {'H(bits)':>8s} {'Count':>7s} {'Top-3 AA'}")
    print(f"  {'-'*9} {'-'*9} {'-'*8} {'-'*7} {'-'*30}")
    
    bin_list = []
    for key, bdata in table['bins'].items():
        pb, qb = [int(x) for x in key.split(',')]
        phi_c, psi_c = bin_to_phi_psi(pb, qb)
        bin_list.append((bdata.get('entropy', 5), phi_c, psi_c, bdata))
    
    bin_list.sort()
    for H, phi_c, psi_c, bdata in bin_list[:20]:
        freqs = bdata['freqs']
        top3 = sorted(freqs.items(), key=lambda x: -x[1])[:3]
        top_str = ", ".join(f"{aa}:{p:.0%}" for aa, p in top3)
        print(f"  {phi_c:+9.1f} {psi_c:+9.1f} {H:8.2f} {bdata['total']:7d} {top_str}")
    
    print()


# ─── Main ───
if __name__ == "__main__":
    args = sys.argv[1:]
    
    if not args:
        print(__doc__)
        sys.exit(0)
    
    cmd = args[0].lower()
    
    if cmd == 'build':
        max_s = int(args[1]) if len(args) > 1 else None
        build_table(max_structures=max_s)
    
    elif cmd == 'predict':
        if len(args) < 2:
            print("  Usage: inverse_fold.py predict <UniProt_ID> [--top N]")
            sys.exit(1)
        uid = args[1]
        top_k = 3
        for i, a in enumerate(args[2:]):
            if a == '--top' and i + 3 < len(args):
                top_k = int(args[i + 3])
        predict_sequence(uid, top_k=top_k)
    
    elif cmd == 'validate':
        max_p = int(args[1]) if len(args) > 1 else 100
        validate(max_proteins=max_p)
    
    elif cmd == 'peak-decode':
        max_p = int(args[1]) if len(args) > 1 else 200
        peak_decode(max_proteins=max_p)
    
    elif cmd == 'stats':
        table_stats()
    
    else:
        print(f"  Unknown command: {cmd}")
        print("  Commands: build, predict, validate, peak-decode, stats")

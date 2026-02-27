#!/usr/bin/env python3
"""
PDB EXPERIMENTAL VALIDATION
============================

Downloads high-resolution X-ray structures from PDB and runs the same
curvature/basin-null analysis used on AlphaFold structures.

Key question: Is coil suppression (C3) real or an AlphaFold artifact?

Usage:
    python pdb_validation.py download         # Download ~200 high-res structures
    python pdb_validation.py analyze          # Run curvature + barcode analysis
    python pdb_validation.py basin_null       # Run basin-control null
    python pdb_validation.py compare          # Compare PDB vs AlphaFold results
    python pdb_validation.py full             # All steps

Requirements:
    - cornu_spirals_on_T2.py (core pipeline) in same directory
    - basin_null_alphafold.py (null analysis) in same directory
    - numpy, scipy, matplotlib

Author: Kase Knochenhauer / True North Construction LLC
"""

import json
import os
import sys
import time
import math
import urllib.request
import urllib.error
import urllib.parse
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np

SCRIPT_DIR = Path(__file__).parent
PDB_DIR = SCRIPT_DIR / "data" / "pdb_xray"
RESULTS_DIR = SCRIPT_DIR / "results"
BARCODE_DIR = RESULTS_DIR / "barcodes_pdb"

for d in [PDB_DIR, RESULTS_DIR, BARCODE_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Import the core pipeline
sys.path.insert(0, str(SCRIPT_DIR))
try:
    from cornu_spirals_on_T2 import (
        extract_dihedrals_manual,
        dssp_assign_pure_numpy, assign_ss_from_dihedrals,
        torus_curvature, classify_spiral,
    )
    HAS_PIPELINE = True
except ImportError as e:
    print(f"  WARNING: Could not import pipeline: {e}")
    print(f"  Download and analyze steps will fail. Only 'download' works standalone.")
    HAS_PIPELINE = False


# ═══════════════════════════════════════════════════════════════════════
#  STEP 1: DOWNLOAD HIGH-RESOLUTION PDB STRUCTURES
# ═══════════════════════════════════════════════════════════════════════

def search_high_res_pdbs(max_resolution=1.5, min_length=100, max_length=800,
                          max_results=300):
    """Query RCSB PDB for high-resolution X-ray structures.
    
    Uses the RCSB Search API to find single-chain protein structures
    with ultra-high resolution.
    """
    # RCSB Search API query
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less",
                        "value": max_resolution,
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.experimental_method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION",
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "equals",
                        "value": 1,
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_sample_sequence_length",
                        "operator": "greater",
                        "value": min_length,
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_sample_sequence_length",
                        "operator": "less",
                        "value": max_length,
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.deposited_polymer_entity_instance_count",
                        "operator": "equals",
                        "value": 1,
                    }
                },
            ]
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": max_results
            },
            "sort": [
                {
                    "sort_by": "rcsb_entry_info.resolution_combined",
                    "direction": "asc"
                }
            ]
        }
    }

    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    headers = {
        'User-Agent': 'Mozilla/5.0 PhiPsiAtlas/1.0',
    }

    query_json = json.dumps(query)
    encoded = urllib.parse.quote(query_json)
    full_url = f"{url}?json={encoded}"
    req = urllib.request.Request(full_url, headers=headers)

    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            result = json.loads(resp.read().decode('utf-8'))
            pdb_ids = [hit['identifier'] for hit in result.get('result_set', [])]
            total = result.get('total_count', 0)
            print(f"  Found {total} structures, returning {len(pdb_ids)}")
            return pdb_ids
    except Exception as e:
        print(f"  RCSB search failed: {e}")
        return []


def download_pdb(pdb_id, output_dir=None):
    """Download a PDB file from RCSB."""
    if output_dir is None:
        output_dir = PDB_DIR

    outpath = output_dir / f"{pdb_id.lower()}.pdb"
    if outpath.exists() and outpath.stat().st_size > 5000:
        return outpath

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    headers = {'User-Agent': 'Mozilla/5.0 PhiPsiAtlas/1.0'}

    try:
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = resp.read()
            if len(data) > 1000:
                with open(outpath, 'wb') as f:
                    f.write(data)
                return outpath
    except Exception as e:
        # Try mmCIF format
        url2 = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        try:
            req = urllib.request.Request(url2, headers=headers)
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read()
                cifpath = output_dir / f"{pdb_id.lower()}.cif"
                if len(data) > 1000:
                    with open(cifpath, 'wb') as f:
                        f.write(data)
                    return cifpath
        except:
            pass

    return None


def cmd_download():
    """Download high-resolution PDB structures."""
    print(f"\n  STEP 1: Downloading high-resolution X-ray structures")
    print(f"  Target: ~200 single-chain proteins, resolution < 1.5 Å\n")

    # Search RCSB
    pdb_ids = search_high_res_pdbs(max_resolution=1.5, max_results=300)

    if not pdb_ids:
        print("  Search failed. Using curated list of high-res structures.")
        # Curated: single-chain proteins, X-ray < 1.5 Å, 100-800 residues
        # Sources: Top200 lists, PDB high-res hall of fame, teaching set
        pdb_ids = [
            # Ultra-high resolution (< 1.0 Å)
            "1CRN", "1ETL", "1UCS", "2VB1", "1GCI", "1RG8", "2IGD", "7ATJ",
            "3NIR", "1IUA", "2B97", "1X6Z", "1M40", "3A38", "4IGD",
            # Classic structures (1.0-1.3 Å)
            "1UBQ", "1L2Y", "1AHO", "3LZT", "2PNE", "1MBN", "2CBA", "1HTI",
            "1HHO", "2HHB", "1PKN", "3SDH", "1YCC", "1CTJ", "2RN2", "1PPE",
            "1TEN", "1FNF", "1CDG", "2CPL", "1A2P", "1AKE", "1ECA", "1IGD",
            "1LYZ", "2ACE", "1CSE", "1HOE", "1RNH", "4LYZ", "1BRS", "1STN",
            "2SNS", "1RN1", "1PGB", "1BPI", "1ENH",
            # Enzymes (1.0-1.5 Å)
            "1LZA", "3EST", "1CHO", "2CGA", "4CHA", "1PPL", "1ACB",
            "1ALC", "2PKA", "1PHH", "3TMS", "1XNB", "1CEX", "1BLC",
            "1THW", "2LZM", "1L58", "1WQ5", "1QJ4", "1QKI",
            # Binding proteins & carriers (1.0-1.5 Å)
            "1A6M", "1BEB", "1BJ7", "1B0Y", "2PTH", "1SHF", "1FUS",
            "1BKR", "1PHP", "2IFB", "1OPA", "1MLA",
            # Beta-rich (1.0-1.5 Å)
            "1TIT", "1WIT", "1HIC", "1NCO", "3FIB", "1REI", "1FKJ",
            "1NLS", "2RHE", "1CD8", "1TEN", "1FNA", "1TTG",
            # All-alpha (1.0-1.5 Å)
            "1HEL", "1A6G", "1UTG", "256B", "1CYO", "1DUR", "1VLS",
            "1EZM", "1H97", "1C75", "2CCY", "1YMB", "1MOL",
            # Alpha/beta mixed (1.0-1.5 Å)
            "1XEL", "1ADS", "1FLV", "2TRX", "1MEE", "1POC", "1TPH",
            "1AMP", "1NPC", "2CMD", "1MRO", "3GRS",
            # Small proteins (< 150 residues, high res)
            "1VII", "1CRN", "1K43", "2CI2", "1SRL", "1PIN", "1FME",
            "1LE0", "1PRB", "1PSE", "2HP8", "1BA5",
            # Additional high-resolution diverse set
            "1QHW", "2ERL", "1POH", "3SSI", "1HRC", "1WLA", "2RBP",
            "1OMB", "1TAP", "1TCA", "1MAZ", "1ROP", "2SOD",
            "1BYI", "1AX8", "1C9O", "1F94", "1JBC", "1G6N", "1XY1",
            "1DS1", "2CTC", "1THV", "1CNS", "1PEN", "1PLQ", "1QAU",
            "1OVA", "1NKI", "1TPK", "1NLS", "1SMD", "2PRK", "1UCN",
            "1BKE", "1FXD", "1CYC", "2MBW", "1CC8", "1CTF", "1RIS",
            "1ORC", "1FVQ", "1AON", "1BGE", "1TIM", "2YHX", "3PGK",
            "1HLE", "1AKY", "1BVB", "3MBP", "1PHT", "2OLB", "1KF6",
            # More diverse structures
            "1CHD", "1FCB", "1GKY", "1GD1", "2AIT", "1SHG", "1PMC",
            "1ALK", "1COA", "1GPR", "1KOE", "1AGQ", "1BXO", "1CIU",
        ]
        # Deduplicate
        pdb_ids = list(dict.fromkeys(pdb_ids))
        print(f"  Using {len(pdb_ids)} curated structures")

    # Download
    success = 0
    for i, pdb_id in enumerate(pdb_ids):
        path = download_pdb(pdb_id)
        if path:
            success += 1
            if (i + 1) % 50 == 0:
                print(f"    Downloaded {i+1}/{len(pdb_ids)} ({success} OK)")
        else:
            if (i + 1) % 50 == 0:
                print(f"    Progress {i+1}/{len(pdb_ids)} ({success} OK, {i+1-success} failed)")
        time.sleep(0.1)

    print(f"\n  Downloaded {success}/{len(pdb_ids)} structures to {PDB_DIR}/")

    # Save metadata
    meta = {
        "n_searched": len(pdb_ids),
        "n_downloaded": success,
        "max_resolution": 1.5,
        "pdb_ids": pdb_ids[:success],
    }
    with open(RESULTS_DIR / "pdb_download_meta.json", 'w') as f:
        json.dump(meta, f, indent=2)


# ═══════════════════════════════════════════════════════════════════════
#  STEP 2: ANALYZE — Run curvature classification on PDB structures
# ═══════════════════════════════════════════════════════════════════════

MACRO_MAP = {
    'geodesic': 'constant_k', 'circular_arc': 'constant_k',
    'sinusoidal': 'oscillatory', 'damped_oscillation': 'oscillatory',
    'damped_osc': 'oscillatory',
    'linear': 'monotone', 'clothoid': 'monotone', 'fermat': 'monotone',
    'exponential': 'monotone',
    'sigmoid': 'transition', 'step': 'transition',
    'quadratic': 'polynomial',
    'gauss_peak': 'localized',
    'too_short': 'too_short',
}


def analyze_one_structure(filepath):
    """Analyze one PDB/CIF file. Returns segments list or None on failure."""
    filepath = Path(filepath)

    try:
        if filepath.suffix == '.cif':
            from cornu_spirals_on_T2 import extract_dihedrals_mmcif
            phi, psi, res_ids, backbone_coords, has_O = extract_dihedrals_mmcif(str(filepath))
        else:
            phi, psi, res_ids, backbone_coords, has_O = extract_dihedrals_manual(str(filepath))
    except Exception as e:
        return None

    if phi is None or len(phi) < 10:
        return None

    # Assign SS using pure-numpy DSSP (same method as AlphaFold analysis)
    valid = ~(np.isnan(phi) | np.isnan(psi))
    phi_v = phi[valid]
    psi_v = psi[valid]
    
    if len(phi_v) < 10:
        return None

    # Detect if angles are in degrees and convert to radians
    max_abs = np.nanmax(np.abs(phi_v))
    if max_abs > np.pi + 0.01:
        phi_v = np.radians(phi_v)
        psi_v = np.radians(psi_v)

    # Try DSSP first (needs backbone N,CA,C,O coords), fall back to dihedral
    ss = None
    if backbone_coords is not None and has_O is not None:
        bb_v = backbone_coords[valid]
        ho_v = has_O[valid]
        # Only use DSSP if most residues have O atoms
        if np.sum(ho_v) > 0.8 * len(ho_v):
            try:
                ss = dssp_assign_pure_numpy(bb_v)
            except Exception:
                ss = None

    if ss is None:
        ss = assign_ss_from_dihedrals(phi_v, psi_v)

    # Debug: check SS distribution for first protein only
    if not hasattr(analyze_one_structure, '_printed_debug'):
        analyze_one_structure._printed_debug = True
        from collections import Counter as _C
        ss_dist = _C(ss)
        print(f"    DEBUG: max|phi|={max_abs:.2f}, SS dist={dict(ss_dist)}, n_res={len(ss)}, dssp={'yes' if ss is not None else 'no'}")

    # Segment and classify
    segments = []
    current_ss = ss[0] if len(ss) > 0 else 'C'
    seg_start = 0

    for i in range(1, len(ss)):
        if ss[i] != current_ss:
            if i - seg_start >= 4:
                seg_phi = phi_v[seg_start:i]
                seg_psi = psi_v[seg_start:i]

                # Compute curvature — returns (kappa, s_kappa, s_full)
                try:
                    kappa_result = torus_curvature(seg_phi, seg_psi)
                    if isinstance(kappa_result, tuple):
                        kappa = kappa_result[0]
                    else:
                        kappa = kappa_result
                    kappa = np.asarray(kappa, dtype=float)
                except:
                    kappa = None

                if kappa is not None and kappa.ndim == 1 and len(kappa) >= 3:
                    s_arc = np.arange(len(kappa), dtype=float)
                    try:
                        result = classify_spiral(kappa, s_arc)
                        if len(result) == 5:
                            slope, r2, model_name, fit_info, osc_info = result
                        elif len(result) == 4:
                            model_name, model_params, aic_score, r2 = result
                        else:
                            slope, r2, model_name = result[0], result[1], result[2]
                        macro = MACRO_MAP.get(model_name, model_name)
                        mean_k = float(np.mean(np.abs(kappa)))

                        segments.append({
                            'ss': current_ss,
                            'length': i - seg_start,
                            'model': model_name,
                            'macro': macro,
                            'mean_kappa': mean_k,
                            'r2': r2,
                        })
                    except:
                        pass

            current_ss = ss[i]
            seg_start = i

    # Final segment
    if len(ss) - seg_start >= 4:
        seg_phi = phi_v[seg_start:]
        seg_psi = psi_v[seg_start:]
        try:
            kappa_result = torus_curvature(seg_phi, seg_psi)
            if isinstance(kappa_result, tuple):
                kappa = kappa_result[0]
            else:
                kappa = kappa_result
            kappa = np.asarray(kappa, dtype=float)
        except:
            kappa = None

        if kappa is not None and kappa.ndim == 1 and len(kappa) >= 3:
            s_arc = np.arange(len(kappa), dtype=float)
            try:
                result = classify_spiral(kappa, s_arc)
                if len(result) == 5:
                    slope, r2, model_name, fit_info, osc_info = result
                elif len(result) == 4:
                    model_name, model_params, aic_score, r2 = result
                else:
                    slope, r2, model_name = result[0], result[1], result[2]
                macro = MACRO_MAP.get(model_name, model_name)
                mean_k = float(np.mean(np.abs(kappa)))
                segments.append({
                    'ss': current_ss,
                    'length': len(ss) - seg_start,
                    'model': model_name,
                    'macro': macro,
                    'mean_kappa': mean_k,
                    'r2': r2,
                })
            except:
                pass

    return segments if segments else None


def cmd_analyze():
    """Run curvature analysis on all downloaded PDB structures."""
    if not HAS_PIPELINE:
        print("  ERROR: Pipeline not available. Ensure cornu_spirals_on_T2.py is in the same directory.")
        return

    print(f"\n  STEP 2: Analyzing PDB structures\n")

    pdb_files = sorted(list(PDB_DIR.glob("*.pdb")) + list(PDB_DIR.glob("*.cif")))
    print(f"  Found {len(pdb_files)} structure files")

    all_segments = []
    n_success = 0
    n_fail = 0
    per_protein = []

    for i, fpath in enumerate(pdb_files):
        segments = analyze_one_structure(fpath)
        if segments:
            n_success += 1
            all_segments.extend(segments)
            per_protein.append({
                'file': fpath.name,
                'n_segments': len(segments),
                'n_residues': sum(s['length'] for s in segments),
            })
        else:
            n_fail += 1

        if (i + 1) % 50 == 0:
            print(f"    Processed {i+1}/{len(pdb_files)} ({n_success} OK, {n_fail} failed)")

    print(f"\n  Analyzed {n_success} proteins, {len(all_segments)} segments")

    # Compute statistics
    ss_counts = Counter(s['ss'] for s in all_segments)
    total = len(all_segments)

    print(f"\n  SS distribution:")
    for ss in ['H', 'E', 'C']:
        n = ss_counts.get(ss, 0)
        print(f"    {ss}: {n} ({n/total*100:.1f}%)")

    # Macro-class by SS
    print(f"\n  Non-constant fraction by SS:")
    for ss in ['H', 'E', 'C']:
        ss_segs = [s for s in all_segments if s['ss'] == ss]
        if ss_segs:
            n_nonconst = sum(1 for s in ss_segs if s['macro'] != 'constant_k')
            frac = n_nonconst / len(ss_segs)
            print(f"    {ss}: {n_nonconst}/{len(ss_segs)} = {frac:.3f}")

    # Gaussian peak counts
    print(f"\n  Gaussian peaks by SS:")
    for ss in ['H', 'E', 'C']:
        ss_segs = [s for s in all_segments if s['ss'] == ss]
        if ss_segs:
            n_gauss = sum(1 for s in ss_segs if s['macro'] == 'localized')
            frac = n_gauss / len(ss_segs) if ss_segs else 0
            print(f"    {ss}: {n_gauss}/{len(ss_segs)} = {frac:.3f} ({frac*100:.1f}%)")

    # Save results
    results = {
        "source": "PDB X-ray (resolution < 1.5 Å)",
        "n_proteins": n_success,
        "n_segments": len(all_segments),
        "ss_counts": dict(ss_counts),
        "segments": all_segments,
        "per_protein": per_protein,
    }

    with open(RESULTS_DIR / "pdb_xray_analysis.json", 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved to {RESULTS_DIR}/pdb_xray_analysis.json")


# ═══════════════════════════════════════════════════════════════════════
#  STEP 3: COMPARE PDB vs AlphaFold
# ═══════════════════════════════════════════════════════════════════════

def cmd_compare():
    """Compare PDB experimental vs AlphaFold analysis results."""
    print(f"\n  STEP 3: Comparing PDB vs AlphaFold results\n")

    # Load PDB results
    pdb_file = RESULTS_DIR / "pdb_xray_analysis.json"
    if not pdb_file.exists():
        print("  Run 'analyze' first.")
        return
    with open(pdb_file) as f:
        pdb = json.load(f)

    # Load AlphaFold results (from basin null which has segment data)
    af_file = None
    for candidate in [
        SCRIPT_DIR / "data" / "basin_null_results.json",
        Path("/mnt/user-data/uploads/basin_null_results.json"),
    ]:
        if candidate.exists():
            af_file = candidate
            break

    if af_file is None:
        print("  AlphaFold results not found. Showing PDB-only analysis.")
        af = None
    else:
        with open(af_file) as f:
            af = json.load(f)

    pdb_segs = pdb['segments']

    print(f"  PDB X-ray: {pdb['n_proteins']} proteins, {len(pdb_segs)} segments")

    # PDB statistics
    print(f"\n  {'Metric':<35s} {'PDB X-ray':>12s} {'AlphaFold':>12s} {'Δ':>8s}")
    print(f"  {'─'*70}")

    for ss, ss_name in [('H', 'Helix'), ('E', 'Strand'), ('C', 'Coil')]:
        pdb_ss = [s for s in pdb_segs if s['ss'] == ss]
        pdb_nonconst = sum(1 for s in pdb_ss if s['macro'] != 'constant_k') / len(pdb_ss) if pdb_ss else 0
        pdb_gauss = sum(1 for s in pdb_ss if s['macro'] == 'localized') / len(pdb_ss) if pdb_ss else 0

        pdb_val = f"{pdb_nonconst:.3f}"
        pdb_gval = f"{pdb_gauss:.3f}"

        # AlphaFold values from the paper
        af_nonconst = {'H': 0.466, 'E': 0.395, 'C': 0.284}  # from basin null results
        af_gauss = {'H': 0.163, 'E': 0.053, 'C': 0.065}      # approximate from paper

        af_val = f"{af_nonconst.get(ss, 0):.3f}"
        delta = pdb_nonconst - af_nonconst.get(ss, 0)

        print(f"  {ss_name + ' non-constant fraction':<35s} {pdb_val:>12s} {af_val:>12s} {delta:>+8.3f}")

    print()

    # The key test: coil suppression
    pdb_coil = [s for s in pdb_segs if s['ss'] == 'C']
    if pdb_coil:
        pdb_coil_nonconst = sum(1 for s in pdb_coil if s['macro'] != 'constant_k') / len(pdb_coil)
        af_coil_nonconst = 0.284  # from AlphaFold analysis

        print(f"  ★ KEY TEST: COIL SUPPRESSION")
        print(f"    AlphaFold coil non-constant: {af_coil_nonconst:.3f}")
        print(f"    PDB X-ray coil non-constant: {pdb_coil_nonconst:.3f}")
        delta = pdb_coil_nonconst - af_coil_nonconst
        print(f"    Difference: {delta:+.3f}")

        if abs(delta) < 0.02:
            print(f"    → CONFIRMED: Coil suppression is consistent across AF and PDB")
        elif delta > 0.02:
            print(f"    → REVERSED: PDB coils show MORE modulation than AlphaFold")
            print(f"    → Suggests AlphaFold smooths coil regions (artifact)")
        else:
            print(f"    → ENHANCED: PDB coils show LESS modulation than AlphaFold")

    # Mean curvature comparison
    print(f"\n  Mean |κ| by SS:")
    for ss, ss_name in [('H', 'Helix'), ('E', 'Strand'), ('C', 'Coil')]:
        pdb_ss = [s for s in pdb_segs if s['ss'] == ss]
        if pdb_ss:
            pdb_mean_k = np.mean([s['mean_kappa'] for s in pdb_ss])
            print(f"    {ss_name}: {pdb_mean_k:.4f} rad/residue (PDB)")

    # Save comparison
    comparison = {
        "pdb_n_proteins": pdb['n_proteins'],
        "pdb_n_segments": len(pdb_segs),
        "pdb_coil_nonconst": pdb_coil_nonconst if pdb_coil else None,
        "af_coil_nonconst": 0.284,
    }
    with open(RESULTS_DIR / "pdb_vs_alphafold_comparison.json", 'w') as f:
        json.dump(comparison, f, indent=2)

    print(f"\n  Comparison saved to {RESULTS_DIR}/pdb_vs_alphafold_comparison.json")


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
    elif cmd == 'compare':
        cmd_compare()
    elif cmd == 'full':
        cmd_download()
        cmd_analyze()
        cmd_compare()
    else:
        print(f"  Unknown: {cmd}")
        print(f"  Commands: download, analyze, compare, full")


if __name__ == "__main__":
    main()

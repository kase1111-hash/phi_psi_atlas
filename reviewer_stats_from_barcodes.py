#!/usr/bin/env python3
"""
REVIEWER STATS FROM BARCODES
==============================

Reads existing barcode JSON files from results/barcodes/ and computes
the three reviewer-requested statistics:

  #7  Median segment length by SS + ΔAICc (if available)
  #8  χ² for Pro+Gly at peaks (from paper values)
  #9  Helix Gaussian peak pair count for exclusion zone

Usage:
    python reviewer_stats_from_barcodes.py

Expects: results/barcodes/*_barcode.json (from cornu_spirals_on_T2.py)
"""

import json
import math
import sys
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
from scipy import stats

SCRIPT_DIR = Path(__file__).parent
BARCODE_DIR = SCRIPT_DIR / "results" / "barcodes"

def load_all_barcodes():
    """Load all barcode JSONs."""
    barcodes = []
    if not BARCODE_DIR.exists():
        print(f"  ⚠ Barcodes directory not found: {BARCODE_DIR}")
        return barcodes
    for f in sorted(BARCODE_DIR.glob("*_barcode.json")):
        try:
            with open(f) as fh:
                bc = json.load(fh)
                bc['_file'] = f.name
                barcodes.append(bc)
        except Exception as e:
            pass
    return barcodes


def main():
    print("\n  REVIEWER STATISTICS FROM BARCODE DATA")
    print("  " + "═" * 55)
    
    barcodes = load_all_barcodes()
    print(f"\n  Loaded {len(barcodes)} protein barcodes")
    
    if not barcodes:
        print("  No barcodes found. Run with paper values only.\n")
        run_chi2_from_paper()
        return
    
    # Collect all segments
    all_segments = []
    for bc in barcodes:
        uid = bc.get('uniprot_id', bc.get('_file', '?'))
        for seg in bc.get('segments', []):
            seg['_protein'] = uid
            all_segments.append(seg)
    
    print(f"  Total segments: {len(all_segments)}")
    
    # Show sample segment to confirm field names
    if all_segments:
        print(f"  Sample segment keys: {list(all_segments[0].keys())}")
    
    # ═══════════════════════════════════════════════════════════════
    #  #7a: SEGMENT LENGTH DISTRIBUTION
    # ═══════════════════════════════════════════════════════════════
    print(f"\n  #7a: SEGMENT LENGTH DISTRIBUTION")
    print("  " + "─" * 50)
    
    for ss, ss_name in [('H', 'Helix'), ('E', 'Strand'), ('C', 'Coil')]:
        lengths = [s['length'] for s in all_segments 
                   if s.get('ss_type', '') == ss and 'length' in s]
        if lengths:
            lengths = np.array(lengths)
            print(f"\n    {ss_name} (n={len(lengths)}):")
            print(f"      Median: {np.median(lengths):.0f} residues")
            print(f"      Mean:   {np.mean(lengths):.1f} ± {np.std(lengths):.1f}")
            print(f"      Range:  {int(np.min(lengths))}–{int(np.max(lengths))}")
            print(f"      Q1/Q3:  {np.percentile(lengths,25):.0f}/{np.percentile(lengths,75):.0f}")
            short = np.sum(lengths <= 5)
            print(f"      ≤5 res: {short} ({short/len(lengths)*100:.1f}%)")
    
    # ═══════════════════════════════════════════════════════════════
    #  #7b: ΔAICc (check if available)
    # ═══════════════════════════════════════════════════════════════
    print(f"\n  #7b: ΔAICc DISTRIBUTION")
    print("  " + "─" * 50)
    
    # Check various possible keys
    aicc_keys = ['delta_aicc', 'daicc', 'aicc_delta', 'all_aicc', 
                 'aicc_values', 'aicc_best', 'aicc_second']
    found_key = None
    for key in aicc_keys:
        if any(key in s for s in all_segments[:20]):
            found_key = key
            break
    
    if found_key:
        print(f"    Found AICc data under key: '{found_key}'")
        deltas = []
        for s in all_segments:
            val = s.get(found_key, None)
            if val is not None:
                try:
                    deltas.append(float(val))
                except (ValueError, TypeError):
                    pass
        
        if deltas:
            deltas = np.array(deltas)
            print(f"    N segments with ΔAICc: {len(deltas)}")
            print(f"    Median ΔAICc: {np.median(deltas):.1f}")
            print(f"    Mean ΔAICc:   {np.mean(deltas):.1f} ± {np.std(deltas):.1f}")
            print(f"    ΔAICc > 10:   {np.sum(deltas > 10)} ({np.sum(deltas > 10)/len(deltas)*100:.1f}%)")
            print(f"    ΔAICc > 4:    {np.sum(deltas > 4)} ({np.sum(deltas > 4)/len(deltas)*100:.1f}%) — strong selection")
            print(f"    ΔAICc > 2:    {np.sum(deltas > 2)} ({np.sum(deltas > 2)/len(deltas)*100:.1f}%) — positive selection")
            print(f"    ΔAICc < 2:    {np.sum(deltas < 2)} ({np.sum(deltas < 2)/len(deltas)*100:.1f}%) — weak/ambiguous")
            
            # By SS
            for ss, ss_name in [('H', 'Helix'), ('E', 'Strand'), ('C', 'Coil')]:
                ss_deltas = [float(s[found_key]) for s in all_segments 
                             if s.get('ss_type','') == ss and s.get(found_key) is not None]
                if ss_deltas:
                    ss_d = np.array(ss_deltas)
                    print(f"    {ss_name}: median ΔAICc = {np.median(ss_d):.1f}, "
                          f">{4}: {np.sum(ss_d>4)/len(ss_d)*100:.0f}%")
    else:
        print(f"    ⚠ No ΔAICc stored in barcode JSONs.")
        print(f"    The current pipeline selects the best model by AICc")
        print(f"    but does not store the runner-up AICc value.")
        print(f"    To add this, modify cornu_spirals_on_T2.py to record")
        print(f"    'delta_aicc' = second_best_aicc - best_aicc per segment.")
        print(f"\n    This is a revision-round item. For initial submission,")
        print(f"    state: 'Model selection uses AICc with penalty for")
        print(f"    parameter count; ΔAICc distribution to be reported in")
        print(f"    supplementary materials.'")
    
    # ═══════════════════════════════════════════════════════════════
    #  #8: χ² TEST FOR PRO+GLY
    # ═══════════════════════════════════════════════════════════════
    run_chi2_from_paper()
    
    # ═══════════════════════════════════════════════════════════════
    #  #9: HELIX GAUSSIAN PEAK PAIRS
    # ═══════════════════════════════════════════════════════════════
    print(f"\n  #9: HELIX GAUSSIAN PEAK PAIRS (EXCLUSION ZONE)")
    print("  " + "─" * 50)
    
    # First, discover what model_class values exist for helix segments
    helix_classes = Counter()
    for seg in all_segments:
        if seg.get('ss_type', '') == 'H':
            helix_classes[seg.get('model_class', '?')] += 1
    
    print(f"\n    Helix model_class values found:")
    for cls, count in helix_classes.most_common():
        print(f"      {cls}: {count}")
    
    # Identify which class is the Gaussian peak
    peak_names = [c for c in helix_classes 
                  if any(kw in c.lower() for kw in ['gauss', 'peak', 'localized', 'gaussian'])]
    if not peak_names:
        print(f"\n    ⚠ Could not identify Gaussian peak class.")
        print(f"    Trying all classes containing 'peak' or 'gauss'...")
        peak_names = [c for c in helix_classes if 'peak' in c.lower() or 'gauss' in c.lower()]
    
    if peak_names:
        print(f"    Using peak class(es): {peak_names}")
    
    total_pairs = 0
    total_peaks = 0
    proteins_with_peaks = 0
    spacing_values = []
    
    for bc in barcodes:
        segs = bc.get('segments', [])
        
        # Find helix Gaussian peaks with their positions
        helix_peaks = []
        for seg in segs:
            if (seg.get('ss_type', '') == 'H' and 
                seg.get('model_class', '') in peak_names):
                # Try various position keys
                pos = None
                for key in ['peak_residue', 'start_idx', 'seg_start', 
                           'start_residue', 'start']:
                    if key in seg:
                        pos = seg[key]
                        break
                if pos is not None:
                    helix_peaks.append(int(pos))
        
        if len(helix_peaks) >= 2:
            proteins_with_peaks += 1
            total_peaks += len(helix_peaks)
            helix_peaks.sort()
            pairs = len(helix_peaks) - 1
            total_pairs += pairs
            
            # Record spacings
            for i in range(len(helix_peaks) - 1):
                spacing = helix_peaks[i+1] - helix_peaks[i]
                if spacing > 0:
                    spacing_values.append(spacing)
    
    print(f"\n    Proteins with ≥2 helix Gaussian peaks: {proteins_with_peaks}")
    print(f"    Total helix Gaussian peaks: {total_peaks}")
    print(f"    Total consecutive pairs: {total_pairs}")
    
    if spacing_values:
        sp = np.array(spacing_values)
        print(f"\n    Spacing statistics:")
        print(f"      Median: {np.median(sp):.0f} residues")
        print(f"      Mean:   {np.mean(sp):.1f} ± {np.std(sp):.1f}")
        excl = np.sum(sp < 10)
        print(f"      Pairs < 10 residues apart: {excl} ({excl/len(sp)*100:.1f}%)")
    
    print(f"\n    FOR PAPER: 'Exclusion zone analysis based on {total_pairs}")
    print(f"    consecutive helix Gaussian peak pairs from")
    print(f"    {proteins_with_peaks} proteins ({total_peaks} peaks total).'")
    
    # ═══════════════════════════════════════════════════════════════
    #  SUMMARY
    # ═══════════════════════════════════════════════════════════════
    print(f"\n  " + "═" * 55)
    print(f"  PAPER-READY SENTENCES:")
    print(f"  " + "═" * 55)
    
    # Segment lengths
    for ss, ss_name in [('H', 'Helix'), ('E', 'Strand'), ('C', 'Coil')]:
        lengths = [s['length'] for s in all_segments 
                   if s.get('ss_type', '') == ss and 'length' in s]
        if lengths:
            med = np.median(lengths)
            print(f"\n    {ss_name}: median segment length = {med:.0f} residues "
                  f"(n = {len(lengths)})")
    
    print(f"\n    Pro+Gly: χ²(1) = 1.3, p = 0.26, 95% CI [8.7%, 13.1%]")
    print(f"\n    Exclusion zone: {total_pairs} helix peak pairs from"
          f" {proteins_with_peaks} proteins")


def run_chi2_from_paper():
    """Run χ² from paper values."""
    print(f"\n  #8: χ² TEST FOR PRO+GLY AT GAUSSIAN PEAKS")
    print("  " + "─" * 50)
    
    n_peaks = 738
    observed_freq = 0.108
    n_pro_gly = round(n_peaks * observed_freq)
    baseline_freq = 0.122
    
    expected_pg = n_peaks * baseline_freq
    expected_other = n_peaks * (1 - baseline_freq)
    observed_other = n_peaks - n_pro_gly
    
    chi2 = ((n_pro_gly - expected_pg)**2 / expected_pg + 
            (observed_other - expected_other)**2 / expected_other)
    p_value = stats.chi2.sf(chi2, df=1)
    
    ci_lo, ci_hi = stats.binom.interval(0.95, n_peaks, observed_freq)
    ci_lo /= n_peaks
    ci_hi /= n_peaks
    
    print(f"\n    n = {n_peaks} peaks")
    print(f"    Observed Pro+Gly: {observed_freq*100:.1f}%")
    print(f"    Baseline:         {baseline_freq*100:.1f}%")
    print(f"    χ²(1) = {chi2:.2f}, p = {p_value:.4f}")
    print(f"    95% CI: [{ci_lo*100:.1f}%, {ci_hi*100:.1f}%]")
    print(f"    → NOT SIGNIFICANT")


if __name__ == "__main__":
    main()

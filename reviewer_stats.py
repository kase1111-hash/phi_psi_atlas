#!/usr/bin/env python3
"""
REVIEWER STATISTICS
====================

Computes statistics requested by reviewers from the main analysis data:

  #7  Median segment length + ΔAICc distribution  
  #8  χ² test for Pro+Gly enrichment at Gaussian peaks
  #9  Number of helix-helix peak pairs for exclusion zone

Requires the main 557-protein analysis JSON and peak-decode data.

Usage:
    python reviewer_stats.py <analysis_json> [peak_decode_json]
    
    analysis_json:   The main output from cornu_spirals_on_T2.py
    peak_decode_json: Optional, peak amino acid identity data

If your data is in a different format, adapt the loading section.

Author: Kase Knochenhauer / True North Construction LLC
"""

import json
import sys
import numpy as np
from collections import Counter
from pathlib import Path


def load_analysis(path):
    """Load the main analysis JSON. Adapt keys as needed."""
    with open(path) as f:
        data = json.load(f)
    return data


# ═══════════════════════════════════════════════════════════════════════
#  #7: SEGMENT LENGTHS + ΔAICc
# ═══════════════════════════════════════════════════════════════════════

def stat_segment_lengths(segments):
    """Report segment length statistics by SS type."""
    print("\n  #7a: SEGMENT LENGTH DISTRIBUTION")
    print("  " + "─" * 50)
    
    for ss in ['H', 'E', 'C']:
        ss_name = {'H': 'Helix', 'E': 'Strand', 'C': 'Coil'}[ss]
        lengths = []
        for seg in segments:
            seg_ss = seg.get('ss', seg.get('ss_type', ''))
            if seg_ss == ss:
                L = seg.get('length', seg.get('n_residues', seg.get('len', 0)))
                if L > 0:
                    lengths.append(L)
        
        if lengths:
            lengths = np.array(lengths)
            print(f"\n    {ss_name} (n={len(lengths)}):")
            print(f"      Median: {np.median(lengths):.0f} residues")
            print(f"      Mean:   {np.mean(lengths):.1f} ± {np.std(lengths):.1f}")
            print(f"      Range:  {np.min(lengths)}–{np.max(lengths)}")
            print(f"      Q1/Q3:  {np.percentile(lengths,25):.0f}/{np.percentile(lengths,75):.0f}")
            # Fraction that are 4-5 residues (minimal fits)
            short = np.sum(lengths <= 5)
            print(f"      ≤5 res: {short} ({short/len(lengths)*100:.1f}%)")


def stat_delta_aicc(segments):
    """Report ΔAICc between best and second-best model."""
    print("\n  #7b: ΔAICc DISTRIBUTION (best vs second-best)")
    print("  " + "─" * 50)
    
    deltas = []
    for seg in segments:
        # Look for AICc values - adapt key names to your data
        aicc_vals = seg.get('aicc_values', seg.get('all_aicc', None))
        
        if aicc_vals is not None:
            if isinstance(aicc_vals, dict):
                vals = sorted(aicc_vals.values())
            elif isinstance(aicc_vals, list):
                vals = sorted(aicc_vals)
            else:
                continue
            
            if len(vals) >= 2:
                deltas.append(vals[1] - vals[0])  # 2nd best - best
        
        # Alternative: look for delta_aicc directly
        da = seg.get('delta_aicc', seg.get('daicc', None))
        if da is not None and aicc_vals is None:
            deltas.append(float(da))
    
    if deltas:
        deltas = np.array(deltas)
        print(f"\n    N segments with AICc data: {len(deltas)}")
        print(f"    Median ΔAICc: {np.median(deltas):.1f}")
        print(f"    Mean ΔAICc:   {np.mean(deltas):.1f} ± {np.std(deltas):.1f}")
        print(f"    ΔAICc > 4:    {np.sum(deltas > 4)} ({np.sum(deltas > 4)/len(deltas)*100:.1f}%) — strong selection")
        print(f"    ΔAICc > 2:    {np.sum(deltas > 2)} ({np.sum(deltas > 2)/len(deltas)*100:.1f}%) — positive selection")
        print(f"    ΔAICc < 2:    {np.sum(deltas < 2)} ({np.sum(deltas < 2)/len(deltas)*100:.1f}%) — weak/ambiguous")
        
        if np.median(deltas) > 4:
            print(f"\n    ★ Median ΔAICc > 4: MODEL SELECTION IS STRONG")
        elif np.median(deltas) > 2:
            print(f"\n    ★ Median ΔAICc > 2: Model selection is adequate")
        else:
            print(f"\n    ⚠ Median ΔAICc < 2: Classification uncertainty exists")
    else:
        print("\n    ⚠ No AICc data found in segments.")
        print("    Expected keys: 'aicc_values', 'all_aicc', or 'delta_aicc'")
        print("    Your segment keys:", list(segments[0].keys()) if segments else "none")


# ═══════════════════════════════════════════════════════════════════════
#  #8: χ² TEST FOR PRO+GLY AT PEAKS
# ═══════════════════════════════════════════════════════════════════════

def stat_peak_chi2(segments=None, peak_decode_path=None):
    """χ² test for Pro+Gly enrichment at Gaussian peaks."""
    print("\n  #8: χ² TEST FOR PRO+GLY AT GAUSSIAN PEAKS")
    print("  " + "─" * 50)
    
    peak_aas = []
    
    # Try loading from peak_decode file
    if peak_decode_path and Path(peak_decode_path).exists():
        with open(peak_decode_path) as f:
            pd = json.load(f)
        # Adapt to your data format
        if isinstance(pd, list):
            peak_aas = [p.get('aa', p.get('residue_name', '')) for p in pd]
        elif isinstance(pd, dict):
            peak_aas = pd.get('peak_amino_acids', pd.get('aas', []))
    
    # Or extract from segments
    if not peak_aas and segments:
        for seg in segments:
            model = seg.get('model_name', seg.get('best_model', seg.get('macro', '')))
            if 'gaussian' in str(model).lower() or model == 'localized':
                aa = seg.get('peak_aa', seg.get('peak_residue_name', None))
                if aa:
                    peak_aas.append(aa)
    
    if not peak_aas:
        print("\n    ⚠ No peak amino acid data found.")
        print("    To compute this, the pipeline needs to record which")
        print("    amino acid is at each Gaussian peak position.")
        print()
        print("    Using paper values for the test:")
        # Use reported values from paper
        n_peaks = 738
        n_pro_gly = round(738 * 0.108)  # 10.8%
        _run_chi2(n_peaks, n_pro_gly)
        return
    
    counts = Counter(peak_aas)
    n_peaks = len(peak_aas)
    n_pro = counts.get('PRO', 0) + counts.get('P', 0)
    n_gly = counts.get('GLY', 0) + counts.get('G', 0)
    n_pro_gly = n_pro + n_gly
    
    print(f"\n    Total peaks: {n_peaks}")
    print(f"    Pro: {n_pro} ({n_pro/n_peaks*100:.1f}%)")
    print(f"    Gly: {n_gly} ({n_gly/n_peaks*100:.1f}%)")
    print(f"    Pro+Gly: {n_pro_gly} ({n_pro_gly/n_peaks*100:.1f}%)")
    
    _run_chi2(n_peaks, n_pro_gly)


def _run_chi2(n_peaks, n_pro_gly):
    """Run χ² test against baseline."""
    from scipy import stats
    
    # Proteome baseline: Pro ~5.1%, Gly ~7.1% → Pro+Gly ~12.2%
    baseline_freq = 0.122
    observed_freq = n_pro_gly / n_peaks
    
    expected_pg = n_peaks * baseline_freq
    expected_other = n_peaks * (1 - baseline_freq)
    observed_other = n_peaks - n_pro_gly
    
    # χ² test (1 df)
    chi2 = ((n_pro_gly - expected_pg)**2 / expected_pg + 
            (observed_other - expected_other)**2 / expected_other)
    p_value = stats.chi2.sf(chi2, df=1)
    
    # Binomial CI for Pro+Gly fraction
    ci_lo, ci_hi = stats.binom.interval(0.95, n_peaks, observed_freq)
    ci_lo /= n_peaks
    ci_hi /= n_peaks
    
    print(f"\n    Baseline (proteome): {baseline_freq*100:.1f}%")
    print(f"    Observed at peaks:   {observed_freq*100:.1f}%")
    print(f"    95% CI:              [{ci_lo*100:.1f}%, {ci_hi*100:.1f}%]")
    print(f"    χ²(1) = {chi2:.2f}, p = {p_value:.4f}")
    
    if p_value > 0.05:
        print(f"\n    ★ NOT SIGNIFICANT (p = {p_value:.3f})")
        print(f"    Pro+Gly is NOT enriched at Gaussian peak positions.")
        print(f"    FOR PAPER: 'No significant enrichment of Pro+Gly at")
        print(f"    peak positions (χ²(1) = {chi2:.1f}, p = {p_value:.2f}; 95% CI")
        print(f"    [{ci_lo*100:.1f}%, {ci_hi*100:.1f}%] includes baseline {baseline_freq*100:.1f}%).'")
    else:
        direction = "enriched" if observed_freq > baseline_freq else "depleted"
        print(f"\n    ⚠ SIGNIFICANT: Pro+Gly is {direction} (p = {p_value:.4f})")


# ═══════════════════════════════════════════════════════════════════════
#  #9: HELIX PEAK PAIR COUNT FOR EXCLUSION ZONE
# ═══════════════════════════════════════════════════════════════════════

def stat_exclusion_pairs(segments):
    """Count helix-helix Gaussian peak pairs used in exclusion zone analysis."""
    print("\n  #9: HELIX GAUSSIAN PEAK PAIRS (EXCLUSION ZONE)")
    print("  " + "─" * 50)
    
    # Group segments by protein
    proteins = {}
    for seg in segments:
        prot = seg.get('protein', seg.get('protein_id', seg.get('uniprot', 'unknown')))
        if prot not in proteins:
            proteins[prot] = []
        proteins[prot].append(seg)
    
    total_pairs = 0
    total_peaks = 0
    proteins_with_peaks = 0
    
    for prot, segs in proteins.items():
        # Find helix Gaussian peaks
        helix_peaks = []
        for seg in segs:
            ss = seg.get('ss', seg.get('ss_type', ''))
            model = seg.get('model_name', seg.get('best_model', seg.get('macro', '')))
            if ss == 'H' and ('gaussian' in str(model).lower() or model == 'localized'):
                # Get peak position
                pos = seg.get('peak_position', seg.get('start', seg.get('seg_start', 0)))
                helix_peaks.append(pos)
        
        if len(helix_peaks) >= 2:
            proteins_with_peaks += 1
            total_peaks += len(helix_peaks)
            # Count consecutive pairs
            helix_peaks.sort()
            n_pairs = len(helix_peaks) - 1
            total_pairs += n_pairs
    
    print(f"\n    Proteins with ≥2 helix peaks: {proteins_with_peaks}")
    print(f"    Total helix Gaussian peaks: {total_peaks}")
    print(f"    Total consecutive pairs: {total_pairs}")
    print(f"\n    FOR PAPER: 'Exclusion zone analysis based on {total_pairs}")
    print(f"    consecutive helix-helix Gaussian peak pairs across")
    print(f"    {proteins_with_peaks} proteins.'")
    
    if total_pairs == 0:
        print("\n    ⚠ No pairs found. Check segment key names.")
        if segments:
            print(f"    Sample segment keys: {list(segments[0].keys())}")


# ═══════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    if len(sys.argv) < 2:
        print("\n  Usage: python reviewer_stats.py <analysis_json> [peak_decode_json]")
        print("  The analysis JSON should contain the segment-level data from")
        print("  the main 557-protein analysis.")
        print()
        print("  Running with paper values only (χ² test)...\n")
        stat_peak_chi2()
        return
    
    analysis_path = sys.argv[1]
    peak_path = sys.argv[2] if len(sys.argv) > 2 else None
    
    print(f"\n  Loading: {analysis_path}")
    data = load_analysis(analysis_path)
    
    # Extract segments - adapt to your JSON structure
    if isinstance(data, list):
        segments = data
    elif isinstance(data, dict):
        segments = data.get('segments', data.get('all_segments', []))
        if not segments:
            # Maybe segments are nested under proteins
            for key, val in data.items():
                if isinstance(val, dict) and 'segments' in val:
                    segments.extend(val['segments'])
                elif isinstance(val, list):
                    segments.extend(val)
    
    print(f"  Loaded {len(segments)} segments")
    
    if not segments:
        print("  ⚠ No segments found. Check JSON structure.")
        print(f"  Top-level keys: {list(data.keys()) if isinstance(data, dict) else 'list'}")
        return
    
    # Run all stats
    stat_segment_lengths(segments)
    stat_delta_aicc(segments)
    stat_peak_chi2(segments, peak_path)
    stat_exclusion_pairs(segments)
    
    print("\n  " + "═" * 50)
    print("  Done. Add these numbers to the paper.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
RIVET CONSERVATION ACROSS EVOLUTION
=====================================

Tests whether curvature signatures are conserved across species.

Six analyses:
  1. Curvature fingerprint by species (table)
  2. Kingdom-level signatures
  3. Divergence-from-human gradient
  4. Gini coefficient conservation (rivet pattern)
  5. Rivet density across the tree of life
  6. Helix-curvature coupling universality

Usage:
    python rivet_conservation.py

Requires: all_species.csv in same directory

Author: Kase Knochenhauer / True North Construction LLC
"""

import csv
import json
import sys
import numpy as np
from pathlib import Path
from collections import defaultdict
from scipy import stats as sp_stats

SCRIPT_DIR = Path(__file__).parent
CSV_PATH = SCRIPT_DIR / "all_species.csv"
RESULTS_DIR = SCRIPT_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

MAMMALS = ['human', 'chimp', 'gorilla', 'mouse', 'rat', 'dog', 'cat', 'cow', 'pig']
VERTEBRATES = MAMMALS + ['chicken', 'zebrafish']
INVERTEBRATES = ['fly', 'worm', 'mosquito', 'honeybee']
FUNGI = ['yeast', 'fission_yeast', 'candida']
PLANTS = ['arabidopsis', 'rice', 'maize', 'soybean']
PROTISTS = ['malaria', 'dictyostelium', 'trypanosoma', 'leishmania']
BACTERIA = ['ecoli', 'bacillus', 'salmonella', 'pseudomonas',
            'tuberculosis', 'staph', 'helicobacter', 'campylobacter']

KINGDOMS = {
    'Mammals': MAMMALS, 'Vertebrates': VERTEBRATES,
    'Invertebrates': INVERTEBRATES, 'Fungi': FUNGI,
    'Plants': PLANTS, 'Protists': PROTISTS, 'Bacteria': BACTERIA,
}

DISTANCE_TIERS = [
    ('Great apes', ['chimp', 'gorilla']),
    ('Rodents', ['mouse', 'rat']),
    ('Other mammals', ['dog', 'cat', 'cow', 'pig']),
    ('Other vertebrates', ['chicken', 'zebrafish']),
    ('Invertebrates', ['fly', 'worm', 'mosquito', 'honeybee']),
    ('Fungi', ['yeast', 'fission_yeast', 'candida']),
    ('Plants', ['arabidopsis', 'rice', 'maize', 'soybean']),
    ('Protists', ['malaria', 'dictyostelium', 'trypanosoma', 'leishmania']),
    ('Bacteria', ['ecoli', 'bacillus', 'salmonella', 'pseudomonas',
                  'tuberculosis', 'staph', 'helicobacter', 'campylobacter']),
]


def load_data():
    species_data = defaultdict(list)
    with open(CSV_PATH, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sp = row['species']
            try:
                protein = {
                    'uniprot_id': row['uniprot_id'],
                    'length': int(row['length']),
                    'n_segments': int(row['n_segments']),
                    'n_gauss_peaks': int(row['n_gauss_peaks']),
                    'mean_kappa': float(row['mean_kappa']),
                    'max_kappa': float(row['max_kappa']),
                    'std_kappa': float(row['std_kappa']),
                    'kappa_gini': float(row['kappa_gini']),
                    'kappa_sq_per_res': float(row['kappa_sq_per_res']),
                    'ss_frac_H': float(row['ss_frac_H']),
                    'ss_frac_E': float(row['ss_frac_E']),
                    'ss_frac_C': float(row['ss_frac_C']),
                    'mean_plddt': float(row['mean_plddt']),
                    'plddt_below_70': float(row['plddt_below_70']),
                    'kappa_sq_q25': float(row['kappa_sq_q25']),
                    'kappa_sq_q50': float(row['kappa_sq_q50']),
                    'kappa_sq_q75': float(row['kappa_sq_q75']),
                    'kappa_sq_q90': float(row['kappa_sq_q90']),
                }
                species_data[sp].append(protein)
            except (ValueError, KeyError):
                continue
    return species_data


def compute_species_stats(proteins):
    good = [p for p in proteins if p['mean_plddt'] > 70 and p['length'] > 50]
    if len(good) < 20:
        return None
    
    means = np.array([p['mean_kappa'] for p in good])
    ginis = np.array([p['kappa_gini'] for p in good])
    stds = np.array([p['std_kappa'] for p in good])
    h_fracs = np.array([p['ss_frac_H'] for p in good])
    e_fracs = np.array([p['ss_frac_E'] for p in good])
    peaks = np.array([p['n_gauss_peaks'] for p in good])
    lengths = np.array([p['length'] for p in good])
    rivet_density = peaks / lengths * 100
    
    return {
        'n_total': len(proteins),
        'n_good': len(good),
        'mean_kappa_mean': float(np.mean(means)),
        'mean_kappa_std': float(np.std(means)),
        'mean_gini': float(np.mean(ginis)),
        'gini_std': float(np.std(ginis)),
        'mean_std_kappa': float(np.mean(stds)),
        'mean_h_frac': float(np.mean(h_fracs)),
        'mean_e_frac': float(np.mean(e_fracs)),
        'median_kappa_sq': float(np.median([p['kappa_sq_q50'] for p in good])),
        'mean_rivet_density': float(np.mean(rivet_density)),
        'rivet_density_std': float(np.std(rivet_density)),
        'mean_length': float(np.mean(lengths)),
        '_means': means, '_ginis': ginis, '_stds': stds,
        '_rivet_density': rivet_density, '_h_fracs': h_fracs,
    }


def main():
    print(f"\n  RIVET CONSERVATION ACROSS EVOLUTION")
    print(f"  {'═'*55}")
    
    data = load_data()
    print(f"  Species loaded: {len(data)}")
    print(f"  Total proteins: {sum(len(v) for v in data.values()):,}")
    
    stats = {}
    for sp, proteins in data.items():
        s = compute_species_stats(proteins)
        if s is not None:
            stats[sp] = s
    
    print(f"  Species with sufficient data: {len(stats)}")
    
    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  1. CURVATURE FINGERPRINT BY SPECIES")
    print(f"  {'─'*55}")
    print(f"    {'Species':20s} {'N':>6s} {'<κ>':>7s} {'Gini':>7s} {'σ(κ)':>7s} {'Rivets':>8s} {'H%':>6s} {'E%':>6s}")
    
    for sp in ['human', 'chimp', 'gorilla', 'mouse', 'rat', 'dog', 'cat', 'cow', 'pig',
               'chicken', 'zebrafish', 'fly', 'worm', 'yeast', 'fission_yeast',
               'arabidopsis', 'rice', 'ecoli', 'bacillus', 'tuberculosis',
               'malaria', 'trypanosoma']:
        if sp not in stats:
            continue
        s = stats[sp]
        print(f"    {sp:20s} {s['n_good']:6d} {s['mean_kappa_mean']:7.2f} "
              f"{s['mean_gini']:7.3f} {s['mean_std_kappa']:7.2f} "
              f"{s['mean_rivet_density']:7.2f}% "
              f"{s['mean_h_frac']*100:5.1f} {s['mean_e_frac']*100:5.1f}")
    
    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  2. KINGDOM-LEVEL CURVATURE SIGNATURES")
    print(f"  {'─'*55}")
    
    kingdom_stats = {}
    for kname, species_list in KINGDOMS.items():
        all_m, all_g, all_r, all_h = [], [], [], []
        n = 0
        for sp in species_list:
            if sp in stats:
                all_m.extend(stats[sp]['_means'].tolist())
                all_g.extend(stats[sp]['_ginis'].tolist())
                all_r.extend(stats[sp]['_rivet_density'].tolist())
                all_h.extend(stats[sp]['_h_fracs'].tolist())
                n += stats[sp]['n_good']
        if n > 50:
            kingdom_stats[kname] = {
                'n': n, '_means': np.array(all_m), '_ginis': np.array(all_g),
                '_rivets': np.array(all_r),
            }
            print(f"    {kname:15s}: n={n:6d}, <κ>={np.mean(all_m):5.2f}, "
                  f"Gini={np.mean(all_g):.3f}, "
                  f"Rivets={np.mean(all_r):.2f}/100res, "
                  f"H%={np.mean(all_h)*100:.1f}")
    
    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  3. CURVATURE DIVERGENCE FROM HUMAN")
    print(f"  {'─'*55}")
    
    if 'human' not in stats:
        print("    ⚠ Human data not found")
    else:
        human = stats['human']
        h_means = human['_means']
        
        print(f"    Human baseline: <κ>={human['mean_kappa_mean']:.2f}, "
              f"Gini={human['mean_gini']:.3f}, "
              f"Rivets={human['mean_rivet_density']:.2f}/100res")
        print()
        print(f"    {'Tier':20s} {'N':>6s} {'Δ<κ>':>8s} {'ΔGini':>8s} {'ΔRivets':>8s} {'MW p(κ)':>10s}")
        
        for tier_name, tier_species in DISTANCE_TIERS:
            t_m, t_g, t_r = [], [], []
            n = 0
            for sp in tier_species:
                if sp in stats:
                    t_m.extend(stats[sp]['_means'].tolist())
                    t_g.extend(stats[sp]['_ginis'].tolist())
                    t_r.extend(stats[sp]['_rivet_density'].tolist())
                    n += stats[sp]['n_good']
            if n < 20:
                continue
            
            t_m, t_g, t_r = np.array(t_m), np.array(t_g), np.array(t_r)
            d_m = np.mean(t_m) - human['mean_kappa_mean']
            d_g = np.mean(t_g) - human['mean_gini']
            d_r = np.mean(t_r) - human['mean_rivet_density']
            
            try:
                _, p_mw = sp_stats.mannwhitneyu(h_means, t_m, alternative='two-sided')
            except:
                p_mw = 1.0
            sig = '★★★' if p_mw < 0.001 else '★★' if p_mw < 0.01 else '★' if p_mw < 0.05 else 'ns'
            print(f"    {tier_name:20s} {n:6d} {d_m:+8.2f} {d_g:+8.3f} "
                  f"{d_r:+8.2f} {p_mw:10.2e} {sig}")
    
    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  4. GINI COEFFICIENT CONSERVATION")
    print(f"  {'─'*55}")
    print(f"    (Gini = curvature inequality: high = peaked/riveted, low = uniform)")
    print()
    
    mammal_ginis = [(sp, stats[sp]['mean_gini']) for sp in MAMMALS if sp in stats]
    bact_ginis = [(sp, stats[sp]['mean_gini']) for sp in BACTERIA if sp in stats]
    
    if mammal_ginis:
        vals = [g for _, g in mammal_ginis]
        cv = np.std(vals) / np.mean(vals) * 100
        print(f"    Mammals: {' '.join(f'{s}={g:.3f}' for s, g in mammal_ginis)}")
        print(f"    Mammal mean: {np.mean(vals):.3f} ± {np.std(vals):.3f} (CV={cv:.1f}%)")
    if bact_ginis:
        vals = [g for _, g in bact_ginis]
        print(f"    Bacteria: {' '.join(f'{s}={g:.3f}' for s, g in bact_ginis[:5])}")
        print(f"    Bacteria mean: {np.mean(vals):.3f} ± {np.std(vals):.3f}")
    
    if mammal_ginis and bact_ginis:
        # Pool individual protein Ginis for proper test
        m_pool = np.concatenate([stats[sp]['_ginis'] for sp in MAMMALS if sp in stats])
        b_pool = np.concatenate([stats[sp]['_ginis'] for sp in BACTERIA if sp in stats])
        _, p_gb = sp_stats.mannwhitneyu(m_pool, b_pool, alternative='two-sided')
        print(f"\n    Mammal vs Bacteria (pooled): MW p = {p_gb:.2e}")
        print(f"    Mammal pool: {np.mean(m_pool):.3f} (n={len(m_pool):,})")
        print(f"    Bacteria pool: {np.mean(b_pool):.3f} (n={len(b_pool):,})")
    
    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  5. RIVET DENSITY ACROSS TREE OF LIFE")
    print(f"  {'─'*55}")
    
    for sp in ['human', 'chimp', 'mouse', 'zebrafish', 'fly', 'worm',
               'yeast', 'arabidopsis', 'ecoli', 'tuberculosis', 'malaria']:
        if sp not in stats:
            continue
        s = stats[sp]
        bar = '█' * int(s['mean_rivet_density'] * 5)
        print(f"    {sp:20s}: {s['mean_rivet_density']:5.2f} ± {s['rivet_density_std']:5.2f} {bar}")
    
    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  6. HELIX-CURVATURE COUPLING IS UNIVERSAL")
    print(f"  {'─'*55}")
    
    all_h_fracs, all_ginis_flat = [], []
    for sp, s in stats.items():
        all_h_fracs.extend(s['_h_fracs'].tolist())
        all_ginis_flat.extend(s['_ginis'].tolist())
    
    r, p = sp_stats.pearsonr(all_h_fracs, all_ginis_flat)
    print(f"    Helix fraction vs Gini (all {len(all_h_fracs):,} proteins):")
    print(f"    Pearson r = {r:.3f}, p = {p:.2e}")
    sig = '★★★' if p < 0.001 else 'ns'
    print(f"    {sig}")
    print()
    
    for kname, species_list in KINGDOMS.items():
        k_h, k_g = [], []
        for sp in species_list:
            if sp in stats:
                k_h.extend(stats[sp]['_h_fracs'].tolist())
                k_g.extend(stats[sp]['_ginis'].tolist())
        if len(k_h) > 50:
            r_k, p_k = sp_stats.pearsonr(k_h, k_g)
            sig = '★★★' if p_k < 0.001 else '★★' if p_k < 0.01 else '★' if p_k < 0.05 else 'ns'
            print(f"    {kname:15s}: r = {r_k:.3f}, p = {p_k:.2e} {sig} (n={len(k_h):,})")
    
    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  VERDICT")
    print(f"  {'═'*55}")
    
    results = {
        'n_species': len(stats),
        'n_proteins_total': sum(s['n_good'] for s in stats.values()),
        'species_stats': {sp: {k: v for k, v in s.items() if not k.startswith('_')}
                         for sp, s in stats.items()},
    }
    out_path = RESULTS_DIR / "rivet_conservation.json"
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved: {out_path}")


if __name__ == "__main__":
    main()

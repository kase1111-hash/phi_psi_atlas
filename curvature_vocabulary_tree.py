#!/usr/bin/env python3
"""
CURVATURE VOCABULARY GRADIENT ACROSS THE TREE OF LIFE
======================================================

The paper shows all-α proteins have richer curvature vocabulary
(entropy 1.87 bits) than all-β (1.19 bits). But is this just a
property of individual folds, or does it scale to PROTEOMES?

Questions:
  1. Do kingdoms have distinct curvature vocabularies?
  2. Does vocabulary complexity correlate with organismal complexity?
  3. Is there a curvature signature that separates domains of life?
  4. After controlling for helix content, is there residual
     kingdom-specific curvature architecture?

Metrics (from all_species.csv, no AICc classification needed):
  - Kappa distribution shape: Gini, Q25/Q50/Q75/Q90 ratios
  - Curvature dynamic range: max_kappa / mean_kappa
  - Curvature concentration: std_kappa / mean_kappa (CV)
  - Inter-quartile curvature ratio: Q90/Q25 (tail heaviness)
  - Per-residue curvature intensity: kappa_sq_per_res

Usage:
    python curvature_vocabulary_tree.py

Author: Kase Knochenhauer / True North Construction LLC
"""

import csv
import json
import numpy as np
from pathlib import Path
from collections import defaultdict
from scipy import stats as sp_stats

SCRIPT_DIR = Path(__file__).parent
CSV_PATH = SCRIPT_DIR / "all_species.csv"
RESULTS_DIR = SCRIPT_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Phylogenetic organization
TREE = {
    'Archaea-like bacteria': ['tuberculosis'],  # Actinobacteria, deep-branching
    'Proteobacteria': ['ecoli', 'salmonella', 'pseudomonas', 'helicobacter', 'campylobacter'],
    'Firmicutes': ['bacillus', 'staph'],
    'Fungi': ['yeast', 'fission_yeast', 'candida'],
    'Plants': ['arabidopsis', 'rice', 'maize', 'soybean'],
    'Protists': ['malaria', 'dictyostelium', 'trypanosoma', 'leishmania'],
    'Ecdysozoa': ['fly', 'worm', 'mosquito', 'honeybee'],
    'Fish': ['zebrafish'],
    'Birds': ['chicken'],
    'Non-primate mammals': ['mouse', 'rat', 'dog', 'cat', 'cow', 'pig'],
    'Primates': ['human', 'chimp', 'gorilla'],
}

# Ordered by approximate divergence from human (Mya)
DIVERGENCE_ORDER = [
    ('Primates', 7),
    ('Non-primate mammals', 90),
    ('Birds', 310),
    ('Fish', 450),
    ('Ecdysozoa', 600),
    ('Fungi', 1000),
    ('Plants', 1500),
    ('Protists', 1500),
    ('Proteobacteria', 3500),
    ('Firmicutes', 3500),
    ('Archaea-like bacteria', 3500),
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
                    'mean_kappa': float(row['mean_kappa']),
                    'max_kappa': float(row['max_kappa']),
                    'std_kappa': float(row['std_kappa']),
                    'kappa_gini': float(row['kappa_gini']),
                    'kappa_sq_per_res': float(row['kappa_sq_per_res']),
                    'ss_frac_H': float(row['ss_frac_H']),
                    'ss_frac_E': float(row['ss_frac_E']),
                    'ss_frac_C': float(row['ss_frac_C']),
                    'mean_plddt': float(row['mean_plddt']),
                    'kappa_sq_q25': float(row['kappa_sq_q25']),
                    'kappa_sq_q50': float(row['kappa_sq_q50']),
                    'kappa_sq_q75': float(row['kappa_sq_q75']),
                    'kappa_sq_q90': float(row['kappa_sq_q90']),
                }
                species_data[sp].append(protein)
            except (ValueError, KeyError):
                continue
    return species_data


def compute_vocabulary(proteins):
    """Compute curvature vocabulary metrics for a set of proteins."""
    good = [p for p in proteins if p['mean_plddt'] > 70 and p['length'] > 50]
    if len(good) < 20:
        return None

    # Per-protein vocabulary metrics
    dynamic_range = []    # max/mean — how extreme are peaks?
    cv_kappa = []         # std/mean — coefficient of variation
    tail_ratio = []       # Q90/Q25 — how heavy-tailed?
    iqr_ratio = []        # Q75/Q25 — inter-quartile spread
    intensity = []        # kappa_sq_per_res — total curvature load
    ginis = []
    h_fracs = []
    e_fracs = []
    c_fracs = []
    means = []

    for p in good:
        mk = p['mean_kappa']
        if mk > 0.5:  # avoid division by ~zero
            dynamic_range.append(p['max_kappa'] / mk)
            cv_kappa.append(p['std_kappa'] / mk)
        
        q25 = p['kappa_sq_q25']
        q50 = p['kappa_sq_q50']
        q75 = p['kappa_sq_q75']
        q90 = p['kappa_sq_q90']
        
        if q25 > 0.1:
            tail_ratio.append(q90 / q25)
            iqr_ratio.append(q75 / q25)
        
        intensity.append(p['kappa_sq_per_res'])
        ginis.append(p['kappa_gini'])
        h_fracs.append(p['ss_frac_H'])
        e_fracs.append(p['ss_frac_E'])
        c_fracs.append(p['ss_frac_C'])
        means.append(mk)

    return {
        'n': len(good),
        'mean_kappa': np.mean(means),
        'mean_dynamic_range': np.mean(dynamic_range) if dynamic_range else 0,
        'mean_cv_kappa': np.mean(cv_kappa) if cv_kappa else 0,
        'mean_tail_ratio': np.mean(tail_ratio) if tail_ratio else 0,
        'mean_iqr_ratio': np.mean(iqr_ratio) if iqr_ratio else 0,
        'mean_intensity': np.mean(intensity),
        'mean_gini': np.mean(ginis),
        'mean_h_frac': np.mean(h_fracs),
        'mean_e_frac': np.mean(e_fracs),
        'mean_c_frac': np.mean(c_fracs),
        # Raw arrays for tests
        '_dynamic_range': np.array(dynamic_range),
        '_cv_kappa': np.array(cv_kappa),
        '_tail_ratio': np.array(tail_ratio) if tail_ratio else np.array([0]),
        '_intensity': np.array(intensity),
        '_ginis': np.array(ginis),
        '_h_fracs': np.array(h_fracs),
        '_means': np.array(means),
    }


def main():
    print(f"\n  CURVATURE VOCABULARY ACROSS THE TREE OF LIFE")
    print(f"  {'═'*55}")

    data = load_data()
    print(f"  Species: {len(data)}, Proteins: {sum(len(v) for v in data.values()):,}")

    # Per-species vocabulary
    sp_vocab = {}
    for sp, proteins in data.items():
        v = compute_vocabulary(proteins)
        if v:
            sp_vocab[sp] = v

    # Per-clade vocabulary
    clade_vocab = {}
    for clade, species_list in TREE.items():
        pooled = []
        for sp in species_list:
            pooled.extend(data.get(sp, []))
        v = compute_vocabulary(pooled)
        if v:
            clade_vocab[clade] = v

    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  1. CURVATURE VOCABULARY BY SPECIES")
    print(f"  {'─'*55}")
    print(f"    {'Species':18s} {'N':>6s} {'<κ>':>6s} {'DR':>6s} {'CV':>6s} {'Tail':>7s} {'Gini':>6s} {'H%':>5s} {'E%':>5s}")

    display_order = ['human', 'chimp', 'gorilla', 'mouse', 'rat', 'dog', 'cow',
                     'chicken', 'zebrafish', 'fly', 'worm',
                     'yeast', 'fission_yeast',
                     'arabidopsis', 'rice',
                     'ecoli', 'bacillus', 'pseudomonas', 'tuberculosis',
                     'malaria', 'dictyostelium']

    for sp in display_order:
        if sp not in sp_vocab:
            continue
        v = sp_vocab[sp]
        print(f"    {sp:18s} {v['n']:6d} {v['mean_kappa']:6.2f} "
              f"{v['mean_dynamic_range']:6.1f} {v['mean_cv_kappa']:6.2f} "
              f"{v['mean_tail_ratio']:7.1f} {v['mean_gini']:6.3f} "
              f"{v['mean_h_frac']*100:4.0f} {v['mean_e_frac']*100:4.0f}")

    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  2. VOCABULARY BY CLADE (evolutionary order)")
    print(f"  {'─'*55}")
    print(f"    {'Clade':25s} {'N':>6s} {'<κ>':>6s} {'DR':>6s} {'CV':>6s} {'Tail':>7s} {'Gini':>6s} {'H%':>5s}")

    for clade_name, div_mya in DIVERGENCE_ORDER:
        if clade_name not in clade_vocab:
            continue
        v = clade_vocab[clade_name]
        print(f"    {clade_name:25s} {v['n']:6d} {v['mean_kappa']:6.2f} "
              f"{v['mean_dynamic_range']:6.1f} {v['mean_cv_kappa']:6.2f} "
              f"{v['mean_tail_ratio']:7.1f} {v['mean_gini']:6.3f} "
              f"{v['mean_h_frac']*100:4.0f}")

    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  3. CURVATURE COMPLEXITY GRADIENT")
    print(f"  {'─'*55}")
    print(f"    Do different domains of life have different")
    print(f"    curvature complexity?")
    print()

    # Pool into 3 domains
    domains = {
        'Eukaryota': ['human', 'chimp', 'gorilla', 'mouse', 'rat', 'dog', 'cat',
                       'cow', 'pig', 'chicken', 'zebrafish', 'fly', 'worm',
                       'mosquito', 'honeybee', 'yeast', 'fission_yeast', 'candida',
                       'arabidopsis', 'rice', 'maize', 'soybean',
                       'malaria', 'dictyostelium', 'trypanosoma', 'leishmania'],
        'Bacteria': ['ecoli', 'bacillus', 'salmonella', 'pseudomonas',
                     'tuberculosis', 'staph', 'helicobacter', 'campylobacter'],
    }

    domain_data = {}
    for dname, species_list in domains.items():
        pooled = []
        for sp in species_list:
            pooled.extend(data.get(sp, []))
        v = compute_vocabulary(pooled)
        if v:
            domain_data[dname] = v
            print(f"    {dname:12s}: n={v['n']:6d}, <κ>={v['mean_kappa']:.2f}, "
                  f"DR={v['mean_dynamic_range']:.1f}, CV={v['mean_cv_kappa']:.2f}, "
                  f"Gini={v['mean_gini']:.3f}, H%={v['mean_h_frac']*100:.1f}")

    # Statistical test
    if 'Eukaryota' in domain_data and 'Bacteria' in domain_data:
        for metric_name, key in [('Dynamic range', '_dynamic_range'),
                                  ('CV(κ)', '_cv_kappa'),
                                  ('Tail ratio', '_tail_ratio'),
                                  ('Curvature intensity', '_intensity'),
                                  ('Gini', '_ginis'),
                                  ('Mean κ', '_means')]:
            euk = domain_data['Eukaryota'][key]
            bac = domain_data['Bacteria'][key]
            if len(euk) > 10 and len(bac) > 10:
                _, p = sp_stats.mannwhitneyu(euk, bac, alternative='two-sided')
                d_cohen = (np.mean(euk) - np.mean(bac)) / np.sqrt((np.var(euk) + np.var(bac)) / 2)
                sig = '★★★' if p < 0.001 else '★★' if p < 0.01 else '★' if p < 0.05 else 'ns'
                print(f"    {metric_name:20s}: Euk={np.mean(euk):.3f} vs Bac={np.mean(bac):.3f}, "
                      f"d={d_cohen:+.3f}, p={p:.2e} {sig}")

    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  4. RESIDUAL VOCABULARY AFTER HELIX CORRECTION")
    print(f"  {'─'*55}")
    print(f"    Does kingdom identity predict curvature BEYOND")
    print(f"    what helix fraction alone explains?")
    print()

    # For each species, compute residual: actual Gini minus Gini
    # predicted by helix fraction (using global regression)
    all_h, all_g, all_sp = [], [], []
    for sp, v in sp_vocab.items():
        for i in range(len(v['_ginis'])):
            all_h.append(v['_h_fracs'][i])
            all_g.append(v['_ginis'][i])
            all_sp.append(sp)

    all_h = np.array(all_h)
    all_g = np.array(all_g)

    # Linear regression: Gini = a * H_frac + b
    slope, intercept, r, p, se = sp_stats.linregress(all_h, all_g)
    predicted = slope * all_h + intercept
    residuals = all_g - predicted

    print(f"    Global regression: Gini = {slope:.4f} × H_frac + {intercept:.4f}")
    print(f"    R² = {r**2:.3f}, p = {p:.2e}")
    print()

    # Compute mean residual per species
    sp_residuals = defaultdict(list)
    for i, sp in enumerate(all_sp):
        sp_residuals[sp].append(residuals[i])

    print(f"    {'Species':18s} {'Mean residual':>14s} {'N':>6s}  Interpretation")
    print(f"    {'-'*60}")
    
    sorted_sp = sorted(sp_residuals.keys(), 
                       key=lambda x: np.mean(sp_residuals[x]))
    
    for sp in sorted_sp:
        res = sp_residuals[sp]
        mr = np.mean(res)
        n = len(res)
        if n < 50:
            continue
        interp = ''
        if mr > 0.003:
            interp = 'MORE peaked than helix predicts'
        elif mr < -0.003:
            interp = 'LESS peaked than helix predicts'
        else:
            interp = '≈ as expected from helix content'
        print(f"    {sp:18s} {mr:+14.4f} {n:6d}  {interp}")

    # Kingdom-level residuals
    print()
    kingdom_groups = {
        'Mammals': ['human', 'chimp', 'gorilla', 'mouse', 'rat', 'dog', 'cat', 'cow', 'pig'],
        'Invertebrates': ['fly', 'worm', 'mosquito', 'honeybee'],
        'Fungi': ['yeast', 'fission_yeast', 'candida'],
        'Plants': ['arabidopsis', 'rice', 'maize', 'soybean'],
        'Bacteria': ['ecoli', 'bacillus', 'salmonella', 'pseudomonas',
                     'tuberculosis', 'staph', 'helicobacter', 'campylobacter'],
    }

    print(f"    {'Kingdom':15s} {'Mean residual':>14s} {'N':>8s} {'p vs 0':>10s}")
    for kname, sp_list in kingdom_groups.items():
        k_res = []
        for sp in sp_list:
            if sp in sp_residuals:
                k_res.extend(sp_residuals[sp])
        if len(k_res) > 50:
            mr = np.mean(k_res)
            _, p_t = sp_stats.ttest_1samp(k_res, 0)
            sig = '★★★' if p_t < 0.001 else '★★' if p_t < 0.01 else '★' if p_t < 0.05 else 'ns'
            print(f"    {kname:15s} {mr:+14.5f} {len(k_res):8d} {p_t:10.2e} {sig}")

    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  5. CURVATURE DISTRIBUTION SHAPE (Q25/Q50/Q75/Q90)")
    print(f"  {'─'*55}")
    print(f"    Normalized quartile profiles reveal distribution shape")
    print(f"    independent of overall magnitude.")
    print()

    # For each species, compute normalized quartile profile
    # Normalize by Q50 (median) to compare shapes
    print(f"    {'Species':18s} {'Q25/Q50':>8s} {'Q75/Q50':>8s} {'Q90/Q50':>8s} {'Shape':>15s}")
    
    for sp in display_order:
        if sp not in sp_vocab:
            continue
        good = [p for p in data[sp] if p['mean_plddt'] > 70 and p['length'] > 50]
        if len(good) < 50:
            continue
        
        q25s = [p['kappa_sq_q25'] for p in good]
        q50s = [p['kappa_sq_q50'] for p in good]
        q75s = [p['kappa_sq_q75'] for p in good]
        q90s = [p['kappa_sq_q90'] for p in good]
        
        # Use median of ratios (only where q50 > 0.5 to avoid div/0)
        r25, r75, r90 = [], [], []
        for i in range(len(good)):
            if q50s[i] > 0.5:
                r25.append(q25s[i] / q50s[i])
                r75.append(q75s[i] / q50s[i])
                r90.append(q90s[i] / q50s[i])
        
        if len(r25) < 50:
            continue
        
        mr25 = np.median(r25)
        mr75 = np.median(r75)
        mr90 = np.median(r90)
        
        # Classify shape
        if mr90 > 50:
            shape = 'very heavy tail'
        elif mr90 > 30:
            shape = 'heavy tail'
        elif mr90 > 15:
            shape = 'moderate tail'
        else:
            shape = 'light tail'
        
        print(f"    {sp:18s} {mr25:8.2f} {mr75:8.2f} {mr90:8.2f} {shape:>15s}")

    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  6. HELIX-CURVATURE COUPLING STRENGTH BY CLADE")
    print(f"  {'─'*55}")
    print(f"    Correlation between helix fraction and curvature")
    print(f"    within each clade. If coupling is universal,")
    print(f"    slopes should be similar.")
    print()

    print(f"    {'Clade':25s} {'r':>7s} {'slope':>8s} {'p':>10s} {'N':>6s}")
    for clade_name, _ in DIVERGENCE_ORDER:
        if clade_name not in clade_vocab:
            continue
        v = clade_vocab[clade_name]
        h = v['_h_fracs']
        g = v['_ginis']
        if len(h) < 30:
            continue
        sl, ic, r_c, p_c, _ = sp_stats.linregress(h, g)
        sig = '★★★' if p_c < 0.001 else '★★' if p_c < 0.01 else '★' if p_c < 0.05 else 'ns'
        print(f"    {clade_name:25s} {r_c:7.3f} {sl:8.4f} {p_c:10.2e} {v['n']:6d} {sig}")

    # ═══════════════════════════════════════════════════════════
    print(f"\n  {'═'*55}")
    print(f"  VERDICT")
    print(f"  {'═'*55}")

    # Save results
    results = {
        'n_species': len(sp_vocab),
        'n_proteins': sum(v['n'] for v in sp_vocab.values()),
        'global_regression': {
            'slope': float(slope),
            'intercept': float(intercept),
            'r_squared': float(r**2),
        },
        'species_vocabulary': {
            sp: {k: v for k, v in vocab.items() if not k.startswith('_')}
            for sp, vocab in sp_vocab.items()
        },
        'clade_vocabulary': {
            cl: {k: v for k, v in vocab.items() if not k.startswith('_')}
            for cl, vocab in clade_vocab.items()
        },
    }

    out_path = RESULTS_DIR / "curvature_vocabulary_tree.json"
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved: {out_path}")


if __name__ == "__main__":
    main()

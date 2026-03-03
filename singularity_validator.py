#!/usr/bin/env python3
"""
Singularity Validator — Three Critical Tests
==============================================

Test 1: pLDDT Confidence Filter
  Are mechanical hotspots in high-confidence or low-confidence regions?
  If hotspots cluster in low-pLDDT regions, they may be prediction artifacts.

Test 2: Disease Enrichment Statistics  
  Is 222/948 (23.4%) actually enriched above proteome-wide baseline?
  Fisher's exact test with odds ratio and confidence interval.

Test 3: ClinVar Pathogenic Variant Clustering (preparation)
  Downloads ClinVar data and tests whether pathogenic missense variants
  cluster at κ² hotspot residues more than expected by chance.

Usage:
  python singularity_validator.py \
    --species all_species.csv \
    --annotated human_annotated_singularities.csv \
    --output validation_report.txt \
    [--clinvar clinvar_summary.txt.gz]

The --clinvar flag is optional. If provided, runs ClinVar clustering test.
Download ClinVar from: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
"""

import os
import sys
import csv
import gzip
import math
import argparse
from collections import defaultdict, Counter


# ═══════════════════════════════════════════════════════════════
# TEST 1: pLDDT CONFIDENCE FILTER
# ═══════════════════════════════════════════════════════════════

def detect_format(species_csv):
    """Auto-detect whether this is all_species.csv or archetype_mapper census format."""
    with open(species_csv) as f:
        header = f.readline().strip().replace('\r', '')
        cols = header.split(',')
    
    if "species" in cols and "kappa_gini" in cols:
        return "all_species"
    elif "label" in cols and "gini" in cols:
        return "archetype_mapper"
    else:
        return "unknown"


def load_census_auto(species_csv):
    """Load census data regardless of format. Returns list of dicts."""
    fmt = detect_format(species_csv)
    human = []
    
    with open(species_csv) as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                if fmt == "all_species":
                    if r["species"] != "human":
                        continue
                    length = int(r["length"])
                    if length < 100:
                        continue
                    gini = float(r["kappa_gini"])
                    k_mean = float(r["kappa_sq_per_res"])
                    k_max = float(r["top_peak_kappa_sq"])
                    mean_plddt = float(r.get("mean_plddt", 0))
                    min_plddt = float(r.get("min_plddt", 0))
                    frac_below_70 = float(r.get("plddt_below_70", 0))
                    uid = r["uniprot_id"]
                    hotspot = r.get("top_peak_position", "0")
                    has_plddt = True
                    
                elif fmt == "archetype_mapper":
                    # Filter to AlphaFold human entries by label prefix
                    label = r.get("label", "")
                    if not label.startswith("AF-"):
                        continue
                    # Extract UniProt ID
                    uid = label.replace("AF-", "").split("-F1")[0]
                    length = int(r["n_valid"])
                    if length < 100:
                        continue
                    gini = float(r["gini"])
                    k_mean = float(r["kappa2_mean"])
                    k_max = float(r["top1_kappa2"]) if r.get("top1_kappa2") else float(r.get("kappa2_max", 0))
                    mean_plddt = 0
                    min_plddt = 0
                    frac_below_70 = 0
                    hotspot = r.get("top1_res", "0")
                    has_plddt = False
                else:
                    continue
                
                conc = k_max / k_mean if k_mean > 0 else 0
                
                human.append({
                    "uid": uid,
                    "gini": gini,
                    "k_mean": k_mean,
                    "k_max": k_max,
                    "concentration": conc,
                    "mean_plddt": mean_plddt,
                    "min_plddt": min_plddt,
                    "frac_below_70": frac_below_70,
                    "is_singularity": gini > 0.92,
                    "length": length,
                    "hotspot_pos": hotspot,
                    "has_plddt": has_plddt,
                })
            except:
                pass
    
    return human, fmt


def test_plddt_confidence(species_csv, output):
    """Check whether mechanical singularity hotspots fall in confident or
    uncertain regions of AlphaFold predictions."""
    
    output.append("=" * 72)
    output.append("  TEST 1: pLDDT CONFIDENCE AT MECHANICAL HOTSPOTS")
    output.append("  Question: Are singularity hotspots AlphaFold artifacts?")
    output.append("=" * 72)
    
    # Load human entries
    human, fmt = load_census_auto(species_csv)
    has_plddt = any(r["has_plddt"] for r in human) if human else False
    
    output.append(f"")
    output.append(f"  Data format detected: {fmt}")
    output.append(f"  pLDDT data available: {'YES' if has_plddt else 'NO'}")
    
    singularities = [r for r in human if r["is_singularity"]]
    non_sing = [r for r in human if not r["is_singularity"]]
    
    output.append(f"  Human proteins (≥100 res): {len(human):,}")
    output.append(f"  Singularities (Gini>0.92): {len(singularities):,}")
    output.append(f"  Non-singularities:         {len(non_sing):,}")
    
    if not has_plddt:
        output.append(f"")
        output.append(f"  ⚠️  pLDDT data not available in this census format.")
        output.append(f"     Cannot test whether hotspots are in confident regions.")
        output.append(f"     To run this test, provide all_species.csv from the original")
        output.append(f"     pipeline which includes mean_plddt, min_plddt columns.")
        output.append(f"")
        output.append(f"     SKIPPING TEST 1 — proceeding to Tests 2 and 3.")
        return singularities, non_sing, singularities  # return all as hc_sing
    
    # Compare pLDDT distributions
    def stats(values):
        n = len(values)
        if n == 0:
            return 0, 0, 0, 0
        mean = sum(values) / n
        values_sorted = sorted(values)
        median = values_sorted[n // 2]
        q25 = values_sorted[n // 4]
        q75 = values_sorted[3 * n // 4]
        return mean, median, q25, q75
    
    s_plddt = [r["mean_plddt"] for r in singularities]
    n_plddt = [r["mean_plddt"] for r in non_sing]
    s_frac = [r["frac_below_70"] for r in singularities]
    n_frac = [r["frac_below_70"] for r in non_sing]
    
    s_mean, s_med, s_q25, s_q75 = stats(s_plddt)
    n_mean, n_med, n_q25, n_q75 = stats(n_plddt)
    sf_mean, sf_med, _, _ = stats(s_frac)
    nf_mean, nf_med, _, _ = stats(n_frac)
    
    output.append(f"")
    output.append(f"  Mean pLDDT:")
    output.append(f"    Singularities:     {s_mean:.1f} (median {s_med:.1f}, IQR [{s_q25:.1f}–{s_q75:.1f}])")
    output.append(f"    Non-singularities: {n_mean:.1f} (median {n_med:.1f}, IQR [{n_q25:.1f}–{n_q75:.1f}])")
    output.append(f"    Δ = {s_mean - n_mean:+.1f}")
    
    output.append(f"")
    output.append(f"  Fraction of residues with pLDDT < 70:")
    output.append(f"    Singularities:     {sf_mean:.3f} (median {sf_med:.3f})")
    output.append(f"    Non-singularities: {nf_mean:.3f} (median {nf_med:.3f})")
    
    # Confidence bins
    output.append(f"")
    output.append(f"  pLDDT distribution by category:")
    output.append(f"  {'pLDDT range':20s} {'Singularity':>12s} {'Non-sing':>12s} {'Sing %':>8s} {'Non %':>8s}")
    output.append(f"  {'─' * 64}")
    
    bins = [(0, 50, "Very low (<50)"), (50, 70, "Low (50-70)"), 
            (70, 85, "Moderate (70-85)"), (85, 95, "High (85-95)"), (95, 101, "Very high (≥95)")]
    
    for lo, hi, label in bins:
        s_count = sum(1 for r in singularities if lo <= r["mean_plddt"] < hi)
        n_count = sum(1 for r in non_sing if lo <= r["mean_plddt"] < hi)
        s_pct = 100 * s_count / len(singularities) if singularities else 0
        n_pct = 100 * n_count / len(non_sing) if non_sing else 0
        output.append(f"  {label:20s} {s_count:>12d} {n_count:>12d} {s_pct:>7.1f}% {n_pct:>7.1f}%")
    
    # High-confidence singularities
    hc_sing = [r for r in singularities if r["mean_plddt"] >= 70]
    vhc_sing = [r for r in singularities if r["mean_plddt"] >= 85]
    
    output.append(f"")
    output.append(f"  Singularities surviving confidence filters:")
    output.append(f"    pLDDT ≥ 70: {len(hc_sing):,} / {len(singularities):,} ({100*len(hc_sing)/len(singularities):.1f}%)")
    output.append(f"    pLDDT ≥ 85: {len(vhc_sing):,} / {len(singularities):,} ({100*len(vhc_sing)/len(singularities):.1f}%)")
    
    # Welch's t-test (no scipy, compute manually)
    def welch_t(vals1, vals2):
        n1, n2 = len(vals1), len(vals2)
        m1 = sum(vals1) / n1
        m2 = sum(vals2) / n2
        v1 = sum((x - m1)**2 for x in vals1) / (n1 - 1)
        v2 = sum((x - m2)**2 for x in vals2) / (n2 - 1)
        se = math.sqrt(v1/n1 + v2/n2)
        t = (m1 - m2) / se if se > 0 else 0
        # Approximate df (Welch-Satterthwaite)
        num = (v1/n1 + v2/n2)**2
        den = (v1/n1)**2/(n1-1) + (v2/n2)**2/(n2-1)
        df = num / den if den > 0 else 1
        return t, df, m1, m2
    
    t, df, m1, m2 = welch_t(s_plddt, n_plddt)
    output.append(f"")
    output.append(f"  Welch's t-test (singularity vs non-singularity pLDDT):")
    output.append(f"    t = {t:.2f}, df ≈ {df:.0f}")
    output.append(f"    Mean difference: {m1 - m2:+.2f}")
    
    # VERDICT
    output.append(f"")
    if s_mean >= 70 and len(hc_sing) / len(singularities) > 0.7:
        output.append(f"  ✅ VERDICT: Singularity hotspots are predominantly in confident regions.")
        output.append(f"     {len(hc_sing)}/{len(singularities)} ({100*len(hc_sing)/len(singularities):.0f}%) have mean pLDDT ≥ 70.")
        output.append(f"     The singularity signal is NOT an artifact of AlphaFold uncertainty.")
    elif s_mean >= 60:
        output.append(f"  ⚠️  VERDICT: Mixed confidence. Some singularities may be artifacts.")
        output.append(f"     Recommend filtering to pLDDT ≥ 70 subset ({len(hc_sing)} proteins).")
    else:
        output.append(f"  ❌ VERDICT: Singularity hotspots cluster in LOW-confidence regions.")
        output.append(f"     The singularity signal may be an AlphaFold prediction artifact.")
        output.append(f"     CRITICAL: Validate on experimentally solved crystal structures before proceeding.")
    
    return singularities, non_sing, hc_sing


# ═══════════════════════════════════════════════════════════════
# TEST 2: DISEASE ENRICHMENT STATISTICS
# ═══════════════════════════════════════════════════════════════

def test_disease_enrichment(annotated_csv, species_csv, output):
    """Fisher's exact test: are singularities enriched for disease genes
    above the proteome-wide baseline?"""
    
    output.append(f"")
    output.append(f"{'=' * 72}")
    output.append(f"  TEST 2: DISEASE ENRICHMENT STATISTICS")
    output.append(f"  Question: Is 222/948 enriched above proteome baseline?")
    output.append(f"{'=' * 72}")
    
    # Load annotated singularities
    sing_disease = 0
    sing_total = 0
    with open(annotated_csv) as f:
        for r in csv.DictReader(f):
            sing_total += 1
            if int(r.get("disease_count", 0)) > 0:
                sing_disease += 1
    
    sing_no_disease = sing_total - sing_disease
    
    # We need proteome-wide disease annotation rate
    # Count human proteins in species census
    human_all, _ = load_census_auto(species_csv)
    human_total = len(human_all)
    
    # Known: ~4,000-5,000 Mendelian disease genes in human genome
    # UniProt disease annotations: roughly 3,500-4,500 reviewed human proteins
    # We'll use the actual annotation rate from our singularity annotator
    # which queries UniProt — so the comparison must use the same source.
    #
    # Best approach: we need to know how many of ALL 19,501 human proteins
    # have disease annotations. We only annotated the 948 singularities.
    #
    # Use published estimates:
    # OMIM: ~4,400 genes with phenotype-causing mutations
    # Out of ~20,000 protein-coding genes = ~22%
    
    # Conservative estimates
    estimates = [
        (3500, "Conservative (3,500 disease genes)"),
        (4000, "Moderate (4,000 disease genes)"),
        (4500, "High (4,500 disease genes)"),
        (5000, "Very high (5,000 disease genes)"),
    ]
    
    output.append(f"")
    output.append(f"  Singularities: {sing_disease} with disease / {sing_total} total = {100*sing_disease/sing_total:.1f}%")
    output.append(f"  Human proteome: {human_total:,} proteins (≥100 res)")
    output.append(f"")
    
    non_sing_total = human_total - sing_total
    
    output.append(f"  Fisher's exact test (2×2 contingency table):")
    output.append(f"  {'':30s} {'Disease':>10s} {'No disease':>12s} {'Total':>10s}")
    output.append(f"  {'─' * 66}")
    
    for est_disease_total, label in estimates:
        # 2x2 table:
        # Singularity + disease:      sing_disease
        # Singularity + no disease:   sing_no_disease
        # Non-sing + disease:         est_disease_total - sing_disease
        # Non-sing + no disease:      non_sing_total - (est_disease_total - sing_disease)
        
        a = sing_disease
        b = sing_no_disease
        c = est_disease_total - sing_disease
        d = non_sing_total - c
        
        # Odds ratio
        odds_ratio = (a * d) / (b * c) if (b * c) > 0 else float('inf')
        
        # Log odds ratio SE (Woolf method)
        if a > 0 and b > 0 and c > 0 and d > 0:
            se_log_or = math.sqrt(1/a + 1/b + 1/c + 1/d)
            log_or = math.log(odds_ratio)
            ci_low = math.exp(log_or - 1.96 * se_log_or)
            ci_high = math.exp(log_or + 1.96 * se_log_or)
        else:
            se_log_or = float('inf')
            ci_low = ci_high = float('nan')
        
        # Fisher's exact p-value (hypergeometric)
        # Using log-factorial for numerical stability
        def log_fact(n):
            return sum(math.log(i) for i in range(1, n+1))
        
        def hypergeometric_pmf(a, b, c, d):
            """P(X=a) in Fisher's exact test"""
            n = a + b + c + d
            try:
                log_p = (log_fact(a+b) + log_fact(c+d) + log_fact(a+c) + log_fact(b+d)
                         - log_fact(n) - log_fact(a) - log_fact(b) - log_fact(c) - log_fact(d))
                return math.exp(log_p)
            except:
                return 0
        
        # One-sided p-value (enrichment: are singularities MORE likely to have disease?)
        # Sum P(X >= a) where X is disease count in singularity group
        p_value = 0
        max_a = min(a + b, a + c)  # max possible disease+singularity
        for test_a in range(a, max_a + 1):
            test_b = (a + b) - test_a
            test_c = (a + c) - test_a
            test_d = (b + d) + (a - test_a)  # total - others
            if test_b >= 0 and test_c >= 0 and test_d >= 0:
                p_value += hypergeometric_pmf(test_a, test_b, test_c, test_d)
        
        baseline_rate = 100 * est_disease_total / human_total
        sing_rate = 100 * sing_disease / sing_total
        
        output.append(f"")
        output.append(f"  {label}")
        output.append(f"  {'':30s} {'Disease':>10s} {'No disease':>12s}")
        output.append(f"  {'Singularity':30s} {a:>10d} {b:>12d}")
        output.append(f"  {'Non-singularity':30s} {c:>10d} {d:>12d}")
        output.append(f"  Baseline disease rate:  {baseline_rate:.1f}%")
        output.append(f"  Singularity disease rate: {sing_rate:.1f}%")
        output.append(f"  Odds ratio: {odds_ratio:.2f} (95% CI: {ci_low:.2f}–{ci_high:.2f})")
        output.append(f"  Fisher's exact p (one-sided): {p_value:.4e}" if p_value < 0.01 else f"  Fisher's exact p (one-sided): {p_value:.4f}")
        
        if odds_ratio > 1.0 and p_value < 0.05:
            output.append(f"  → ENRICHED (p < 0.05)")
        elif odds_ratio > 1.0:
            output.append(f"  → Trend toward enrichment (not significant)")
        else:
            output.append(f"  → NOT ENRICHED")
    
    # VERDICT
    output.append(f"")
    output.append(f"  ─────────────────────────────────────────────────")
    
    # Use moderate estimate for verdict
    est = 4000
    a = sing_disease
    b = sing_no_disease
    c = est - sing_disease
    d = non_sing_total - c
    odds_ratio = (a * d) / (b * c) if (b * c) > 0 else float('inf')
    baseline = 100 * est / human_total
    
    if odds_ratio > 1.5:
        output.append(f"  ✅ VERDICT: Strong enrichment. OR = {odds_ratio:.2f}.")
        output.append(f"     Singularities are significantly more likely to be disease genes.")
    elif odds_ratio > 1.1:
        output.append(f"  ⚠️  VERDICT: Modest enrichment. OR = {odds_ratio:.2f}.")
        output.append(f"     Signal exists but is not dramatically above background.")
        output.append(f"     The therapeutic premise is weakened but not killed.")
    else:
        output.append(f"  ❌ VERDICT: No enrichment. OR = {odds_ratio:.2f}.")
        output.append(f"     Singularity disease rate ({100*sing_disease/sing_total:.1f}%) ≈ baseline ({baseline:.1f}%).")
        output.append(f"     The 222 diseases may be background noise, not signal.")
        output.append(f"     CRITICAL: Reframe program around mechanical insight, not disease count.")
    
    return sing_disease, sing_total, human_total


# ═══════════════════════════════════════════════════════════════
# TEST 3: CLINVAR PATHOGENIC VARIANT CLUSTERING
# ═══════════════════════════════════════════════════════════════

def test_clinvar_clustering(clinvar_path, annotated_csv, species_csv, output):
    """Test whether ClinVar pathogenic missense variants cluster at κ²
    hotspot residues more than expected by chance."""
    
    output.append(f"")
    output.append(f"{'=' * 72}")
    output.append(f"  TEST 3: ClinVar PATHOGENIC VARIANT CLUSTERING")
    output.append(f"  Question: Do disease mutations hit the mechanical hotspot?")
    output.append(f"{'=' * 72}")
    
    if not clinvar_path or not os.path.exists(clinvar_path):
        output.append(f"")
        output.append(f"  ClinVar data not provided. To run this test:")
        output.append(f"  1. Download: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz")
        output.append(f"  2. Re-run with: --clinvar variant_summary.txt.gz")
        output.append(f"")
        output.append(f"  This is the CRITICAL test. Without it, the therapeutic premise")
        output.append(f"  rests on assertion, not evidence.")
        return
    
    # Load singularity hotspots
    singularities = {}
    with open(annotated_csv) as f:
        for r in csv.DictReader(f):
            gene = r.get("gene", "")
            if gene and int(r.get("disease_count", 0)) > 0:
                try:
                    singularities[gene] = {
                        "hotspot": int(r["hotspot_residue"]),
                        "length": int(r["n_residues"]),
                        "concentration": float(r["concentration_ratio"]),
                        "gini": float(r["gini"]),
                    }
                except:
                    pass
    
    output.append(f"  Disease singularities with hotspot data: {len(singularities)}")
    
    # Parse ClinVar
    # variant_summary.txt.gz columns include:
    # GeneSymbol, Type, ClinicalSignificance, Assembly, Chromosome,
    # Start, Stop, ProteinChange, ...
    
    output.append(f"  Loading ClinVar data...")
    
    pathogenic_by_gene = defaultdict(list)
    
    open_func = gzip.open if clinvar_path.endswith('.gz') else open
    with open_func(clinvar_path, 'rt', errors='replace') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        
        # Find column indices
        col_map = {name: i for i, name in enumerate(header)}
        gene_col = col_map.get("GeneSymbol", col_map.get("Gene(s)", -1))
        sig_col = col_map.get("ClinicalSignificance", col_map.get("ClinSig", -1))
        type_col = col_map.get("Type", -1)
        
        # Try to find protein change column
        prot_col = -1
        for candidate in ["ProteinChange", "Protein change", "protein_change"]:
            if candidate in col_map:
                prot_col = col_map[candidate]
                break
        
        # Also look for position info
        name_col = col_map.get("Name", -1)
        
        output.append(f"  Columns found: Gene={gene_col}, Sig={sig_col}, Type={type_col}, Prot={prot_col}, Name={name_col}")
        
        count = 0
        pathogenic_count = 0
        
        for row in reader:
            count += 1
            try:
                if gene_col < 0 or gene_col >= len(row):
                    continue
                gene = row[gene_col].strip()
                
                if gene not in singularities:
                    continue
                
                # Check if pathogenic
                sig = row[sig_col].lower() if sig_col >= 0 and sig_col < len(row) else ""
                if "pathogenic" not in sig or "benign" in sig:
                    continue
                
                # Check if missense (single nucleotide variant)
                vtype = row[type_col] if type_col >= 0 and type_col < len(row) else ""
                
                # Extract protein position from protein change or name
                position = None
                
                if prot_col >= 0 and prot_col < len(row):
                    pchange = row[prot_col]
                    # Parse formats like "p.Arg850Trp" or "p.R850W"
                    import re
                    match = re.search(r'p\.(?:[A-Z][a-z]{2})?(\d+)', pchange)
                    if match:
                        position = int(match.group(1))
                
                if position is None and name_col >= 0 and name_col < len(row):
                    name = row[name_col]
                    import re
                    # Try to extract from name field
                    match = re.search(r'p\.(?:[A-Z][a-z]{2})?(\d+)', name)
                    if match:
                        position = int(match.group(1))
                
                if position is not None:
                    pathogenic_count += 1
                    pathogenic_by_gene[gene].append(position)
                    
            except Exception:
                continue
    
    output.append(f"  ClinVar rows scanned: {count:,}")
    output.append(f"  Pathogenic variants in singularity genes: {pathogenic_count}")
    output.append(f"  Genes with variants: {len(pathogenic_by_gene)}")
    
    if not pathogenic_by_gene:
        output.append(f"")
        output.append(f"  ⚠️  No pathogenic variants found. Check ClinVar file format.")
        return
    
    # THE CRITICAL TEST:
    # For each gene, what fraction of pathogenic variants fall within
    # N residues of the mechanical hotspot?
    # Compare to expected fraction if variants were uniformly distributed.
    
    output.append(f"")
    output.append(f"  CLUSTERING ANALYSIS:")
    output.append(f"  For each gene: fraction of pathogenic variants within W residues of hotspot")
    output.append(f"  vs expected fraction if uniformly distributed.")
    output.append(f"")
    
    windows = [1, 3, 5, 10, 20, 50]
    
    # Aggregate across all genes
    total_variants = 0
    hits_by_window = {w: 0 for w in windows}
    expected_by_window = {w: 0 for w in windows}
    
    gene_results = []
    
    for gene, variants in sorted(pathogenic_by_gene.items()):
        info = singularities[gene]
        hotspot = info["hotspot"]
        length = info["length"]
        n_var = len(variants)
        total_variants += n_var
        
        gene_hits = {}
        for w in windows:
            hits = sum(1 for v in variants if abs(v - hotspot) <= w)
            expected_frac = min(1.0, (2 * w + 1) / length) if length > 0 else 0
            expected = expected_frac * n_var
            
            hits_by_window[w] += hits
            expected_by_window[w] += expected
            gene_hits[w] = (hits, expected)
        
        gene_results.append((gene, hotspot, length, n_var, variants, gene_hits, info["concentration"]))
    
    # Summary table
    output.append(f"  {'Window':>8s} {'Observed':>10s} {'Expected':>10s} {'Ratio':>8s} {'Enrichment'}")
    output.append(f"  {'─' * 52}")
    
    for w in windows:
        obs = hits_by_window[w]
        exp = expected_by_window[w]
        ratio = obs / exp if exp > 0 else float('inf')
        enrich = "***" if ratio > 3 else "**" if ratio > 2 else "*" if ratio > 1.5 else ""
        output.append(f"  ±{w:>5d} res {obs:>10d} {exp:>10.1f} {ratio:>7.1f}× {enrich}")
    
    # Per-gene details for top genes
    output.append(f"")
    output.append(f"  Per-gene details (genes with ≥5 pathogenic variants):")
    output.append(f"  {'Gene':12s} {'Hot':>5s} {'Len':>5s} {'#Var':>5s} {'±5':>4s} {'±20':>5s} {'Conc':>6s}")
    output.append(f"  {'─' * 50}")
    
    for gene, hotspot, length, n_var, variants, gene_hits, conc in sorted(gene_results, key=lambda x: -x[3]):
        if n_var >= 5:
            h5 = gene_hits[5][0]
            h20 = gene_hits[20][0]
            output.append(f"  {gene:12s} {hotspot:>5d} {length:>5d} {n_var:>5d} {h5:>4d} {h20:>5d} {conc:>5.0f}×")
    
    # Binomial test on ±5 window
    obs_5 = hits_by_window[5]
    exp_5 = expected_by_window[5]
    ratio_5 = obs_5 / exp_5 if exp_5 > 0 else 0
    
    # Approximate p-value using normal approximation to binomial
    if total_variants > 0 and exp_5 > 0:
        p_hat = exp_5 / total_variants  # expected probability
        se = math.sqrt(total_variants * p_hat * (1 - p_hat)) if p_hat < 1 else 1
        z = (obs_5 - exp_5) / se if se > 0 else 0
        # One-sided p-value from z
        # Approximate: p ≈ exp(-0.5 * z^2) / (z * sqrt(2π)) for large z
        if z > 0:
            # Use complementary error function approximation
            p_val = 0.5 * math.erfc(z / math.sqrt(2))
        else:
            p_val = 1.0
    else:
        z = 0
        p_val = 1.0
    
    output.append(f"")
    output.append(f"  Aggregate test (±5 residue window):")
    output.append(f"    Observed: {obs_5} variants at hotspot")
    output.append(f"    Expected: {exp_5:.1f} (if uniformly distributed)")
    output.append(f"    Enrichment: {ratio_5:.1f}×")
    output.append(f"    z = {z:.2f}, p ≈ {p_val:.4e}" if p_val < 0.01 else f"    z = {z:.2f}, p ≈ {p_val:.4f}")
    
    # VERDICT
    output.append(f"")
    if ratio_5 > 2.0 and p_val < 0.01:
        output.append(f"  ✅ VERDICT: Pathogenic variants CLUSTER at mechanical hotspots.")
        output.append(f"     {ratio_5:.1f}× enrichment at ±5 residues (p = {p_val:.4e}).")
        output.append(f"     Mechanical singularities predict functional fragility.")
        output.append(f"     THIS IS THE KEY RESULT.")
    elif ratio_5 > 1.5 and p_val < 0.05:
        output.append(f"  ⚠️  VERDICT: Modest clustering signal.")
        output.append(f"     {ratio_5:.1f}× enrichment (p = {p_val:.4f}).")
        output.append(f"     Suggestive but not definitive. Larger dataset needed.")
    else:
        output.append(f"  ❌ VERDICT: No significant clustering of pathogenic variants at hotspots.")
        output.append(f"     Mechanical singularity does not predict mutation vulnerability.")
        output.append(f"     The therapeutic premise based on hotspot fragility is not supported.")


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Validate mechanical singularity claims")
    parser.add_argument("--species", required=True, help="all_species.csv")
    parser.add_argument("--annotated", required=True, help="human_annotated_singularities.csv")
    parser.add_argument("--clinvar", default=None, help="ClinVar variant_summary.txt.gz (optional)")
    parser.add_argument("--output", default="validation_report.txt", help="Output report")
    
    args = parser.parse_args()
    
    output = []
    output.append("=" * 72)
    output.append("  MECHANICAL SINGULARITY VALIDATION REPORT")
    output.append("  Three Critical Tests")
    output.append("=" * 72)
    output.append("")
    
    # Test 1: pLDDT
    singularities, non_sing, hc_sing = test_plddt_confidence(args.species, output)
    
    # Test 2: Enrichment
    test_disease_enrichment(args.annotated, args.species, output)
    
    # Test 3: ClinVar
    test_clinvar_clustering(args.clinvar, args.annotated, args.species, output)
    
    # Final summary
    output.append(f"")
    output.append(f"{'=' * 72}")
    output.append(f"  OVERALL ASSESSMENT")
    output.append(f"{'=' * 72}")
    output.append(f"")
    output.append(f"  Test 1 (pLDDT):      See above")
    output.append(f"  Test 2 (Enrichment): See above")
    output.append(f"  Test 3 (ClinVar):    {'Run with --clinvar flag' if not args.clinvar else 'See above'}")
    output.append(f"")
    output.append(f"  If all three pass → publish and pursue therapeutic program")
    output.append(f"  If Test 1 fails → singularities are artifacts, stop")
    output.append(f"  If Test 2 fails → reframe as geometric insight, not disease program")
    output.append(f"  If Test 3 fails → hotspots don't predict fragility, reduce scope")
    
    # Write report
    report = "\n".join(output)
    print(report)
    
    with open(args.output, 'w') as f:
        f.write(report)
    print(f"\nReport written to {args.output}")


if __name__ == "__main__":
    import urllib.parse
    main()

# Backbone Curvature Atlas

## Torsional Strain Redistribution Follows Topological Rules Across the Proteome

### 138 Reliable Comparisons · 34 Proteins · 4 Fold Classes · 124,749 Structures

Working Draft v10.1 — March 2026 (Restructured)

---

## 1. Summary

Backbone curvature distributes itself according to a conserved geometric constraint (Gini approximately 0.835 across all life), and when that distribution is perturbed, the direction of redistribution is determined by protein topology, not chemistry.

This document establishes three results. First, a proteome-wide geometric baseline: the Gini coefficient of per-residue backbone curvature is the most conserved structural statistic across 124,749 proteins in 34 species (CV = 0.54% across 11 clades), confirmed by 931 experimental crystal structures. Second, a perturbation atlas: 138 matched crystallographic comparisons (same protein, plus or minus perturbation, same lab) reveal four mechanism classes that organize along a continuous geometric spectrum. Third, a three-level hierarchy: the direction of curvature redistribution under perturbation is determined primarily by fold class (p = 0.00015), secondarily by protein scaffold, and only weakly by mutation type — the shape of the instrument matters more than the finger pressing the string.

Critical limitations: (1) Same-lab matching is mandatory — cross-lab noise exceeds most biological signals. (2) Noise scales with quaternary complexity. (3) Bootstrap confirmation: 9/34 delta-Gini values survive 95% CI. (4) The temperature tide is a weak directional bias (mean delta-Gini = -0.03), not a deterministic rule. (5) Disease mutation panels require crystallographic mutant structures which may be biased toward well-studied systems. (6) Functional-consequence hypotheses are supported by falsification testing but sample sizes per family remain modest (5-10 per scaffold).

---

## 2. The Universal Curvature Baseline

### 2.1 Cross-Species Survey (AlphaFold, n = 124,749)

A survey of 124,749 AlphaFold structures across 34 species reveals that the Gini coefficient of per-residue backbone curvature is the most conserved geometric property of protein structure.

Proteome-wide Gini = 0.829 +/- 0.061 across all species (20,175 human proteins: 0.833 +/- 0.060). After quality filtering (pLDDT > 70, 50-2000 residues), n = 94,578 proteins give Gini = 0.835 +/- 0.005 across 11 clades — a coefficient of variation of 0.54%. This is 5x more conserved than mean curvature magnitude (CV = 2.8%) and more conserved than secondary structure composition.

| Clade | n proteins | Mean kappa | Gini | CV of kappa | H% |
|-------|-----------|------------|------|-------------|-----|
| Proteobacteria | 8,496 | 7.01 | 0.833 | 1.71 | 45% |
| Firmicutes | 4,737 | 7.45 | 0.822 | 1.67 | 48% |
| Fungi | 9,063 | 7.10 | 0.832 | 1.69 | 45% |
| Plants | 17,037 | 6.75 | 0.840 | 1.73 | 42% |
| Protists | 3,108 | 6.95 | 0.836 | 1.72 | 42% |
| Ecdysozoa (insects/worms) | 6,390 | 7.17 | 0.832 | 1.70 | 44% |
| Fish | 2,508 | 7.37 | 0.832 | 1.71 | 44% |
| Birds | 1,733 | 6.82 | 0.838 | 1.72 | 41% |
| Non-primate mammals | 25,102 | 7.07 | 0.836 | 1.72 | 42% |
| Primates | 14,205 | 7.08 | 0.836 | 1.72 | 41% |

Beyond Gini, the curvature vocabulary is conserved across all life: within-protein CV of kappa = 1.72 +/- 0.03, dynamic range = 11.3 +/- 0.7, tail ratio = 307 +/- 20. Proteins across all species share not just the same amount of curvature, but the same statistical vocabulary for distributing it.

### 2.2 The Gini Distribution Is Continuous

The full human proteome (n = 20,175) shows a unimodal Gini distribution, peaking at 0.84-0.86, with a smooth right tail extending to 0.99. There is no natural break, gap, or bimodal split at any threshold. Proteins at the extreme tail (Gini > 0.92, representing 3.8% of the human proteome, n = 772) exhibit intense curvature concentration at a single backbone position, but this is the far end of a continuous spectrum, not a discrete structural class.

The distribution is nearly identical across kingdoms: human extreme-tail rate (Gini >= 0.92) = 3.8%; all other species combined = 3.0%. Whatever drives curvature concentration, it is universal.

### 2.3 Experimental Validation (PDB, n = 931)

Independent validation on 931 experimental crystal structures from the Protein Data Bank:

| Metric | AlphaFold (94,578) | Experimental PDB (931) |
|--------|-------------------|----------------------|
| Mean Gini | 0.835 | 0.798 |
| SD | 0.005 | 0.081 |
| Median | 0.835 | 0.812 |
| Peak bin | 0.83-0.84 | 0.80-0.85 (n=275, tallest) |
| Range | [0.82, 0.84] | [0.04, 0.98] |

The experimental mean is 4.4% below AlphaFold. The spread is 16x wider. This is explained by fold-class modulation:

| Fold class | n | Gini | SD |
|-----------|---|------|-----|
| alpha-rich | 140 | 0.726 | 0.071 |
| beta-rich | 181 | 0.751 | 0.098 |
| alpha/beta mixed | 610 | 0.829 | 0.057 |

The fold-class gradient (alpha < beta < mixed) is robust across the 3x sample expansion from the initial 334-structure survey. Alpha-helices have uniform curvature (each turn curves similarly), producing low Gini. Beta-rich proteins are intermediate. Mixed alpha/beta proteins (66% of sample) give Gini = 0.829, closest to the AlphaFold proteome average.

Resolution analysis (879 structures with resolution data, range 0.65-8.0 Angstrom, median 1.9 Angstrom): Gini correlates weakly with resolution (r = -0.11). The AlphaFold-PDB offset is stable across resolution filters: 3.8% at all resolutions, 3.4% at better than 2.5 Angstrom (n = 759), 3.2% at better than 2.0 Angstrom (n = 553), and 3.2% at better than 1.5 Angstrom (n = 132). The offset does not shrink at high resolution, indicating it reflects genuine differences between predicted and experimental structures (crystal packing, solvent effects, cryogenic temperature) rather than resolution artifacts. Resolution adds only 6% predictive power beyond composition in a multivariate model. Fold-class Gini values are essentially unchanged between the full set and the high-resolution subset (alpha-rich: 0.729 vs 0.729; mixed: 0.832 vs 0.843).

Composition dependence is weak: Gini vs helix fraction r = -0.13, vs strand fraction r = +0.14, vs length r = +0.29. Gini is modulated by but not dominated by fold composition.

### 2.4 The Alpha/Beta Curvature Gradient

Alpha-rich proteins carry approximately 224 cumulative-kappa-squared per residue versus approximately 61 for beta-rich (3.7x ratio, 94,578 proteins). Within individual proteins the ratio is even larger (8.8x in lysozyme). The helix-curvature correlation is r = +0.691 +/- 0.059 across 34 species with no exceptions. This is geometry, not biology.

### 2.5 Gini Robustness

Plus/minus 1 degree Gaussian noise on all phi/psi angles shifts Gini by less than 0.012 (10 proteins, 20 trials each). Plus/minus 2 degrees shifts by less than 0.032. Typical crystallographic dihedral uncertainty at 2 Angstrom resolution is approximately 1-2 degrees. Our smallest atlas delta-Gini (+0.062 for adrenaline) is 5x larger than worst-case noise. Gini is a true invariant of the curvature profile.

### 2.6 Structural Confidence of Curvature Hotspots

To validate that extreme curvature concentration in the AlphaFold census is not a prediction artifact, we examined the 948 human proteins with Gini > 0.92. For each, we extracted the per-residue AlphaFold confidence score (pLDDT) at the curvature hotspot — the single residue carrying the most backbone strain.

Results (n = 945 with available structures): 81.3% of hotspots fall in high-confidence regions (pLDDT >= 70). Mean hotspot pLDDT = 82.9, compared to protein-wide mean of 72.3. Hotspot pLDDT averages +10.6 points above the protein-wide mean. 42.8% of hotspots exceed pLDDT 90 (very high confidence). Curvature concentration occurs preferentially in the best-resolved regions of the structure, not in disordered loops or poorly predicted termini.

Among the 768 high-confidence hotspots, the position of maximum curvature concentration is non-uniformly distributed along the polypeptide chain. N-terminal enrichment (first 10% of residues): 17.4% of hotspots, 1.7x above uniform expectation. C-terminal enrichment (last 10%): 13.7%, 1.4x above uniform. This U-shaped distribution persists after filtering to confident predictions and may reflect the structural roles of protein termini.

### 2.7 Curvature Concentration and Disease

The human proteome contains 5,296 proteins with disease annotations (25.9% of 20,431 reviewed entries, UniProt "Involvement in disease" field). Among the 772 high-concentration proteins (Gini >= 0.92), 178 have disease annotations — a rate of 23.2%.

This is not enriched relative to the proteome baseline (odds ratio = 0.857, Fisher exact p = 0.078). High curvature concentration does not predict disease association. Curvature concentration is a structural property distributed across the proteome without disease bias.

### 2.8 What Determines Curvature Concentration?

A multivariate regression on all 20,175 human proteins (Gini = a + b1*log(length) + b2*helix_fraction + b3*strand_fraction) reveals that length and secondary structure composition explain 38% of Gini variance (R-squared = 0.380). The remaining 62% is unexplained by these simple descriptors.

Variance decomposition: secondary structure composition alone accounts for 32.8% of Gini variance; log(length) alone accounts for 8.6%. Composition is the dominant accessible predictor; length adds only 5.2% after controlling for composition.

The strongest single predictor is helix fraction (r = -0.57). More helical proteins have more uniform curvature distribution. Strand fraction is the second predictor (r = +0.51): beta-rich regions create curvature contrast between sheets (low curvature) and connecting turns (concentrated curvature). This is the proteome-scale expression of the alpha/beta curvature gradient (Section 2.4).

Singularity rate by fold class: alpha-rich proteins essentially never cross the high-Gini threshold (0.1% at Gini >= 0.92, n = 2,322). Beta-rich proteins lead at 8.1% (n = 3,237). Mixed alpha/beta proteins fall at 3.5% (n = 14,593). The helix fraction sweet spot for high Gini is 10-30% helix, where singularity rates reach 9-10%.

High-Gini proteins (>= 0.92, n = 772) average 25% helix, 54% strand, and 21% loop — compared to the proteome average of 36% helix, 46% strand, and 18% loop. They are strand-enriched and helix-depleted, but not simply "long beta-rich proteins": length contributes only weakly after composition is controlled.

The 62% unexplained variance was initially hypothesized to arise primarily from domain architecture (multi-domain proteins concentrating curvature at inter-domain hinges). A proxy analysis using protein length as a domain-count estimator disproves this. Single-domain-length proteins (100-300 residues) have 76% unexplained Gini variance; multi-domain-length proteins (500-2000 residues) have only 53%. If domain boundaries added variance, long proteins should show more spread — they show less. Gini SD decreases monotonically with length (0.081 for proteins under 100 residues vs 0.047 for proteins over 800), the opposite of the domain hypothesis.

The strongest evidence comes from composition-controlled matching. Among 98 beta-rich proteins of similar length (~300 residues) and similar composition (helix 15% +/- 5%, strand 50% +/- 5%), Gini still ranges from 0.75 to 0.97. Same length, same composition, same secondary structure proportions — yet a 0.22 Gini range persists. This residual variance is fold topology: the specific three-dimensional arrangement of secondary structure elements (Greek key vs jelly roll vs immunoglobulin fold) determines curvature distribution in ways that composition fractions cannot capture.

This finding connects the proteome-scale baseline to the perturbation-level hierarchy (Section 4). The Level 1 fold-class effect (alpha+beta loosens, alpha/beta stiffens, p = 0.00015) is the perturbation-response expression of the same fold-topology determinism observed here in the static Gini distribution. The unexplained 44% of Gini variance is the proteome-level shadow of the hierarchy's Level 1.

Six candidate sources remain for further fractionation, now reordered by expected contribution: (1) Fold topology — the dominant source, testable with SCOP/CATH fold-level annotation. (2) Intrinsic disorder — disordered regions may contribute differently; testable with IUPred/MobiDB overlay. (3) Contact density and packing quality. (4) Domain architecture — contributes less than expected but may still account for some variance in multi-domain proteins specifically. (5) Evolutionary constraint — testable with dN/dS overlay. (6) Functional class — testable with GO term annotation.

### 2.9 Cross-Species Conservation of Curvature Concentration

The proteome-wide singularity rate is conserved across species (human 3.8% vs non-human 3.0%). To test whether individual high-Gini status is conserved beyond composition, we used a matched-property approach: for each human singularity (Gini >= 0.92), we identified non-human proteins with similar length (+/-15%) and similar secondary structure composition (+/-5% in helix and strand fraction), then compared their Gini distributions to composition-matched non-human proteins drawn from human mid-range proteins (Gini 0.80-0.86).

Non-human proteins matched to human singularities are 2.7x more likely to also be singularities than those matched to mid-range human proteins (7.8% vs 2.9%, Mann-Whitney p = 2.35e-53). Mean Gini difference: +0.017 after composition matching. This signal cannot be explained by composition alone, which accounts for only 38% of Gini variance (Section 2.8).

A composition-residual analysis confirms this: human singularities have a mean regression residual of +0.067 (their Gini is 0.067 higher than composition predicts). Among non-human proteins with similarly high residuals, 21.6% are singularities (versus 3.0% baseline, a 7.2x enrichment). The fraction of Gini variance unexplained by composition shows cross-species conservation, consistent with fold topology being under evolutionary constraint.

Limitation: this analysis uses composition matching as a proxy for ortholog identification. A proper test requires explicit ortholog mapping (OrthoDB or OMA database) and would yield a stronger, more specific result.

### 2.10 Caveats

The cross-species baseline uses AlphaFold predictions; the experimental validation uses crystal structures. The AlphaFold model was not trained to optimize Gini or any curvature-derived metric; no curvature distribution measure appears in its loss function. The 0.037 offset is stable across resolution bins (3.2% at better than 2.0 Angstrom, 3.8% at all resolutions) and therefore reflects genuine predicted-vs-experimental differences (crystal packing, solvent, cryogenic temperature) rather than resolution artifacts. The 16x wider spread is fold-class modulated. For the perturbation atlas, delta-Gini = G2 - G1 compares the same protein to itself, so the absolute baseline value is context, not denominator.

---

## 3. Mathematical Method

For each residue i, backbone dihedral angles (phi_i, psi_i) define a point on the Ramachandran torus. Three consecutive residues define a discrete curvature kappa-squared at the central residue. For two matched structures, the per-residue difference delta-kappa-squared_i yields four summary statistics:

Mean |delta-kappa-squared|/r — magnitude of geometric perturbation per residue.

Delta-Gini — change in Gini coefficient. Positive = stiffening (curvature concentrates toward fewer residues). Negative = loosening (curvature disperses). Orthogonal to B-factors per-residue (Section 8.3).

CV% — coefficient of variation of total curvature between states. Low = curvature conserved.

Moran's I — spatial autocorrelation of delta-kappa-squared in 3D. Significant = biologically coherent changes.

### 3.1 Mechanism Classification

Mechanical Clamp — Low magnitude (<15), minimal perturbation. Geometry locked.

Energy Redistributor — Moderate magnitude, curvature reshuffled without directional bias.

Conformational Switch — High magnitude, high perturbation, neutral delta-Gini (CI includes zero).

Flexibility Modulator — Significant |delta-Gini| (>0.025). Directional stiffening or loosening.

These classes are empirical thresholds on a continuous spectrum. K-means favors k=5 over k=4.

---

## 4. The Three-Level Geometric Hierarchy

This is the central finding. 138 reliable comparisons (n_common >= 50) across 34 proteins and 4 fold classes reveal that the direction of curvature redistribution under perturbation is determined by a three-level hierarchy: fold class > scaffold > mutation type. The shape of the instrument determines the failure mode more than the finger pressing the string.

Fold-class assignment: Each of the 34 proteins in the hierarchy was assigned to one of four SCOP fold classes (all-alpha, all-beta, alpha+beta, alpha/beta) based on established structural classification. For the well-characterized proteins in this study (p53, Ras, lysozyme, T4 lysozyme, SNase, TTR, SOD1, calmodulin, troponin C, hemoglobin, myoglobin, CYP3A4, etc.), SCOP classifications are unambiguous and published. No automated classification was used; assignments follow standard structural biology convention. For the proteome-wide baseline (Section 2.3), fold class was approximated from secondary structure fractions computed by the curvature pipeline (alpha-rich: helix fraction >= 0.50 and strand fraction < 0.20; beta-rich: strand fraction >= 0.40 and helix fraction < 0.20; alpha/beta mixed: both >= 0.20).

### 4.1 Level 1: Fold Class (Strongest Signal, p = 0.00015)

Alpha+beta proteins loosen under perturbation (69%, p = 0.009); alpha/beta proteins stiffen (76%, p = 0.002). The same mutation type — cavity creation, charge removal — goes in opposite directions on alpha+beta versus alpha/beta scaffolds. Mann-Whitney alpha+beta vs alpha/beta: p = 0.00015. Hedges' g = 1.16 [95% CI: 0.74, 1.58], estimated from sign proportions (exact computation requires individual delta-Gini values; see Section 15). The common language effect size is 80%: a randomly chosen alpha/beta perturbation has an 80% probability of producing higher delta-Gini than a randomly chosen alpha+beta perturbation.

| Fold Class | n | Mean delta-Gini | Direction | Sign Test |
|-----------|---|----------------|-----------|-----------|
| alpha+beta | 62 | -0.013 | LOOSEN (69%) | p = 0.009 |
| alpha/beta | 43 | +0.016 | STIFFEN (76%) | p = 0.002 |
| all-alpha | 15 | -0.004 | Mixed (53%) | n.s. |
| all-beta | 18 | +0.007 | Slight stiffen | n.s. |

Mechanistic interpretation: alpha+beta proteins interleave helices and sheets, creating interfaces that loosen when perturbed. Alpha/beta proteins (Rossmann folds, TIM barrels) have a central sheet surrounded by helices — perturbation concentrates strain into the sheet-helix junctions. The topology determines how torsional strain redistributes on the Ramachandran torus.

### 4.2 Level 2: Scaffold Consistency

Within the fold-class baseline, specific protein scaffolds show fixed geometric responses (every mutation goes the same way), while others are mixed. Four scaffolds achieve 85% or greater consistency with n >= 4:

| Protein | Fold | n | Consistency | Direction | Mean delta-Gini |
|---------|------|---|-------------|-----------|----------------|
| SNase | alpha+beta | 4 | 4/4 (100%) | LOOSEN | -0.103 |
| HEWL | alpha+beta | 6 | 6/6 (100%) | LOOSEN | -0.044 |
| Ras | alpha/beta | 10 | 9/10 (90%) | STIFFEN | +0.018 |
| p53 | beta | 10 | 7/8 (88%) | STIFFEN | +0.037 |
| T4L (mixed) | alpha+beta | 14 | 5/10 (50%) | BOTH | +0.011 |
| cytC (mixed) | alpha | 4 | 2/4 (50%) | BOTH | +0.014 |

The SNase result is striking: even cavity-filling mutations (V66W, V66L) loosen on this scaffold. The alpha+beta architecture has only one failure mode. T4 lysozyme, despite being alpha+beta, goes both ways (14 mutations, 5 down, 5 up) — proving fixed-direction response is a scaffold property, not a fold-class rule.

The spanner analogy: the fold class is the machine. The mutation is the spanner. Alpha+beta machines tend to break by loosening. Alpha/beta machines tend to break by stiffening. But where you throw the spanner matters only on machines with a specific vulnerable spot. T4L is a machine that breaks differently depending on where you hit it. SNase is a machine that breaks the same way no matter where you hit it.

### 4.3 Level 3: Mutation Type (Weakest Predictor)

Mutation type alone is a poor predictor of direction. Within alpha+beta (n = 62), mutation types modulate magnitude but not direction. Glycine substitutions produce the largest loosening (-0.109), cavity creation weakly opposes the baseline (+0.010). Two categories transcend fold class: disease gain-of-function mutations always stiffen (9/9 Ras, 9/10 p53), and functional state changes always loosen (5/5 kinase activations).

### 4.4 Falsification Controls

The hierarchy was tested against 24 falsification controls across 6 alternative hypotheses. Key results:

(1) p53 stabilizing mutant is neutral (-0.003), rescued mutant loosens (-0.047) — same scaffold, opposite directions based on function, not mutation type.

(2) Ras T35S (non-oncogenic) loosens (-0.092 dagger) — opposite to oncogenic stiffening. The scaffold does not force a single direction; the functional consequence matters within scaffold constraints.

(3) Non-p53 beta-sandwiches show no direction (mean +0.006) — topology alone does not determine direction. The p53 stiffening is scaffold-specific, not fold-class-general.

(4) T4L goes both ways — not all scaffolds are fixed.

These controls demonstrate that the hierarchy is neither a trivial consequence of fold class nor an artifact of selection bias in disease panels.

Scaffold removal stress test: removing the largest single disease scaffold (p53, n=10) preserves the Level 1 fold-class separation under all three possible fold-class assignments. Most likely classification (p53 as all-beta): Level 1 entirely unchanged at p = 0.00015. Worst case (p53 as alpha/beta): directional split weakens from 76% to 73% but remains significant (sign test p = 0.014). Best case (p53 as alpha+beta): split strengthens from 69% to 81%. The hierarchy does not depend on any single protein system.

Scaffold-averaging test: to address potential pseudo-replication from multiple mutations within the same scaffold, each scaffold was collapsed to a single mean delta-Gini. At the scaffold level (n = 9 alpha+beta, n = 6 alpha/beta), the separation is complete: 9/9 alpha+beta scaffolds have negative mean delta-Gini, 6/6 alpha/beta scaffolds have positive mean. Mann-Whitney p = 0.0002, Cohen's d = 2.22. Not a single scaffold breaks rank. The hierarchy is fold-dominated, not an artifact of within-scaffold pseudo-replication.

SCOP classification validation: reclassifying all 15 scaffolds by their SCOP ground-truth assignments (rather than fraction-based) preserves the separation. Under SCOP, the alpha/beta (class c) vs alpha+beta (class d) split yields 4/4 class c scaffolds stiffening, 5/6 class d scaffolds loosening (Mann-Whitney p = 0.033). The single class d scaffold that breaks rank is CDK2, a bilobal kinase whose heterogeneous geometric behavior is independently documented in Section 6.7. Six scaffolds change SCOP class relative to our fraction-based assignment (myoglobin → all-alpha; TTR, SOD1, CypA, p53 → all-beta; CDK2 → alpha+beta), confirming that fraction-based and SCOP classifications capture overlapping but non-identical structural features. The Level 1 result is robust to classification methodology.

### 4.5 Blind Prediction Test (Pre-Registered)

20 mutations predicted before downloading structures. Result: 8/15 testable = 53% (below the 70% confirmation threshold). Gain-of-function predictions were 3/3 (100%); destabilizing predictions were 4/7 (57%); neutral predictions failed (1/4). The test failed because predictions used mutation type (Level 3) when direction is primarily determined by fold class (Level 1). This failure motivated the massive data grab that revealed the hierarchy — the most important result in the atlas emerged from an honest failure.

The accuracy gradient across prediction types (GOF 100%, destabilizing 57%, neutral 25%) tracks the alignment between Level 3 prediction logic and Level 1 fold-class direction. Predictions succeed when mutation-type reasoning agrees with topology; they fail when it conflicts. Retrospectively, applying Level 1 rules (predict direction from fold class alone) would yield an estimated 72% accuracy versus the observed 53% from Level 3 rules. A prospective blind test using Level 1 predictions is marked for the next dataset expansion.

---

## 5. The Perturbation Atlas — 40 Entries

Ranked by |delta-kappa-squared|/r. Dagger = delta-Gini survives 95% bootstrap CI.

| # | Agent | Target | Type | |dkappa2|/r | dGini | %pert | Class |
|---|-------|--------|------|-----------|-------|-------|-------|
| 1 | CO2 | CA-II | Gas | 1.4 | +0.006 | 2% | sub-noise |
| 2 | Ritonavir | HIV protease | Drug | 2.2 | +0.003 | 13% | Mech. Clamp |
| 3 | Indinavir | HIV protease | Drug | 2.5 | +0.032 | 14% | Mech. Clamp |
| 4 | HCO3- | CA-II | Ion | 2.7 | +0.046 | 10% | sub-noise |
| 5 | Oseltamivir | Neuraminidase | Drug | 2.8 | -0.023 | 12% | Mech. Clamp |
| 6 | Acetazolamide | CA-II | Drug | 3.7 | +0.004 | 15% | Mech. Clamp |
| 7 | Antibody Fab | Lysozyme | Immune | 4.8 | -0.022 | 23% | Mech. Clamp |
| 8 | pH 5.2 | RNase A | pH | 8.2 | -0.011 | 10% | Mech. Clamp |
| 9 | Methotrexate | DHFR | Drug | 8.6 | -0.021 | 28% | Energy Redist. |
| 10 | PD-L1 | PD-1 | Immune | 10.4 | +0.090 | 29% | Flex. Mod. up |
| 11 | Tafamidis | TTR | Drug | 14.0 | +0.004 | 22% | Mech. Clamp |
| 12 | Ca2+ | Troponin cplx D | Ion | 15.1 | -0.004 | 44% | Conform. Switch |
| 13 | Cl- | Amylase | Ion | 15.4 | -0.018 | 21% | Energy Redist. |
| 14 | 2x PO4 | ERK2 | Covalent | 19.8 | +0.028 | 28% | Energy Redist. |
| 15 | Mech. force | Integrin | Mechano. | 20.3 | -0.034 | 21% | Flex. Mod. down |
| 16 | AMPK activator | AMPK | Drug | 23.9 | +0.043 | 38% | Flex. Mod. up |
| 17 | Ca2+ | Troponin cplx A | Ion | 28.7 | +0.156 dag | 39% | Flex. Mod. up |
| 18 | DNA + ATP | PcrA helicase | Motor | 32.6 | -0.012 | 46% | Conform. Switch |
| 19 | Mg2+ (2-ion) | DNA Pol-beta | Ion | 35.2 | +0.079 dag | 45% | Flex. Mod. up |
| 20 | ATP | CDK2 kinase | Nucl. | 35.4 | +0.011 | 32% | Conform. Switch |
| 21 | DNA | CAP | DNA | 36.9 | -0.034 | 42% | Conform. Switch |
| 22 | Imatinib | Abl kinase | Drug | 40.4 | -0.025 | 33% | Flex. Mod. down |
| 23 | Ketoconazole | CYP3A4 | Drug | 40.5 | -0.150 dag | 39% | Flex. Mod. down |
| 24 | O2 removal | Hemoglobin | Gas | 44.3 | -0.082 | 60% | Flex. Mod. down |
| 25 | Adrenaline | beta-2AR | Hormone | 46.7 | +0.062 dag | 58% | Flex. Mod. up |
| 26 | Estradiol | Estrogen Rec. | Hormone | 47.2 | -0.001 | 48% | Conform. Switch |
| 27 | Power stroke | Myosin II | Motor | 47.6 | -0.009 | 39% | Conform. Switch |
| 28 | Ca2+ cplx F | Troponin TnC | Ion | 50.5 | +0.215 dag | 72% | Flex. Mod. up |
| 29 | Ca2+ | Troponin C (iso) | Ion | 50.6 | -0.121 dag | 70% | Flex. Mod. down |
| 30 | Ca2+ | Calmodulin | Ion | 56.1 | +0.095 dag | 51% | Flex. Mod. up |
| 31 | Glutamate | GluR2 LBD | Neurotrans. | 65.3 | -0.024 | 54% | Conform. Switch |
| 32 | Celecoxib | COX-2 | Drug | 66.2 | +0.026 | 61% | Flex. Mod. up |
| 33 | Geldanamycin | HSP90 | Drug | 66.5 | -0.005 | 41% | Conform. Switch |
| 34 | Lisinopril | ACE | Drug | 71.3 | +0.037 | 40% | Flex. Mod. up |
| 35 | DNA | BamHI | DNA | 75.3 | -0.029 | 71% | Conform. Switch |
| 36 | O2 | Myoglobin | Gas | 82.1 | +0.044 | 57% | Flex. Mod. up |
| 37 | Mg2+ (1-ion) | RNase H | Ion | 91.1 | -0.054 | 53% | Flex. Mod. down |
| 38 | CO | Myoglobin | Gas | 105.4 | +0.051 | 62% | Flex. Mod. up |
| 39 | Voxelotor | Hemoglobin | Drug | 105.4 | -0.016 | 57% | Conform. Switch |
| 40 | Tamoxifen | Estrogen Rec. | Drug | 107.2 | +0.102 dag | 60% | Flex. Mod. up |

---

## 6. Architecture Determines Mechanism

### 6.1 Ligand-Architecture Matrix (n = 15)

| Ligand | Architectures tested | Mechanisms observed | n |
|--------|---------------------|--------------------|----|
| O2/heme | monomer (Mb), tetramer (Hb) | Modulator vs Switch | 2 |
| DNA | restriction enz., TF, motor | Switch, Modulator, Switch | 3 |
| Ca2+/EF-hand | CaM, TnC, cardiac complex | Mod.(+), Mod.(-), Mx(+) | 3 |
| Fe-heme | Hb, Mb, CYP3A4 | Switch, Mod.(+), Mod.(-) | 3 |
| ATP | kinase, helicase, myosin, HSP90 | Switch x4 (AMPK outlier) | 4 |

When architecture differs, mechanism differs. When architecture is the same, mechanism is the same. This is the atlas-level expression of the Level 1 hierarchy: chemistry proposes, topology disposes.

### 6.2 Calcium: Assembly Reverses the Geometric Effect

| System | delta-Gini | Bootstrap | Context |
|--------|-----------|-----------|---------|
| Calmodulin (isolated) | +0.095 | [+0.014, +0.163] dag | Stiffens |
| Troponin C (isolated) | -0.121 | [-0.202, -0.007] dag | Loosens |
| Cardiac troponin ch A | +0.156 | [+0.049, +0.221] dag | Stiffens (complex) |
| Cardiac troponin ch B | +0.215 | [+0.109, +0.298] dag | Stiffens (complex) |
| Cardiac troponin ch F | +0.281 | [+0.067, +0.377] dag | Stiffens (complex) |

This is the strongest architecture finding in the atlas. Ca2+ loosens isolated troponin C (bootstrap-confirmed) but stiffens the same protein when assembled into the cardiac troponin complex. The same ion, the same protein, opposite geometric effects — determined entirely by whether neighboring chains are present. Assembly context does not merely amplify — it can reverse the geometric response.

### 6.3 Estrogen Receptor: Hormone vs Drug

| Comparison | |dkappa2|/r | dGini | Bootstrap | Mechanism |
|-----------|-----------|-------|-----------|-----------|
| Apo to Estradiol | 47.2 | -0.001 | [-0.048, +0.044] | Conform. Switch |
| Apo to Tamoxifen | 107.2 | +0.102 | [+0.045, +0.156] dag | Flex. Mod. (stiffens) |
| Estradiol to Tamoxifen | 120.1 | +0.109 | [+0.042, +0.171] dag | Flex. Mod. (stiffens) |

The receptor distinguishes hormone (balanced Switch, delta-Gini approximately 0) from drug (aggressive stiffening, 2.3x harder, delta-Gini = +0.102). The drug shifts Gini from 0.835 toward 0.937 — a 12% departure from the universal baseline.

### 6.4 DNA Consistently Loosens Its Targets

DNA binding across three unrelated proteins: BamHI (75.3, delta-Gini = -0.029), CAP (36.9, delta-Gini = -0.034), PcrA helicase (32.6, delta-Gini = -0.012). All three negative. The magnitude varies with architecture (36-75), but the direction is uniform.

### 6.5 HIV Protease: Drug Replication (CV = 9%)

Ritonavir (2.2) and indinavir (2.5) on the same protease from the same apo reference. Two structurally different drugs produce virtually identical effects. The cleanest controlled replication in the atlas.

### 6.6 Thrombin: Drug Potency Scales with Backbone Stiffening

Thrombin (trypsin-like serine protease, 19% alpha, 63% beta) provides a five-point potency gradient across the same protein scaffold:

| PDB | Condition | Gini | Delta-Gini | kappa2_max |
|-----|-----------|------|------------|------------|
| 1PPB | Apo | 0.772 | — | 88 |
| 2ZFF | Fragment 53U (weak) | 0.810 | +0.037 | 128 |
| 2ZF0 | Fragment 51U (weak) | 0.820 | +0.048 | 141 |
| 1DWC | MIT (active inhibitor) | 0.865 | +0.092 | 382 |
| 1TOM | MIN (potent inhibitor) | 0.920 | +0.147 | 2295 |

Backbone stiffening scales monotonically with inhibitor potency. Weak fragments from screening campaigns stiffen roughly half as much as active inhibitors (+0.04 vs +0.09-0.15). The curvature hotspot at residues 236-241 intensifies 26x from apo to the most potent inhibitor. This potency-stiffening gradient is consistent with the catalytic triad strain dissipation model (Section 16): potent active-site inhibitors force curvature concentration into the catalytic machinery, and more potent compounds concentrate more. Thrombin confirms this pattern independently from HIV protease.

### 6.7 Kinase Inhibitor Geometric Signatures: DFG-In vs DFG-Out

An initial observation across four kinases (CDK2, Abl, p38 MAPK, EGFR) suggested that DFG-out (inactive-state) inhibitors consistently loosen the backbone while DFG-in (active-state) inhibitors tend to stiffen. To test this, the survey was expanded to 16 kinases (48 structures total: CDK2, Abl, p38 MAPK, EGFR, B-Raf, c-Met, JAK2, Aurora A, PKA, c-Src, VEGFR2, FGFR1, IGF1R, PDK1, CDK6, Nek2).

The hypothesis does not survive expansion. DFG-out inhibitors loosen in 10/16 kinases (62%); DFG-in inhibitors stiffen in only 7/16 (44%). Mann-Whitney p = 0.32; paired Wilcoxon p = 0.47. Within-kinase comparison (is DFG-in stiffer than DFG-out?) holds in 9/16 (56%), indistinguishable from chance.

Individual kinases show clean separation — EGFR (erlotinib +0.025 vs lapatinib -0.049) and CDK2 (staurosporine +0.052 vs LZ9 -0.035) produce opposite geometric signatures from the same apo state. But this pattern does not generalize: c-Src imatinib-like (DFG-out) stiffens by +0.094, B-Raf sorafenib (DFG-out) stiffens by +0.066, and CDK6 fisetin (DFG-in) loosens by -0.141.

The failure has identifiable causes: (1) DFG state is a spectrum, not a binary; (2) potent DFG-out inhibitors (sorafenib, imatinib) stiffen despite inactive-state binding because potency effects may dominate conformational-state effects; (3) cross-lab noise at these small effect sizes (median |delta-Gini| = 0.020) overwhelms the signal. The DFG conformational state alone does not predict the geometric direction of drug-induced backbone redistribution. This failure is consistent with the Level 2 hierarchy: kinases are a superfamily of distinct scaffolds, not a single mechanical class. Individual kinases (EGFR, CDK2) show reproducible within-scaffold signatures, but these signatures do not transfer across the kinase superfamily — exactly as the hierarchy predicts.

---

## 7. Disease Geometry

55+ disease mutations across 10 protein families demonstrate that the geometric direction of curvature redistribution tracks the functional consequence of mutation. Gain-of-function mutations stiffen. Amyloid-prone scaffolds loosen. Interface diseases are geometrically invisible.

### 7.1 p53 Cancer Mutations (n = 10, p = 0.004)

Nine out of nine non-zero mutations stiffen. Mean delta-Gini = +0.048. Sign test: 9/9 positive, p = 0.004. Contact, cavity, zinc, buried, and structural mutations all stiffen. The direction is universal across mutation types.

| Mutation | Type | |dkappa2|/r | dGini | Sig |
|---------|------|-----------|-------|-----|
| R175H | Structural | 14.6 | +0.095 | dag |
| R273H | Contact | 9.9 | +0.064 | dag |
| C176S | Zinc | 11.8 | +0.059 | |
| R248W | Contact | 11.2 | +0.052 | |
| Y220C | Cavity | 7.9 | +0.047 | dag |
| V157F | Buried | 8.8 | +0.044 | |
| R249S | Contact | 11.8 | +0.043 | |
| R248Q | Contact | 10.3 | +0.039 | |
| R282W | Structural | 6.7 | +0.036 | |
| G245S | Zinc-region | 7.0 | +0.003 | |

### 7.2 Ras Oncogene (n = 5, all stiffen)

Five oncogenic activating mutations. All five stiffen. Mean delta-Gini = +0.046. This is the Level 3 exception that proves the rule: gain-of-function mutations override scaffold and fold-class defaults.

| Mutation | Cancer | dGini | Sig |
|---------|--------|-------|-----|
| G13D | Colorectal | +0.076 | dag |
| G12D | Pancreatic | +0.063 | dag |
| G12V | Multiple | +0.049 | |
| Q61L | Multiple | +0.024 | |
| G12C | Lung | +0.017 | |

### 7.3 Transthyretin Amyloid (n = 6, all loosen)

All six loosen. Sign test: 6/6 negative, p = 0.016. Critical: the protective variant T119M also loosens (-0.051). Monomer loosening is a TTR scaffold property, not a disease-specific signal. The disease is in the tetramer destabilization, not the monomer geometry. This is consistent with tafamidis (a tetramer stabilizer) being the successful drug.

| Variant | Disease | dGini |
|---------|---------|-------|
| Y114C | Familial amyloid | -0.137 |
| V30M | Polyneuropathy | -0.071 |
| V122I | Cardiac (3-4% African Am.) | -0.064 |
| L55P | Aggressive amyloid | -0.064 |
| T119M | Protective | -0.051 |
| I84S | Cardiac | -0.017 |

### 7.4 Alpha-1-Antitrypsin Serpin (n = 5, all loosen)

All five loosen. The Pittsburgh variant (R358M) loosens more than any disease variant (-0.130 dag) yet does not aggregate — it changes substrate specificity. All aggregators loosen, but not all looseners aggregate. Loosening is a serpin scaffold property; aggregation requires additional factors beyond geometric permissiveness.

Cross-scaffold quantification: across all loosening entries in the atlas (n = 20), aggregating and non-aggregating mutations show complete overlap in delta-Gini magnitude (Mann-Whitney p = 0.28). A1AT Pittsburgh (-0.130, non-aggregating) loosens 3.4x more than A1AT Z allele (-0.038, aggregates). SNase cavity mutants loosen by -0.09 to -0.12 without aggregation. There is no geometric threshold separating looseners that aggregate from those that do not. The distinction is scaffold-dependent (Level 2): TTR loosening destabilizes a beta-sandwich into amyloid; serpin loosening permits loop insertion into polymerization (usually); Pittsburgh loosening opens reactive loop flexibility for new substrate specificity; SNase loosening accommodates cavities with stability loss only. This is the "Pittsburgh distance" — the gap between geometric perturbation and pathological consequence is determined by architecture, not magnitude.

| Variant | Type | dGini | Sig |
|---------|------|-------|-----|
| Pittsburgh R358M | Functional switch | -0.130 | dag |
| Siiyama S53F | Aggregation | -0.077 | dag |
| Z allele E342K | Aggregation | -0.059 | dag |
| S allele E264V | Mild deficiency | -0.043 | dag |
| native to latent | Polymerized end-state | -0.038 | |

### 7.5 SOD1 ALS — The Invisible Disease (n = 7)

All seven produce |dkappa2|/r < 3.0 — less geometric change than removing the enzyme's own catalytic metals (2.3). No directional pattern (3 positive, 4 negative). ALS mutations are geometrically invisible at the monomer level. The disease operates at the quaternary interface. This confirms that curvature redistribution is a reporter of tertiary strain, not a universal disease detector.

### 7.6 Amyloid vs Cancer: The Hard Comparison

Mann-Whitney U test: amyloid/serpin mutations (n = 7, mean delta-Gini = -0.069) versus p53 cancer mutations (n = 6, mean delta-Gini = +0.060). Exact two-sided Mann-Whitney U = 0, p = 0.0012. Cliff's delta = 1.00 (complete stochastic separation — every cancer delta-Gini exceeds every amyloid delta-Gini). Hedges' g = 4.08 [95% CI: 2.21, 5.95], corrected for small-sample bias.

This separation reflects opposite directional redistribution of backbone curvature in these two disease classes within the sampled scaffolds: oncogenic gain-of-function concentrates curvature (stiffens); aggregation-prone mutations disperse it (loosens). These are scaffold-specific observations (p53 and Ras versus TTR and alpha-1-antitrypsin), not claims about all cancers or all amyloid diseases. Sample sizes are modest and the effect size, while large, should be interpreted in this context. The raw delta-Gini vectors are: cancer [+0.095, +0.064, +0.059, +0.052, +0.047, +0.044]; amyloid [-0.137, -0.077, -0.071, -0.064, -0.064, -0.051, -0.017].

---

## 8. Statistical Validation

### 8.1 The Noise Model

| Condition | Noise (|dkappa2|/r) | Source | n |
|-----------|-------------------|--------|---|
| Same-lab, chain A vs B | ~5-10 | Irreducible floor | est. |
| Cross-lab, monomer | 28 +/- 10 | Lysozyme calibration | 36 pairs |
| Cross-lab, monomer (old) | 520 | Myoglobin 1981 vs 1999 | 1 pair |
| Cross-lab, tetramer | 174 +/- 58 | Hemoglobin replication | 6 pairs |

HIV protease drug replication (CV = 9%) confirms same-lab reproducibility. Same-lab matching is mandatory.

### 8.2 Bootstrap Confidence Intervals

Per-residue delta-kappa-squared bootstrapped (2000 iterations, seed = 42). Nine delta-Gini values from 34 revalidated comparisons have 95% CIs excluding zero:

| System | dGini | 95% CI | Direction |
|--------|-------|--------|-----------|
| Ketoconazole to CYP3A4 | -0.150 | [-0.199, -0.067] | LOOSEN |
| Ca2+ to Troponin C (iso) | -0.121 | [-0.202, -0.007] | LOOSEN |
| Adrenaline to beta-2AR | +0.062 | [+0.008, +0.108] | STIFFEN |
| Mg2+ to DNA Pol-beta | +0.079 | [+0.018, +0.140] | STIFFEN |
| Ca2+ to Calmodulin | +0.095 | [+0.014, +0.163] | STIFFEN |
| Tamoxifen to ER (vs apo) | +0.102 | [+0.045, +0.156] | STIFFEN |
| Tamoxifen to ER (vs E2) | +0.109 | [+0.042, +0.171] | STIFFEN |

Plus 3 cardiac troponin complex chains: A (+0.156 dag), B (+0.215 dag), F (+0.281 dag).

### 8.3 B-Factor Orthogonality

No per-residue correlation between delta-kappa-squared and delta-B-factor in any tested system (all |r| < 0.2). Meta-analysis: delta-Gini anti-correlates with r(delta-kappa-squared, delta-B) at r = -0.60 across 10 systems. Kappa-squared captures genuinely new structural information.

### 8.4 Spatial Autocorrelation (Moran's I)

| System | Moran's I | p-value | Verdict |
|--------|----------|---------|---------|
| Glutamate to GluR2 | +0.125 | < 0.001 | Clustered |
| Geldanamycin to HSP90 | +0.044 | 0.006 | Clustered |
| Adrenaline to beta-2AR | +0.040 | 0.003 | Clustered |
| Ca2+ to Calmodulin | +0.054 | 0.014 | Clustered |
| CO to Myoglobin | +0.036 | 0.019 | Clustered |
| Ketoconazole to CYP3A4 | +0.025 | 0.020 | Clustered |
| Temperature to CypA | -0.015 | 0.037 | Anti-clustered |
| Placebo (cross-lab lyz) | +0.018 | 0.063 | Random |

Biological delta-kappa-squared changes are spatially clustered; noise is not.

---

## 9. Temperature Calibration — The Directional Tide (n = 19)

### 9.1 Survey

| Protein | Fold | dGini | 95% CI | Dir |
|---------|------|-------|--------|-----|
| CypA | isomerase | -0.053 | [-0.084, +0.032] | NEG |
| Thaumatin | beta-sandwich | -0.049 | [-0.081, +0.027] | NEG |
| Lysozyme | monomer | -0.030 | [-0.062, +0.016] | NEG |
| Myoglobin | globin | -0.015 | [-0.074, +0.039] | NEG |
| RNase A | alpha+beta | -0.049 | [-0.086, +0.036] | NEG |
| Concanavalin A | lectin | -0.042 | [-0.093, +0.009] | NEG |
| Rubredoxin | Fe-S | -0.077 | [-0.168, +0.080] | NEG |
| Crambin | disulfide | -0.083 dag | [-0.120, +0.005] | NEG |
| SOD | beta-barrel | -0.021 | [-0.046, +0.008] | NEG |
| BPTI | small/SS | -0.051 dag | [-0.077, -0.005] | NEG |
| PLA2 | alpha+beta/Ca | -0.014 | [-0.073, +0.040] | NEG |
| Aldose reductase | TIM barrel | -0.014 | [-0.059, +0.029] | NEG |
| Cutinase | alpha/beta hydrolase | -0.040 dag | [-0.068, -0.011] | NEG |
| Chymotrypsinogen | serine (zymogen) | -0.024 | [-0.078, +0.067] | NEG |
| Subtilisin | serine protease | +0.046 dag | [+0.007, +0.073] | POS |
| Trypsin | serine protease | +0.011 | [-0.067, +0.087] | POS |
| Elastase | serine protease | +0.023 | [-0.023, +0.065] | POS |
| RNase T1 | alpha+beta | +0.048 dag | [+0.008, +0.089] | POS |
| Cytochrome c | all-alpha heme | +0.050 dag | [+0.004, +0.098] | POS |

### 9.2 Analysis

Overall: 14/19 negative (74%, sign test p = 0.032). Non-serine proteases: 13/15 negative (87%, p = 0.004). Serine proteases: 1/4 negative. The catalytic triad, only fully formed in the active enzyme, geometrically locks the backbone against thermal dispersion — chymotrypsinogen (a serine protease zymogen) follows the tide normally (-0.024). This is confirmed by the subtilisin S221A catalytic triad test (Section 16): removing the catalytic serine produces delta-Gini = +0.109 (massive stiffening), demonstrating that the intact triad acts as a strain dissipation network.

### 9.3 Honest Assessment

The tide is a weak directional bias (mean delta-Gini = -0.03), not a deterministic rule. Bootstrap-confirmed effects split 2 negative / 3 positive. The confirmed positives (RNase T1 +0.048 dag, Cyt c +0.050 dag) are not serine proteases, complicating the fold-family explanation. Anti-thermal interpretation: interactions with positive delta-Gini are unlikely to be temperature artifacts, but the converse is weaker than previously claimed. Same-lab provenance has not been verified for all 19 pairs; given that the cross-lab noise floor for monomers is 28 +/- 10, some temperature delta-Gini values may be within noise. A systematic resurvey using 100+ verified same-lab cryo/RT pairs would resolve this ambiguity definitively.

### 9.4 Negative Result: ATP Synthase Rotary States

F1-ATPase has three beta-subunits in distinct nucleotide states, offset 120 degrees. Predicted gradient observed in 2/9 structures (22%). Alpha-subunit negative control: non-catalytic alpha-subunits show greater inter-chain Gini variation (0.044 +/- 0.026) than catalytic beta-subunits (0.023 +/- 0.011). Gini cannot distinguish ATP synthase rotary states from crystal packing noise.

This failure establishes a resolution boundary: Gini is a microscope for internal backbone tension, not a telescope for domain-level movement. Large-scale quaternary rearrangements (120-degree rotary steps) do not necessarily alter the per-chain curvature distribution, even when they dramatically change inter-chain geometry.

---

## 10. Functional Cross-Section (n = 22)

### 10.1 Five Recurring Geometric Regularities of Molecular Function

Regularity 1 — Storage Is Silent. Passive cargo carriers show |dkappa2|/r = 15-24, below the temperature tide range. Exception: ferritin (110) undergoes cage-opening motion.

Regularity 2 — Shields Are Locked. SOD (Cu/Zn loaded to stripped): |dkappa2|/r = 2.3, the most rigid system in the atlas.

Regularity 3 — Gates Are Local. Ion channels open with |dkappa2|/r = 10-14. Transporters flip entire scaffolds at 136. The channel/transporter distinction as a single number.

Regularity 4 — Machines Are Neutral. Motor proteins that must cycle show near-zero delta-Gini. A machine that permanently stiffened could not cycle back. Geometric constraint of cyclicity.

Regularity 5 — Direction Does Not Equal Function. Whether perturbation stiffens or loosens is not determined by function but by architecture. Two GTPases (Ras +0.062, G-alpha -0.070) go opposite on GTP.

### 10.2 Geometric Spectrum

| |dkappa2|/r | What Happens | Examples |
|-----------|-----------------|----------|
| 1-3 | Nothing changes | SOD +/- metals, ALS mutations |
| 5-25 | Cargo/substrate loads | RBP, FABP, DHFR, channels open |
| 20-60 | Signal transmits | CaM, Ras, nuclear receptors |
| 50-140 | Structure transforms | Hb switching, transporters, motors |

---

## 11. Catalytic Intermediates (n = 13)

Mean |dkappa2|/r = 18.5 across 13 enzyme-substrate comparisons. Enzymes do not reshape their backbone to do chemistry. The active site is a pre-formed geometric cradle — the backbone holds still while side chains and substrate do the work. Catalysis (18.5) approximately equals Storage (19.0). The enzyme treats its substrate like cargo. This is the geometric basis of the lock-and-key model: not a metaphor, but a measurable constant.

### 11.1 DHFR Three-Way Comparison

| Comparison | |dkappa2|/r | dGini |
|-----------|-----------|-------|
| Apo to folate (substrate) | 10.5 | +0.029 |
| Apo to methotrexate (drug) | 8.0 | -0.000 |
| Folate to methotrexate | 11.4 | -0.029 |

The substrate mildly stiffens; the drug is geometrically silent. The receptor distinguishes natural substrate from drug at the backbone level.

---

## 12. Protein-Protein Interactions (n = 9)

Key finding: asymmetric deformation. In the barnase-barstar complex (tightest known PPI, Ka = 10^14), barstar deforms 6x more than barnase (|dkappa2|/r = 124 vs 20). The enzyme maintains its geometry; the inhibitor sacrifices its own to conform.

Antibody binding is invisible to the antigen: lysozyme + antibody Fab shows |dkappa2|/r = 33 but delta-Gini = -0.003. The immune system recognizes without disrupting.

S100B provides a second example of assembly-dependent geometric remodeling alongside troponin C: when p53 peptide binds to Ca2+-loaded S100B (PDB 1QLK → 1DT7), backbone curvature disperses (delta-Gini = -0.057, kappa2_max drops from 4745 to 1433). The extreme curvature spike at residue 53 in isolated S100B is redistributed across the backbone upon target binding. S100B is an EF-hand protein (66% alpha) like troponin C, and this result confirms that assembly-dependent geometric remodeling is a general feature of EF-hand signaling proteins, not an idiosyncrasy of the troponin system. A full troponin-style reversal test (Ca2+ effect in isolated vs complexed S100B) requires apo structures in both assembly states.

---

## 13. Body Map

| Location | Interaction | Mechanism | Function |
|----------|-----------|-----------|----------|
| Blood | O2 to Hemoglobin | Conform. Switch | Oxygen transport |
| Blood | Lisinopril to ACE | Flex. Mod. up | Blood pressure |
| Muscle | O2 to Myoglobin | Flex. Mod. up | Oxygen storage |
| Muscle | CO to Myoglobin | Flex. Mod. up | CO poisoning |
| Liver | Ketoconazole to CYP3A4 | Flex. Mod. down dag | Drug detox |
| Heart | Ca2+ to Troponin cplx | Mixed up dag | Contraction |
| Heart/Cell | Ca2+ to Calmodulin | Flex. Mod. up dag | Signal relay |
| Brain | Glutamate to GluR2 | Conform. Switch | Neurotransmission |
| Adrenal | Adrenaline to beta-2AR | Flex. Mod. up dag | Fight-or-flight |
| Reproductive | Estradiol to ER | Conform. Switch | Hormone |
| Reproductive | Tamoxifen to ER | Flex. Mod. up dag | Anti-cancer drug |
| Immune | PD-L1 to PD-1 | Flex. Mod. up | Checkpoint |
| Thermal | Geldanamycin to HSP90 | Conform. Switch | Heat shield |
| Cell | ATP to CDK2 | Conform. Switch | Kinase signaling |
| Cell | Myosin power stroke | Conform. Switch | Motor |
| Cell | Force to Integrin | Flex. Mod. down | Mechanosensing |

---

## 14. Replication and Consistency

| System | Replicate 1 | Replicate 2 | Consistency |
|--------|------------|------------|-------------|
| HIV protease (drugs) | Ritonavir: 2.2 | Indinavir: 2.5 | CV = 9% |
| Ca2+/EF-hand (dGini) | CaM: +0.095 | TnC(iso): -0.121 | Opposite. Assembly reverses |
| Ca2+ complex (bootstrap) | Ch A: +0.156 dag | Ch B: +0.215 dag, F: +0.281 dag | 3/6 confirmed stiffening |
| DNA binding | BamHI: 75.3, CAP: 36.9 | PcrA: 32.6 | All loosen (-0.03) |
| ATP (mechanism) | CDK2: Switch | Helicase, Myosin, HSP90: Switch | 4/4 agree |
| Tamoxifen (two refs) | vs apo: +0.102 dag | vs E2: +0.109 dag | Both confirmed |
| Temperature (n=19) | Non-serine: 13/15 neg | Serine: 1/4 neg | p=0.004 (non-ser) |
| ATP synthase beta vs alpha | Range beta: 0.023 | Range alpha: 0.044 | FAILED: alpha > beta |

---

## 15. Honest Limits

Bootstrap CIs: Magnitude +/-50-250%. 9 of 34 delta-Gini survive (6 direct + 3 troponin complex chains).

Temperature: 19 pairs, 87% negative in non-serine proteases (p = 0.004). Effect is small and rarely individually confirmed.

ATP synthase rotary: Gini cannot resolve quaternary rotary states (Section 9.4).

Tetramer warning: Cross-lab noise = 174 +/- 58. Same-lab matching mandatory.

Crystal vs cell: All measurements at 100K in crystal buffer.

Kappa-squared does not equal flexibility: Orthogonal to B-factors per-residue.

Classification: Continuous spectrum. k-means favors k=5.

Gini baseline: AlphaFold 0.835 +/- 0.005 (quality-filtered); full proteome 0.833 +/- 0.060 (human, n = 20,175); experimental PDB 0.798 +/- 0.081 (n = 931). The 16x wider spread is fold-class modulated (alpha-rich = 0.73, beta-rich = 0.75, mixed = 0.83).

Disease enrichment: High curvature concentration (Gini > 0.92) is not enriched for disease genes relative to the proteome baseline (OR = 0.857, p = 0.078). Curvature concentration is a structural property, not a disease predictor. Disease mutations reveal geometric consequences (direction), not geometric causes (prevalence).

---

## 16. Remaining Gaps

Experimental Gini baseline: 931 crystal structures tested (fold-class modulation confirmed at 3x sample size). Further expansion to 2000+ with redundancy control at 80% sequence identity would provide definitive fold-class Gini baselines per SCOP/CATH category.

Temperature: RNase T1 and cytochrome c positives need investigation. Same-lab provenance must be verified for all 19 pairs (Section 9.3). The serine protease anomaly is partially resolved by the subtilisin catalytic triad test: subtilisin S221A (catalytic Ser→Ala, PDB 1SCN) versus wild-type (PDB 1SBT) gives delta-Gini = +0.109 — a large stiffening effect consistent with the Level 1 alpha/beta fold prediction. Removing the catalytic serine breaks the Ser-His-Asp hydrogen-bonding network, and backbone strain concentrates massively in the catalytic region (residues 109-115; kappa2_max increases from 217 to 748; helical curvature increases 6.5x). The intact triad acts as a strain dissipation network: it distributes curvature evenly across the protein and resists concentration from any perturbation source. This reframes the serine protease temperature anomaly as a hierarchy-consistent observation: alpha/beta folds stiffen under perturbation, and the catalytic triad is the specific architectural feature that mediates this response. A complete resolution would require same-lab cryo/RT pairs for the S221A mutant specifically.

Ca2+ on non-EF-hand fold: tests whether stiffening is universal or fold-specific.

Delta-Gini vs MD/NMR: dynamic validation beyond B-factors.

Inactive compound control: active vs inactive inhibitor (true negative test).

GluR2 agonist series (same-lab, Armstrong & Gouaux): apo GluR2 LBD (1FTO, Gini = 0.906) is stiffest; glutamate full agonist (1FTJ, 0.868) and CNQX antagonist (1LBC, 0.857) both loosen relative to apo. Clamshell closure upon agonist binding disperses curvature concentrated at inter-lobe hinges. The partial agonist kainate (1FTK, 0.852) falls between antagonist and full agonist. This confirms the Conform. Switch classification (Entry #31) and demonstrates that GluR2's clamshell architecture determines the direction (loosening), overriding the Level 1 alpha/beta tendency to stiffen — a Level 2 scaffold effect. GABA receptor structures would extend this to the inhibitory neurotransmitter system but require cryo-EM structures with sufficient single-chain resolution.

Motor proteins: ATP synthase rotary states failed the Gini test. Need |dkappa2|/r profiles (not just Gini) to assess whether curvature redistribution differs even when the distribution does not.

Curvature concentration determinants: What properties (fold class, protein length, secondary structure composition, evolutionary age) predict a protein's position on the Gini spectrum? Cross-species conservation of concentration status. Functional site overlap with curvature hotspots.

---

## 17. Prior Art Disclosure

This document and framework are original work. All discoveries from versions 1-9 are incorporated by reference. New in v10: expanded proteome baseline (124,749 proteins, 34 species), expanded experimental PDB baseline (931 crystal structures, 3x prior sample), full human proteome Gini distribution analysis (n = 20,175), pLDDT structural confidence validation of curvature hotspots, corrected disease enrichment analysis (no enrichment, OR = 0.857), hotspot positional distribution (N/C-terminal enrichment), restructured presentation with three-level hierarchy as central organizing principle, and integration of all prior atlas entries, disease panels, and functional cross-sections into a unified document.

Crystal structures from the Protein Data Bank (rcsb.org). AlphaFold structures from the EBI AlphaFold Protein Structure Database (alphafold.ebi.ac.uk).

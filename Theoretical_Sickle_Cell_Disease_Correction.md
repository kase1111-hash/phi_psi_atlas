OPEN SOURCE PRIOR ART DISCLOSURE
Curvature Fix Maps for Structure-Based Drug Design:
Theoretical Candidates for Sickle Cell Disease Geometric Correction
Kase Branham
February 28, 2026
License: Public Domain (CC0 1.0 Universal)
1. Purpose of This Disclosure
This document constitutes a public disclosure of theoretical drug candidates, computational methods, and screening frameworks for the geometric correction of sickle cell hemoglobin. By publishing these concepts with full technical detail, this disclosure establishes prior art that prevents any party from obtaining patent protection on the described methods, candidates, or frameworks. All content is released into the public domain under CC0 1.0 Universal license.
2. Technical Background
2.1 The curvature fix map concept
We have developed a method for computing the complete residue-by-residue dihedral angle corrections needed to restore a disease-mutant protein to its wild-type backbone geometry. The method uses torus curvature (κ²) computed from the discrete Frenet frame on the Ramachandran torus T², comparing crystal structures of wild-type and mutant proteins to generate a “damage map” of curvature perturbations. The negative of this damage map defines the “fix map”—the specific dihedral corrections at each residue that would restore wild-type geometry.
2.2 Application to sickle cell disease
Sickle cell disease is caused by a single point mutation (E6V) in hemoglobin β-globin. Comparison of crystal structures 2DN2 (wild-type T-state, PDB) and 2HBS (sickle T-state, PDB) reveals that this single mutation causes distributed curvature perturbations across 85 of 146 residues. The three largest perturbations occur at residues 128 (Δκ² = +1859), 27 (Δκ² = −1128), and 136 (Δκ² = −1054)—far from the mutation site at position 6. The average dihedral correction needed is 4.1° in φ and 4.6° in ψ per perturbed residue.
2.3 Existing treatment baseline
Voxelotor (Oxbryta), the only FDA-approved curvature-modifying drug for sickle cell, achieves a Pearson correlation of r = −0.12 between its curvature effect profile and the theoretical fix map. This indicates weak geometric correction—the drug modifies hemoglobin geometry in the right general direction but with poor targeting precision.
3. Pharmacophore Definition
The fix map identifies three spatially distinct clusters of high-priority correction targets:
Cluster	Residues	Span (Å)	Required Action	Druggability
A	106, 111, 113, 115	10.3	Stiffen 111; relax 106, 113, 115	Small molecule
B	128, 136, 138	15.4	Stiffen 128; relax 136, 138	Macrocycle
C	7, 27	22.0	Stiffen 7; relax 27	Small polar molecule

3D coordinates for all pharmacophore points are derived from PDB structure 2HBS chain B (Cα positions). The complete coordinate set and pairwise distance matrix are provided in the supplementary data.
4. Theoretical Drug Candidates
4.1 TN-A1: Bifunctional helix G corrector
Scaffold: Pyridine-carboxamide with flexible ethylamine arms. Molecular weight: approximately 250–300 Da. The rigid pyridine core provides steric constraint at residue 111 (stiffening), while flexible amine arms extending toward residues 106, 113, and 115 allow backbone freedom (relaxation). Amide groups form hydrogen bonds to lock backbone geometry at the stiffening site. Structural analogy: imatinib core scaffold (pyridine-carboxamide kinase inhibitor).
4.2 TN-B1: Macrocyclic helix H corrector
Scaffold: 16-membered macrolactam with pendant hydrogen-bond donor. Molecular weight: approximately 400–500 Da. An aromatic group with primary amine at one pole provides rigid H-bond donation at residue 128 (stiffening). A flexible aliphatic chain at the opposite pole allows backbone relaxation at residues 136 and 138. Ring diameter of approximately 15 Å matches the cluster B span. Structural analogy: erythromycin (macrolide antibiotic) and lorlatinib (macrocyclic kinase inhibitor).
4.3 TN-C1: N-terminal backbone competitor
Scaffold: Salicylate-acetamide derivative. Molecular weight: approximately 195 Da. A carboxylate group competes for backbone hydrogen bonds at residue 27 (relaxation). A phenol provides secondary H-bond contact, and an acetamide group mimics backbone geometry as a φ/ψ competitor. Structural analogy: aspirin, diflunisal, 5-hydroxymethylfurfural (5-HMF, which is already known to modify hemoglobin at the N-terminal valine).
4.4 TN-SCD-1: Three-component cocktail
A combination therapy consisting of TN-A1 + TN-B1 + TN-C1, administered together to achieve multi-cluster geometric correction. Predicted combined correlation with the fix map: r ≈ −0.25 to −0.35, representing a 2–3× improvement over voxelotor monotherapy (r = −0.12). Each component targets a distinct spatial cluster, and each can be independently dosed and optimized.
5. Drug Screening Method
We disclose a general method for scoring drug candidates by their geometric correction profile:
Step 1: Obtain crystal structures of the wild-type and disease-mutant protein.
Step 2: Compute the torus curvature fix map (Δκ² per residue).
Step 3: For each candidate drug, compute the drug-induced curvature change profile (via molecular dynamics simulation of the drug-protein complex).
Step 4: Score the drug by Pearson correlation between its curvature effect profile and the fix map. A score of r = −1.0 represents perfect geometric correction; r = 0 represents no correction; r = +1.0 represents perfect amplification of the disease effect.
Step 5: Assess off-target risk using the two-axis vulnerability model (Section 6).
6. Proteome Vulnerability Model
We disclose a two-axis model for predicting off-target drug effects from curvature data alone:
Axis 1: Exposure. Structural similarity to the target protein (hemoglobin), computed as a weighted distance in curvature feature space (helix fraction, Σκ²/residue, Gini coefficient, chain length). Higher exposure indicates greater likelihood of drug interaction.
Axis 2: Fragility. Mechanical sensitivity to curvature perturbation, computed as the Gini coefficient divided by the logarithm of the total curvature budget. Higher fragility indicates that small perturbations could disrupt native geometry.
In the human proteome (n = 13,428 high-quality proteins), we find zero overlap between the top 5% most exposed and top 5% most fragile proteins. This indicates a wide therapeutic window for hemoglobin-targeting curvature correctors: the proteins most likely to interact with the drug are mechanically robust, and the proteins most vulnerable to perturbation are structurally dissimilar to hemoglobin.
7. Repurposing Candidates
The following existing FDA-approved or investigational drugs are identified as priority candidates for curvature-based screening against the hemoglobin fix map:
Hemoglobin-interacting: Voxelotor (baseline), efaproxiral (RSR13), 5-hydroxymethylfurfural derivatives, nitroglycerin.
Helix bundle modulators: Trifluoperazine, W-7, calmidazolium (calmodulin inhibitors); dantrolene (ryanodine receptor); ivacaftor (CFTR corrector—already FDA-approved as a helix bundle geometry corrector for cystic fibrosis).
Macrocyclic scaffolds: Lorlatinib, erythromycin, cyclosporine, vancomycin (as structural templates for cluster B targeting).
Backbone competitors: Aspirin, diflunisal, probenecid (as structural templates for cluster C targeting).
8. Generalization to Other Diseases
The curvature fix map method generalizes to any disease caused by a structural mutation with known wild-type and mutant crystal structures. Priority targets for future fix maps include: cystic fibrosis (CFTR ΔF508), Alzheimer’s disease (Aβ misfolding variants), Parkinson’s disease (α-synuclein aggregation variants), Marfan syndrome (fibrillin-1 mutations), and familial hypertrophic cardiomyopathy (myosin heavy chain mutations).
9. Legal Notice
This document is a public disclosure of intellectual content intended to establish prior art under 35 U.S.C. § 102. All methods, candidates, frameworks, and data described herein are released into the public domain under CC0 1.0 Universal (CC0 1.0) Public Domain Dedication. No rights are reserved. Any party may use, reproduce, modify, distribute, or build upon this work for any purpose without restriction or attribution requirement. This disclosure is intended to prevent patenting of the described concepts by any party.

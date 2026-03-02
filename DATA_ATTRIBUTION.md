# Data Attribution

This project uses protein structure data from the following external sources.
No CIF or PDB coordinate files are distributed with this repository; they are
downloaded at runtime by the scripts listed below. This document credits the
databases, depositors, and publications whose data this project relies on.

---

## 1. AlphaFold Protein Structure Database

**Provider:** DeepMind / EMBL-EBI
**URL:** https://alphafold.ebi.ac.uk
**License:** CC-BY 4.0
**Citation:**

> Jumper, J. et al. "Highly accurate protein structure prediction with
> AlphaFold." *Nature* 596, 583–589 (2021).
> https://doi.org/10.1038/s41586-021-03819-2
>
> Varadi, M. et al. "AlphaFold Protein Structure Database: massively expanding
> the structural coverage of protein-sequence space with high-accuracy models."
> *Nucleic Acids Research* 50, D439–D444 (2022).
> https://doi.org/10.1093/nar/gkab1061

**File format:** mmCIF (`.cif`)
**Naming convention:** `AF-{UniProtID}-F1-model_v{N}.cif`

### Structures used

The following AlphaFold predicted structures are downloaded and used across
several scripts in this project:

#### Core analysis dataset (`cornu_spirals_on_T2.py`)

Downloaded as PDB from RCSB (see Section 2), but AlphaFold structures for
these UniProt IDs are also fetched when running in `alphafold` mode.

#### DALI validation pairs (`dali_batch.py`, `cif_to_pdb.py`, `download_pdb_pairs.py`)

| Pair | Query (UniProt) | Gene | Hit (UniProt) | Gene |
|------|-----------------|------|---------------|------|
| 1 | P00918 | CA2 | P59991 | KRTAP12-2 |
| 2 | P60174 | TPI | P35270 | SPR |
| 3 | P62258 | 14-3-3e | O43808 | SLC25A17 |
| 4 | P11166 | GLUT1 | P09914 | IFIT1 |
| 5 | P02751 | FN1 | P21333 | FLNA |
| 6 | P01857 | IgG1 | P32969 | RPL9 |
| 7 | P38646 | GRP75 | O15344 | MID1 |
| 8 | P07900 | HSP90a | O14829 | PPEF1 |
| 9 | P04626 | HER2 | O95302 | FKBP9 |
| 10 | P21860 | ErbB3 | P28799 | GRN |

#### True-positive validation pairs (`download_tp_pairs.py`)

| Pair | Query (UniProt) | Gene | Hit (UniProt) | Gene |
|------|-----------------|------|---------------|------|
| 11 | P00338 | LDH-A | P07195 | LDH-B |
| 12 | P68871 | Hb-beta | P69891 | Hb-G1 |
| 13 | P00338 | LDH-A | P40926 | MDH2 |
| 14 | P62258 | 14-3-3e | P31946 | YWHAB |
| 15 | P69905 | Hb-alpha | P02144 | Myoglobin |
| 16 | P00918 | CA2 | P23280 | CA6 |

#### Benchmark query panel (`fetch_benchmark_proteins.py`)

| UniProt ID | Description |
|------------|-------------|
| P69905 | Hemoglobin alpha |
| P68871 | Hemoglobin beta |
| P02144 | Myoglobin |
| P00918 | Carbonic anhydrase 2 |
| P00533 | EGFR |
| P04626 | ErbB2/HER2 |
| P21860 | ErbB3 |
| P01857 | IgG1 Fc |
| P01834 | Ig kappa chain C |
| P02751 | Fibronectin |
| P00338 | L-lactate dehydrogenase A |
| P04406 | GAPDH |
| P60174 | Triosephosphate isomerase |
| P01308 | Insulin |
| P62258 | 14-3-3 epsilon |
| P11166 | GLUT1 |
| P06748 | Nucleophosmin |
| P10636 | Tau |
| P07900 | HSP90-alpha |
| P38646 | Mortalin/GRP75 |

#### Proteome-scale surveys (`proteome_curvature_survey.py`, `cross_species_survey.py`, `phyllotaxis_protein.py`, `helix_geometry.py`)

These scripts process bulk AlphaFold proteome downloads (e.g. the full human
proteome). Users supply their own AlphaFold cache directory. All structures
carry the CC-BY 4.0 license from the AlphaFold Protein Structure Database.

---

## 2. RCSB Protein Data Bank (PDB)

**Provider:** Research Collaboratory for Structural Bioinformatics
**URL:** https://www.rcsb.org
**License:** CC0 1.0 (PDB data are in the public domain)
**Citation:**

> Berman, H.M. et al. "The Protein Data Bank." *Nucleic Acids Research* 28,
> 235–242 (2000). https://doi.org/10.1093/nar/28.1.235
>
> Burley, S.K. et al. "RCSB Protein Data Bank (RCSB.org): delivery of
> experimentally-determined PDB structures alongside one million computed
> structure models of proteins from artificial intelligence/machine learning."
> *Nucleic Acids Research* 51, D483–D489 (2023).
> https://doi.org/10.1093/nar/gkac1052

**File formats:** PDB (`.pdb`) and mmCIF (`.cif`)

### Core folding-kinetics dataset (`cornu_spirals_on_T2.py`)

These experimentally determined structures are downloaded from RCSB PDB:

| PDB ID | Protein | Chain | Reference context |
|--------|---------|-------|-------------------|
| 2CI2 | CI2 | I | Jackson & Fersht 1991 |
| 1SRL | SH3-src | A | Grantcharova & Baker 1997 |
| 2PTL | Protein-L | A | Scalley et al. 1997 |
| 1CSP | CspB | A | Schindler et al. 1995 |
| 1UBQ | Ubiquitin | A | Khorasanizadeh et al. 1993 |
| 1LMB | Lambda repressor | 3 | Burton et al. 1997 |
| 1VII | Villin HP | A | Kubelka et al. 2003 |
| 1ENH | Engrailed HD | A | Mayor et al. 2003 |
| 1PGA | Protein G | A | Park et al. 1999 |
| 1FKB | FKBP | A | Main et al. 1999 |
| 1BNI | Barnase | A | Matouschek et al. 1990 |
| 1BTA | Barstar | A | Schreiber & Fersht 1993 |
| 3CHY | CheY | A | Lopez-Hernandez & Serrano 1996 |
| 2ACY | Acylphosphatase | A | Chiti et al. 1999 |
| 1HRC | Cytochrome c | A | Mines et al. 1996 |
| 7RSA | RNase A | A | Houry et al. 1999 |
| 2LZM | Lysozyme | A | Khorasanizadeh et al. 1993 |
| 1RX2 | DHFR | A | Ionescu et al. 2000 |
| 2TRX | Thioredoxin | A | Georgescu et al. 1998 |
| 1SHG | SH3-spectrin | A | Viguera et al. 1994 |
| 1TEN | TNfn3 | A | Clarke et al. 1997 |
| 1ARR | Arc repressor | A | Milla & Sauer 1994 |
| 1MBD | Myoglobin | A | Jennings & Wright 1993 |

### Coordinate-based validation (`validation_from_coordinates.py`)

| PDB ID | Protein |
|--------|---------|
| 1MBN | Myoglobin (deoxy) |
| 1A6M | Myoglobin (met) |
| 1MBD | Myoglobin (CO) |
| 1LYZ | Lysozyme (hen) |
| 1HEL | Lysozyme (human) |
| 1UBQ | Ubiquitin |
| 1UBI | Ubiquitin (2) |

### Allosteric analysis (`allosteric_barcodes.py`)

| PDB ID | State | Resolution | Depositor |
|--------|-------|------------|-----------|
| 2DN2 | Human hemoglobin T-state (deoxy) | 1.25 A | Park et al. 2006 |
| 2DN1 | Human hemoglobin R-state (oxy) | 1.25 A | Park et al. 2006 |

### Reviewer validation — Test A: 500+ diverse structures (`reviewer_validation.py`)

A curated set of 530 PDB entries spanning alpha-rich, beta-rich, alpha/beta,
large, small, enzyme, and binding-protein fold classes. Selected from CATH
representatives, the PDB Top 500, and manual curation. Full list of PDB IDs:

1MBN, 1HHO, 1A6N, 1CLL, 1CFD, 2CCL, 1GKP, 1MBA, 1UTG, 1YCC, 1CYC, 1BBH,
1HRC, 1AKR, 1ALK, 256B, 1CPQ, 1ECM, 1FLP, 1GPR, 1HEL, 1LGJ, 1MHR, 1PBN,
1RVS, 1VLS, 3INS, 1YPA, 1ASZ, 2MHR, 1ECA, 1ISA, 1HBG, 2DHB, 1FLG, 1MOF,
1BAB, 2IFB, 2LHB, 1ITH, 1CPE, 1FHA, 1AHL, 4CPA, 1BJ7, 1BGC, 1DPS, 1A62,
1HCN, 1NKR, 1REI, 1FKN, 1PGB, 1IGD, 1TEN, 1FNA, 1TIT, 1WIT, 2IG2, 1CFB,
1AXJ, 1CQE, 1DAB, 1EKP, 1FLT, 1GCN, 1HFI, 1JPC, 1KAP, 1LCD, 1LCL, 1MAI,
1MUP, 1NLS, 1OAA, 1PLC, 1QCF, 1RAF, 1RKI, 1SEM, 1TEA, 1UBI, 2CDV, 2SNS,
3CHY, 2AIT, 1PLQ, 1CSK, 1BAM, 1BPI, 1QJP, 1BEO, 7RSA, 1AUA, 1PMC, 1NPC,
1QDD, 1AOE, 1HKA, 1GPB, 1LYZ, 1AKE, 1HSB, 1PHT, 1TPH, 1TIM, 3DFR, 1DRF,
2ACE, 1CHD, 1CTS, 1DOR, 1ENH, 1FRD, 1GCA, 1HOE, 1LDM, 1MJC, 1NXB, 1PHP,
1RNS, 1SBP, 1THV, 1XEL, 2CBA, 2GST, 2PTC, 3ADK, 4ENL, 5RXN, 1PPE, 1OCA,
1EZM, 2IHL, 1ALF, 1AOF, 1AYL, 1BGE, 1BIF, 1BMD, 1CDG, 1CRN, 1CTJ, 1DPE,
1DUB, 1EMD, 1EOK, 1FCB, 1FIP, 1GKY, 1SMD, 1TQN, 1O8A, 1W0F, 1FMV, 1VOM,
1EVE, 1EA5, 1GOS, 1HXW, 1PKD, 1AO0, 1B41, 1CDL, 1CLC, 1DL2, 1E2X, 1F13,
1GTM, 1HDC, 1IEP, 1JR1, 1KWF, 1LFO, 1MML, 1NAR, 1OHO, 1PMR, 1QPQ, 1RCO,
1RTQ, 1SFP, 1TCA, 1URF, 2AAA, 2RN2, 3EST, 4CTS, 1ACB, 1AOZ, 1AUO, 1B0Y,
1B3A, 1BKR, 1BQC, 1BTE, 1BXO, 1CEX, 1CMB, 1CVO, 1HHP, 1PTF, 1SH1, 2HPR,
1ABE, 1AHO, 1FKB, 1MCI, 1ROP, 1UBQ, 1VII, 2ABD, 2RHE, 3ICB, 1EDN, 1ERY,
1FAS, 1FWP, 1GPE, 1HMV, 1JLI, 1L58, 1MBG, 1N0W, 1PEN, 1Q2S, 1SPH, 2CI2,
1RGS, 1AMP, 1APR, 1CTF, 1DIG, 1EGF, 1FDD, 1FNF, 1GCI, 1H4X, 1IDC, 1KIV,
1LIL, 1MEE, 1MIN, 1NUC, 1ONC, 1PGA, 1QGV, 1RIS, 1STN, 1ERK, 4EK3, 1HCK,
1CDK, 1ATP, 1JBP, 2PTH, 1CSE, 2EST, 1PPB, 1SGT, 1TRY, 3TEC, 1ACH, 1GEN,
1R1R, 2DRI, 3RUB, 1AGT, 1BAN, 1BDM, 1BVU, 1CAH, 1CEL, 1CLE, 1COX, 1DDT,
1DRE, 1EBH, 1FAX, 1GHR, 1GLM, 1HEW, 1HUW, 1JAG, 1KFN, 1LAM, 1M6T, 1NAG,
1OVA, 1PDO, 1POC, 1QBI, 1RBO, 1SHK, 1THG, 1XYZ, 2ARC, 2CMD, 2TMA, 1CGP,
2CGP, 1A52, 1ERE, 3ERT, 1FTJ, 1FTK, 1YET, 1YER, 1O86, 3PJR, 1FKG, 1MFB,
1RTB, 1SKH, 1WHI, 1XIA, 2BBK, 2MSB, 3MBP, 1AAC, 1ACI, 1ARB, 1BAK, 1CAA,
1DIL, 1DOT, 1EBZ, 1FDH, 1GAP, 1HFA, 1IOB, 1JBE, 1KRN, 1LIB, 1MCT, 1NGR,
1OSA, 1PDA, 1QNF, 1RPO, 1STP, 1TML, 2ACQ, 2CAB, 2FBJ, 2LIV, 1THW, 2VHK,
3K0M, 3K0N, 193L, 1BWH, 1BZR, 1BZP, 1A6M, 1A6G, 1J1E, 1J1D, 3TCT, 2TN4,
5E6E, 5E83, 4HHB, 2HTQ, 3RRQ, 4TMN, 8TLN, 5IKR, 1HSG, 1TTC, 1RKP, 1BPD,
1BHM, 3ZDX, 4ZQK, 4TLT

### Reviewer validation — Test B: cryo/RT temperature pairs (`reviewer_validation.py`)

| Label | Cryo PDB | RT PDB | Chain | Source publication |
|-------|----------|--------|-------|--------------------|
| CypA | 3K0N | 3K0M | A | Fraser 2009 |
| Thaumatin | 1THW | 2VHK | A | Ko 2003 |
| Myoglobin | 1BZR | 1BZP | A | Ostermann 2000 |
| Lysozyme | 193L | 1BWH | A | Walsh 1998 |
| RNase A | 1AFU | 1RAT | A | Rasmussen 1992 |
| Concanavalin A | 1NLS | 1GKB | A | Parkin 1996 |
| Subtilisin | 1SCJ | 1SCA | A | Davail 1994 |
| Trypsin | 1TRY | 1TPO | A | Walter 1982 |
| Insulin | 3INS | 4INS | A | Smith 1984 |
| Ribonuclease S | 1RCA | 1RNS | A | Tilton 1992 |
| Elastase | 3EST | 1EST | A | Meyer 1988 |
| Chymotrypsin | 1ACB | 4CHA | A | Tsukada 1977 |
| Rubredoxin | 1BRF | 5RXN | A | Day 1992 |
| Crambin | 1CRN | 1CBN | A | Jelsch 2000 |
| Superoxide dismutase | 1SXA | 1SXC | A | Carugo 1999 |

### Reviewer validation — Test C: ATP synthase rotary states (`reviewer_validation.py`)

| Label | PDB ID | Chains | Depositor/Publication |
|-------|--------|--------|-----------------------|
| F1-ATPase ground state | 1BMF | D, E, F | Abrahams et al. 1994 |
| F1-ATPase transition (AlF4) | 1E1R | D, E, F | Braig et al. 2000 |
| F1-ATPase AMPPNP | 1H8E | D, E, F | Menz et al. 2001 |

### Extended validation — additional temperature pairs (`extended_validation.py`)

| Label | Cryo PDB | RT PDB | Chain | Fold class |
|-------|----------|--------|-------|------------|
| Ribonuclease T1 | 9RNT | 1RGE | A | alpha+beta |
| Cytochrome c | 1CRC | 1AKK | A | all-alpha |
| BPTI | 4PTI | 5PTI | A | small/disulfide |
| Thermolysin | 8TLN | 4TMN | A | metalloprotease |
| Phospholipase A2 | 1POA | 1PSJ | A | alpha+beta/Ca2+ |
| Chymotrypsinogen | 1CHG | 2CGA | A | serine protease |
| Aldose reductase | 1US0 | 2ACR | A | TIM barrel |
| Cutinase | 1CEX | 1CUS | A | alpha/beta hydrolase |

### Extended validation — additional ATP synthase structures (`extended_validation.py`)

| Label | PDB ID | Chains | Source |
|-------|--------|--------|--------|
| Ground state | 1BMF | D, E, F (beta) + A, B, C (alpha) | Abrahams 1994 |
| Transition (AlF4) | 1E1R | D, E, F + A, B, C | Braig 2000 |
| AMPPNP bound | 1H8E | D, E, F + A, B, C | Menz 2001 |
| ADP crystal form 2 | 1E79 | D, E, F + A, B, C | Braig 2000 |
| TNP-AMP bound | 1NBM | D, E, F + A, B, C | Menz 2001 |
| Mixed occupancy | 1OHH | D, E, F + A, B, C | Kagawa 2003 |
| Pre-hydrolysis | 2JDI | D, E, F + A, B, C | Bowler 2006 |
| Post-hydrolysis | 2JJ2 | D, E, F + A, B, C | Bowler 2007 |
| Yeast F1 | 2V7Q | D, E, F + A, B, C | Kabaleeswaran 2006 |

### High-resolution validation set (`pdb_validation.py`)

A curated set of ~200 high-resolution (< 1.5 A) X-ray crystal structures used
for independent curvature validation. Sourced from the PDB Top 200 lists, the
PDB high-resolution hall of fame, and the PDB teaching set. The full list is
defined in `pdb_validation.py:cmd_download()`.

### Drug similarity reference (`drug_similarity_search.py`)

| PDB ID | Usage |
|--------|-------|
| 2DN1 | Beta-globin R-state reference for drug-curvature analysis |

---

## 3. DALI Structure Comparison Server

**Provider:** Holm group, University of Helsinki
**URL:** http://ekhidna2.biocenter.helsinki.fi/dali/
**Citation:**

> Holm, L. "DALI and the persistence of protein shape." *Protein Science* 29,
> 128–140 (2020). https://doi.org/10.1002/pro.3749

The DALI server is used for independent structural superposition validation
of geometry-only hits (`dali_batch.py`, `dali_validation.py`,
`dali_pair_guide.py`). Structures are submitted to the DALI web server and
Z-scores are recorded manually.

---

## 4. UniProt

**Provider:** UniProt Consortium
**URL:** https://www.uniprot.org
**License:** CC-BY 4.0
**Citation:**

> UniProt Consortium. "UniProt: the Universal Protein Knowledgebase in 2023."
> *Nucleic Acids Research* 51, D523–D531 (2023).
> https://doi.org/10.1093/nar/gkac1052

UniProt identifiers are used throughout this project to identify proteins.
Functional annotations from UniProt are cached in `data/uniprot_features/`
for peak-decode and inverse-fold analyses.

---

## 5. CATH Protein Structure Classification

**Provider:** Orengo group, University College London
**URL:** https://www.cathdb.info
**Citation:**

> Sillitoe, I. et al. "CATH: increased structural coverage of functional
> space." *Nucleic Acids Research* 49, D266–D273 (2021).
> https://doi.org/10.1093/nar/gkaa1079

CATH domain classifications are used for fold-diversity benchmarking in
`cath_conservation_test.py` and `cath_fast.py`.

---

## 6. Embedded Code Attribution

### Pure-Python DSSP implementation (`cornu_spirals_on_T2.py`)

The DSSP secondary structure assignment code (Section 1 of
`cornu_spirals_on_T2.py`) is adapted from:

- **MDAnalysis** pydssp_numpy.py (MIT License)
- Which is itself a reimplementation of **PyDSSP** v0.9.0 by Shintaro Minami
  (MIT License)
- Based on the Kabsch-Sander hydrogen bond energy criterion from the original
  DSSP algorithm by Kabsch & Sander 1983.

---

## Summary

| Source | License | Data type | How obtained |
|--------|---------|-----------|--------------|
| AlphaFold DB | CC-BY 4.0 | Predicted structures (mmCIF) | Downloaded at runtime |
| RCSB PDB | CC0 1.0 | Experimental structures (PDB/mmCIF) | Downloaded at runtime |
| DALI Server | Academic use | Z-scores | Web submission |
| UniProt | CC-BY 4.0 | Protein identifiers & annotations | API queries at runtime |
| CATH | Academic use | Domain classifications | Downloaded at runtime |

No structure coordinate files are committed to this repository. All data is
fetched on demand from the above public databases.

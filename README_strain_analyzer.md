# Backbone Strain Analyzer

**Interactive walkthrough of torsional curvature analysis on any protein pair.**

Computes discrete torsional curvature (κ²), Gini coefficient, and mechanical mode classification from two crystal structures. Based on the Backbone Strain Atlas v15 (Knochenhauer 2026).

---

## Requirements

- Python 3.8+
- numpy
- scipy

```bash
pip install numpy scipy
```

No other dependencies. No internet connection needed. No force fields or simulations.

---

## Quick Start

```bash
# Download two structures from RCSB (mmCIF format)
# Example: HIV-PR apo (1HHP) vs indinavir-bound (1HSG)

python strain_analyzer.py 1HHP.cif 1HSG.cif A
```

That's it. The script walks you through 8 steps interactively.

---

## Usage

```
python strain_analyzer.py <apo.cif> <drug.cif> <chain> [--pocket-ref <ref.cif>]
```

| Argument | Description |
|----------|-------------|
| `apo.cif` | Reference structure (apo or first state) in mmCIF format |
| `drug.cif` | Perturbed structure (drug-bound or second state) in mmCIF format |
| `chain` | Chain identifier (e.g., `A`, `B`, `E`) |
| `--pocket-ref` | Optional: define binding pocket from a different structure |

### Getting structures

1. Go to [rcsb.org](https://www.rcsb.org)
2. Search for your PDB ID (e.g., `1HSG`)
3. Download → mmCIF format
4. Save as `1HSG.cif`

---

## The 8 Steps

Each step prints results and offers **[A/B/C/D]** options for deeper exploration. Press **Enter** to skip options and continue.

| Step | What it does | Key output |
|------|-------------|------------|
| **1. Load** | Parse mmCIF, extract backbone + ligands | Residue count, ligand ID |
| **2. Torsion angles** | Compute φ/ψ per residue | Secondary structure composition |
| **3. κ² (curvature)** | Frenet-Serret normalized curvature | Per-residue bending energy proxy |
| **4. Gini coefficient** | Inequality of κ² distribution | Single number per structure |
| **5. ΔGini** | Direction (stiffen/loosen) + magnitude | The mechanical fingerprint |
| **6. Hotspots** | Top 20 most-perturbed residues | Ranked by |Δκ²| |
| **7. Local vs distal** | Proximity to ligand | % of hotspots near binding site |
| **8. Mechanical mode** | Uniform / concentrator / exporter | How strain redistributes |

---

## Interpreting Results

### ΔGini (Step 5)

| Value | Meaning |
|-------|---------|
| **> 0** | Drug **stiffened** the backbone (concentrated bending energy) |
| **< 0** | Drug **loosened** the backbone (distributed bending energy) |
| **|ΔGini| < 0.012** | Below noise floor — not significant |
| **|ΔGini| > 0.03** | Substantial mechanical effect |

### Mechanical Mode (Step 8)

| Mode | Local ΔGini | Distal ΔGini | What it means |
|------|------------|-------------|---------------|
| **Uniform** | + or − | Same direction | Whole protein shifts together |
| **Concentrator** | + | − | Drug clamps pocket, relaxes periphery |
| **Exporter** | − | + | Drug relaxes pocket, stiffens periphery |

### Known reference values

| Target | Typical ΔGini | Direction |
|--------|--------------|-----------|
| CA-II | +0.004 to +0.112 | Always stiffens |
| HIV-PR | +0.003 to +0.151 | 92% stiffen |
| CDK2 | −0.008 to −0.020 | 71% loosen |
| HSP90 | −0.039 to +0.005 | 71% loosen |

---

## Caveats

- **Sphere size**: If >30% of residues are within 10 Å of the ligand (small proteins, large cofactors), the local/distal classification loses power. The script warns you.
- **Crystal form**: ΔGini includes crystallographic noise. For definitive results, compare multiple drug-bound structures to the same apo reference.
- **Not a potency predictor**: ΔGini does not correlate with binding affinity (ρ = +0.15, not significant). It measures mechanical redistribution, not thermodynamics.
- **Direction requires a panel**: A single comparison tells you direction for that drug. To determine the target's consensus direction, you need ≥15 comparisons.

---

## Examples

```bash
# CA-II: acetazolamide binding
python strain_analyzer.py 2CBA.cif 1YDA.cif A

# TTR: tafamidis (concentrator mode)
python strain_analyzer.py 1TTC.cif 3TCT.cif A

# CDK2: indirubin (loosening)
python strain_analyzer.py 1HCK.cif 1KF7.cif A

# Hemoglobin: T-state to R-state (induced-fit)
python strain_analyzer.py 4HHB.cif 1HHO.cif A

# ERK: phosphorylation (not a drug — post-translational modification)
python strain_analyzer.py 1ERK.cif 2ERK.cif A
```

---

## Non-interactive mode

For scripting, pipe empty lines to skip all options:

```bash
yes '' | python strain_analyzer.py 2CBA.cif 1YDA.cif A
```

---

## What this measures

κ² approximates the geometric component of **backbone bending energy density** in the Kirchhoff elastic rod framework, computed in conformational space (φ, ψ) rather than Cartesian space (x, y, z). The Gini coefficient measures how unequally that energy is distributed. ΔGini measures whether a perturbation concentrated or spread the bending energy.

For the full theoretical foundation, artifact controls, and 39-comparison validation dataset, see **The Backbone Strain Atlas v15**.

---

*Knochenhauer 2026. Method validated across 11 targets, 39 drug comparisons, 6 apo-apo controls, 10 energy-weighting tests, 6 artifact controls.*

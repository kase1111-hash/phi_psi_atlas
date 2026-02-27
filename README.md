# PhiPsi_Atlas

**Differential geometry of protein backbones on the Ramachandran torus**

PhiPsi_Atlas treats the (φ, ψ) dihedral trajectory of every protein as a curve on the flat torus T² = S¹ × S¹ and extracts geometric invariants — curvature classifications, winding numbers, and topological barcodes — that characterize backbone structure without 3D coordinates.

---

## What it does

Every protein backbone traces a path through Ramachandran space. PhiPsi_Atlas computes:

- **Torus curvature κ(s)** — geodesic curvature of the (φ,ψ) curve on T²
- **10-model AICc classification** — each secondary structure segment classified as geodesic, circular arc, clothoid, quadratic, Gaussian peak, sinusoidal, damped oscillation, exponential, sigmoid, or step
- **Winding numbers (p, q)** — topological winding of the backbone around each torus cycle
- **Geometric barcodes** — full per-segment descriptor vectors saved as JSON, enabling geometric homology search
- **Basin-control null** — Markov-walk null model preserving Ramachandran basin identity to test whether curvature patterns exceed basin geometry
- **Superpotential landscape** — energy-like surface on T² derived from curvature accumulation

## Key findings (n = 200 human proteins)

| SS Type | Non-constant κ | Basin null | p-value | Direction |
|---------|---------------|------------|---------|-----------|
| α-Helix | 65.0% | 45.5% | 1.4 × 10⁻⁵⁸ | **Enriched** |
| β-Strand | 30.0% | 11.4% | 1.7 × 10⁻⁴⁴ | **Enriched** |
| Coil | 29.5% | 34.6% | 5.5 × 10⁻⁷ | **Suppressed** |

Helices are the most curvature-modulated secondary structure on T². Coil is the only type smoother than its basin null. ~90% of proteins wind into negative ψ.

## Installation

```bash
git clone https://github.com/[you]/PhiPsi_Atlas.git
cd PhiPsi_Atlas
pip install numpy scipy
```

Optional but recommended:
```bash
pip install requests    # AlphaFold downloads
pip install biopython   # mmCIF/PDB parsing (fallback: built-in parser)
```

No GPU required. Runs on CPU. Tested on Python 3.9+.

## Quick start

### Analyze a single PDB structure
```bash
python cornu_spirals_on_T2.py            # runs on bundled test proteins
python cornu_spirals_on_T2.py pdb 1UBQ   # fetch and analyze ubiquitin
```

### Analyze AlphaFold structures
```bash
# Single protein by UniProt ID
python cornu_spirals_on_T2.py alphafold P00533

# Batch from file (one UniProt ID per line)
python cornu_spirals_on_T2.py alphafold @human_2000.txt
```

If you have a local AlphaFold cache (bulk download), place it at `alphafold_cache/` and the pipeline will find structures there before attempting any network download.

### Run the basin-control null
```bash
# Quick test (~5 min, 200 proteins, 2 null reps)
python basin_null_alphafold.py --quick

# Full power (~hours, all cached structures, 5 reps)
python basin_null_alphafold.py --reps 5

# Custom
python basin_null_alphafold.py --max-proteins 500 --reps 3
```

### Inspect results
```bash
# Dashboard with distributions, winding maps, cluster analysis
python analyze_barcodes.py

# Specific queries
python analyze_barcodes.py --find P00533          # look up a protein
python analyze_barcodes.py --top-defects 20       # most defective proteins
python analyze_barcodes.py --top-geodesic 20      # smoothest proteins
python analyze_barcodes.py --winding-map           # ASCII (p,q) scatter
python analyze_barcodes.py --clusters              # group by SS composition
python analyze_barcodes.py --defect-catalog        # all Gaussian peaks ranked
```

## Project structure

```
PhiPsi_Atlas/
├── cornu_spirals_on_T2.py      # Core library + CLI entry point
│                                  Sections 1–15: math, curvature, classification,
│                                  winding, barcodes, AlphaFold ingestion, reporting
├── basin_null_alphafold.py      # Basin-control Markov null analysis
├── analyze_barcodes.py          # Barcode database dashboard + queries
├── human_2000.txt               # Example UniProt ID list (human proteome sample)
├── alphafold_cache/             # Local AlphaFold structure cache (user-provided)
├── data/
│   └── alphafold/               # Downloaded structures (auto-created)
└── results/
    ├── barcodes/                # Per-protein JSON barcodes
    ├── basin_null_results.json  # Null analysis output
    └── reports/                 # Per-protein text reports
```

## Geometric barcode format

Each protein produces a JSON barcode in `results/barcodes/`:

```json
{
  "uniprot_id": "P68871",
  "chain_length": 147,
  "p_winding": 3.21,
  "q_winding": -8.44,
  "Q_magnitude": 9.03,
  "ss_composition": {"H": 0.76, "E": 0.0, "C": 0.24},
  "geodesic_fraction": 0.32,
  "defect_fraction": 0.18,
  "segments": [
    {
      "ss_type": "H",
      "spiral_class": "gauss_peak",
      "length": 12,
      "mean_kappa": 0.041,
      "dkappa_ds": 0.003,
      "r_squared": 0.87,
      "start_residue": 4,
      "end_residue": 15
    }
  ]
}
```

## The ten curvature models

| # | Model | κ(s) | Macro-category |
|---|-------|------|----------------|
| 1 | Geodesic | ≈ 0 | constant_k |
| 2 | Circular arc | c | constant_k |
| 3 | Clothoid | a + bs | monotone |
| 4 | Quadratic | a + bs + cs² | polynomial |
| 5 | Gaussian peak | A·exp(−(s−s₀)²/2σ²) + c | localized |
| 6 | Sinusoidal | a + b·sin(ωs + φ) | oscillatory |
| 7 | Damped oscillation | a + be⁻ᶜˢ·sin(ωs + φ) | oscillatory |
| 8 | Exponential | aeᵇˢ + c | monotone |
| 9 | Sigmoid | L/(1+e⁻ᵏ⁽ˢ⁻ˢ⁰⁾) + c | transition |
| 10 | Step | piecewise constant | transition |

Best model selected by corrected Akaike Information Criterion (AICc). For null comparison, models collapse to six macro-categories.

## How the basin-control null works

Standard residue-shuffle nulls destroy all sequential correlation and trivially reject. The basin-control null asks a harder question: **does the observed curvature exceed what Ramachandran basin geometry alone predicts?**

For each real segment of length L in SS type S:
1. Build a step pool of (Δφ, Δψ) from all real segments of the same SS type
2. Generate surrogate paths by Markov walk from random SS-matched starting points
3. Classify surrogates identically to real data
4. Compare distributions by χ² and binary (constant vs non-constant) contingency

This preserves basin identity, local step-size statistics, and segment lengths while randomizing path ordering.

## Performance

- **Single protein**: ~2–5 seconds (depending on chain length)
- **Batch 200 proteins**: ~10 minutes
- **Basin null (200 proteins, 2 reps)**: ~15 minutes
- **Full proteome (20k proteins)**: ~12 hours for barcodes, ~24 hours for null


## Citation

If you use PhiPsi_Atlas in your work:

```
@software{phipsi_atlas,
  title = {PhiPsi_Atlas: Differential geometry of protein backbones on the Ramachandran torus},
  year = {2026},
  url = {https://github.com/kase1111-hash/PhiPsi_Atlas}
}
```

Preprint: *Curvature Modulation on the Ramachandran Torus: Secondary Structure Signatures at Proteome Scale* (in preparation)

## License

MIT

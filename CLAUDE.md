# CLAUDE.md

## Project

Monte Carlo framework proving the Standard Genetic Code (SGC) is statistically
optimized. Two metrics — Dirichlet energy on the codon Hamming graph and noise
distortion from single-nucleotide misdecoding — show the SGC outperforms
>99.9999% of 1M random codes with identical degeneracy structure.

Companion code for: "Triplet architecture enables deep error-minimization in
the genetic code" (JME submission).

## Commands

```bash
# build the 22-property parquet (prerequisite for phase-2 runs)
python build_aa_props.py

# phase-1 sanity check (~2 min, single-process)
python run_phase1_demo.py

# phase-2 triplet (1M samples, ~2 min on 64 cores)
python run_phase2_batch.py --n 1000000 --workers 32 --out results/phase2_triplet.auto.json

# phase-2 doublet (n=2, reduced alphabet)
python run_phase2_n2n4.py --n-bases 2 --n 1000000 --workers 32 --raa raa10 --out results/phase2_doublet.raa10.json

# 2x2 factorial + sensitivity (generates all manuscript numbers + raw distributions)
python run_publication_controls.py --n-null 1000000 --workers 32

# verify synonymy numbers
python verify_synonymy.py

# run tests
python -m pytest tests/ -v

# generate figures
python figures/fig1_optimization_factorial.py
python figures/fig2_wobble_mechanism.py
python figures/figS1_raw_distributions.py
```

SLURM scripts (`slurm_triplet.sbatch`, `slurm_doublet.sbatch`,
`slurm_quadruplet.sbatch`) request 32 CPUs / 64G RAM.

## Architecture

### Module layout

- **`src/io/`** — data loading. `codon_io.py` loads the SGC table; `aa_props_io.py`
  handles feature loading, correlation pruning, PCA; `aa_props_lib.py` is the raw
  22-property database.
- **`src/sims/`** — simulation primitives. `codon_graph.py` builds the Hamming graph;
  `random_codes.py` has degeneracy-preserving shuffle generators.
- **`src/metrics/`** — `dirichlet.py` computes E = Σ_edges ||p(u)-p(v)||²;
  `mi.py` computes mutual information by codon position.
- **`src/receiver/`** — `thermo_noise.py` computes noise distortion via Boltzmann
  softmax over Hamming-1 neighbors.

### Data flow

```
aa_props_lib.py (22 properties)
  → build_aa_props.py (auto-select, PCA)
  → data/processed/aa_props.parquet (8 PCs, 97.1% variance)

codon_table.csv (64 codons → AA)
  → codon_graph.py (Hamming graph: 4^n nodes, edges at Hamming distance 1)

parquet + codon_table.csv
  → run_publication_controls.py
  → 1M degeneracy-preserving shuffles per condition
  → results/publication_controls.json + results/null_distributions.npz
```

## Key design details

- **Dirichlet energy** uses squared L2 distance; **noise distortion** uses unsquared L2.
- The noise distortion in `run_publication_controls.py` uses **uniform weighting**
  (all error pathways weighted equally). The separate Boltzmann-weighted noise model
  in `thermo_noise.py` is used by `run_phase2_batch.py`.
- The null model permutes AA labels among **individual sense codons** (not codon
  blocks), preserving per-AA codon counts. Block structure is NOT preserved.

## Null model

Degeneracy-preserving shuffle: permutes AA labels among 61 sense codons while
keeping per-AA codon counts fixed (Leu=6, Ser=6, ..., Met=1, Trp=1). Stops
remain stops. Percentile formula: `(count + 1) / (n + 1)`.

## Canonical artifacts

- `results/phase1_demo.json`
- `results/phase2_triplet.auto.json`
- `results/phase2_doublet.raa10.json`
- `results/phase2_quadruplet.auto.json`
- `results/publication_controls.json`
- `results/sgc_baselines.json`

## Import pattern

All scripts run from repo root. `conftest.py` adds `src/` to `sys.path`.
Imports use `from src.io.codon_io import ...` style.

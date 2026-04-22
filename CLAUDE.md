# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Monte Carlo framework proving the Standard Genetic Code (SGC) is statistically
optimized. Two metrics — Dirichlet energy on the codon Hamming graph and noise
distortion from single-nucleotide misdecoding — show the SGC outperforms
>99.9999% of 1M random codes with identical degeneracy structure.

Companion code for: "Triplet architecture enables deep error-minimization in
the genetic code" (JME submission).

## Commands

```bash
# prerequisite: build the 22-property parquet (must run before phase-2)
python build_aa_props.py

# phase-1 sanity check (~2 min, single-process)
python run_phase1_demo.py

# phase-2 triplet (1M samples, ~2 min on 64 cores) — Boltzmann-weighted noise via thermo_noise.py
python run_phase2_batch.py --n 1000000 --workers 32 --out results/phase2_triplet.auto.json

# phase-2 doublet (n=2, reduced alphabet)
python run_phase2_n2n4.py --n-bases 2 --n 1000000 --workers 32 --raa raa10 --out results/phase2_doublet.raa10.json

# 2x2 factorial + sensitivity (generates all manuscript numbers + raw distributions)
python run_publication_controls.py --n-null 1000000 --workers 32

# 2 extra random seeds (for Table S1 mean ± SD over 3 seeds)
python run_sensitivity_extra_seeds.py

# preemptive defense: f^T L^2 f + Boltzmann-weighted D for all 3 conditions
python run_preempt_metrics.py --n-null 100000 --workers 32

# verify synonymy numbers match manuscript
python verify_synonymy.py

# build Table S1 (CSV + LaTeX, aggregates closest/farthest/3-seed random)
python build_table_s1.py

# run all tests
python -m pytest tests/ -v

# run a single test file or test function
python -m pytest tests/test_synonymy.py -v
python -m pytest tests/test_graph.py::test_hamming -v

# generate all manuscript figures (fig1–fig4 + figS1)
python figures/fig1.py
python figures/fig2.py
python figures/fig3.py
python figures/fig4.py
python figures/figS1_raw_distributions.py
```

SLURM scripts (`slurm_triplet.sbatch`, `slurm_doublet.sbatch`,
`slurm_quadruplet.sbatch`) request 32 CPUs / 64G RAM.

## Architecture

### Module layout

- **`src/io/`** — data loading. `codon_io.py` loads the SGC table; `aa_props_io.py`
  handles feature loading, correlation pruning, PCA; `aa_props_lib.py` is the raw
  22-property database (hardcoded canonical values for 20 amino acids).
- **`src/sims/`** — simulation primitives. `codon_graph.py` builds the Hamming graph
  (`build_graph_n(n)` creates 4^n nodes with edges at Hamming distance 1);
  `random_codes.py` has degeneracy-preserving shuffle generators.
- **`src/metrics/`** — `dirichlet.py` computes E = sum_edges ||p(u)-p(v)||^2;
  `mi.py` computes mutual information by codon position (used in phase-1 only).
- **`src/receiver/`** — `thermo_noise.py` computes noise distortion via Boltzmann
  softmax over Hamming-1 neighbors. Has both a triplet-specific version and a
  general version (`_multi_general`) for arbitrary codon lengths.

### Data flow

```
aa_props_lib.py (22 properties, hardcoded)
  -> build_aa_props.py (auto-select, PCA)
  -> data/processed/aa_props.parquet (8 PCs, 97.1% variance)

codon_table.csv (64 codons -> AA)
  -> codon_graph.py (Hamming graph: 4^n nodes, edges at Hamming distance 1)

parquet + codon_table.csv
  -> run_publication_controls.py
  -> 1M degeneracy-preserving shuffles per condition
  -> results/publication_controls.json + results/null_distributions.npz
```

### Import pattern

All scripts run from repo root. `conftest.py` adds `src/` to `sys.path`.
Imports use `from src.io.codon_io import ...` style.

## Key design details

- **Dirichlet energy** uses **squared** L2 distance; **noise distortion** uses
  **unsquared** L2. Do not change one to match the other.
- Two noise distortion implementations exist for different purposes:
  - `run_publication_controls.py` computes distortion inline with **uniform weighting**
    (all error pathways weighted equally). This is the canonical version for manuscript
    numbers.
  - `thermo_noise.py` uses **Boltzmann-weighted** probabilities with wobble/asymmetry
    parameters. Used by `run_phase2_batch.py`.
- The null model permutes AA labels among **individual sense codons** (not codon
  blocks), preserving per-AA codon counts. Block structure is NOT preserved.
  AUG->Met is pinned. Stops remain stops.
- Percentile formula uses Laplace smoothing: `(count + 1) / (n + 1)`.

### RAA10 reduced alphabet

The 2x2 factorial design uses a 10-class Reduced Amino Acid alphabet (RAA10)
defined in both `run_publication_controls.py` and `verify_synonymy.py`. The
mapping groups 20 AAs into 10 classes (e.g., {A,G}->X1, {V,L,I,M}->X2, etc.).
Representative selection for each class (closest/farthest/random to centroid)
is the subject of the sensitivity analysis.

## Canonical artifacts

### Numerical results

- `results/phase1_demo.json`
- `results/phase2_triplet.auto.json` — Boltzmann-weighted noise, n=3 (for
  cross-check against the uniform-weighted manuscript numbers)
- `results/phase2_doublet.raa10.json`
- `results/phase2_quadruplet.auto.json`
- `results/publication_controls.json` — **the definitive manuscript numbers**
  (factorial + sensitivity at 1M shuffles per condition)
- `results/sensitivity_extra_seeds.json` — random-strategy seeds 1 and 2
  (seed 42 is in `publication_controls.json`); feeds Table S1's 3-seed mean
- `results/preempt_metrics.json` — preemptive defense: f^T L^2 f and
  Boltzmann-weighted D vs the manuscript's f^T L f and uniform D
- `results/sgc_baselines.json`
- `results/null_distributions.npz` — gitignored (46 MB), regenerate with
  `run_publication_controls.py`
- `supplementary/table_S1.{csv,tex}` and `table_S1_summary.json`

### Figures (canonical manuscript set)

| Manuscript | Script | Deliverables |
|:-:|---|---|
| Figure 1 | `figures/fig1.py` | framework: Hamming graph, PCA space, GCU schematic |
| Figure 2 | `figures/fig2.py` | SGC null: E histogram, D histogram, joint D-E scatter |
| Figure 3 | `figures/fig3.py` | factorial overlay + position-specific synonymy bars |
| Figure 4 | `figures/fig4.py` | wobble mechanism: GCU triplet + GC doublet neighborhoods |
| Figure S1 | `figures/figS1_raw_distributions.py` | per-condition raw null distributions |

Each script outputs `.pdf` (vector, fonttype=42 so text stays editable) and
`.png` (300 DPI). Each has a matching `fig{N}_caption.txt` draft.

### Manuscript and supplementary prose

- `manuscript_JME.md`
- `cover_letter_JME.md`
- `submission_checklist_JME.md`
- `supplementary_tables_JME.md`
- `supplementary_materials_JME.md` — combined, submission-ready
- `literature_comparison.md` — verbatim-equation comparison to Tlusty
  2007/2010 and Radványi & Kun 2021; drives the "what is novel" framing
- `rebuild_notes.md` — decisions, discrepancies, preemptive-defense numbers
- `HISTORY.md` — append-only log of Q&A in this project

## Tests

8 tests across 3 files:
- `test_graph.py` — Hamming graph construction and distance
- `test_phase1_compat.py` — codon loader preserves base columns; MI computation runs
- `test_synonymy.py` — position-specific synonymy regression values for all 3
  conditions (A/B/C) plus doublet class count assertion

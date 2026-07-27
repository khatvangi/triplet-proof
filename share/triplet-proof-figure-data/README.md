# Figure data for redrawing in Excel / Origin

This bundle contains the numerical data behind the figures in the paper
"Triplet architecture enables deep error-minimization in the genetic code"
(Khatvani et al., submitted to *Journal of Molecular Evolution*).

All files are plain CSV with a header row. Decimal separator: period (.).
Field separator: comma (,). Encoding: UTF-8.

## Files

### Input data (what the pipeline consumes)

- **codon_table.csv** — the standard genetic code: 64 codons → amino acid
  (column `aa`), plus `is_stop` and `is_start` flags. 64 rows.

- **aa_props_raw.csv** — 22 physicochemical descriptors for 20 amino acids
  (hydropathy, volume, polarity, pI, charge, hydrogen-bond donors /
  acceptors, aromatic flag, rotatable bonds, helix / sheet / turn
  propensities, etc.). 20 rows × 23 columns (aa + 22 descriptors).

- **aa_props_PCs.csv** — the 8 principal components derived from
  `aa_props_raw.csv` after correlation pruning and z-scoring. These are
  the actual property vectors used in all distance calculations in the
  paper (explains 97.1 % of the variance of the 22 descriptors). 20 rows
  × 9 columns (aa + pca1..pca8).

### Null distributions (100K-sample subsamples of the 1M-sample runs)

Each of these is a random subsample of 100 000 rows from the full
1 000 000-sample null distributions in `results/null_distributions.npz`
of the code repository. Subsample uses numpy RNG with seed
42 so the same 100 000 rows appear in every re-run of the
export script. For figure redrawing this is visually indistinguishable
from the full 1M; for exact z-score reproduction use the full 1M in the
repo.

Each null file has two columns, D and E:
- `D` = noise distortion (per-codon average of ‖Δf‖ over non-synonymous
  Hamming-1 neighbours; unsquared L2 distance in the 8-PC space)
- `E` = Dirichlet energy (sum over graph edges of ‖Δf‖² in the 8-PC
  space)

- **null_A_doublet.csv** — Condition A: 10-AA doublet projection, n=2.
  100 000 rows × 2 columns.
- **null_B_SGC.csv** — Condition B: 20-AA SGC, n=3.
  100 000 rows × 2 columns.
- **null_C_10AA_triplet.csv** — Condition C: 10-AA RAA10 triplet, n=3.
  100 000 rows × 2 columns.

### SGC observed values and statistics

- **sgc_baselines.csv** — the SGC's observed D and E for each condition
  (the red vertical line / red star in Figures 2, 3, S1). 3 rows × 4
  columns.

- **z_scores_summary.csv** — z-scores, p-values, null means and standard
  deviations, count of random codes with value ≤ SGC, and D–E correlations
  for each condition. Everything needed to reproduce the annotations in
  the figures. 3 rows × 15 columns.

- **synonymy_fractions.csv** — position-specific synonymous-substitution
  fraction (%) for each condition. Drives the bar chart in Figure 3b.
  Empty cell for Condition A position 3 because doublet codons have no
  third position. 3 rows × 4 columns.

- **table_S1.csv** — sensitivity of Condition C's z-score to the choice
  of representative amino acid per RAA10 class (closest / farthest /
  random-3-seed-mean). Drives Table S1 in the supplementary materials.
  3 rows × 4 columns.

## Which file feeds which figure

- **Figure 1a** (Hamming graph) — built from `codon_table.csv` alone,
  plus a layout algorithm. Needs a graph / network tool, not Excel /
  Origin.
- **Figure 1b** (PCA space of 20 AAs) — `aa_props_PCs.csv` columns
  `pca1` and `pca2`; colour-code by amino-acid class (see paper).
- **Figure 1c** (GCU neighborhood schematic) — `codon_table.csv` +
  `aa_props_PCs.csv`; computed edge-wise distances reported in the
  figure text. Hand-drawn layout.

- **Figure 2a** (Dirichlet energy null) — `null_B_SGC.csv` column `E`
  as a histogram; vertical line at `E_SGC` from `sgc_baselines.csv`
  condition B; annotation from `z_scores_summary.csv` column
  `E_z_score`.
- **Figure 2b** (Noise distortion null) — `null_B_SGC.csv` column `D`
  as a histogram; vertical line at `D_SGC`; annotation from
  `D_z_score`.
- **Figure 2c** (Joint D–E scatter) — `null_B_SGC.csv` as (D, E)
  scatter; red star at (`D_SGC`, `E_SGC`) for condition B; correlation
  `r` = `DE_correlation` from `z_scores_summary.csv`.

- **Figure 3a** (three-condition D-histogram overlay) — all three
  `null_*.csv` files, column `D` from each, plotted on a single x-axis
  with transparency, with three vertical lines at the respective
  `D_SGC` values from `sgc_baselines.csv`.
- **Figure 3b** (synonymy bars) — `synonymy_fractions.csv`.

- **Figure 4** (wobble schematic) — conceptual; built from
  `codon_table.csv` neighbour relations plus hand-drawn layout.

- **Figure S1** (raw distributions per condition) — the three
  `null_*.csv` files as separate histograms, each with its `D_SGC` red
  vertical line. For an honest cross-condition comparison plot all
  three on the same x-axis range, e.g. D ∈ [17, 36].

## If you need the full 1M-sample null distributions

Clone the code repository at https://github.com/khatvangi/triplet-proof
and look in `results/null_distributions.npz` — that's a NumPy zip file
with six 1 000 000-element arrays (A_D, A_E, B_D, B_E, C_D, C_E). The
script that regenerates this file is `run_publication_controls.py`.

## Contact

If anything is unclear or you find a discrepancy between a file here
and the corresponding figure in the paper, email khatvangi@...

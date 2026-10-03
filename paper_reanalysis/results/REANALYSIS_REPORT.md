# Reanalysis of the submitted JME manuscript: the null model, not the chemistry

Run: `scripts/run_reanalysis.py --n 1000000`, seed `20260923`, 1007 s.
Plan frozen before any cell was computed: `ANALYSIS_PLAN.md` (commit `da61661`).
Implementation: commit `d9b657b`. Grid: `results/grid.csv`. Figure:
`results/fig_null_vs_chemistry.{png,pdf}`.

This is a re-scoring of the submitted analysis, not a new experiment. No claim
below rests on a number that was not produced by the run recorded here.


## 1. The gate

Per the plan, nothing was scored until the local implementation reproduced the
manuscript's own SGC values through the same code path used for every null draw:

| quantity | published | reproduced here | absolute error |
|---|---:|---:|---:|
| `E` | 2984.6754069246 | 2984.6754069246 | 4.5e-13 |
| `D` | 23.9223597432532 | 23.9223597432532 | 2.5e-14 |

The gate passed. Had it failed the run was to stop and report nothing.


## 2. The grid

Two nulls x five property representations x two objectives, all through one
scoring function, 10^6 draws per cell.

**z (negative = the standard code scores better than the null)**

| properties | E, codon null | D, codon null | E, block null | D, block null |
|---|---:|---:|---:|---:|
| 8 PCs, whitened — *what the paper computed* | **-14.60** | **-16.56** | -1.92 | -2.54 |
| 8 PCs, unwhitened | -12.81 | -15.55 | -2.17 | -3.01 |
| 20 descriptors, pre-PCA | -12.97 | -15.85 | -2.19 | -2.97 |
| polar requirement, 1-D | -7.85 | -11.59 | -2.75 | -3.40 |
| **chemistry-blind**, all AA pairs at distance 1 | **-17.07** | **-17.07** | *invariant* | *invariant* |

The two nulls:

* **codon null** — permute amino-acid labels among the 61 individual sense
  codons, preserving each amino acid's codon count. This is what
  `run_publication_controls.py` actually samples.
* **block null** — permute the 20 amino acids among the 20 synonymous blocks
  (`L = pi[L_sgc]`; Haig & Hurst, Freeland & Hurst). This is what the submitted
  Methods paragraph describes.


## 3. Three findings

### 3.1 The published z-scores do not measure chemistry

Under the codon null, the **chemistry-blind** representation — in which every
pair of amino acids is placed at distance 1, so the score contains no
physicochemical information whatsoever — returns `z = -17.07` on both
objectives. That is **more extreme than the published -14.60 and -16.56**. It is
the most extreme cell in the entire grid.

A statistic that a chemistry-free representation maximises is not a measurement
of chemical optimisation. What the codon null actually detects is the standard
code's block structure: the SGC places 67 of its 263 sense-sense Hamming-1 edges
inside synonymous blocks, while a codon-level permutation places a mean of 12.51
(263 - 250.4945, read off this run's own blind codon-null mean). Any distance function at all, chemical or not, scores those 67 zero-cost
edges as savings.

This single control is, in my judgement, sufficient on its own to retract the
optimisation reading of the headline numbers. It requires no argument about
which null is correct.

### 3.2 Under the block null the blind control goes to exactly nothing — as it must

`E` under the block null is **exactly invariant**: 1,000,000 of 1,000,000 draws
tie with the standard code, `p = 1`. The mechanism is transparent — with all
amino-acid pairs at distance 1, `E` degenerates into a count of *non-synonymous*
edges, and the block null holds the block partition fixed, so that count cannot
move. Its value is 196, which is exactly the cross-block edge count independently
reported in `code_evolution/results/block_null/REPORT.md` §4.1 (263 = 196
cross-block + 67 synonymous). Two independently written code paths agree on an
integer.

This is the behaviour a correct null must show: the negative control must return
nothing. The codon null fails that test; the block null passes it.

*(Note on the table: `null_sd` for this cell prints as 7.1e-15 rather than 0.
That is rounding in the mean of 10^6 identical float64 values, not variance.
`zstats` now flags invariance against the scale of the statistic and counts ties
against the standard code; before that fix this cell printed `z = 1.0` and
`p = 1e-6`, which would have appeared in the table as the strongest result in
the grid when it is the weakest possible one.)*

### 3.3 Chemistry survives, at 2-3.4 z rather than 14-17

Under the block null the chemical signal is real but an order of magnitude
smaller than published:

| representation | objective | z | p (Laplace) | draws better than the SGC |
|---|---|---:|---:|---:|
| 8 PCs, whitened | `E` | -1.92 | 0.0225 | 22,535 / 10^6 |
| 8 PCs, whitened | `D` | -2.54 | 0.0053 | 5,293 / 10^6 |
| 8 PCs, unwhitened | `D` | -3.01 | 0.0011 | 1,107 / 10^6 |
| polar requirement, 1-D | `E` | -2.75 | 1.5e-4 | 152 / 10^6 |
| polar requirement, 1-D | `D` | -3.40 | 2.1e-5 | 19 / 10^6 |

The standard code is beatable by direct sampling in every one of these cells.
The strongest honest signal in the grid comes from polar requirement — a single
one-dimensional scale — not from the eight-component PCA the paper uses.


## 4. The undisclosed whitening cuts against the paper in both directions

`aa_prop_vector_map(..., standardize=True)` z-scores the PC scores, which
re-weights the eight components to equal variance. The manuscript does not
disclose this step. Its effect is not neutral:

| null | whitened | unwhitened | effect of whitening |
|---|---:|---:|---:|
| codon (`E`) | -14.60 | -12.81 | **inflates by 1.79 z** |
| block (`E`) | -1.92 | -2.17 | **deflates by 0.25 z** |

So the undisclosed transformation makes the defective statistic look stronger
and the defensible statistic look weaker. I have no evidence this was
deliberate, and the analysis plan did not predict the direction; it is reported
because it is what the grid shows.


## 5. Cross-validation

The `pc8_whitened` block-null cell reproduces the independently preregistered
10^7-draw run (`code_evolution/results/block_null/`, seed `20260903`):

| objective | this run (10^6) | preregistered (10^7) | difference |
|---|---:|---:|---:|
| `E` | -1.9164 | -1.9153 | 0.0011 |
| `D` | -2.5385 | -2.5369 | 0.0016 |

Different seed, different sample size, different implementation, same answer.

**Which matrix this belongs to.** Both figures use the hand-compiled descriptor
matrix in `data/processed/aa_props.parquet`. The manuscript's descriptor matrix
was later rebuilt from named AAindex accessions, which moves the block-null
values to `E` z = -1.90 (p = 0.027) and `D` z = -2.34 (p = 0.010). This
cross-validation therefore certifies that two implementations agree on one
matrix; it is not a check on the rebuilt numbers. The chemistry-blind cell is
the part that carries across, because it contains no chemistry and cannot move
when the matrix changes.

**Not cross-validated:** the `prs1d` rows here are `E` and `D` computed on a
one-dimensional polar-requirement feature. They are *not* the `PRS` objective of
that report (`z = -3.7397`, φ6 transition/transversion weighting). Same sign and
same order of magnitude; not a replication. Do not quote them as one.


## 6. What the paper may and may not now claim

**May not:**

* that the standard code beats >99.9999% of random codes on chemical distance —
  the codon-level figure is a synonymy measurement, and a chemistry-blind
  control beats it;
* `z = -14.6` / `z = -16.6` as evidence of chemical optimisation;
* the printed `D` formula, which carries a `1/|N(c)|` the code never applies
  (`corr(D, H1) = 1.0000000000`).

**May, with the block null and full disclosure:**

* that the standard code's assignment of amino acids to its *existing* block
  structure is better than chance on chemical distance, at `p = 0.023` (`E`) and
  `p = 0.0053` (`D`), and `p = 2.1e-5` on polar requirement;
* that the effect is modest and the code is beatable by direct sampling — tens
  of thousands of random block permutations per million score better;
* that the large published effect is attributable to the block structure itself,
  which the codon-level null destroys and therefore credits to chemistry.


## 7. What this run does not settle

* It does not establish *why* the block structure exists. Nothing here separates
  selection on synonymy from a historical or biosynthetic origin for it.
* It does not test the manuscript's `n = 2` / `n = 3` architecture comparison.
  That is a separate rerun and is not in this grid.
* It does not address the reviewer's objections outside the null-model issue.
* 10^6 draws bound the codon-null cells at `p < 1e-6` but cannot resolve them
  further; the preregistered 10^7 run is the reference for the block cells.


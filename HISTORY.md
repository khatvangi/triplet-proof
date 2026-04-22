## 2026-04-22

**Q:** Pre-compute the two defensive metrics: f⁺ L² f (Tlusty-style squared
Laplacian) and Boltzmann-weighted D (thermo_noise.py). Verify qualitative
agreement with the manuscript's f⁺ L f and uniform-D.

**A:** Ran `run_preempt_metrics.py` at 100K shuffles × 3 conditions. Output
at `results/preempt_metrics.json`. Details in `rebuild_notes.md`. Summary:
- **f⁺ L² f tracks f⁺ L f tightly** (null r = 0.987; z-scores differ by
  ≤1.1 units across all conditions). Clean defense against a Tlusty-style
  reviewer.
- **Boltzmann-weighted D is NOT qualitatively identical to uniform D** —
  it gives z = −24.9 for the SGC (vs −16.6 uniform) and z = −20.1 for the
  10-AA triplet control (vs −14.8). The realistic model STRENGTHENS the
  SGC's advantage by 5–8 z-units on triplets. For doublet (A) they
  coincide (no wobble position). This is actually better than "identical" —
  reporting the uniform version is genuinely conservative; a mechanism-
  minded reviewer would compute a *stronger* result, not a weaker one.

## 2026-04-21 (later)

**Q:** Pre-submission supplementary cleanup: reconcile Table S1 between the
stale `supplementary_tables_JME.md` and the 3-seed CSV/TeX; write Fig S1
caption; produce one submission-ready `supplementary_materials_JME.md`.

**A:** Three files produced / updated:
- `supplementary_tables_JME.md` — Table S1 replaced with 3-seed protocol
  values matching `table_S1.csv/.tex`. Table S2 untouched.
- `figures/figS1_caption.txt` — new caption (CV values confirmed current).
- `supplementary_materials_JME.md` — single submission-ready file combining
  Table S1, Table S2, Figure S1 reference + caption, and a file manifest.

Discrepancies found and resolved:
- **Stale Table S1 z-scores**: the old file had closest (−14.87, −13.41),
  farthest (−14.57, −13.12), random (−15.08, −13.73). Current CSV/TeX /
  `publication_controls.json` have closest (−14.79, −13.34), farthest
  (−14.61, −13.18). Differences of ~0.08 z-units suggest the old file
  came from a prior Monte Carlo run; new file uses current values.
- **Stale D/E ranges**: old said D-range = 0.51, E-range = 0.61 (single-seed
  random). With 3-seed random-mean: D-range = 0.56, E-range = 0.92.
- **DE correlations transferred cleanly**: closest 0.958, farthest 0.965,
  and single-seed random 0.955 (stale file) all match current
  `publication_controls.json` to 3 decimals. Added the random 3-seed mean
  DE corr = 0.961 ± 0.008 (computed from seed 42 = 0.955, seed 1 = 0.969,
  seed 2 = 0.958).
- **CV values in Fig S1 verified against `publication_controls.json`**:
  A = 0.0343 (shown 3.4 %), C = 0.0211 (shown 2.1 %), B = 0.0164 (shown
  1.6 %). All current. No regeneration of the figure needed.

## 2026-04-21

**Q:** Verify whether our Dirichlet energy and noise distortion metrics have
prior art in the genetic code literature. Compare to Tlusty 2007, Tlusty 2010,
and Radványi & Kun 2021. Don't make citation decisions — just report.

**A:** Full comparison written to `literature_comparison.md`. Headline findings:
- Our **Dirichlet energy** E = Σ ‖f(u)−f(v)‖² = fᵀ L f is the standard
  unweighted graph Laplacian quadratic form. Tlusty 2007 Eq. B1 defines
  exactly the same weighted Laplacian on the codon error-graph; our E is
  the unweighted / uniform-misreading specialization. Not novel.
- Our **noise distortion** D = (1/|C|) Σ_c Σ_{c'∈N(c)} ‖f(c)−f(c')‖
  (unsquared) has the same functional form as Tlusty's "error-load"
  (2007 Eq. 1, 2010 Eq. 1) and Radványi & Kun's 2021 "distortion" (Eq. 1).
  Tlusty uses weighted misreading matrix + probabilistic code matrix;
  Radványi & Kun use codon usage + Kimura mutation model. Our version
  specializes to uniform codon prior + unit-weight Hamming-1 + multi-D PCA
  distance. Not novel as a metric; specific parameter choices differ.
- **What IS new in our work**: degeneracy-preserving null model (empirical
  z-scores and p-values), 2×2 factorial alphabet-vs-architecture design,
  10⁶-sample Monte Carlo, PCA-reduced multi-dim property vector.
- Recommended positioning: cite Tlusty as framework source, Radványi & Kun
  as recent related application; frame contribution as factorial design +
  empirical null, not introduction of the metrics.

## 2026-04-21 (earlier)

**Q:** Full rebuild of all four manuscript figures + Table S1 per Paudyal's
review. Keep existing color palette. Three random seeds for Table S1. Verify
z-scores match manuscript. Flag any discrepancies.

**A:** Full rebuild delivered:
- `figures/fig{1,2,3,4}.{pdf,png}` + `fig{1,2,3,4}_caption.txt`
- `supplementary/table_S1.{csv,tex}` (closest, farthest, random mean ± SD
  over 3 seeds)
- `rebuild_notes.md` summarizing decisions and discrepancies
- Verified all manuscript z-scores match `publication_controls.json` ✓
- **Flagged two real discrepancies**: (1) Fig 1c prose claimed ‖Δf‖ = 8.7
  for "radical" but max pairwise distance in PC space is 6.23 — used actual
  values 0 and 3.5–4.9 for Ala neighborhood; (2) Methods "<0.5 units"
  sensitivity claim is false for both metrics once the random strategy is
  sampled with 3 seeds — D range 0.56, E range 0.92. Table S1 reports
  honest numbers.
- Design changes: Fig 3 now single-x-axis overlay (not stacked mini-hists)
  without Δz brackets; Fig 4 uses colored+bold substituted bases instead
  of underlines (matplotlib mathtext doesn't support `\underline`).

## 2026-04-14

**Q:** Generate all manuscript figures for MBE Letter submission — extract .npy data from .npz archive, fix and run existing figure scripts, create joint D-E scatter, create PCA amino acid space panel, and assemble final manuscript-numbered figure files.

**A:** Completed all figure generation:
- Extracted 3 `.npy` files from `null_distributions.npz` (each 1M x 2 with [D, E] columns)
- Verified 3 existing scripts run cleanly: `fig1_optimization_factorial.py`, `fig2_wobble_mechanism.py`, `figS1_raw_distributions.py`
- Created `fig_joint_scatter.py` — scatter of D vs E for Condition B with marginals and SGC star
- Created `fig1b_pca_space.py` — PC1 vs PC2 of 20 AAs colored by RAA10 grouping
- Created manuscript-numbered composites:
  - `fig2_sgc_optimization.py` → panels (a) histogram + (b) joint scatter
  - `fig3_factorial_synonymy.py` → panels (a) stacked mini-hists with brackets + (b) synonymy bars
  - `fig4_wobble_mechanism.py` → wrapper calling existing wobble schematic
- Figure 1 panels (a) and (c) left for hand-drawing in Illustrator; only data-driven panel (b) generated
- All tests still pass (8/8)

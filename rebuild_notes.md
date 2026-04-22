# Figure rebuild — notes

Date: 2026-04-21
Scope: full rebuild of Figures 1–4 and Table S1 for the JME submission.
Previous state: the pre-rebuild `fig1_optimization_factorial.py`,
`fig2_wobble_mechanism.py`, `figS1_raw_distributions.py` still exist and run
cleanly. Manuscript-numbered scripts (`fig1.py`, `fig2.py`, `fig3.py`,
`fig4.py`) are the canonical deliverables for the current submission.

## Discrepancies between the manuscript and the data

1. **Fig 1 panel (c), ‖Δf‖ values.** The rebuild instructions specified
   "‖Δf‖ = 0 for synonymous, 1.2 for conservative, 8.7 for radical". These
   numbers do not match the data.

   - Maximum pairwise ‖Δf‖ across all 20 AAs in the 8-PC space is 6.23
     (G – P pair). A value of 8.7 is outside the achievable range of this
     property space.
   - For the actual Ala (GCU) Hamming-1 neighborhood shown in the figure:
     synonymous = 0 ✓, non-synonymous = 3.5–4.9 (not 1.2 or 8.7).
   - Action: figure uses the real computed distances (0, 3.5–4.9). Manuscript
     prose referencing 0/1.2/8.7 needs to be updated to match — or, if the
     prose intended illustrative values across any AA pair (not constrained
     to an Ala neighborhood), the "radical" example should be revised to a
     real radical pair such as G – P (6.23) or E – P (5.81).

2. **Methods claim "sensitivity <0.5 units".** False for both metrics once
   the random strategy is properly sampled with 3 seeds instead of 1.

   | Metric | Range with 1 random seed | Range with 3 random seeds (honest) |
   |:-:|:-:|:-:|
   | Noise D | 0.48 | **0.56** |
   | Dirichlet E | 0.56 | **0.92** |

   The E range was already above 0.5 with the original single-seed report;
   the manuscript never had a correct "<0.5" claim on E. With 3 seeds it is
   also false for D. Suggested Methods wording:

   > "Sensitivity to representative choice: z-score ranges across strategies
   > (closest to centroid, farthest from centroid, mean over three random
   > draws) are 0.56 units for noise distortion and 0.92 units for Dirichlet
   > energy — small compared to the SGC's ~15 standard-deviation displacement
   > from the null."

   Table S1 (`supplementary/table_S1.csv` and `.tex`) reports the numbers
   honestly.

3. **All z-scores in the manuscript text match the data.** Verified against
   `publication_controls.json`:
   - z_D(A) = −3.80 ✓, z_D(B) = −16.56 ✓, z_D(C) = −14.81 ✓
   - z_E(B) = −14.60 ✓, z_E(C) = −13.35 ✓
   - DE-correlation for Condition B: 0.959 ✓

## Layout and design decisions

1. **Figure 1 panel (a) — Hamming graph layout.** A 61-node graph with
   avg-degree 8.6 is inherently visually dense at any 2 × 2 inch panel size.
   Spring layout (seed = 7, k = 0.35, 300 iterations) gave the best balance
   of cluster separation and aspect ratio. Tried spectral and kamada_kawai:
   both produced more tangled layouts for this graph. Central nodes are
   unlabelled; five peripheral anchors (UUU Phe, GGG Gly, CAG Gln, AAA Lys,
   UCC Ser) provide compass orientation. AUG was dropped from anchors
   because the spring layout places it at the graph centroid (radius 0.14)
   where the label would collide with other nodes.

2. **Figure 2 layout.** Three panels in a row (double-column, 7.5 × 3.0 in).
   Slightly wider than the prompt's 7.2-inch cap so that panel (a)'s ylabel
   doesn't clip; reviewer can shrink to 7.0 in without loss.

3. **Figure 3 — overlay vs. stacked.** Chose **single-x-axis overlay** with
   alpha = 0.5 per histogram over stacked mini-histograms. Rationale: the B
   and C null distributions (triplet conditions) substantially overlap each
   other in D-space (both ≈ 28–34), which is the visual evidence that the
   factorial architecture effect (n = 2 → n = 3) dominates the alphabet
   effect (10 → 20 AAs). Stacked mini-histograms lose this co-location
   property.

   Dropped the (a)/(b)/(c) three-panel structure from the prompt; this is a
   two-panel figure (overlay + synonymy bars). The prompt explicitly allowed
   this: "drop panel (c) and make this a 2-panel figure".

   **Did not add Δz "architecture" and "alphabet" brackets** as requested —
   the previous version (`fig1_optimization_factorial.py`) had them; the new
   `fig3.py` omits them per the prompt's "DO NOT add Δz brackets claiming
   'architecture effect = −11.0 z-units'" directive. The factorial
   decomposition is now left to the prose.

4. **Figure 4 — substituted-base highlighting.** The prompt asked for
   underlined substituted bases. Matplotlib's mathtext backend does not
   support `\underline`, and Unicode combining low line (U+0332) renders
   blank in most matplotlib fontsets. **Used per-character coloring
   instead**: the substituted base is drawn bold, larger, and in the edge
   color (red for non-synonymous, blue for synonymous). This is at least as
   readable as an underline at the print size and matches the color legend
   of the figure. Noted in the figure's bottom text line.

5. **Color palette kept unchanged.** After the earlier pushback on the
   SGC-vs-20-AA-triplet color split, the existing Paul Tol palette was
   retained:
   - Red (#CC6677) — Condition B / SGC / 20-AA triplet
   - Blue (#4477AA) — Condition A / doublet / synonymous edge (Fig 4)
   - Green (#228833) — Condition C / 10-AA triplet
   - Grey (#BBBBBB) — null distributions
   - Paul Tol qualitative 10-group palette for the AA class encoding

## Deliverables

Generated:

- `figures/fig1.pdf`, `fig1.png`, `fig1_caption.txt`
- `figures/fig2.pdf`, `fig2.png`, `fig2_caption.txt`
- `figures/fig3.pdf`, `fig3.png`, `fig3_caption.txt`
- `figures/fig4.pdf`, `fig4.png`, `fig4_caption.txt`
- `supplementary/table_S1.csv`
- `supplementary/table_S1.tex`
- `supplementary/table_S1_summary.json`
- `run_sensitivity_extra_seeds.py` — script that generates the extra random
  seeds (1 and 2; seed 42 was already in `publication_controls.json`)
- `results/sensitivity_extra_seeds.json` — raw output of that run
- `build_table_s1.py` — assembles the table from both sources

The existing non-manuscript-numbered scripts remain in place:

- `figures/fig1_optimization_factorial.py`, `fig1b_pca_space.py`,
  `fig2_sgc_optimization.py`, `fig2_wobble_mechanism.py`,
  `fig3_factorial_synonymy.py`, `fig4_wobble_mechanism.py`,
  `fig_joint_scatter.py`, `figS1_raw_distributions.py`

These can be deleted once the manuscript is frozen; they are not referenced
by the manuscript figures but may be useful for exploratory re-plots.

## Items left alone but worth flagging

1. **Journal confirmation.** Repo and filenames say JME (Journal of Molecular
   Evolution). An earlier prompt in this workflow referenced MBE. All new
   output uses JME conventions. If the target is MBE instead, captions are
   generic enough to carry over; no reformatting needed besides column-width
   adjustment of Table S1's LaTeX.

2. **Arial vs. Liberation Sans.** Arial is not installed on boron. Used
   Liberation Sans which is metric-compatible (same glyph widths at the same
   point size). Publishers typically accept either. If the submission
   system demands Arial, the PDFs can be re-rendered on a system with Arial
   without changing any layout.

3. **Paudyal Comments 72 / 83** (whether Dirichlet energy and noise
   distortion are prior art). Cannot answer from the code or data — this is
   a literature question. A focused search of Koonin–Novozhilov 2017, Freeland
   & Hurst 1998, Archetti 2004 et al. for graph-Laplacian-style metrics on
   codon graphs would resolve it. Do not fabricate citations.

4. **Supplementary figure S1** (raw distributions per condition) was **not**
   rebuilt — the existing `figS1_raw_distributions.py` still renders cleanly
   and uses the same color scheme, so it can be kept as-is. Re-run if
   needed:

   ```bash
   python3 figures/figS1_raw_distributions.py
   ```

## Pre-emptive defense computations (2026-04-22)

Ran `run_preempt_metrics.py` at 100 000 shuffles per condition to pre-compute
two additional metrics we might be asked about in review:

### f⁺ L² f (squared graph Laplacian; Tlusty 2007 Eq. 3 analog)

| condition | f⁺ L f (manuscript) | f⁺ L² f |
|:-:|:-:|:-:|
| A: 10 AA, n=2 | −3.31 | −3.31 |
| C: 10 AA, n=3 | −13.30 | −12.24 |
| B: 20 AA, n=3 (SGC) | −14.59 | −13.52 |

Null-distribution correlation r(f⁺ L f, f⁺ L² f) = 0.987 across all three
conditions. The two metrics are nearly redundant on the Hamming-1 codon
graph, so the choice between them does not change any qualitative
conclusion. Keep as a one-line footnote if a reviewer asks.

### Boltzmann-weighted D (thermo_noise.py; K=3.0, wobble=0.4, RT=1.0)

| condition | D uniform (manuscript) | D Boltzmann |
|:-:|:-:|:-:|
| A: 10 AA, n=2 | −3.79 | −3.79 |
| C: 10 AA, n=3 | −14.76 | −20.12 |
| B: 20 AA, n=3 (SGC) | −16.56 | −24.89 |

Null-distribution correlation r(D_uniform, D_boltz) = 0.62. These are
*not* the same object — the Boltzmann weighting redirects probability
mass onto pos-3 errors, which are ~69 % synonymous in the SGC. The
realistic model therefore *strengthens* the SGC's z-score by 5–8 units
on the triplet conditions. For n = 2 (no distinct wobble position) the
two coincide.

**Framing for the manuscript:** the uniform-weighted D we report is the
conservative choice. A reviewer demanding a mechanistic weighting should
take the Boltzmann-weighted result (z = −24.9 for the SGC at 10⁵
samples) as *stronger* evidence, not weaker. Adding a brief Methods
paragraph that says this inoculates against either critique.

Raw numerical output: `results/preempt_metrics.json`.
Script: `run_preempt_metrics.py`.

# Codon-position-resolved translational misreading rates: extraction from primary sources

Date: 2026-09-23. Scope: data extraction only. No fitting or optimisation.
Machine-readable version: `POSITION_RATES.json` (built by `scripts/build_position_rates_json.py`).

**Definitions.** ε_k is misreading caused by a single codon–anticodon mismatch at
**mRNA codon position k**. Position 1 pairs with anticodon 36; position 3 is the
wobble position and pairs with anticodon 34. All ratios are normalised to ε₂ = 1.
Every assay below sees only **amino-acid-changing** misreading. Synonymous
misreading, which is mostly at position 3, is invisible to all of them, so
ε₃ means "non-synonymous misreading per non-synonymous position-3 opportunity".

## Main table

| Source | Organism / system | ε₁ : ε₂ : ε₃ (ε₂ = 1) | 95% interval on ratio | n | How computed | Caveats |
|---|---|---|---|---|---|---|
| Katoh & Suga 2023 NAR, Supp Fig 8 (strict) | E. coli, reconstituted in-vitro translation, cognate aa withdrawn | **0.98 : 1 : 1.39** | ε₁/ε₂ 0.69–1.40; ε₃/ε₂ 0.95–2.05 (Katz) | single-mismatch combinations tested 64 / 71 / 21; with misincorporation 30 / 34 / 14 | per-opportunity frequency = (combinations with intensity > 0) / (combinations tested), restricted to codon–aa combinations with exactly one '+' mismatch row and no '±' | no cognate competitor, so this measures what near-cognates *can* decode, not in-vivo rates. Intensity is a 0–5 semi-quantitative score. Position-3 n is small. Position-3 "mismatch" is judged after modified-wobble rules |
| same, authors' "1 mismatch" class | same | 0.94 : 1 : 1.39 | ε₁/ε₂ 0.65–1.36; ε₃/ε₂ 0.95–2.05 | 62 / 71 / 21; positive 28 / 34 / 14 | as above, but using the authors' `total == 1` class | differs from strict only by UCC→Ala and UCU→Ala, which the figure labels total "2" despite a +/−/− pattern (see inconsistencies) |
| same, weighted by intensity | same | 0.99 : 1 : 1.79 | — | as strict | mean intensity (0–5, zeros included) per tested combination | intensity is ordinal, so treat as indicative only |
| Landerer et al. 2024 MBE, Data S4 | **E. coli**, in vivo, 80 public proteomes (open-search MS) | **1.37 : 1 : 0.71** | dataset bootstrap: ε₁/ε₂ 1.10–2.55; ε₃/ε₂ 0.32–2.25 | single-mismatch PSMs 13,029 / 10,505 / 2,157 (of 58,906) | our assignment: an event (codon c → aa d) goes to position k if a codon of d differs from c only at k. ε_k = N_k / Σ_c B_c·m_k(c), with B_c = codon PSM exposure and m_k = detectable non-synonymous single changes at k (stop and I↔L excluded) | one dataset (PXD014341) holds 47% of single-mismatch PSMs. **Leaving it out gives 1.91 : 1 : 1.38.** 28% of PSMs need ≥2 mismatches and 29% are position-ambiguous (1/2 or 1/3), so both are excluded. Cys-origin events dominate (C→S, C→A) and may be artefacts |
| Landerer et al. 2024, Data S4 | **S. cerevisiae**, 72 public proteomes | **1.40 : 1 : 0.71** | dataset bootstrap: ε₁/ε₂ 0.76–2.54; ε₃/ε₂ 0.35–1.19 | single-mismatch PSMs 3,106 / 2,470 / 636 (of 17,765) | as above | 59% of yeast PSMs need ≥2 mismatches, so most events are not near-cognate. Leaving out the largest dataset (PXD026077, 11%) gives 1.16 : 1 : 0.57. Position 1 is dominated by GCT→Pro (1,348 of 3,106) |
| Landerer et al. 2024, Fig. S5c,d | E. coli, yeast | **figure only** | — | — | not digitised | authors' text: "Misincorporations caused by a single mismatch are most frequently observed at the second codon position". This is a share of events, not per opportunity |
| Kramer & Farabaugh 2007 RNA, Table 1 | E. coli, in vivo, luciferase K529 reporter, tRNA-Lys(UUU) only | **not computable** (Table 1 partial; full set is figure only) | — | 14 near-cognate codons tested (6 at position 1, 6 at position 2, 2 at position 3) | Table 1 lists only the codons above background (below) | one tRNA, one site. The sense position-1 near-cognates (CAA, CAG, GAA, GAG) appear only in the Fig. 2 bar chart |
| Mordret et al. 2019 Mol Cell | E. coli, yeast | **not accessed; figure only** | — | — | — | see "Could not access" |
| Manickam et al. 2016 NAR | E. coli, in vivo reporter | **figure only** | — | — | — | rates appear only in bar-chart figures |

### Share of detected single-mismatch events by position (raw, not per opportunity)

| Source | pos 1 | pos 2 | pos 3 | Unit |
|---|---|---|---|---|
| Katoh 2023 (strict) | 38.5% | 43.6% | 17.9% | 78 positive combinations |
| Landerer E. coli | 50.7% | 40.9% | 8.4% | 25,691 PSMs |
| Landerer E. coli | 51.3% | 41.7% | 7.0% | 6,559 unique sites |
| Landerer yeast | 50.0% | 39.8% | 10.2% | 6,212 PSMs |
| Landerer yeast | 49.3% | 43.5% | 7.1% | 1,923 unique sites |

The low position-3 *share* comes mostly from the code itself: most position-3
changes are synonymous. Per opportunity (main table), position 3 is not
depleted by a comparable factor.

### Kramer & Farabaugh 2007, Table 1, WT column (misreading frequency ×10⁻⁴, as printed)

Transcribed from the typeset table image `raw/kramer2007/87tbl1.jpg`
(sha256 `10d19b26…58d0`), column "WT":

| Codon | aa | Position misread | WT ×10⁻⁴ |
|---|---|---|---|
| UAA | Ter | 1 | 4.1 ± 0.51 |
| UAG | Ter | 1 | 14 ± 0.66 |
| AUA | Ile | 2 | 3.5 ± 0.44 |
| AGA | Arg | 2 | 36 ± 2.7 |
| AGG | Arg | 2 | 31 ± 1.5 |
| AAU | Asn | 3 | 16 ± 0.78 |
| AAC | Asn | 3 | 3.8 ± 0.52 |
| Others | — | — | 3.1 ± 0.74 (footnote: "Average of all remaining tested codons") |

The text says the median over the 24 codons tested was "no more than 3.4 × 10⁻⁴"
and GAA (position 1) was "no more than 2.0 × 10⁻⁴". "Others" is at the
UUU-control floor. A mean ε per position would require assigning values to the
unlisted codons, so we did not compute one.

## Provenance per source

**Katoh T, Suga H. 2023.** *Nucleic Acids Res* 51(19):10642–10652.
doi:10.1093/nar/gkad707, PMC10602915.
- Files: `raw/katoh2023_suppl.zip` (Europe PMC) → `raw/katoh2023/suppl/Supplementary Figures.pdf`; full-text XML `raw/katoh2023_fulltext.xml`.
- Table used: Supplementary Figure 8, pdf pages 12–15. The panels are NUN, NCN, NAN, NGN. For each of 61 codons × 20 aa there are rows `intensity`, `1st`, `2nd`, `3rd`, `total`.
- `scripts/parse_katoh_s8.py` read each cell as the literal vector-text glyph (PyMuPDF word boxes), not as pixels. Output: `katoh2023_s8_parsed.csv`, 1,159 rows.
- **Validation:** the parse reproduces main-text Figure 4 exactly. Class sizes 1 / 2 / 3 / 1/2 / 2/3 / 1/2/3 = 154 / 348 / 211 / 160 / 244 / 42 (sum 1,159). Misdecoding frequencies are 49.4 / 18.7 / 14.7 / 65.6 / 29.5 / 40.5%. The per-intensity counts in every bar also match.
- Important: Figure 4's classes are the **number** of mismatches (1, 2, 3), not their position. The 49.4 / 18.7 / 14.7% figures must not be read as ε₁ / ε₂ / ε₃. Position comes only from the Supp Fig 8 `1st/2nd/3rd` rows.
- **Source inconsistencies** (printed total disagrees with the printed +/− rows): AGC→Q and AGU→Q (+/+/+, total "2"); UCC→A and UCU→A (+/−/−, total "2"). Another 42 rows differ only because '±' is correlated across isoacceptors, which is expected.
- Computation: `scripts/katoh_position_rates.py` → `katoh2023_position_summary.json`. This file also has the per-position split by codon second base.

**Landerer C, Poehls J, Toth-Petroczy A. 2024.** *Mol Biol Evol* 41(3):msae048.
doi:10.1093/molbev/msae048, PMC10939442.
- Files: `raw/landerer2024/` (Europe PMC supplementary zip; Data S1–S4, supplementary PDF, full-text XML).
- Table used: Data S4 `{ecoli,yeast}_filtered/*_substitution_errors.csv` and `*_codon_counts.csv`.
- Computation: `scripts/landerer_position_rates.py` → `landerer2024_position_summary.json`; `scripts/landerer_dataset_bootstrap.py` → `landerer2024_bootstrap.json` (2,000 dataset resamples, seed 0).
- **The position assignment is ours**, from codon sequence (Hamming-1 near-cognate). It is not a number the authors report. It does not consult actual tRNA anticodons or modifications.
- Sensitivities, E. coli / yeast ε₁:ε₂:ε₃:
  - Cys-origin events and Cys codons removed: 1.22 : 1 : 0.71 / 1.39 : 1 : 0.71.
  - Opportunities limited to (origin, destination) pairs observed at least once: 1.53 : 1 : 0.78 / 1.52 : 1 : 0.61.

**Kramer EB, Farabaugh PJ. 2007.** *RNA* 13(1):87–96. doi:10.1261/rna.294907, PMC1705757.
- Files: `raw/kramer2007/kramer2007_pmc.html` (PMC article page), `raw/kramer2007/87tbl1.jpg` (Table 1 image from the PMC CDN).
- Table used: Table 1, WT column. Per-codon values for all 14 near-cognates are in Figure 2 (bar chart), which was not digitised.

**Mordret E et al. 2019.** *Mol Cell* 75(3):427–441.e5. doi:10.1016/j.molcel.2019.06.041, PMID 31353208.
- Could not access. See below.

**Manickam N, Joshi K, Bhatt MJ, Farabaugh PJ. 2016.** *Nucleic Acids Res* 44(4):1871.
doi:10.1093/nar/gkv1506, PMC4770228.
- Files: `raw/manickam2016/`.
- Figure only. Checked because it extends the Kramer reporter to more tRNAs. It was not one of the four priority sources.

Checksums of every saved raw file: `raw/SHA256SUMS.txt` and `raw/landerer2024/S4_SHA256SUMS.txt`.

## Could not access

- **Mordret et al. 2019 (priority a).** The paper is marked bronze open access (Unpaywall), but:
  - cell.com and sciencedirect.com returned HTTP 403 (Cloudflare challenge) to curl, WebFetch and a headless Playwright browser;
  - bioRxiv preprint 10.1101/255943 returned 429;
  - there is no PMC copy;
  - github.com/ernestmordret/substitutions contains pipeline code only, no result tables.

  Full-text excerpts via scite show that position data are in **Figure 3A**: three 4×4 "confusion matrices" per organism, "shade intensity reflects the logarithm (base 2) of the number of peptides". That is **figure only**. No numbers were taken from it. Quoted qualitative statements (scite excerpts): "The most frequently observed substitution types involve mismatches between uracil and guanine in the first or the second position of the codon"; for yeast, "predominantly feature G:U mismatches at the first or second positions ... elevated level of G:G mispairing in the second position". The supplementary tables should be downloadable by hand in an ordinary browser session.
- **Landerer Fig. S5c,d.** Accessible, but figure only.

## Caveats that apply across sources

1. **Different quantities.** Katoh is a per-combination *capability* in vitro with the cognate tRNA absent. Landerer is in-vivo *detected PSMs*, which depend on MS detectability. Kramer is an in-vivo per-codon *frequency* for one tRNA. They are not interchangeable.
2. **Only non-synonymous misreading is observable.** ε₃ is conditional on the position-3 change altering the amino acid.
3. **Mismatch position vs codon difference.** Katoh assigns position from actual anticodon pairing, including modified wobble bases. Our Landerer assignment uses codon Hamming distance. A "position-3" event in Landerer may be a wobble-tolerated pairing in reality, and the reverse.
4. **Concentration.** A handful of codon→aa pairs carry most events at each position (listed under `top_pairs_by_position` in the JSON). The E. coli ratio depends strongly on one dataset.
5. **Artefacts.** Landerer removed known-PTM mass matches, but the most frequent substitutions (C→S, C→A, and in yeast A→P) are candidate sample-preparation artefacts. Removing Cys changes ε₁/ε₂ in E. coli from 1.37 to 1.22.

## Plain reading of the numbers

- No source supports a strongly suppressed ε₂ relative to ε₁. ε₁/ε₂ is 0.94–0.98 in vitro (Katoh) and 1.2–1.5 in vivo (Landerer), with intervals from about 0.65 to 2.5.
- No source supports ε₃ ≫ ε₂ per non-synonymous opportunity. ε₃/ε₂ is 1.39 in Katoh (95% 0.95–2.05) and about 0.71 in Landerer. The E. coli bootstrap range is 0.32–2.25, and leaving out one dataset turns it into 1.38.
- The position-3 *share* of detected events (7–18%) is low mainly because most position-3 changes are synonymous.
- None of these sources pins ε₃/ε₂ to better than roughly a factor of 2.

# INVENTORY — Landerer, Poehls & Toth-Petroczy 2024, Data S2

File inspected: `/storage/kiran-stuff/codon-deployment/data/raw/Data_S2_error_detection_rate.xlsx`
(copied unchanged to `raw/Data_S2_error_detection_rate.xlsx`).

sha256 `08ef9455c90b937bf622a93a459e4106b65484133a9a34d280df627396dc266b`. This is
byte-identical to `Data_S2_error_detection_rate.xlsx` inside the journal's own
supplementary zip (`raw/landerer2024/msae048_supplementary_data.zip`, downloaded
from Europe PMC, PMC10939442). So the local copy is the published file.

Source: Landerer C, Poehls J, Toth-Petroczy A. 2024. *Mol Biol Evol* 41(3):msae048.
doi:10.1093/molbev/msae048.

## Structure

The workbook has no merged cells, no cell comments and no defined names.

| Sheet | Range | Data rows | Columns |
|---|---|---|---|
| `S. cerevisiae` | A1:J62 | 61 (one per sense codon, DNA alphabet AAA…TTT) | 10 |
| `E. coli` | A1:J62 | 61 | 10 |

Both sheets have the same columns (header row 1, verbatim):

| Column | Header | Meaning (from the paper and the S4 files) |
|---|---|---|
| A | `Codon` | source (mRNA) sense codon, T not U; stop codons absent |
| B | `mean` | mean over PRIDE datasets of the per-codon detection rate = `error_count / base_count` (see Data S4 `*_codon_counts.csv`) |
| C | `median` | median of the same |
| D | `sd` | standard deviation |
| E | `var` | variance (displays as 0.0 at default precision; the true values run from 1.6e-10 to 3.4e-3) |
| F | `range` | max − min |
| G | `iqr` | interquartile range |
| H | `se` | standard error |
| I | `cv` | coefficient of variation |
| J | `n_unique_aa_misincroporations` (typo in source) | number of distinct destination amino acids seen from that codon (yeast 1–14, E. coli 1–18) |

Value ranges of `mean`: yeast 1.0e-5 to 4.7e-3; E. coli 3.3e-5 to 2.0e-2. Neither sheet
has any missing cells.

## Does Data S2 resolve errors by mismatch position or by near-cognate codon?

**No.** Each row is one *source* codon with summary statistics of its total
detection rate, pooled over all destination amino acids. The table has no
destination amino acid, no near-cognate codon, no mismatch position and no
mismatch type. Codon-position rates **cannot** be derived from Data S2.

(`/storage/kiran-stuff/codon-deployment/data/raw/make_codon_error_rates_from_landerer.py`
reads only `Codon` and `mean` from the `E. coli` sheet. It was only read, not changed.)

## The other Landerer supplementary files (downloaded to `raw/landerer2024/`)

| File | Content | Position-resolved? |
|---|---|---|
| `Data_S1_datasets.xlsx` | 6 sheets (`S.cerevisiae_all/_curated/_used`, `E.coli_all/_curated/_used`): PRIDE IDs, n misincorporations, n peptides, curation flags | no |
| `Data_S3_ptm_masking.xlsx` | 130 rows: PTMs/artefacts (`delta_m`, `site`, `position`, `title`) whose mass mimics a substitution. These were removed upstream | no |
| `Data_S4_detected_misincorporations.zip` | `ecoli_filtered/` (80 datasets) and `yeast_filtered/` (72): per dataset, `*_substitution_errors.csv` (one row per PSM: protein, site, **source codon**, origin aa, **destination aa**), `*_codon_counts.csv` (codon, base_count, error_count, detection_rate) and `*_peptide_counts.csv` | **yes, indirectly.** Source codon + destination aa lets each event be assigned to the mismatched codon position. See POSITION_RATES.md |
| `Landerer_MBE_supplementary_material_accepted.pdf` | methods; Fig. S5c,d = "Percentage of observed amino acid misincorporation caused by a given codon position" | Fig. S5 is a figure only (not digitised) |

Notes on Data S4:
- In every dataset, the `error_count` summed over `*_codon_counts.csv` equals the number of rows in `*_substitution_errors.csv` (checked by assertion), so detection rates are in PSM units.
- The `Data_S1` `*_used` sheets list exactly the 80 E. coli / 72 yeast datasets present in S4. The other rows in those sheets are `sum`/`min`/`max` summary rows.
- 84 E. coli PSM rows (0.14%) give a codon that does not encode the stated origin amino acid (73 of them `ATC`, origin `A`). These were dropped from our analysis. Yeast has none.
- E. coli `codon_counts` contains one `TGA` occurrence (base_count 1, error_count 0). It was excluded from exposure.
- The macOS `__MACOSX/` resource-fork folder from the zip was deleted. Checksums of the S4 files are in `raw/landerer2024/S4_SHA256SUMS.txt`.

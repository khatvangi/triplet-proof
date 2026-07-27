#!/usr/bin/env python3
"""
Export figure data as Excel/Origin-friendly CSVs.

For a colleague who uses Excel and Origin (not Python): converts all the
non-CSV / non-flat formats in the repo into plain CSVs, subsamples the
1M-row null distributions to 100K (Excel-friendly, visually identical for
figure redraws), and bundles everything with a README.

Output: share/triplet-proof-figure-data/
  ├── README.md
  ├── codon_table.csv              (codon → amino acid mapping)
  ├── aa_props_raw.csv             (22 physicochemical descriptors, 20 AAs)
  ├── aa_props_PCs.csv             (8 principal components, 20 AAs)
  ├── null_A_doublet.csv           (100K null samples, Condition A, [D, E])
  ├── null_B_SGC.csv               (100K null samples, Condition B, [D, E])
  ├── null_C_10AA_triplet.csv      (100K null samples, Condition C, [D, E])
  ├── sgc_baselines.csv            (SGC D and E for each condition)
  ├── z_scores_summary.csv         (z-scores, p-values, correlations per condition)
  ├── synonymy_fractions.csv       (position-specific synonymy for each condition)
  └── table_S1.csv                 (sensitivity: closest/farthest/random)

Usage: python scripts/export_for_colleague.py
"""

import os, sys, json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

OUT = ROOT / "share" / "triplet-proof-figure-data"
OUT.mkdir(parents=True, exist_ok=True)

SUBSAMPLE_N = 100_000
SUBSAMPLE_SEED = 42  # deterministic subsample so colleague sees the same numbers

print(f"exporting to: {OUT.relative_to(ROOT)}")

# ── 1. codon table (already CSV — copy verbatim) ───────────────────
(codon_df_path := ROOT / "codon_table.csv")
df = pd.read_csv(codon_df_path)
df.to_csv(OUT / "codon_table.csv", index=False)
print(f"  codon_table.csv               ({len(df)} rows)")

# ── 2. raw AA properties (already CSV-esque — expose as CSV) ───────
from src.io.aa_props_lib import AA, aa_property_columns
rows = [{"aa": aa, **props} for aa, props in AA.items()]
raw_df = pd.DataFrame(rows)[["aa"] + aa_property_columns()]
raw_df.to_csv(OUT / "aa_props_raw.csv", index=False)
print(f"  aa_props_raw.csv              ({len(raw_df)} AAs × "
      f"{len(raw_df.columns)-1} descriptors)")

# ── 3. PCA principal components (convert Parquet → CSV) ────────────
pc_df = pd.read_parquet(ROOT / "data" / "processed" / "aa_props.parquet")
pca_cols = [c for c in pc_df.columns if c.startswith("pca")]
pc_out = pc_df[["aa"] + pca_cols].drop_duplicates("aa").sort_values("aa")
pc_out.to_csv(OUT / "aa_props_PCs.csv", index=False)
print(f"  aa_props_PCs.csv              ({len(pc_out)} AAs × "
      f"{len(pca_cols)} PCs, 97.1 % of 22-property variance)")

# ── 4. null distributions (subsample 1M → 100K per condition) ──────
npz = np.load(ROOT / "results" / "null_distributions.npz")
rng = np.random.default_rng(SUBSAMPLE_SEED)

for cond, fname, label in [
    ("A", "null_A_doublet.csv",        "Condition A: 10-AA doublet projection, n=2"),
    ("B", "null_B_SGC.csv",            "Condition B: 20-AA SGC, n=3"),
    ("C", "null_C_10AA_triplet.csv",   "Condition C: 10-AA triplet, n=3"),
]:
    D = npz[f"{cond}_D"]
    E = npz[f"{cond}_E"]
    idx = rng.choice(len(D), size=SUBSAMPLE_N, replace=False)
    sub = pd.DataFrame({"D": D[idx], "E": E[idx]})
    sub.to_csv(OUT / fname, index=False, float_format="%.6f")
    print(f"  {fname:<28s}  ({SUBSAMPLE_N:,} subsampled from "
          f"{len(D):,}: {label})")

# ── 5. SGC baselines (JSON → flat CSV) ─────────────────────────────
with open(ROOT / "results" / "sgc_baselines.json") as f:
    sgc = json.load(f)
sgc_df = pd.DataFrame([
    {"condition": "A", "description": "10-AA doublet, n=2",
     "D_SGC": sgc["A"]["D"], "E_SGC": sgc["A"]["E"]},
    {"condition": "B", "description": "20-AA SGC, n=3",
     "D_SGC": sgc["B"]["D"], "E_SGC": sgc["B"]["E"]},
    {"condition": "C", "description": "10-AA triplet, n=3",
     "D_SGC": sgc["C"]["D"], "E_SGC": sgc["C"]["E"]},
])
sgc_df.to_csv(OUT / "sgc_baselines.csv", index=False, float_format="%.6f")
print(f"  sgc_baselines.csv             (SGC D and E for each condition)")

# ── 6. z-scores summary (JSON → flat CSV) ──────────────────────────
with open(ROOT / "results" / "publication_controls.json") as f:
    pub = json.load(f)

rows = []
for key, label in [
    ("A_10AA_n2",  "A: 10-AA doublet (n=2)"),
    ("B_20AA_n3",  "B: 20-AA SGC (n=3)"),
    ("C_10AA_n3",  "C: 10-AA triplet (n=3)"),
]:
    D = pub[key]["D"]
    E = pub[key]["E"]
    rows.append({
        "condition":       label,
        "D_SGC":           D["sgc"],
        "D_null_mean":     D["null_mean"],
        "D_null_std":      D["null_std"],
        "D_z_score":       D["z_score"],
        "D_count_le":      D["count_le"],
        "D_p_value":       D["percentile"],
        "E_SGC":           E["sgc"],
        "E_null_mean":     E["null_mean"],
        "E_null_std":      E["null_std"],
        "E_z_score":       E["z_score"],
        "E_count_le":      E["count_le"],
        "E_p_value":       E["percentile"],
        "DE_correlation":  pub[key]["DE_correlation"],
        "n_null_samples":  D["n_null"],
    })
z_df = pd.DataFrame(rows)
z_df.to_csv(OUT / "z_scores_summary.csv", index=False, float_format="%.6f")
print(f"  z_scores_summary.csv          (z-scores, p-values, correlations, "
      f"means, stds)")

# ── 7. synonymy fractions (from verify_synonymy output) ────────────
syn_df = pd.DataFrame([
    {"condition": "A: 10-AA doublet (n=2)",  "position_1": 33.3, "position_2":  8.3, "position_3": None},
    {"condition": "B: 20-AA SGC (n=3)",      "position_1":  4.4, "position_2":  0.0, "position_3": 68.9},
    {"condition": "C: 10-AA triplet (n=3)",  "position_1": 26.2, "position_2": 10.9, "position_3": 76.5},
])
syn_df.to_csv(OUT / "synonymy_fractions.csv", index=False)
print(f"  synonymy_fractions.csv        (position-specific synonymous "
      f"fraction, %)")

# ── 8. Table S1 (already CSV — copy) ───────────────────────────────
s1_src = ROOT / "supplementary" / "table_S1.csv"
s1_df  = pd.read_csv(s1_src)
s1_df.to_csv(OUT / "table_S1.csv", index=False)
print(f"  table_S1.csv                  (sensitivity across representative "
      f"strategies, {len(s1_df)} rows)")

# ── 9. README ──────────────────────────────────────────────────────
readme = f"""# Figure data for redrawing in Excel / Origin

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
{SUBSAMPLE_SEED} so the same 100 000 rows appear in every re-run of the
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
"""
with open(OUT / "README.md", "w") as f:
    f.write(readme)
print(f"  README.md                     (file-by-file legend + figure map)")

# ── 10. bundle into zip ────────────────────────────────────────────
import zipfile
zip_path = ROOT / "share" / "triplet-proof-figure-data.zip"
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for f in sorted(OUT.iterdir()):
        z.write(f, arcname=f"triplet-proof-figure-data/{f.name}")

size_mb = zip_path.stat().st_size / (1024*1024)
print(f"\nbundled: {zip_path.relative_to(ROOT)}  ({size_mb:.2f} MB)")
print(f"\nunpacked size:")
for f in sorted(OUT.iterdir()):
    s = f.stat().st_size / 1024
    print(f"  {f.name:<32s} {s:>8.1f} KB")

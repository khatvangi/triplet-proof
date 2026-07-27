# Supplementary Tables for Journal of Molecular Evolution Submission

## Table S1. Sensitivity of Condition C (10 amino acids on triplet codons) to representative-set choice

Condition C was re-evaluated using three representative-selection rules for each reduced-alphabet (RAA10) class: closest to centroid, farthest from centroid, and random class member. Each entry is computed against its own 100 000-shuffle null distribution. For the random strategy, three independent seeds (1, 2, 42) were run and the reported value is the mean across seeds; SD is the sample standard deviation across the three seeds. Values below match `supplementary/table_S1.csv` and `supplementary/table_S1.tex`.

| Representative set | Distortion z (D) | Dirichlet z (E) | D–E correlation |
|---|---:|---:|---:|
| Closest to centroid (main text) | −14.79 | −13.34 | 0.958 |
| Farthest from centroid | −14.61 | −13.18 | 0.965 |
| Random class member (mean ± SD, 3 seeds) | −15.17 ± 0.16 | −14.09 ± 0.54 | 0.961 ± 0.008 |

Individual random-seed values:

- D z-scores: {seed 42: −15.09, seed 1: −15.36, seed 2: −15.07}
- E z-scores: {seed 42: −13.74, seed 1: −14.71, seed 2: −13.83}
- D–E correlations: {seed 42: 0.955, seed 1: 0.969, seed 2: 0.958}

Ranges across representative strategies (using the random-strategy mean for the random row):

- Distortion z-score range: **0.56**
- Dirichlet z-score range: **0.92**

Note on the "<0.5 units" wording: an earlier draft of the Methods stated that the sensitivity z-score range was below 0.5 units. That claim held only under a single random-seed comparison. With three random seeds, the random-strategy mean shifts, and both ranges exceed 0.5. The values above are the honest numbers; recommended Methods wording is "z-score ranges are 0.56 (D) and 0.92 (E), small compared to the ~15 SD displacement of the SGC from its null".

Sources: `results/publication_controls.json` (closest, farthest, random seed 42) and `results/sensitivity_extra_seeds.json` (random seeds 1 and 2). See `build_table_s1.py` for the aggregation script and `supplementary/table_S1_summary.json` for the machine-readable numeric summary.


## Table S2. Physicochemical descriptor set used for PCA

Amino acids were represented in an 8-dimensional physicochemical space derived by principal component analysis (PCA) of a curated descriptor library. The library begins with **22 descriptors** spanning hydrophobicity, size/volume, polarity and charge, hydrogen-bonding capacity, aromaticity, conformational (secondary-structure) propensity, flexibility, and solvent accessibility. Descriptors were z-scored, then **two were removed by correlation pruning** (Pearson |r| > 0.95 with a retained descriptor), leaving **20 retained descriptors**. PCA on the retained set yields **8 principal components explaining 97.10% of total variance**; each amino acid is represented by its 8-dimensional PC score vector (Fig. 1b).

Descriptor sources follow canonical published scales, including AAindex-curated measures (Kawashima and Kanehisa 2000). The full 20 × 22 numeric descriptor matrix is provided in machine-readable form as `supplementary/table_S2.csv`.

| # | Descriptor code | Property | Source / scale | Status |
|---:|---|---|---|---|
| 1 | `hydro_kd` | Hydropathy | Kyte–Doolittle hydropathy index | Retained |
| 2 | `vol_zam` | Residue volume (Å³) | Zamyatnin residue volumes | Retained |
| 3 | `polarity_gr` | Polarity | Grantham polarity | Retained |
| 4 | `pI` | Isoelectric point | Standard pI | Retained |
| 5 | `sc_pka` | Side-chain pKa | Side-chain pKa (0 if none) | Retained |
| 6 | `charge_pH7` | Net charge at pH 7 | Formal side-chain charge (His ≈ +0.1) | Retained |
| 7 | `donors` | H-bond donors | Side-chain donor count | Retained |
| 8 | `acceptors` | H-bond acceptors | Side-chain acceptor count | Retained |
| 9 | `arom` | Aromaticity | Aromatic side chain (0/1) | Retained |
| 10 | `ring_count` | Ring count | Number of ring systems in side chain | Retained |
| 11 | `rotb` | Rotatable bonds | Side-chain rotatable-bond count | Retained |
| 12 | `refractivity` | Molar refractivity | Amino-acid molar refractivity | Pruned (redundant with vol_zam (r=0.99)) |
| 13 | `bulkiness` | Bulkiness | Zimmerman bulkiness | Retained |
| 14 | `cf_alpha` | α-helix propensity | Chou–Fasman P(α) | Retained |
| 15 | `cf_beta` | β-sheet propensity | Chou–Fasman P(β) | Retained |
| 16 | `cf_turn` | Turn propensity | Chou–Fasman P(turn) | Retained |
| 17 | `acc_pref` | Solvent accessibility | Relative accessibility preference | Retained |
| 18 | `flex` | Flexibility | Normalized B-factor flexibility index | Retained |
| 19 | `helical_moment` | Helical moment | Helix hydrophobic-moment proxy | Retained |
| 20 | `beta_moment` | β-strand moment | β-strand moment proxy | Retained |
| 21 | `hphob_fr` | Hydrophobic fraction | Hydrophobic surface-tendency proxy | Pruned (redundant with helical_moment (r=0.97)) |
| 22 | `sc_mass` | Side-chain mass (Da) | Approx. side-chain mass | Retained |

*Two descriptors were dropped in correlation pruning: `refractivity` (r = 0.99 with `vol_zam`) and `hphob_fr` (r = 0.97 with `helical_moment`; r = 0.97 with `acc_pref`). Removing them eliminates redundant variance without discarding a distinct physicochemical axis.*

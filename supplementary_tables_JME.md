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

The local property library begins with 22 amino acid descriptors and retains 20 after correlation pruning. The retained standardized descriptors are listed below. PCA on this retained set yields 8 principal components explaining 97.1% of total variance.

Retained descriptors:

1. `hydro_kd_z`
2. `vol_zam_z`
3. `polarity_gr_z`
4. `pI_z`
5. `sc_pka_z`
6. `charge_pH7_z`
7. `donors_z`
8. `acceptors_z`
9. `arom_z`
10. `ring_count_z`
11. `rotb_z`
12. `bulkiness_z`
13. `cf_alpha_z`
14. `cf_beta_z`
15. `cf_turn_z`
16. `acc_pref_z`
17. `flex_z`
18. `helical_moment_z`
19. `beta_moment_z`
20. `sc_mass_z`

Original descriptor pool before pruning:

1. `hydro_kd`
2. `vol_zam`
3. `polarity_gr`
4. `pI`
5. `sc_pka`
6. `charge_pH7`
7. `donors`
8. `acceptors`
9. `arom`
10. `ring_count`
11. `rotb`
12. `refractivity`
13. `bulkiness`
14. `cf_alpha`
15. `cf_beta`
16. `cf_turn`
17. `acc_pref`
18. `flex`
19. `helical_moment`
20. `beta_moment`
21. `hphob_fr`
22. `sc_mass`

PCA summary:

- Number of retained descriptors: 20
- Number of PCs retained: 8
- Variance explained: 0.9710

# Supplementary Tables for Journal of Molecular Evolution Submission

## Table S1. Sensitivity of Condition C (10 amino acids on triplet codons) to representative-set choice

Condition C was re-evaluated using three representative-selection rules for each reduced-alphabet class: closest to centroid, farthest from centroid, and random class member. Values below are taken from `results/publication_controls.json`.

| Representative set | Distortion z | Dirichlet z | D-E correlation |
|--------------------|-------------:|------------:|----------------:|
| Closest to centroid | -14.87 | -13.41 | 0.958 |
| Farthest from centroid | -14.57 | -13.12 | 0.965 |
| Random class member | -15.08 | -13.73 | 0.955 |

Range across representative sets:

- Distortion z range: 0.51
- Dirichlet z range: 0.61

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

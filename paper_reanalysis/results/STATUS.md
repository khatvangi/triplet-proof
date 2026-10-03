# Reanalysis status

grid completed: 2026-09-23 14:40:54
report built:   2026-09-23 14:41:08

## z grid

|                                           |   ('D', 'block') |   ('D', 'codon') |   ('E', 'block') |   ('E', 'codon') |
|:------------------------------------------|-----------------:|-----------------:|-----------------:|-----------------:|
| 8 PCs, whitened (what the paper computed) |            -2.54 |           -16.56 |            -1.92 |           -14.6  |
| 8 PCs, unwhitened                         |            -3.01 |           -15.55 |            -2.17 |           -12.81 |
| 20 descriptors (pre-PCA)                  |            -2.97 |           -15.85 |            -2.19 |           -12.97 |
| polar requirement (1-D)                   |            -3.4  |           -11.59 |            -2.75 |            -7.85 |
| CHEMISTRY-BLIND (all AAs equidistant)     |             1    |           -17.07 |           nan    |           -17.07 |

## cross-validation vs the preregistered 10^7 block null

{
  "E": {
    "this_run": -1.9163700688569454,
    "block_null_D6_10e7": -1.9153,
    "abs_diff": 0.0010700688569453831
  },
  "D": {
    "this_run": -2.5385027843551877,
    "block_null_D6_10e7": -2.5369,
    "abs_diff": 0.0016027843551875698
  }
}

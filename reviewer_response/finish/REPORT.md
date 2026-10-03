# finish: R4 stop-codon variant and missing position_model P cells

## (1) R4: n=3 on 61 sense codons (stops excluded), uniform weights, k=10

| usage | best Dw (this run) | hits random / warm | usage_cert R4 | best over both | n=2 optimum | gap | gap % | syn p1 / p2 / p3 (raw) |
|:-:|--:|:-:|--:|--:|--:|--:|--:|:-:|
| U | 18.866276960 | 4/512 / 0/64 | 18.857417193 | 18.857417193 | 20.017516292 | 1.151239 | 5.75 | 0.000 / 0.352 / 1.000 |
| P | 17.383593506 | 5/512 / 0/64 | 17.383593506 | 17.383593506 | 18.835846661 | 1.452253 | 7.71 | 0.575 / 0.830 / 0.000 |

gap = n=2 optimum − best n=3 (61-sense) Dw; positive = the triplet code does better. synonymy = fraction of sense→sense arcs at that position whose endpoints share a class, in the best code from this run. Usage-weighted synonymy is in RESULTS.json.

caveats: this run uses usage_cert's kernels and schedule on disjoint seeds. For U it did NOT reproduce usage_cert's best (18.866 vs 18.857; hit rates 4/512 and 1/512), so the U optimum is not converged and 18.857 is only an upper bound. For P both runs landed on the same value. No warm start from the n=2 lift reached the best in either usage. The synonymy profile belongs to one best code. Other codes with equal Dw may spread synonymy across the positions differently, so do not read the p1-vs-p3 split as unique.

## (2) position_model, usage P, n=3 cells that were missing

missing before this run: [['P', 3.0, 5.0], ['P', 3.0, 10.0], ['P', 3.0, 30.0], ['P', 3.0, 100.0]]. reproduction gate (re-running cell (3,1,2) with its original seeds): **PASS** (best 25.98971561000055 vs 25.98971561000055).

| r1 | r3 | D3 best | hits random | warm hits | A (H-inv) | A (H-term) | certified | syn3 raw | SGC excess % |
|:-:|:-:|--:|:-:|:-:|--:|--:|:-:|--:|--:|
| 3 | 5 | 33.159572866 | 0/128 | 8/11, 0/11, 0/11 | -0.00000 | +0.53185 | True | 1.000 | 39.7 |
| 3 | 10 | 33.159572866 | 0/128 | 10/11, 0/11, 9/11 | -0.00000 | +0.68905 | True | 1.000 | 72.4 |
| 3 | 30 | 33.159572866 | 0/128 | 10/11, 0/11, 10/11 | -0.00000 | +0.86604 | True | 1.000 | 203.4 |
| 3 | 100 | 33.159572866 | 0/128 | 10/11, 0/11, 10/11 | -0.00000 | +0.95523 | True | 1.000 | 661.8 |

warm hits are in the order lift H-inv, lift H-term, previous n=3 cell.

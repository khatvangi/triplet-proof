# n=3 vs n=4 optima at k = 20: OLD vs rebuilt descriptor matrix

NEW = `descriptor_rebuild/NEW_VECTORS.csv` (8 PCs, each z-scored). OLD = `n4/RESULTS.json`. Kernels, graph, usage, search and lift/certification logic are imported unchanged from `n4/n4.py`. G = 1 - D4*/D3*. D4* is pooled over searches of the same (usage, w4), as in n4.

## Gates (old vectors)

- SGC D = 23.922359743253192 vs n4 gate 1 23.922359743253192, |diff| = 0.0e+00 -> **PASS**
- n=3 U UNIF, n4's seed block, 256 restarts: best 23.067734683932564 vs n4 23.067734683932564 (|diff| 0.0e+00; n4 spread best..p50 = 23.067735..23.348328), hits 1/256 -> **PASS**
- info: SGC D with NEW vectors = 23.839669

Search: n=3 two independent batches of 256 random restarts (SA 500*N^2); n=4 128 random + 40 warm (SA 150*N^2). Restart reductions: ['n4 U UNIF Q-mon: 96+16', 'n4 U UNIF Q-inv: 96+16', 'n4 U K-str Q-mon: 96+16', 'n4 U K-str Q-inv: 96+16', 'n4 U K-int Q-mon: 96+16', 'n4 U K-int Q-inv: 96+16', 'n4 U w3=2 Q-mon: 96+16', 'n4 U w3=2 Q-inv: 96+16', 'n4 U w3=5 Q-mon: 96+16', 'n4 U w3=5 Q-inv: 96+16', 'n4 U w3=10 Q-mon: 96+16', 'n4 U w3=10 Q-inv: 96+16'].

## Usage U

| weights | hyp | OLD D3* | OLD D4* | OLD G | NEW D3* | NEW D4* | NEW G | NEW cert | NEW rand hits |
|---|---|---|---|---|---|---|---|---|---|
| UNIF (1, 1, 1) | Q-mon | 23.067735 | 22.766739 | +0.0130 | 23.365725 | 23.365725 | +0.0000 | yes | 0/96 |
| UNIF (1, 1, 1) | Q-inv | 23.067735 | 22.766739 | +0.0130 | 23.365725 | 23.365725 | +0.0000 | yes | 0/96 |
| K-str (0.98, 1, 1.39) | Q-mon | 23.456460 | 22.783332 | +0.0287 | 23.853647 | 23.143970 | +0.0298 | yes | 0/96 |
| K-str (0.98, 1, 1.39) | Q-inv | 23.456460 | 23.456460 | +0.0000 | 23.853647 | 23.853647 | +0.0000 | yes | 0/96 |
| K-int (0.99, 1, 1.79) | Q-mon | 24.007356 | 23.021302 | +0.0411 | 24.630416 | 23.259405 | +0.0557 | yes | 0/96 |
| K-int (0.99, 1, 1.79) | Q-inv | 24.007356 | 24.007356 | -0.0000 | 24.630416 | 24.630416 | +0.0000 | yes | 0/96 |
| w3=2 (1, 1, 2) | Q-mon | 24.466055 | 23.067735 | +0.0572 | 24.939213 | 23.365725 | +0.0631 | yes | 0/96 |
| w3=2 (1, 1, 2) | Q-inv | 24.466055 | 24.466055 | +0.0000 | 24.939213 | 24.939213 | +0.0000 | yes | 0/96 |
| w3=5 (1, 1, 5) | Q-mon | 28.101396 | 23.067735 | +0.1791 | 29.188913 | 23.365725 | +0.1995 | yes | 0/96 |
| w3=5 (1, 1, 5) | Q-inv | 28.101396 | 28.101396 | +0.0000 | 29.188913 | 29.188913 | -0.0000 | yes | 0/96 |
| w3=10 (1, 1, 10) | Q-mon | 32.905096 | 23.067735 | +0.2990 | 35.243498 | 23.365725 | +0.3370 | yes | 0/96 |
| w3=10 (1, 1, 10) | Q-inv | 32.905096 | 32.905096 | +0.0000 | 35.243498 | 35.243498 | +0.0000 | yes | 0/96 |

## Usage P

| weights | hyp | OLD D3* | OLD D4* | OLD G | NEW D3* | NEW D4* | NEW G | NEW cert | NEW rand hits |
|---|---|---|---|---|---|---|---|---|---|
| (not run) | | | | | | | | | |

## n=3 batch agreement (NEW vectors)

| usage | w sorted | batch A | batch B | agree (1e-6) | hits A / B | D3* |
|---|---|---|---|---|---|---|
| U | (1, 1, 1) | 23.365725312 | 23.397184312 | NO | 1/0 | 23.365725312 |
| U | (0.98, 1, 1.39) | 23.853646961 | 23.914448991 | NO | 1/0 | 23.853646961 |
| U | (1, 1, 1.39) | 24.058577914 | 24.116199041 | NO | 1/0 | 24.058577914 |
| U | (0.98, 1, 1) | 23.143969989 | 23.194399373 | NO | 1/0 | 23.143969989 |
| U | (1, 1.39, 1.39) | 28.253969445 | 28.227886530 | NO | 0/1 | 28.227886530 |
| U | (0.98, 1.39, 1.39) | 27.945979840 | 27.934039781 | NO | 0/1 | 27.934039781 |
| U | (0.99, 1, 1.79) | 24.652945526 | 24.630416249 | NO | 0/1 | 24.630416249 |
| U | (1, 1, 1.79) | 24.701193247 | 24.641889134 | NO | 0/1 | 24.641889134 |
| U | (0.99, 1, 1) | 23.259405067 | 23.289723908 | NO | 1/0 | 23.259405067 |
| U | (1, 1.79, 1.79) | 33.160783686 | 33.109263061 | NO | 0/1 | 33.109263061 |
| U | (0.99, 1.79, 1.79) | 33.027332898 | 32.833877688 | NO | 0/1 | 32.833877688 |
| U | (1, 1, 2) | 24.958768474 | 24.939213366 | NO | 0/1 | 24.939213366 |
| U | (1, 2, 2) | 35.540264563 | 35.632421604 | NO | 1/0 | 35.540264563 |
| U | (1, 1, 5) | 29.188912867 | 29.326995719 | NO | 1/0 | 29.188912867 |
| U | (1, 5, 5) | 71.539808773 | 71.502877545 | NO | 0/1 | 71.502877545 |
| U | (1, 1, 10) | 35.478730360 | 35.243498281 | NO | 0/1 | 35.243498281 |
| U | (1, 10, 10) | 130.918729415 | 130.948832542 | NO | 1/0 | 130.918729415 |

0/17 n=3 searches agree across batches. Hits are counted against the pooled best, so a batch that missed it shows 0.

## Conclusion (usage U only; P was not run)

- **Q-inv: unchanged.** G = 0 in all 6 NEW cells. In every one the best n=4 code is the silent position-4 lift of the best triplet code. OLD was also 0 in all non-UNIF cells.
- **Q-mon at in vitro rates: unchanged in kind, slightly larger.** K-str 2.9% -> 3.0%, K-int 4.1% -> 5.6%. That is still 'a few percent'. Graded weights: w3=2 5.7% -> 6.3%, w3=5 17.9% -> 20.0%, w3=10 29.9% -> 33.7%. In every Q-mon cell D4* equals min_p D3*(w\p). So G here is a ratio of two n=3 optima: the quadruplet gain is just deleting the heavy third position.
- **UNIF: OLD G = +1.3%, NEW G = 0.** OLD found an n=4 code 0.30 below every lift. NEW found none, with only 96+16 restarts. This cell is search-limited, not a matrix effect that can be claimed.
- Absolute D values are not comparable across matrices (the distance matrix changed; SGC D is 23.92 OLD vs 23.84 NEW). Only G is comparable.

## Caveats

- **The n=3 search is not converged.** 0/17 two-batch pairs agree to 1e-6. Every batch best is hit by exactly 1 of 256 restarts, and the batch-to-batch gap is 0.01-0.23 (up to 0.7% of D3). Every D3* and D4* is a search upper bound, and G carries an uncertainty of that order. This is enough to move the UNIF cell and the smaller Q-mon gains by about +/-1 percentage point. It is not enough to turn Q-inv's 0 into a gain, or the w3=5/10 gains into nothing.
- n=4 restarts were cut to the floor (96 random + 16 warm) in all 12 cells because of the budget. No random n=4 restart reached D4*; every best came from a warm start.
- Usage P was not run: the U block took 1074 s on a machine shared with other jobs (load average ~87 on 64 cores), and the run was stopped before P to stay near the 28-minute budget.

elapsed 1074.1 s; finished = False

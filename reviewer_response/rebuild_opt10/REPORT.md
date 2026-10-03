# Ten-class optima under the rebuilt descriptor matrix

Vectors: `reviewer_response/descriptor_rebuild/NEW_VECTORS.csv` (pc1..pc8, ddof=0 across 20 AAs), representatives X1=A, X2=L, X3=W, X4=S, X5=N, X6=D, X7=R, X8=H, X9=P, X10=C. Proteome usage P: Swiss-Prot 2026_03. Kernels, annealing + steepest descent, lift certification and the position grid are imported unchanged from `usage_cert/usage_cert.py` and `position_model/position_model.py` (the module `position_empirical/` imports); only the class distance matrix changes. OLD = `usage_cert/RESULTS.json`, `position_empirical/RESULTS.json`.

## Gate

Old vectors (aa_props_export.csv, pca1..8 z-scored), usage_cert seed list for R2 U n=2 (seed0 = 401000). Seed 401000 gave 20.031754609; seed 401001 gave 20.017516292, |diff| from 20.017516292 = 3.3e-11 → **PASS**. SGC projection (old vectors, U) reproduces usage_cert to 0.0e+00.

## (a) usage_cert R2/R3: equal position weights, free class sizes, all 4^n codons sense

| quantity | usage | OLD | NEW | conclusion unchanged? |
|---|:-:|--:|--:|---|
| D2* | U | 20.017516292 | 21.344573500 | value changes (new metric) |
| D3* | U | 20.017516292 | 21.344573500 | value changes (new metric) |
| triplet advantage D2*−D3* | U | 0.000000 (0.000%) | 0.000000 (0.000%) | yes (OLD sign 0, NEW sign 0) |
| hits n=2 random | U | 272/512 | 458/512 | — |
| hits n=3 random / warm | U | 68/512 / 61/64 | 14/512 / 58/64 | — |
| certified D3* ≤ lift | U | True | True | yes |
| n=3 syn p1/p2/p3 (best code) | U | — | 0.250/1.000/0.042 | — |
| D2* | P | 18.835846661 | 19.187429770 | value changes (new metric) |
| D3* | P | 18.766672019 | 18.885504832 | value changes (new metric) |
| triplet advantage D2*−D3* | P | 0.069175 (0.367%) | 0.301925 (1.574%) | yes (OLD sign +, NEW sign +) |
| hits n=2 random | P | 211/512 | 395/512 | — |
| hits n=3 random / warm | P | 68/512 / 34/64 | 19/256 / 9/64 | — |
| certified D3* ≤ lift | P | True | True | yes |
| n=3 syn p1/p2/p3 (best code) | P | — | 0.969/0.542/0.000 | — |

Manuscript statements under test: "D2* = D3* = 20.0175 (U)" and "18.836 vs 18.767 (P)".

## (c) reference codes, equal weights

| quantity | usage | OLD | NEW | conclusion unchanged? |
|---|:-:|--:|--:|---|
| SGC 10-class projection (61 sense) | U | 24.021366 | 25.382848 | — |
| doublet projection (9 classes) | U | 21.635431 | 22.070411 | — |
| SGC − doublet projection | U | +2.385935 | +3.312437 | yes |
| SGC excess over n=3 optimum | U | 4.003850 (20.00%) | 4.038274 (18.92%) | yes (SGC above optimum) |
| SGC 10-class projection (61 sense) | P | 21.019747 | 21.386877 | — |
| doublet projection (9 classes) | P | 19.685652 | 19.572164 | — |
| SGC − doublet projection | P | +1.334095 | +1.814713 | yes |
| SGC excess over n=3 optimum | P | 2.253075 (12.01%) | 2.501372 (13.24%) | yes (SGC above optimum) |

## (b) position_empirical grid

H-inv doublet = (w1, w2); H-term = (w1, w3); A = 1 − D3*/D2*. Conjecture: D3* = min_p D2*(w minus p) to 1e-9. "Unchanged" means the same silenced position, the same conjecture verdict, and the same sign of A_inv and A_term.

| set | usage | | D2*(H-inv) | D2*(H-term) | D3* | silenced | A_inv | A_term | min_p D2* − D3* | conj. | cert. | SGC Dpos | unchanged? |
|---|:-:|:-:|--:|--:|--:|:-:|--:|--:|--:|:-:|:-:|--:|---|
| L-Ec (1.37, 1.0, 0.71) | U | OLD | 23.073322 | 19.663333 | 16.607527 | 1 | +0.28023 | +0.15541 | 0.00e+00 | EQUAL | yes | 26.755029 |  |
| L-Ec (1.37, 1.0, 0.71) | U | NEW | 24.935313 | 21.362076 | 17.963205 | 1 | +0.27961 | +0.15911 | 0.00e+00 | EQUAL | yes | 28.228693 | yes |
| L-Ec-LOO (1.91, 1.0, 1.38) | U | OLD | 27.533146 | 32.001407 | 23.155911 | 1 | +0.15898 | +0.27641 | -7.11e-15 | EQUAL | yes | 34.218731 |  |
| L-Ec-LOO (1.91, 1.0, 1.38) | U | NEW | 29.911784 | 34.594878 | 25.030234 | 1 | +0.16320 | +0.27648 | 0.00e+00 | EQUAL | yes | 36.096566 | yes |
| L-Sc (1.4, 1.0, 0.71) | U | OLD | 23.321090 | 19.911101 | 16.607527 | 1 | +0.28788 | +0.16592 | 0.00e+00 | EQUAL | yes | 27.051397 |  |
| L-Sc (1.4, 1.0, 0.71) | U | NEW | 25.220074 | 21.631234 | 17.963205 | 1 | +0.28774 | +0.16957 | 0.00e+00 | EQUAL | yes | 28.539329 | yes |
| L-Sc-LOO (1.16, 1.0, 0.57) | U | OLD | 21.338946 | 16.282755 | 14.961326 | 1 | +0.29887 | +0.08116 | 0.00e+00 | EQUAL | yes | 24.235572 |  |
| L-Sc-LOO (1.16, 1.0, 0.57) | U | NEW | 22.907684 | 17.689415 | 16.239611 | 1 | +0.29108 | +0.08196 | 3.55e-15 | EQUAL | yes | 25.578573 | yes |
| K-str (0.98, 1.0, 1.39) | U | OLD | 19.782345 | 23.003329 | 19.782345 | 3 | +0.00000 | +0.14002 | 3.55e-15 | EQUAL | yes | 25.063102 |  |
| K-str (0.98, 1.0, 1.39) | U | NEW | 21.113277 | 24.886529 | 21.113277 | 3 | -0.00000 | +0.15162 | -7.11e-15 | EQUAL | yes | 26.500848 | yes |
| K-int (0.99, 1.0, 1.79) | U | OLD | 19.899930 | 26.424488 | 19.899930 | 3 | +0.00000 | +0.24691 | 0.00e+00 | EQUAL | yes | 26.432983 |  |
| K-int (0.99, 1.0, 1.79) | U | NEW | 21.228925 | 28.693784 | 21.228925 | 3 | -0.00000 | +0.26016 | -3.55e-15 | EQUAL | yes | 27.963460 | yes |
| UNIF (1.0, 1.0, 1.0) | U | OLD | 20.017516 | 20.017516 | 20.017516 | 3 | +0.00000 | +0.00000 | 0.00e+00 | EQUAL | yes | 24.021366 |  |
| UNIF (1.0, 1.0, 1.0) | U | NEW | 21.344573 | 21.344573 | 21.344573 | 1 | +0.00000 | +0.00000 | 0.00e+00 | EQUAL | yes | 25.382848 | NO |
| L-Ec (1.37, 1.0, 0.71) | P | OLD | 21.485736 | 18.100281 | 15.432510 | none | +0.28173 | +0.14739 | 1.79e-02 | FALSE | yes | 23.461384 |  |
| L-Ec (1.37, 1.0, 0.71) | P | NEW | 21.789980 | 18.265461 | 15.522658 | none | +0.28762 | +0.15016 | 1.40e-01 | FALSE | yes | 23.826412 | yes |
| L-Ec-LOO (1.91, 1.0, 1.38) | P | OLD | 25.353142 | 29.789256 | 21.536723 | none | +0.15053 | +0.27703 | 2.06e-02 | FALSE | yes | 29.398797 |  |
| L-Ec-LOO (1.91, 1.0, 1.38) | P | NEW | 25.588295 | 30.206630 | 21.670711 | none | +0.15310 | +0.28258 | 1.90e-01 | FALSE | yes | 29.931069 | yes |
| L-Sc (1.4, 1.0, 0.71) | P | OLD | 21.700592 | 18.315137 | 15.437569 | none | +0.28861 | +0.15711 | 1.28e-02 | FALSE | yes | 23.710413 |  |
| L-Sc (1.4, 1.0, 0.71) | P | NEW | 22.000997 | 18.476479 | 15.532755 | none | +0.29400 | +0.15932 | 1.30e-01 | FALSE | yes | 24.078934 | yes |
| L-Sc-LOO (1.16, 1.0, 0.57) | P | OLD | 19.981745 | 14.961932 | 13.757382 | none | +0.31150 | +0.08051 | 5.87e-02 | FALSE | yes | 21.414172 |  |
| L-Sc-LOO (1.16, 1.0, 0.57) | P | NEW | 20.312857 | 15.086847 | 13.768421 | none | +0.32218 | +0.08739 | 1.93e-01 | FALSE | yes | 21.732941 | yes |
| K-str (0.98, 1.0, 1.39) | P | OLD | 18.602367 | 21.395494 | 18.598192 | none | +0.00022 | +0.13074 | 4.17e-03 | FALSE | yes | 21.700601 |  |
| K-str (0.98, 1.0, 1.39) | P | NEW | 18.944360 | 21.687588 | 18.776254 | none | +0.00887 | +0.13424 | 1.68e-01 | FALSE | yes | 22.126158 | yes |
| K-int (0.99, 1.0, 1.79) | P | OLD | 18.719107 | 24.376979 | 18.719107 | 3 | -0.00000 | +0.23210 | -3.55e-15 | EQUAL | yes | 22.652198 |  |
| K-int (0.99, 1.0, 1.79) | P | NEW | 19.065895 | 24.622690 | 19.059559 | none | +0.00033 | +0.22594 | 6.34e-03 | FALSE | yes | 23.141234 | NO |
| UNIF (1.0, 1.0, 1.0) | P | OLD | 18.835847 | 18.835847 | 18.766672 | none | +0.00367 | +0.00367 | 6.92e-02 | FALSE | yes | 21.019747 |  |
| UNIF (1.0, 1.0, 1.0) | P | NEW | 19.187430 | 19.187430 | 18.885505 | none | +0.01574 | +0.01574 | 3.02e-01 | FALSE | yes | 21.386877 | yes |

Hits (NEW, n=3): U L-Ec rand 0/128 p1:5/6 p2:0/6 p3:0/6; U L-Ec-LOO rand 0/128 p1:5/6 p2:0/6 p3:0/6; U L-Sc rand 0/128 p1:5/6 p2:0/6 p3:0/6; U L-Sc-LOO rand 1/128 p1:5/6 p2:0/6 p3:0/6; U K-str rand 0/128 p1:0/6 p2:0/6 p3:5/6; U K-int rand 0/128 p1:0/6 p2:0/6 p3:5/6; U UNIF rand 3/128 p1:5/6 p2:5/6 p3:5/6; P L-Ec rand 4/128 p1:0/6 p2:0/6 p3:0/6; P L-Ec-LOO rand 2/128 p1:0/6 p2:0/6 p3:0/6; P L-Sc rand 5/128 p1:0/6 p2:0/6 p3:0/6; P L-Sc-LOO rand 6/128 p1:0/6 p2:0/6 p3:0/6; P K-str rand 4/128 p1:0/6 p2:0/6 p3:1/6; P K-int rand 1/128 p1:0/6 p2:0/6 p3:0/6; P UNIF rand 13/128 p1:0/6 p2:0/6 p3:1/6

## Restarts and reductions

(a): n=2 512 random; n=3 512 random + 64 silent-lift warm (usage_cert counts). (b): n=2 128 random + 16 warm from the previous n=2 optimum; n=3 128 random + 16 warm split over the three silent lifts (ceil → 6 per lift, 18 total). This is the stated floor, **below position_empirical's 256 + 32 (n=2: 256)**, chosen up front because the machine was shared (load ≈ 88 on 64 cores). Further reductions: ['(a) P n=3 random 512->256']. Elapsed 1111.0 s.

Manuscript Table 2 could not be located in `manuscript_JME.md` (no D2*/D3* table there); the OLD column for (b) is `position_empirical/RESULTS.json`, which is the source those values were drawn from.
## Reading the two NO flags and the search caveats

1. **UNIF U, silenced position 3 → 1: not a change of conclusion.** Under equal weights and uniform usage the three silent lifts are symmetric and all score D2* exactly, so which position ends up silenced depends only on which tied optimum the search returns. In this run the (a) U best code silences p2 (syn 0.250/1.000/0.042), while the (b) UNIF U best code silences p1. The substantive statement, D3* = D2* with no triplet advantage, is unchanged.
2. **K-int P: conclusion CHANGES.** OLD: D3* = min_p D2* (p3 silent, gap −3.6e-15). NEW: a non-silent triplet code reaches D3 = 19.059559396, which is 6.3e-3 below min_p D2* = 19.065894655. That code's best-known value is exact to 1e-14 (the lift check), so the conjecture is *falsified* for this cell under the new vectors, whatever the true D3* is. The best code was hit by 1/128 random starts and by 0/18 warm starts, so D3* may lie lower still. The n=2 side is converged (62/128 random hits at (0.99, 1.0)). Under NEW vectors the P row now fails the conjecture in **all seven** weight sets, where OLD had six of seven.
3. **U rows: EQUAL is a search result, not a proof.** In every U cell of (b), random restarts reached the n=3 optimum in only 0–3 of 128 runs, and the optimum came from silent-lift warm starts. A non-silent code below the lift that SA+descent cannot reach would go undetected. The same limitation held for OLD (the conjecture is only certified in the direction D3* ≤ min_p D2*).
4. **(a) P: the triplet advantage grows from 0.069 (0.37%) to 0.302 (1.57%),** so its sign and the certification are unchanged. The n=3 P search was cut to 256 random starts. The best code was hit by 19/256 random and 9/64 warm starts, so it is reproducible within this run. It has not been reproduced on disjoint seeds.
5. The HISTORY.md files were not updated, because the task restricts writes to this directory.

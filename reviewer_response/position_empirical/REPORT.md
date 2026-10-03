# Position-weighted optima at empirical misreading-rate vectors

Objective, kernels, SA + steepest descent, and gates are imported unchanged from `reviewer_response/position_model/position_model.py`. k = 10 (RAA10 reps), class sizes free, no stop codons, 64 triplet / 16 doublet codons. Weights from `position_rates/POSITION_RATES.md`.

H-inv doublet = (w1, w2); H-term doublet = (w1, w3). A = 1 - D3*/D2*. Conjecture: D3* = min_p D2*(w without p), tolerance 1e-9. Syn(p1/p2/p3) = raw fraction of arcs at each position that are synonymous in the best triplet code. SGC = 10-class projection of the standard code, 61 sense codons, same weights (not an optimum).

| set | w | usage | D2*(H-inv) | D2*(H-term) | D2* silence p1 / p2 / p3 | D3* | silenced | syn p1/p2/p3 | A_inv | A_term | min_p D2* - D3* | conj. | cert. | SGC Dpos |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L-Ec | (1.37, 1.0, 0.71) | U | 23.073322 | 19.663333 | 16.607527 / 19.663333 / 23.073322 | 16.607527 | 1 | 1.000/0.333/0.000 | +0.28023 | +0.15541 | 0.00e+00 | EQUAL | yes | 26.755029 |
| L-Ec-LOO | (1.91, 1.0, 1.38) | U | 27.533146 | 32.001407 | 23.155911 / 32.001407 / 27.533146 | 23.155911 | 1 | 1.000/0.000/0.333 | +0.15898 | +0.27641 | -7.11e-15 | EQUAL | yes | 34.218731 |
| L-Sc | (1.4, 1.0, 0.71) | U | 23.321090 | 19.911101 | 16.607527 / 19.911101 / 23.321090 | 16.607527 | 1 | 1.000/0.333/0.000 | +0.28788 | +0.16592 | 0.00e+00 | EQUAL | yes | 27.051397 |
| L-Sc-LOO | (1.16, 1.0, 0.57) | U | 21.338946 | 16.282755 | 14.961326 / 16.282755 / 21.338946 | 14.961326 | 1 | 1.000/0.333/0.000 | +0.29887 | +0.08116 | 0.00e+00 | EQUAL | yes | 24.235572 |
| K-str | (0.98, 1.0, 1.39) | U | 19.782345 | 23.003329 | 23.238500 / 23.003329 / 19.782345 | 19.782345 | 3 | 0.000/0.333/1.000 | +0.00000 | +0.14002 | 3.55e-15 | EQUAL | yes | 25.063102 |
| K-int | (0.99, 1.0, 1.79) | U | 19.899930 | 26.424488 | 26.542074 / 26.424488 / 19.899930 | 19.899930 | 3 | 0.000/0.333/1.000 | +0.00000 | +0.24691 | 0.00e+00 | EQUAL | yes | 26.432983 |
| UNIF | (1.0, 1.0, 1.0) | U | 20.017516 | 20.017516 | 20.017516 / 20.017516 / 20.017516 | 20.017516 | 3 | 0.333/0.000/1.000 | +0.00000 | +0.00000 | 0.00e+00 | EQUAL | yes | 24.021366 |
| L-Ec | (1.37, 1.0, 0.71) | P | 21.485736 | 18.100281 | 15.450391 / 18.100281 / 21.485736 | 15.432510 | none | 0.969/0.479/0.000 | +0.28173 | +0.14739 | 1.79e-02 | FALSE | yes | 23.461384 |
| L-Ec-LOO | (1.91, 1.0, 1.38) | P | 25.353142 | 29.789256 | 21.557355 / 29.789256 / 25.353142 | 21.536723 | none | 0.969/0.000/0.479 | +0.15053 | +0.27703 | 2.06e-02 | FALSE | yes | 29.398797 |
| L-Sc | (1.4, 1.0, 0.71) | P | 21.700592 | 18.315137 | 15.450391 / 18.315137 / 21.700592 | 15.437569 | none | 0.969/0.479/0.000 | +0.28861 | +0.15711 | 1.28e-02 | FALSE | yes | 23.710413 |
| L-Sc-LOO | (1.16, 1.0, 0.57) | P | 19.981745 | 14.961932 | 13.816034 / 14.961932 / 19.981745 | 13.757382 | none | 0.969/0.479/0.000 | +0.31150 | +0.08051 | 5.87e-02 | FALSE | yes | 21.414172 |
| K-str | (0.98, 1.0, 1.39) | P | 18.602367 | 21.395494 | 21.628973 / 21.395494 / 18.602367 | 18.598192 | none | 0.000/0.479/0.969 | +0.00022 | +0.13074 | 4.17e-03 | FALSE | yes | 21.700601 |
| K-int | (0.99, 1.0, 1.79) | P | 18.719107 | 24.376979 | 24.493719 / 24.376979 / 18.719107 | 18.719107 | 3 | 0.000/0.417/1.000 | -0.00000 | +0.23210 | -3.55e-15 | EQUAL | yes | 22.652198 |
| UNIF | (1.0, 1.0, 1.0) | P | 18.835847 | 18.835847 | 18.835847 / 18.835847 / 18.835847 | 18.766672 | none | 0.969/0.479/0.000 | +0.00367 | +0.00367 | 6.92e-02 | FALSE | yes | 21.019747 |

Gates: A True, B True, C (all three silent positions) True. Restarts: n=2 256 random (+16 warm from the previous n=2 optimum when one exists); n=3 256 random + 33 warm (11 per silent lift). Reductions: none. Elapsed 1041.5 s.
Notes (all from this run's RESULTS.json):

- Positive controls: UNIF reproduces position_model exactly (U 20.017516292; P 18.766672019, which is position_model's n=3 P w=(1,1,1) value).
- Certification: in all 14 cells D3* <= every silent-lift value, and the random restarts alone reach the certified optimum (best_from = random everywhere). Lift exactness |D2 - D3(lift)| <= 1.1e-14.
- Conjecture D3* = min_p D2*: holds in all 7 U cells, fails in 6 of 7 P cells (L-Ec, L-Ec-LOO, L-Sc, L-Sc-LOO, K-str, UNIF fail; K-int holds). Gaps are 0.004-0.069 (0.02-0.4% of D3*). In each failing cell the optimum is a near-silent code that breaks one synonymous arc pair at the silenced position (syn 0.969, and one class holds a single codon).
- Which position is silenced follows the largest weight: Landerer vectors (w1 largest) silence position 1; Katoh vectors (w3 largest) silence position 3. With the Landerer weights the optimal triplet is not wobble-like.
- A_inv = 0 whenever position 3 is silenced (the H-inv lift is optimal), up to the P tie-breaking gaps; A_term > 0 in every non-uniform cell because H-term puts the doublet's second weight on the high-error position.

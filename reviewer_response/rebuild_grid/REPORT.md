# rebuild_grid: ten-class position-weight grid, OLD vs NEW descriptors

NEW vectors: `reviewer_response/descriptor_rebuild/NEW_VECTORS.csv` (z-scored: pc1..pc8, ddof=0 across 20 AAs), reps {'X1': 'A', 'X2': 'L', 'X3': 'W', 'X4': 'S', 'X5': 'N', 'X6': 'D', 'X7': 'R', 'X8': 'H', 'X9': 'P', 'X10': 'C'}. Machinery imported read-only from `position_model/` (kernels, runSearch, 128 random + max(32, 2·sources) warm restarts, warm sources, seed schedule 7000+100000·cell, certification). Only the 10×10 distance matrix changed.

OLD = `position_model/RESULTS.json`, plus `finish/RESULTS.json` for the four P cells at r1=3, r3∈{5,10,30,100}.

**Gate** (OLD vectors, U n=2 w=(1,5) (H-term of r1=1, r3=5), original seeds and warm start): best 53.053250863750065 vs original 53.053250863750065, |Δ| = 0.0e+00; random hits 48/128 vs 48/128 → **PASS**.

Restart reductions: none (every cell ran the full 128 + 32).

## Findings

- Every qualitative feature of the grid survives the descriptor rebuild. At r1=1, H-inv gives A=0 in 11 of the 12 cells (all except P (1,1,1)): the best triplet code is the p3-silenced lift. H-term A rises with r3 (U: 0.305→0.976; P: 0.268→0.973). At r1=3, r3≥5, A_inv=0 with p3 silenced. The only cells with A_inv>0 are r1=3, r3∈{1,2}, where p1 is silenced (or, in P (3,2), no position is), and P (1,1,1).
- The only cells where the triplet beats every silent lift ('below') are P (1,1,1) (A=+0.0157 new vs +0.0037 old) and P (3,1,2) (26.082 vs min D2 26.221). The same two cells were 'below' with the old vectors.
- NEW search failures: in 3 cells (U (3,1,1), U (3,1,2), P (3,1,1)) the search ended ABOVE the drop-p1 silent lift. position_model's H-inv certificate passes all three anyway, because it never checks the p1 lift. With the old vectors these cells reached the lift value. Their true D3* equals min D2 or lies below it; the corrected A values are listed below. The r1=3, r3=1/2 cells need a rerun with a drop-p1 lift warm start before anyone quotes them.
- All numbers shift by only a few percent (D2*/D3* up 6.6–8.5% in U; in P they move between −0.6% and +1.9%). The SGC-excess column changes by at most about 60 percentage points, at r3=100 where the excess is already about 700–1600%.

Weights: triplet (r1, 1, r3). H-inv doublet (r1, 1), which drops p3. H-term doublet (r1, r3), which drops p2. A = 1 − D3*/D2*. min D2 = min over p of D2*(w minus p), where dropping p1 gives (1, r3). '=min' means D3* equals min D2 to 1e-9 relative. 'sil' lists the codon positions whose arcs are all synonymous in the best n=3 code. 'cert' is position_model's certificate: D3* ≤ value of the H-inv silent lift. In the min-D2 column, 'yes' means D3* equals min D2. 'below' means the triplet code beats every silent lift. 'ABOVE' means the search ended higher than a known silent lift: a **search failure**. That lift is a better n=3 code, so the true D3* is ≤ min D2. position_model's certificate misses this because it checks only the H-inv (drop-p3) lift.

## usage U

| r1 | r3 | D2*(H-inv) old → new | D2*(H-term) old → new | D3* old → new | A_inv old → new | A_term old → new | sil old → new | D3*=min D2 (argmin p) old → new | cert old / new | rand hits new |
|:-:|:-:|--:|--:|--:|--:|--:|:-:|:-:|:-:|:-:|
| 1 | 1 | 20.017516 → 21.344573 | 20.017516 → 21.344573 | 20.017516 → 21.344573 | +0.0000 → +0.0000 | +0.0000 → +0.0000 | 2 → 3 | yes (1/2/3) → yes (1/2/3) | True / True | 2/128 |
| 1 | 2 | 20.017516 → 21.344573 | 28.276450 → 30.719257 | 20.017516 → 21.344573 | -0.0000 → +0.0000 | +0.2921 → +0.3052 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 1 | 5 | 20.017516 → 21.344573 | 53.053251 → 57.317651 | 20.017516 → 21.344573 | +0.0000 → +0.0000 | +0.6227 → +0.6276 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 1 | 10 | 20.017516 → 21.344573 | 93.895186 → 100.987520 | 20.017516 → 21.344573 | +0.0000 → +0.0000 | +0.7868 → +0.7886 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 1 | 30 | 20.017516 → 21.344573 | 256.442267 → 275.666997 | 20.017516 → 21.344573 | +0.0000 → +0.0000 | +0.9219 → +0.9226 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 1 | 100 | 20.017516 → 21.344573 | 825.357048 → 887.045165 | 20.017516 → 21.344573 | +0.0000 → +0.0000 | +0.9757 → +0.9759 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 3 | 1 | 36.535384 → 39.638479 | 36.535384 → 39.638479 | 20.017516 → 21.550588 | +0.4521 → +0.4563 | +0.4521 → +0.4563 | 1 → 1 | yes (1) → ABOVE (1) | True / True | 0/128 |
| 3 | 2 | 36.535384 → 39.638479 | 48.293966 → 52.309082 | 28.276450 → 30.769816 | +0.2261 → +0.2237 | +0.4145 → +0.4118 | 1 → 1 | yes (1) → ABOVE (1) | True / True | 1/128 |
| 3 | 5 | 36.535384 → 39.638479 | 76.570416 → 83.055080 | 36.535384 → 39.638479 | -0.0000 → +0.0000 | +0.5229 → +0.5227 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 3 | 10 | 36.535384 → 39.638479 | 117.865084 → 127.760524 | 36.535384 → 39.638479 | -0.0000 → +0.0000 | +0.6900 → +0.6897 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 3 | 30 | 36.535384 → 39.638479 | 281.685559 → 302.962561 | 36.535384 → 39.638479 | -0.0000 → +0.0000 | +0.8703 → +0.8692 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 3 | 100 | 36.535384 → 39.638479 | 850.600340 → 914.340729 | 36.535384 → 39.638479 | -0.0000 → +0.0000 | +0.9570 → +0.9566 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |

## usage P

| r1 | r3 | D2*(H-inv) old → new | D2*(H-term) old → new | D3* old → new | A_inv old → new | A_term old → new | sil old → new | D3*=min D2 (argmin p) old → new | cert old / new | rand hits new |
|:-:|:-:|--:|--:|--:|--:|--:|:-:|:-:|:-:|:-:|
| 1 | 1 | 18.835847 → 19.187430 | 18.835847 → 19.187430 | 18.766672 → 18.885505 | +0.0037 → +0.0157 | +0.0037 → +0.0157 | none → none | below (1/2/3) → below (1/2/3) | True / True | 10/128 |
| 1 | 2 | 18.835847 → 19.187430 | 25.997710 → 26.221348 | 18.835847 → 19.187430 | +0.0000 → +0.0000 | +0.2755 → +0.2683 | 3 → 3 | yes (3) → yes (3) | True / True | 1/128 |
| 1 | 5 | 18.835847 → 19.187430 | 47.295310 → 47.323103 | 18.835847 → 19.187430 | +0.0000 → +0.0000 | +0.6017 → +0.5945 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 1 | 10 | 18.835847 → 19.187430 | 82.513731 → 82.492694 | 18.835847 → 19.187430 | +0.0000 → +0.0000 | +0.7717 → +0.7674 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 1 | 30 | 18.835847 → 19.187430 | 223.387412 → 222.493493 | 18.835847 → 19.187430 | +0.0000 → +0.0000 | +0.9157 → +0.9138 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 1 | 100 | 18.835847 → 19.187430 | 716.445299 → 711.944786 | 18.835847 → 19.187430 | +0.0000 → +0.0000 | +0.9737 → +0.9730 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 3 | 1 | 33.159573 → 33.255266 | 33.159573 → 33.255266 | 18.835847 → 19.313355 | +0.4320 → +0.4192 | +0.4320 → +0.4192 | 1 → 1 | yes (1) → ABOVE (1) | True / True | 0/128 |
| 3 | 2 | 33.159573 → 33.255266 | 44.833556 → 45.408778 | 25.989716 → 26.082142 | +0.2162 → +0.2157 | +0.4203 → +0.4256 | none → none | below (1) → below (1) | True / True | 4/128 |
| 3 | 5 | 33.159573 → 33.255266 | 70.831266 → 71.630126 | 33.159573 → 33.255266 | -0.0000 → +0.0000 | +0.5319 → +0.5357 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 3 | 10 | 33.159573 → 33.255266 | 106.640582 → 106.799717 | 33.159573 → 33.255266 | -0.0000 → +0.0000 | +0.6891 → +0.6886 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 3 | 30 | 33.159573 → 33.255266 | 247.541192 → 247.478082 | 33.159573 → 33.255266 | -0.0000 → +0.0000 | +0.8660 → +0.8656 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |
| 3 | 100 | 33.159573 → 33.255266 | 740.599078 → 737.402093 | 33.159573 → 33.255266 | -0.0000 → +0.0000 | +0.9552 → +0.9549 | 3 → 3 | yes (3) → yes (3) | True / True | 0/128 |

## NEW cells where the search ended above min D2 (search failures)

- U (r1=3, r3=1): search D3* = 21.550588406 > min D2 = 21.344573500 (drop p1). Silent-lift values under w for dropping p1/p2/p3: 21.344573500, 39.638479013, 39.638479013. Using the best known n=3 value (the lift): A_inv = +0.46152, A_term = +0.46152. Treat these as bounds; this cell is not certified.
- U (r1=3, r3=2): search D3* = 30.769815835 > min D2 = 30.719256792 (drop p1). Silent-lift values under w for dropping p1/p2/p3: 30.719256792, 52.309082165, 39.638479013. Using the best known n=3 value (the lift): A_inv = +0.22501, A_term = +0.41274. Treat these as bounds; this cell is not certified.
- P (r1=3, r3=1): search D3* = 19.313354749 > min D2 = 19.187429770 (drop p1). Silent-lift values under w for dropping p1/p2/p3: 19.187429770, 33.255266237, 33.255266237. Using the best known n=3 value (the lift): A_inv = +0.42303, A_term = +0.42303. Treat these as bounds; this cell is not certified.

## SGC 10-class projection (61 sense codons) under each triplet weight

| usage | r1 | r3 | SGC Dpos old | SGC Dpos new | excess over D3* old % | excess over D3* new % | new, over best known (min of D3*, min D2) % |
|:-:|:-:|:-:|--:|--:|--:|--:|--:|
| U | 1 | 1 | 24.021366 | 25.382848 | 20.0 | 18.9 | 18.9 |
| U | 1 | 2 | 27.199096 | 28.780515 | 35.9 | 34.8 | 34.8 |
| U | 1 | 5 | 36.732284 | 38.973515 | 83.5 | 82.6 | 82.6 |
| U | 1 | 10 | 52.620931 | 55.961849 | 162.9 | 162.2 | 162.2 |
| U | 1 | 30 | 116.175520 | 123.915185 | 480.4 | 480.5 | 480.5 |
| U | 1 | 100 | 338.616580 | 361.751860 | 1591.6 | 1594.8 | 1594.8 |
| U | 3 | 1 | 43.779230 | 46.091868 | 118.7 | 113.9 | 115.9 |
| U | 3 | 2 | 46.956959 | 49.489535 | 66.1 | 60.8 | 61.1 |
| U | 3 | 5 | 56.490148 | 59.682536 | 54.6 | 50.6 | 50.6 |
| U | 3 | 10 | 72.378795 | 76.670870 | 98.1 | 93.4 | 93.4 |
| U | 3 | 30 | 135.933383 | 144.624205 | 272.1 | 264.9 | 264.9 |
| U | 3 | 100 | 358.374443 | 382.460880 | 880.9 | 864.9 | 864.9 |
| P | 1 | 1 | 21.019747 | 21.386877 | 12.0 | 13.2 | 13.2 |
| P | 1 | 2 | 23.191216 | 23.714131 | 23.1 | 23.6 | 23.6 |
| P | 1 | 5 | 29.705623 | 30.695895 | 57.7 | 60.0 | 60.0 |
| P | 1 | 10 | 40.562968 | 42.332168 | 115.3 | 120.6 | 120.6 |
| P | 1 | 30 | 83.992346 | 88.877260 | 345.9 | 363.2 | 363.2 |
| P | 1 | 100 | 235.995172 | 251.785082 | 1152.9 | 1212.2 | 1212.2 |
| P | 3 | 1 | 37.621707 | 38.221681 | 99.7 | 97.9 | 99.2 |
| P | 3 | 2 | 39.793176 | 40.548936 | 53.1 | 55.5 | 55.5 |
| P | 3 | 5 | 46.307583 | 47.530699 | 39.7 | 42.9 | 42.9 |
| P | 3 | 10 | 57.164927 | 59.166972 | 72.4 | 77.9 | 77.9 |
| P | 3 | 30 | 100.594306 | 105.712064 | 203.4 | 217.9 | 217.9 |
| P | 3 | 100 | 252.597132 | 268.619886 | 661.8 | 707.8 | 707.8 |

The SGC Dpos with OLD vectors recomputed here matches position_model's stored values to 1e-9 (asserted).


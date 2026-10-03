# rebuild_grid2: completed rebuilt-matrix grids (fixed certificate) + rebuild_opt20 usage P

Read-only reuse: `rebuild_grid/rebuild_grid.py` (imported by `grid_fixed.py`) and `rebuild_opt20/run.py` (copied to `opt20P.py`, edits marked `rebuild_grid2`). NEW vectors = `descriptor_rebuild/NEW_VECTORS.csv`.

## (A) ten-class position-weight grid, NEW vectors

**Provenance.** The `rebuild_grid` run kept going after its REPORT.md was written. It finished all 48 cells (24 n=2, 24 n=3) with no restart cuts (`rebuild_grid/RESULTS.json`, finished = True). The n=2 optima used here come from that run. Each doublet label was re-scored and matched its stored D2* to 1e-9. The n=3 cells listed below as `rerun` were searched again in `grid_fixed.py`, with two fixes:

1. **Certificate:** D3* is checked against all three silent lifts (drop p1, drop p2, drop p3). Each lift is built at its own dropped position, so its value under w equals its doublet optimum exactly (asserted).
2. **Warm starts:** the three silent lifts, plus the previous n=3 cell (rebuild_grid's neighbour rule). rebuild_grid had only two lift warm starts. Its `lift_Hinv` is the drop-p3 lift. Its `lift_Hterm` came from `UC.liftLabels`, which is silent at p3, not p2, so it was never the H-term lift. The drop-p1 lift was missing, and that is the lift that wins whenever r1 = 3 and r3 ≤ 2.

Random restarts reuse rebuild_grid's seed schedule (7000 + 100000·cell_index). With 128 random restarts, the random part reproduces rebuild_grid's random restarts exactly. U r1=1 rows were not rerun; they are rebuild_grid's values, re-checked here against all three lift values it stored.

Columns: weights (r1, 1, r3). D2*(H-inv) = D2*(r1, 1), which drops p3. D2*(H-term) = D2*(r1, r3), which drops p2. A = 1 − D3*/D2*. 'vs min D2' compares D3* with min_p D2*(w∖p). 'yes' = equal (1e-9 rel), 'below' = the triplet code beats every silent lift, 'ABOVE' = search failure. argmin p is in parentheses. 'sil' = codon position(s) whose arcs are all synonymous in the best n=3 code. 'cert3' = D3*(search) ≤ min of the three lift values.

### usage U

| r1 | r3 | D2*(H-inv) | D2*(H-term) | D3* | A_inv | A_term | sil | D3* vs min D2 (argmin p) | cert3 | rand hits | restarts | rebuild_grid D3* (vs min D2) | source |
|:-:|:-:|--:|--:|--:|--:|--:|:-:|:-:|:-:|:-:|:-:|--:|:-:|
| 1 | 1 | 21.344573 | 21.344573 | 21.344573500 | +0.0000 | +0.0000 | 3 | yes (1/2/3) | yes | 2/128 | 128+32 | 21.344573500 (yes) | rebuild_grid |
| 1 | 2 | 21.344573 | 30.719257 | 21.344573500 | +0.0000 | +0.3052 | 3 | yes (3) | yes | 0/128 | 128+32 | 21.344573500 (yes) | rebuild_grid |
| 1 | 5 | 21.344573 | 57.317651 | 21.344573500 | +0.0000 | +0.6276 | 3 | yes (3) | yes | 0/128 | 128+32 | 21.344573500 (yes) | rebuild_grid |
| 1 | 10 | 21.344573 | 100.987520 | 21.344573500 | +0.0000 | +0.7886 | 3 | yes (3) | yes | 0/128 | 128+32 | 21.344573500 (yes) | rebuild_grid |
| 1 | 30 | 21.344573 | 275.666997 | 21.344573500 | +0.0000 | +0.9226 | 3 | yes (3) | yes | 0/128 | 128+32 | 21.344573500 (yes) | rebuild_grid |
| 1 | 100 | 21.344573 | 887.045165 | 21.344573500 | +0.0000 | +0.9759 | 3 | yes (3) | yes | 0/128 | 128+32 | 21.344573500 (yes) | rebuild_grid |
| 3 | 1 | 39.638479 | 39.638479 | 21.344573500 | +0.4615 | +0.4615 | 1 | yes (1) | yes | 0/128 | 128+32 | 21.550588406 (ABOVE) | rerun |
| 3 | 2 | 39.638479 | 52.309082 | 30.719256792 | +0.2250 | +0.4127 | 1 | yes (1) | yes | 0/128 | 128+32 | 30.769815835 (ABOVE) | rerun |
| 3 | 5 | 39.638479 | 83.055080 | 39.638479013 | +0.0000 | +0.5227 | 3 | yes (3) | yes | 0/128 | 128+32 | 39.638479013 (yes) | rerun |
| 3 | 10 | 39.638479 | 127.760524 | 39.638479013 | +0.0000 | +0.6897 | 3 | yes (3) | yes | 0/128 | 128+32 | 39.638479013 (yes) | rerun |
| 3 | 30 | 39.638479 | 302.962561 | 39.638479013 | +0.0000 | +0.8692 | 3 | yes (3) | yes | 0/128 | 128+32 | 39.638479013 (yes) | rerun |
| 3 | 100 | 39.638479 | 914.340729 | 39.638479013 | +0.0000 | +0.9566 | 3 | yes (3) | yes | 0/128 | 128+32 | 39.638479013 (yes) | rerun |

### usage P

| r1 | r3 | D2*(H-inv) | D2*(H-term) | D3* | A_inv | A_term | sil | D3* vs min D2 (argmin p) | cert3 | rand hits | restarts | rebuild_grid D3* (vs min D2) | source |
|:-:|:-:|--:|--:|--:|--:|--:|:-:|:-:|:-:|:-:|:-:|--:|:-:|
| 1 | 1 | 19.187430 | 19.187430 | 18.885504832 | +0.0157 | +0.0157 | none | below (1/2/3) | yes | 8/96 | 96+16 | 18.885504832 (below) | rerun |
| 1 | 2 | 19.187430 | 26.221348 | 19.187429770 | +0.0000 | +0.2683 | 3 | yes (3) | yes | 1/96 | 96+16 | 19.187429770 (yes) | rerun |
| 1 | 5 | 19.187430 | 47.323103 | 19.187429770 | +0.0000 | +0.5945 | 3 | yes (3) | yes | 0/96 | 96+16 | 19.187429770 (yes) | rerun |
| 1 | 10 | 19.187430 | 82.492694 | 19.187429770 | +0.0000 | +0.7674 | 3 | yes (3) | yes | 0/96 | 96+16 | 19.187429770 (yes) | rerun |
| 1 | 30 | 19.187430 | 222.493493 | 19.187429770 | +0.0000 | +0.9138 | 3 | yes (3) | yes | 0/96 | 96+16 | 19.187429770 (yes) | rerun |
| 1 | 100 | 19.187430 | 711.944786 | 19.187429770 | +0.0000 | +0.9730 | 3 | yes (3) | yes | 0/96 | 96+16 | 19.187429770 (yes) | rerun |
| 3 | 1 | 33.255266 | 33.255266 | 19.187429770 | +0.4230 | +0.4230 | 1 | yes (1) | yes | 0/96 | 96+16 | 19.313354749 (ABOVE) | rerun |
| 3 | 2 | 33.255266 | 45.408778 | 26.082142011 | +0.2157 | +0.4256 | none | below (1) | yes | 4/96 | 96+16 | 26.082142011 (below) | rerun |
| 3 | 5 | 33.255266 | 71.630126 | 33.255266237 | +0.0000 | +0.5357 | 3 | yes (3) | yes | 0/96 | 96+16 | 33.255266237 (yes) | rerun |
| 3 | 10 | 33.255266 | 106.799717 | 33.255266237 | +0.0000 | +0.6886 | 3 | yes (3) | yes | 0/96 | 96+16 | 33.255266237 (yes) | rerun |
| 3 | 30 | 33.255266 | 247.478082 | 33.255266237 | +0.0000 | +0.8656 | 3 | yes (3) | yes | 0/96 | 96+16 | 33.255266237 (yes) | rerun |
| 3 | 100 | 33.255266 | 737.402093 | 33.255266237 | +0.0000 | +0.9549 | 3 | yes (3) | yes | 0/96 | 96+16 | 33.255266237 (yes) | rerun |

### OLD vectors (position_model / finish), for comparison

See `rebuild_grid/REPORT.md` (the 'old' side of each 'old → new' column). The OLD D3* values are unchanged by this rerun.

## (B) rebuild_opt20, usage P: n=3 vs n=4 optima at k = 20, OLD vs NEW matrix

Gate A (old vectors, SGC D): 23.922359743253192 vs n4 23.922359743253192, |diff| 0.0e+00 -> **PASS**. Gate B was not rerun; rebuild_opt20 passed it with identical code. Search: n=3 two independent batches of 256 random restarts (SA 500*N^2). n=4 random + warm from the triplet lift and all four deletion lifts (SA 150*N^2). Restart reductions: ['n4 P UNIF Q-mon: 96+40', 'n4 P UNIF Q-inv: 96+40', 'n4 P K-str Q-mon: 96+40', 'n4 P K-str Q-inv: 96+40', 'n4 P K-int Q-mon: 96+40', 'n4 P K-int Q-inv: 96+40', 'n4 P w3=2 Q-mon: 96+40', 'n4 P w3=2 Q-inv: 96+40', 'n4 P w3=5 Q-mon: 96+40', 'n4 P w3=5 Q-inv: 96+40', 'n4 P w3=10 Q-mon: 96+40', 'n4 P w3=10 Q-inv: 96+40']. The seed counter resumes where rebuild_opt20's U block stopped. Finished: True; elapsed 750.7 s.

| weights | hyp | OLD D3* | OLD D4* | OLD G | NEW D3* | NEW D4* | NEW G | NEW cert | NEW rand hits | D4*=min del |
|---|---|---|---|---|---|---|---|---|---|---|
| UNIF (1, 1, 1) | Q-mon | 22.109961 | 22.006804 | +0.0047 | 21.989661 | 21.987900 | +0.0001 | yes | 0/96 | no |
| UNIF (1, 1, 1) | Q-inv | 22.109961 | 22.006804 | +0.0047 | 21.989661 | 21.987900 | +0.0001 | yes | 0/96 | no |
| K-str (0.98, 1, 1.39) | Q-mon | 22.631766 | 21.909669 | +0.0319 | 22.300221 | 21.817840 | +0.0216 | yes | 0/96 | yes |
| K-str (0.98, 1, 1.39) | Q-inv | 22.631766 | 22.613483 | +0.0008 | 22.300221 | 22.259774 | +0.0018 | yes | 0/96 | no |
| K-int (0.99, 1, 1.79) | Q-mon | 23.223370 | 22.003514 | +0.0525 | 22.731009 | 21.913324 | +0.0360 | yes | 0/96 | yes |
| K-int (0.99, 1, 1.79) | Q-inv | 23.223370 | 23.176235 | +0.0020 | 22.731009 | 22.699812 | +0.0014 | yes | 0/96 | no |
| w3=2 (1, 1, 2) | Q-mon | 23.573031 | 22.109961 | +0.0621 | 23.067391 | 21.989661 | +0.0467 | yes | 0/96 | yes |
| w3=2 (1, 1, 2) | Q-inv | 23.573031 | 23.541056 | +0.0014 | 23.067391 | 23.040627 | +0.0012 | yes | 0/96 | no |
| w3=5 (1, 1, 5) | Q-mon | 26.753633 | 22.109961 | +0.1736 | 26.143160 | 21.989661 | +0.1589 | yes | 0/96 | yes |
| w3=5 (1, 1, 5) | Q-inv | 26.753633 | 26.644814 | +0.0041 | 26.143160 | 26.143160 | -0.0000 | yes | 0/96 | yes |
| w3=10 (1, 1, 10) | Q-mon | 30.956656 | 22.109961 | +0.2858 | 31.330458 | 21.989661 | +0.2981 | yes | 0/96 | yes |
| w3=10 (1, 1, 10) | Q-inv | 30.956656 | 30.956656 | +0.0000 | 31.330458 | 31.330458 | +0.0000 | yes | 0/96 | yes |

### n=3 batch agreement (NEW vectors, usage P)

| usage | w sorted | batch A | batch B | agree (1e-6) | hits A / B | D3* |
|---|---|---|---|---|---|---|
| P | (1, 1, 1) | 22.008808262 | 21.989661128 | NO | 0/1 | 21.989661128 |
| P | (0.98, 1, 1.39) | 22.300220777 | 22.300500205 | NO | 1/0 | 22.300220777 |
| P | (1, 1, 1.39) | 22.464458929 | 22.399866473 | NO | 0/1 | 22.399866473 |
| P | (0.98, 1, 1) | 21.817840465 | 21.839450570 | NO | 1/0 | 21.817840465 |
| P | (1, 1.39, 1.39) | 26.432934146 | 26.370418965 | NO | 0/1 | 26.370418965 |
| P | (0.98, 1.39, 1.39) | 26.085663856 | 26.212858153 | NO | 1/0 | 26.085663856 |
| P | (0.99, 1, 1.79) | 22.916959372 | 22.731009432 | NO | 0/1 | 22.731009432 |
| P | (1, 1, 1.79) | 23.032035981 | 22.962334596 | NO | 0/1 | 22.962334596 |
| P | (0.99, 1, 1) | 21.913324364 | 21.930276662 | NO | 1/0 | 21.913324364 |
| P | (1, 1.79, 1.79) | 30.710901630 | 30.802306776 | NO | 1/0 | 30.710901630 |
| P | (0.99, 1.79, 1.79) | 30.621615003 | 30.703960171 | NO | 1/0 | 30.621615003 |
| P | (1, 1, 2) | 23.231914552 | 23.067391278 | NO | 0/1 | 23.067391278 |
| P | (1, 2, 2) | 32.989655030 | 32.989655030 | yes | 1/1 | 32.989655030 |
| P | (1, 1, 5) | 26.415815419 | 26.143159901 | NO | 0/1 | 26.143159901 |
| P | (1, 5, 5) | 65.290997603 | 65.434415149 | NO | 1/0 | 65.290997603 |
| P | (1, 1, 10) | 31.348471828 | 31.330458487 | NO | 0/1 | 31.330458487 |
| P | (1, 10, 10) | 118.973011007 | 118.914027490 | NO | 0/1 | 118.914027490 |

1/17 n=3 searches agree across batches.

## Conclusions

- **(A) The fixed certificate removes every ABOVE cell.**
  - U (3,1): certified at 21.344573500 (A_inv = A_term = +0.4615).
  - U (3,2): certified at 30.719256792 (A_inv = +0.2250, A_term = +0.4127).
  - P (3,1): also ABOVE in rebuild_grid (19.313354749). rebuild_grid finished it after its REPORT.md was written, so the report never flagged it. Now certified at 19.187429770 (A = +0.4230).
  - All three equal D2* at dropped position p = 1, with position 1 silenced. The winning restarts came from the drop-p1 lift, which rebuild_grid never used as a warm start.
  - The random part reproduces rebuild_grid's random best exactly in all six U reruns (128 random restarts). The value changed only through the warm starts.
- **(A) NEW grid, all 24 cells certified against all three lifts.**
  - 'below' survives only at P (1,1) (18.885504832 vs 19.187429770) and P (3,2) (26.082142011 vs 26.221348003). These are the same two 'below' cells as OLD.
  - Every other cell is 'yes'. Position 3 is silenced when r1 = 1 or r3 ≥ 5. Position 1 is silenced when r1 = 3 and r3 ≤ 2. This matches OLD's pattern.
- **(B) Usage P reproduces the U pattern.**
  - Q-mon: G = 2.2% (K-str), 3.6% (K-int), and 4.7% / 15.9% / 29.8% (w3 = 2 / 5 / 10). In every non-UNIF Q-mon cell, D4* = min_p D3*(w∖p): the quadruplet gain is just deleting the heaviest position. OLD P gave 3.2 / 5.3 / 6.2 / 17.4 / 28.6%, so NEW is the same kind of effect, slightly smaller at the in vitro rates.
  - Q-inv: G ≤ +0.0018 everywhere. In the UNIF, K-str, K-int and w3=2 cells, the n=4 search went below the lift of the triplet best by 0.002 to 0.04. That is inside the n=3 search noise: the two batches differ by 0.0003 to 0.27 (median about 0.07), and only 1/17 agree. A better, still unfound n=3 code, whose lift would close the gap, is the plausible explanation. This is not evidence of a quadruplet gain. OLD P Q-inv gave the same size (+0.0008 to +0.0047).
  - UNIF P: G = +0.0001 NEW vs +0.0047 OLD. Both are below the noise.

## Caveats

- (B) n=4 used 96 random + 40 warm in all 12 cells. I set this for time; the budget planner did not choose it, and that is why every cell is listed under 'Restart reductions'. No random n=4 restart reached D4*; every best came from a warm start. Gate B was copied from rebuild_opt20, not rerun. Gate A passed.
- (A) P reruns used 96 random + 16 warm (the floor), for time. U r1 = 3 reruns used 128 + 32. U r1 = 1 cells are rebuild_grid's values (not rerun); all are certified against its three stored lift values.
- The n=2 optima come from rebuild_grid (128 + 32 restarts, not rerun). If a doublet search failed, the error would carry into D2*, min D2 and the lifts.

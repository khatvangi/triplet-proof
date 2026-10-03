# Doublet vs triplet per-residue optima of D (arch_cert)

Script: `arch_cert.py` (run from repo root: `python reviewer_response/arch_cert/arch_cert.py`).
Raw numbers, per-scenario settings and the best labelings (codon -> class) are in `RESULTS.json`;
console log in `run.log`.

D = (1/N_sense) * sum over sense codons of sum over sense Hamming-1 neighbours of dist(label, label')
= 2*H1/N_sense. Class vectors: pca1..pca8 of `reviewer_response/aa_props_export.csv`, z-scored
across the 20 AAs (ddof=0), representatives X1=A X2=L X3=W X4=S X5=N X6=D X7=R X8=H X9=P X10=C.

## Gates (both passed)

| gate | computed D | expected | pass |
|---|---|---|---|
| SGC RAA10 triplet projection, UAA/UAG/UGA excluded (N=61) | 21.268639884706705 | 21.268639884706705 | yes |
| majority-vote doublet projection (`build_doublet_table`, imported from `run_publication_controls.py`) (N=16) | 18.937298834660997 | 18.937298834660997 | yes |

The doublet projection has no stops and 9 classes: X8 (His) is lost to the majority vote.

## Results

Ratio = D(n=3) / D(n=2) for the same scenario and row pair. hits = restarts whose final D is within
1e-9 of the best value found in that scenario.

| scenario | k | n | best D | hits/restarts | p5 | p50 | D(n=3)/D(n=2) | 3·log4(k) |
|---|---|---|---|---|---|---|---|---|
| S1 blind | 4 | 2 | 3.000000000 | 256/256 | 3.000000 | 3.000000 | 1.000000 | 3 (equal) |
| S1 blind | 4 | 3 | 3.000000000 | 256/256 | 3.000000 | 3.000000 | | 3 (equal) |
| S1 blind | 8 | 2 | 5.000000000 | 256/256 | 5.000000 | 5.000000 | 1.000000 | — |
| S1 blind | 8 | 3 | 5.000000000 | 254/256 | 5.000000 | 5.000000 | | — |
| S1 blind | 10 | 2 | 5.250000000 | 256/256 | 5.250000 | 5.250000 | 1.107143 | — |
| S1 blind | 10 | 3 | 5.812500000 | 210/256 | 5.812500 | 5.812500 | | — |
| S1 blind | 16 | 2 | 6.000000000 | 256/256 | 6.000000 | 6.000000 | 1.000000 | 6 (equal) |
| S1 blind | 16 | 3 | 6.000000000 | 256/256 | 6.000000 | 6.000000 | | 6 (equal) |
| S2 chem, balanced | 10 | 2 | 19.957890936 | 1494/2048 | 19.957891 | 19.957891 | 1.108521 | |
| S2 chem, balanced | 10 | 3 | 22.123736485 | **33/1024** | 22.161551 | 22.343816 | | |
| S3 chem, projection sizes | 9 | 2 | 18.036971267 | 2048/2048 | 18.036971 | 18.036971 | 1.033836 | |
| S3 chem, projection sizes | 10 | 3 | 18.647263938 | **5/1024** | 18.759861 | 19.078404 | | |
| S4 chem, free sizes | 10 | 2 | 16.941837482 | 676/2048 | 16.941837 | 16.944911 | 0.386151 | |
| S4 chem, free sizes | 10 | 3 | 6.542110309 | 50/512 | 6.542110 | 6.635373 | | |

N_sense: 16 (n=2, all scenarios); 64 (n=3, S1/S2); 61 (n=3, S3/S4, SGC stops removed from the graph).
S3 class sizes: doublet X1 2, X2 4, X3 1, X4 3, X5 2, X6 1, X7 1, X9 1, X10 1; triplet X1 8, X2 14,
X3 5, X4 10, X5 4, X6 4, X7 8, X8 2, X9 4, X10 2. S3 references: the projection codes themselves
score 18.937 (doublet) and 21.269 (triplet), i.e. 5.0% and 14.1% above the size-matched optima found.

## How far each number can be trusted

- **Certified in practice (every or nearly every restart agrees):** all of S1, S2 n=2, S3 n=2.
  S1 at k=4 and k=16 equals the analytic 3·log4(k) at both lengths (subcube partitions).
- **Well supported:** S4 n=2 (676/2048, p5 = best) and S4 n=3 (50/512, p5 = best).
- **Upper bounds only, not certified:** S2 n=3 (33/1024) and especially S3 n=3 (5/1024). An earlier
  development run with a 5x shorter schedule reported 18.647757 for S3 n=3, above the value here, so
  that landscape has near-degenerate optima that restarts only just resolve. The runs at 512 and
  1024 restarts shared seeds, so their agreement on 18.647263938 is not independent confirmation.
  Treat these two as "best found" and the corresponding ratios (1.1085, 1.0338) as resting on them.

## Observations (no interpretation beyond the numbers)

1. Chemistry-blind, k a power of 4: the per-residue optimum is identical at n=2 and n=3 (ratio 1).
   The k=8 case also gives ratio 1. The only S1 gap is k=10, where balanced sizes are forced to be
   uneven differently (16 = 6x2 + 4x1, 64 = 4x7 + 6x6).
2. With chemistry and balanced sizes (S2), the ratio 1.1085 is nearly the same as the blind k=10
   ratio 1.1071, so most of the S2 gap already appears without chemistry.
3. With projection-code sizes (S3), the triplet optimum is 3.4% above the doublet optimum. The
   comparison is k=10 vs k=9 on different size vectors, so the two are not the same optimization problem.
4. With free sizes (S4) both optima collapse to one large class (X8, His) plus nine single-codon
   classes (7+9x1 at n=2; 52+9x1 at n=3). The 0.386 ratio reflects that nine forced singletons are
   9/16 of the doublet space but 9/61 of the triplet space. It is a normalisation effect, not an
   architectural advantage.

## Method and deviations from the spec

- Simulated annealing with incremental delta evaluation. Initial T = mean |delta| over 2000 random
  moves, geometric cooling by 1e4 over the run. Each restart ends with steepest descent over all
  pairwise swaps (plus all legal recolors for S4) until no move improves by more than 1e-12. The
  reported D is recomputed from scratch on the final labeling. 48 worker processes (the machine was
  shared: load average 85-145 on 64 cores).
- Schedule length (moves per restart): S1 1000·N², S2/S3 5000·N², S4 10000·N².
- Restarts: 256 for S1; 2048 for n=2 chemistry scenarios; 1024 for S2/S3 at n=3; 512 for S4 at n=3.
  This is more than the 256 minimum because a probe on S2 n=3 showed hits saturate with schedule
  length (0/48 at 1000·N², 5/48 at 5000·N², 4/48 at 20000·N²), so only more restarts help.
  The 64-restart budget fallback never fired.
- **Deviation, S4 move set:** the spec asks for single-codon recolors. With recolors alone, a
  single-codon class can never move: recoloring its codon would empty the class, and relocating it
  needs a two-step uphill path. The optimum is mostly such classes, so recolor-only annealing froze
  (development run: S4 n=2 1/2048 hits at 16.9418; S4 n=3 1/512 at 6.5910). S4 therefore mixes 50%
  recolors with 50% size-preserving swaps, and steepest descent covers both. That gave 676/2048
  and 50/512 hits, and a lower S4 n=3 optimum (6.5421).
- Wall time: the final run took 12.6 min (`elapsed_seconds` 754 in RESULTS.json). Two earlier
  development runs (1.6 min, 7.9 min) and the probe (about 1 min) are superseded and their numbers are
  not reported except where quoted above.

# Literature comparison — Dirichlet energy and noise distortion in the genetic code

Date: 2026-04-21
Purpose: establish whether the two metrics used in our pipeline have prior art,
in response to Paudyal's Comments 72 and 83.

## Our definitions (what the code actually computes)

From `src/metrics/dirichlet.py` and `run_publication_controls.py::compute_sgc_metrics`:

**E (Dirichlet energy):**

    E = Σ_{(u,v) ∈ edges(G)} ‖f(u) − f(v)‖²

- G = Hamming-1 graph on 4³ = 64 codons (261 edges; 263 sense-only)
- f(c) = property vector of the amino acid encoded by codon c (8-dim, PCA-reduced, z-scored)
- Sum over **undirected** edges (each pair counted once)
- Distance is **squared** L2
- Edges touching stops are skipped

This is the standard graph Laplacian quadratic form: E = fᵀ L f where L is
the unweighted Hamming-1 graph Laplacian.

**D (noise distortion, canonical version used in the manuscript):**

    D = (1/|C|) × Σ_{c ∈ C} Σ_{c' ∈ N(c), nonsyn, sense} ‖f(c) − f(c')‖

- |C| = number of sense codons (61 for SGC)
- N(c) = Hamming-1 neighbors of c
- Distance is **unsquared** L2
- No `1/|N(c)|` per-codon averaging (Kiran's prompt referenced `1/|N(c)|` but our code does not normalize by neighbor count)
- Each undirected edge is counted **twice** (once from each endpoint)
- Synonymous and stop neighbors contribute 0 (filtered explicitly, but would be 0 anyway)
- No weighting by codon usage or mutation probability

## Paper 1 — Tlusty 2007, *J. Theor. Biol.* 249:331–342 (arXiv:1007.4122)

**Distortion (Eq. 1):**

    H_ED = Σ_{α→i→j→β} P_αijβ C_αβ = Σ_{α,i,j,β} P_α E_αi R_ij D_jβ C_αβ

where
- P_α = prior probability amino-acid α is required
- E_αi = encoder matrix (α → codon i)
- R_ij = reading/misreading matrix (codon i read as j)
- D_jβ = decoder matrix (codon j → amino-acid β)
- C_αβ = chemical distance between amino-acids α and β (**unsquared**)

For a deterministic SGC (E and D delta-functions), uniform priors, and
uniform Hamming-1 misreading weights, this reduces to

    H_ED ∝ Σ_{(i,j) ∈ Hamming-1 edges} C_{a(i), a(j)}

which is our D up to (i) overall normalization and (ii) directed- vs undirected-
edge convention. **Same functional form, same distance (unsquared), same
summation pattern.**

**Graph Laplacian (Appendix B, Eq. B1):**

    Δ_ij = −R_ij        (i ≠ j)
    Δ_ii = Σ_{j ≠ i} R_ij

This is the **weighted** graph Laplacian with edge weights equal to the
misreading probabilities R_ij. Our unweighted Hamming-1 Laplacian is the
special case R_ij = const for all Hamming-1 pairs and 0 otherwise.

**Dirichlet-energy-like form (Eq. 3, second-order expansion near non-coding
state):**

    F ~ Σ_{α,β,i,j} (T δ_ij δ_αβ − r²_ij c_αβ) e_αi e_βj

where r = −(Δ − Δ_C) and r² is the **squared** normalized misreading matrix.
This is a quadratic form on the order parameter e_αi that involves the
Laplacian (through r²). In standard notation, fᵀ L f (unweighted Hamming-1)
is the **R_ij = const, no-code-matrix, single-distance-channel** limit of
Tlusty's construction.

## Paper 2 — Tlusty 2010, *Phys. Life Rev.* 7:362–376 (arXiv:1007.3906)

**Error-load (Eq. 1):**

    error-load = Σ_{i,j,α,β} r_ij p_iα p_jβ c_αβ

Same as Tlusty 2007 Eq. 1 but with p_iα as a probabilistic code matrix (p_iα
is the probability that codon i encodes amino-acid α). The distance c_αβ is
unsquared. The accompanying text explicitly says "The error-load measures
the smoothness of the code" (line 216) — confirming this is the distortion
measure and that "smoothness" is the intended semantic interpretation.

**Fitness (Eq. 4):**

    fitness = −error-load + w_D × diversity − w_C × cost

The error-load enters with a negative sign; fitness maximization = error-load
minimization, which is the same optimization our null model tests the SGC
against (we test whether the SGC's error-load is lower than random codes).

**Graph Laplacian:** defined identically to Tlusty 2007 Eq. B1.

## Paper 3 — Radványi & Kun 2021, *Sci. Rep.* 11:10963

**Distortion (Eq. 1):**

    D = Σ_{i,j} P(c_i) × P(Y = c_j | X = c_i) × d(aa_i, aa_j)

where
- P(c_i) = **codon usage** probability in a given organism
- P(Y = c_j | X = c_i) = **Kimura 2-parameter mutation** probability
  (transitions vs transversions, κ = ti/tv ratio = 2.5 in their analysis)
- d(aa_i, aa_j) = physicochemical distance, **unsquared**, with
  d(aa, aa) = 0 (synonymous mutations contribute 0)
- Single-property distances: d is computed separately for hydropathy (D_Hyd),
  polar requirement (D_Pol), molecular volume (D_Vol), isoelectric point
  (D_pI)

Citations 33, 34 in the paper are Tlusty's work — they explicitly say
"the information theoretic concept of distortion... is used to estimate the
average effect" of mutations. **Same object as Tlusty's Eq. 1**, specialized
to (i) organism-specific codon usage priors rather than uniform, (ii) Kimura
mutation model rather than uniform Hamming-1 misreading, (iii) scalar
property distances rather than vector.

They **do not** compute a separate Dirichlet energy or Laplacian quadratic
form. Their analysis is entirely in distortion-space.

## Head-to-head comparison

| Object | Tlusty 2007/2010 | Radványi & Kun 2021 | Our code |
|:-:|:-:|:-:|:-:|
| Distortion formula | Σ r_ij p_iα p_jβ c_αβ | Σ P(c_i) P(c_j\|c_i) d(aa_i, aa_j) | (1/\|C\|) Σ Σ ‖f(c)−f(c')‖ |
| Distance type | Unsquared | Unsquared | Unsquared |
| Distance dimension | 1D (single property) | 1D (single property, 4 analyses) | Multi-D (8-PC) |
| Codon prior | Uniform (implicit) | Empirical codon usage | Uniform |
| Mutation weights | r_ij (weighted misreading) | Kimura ti/tv ratio | Unit weights on Hamming-1 |
| Graph Laplacian formula | Δ_ij = −R_ij (off-diag) | Not used | Δ_ij = −1 if Hamming-1, else 0 |
| Quadratic smoothness form | f^T L² f (in phase-transition analysis) | Not used | f^T L f (direct, as our E) |

## Plain-language verdict

1. **Our "Dirichlet energy" E** is the standard graph Laplacian quadratic
   form f^T L f on the unweighted Hamming-1 codon graph. This is the
   simplest, unweighted instance of the construction Tlusty uses. Tlusty's
   analysis works with the *weighted* Laplacian (edge weights = misreading
   probabilities) and in the phase-transition analysis uses the *squared*
   Laplacian. Our E is not a new metric — it is the direct, unweighted
   specialization of an object in common use in this literature since at
   least Tlusty 2007.

2. **Our "noise distortion" D** is the same object as Tlusty's "error-load"
   (Eq. 1 in both 2007 and 2010) and Radványi & Kun's "distortion" (Eq. 1
   in 2021), specialized to: (i) deterministic codon-to-amino-acid mapping,
   (ii) uniform codon priors, (iii) Hamming-1 unit-weight neighbor edges,
   (iv) multi-dimensional property vector with unsquared L2 distance, (v)
   per-codon average rather than per-pair sum. The functional form and
   semantic interpretation are Tlusty's; our specific parameter choices
   and normalization are different but not conceptually novel.

3. **What is new in our work (relative to these three papers):**
   - The null model: degeneracy-preserving shuffle of amino-acid labels
     across sense codons, preserving per-amino-acid codon counts. Tlusty's
     null is a thermodynamic ensemble of code matrices at temperature T;
     Radványi & Kun has no explicit null — they test correlations against
     environmental variables. Our null directly tests "is the SGC more
     error-robust than random codes with matched degeneracy?"
   - The 2 × 2 factorial (alphabet × architecture) with RAA10-projected
     doublet codes as the `n = 2` control. This is not in any of the three.
   - Running 10⁶ samples and computing empirical z-scores rather than
     asymptotic / analytical phase-transition locations.
   - Using a PCA-reduced 8-dimensional property vector (from 22 descriptors)
     for the distance. Tlusty uses abstract c_αβ; Radványi uses four
     individual scalar properties. Neither uses a multi-dim PCA space.

## Caveats

- I am reporting what the published papers **say**, not evaluating whether
  Tlusty's formulation is the first in the literature. Earlier work
  (Freeland & Hurst 1998, Freeland et al. 2000, Archetti 2004) may have
  applied related distortion-style objects. I did not read those.

- Radványi & Kun cite Tlusty (references 33, 34 in their paper) as the
  source of "the information theoretic concept of distortion." Their
  framing is that they are **applying** Tlusty's framework, not inventing a
  new one. This is the correct pattern for our paper to follow.

- The distinction between "squared" (as in the Dirichlet energy) and
  "unsquared" (as in the distortion) is a real and legitimate design
  choice, not a novelty. Both appear in Tlusty's papers — squared in the
  phase-transition analysis, unsquared in the primary distortion definition.

## Suggested action

Adjust the Introduction / Methods to cite Tlusty 2007 and/or 2010 as the
source of the rate-distortion framework for genetic code analysis and
Radványi & Kun 2021 as a recent application of the same framework to a
different biological question (mutational robustness vs. environment).
Position the contribution of this paper as: (i) a degeneracy-preserving
null model that lets us compute empirical z-scores and p-values, (ii) a
2 × 2 factorial that isolates codon architecture from alphabet size, and
(iii) a 10⁶-sample Monte Carlo implementation with a multi-dimensional
PCA-based property distance. Do **not** claim to have introduced Dirichlet
energy or noise distortion to this problem.

I have not made any citation decisions beyond listing what the three
papers contain. The Introduction rewrite is yours to do.

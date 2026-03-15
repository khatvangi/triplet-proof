# Article type: Original Article

# Title

Triplet architecture enables deep error-minimization in the genetic code

# Running title

Triplet architecture and genetic code robustness

# Authors

[Author names to be inserted]

# Affiliations

1. [Affiliation 1]
2. [Affiliation 2]

# Corresponding author

[Corresponding author name]  
[Institution]  
[Postal address]  
[Email]

# Keywords

Genetic code; Molecular evolution; Error minimization; Codon architecture; Wobble position; Theoretical evolution

## Abstract

The standard genetic code (SGC) assigns amino acids to trinucleotide codons in a pattern that reduces the impact of single-nucleotide errors, but the relative contributions of amino acid assignment and codon architecture have remained difficult to separate. We quantified error-minimization on the codon Hamming graph using two complementary graph-theoretic metrics: Dirichlet energy, which measures local smoothness of amino acid properties across neighboring codons, and single-step distortion, which measures cumulative property change across one-nucleotide sense-neighbor substitutions. Amino acids were represented in an 8-dimensional physicochemical space derived by principal component analysis of a broad descriptor set. Against 1,000,000 random codes with matched codon-count degeneracy, the SGC outperformed all sampled alternatives on both metrics (Dirichlet z = -14.61; distortion z = -16.57; empirical p < 10^-6 for both). To distinguish codon length from alphabet-size effects, we analyzed a 2 x 2 design comparing a 10-amino-acid doublet projection, the 20-amino-acid triplet code, and a 10-amino-acid triplet control. Moving from doublet to triplet architecture shifted the z-score by 11.0 units for distortion and 10.0 units for Dirichlet energy, whereas expanding from 10 to 20 amino acids at fixed triplet length shifted the z-score by only 1.8 and 1.3 units, respectively. In the SGC, 68.9% of position-3 substitutions are synonymous, compared with 4.4% at position 1 and 0% at position 2; in the 10-amino-acid triplet control, position-3 synonymy rises to 76.5%. These results indicate that triplet architecture provides a high-synonymy dimension that permits substantially deeper error-minimization than doublet projections. More broadly, they support the view that triplet coding is not only sufficient for encoding the modern amino acid repertoire, but also especially favorable for robust error-buffering.

## Introduction

Single-nucleotide errors during translation are an unavoidable cost of decoding mRNA into protein. Their phenotypic severity depends strongly on which amino acid replaces which: some substitutions are nearly neutral, whereas others are severely disruptive. The arrangement of amino acids in the standard genetic code (SGC) places chemically similar amino acids on neighboring codons, so that common translation errors tend to produce less damaging replacements (Woese 1965; Sonneborn 1965).

Whether this pattern reflects direct selection for error-minimization or arises largely as a by-product of the code's evolutionary history has been debated for decades. Haig and Hurst (1991) introduced Monte Carlo randomization of amino acid assignments as a quantitative test of code robustness. Freeland and Hurst (1998) estimated that the SGC outperforms approximately 99.99% of random codes under a polar-requirement metric, and later work confirmed strong but metric-dependent optimization under alternative physicochemical measures (Gilis et al. 2001; Goodarzi et al. 2004; Itzkovitz and Alon 2007). Novozhilov et al. (2007) argued that partial local optimization of a random starting code can yield robustness levels comparable to the SGC, supporting a mixed picture of contingency and selection. Koonin and Novozhilov (2009) further emphasized that the apparent strength of optimization depends on the chosen amino acid property scale.

Two unresolved issues follow from this literature. First, any result based on one or a few hand-picked amino acid property scales remains vulnerable to property-selection bias. Second, nearly all prior tests evaluate optimization within the triplet code itself, rather than asking whether triplet architecture is qualitatively more favorable to robust coding than shorter architectures. Doublet codes provide only 16 codons, whereas triplet codes provide 64 codons and a highly synonymous third position. The evolutionary significance of this additional coding dimension has not been cleanly quantified.

Here we address both issues. We represent amino acids in an 8-dimensional orthogonal property space derived by principal component analysis (PCA) of a broad descriptor set, reducing the degree of freedom associated with manual property choice. We then evaluate the code using two graph-theoretic metrics on the codon Hamming graph: Dirichlet energy, which captures smoothness of the property landscape, and single-step distortion, which captures cumulative change induced by local nucleotide errors. Finally, we compare doublet and triplet constructions in a 2 x 2 factorial design that separates architecture from alphabet size. This analysis shows that triplet architecture supports a much deeper level of error-minimization than doublet projections and that the mechanistic basis of this difference is the emergence of a highly synonymous wobble-like third position, consistent with the classical wobble framework (Crick 1966).

## Materials and methods

### Amino acid property space

Amino acid physicochemical descriptors were compiled from canonical published scales, including AAindex-curated measures (Kawashima and Kanehisa 2000), and implemented in a local property library. We began with 22 descriptors spanning hydrophobicity, volume, polarity, charge, hydrogen-bonding capacity, flexibility, and secondary-structure propensity. After correlation pruning, 20 standardized descriptors were retained. Principal component analysis extracted 8 orthogonal components explaining 97.1% of total variance. Each amino acid was represented by its 8-dimensional PC score vector.

### Codon Hamming graph

For a code of length n, the Hamming graph H(4, n) contains 4^n codons as nodes, with edges connecting codons that differ at exactly one position. Sense codons inherit the property vector of their assigned amino acid; stop codons are excluded from graph-based metric calculations.

### Metrics

Dirichlet energy was computed as

E = sum_{(u,v) in edges} ||f(u) - f(v)||^2

where f(c) is the amino acid property vector assigned to codon c. Lower values indicate a smoother local property landscape.

Single-step distortion was computed as

D = (1 / |C|) sum_{c in C} sum_{c' in N_sense(c)} ||f(c) - f(c')||

where C is the set of sense codons and N_sense(c) is the set of sense codons reachable from c by a single-nucleotide substitution. Synonymous substitutions contribute zero by construction, and substitutions to stop codons are omitted.

### Null model

Random codes were generated by permuting amino acid labels among sense codons while preserving the number of codons assigned to each label in the reference code. Thus, the null model preserves codon-count degeneracy but not codon-block topology. For each permutation, both metrics were computed. Main analyses used 1,000,000 permutations per condition; sensitivity analyses used 100,000 permutations. Empirical one-sided tail probabilities were estimated as (count + 1) / (n + 1), where count is the number of random codes with metric value less than or equal to the observed value.

### Cross-length factorial design

Three conditions were analyzed.

Condition A (10 amino acids, n = 2): a doublet code derived from the SGC position-1/2 structure using a reduced 10-class alphabet following Murphy et al. (2000). For each reduced class, the representative amino acid was chosen as the member closest to the class centroid in 8-PC space.

Condition B (20 amino acids, n = 3): the standard genetic code.

Condition C (10 amino acids, n = 3): a triplet control using the same reduced 10-amino-acid alphabet projected onto the triplet SGC structure.

For all conditions, null permutations preserved codon counts per label, and both D and E were computed. Architecture effects were defined as z_C - z_A. Alphabet effects were defined as z_B - z_C.

### Sensitivity analysis

For the 10-amino-acid triplet control, two alternative representative sets were also tested: the amino acid farthest from the class centroid and a random class member. The z-score range across representative-set choices remained below 0.7 for both metrics.

### Synonymy analysis

Position-specific synonymy was defined as the fraction of all possible single-nucleotide substitutions at a given codon position that preserved the same amino acid identity, or the same reduced-alphabet class in reduced-alphabet analyses.

## Results

### The standard code is extreme under both graph-theoretic metrics

Across 1,000,000 random codes with matched codon-count degeneracy, the SGC outperformed all sampled alternatives on both metrics. For the 20-amino-acid triplet code, single-step distortion yielded z = -16.57 and Dirichlet energy yielded z = -14.61, with zero random codes matching or beating the SGC on either metric. The two metrics were strongly correlated across the null ensemble (r = 0.958), indicating that they capture related, though not identical, aspects of code robustness.

These values place the SGC far beyond the range observed in the random ensemble and support a much stronger signal than the classic one-property analyses of code optimization (Haig and Hurst 1991; Freeland and Hurst 1998).

### Doublet projections are optimized, but only shallowly

The SGC-derived doublet projection, evaluated under the 10-amino-acid reduced alphabet, remained better than most random codes but not exceptionally so. For Condition A, single-step distortion yielded z = -3.80 and Dirichlet energy yielded z = -3.32. The empirical tail probabilities were 1.08 x 10^-4 and 4.66 x 10^-4, respectively. Thus, the doublet projection lies in the null tail, but not beyond the range sampled in one million shuffles.

### Triplet architecture dominates alphabet size

The critical control was the 10-amino-acid triplet construction. Under this condition, single-step distortion yielded z = -14.81 and Dirichlet energy yielded z = -13.34, with zero random codes matching or beating the observed values on either metric.

The factorial decomposition was therefore:

- architecture effect (n = 2 to n = 3, same 10 amino acids): -11.01 z-units for distortion and -10.02 z-units for Dirichlet energy
- alphabet effect (10 to 20 amino acids, same triplet length): -1.76 z-units for distortion and -1.27 z-units for Dirichlet energy
- architecture:alphabet ratio: 6.3:1 for distortion and 7.9:1 for Dirichlet energy

Under both metrics, triplet architecture contributed roughly six- to eight-fold more to optimization depth than alphabet expansion.

### The wobble-like third position is the mechanistic basis

The SGC exhibits strong positional asymmetry in synonymy. In the full 20-amino-acid triplet code, 68.9% of all position-3 substitutions are synonymous, compared with 4.4% at position 1 and 0% at position 2. Across all three positions combined, 24.4% of one-nucleotide substitutions are synonymous. In the 10-amino-acid triplet control, position-3 synonymy rises to 76.5%.

By contrast, in the reduced-alphabet doublet projection, the two positions show only 33.3% and 8.3% synonymy, respectively, with an overall synonymous fraction of 20.8%. Neither position functions as a wobble-like buffer comparable to triplet position 3.

This asymmetry provides a structural explanation for the difference in optimization depth. In triplets, positions 1 and 2 can encode amino acid family structure while position 3 absorbs many local errors as silent substitutions. In doublets, both positions continue to carry substantial discriminating information, leaving no comparably silent axis in codon space.

## Discussion

Three conclusions follow from these results.

First, the standard genetic code is extremely optimized relative to random alternatives when evaluated in a multivariate physicochemical space rather than along a single hand-selected property axis. Because the PCA representation reduces researcher discretion in property choice, the resulting signal is less vulnerable to the criticism that code optimality depends on arbitrary scale selection (Koonin and Novozhilov 2009).

Second, the depth of error-minimization depends strongly on coding architecture. The reduced-alphabet doublet projection is optimized, but only shallowly. The triplet architecture, by contrast, enables access to a part of code space that the random null does not reach in one million samples. This extends the assignment-level perspective of previous work by showing that architectural degrees of freedom matter at least as much as, and here much more than, the number of distinct amino acid labels.

Third, the mechanistic basis of this architectural effect is positional synonymy. The third codon position provides a highly tolerant dimension that converts a large fraction of local mutational or translational perturbations into silent substitutions. This creates a directional structure in the triplet Hamming graph that has no analog in the doublet projection analyzed here.

This manuscript does not claim that the present reduced-alphabet doublet projection is a historical reconstruction of a primordial code, nor that the null model captures every biologically plausible constraint. It also does not establish why triplets evolved historically. Rather, it shows that when matched random ensembles are compared under a consistent property representation, triplet architecture supports a qualitatively deeper level of error-minimization than the doublet projection examined here.

Longer codes were not optimized in this study. Although longer codons could in principle create additional synonymous directions, evaluating that possibility would require an explicit model of longer-code assignment, cellular cost, and evolutionary accessibility. The present result is therefore narrower and cleaner: relative to doublets, triplets cross a threshold in achievable error-minimization depth.

## Acknowledgements

[Add acknowledgements here, if any.]

## Funding

[Add funding information here. If none: "This work received no specific funding."]

## Author contributions

[Add author contribution statement here.]

## Data availability

All analysis outputs supporting this manuscript are contained in the accompanying repository. A public archival link should be inserted here before submission.

## Code availability

All code used for the analyses is contained in the accompanying repository. A public archival link should be inserted here before submission.

## Conflict of interest

The authors declare no conflict of interest.

## References

Crick FHC (1966) Codon-anticodon pairing: the wobble hypothesis. J Mol Biol 19:548-555. https://doi.org/10.1016/S0022-2836(66)80022-0

Freeland SJ, Hurst LD (1998) The genetic code is one in a million. J Mol Evol 47:238-248. https://doi.org/10.1007/PL00006381

Gilis D, Massar S, Cerf NJ, Rooman M (2001) Optimality of the genetic code with respect to protein stability and amino-acid frequencies. Genome Biol 2:research0049.1-research0049.12. https://doi.org/10.1186/gb-2001-2-11-research0049

Goodarzi H, Nejad HA, Torabi N (2004) On the optimality of the genetic code, with the consideration of termination codons. Biosystems 77:163-173. https://doi.org/10.1016/j.biosystems.2004.04.002

Haig D, Hurst LD (1991) A quantitative measure of error minimization in the genetic code. J Mol Evol 33:412-417. https://doi.org/10.1007/BF02103132

Itzkovitz S, Alon U (2007) The genetic code is nearly optimal for allowing additional information within protein-coding sequences. Genome Res 17:405-412. https://doi.org/10.1101/gr.5987307

Kawashima S, Kanehisa M (2000) AAindex: amino acid index database. Nucleic Acids Res 28:374. https://doi.org/10.1093/nar/28.1.374

Koonin EV, Novozhilov AS (2009) Origin and evolution of the genetic code: the universal enigma. IUBMB Life 61:99-111. https://doi.org/10.1002/iub.146

Murphy LR, Wallqvist A, Levy RM (2000) Simplified amino acid alphabets for protein fold recognition and implications for folding. Protein Eng 13:149-152. https://doi.org/10.1093/protein/13.3.149

Novozhilov AS, Wolf YI, Koonin EV (2007) Evolution of the genetic code: partial optimization of a random code for robustness to translation error in a rugged fitness landscape. Biol Direct 2:24. https://doi.org/10.1186/1745-6150-2-24

Sonneborn TM (1965) Degeneracy of the genetic code: extent, nature, and genetic implications. In: Bryson V, Vogel HJ (eds) Evolving genes and proteins. Academic Press, New York, pp 377-397

Woese CR (1965) On the evolution of the genetic code. Proc Natl Acad Sci USA 54:1546-1552. https://doi.org/10.1073/pnas.54.6.1546

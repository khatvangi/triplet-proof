# RULES (frozen before any download or computation) — 2026-09-23

## Descriptor sourcing
Every descriptor must be one of:
(a) an AAindex1 entry, identified by accession, values read programmatically from the
    official AAindex1 file (https://www.genome.jp/ftp/db/community/aaindex/aaindex1),
    saved under sources/ with sha256 and download date;
(b) computed deterministically from side-chain structure by a stated rule (RDKit on
    canonical SMILES of the 20 free amino acids; SMILES and RDKit version recorded);
(c) a value from a named primary table saved under sources/ with DOI and table location.
Nothing hand-entered without (c).

## Concepts, fixed order (the 22 of src/io/aa_props_lib.py)
1 hydropathy (Kyte-Doolittle); 2 side-chain volume (Zamyatnin); 3 polarity (Grantham);
4 isoelectric point; 5 side-chain pKa; 6 net charge at pH 7; 7 side-chain H-bond donors;
8 acceptors; 9 aromaticity; 10 ring count; 11 rotatable bonds; 12 molar refractivity;
13 bulkiness (Zimmerman); 14 Chou-Fasman alpha; 15 beta; 16 turn; 17 solvent accessibility;
18 flexibility; 19 helical moment; 20 strand moment; 21 hydrophobic surface; 22 side-chain mass.

## Selection rule
For each concept choose the AAindex entry whose title/author matches the named source.
If the named source is not in AAindex, choose the closest AAindex entry measuring the same
concept and record why. If none exists and no (b)/(c) route exists, DROP the concept and
record it. Choices are made from AAindex titles/authors and chemistry only — never by
looking at downstream results. Choices are written to table_S2_new.csv before any
statistic is computed with the new matrix.
Structure-derived (b) counts are computed on the SIDE CHAIN (atoms other than the
backbone N, CA-H excluded, C(=O)O) using a stated SMARTS/RDKit rule recorded in the code.

## Pipeline (identical to the original otherwise)
- z-score each descriptor across the 20 AAs (ddof=0);
- walk concepts in the fixed order; drop a descriptor if |r| >= 0.95 with an EARLIER kept
  descriptor (record which);
- PCA by SVD of the z-matrix;
- primary: smallest k with cumulative variance >= 97%; also fixed k = 8;
- primary matrix = each retained PC z-scored (ddof=0) ("rescaled"); secondary = unscaled PCs.
- d(a,b) = Euclidean distance between AA vectors.

## Statistics (fixed before running)
Gate first: with OLD vectors (reviewer_response/aa_props_export.csv) the code must reproduce
SGC D = 23.9223597 and doublet-projection D = 18.9372988 to 1e-9 before any new-matrix run.
Then (1) codon-shuffle null 1e6 seed 20260924 (A, B, C; plus equidistant d=1);
(2) block null: A exhaustive 9!, C exhaustive 10!, B 1e7 seed 20260924 (Ser one block),
z and p = fraction strictly better; old-vector C must give p(D)=0.34, p(E)=0.32, A p(D)=0.0068;
(3) mean d over 190 pairs, SGC mean d per non-synonymous single-nt misreading (uniform, Katoh
0.98/1/1.39), medians under the two nulls (1e6), best block code in 1e6 + pairwise-swap descent;
(4) Landerer Data S4 split per decode_split.py; (5) usage-weighted Dw (uniform + proteome).
Qualitative verdict per statistic is stated as one line and judged Y/N on sign/tail side and
the conventional 0.05 threshold, not on exact magnitude.

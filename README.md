# triplet-proof

Monte Carlo framework demonstrating that the standard genetic code (SGC) is
statistically optimized for error-minimization, and that triplet codon
architecture — not amino acid alphabet size — is the primary enabler.

Companion code for:
**"Triplet architecture enables deep error-minimization of the genetic code"**
(submitted to *Journal of Molecular Evolution*)

## Key findings

1. **SGC optimization** (p < 10⁻⁶): the SGC outperforms all 1,000,000 random
   codes with matched degeneracy on both metrics — noise distortion
   (z = −16.6) and Dirichlet energy (z = −14.6).
2. **Triplet architecture crosses an optimization threshold that doublet
   architecture does not**, at either alphabet size. The 10-amino-acid doublet
   projection stays inside its null tail (z = −3.8; the best random doublet
   reaches z = −4.0, within one unit). Both triplet conditions lie beyond the
   reach of 10⁶ random codes: 10 amino acids z = −14.8, 20 amino acids
   z = −16.6.
3. **Mechanism**: position 3 (wobble) absorbs 69–77% of single-nucleotide
   errors as synonymous substitutions, creating error-buffering corridors that
   doublet codes structurally cannot support.

> **A note on comparing z-scores across conditions.** The null distributions at
> n = 2 and n = 3 operate over code spaces of different size and degeneracy, so
> their variance structure differs and Δz values are **not strictly
> commensurable across architectures**. `publication_controls.json` does record
> a `factorial_summary` block with differences and ratios; treat those as
> directional indicators only. The manuscript interprets the result in terms of
> whether each condition crosses its null envelope, not by z-score subtraction.

## Reproduce

```bash
# prerequisites: python 3.10+, numpy, pandas, scipy, scikit-learn, networkx, matplotlib

# 1. build amino acid property space (22 descriptors → 8 PCs, 97.1% variance)
python build_aa_props.py

# 2. phase-1 sanity check (~2 min, single process)
python run_phase1_demo.py

# 3. full 1M triplet run (~2 min on 64 cores)
python run_phase2_batch.py --n 1000000 --workers 32 \
  --out results/phase2_triplet.auto.json

# 4. doublet (n=2) with 10-class reduced alphabet
python run_phase2_n2n4.py --n-bases 2 --n 1000000 --workers 32 --raa raa10 \
  --out results/phase2_doublet.raa10.json

# 5. 2×2 factorial + sensitivity analysis (generates all manuscript numbers)
python run_publication_controls.py --n-null 1000000 --workers 32

# 6. verify synonymy numbers
python verify_synonymy.py

# 7. run tests
python -m pytest tests/ -v

# 8. generate figures (each emits .pdf + .png)
python figures/fig1.py
python figures/fig2.py
python figures/fig3.py
python figures/fig4.py
python figures/figS1_raw_distributions.py
```

Two further scripts reproduce the supplementary analyses:

```bash
# extra random seeds for the Table S1 3-seed sensitivity protocol
python run_sensitivity_extra_seeds.py

# preemptive defense: squared Laplacian fᵀL²f and Boltzmann-weighted distortion
python run_preempt_metrics.py --n-null 100000 --workers 32

# build Table S1 (CSV + LaTeX)
python build_table_s1.py
```

SLURM batch scripts (`slurm_triplet.sbatch`, `slurm_doublet.sbatch`,
`slurm_quadruplet.sbatch`) request 32 CPUs / 64 GB RAM.

## Directory structure

```
triplet-proof/
├── src/                            # core library
│   ├── io/
│   │   ├── codon_io.py             # load codon_table.csv → DataFrame
│   │   ├── aa_props_io.py          # feature loading, correlation pruning, PCA
│   │   └── aa_props_lib.py         # 22 amino acid properties (canonical values)
│   ├── sims/
│   │   ├── codon_graph.py          # Hamming graph on 4^n nodes
│   │   └── random_codes.py         # degeneracy-preserving shuffle
│   ├── metrics/
│   │   ├── dirichlet.py            # E = Σ_edges ||p(u) - p(v)||²
│   │   └── mi.py                   # mutual information by codon position
│   └── receiver/
│       └── thermo_noise.py         # noise distortion (Boltzmann softmax)
├── tests/
│   ├── test_graph.py               # Hamming graph construction
│   ├── test_phase1_compat.py       # codon loader + phase-1 compatibility
│   └── test_synonymy.py            # position-specific synonymy regression
├── figures/
│   ├── fig1.py                     # framework: Hamming graph, PCA space, schematic
│   ├── fig2.py                     # SGC vs null: E, D histograms + joint scatter
│   ├── fig3.py                     # factorial overlay + position-specific synonymy
│   ├── fig4.py                     # wobble mechanism: triplet vs doublet neighborhoods
│   └── figS1_raw_distributions.py  # per-condition raw nulls (supplementary)
├── results/                        # canonical outputs (see below)
├── data/processed/
│   ├── aa_props.parquet            # 8 PCA columns (97.1% variance)
│   └── aa_props_meta.json          # PCA metadata
├── codon_table.csv                 # standard genetic code (64 codons)
├── aa_props.csv                    # 6 raw properties (legacy, phase-1 only)
├── build_aa_props.py               # build parquet from 22-property library
├── run_phase1_demo.py              # phase-1 sanity check
├── run_phase2_batch.py             # 1M triplet Monte Carlo
├── run_phase2_n2n4.py              # doublet/quadruplet Monte Carlo
├── run_publication_controls.py     # 2×2 factorial + sensitivity
├── verify_synonymy.py              # position-specific synonymy computation
├── supplementary/                  # Table S1, Table S2 (CSV + LaTeX)
├── share/                          # portable CSV export of the raw nulls
│   └── triplet-proof-figure-data/  # 3 x 1M (D, E) rows + summary tables
├── scripts/
│   └── export_for_colleague.py     # regenerates the share/ export
├── manuscript_JME.md               # manuscript draft
├── supplementary_tables_JME.md     # Tables S1, S2
└── conftest.py                     # pytest path setup
```

## Metrics

**Dirichlet energy** measures smoothness of amino acid properties across the
codon Hamming graph:

    E = Σ_{(u,v) ∈ edges} ||f(u) − f(v)||²

**Noise distortion** measures total property change across all non-synonymous
single-nucleotide substitutions, averaged per sense codon:

    D = (1/|C|) Σ_{c ∈ C} Σ_{c' ∈ N_sense(c)} ||f(c) − f(c')||

Lower values = more robust code. The SGC scores lower than all 1M random codes
on both metrics.

## Null model

Random codes are generated by permuting amino acid labels among individual
sense codons while preserving the number of codons assigned to each amino acid
(degeneracy structure). Stop codons remain fixed. This null preserves degeneracy
but not codon-block topology.

Percentiles use Laplace smoothing: `(count + 1) / (n + 1)`.

## Canonical results

| File | Condition | Key result |
|------|-----------|------------|
| `phase1_demo.json` | triplet, 6 raw features, 100K | sanity check |
| `phase2_triplet.auto.json` | triplet, 8 PCs, 1M | 0/1M on E and N |
| `phase2_doublet.raa10.json` | doublet, 10 classes, 1M | E: 0.13%, N: 0.018% |
| `phase2_quadruplet.auto.json` | quadruplet, 8 PCs, 1M | 0/1M on E and N |
| `publication_controls.json` | 2×2 factorial, 1M per condition | **definitive manuscript numbers** |
| `sensitivity_extra_seeds.json` | random-representative seeds 1, 2 | feeds Table S1's 3-seed mean |
| `preempt_metrics.json` | fᵀL²f and Boltzmann-weighted D | robustness to metric convention |
| `sgc_baselines.json` | SGC metric values for figures | — |

### Raw null distributions

`results/null_distributions.npz` (46 MB) and the three extracted `.npy` arrays
are **gitignored**, so they are not part of this repository or its Zenodo
archive. They are fully reproducible with `run_publication_controls.py`.

For inspection without rerunning the Monte Carlo, `share/triplet-proof-figure-data/`
contains the same distributions as plain CSV (three files of 1,000,000 rows,
each row a `(D, E)` pair), together with the PCA property table, the SGC
baselines, the synonymy fractions, and the z-score summary. Every figure and
every z-score in the manuscript can be regenerated from that directory alone.

## Tests

```bash
python -m pytest tests/ -v
# 8 tests: graph construction, codon loading, synonymy regression, doublet class count
```

## Citation

If you use this code or data, please cite the archived release. Machine-readable
metadata is in `CITATION.cff`; GitHub renders it as a "Cite this repository"
button in the sidebar.

<!-- ZENODO-DOI: replace this block after the first release is archived -->
> **DOI pending.** This repository is archived on Zenodo; the concept DOI
> (which always resolves to the most recent version) will be inserted here
> once the first release is deposited.

## License

Dual-licensed by content type:

| Content | License |
|---------|---------|
| Source code (`src/`, `tests/`, `figures/*.py`, `scripts/`, `run_*.py`, `build_*.py`, `*.sbatch`) | [MIT](LICENSE) |
| Manuscript, figures, tables, derived data | [CC BY 4.0](LICENSE-CC-BY-4.0.md) |

The amino acid descriptors compiled in `src/io/aa_props_lib.py` derive from the
[AAindex database](https://www.genome.jp/aaindex/) and remain subject to its
terms of use.

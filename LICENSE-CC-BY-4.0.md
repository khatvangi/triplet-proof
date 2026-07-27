# Creative Commons Attribution 4.0 International (CC BY 4.0)

SPDX-License-Identifier: CC-BY-4.0

Copyright (c) 2026 Hannah E. Rebbeck, J. Paudyal, Boggavarapu Kiran

## Scope

This repository is dual-licensed by content type. The `LICENSE` file contains
the MIT License, which governs all **software**: Python source under `src/`,
`tests/`, `figures/`, and `scripts/`, the top-level `run_*.py` and `build_*.py`
drivers, and the SLURM batch scripts. `LICENSE` is kept as verbatim MIT text so
that automated license detectors identify it correctly; this file records the
scope split.

The CC BY 4.0 terms below cover the non-software content of this repository:

- `manuscript_JME.md`, `supplementary_materials_JME.md`,
  `supplementary_tables_JME.md`, `supplementary_materials_JME.pdf`
- `cover_letter_JME.md`, `submission_checklist_JME.md`
- `literature_comparison.md`, `rebuild_notes.md`
- All rendered figures: `figures/*.pdf`, `figures/*.png`, `figures/*.svg`,
  and the accompanying `figures/*_caption.txt` files
- All tabular and derived data: `results/*.json`, `supplementary/table_S*.csv`,
  `supplementary/table_S*.tex`, `share/**`, `codon_table.csv`, `aa_props.csv`,
  and `data/processed/**`

Source code is licensed separately under the MIT License; see `LICENSE`.

## Terms

You are free to:

- **Share** — copy and redistribute the material in any medium or format
- **Adapt** — remix, transform, and build upon the material for any purpose,
  including commercially

under the following terms:

- **Attribution** — You must give appropriate credit, provide a link to the
  license, and indicate if changes were made. You may do so in any reasonable
  manner, but not in any way that suggests the licensor endorses you or your
  use.
- **No additional restrictions** — You may not apply legal terms or
  technological measures that legally restrict others from doing anything the
  license permits.

## Attribution

Please cite the associated work. See `CITATION.cff` for machine-readable
citation metadata.

## Full legal code

The complete license text is available at:

- Human-readable summary: https://creativecommons.org/licenses/by/4.0/
- Full legal code: https://creativecommons.org/licenses/by/4.0/legalcode

## Third-party data

The amino acid physicochemical descriptors in `src/io/aa_props_lib.py` are
canonical published values compiled from the AAindex database
(https://www.genome.jp/aaindex/). AAindex is subject to its own terms of use;
attribute AAindex when redistributing those descriptor values.

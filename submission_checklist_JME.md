# Journal of Molecular Evolution Submission Checklist

## Manuscript package

- [ ] Final manuscript file completed: [manuscript_JME.md](/storage/kiran-stuff/triplet-proof/manuscript_JME.md)
- [ ] Title finalized
- [ ] Running title finalized
- [ ] Author names inserted
- [ ] Affiliations inserted
- [ ] Corresponding author details inserted
- [ ] Keywords reviewed
- [ ] Abstract checked against journal word limits

## Back matter

- [ ] Acknowledgements completed or removed
- [ ] Funding statement completed
- [ ] Author contributions completed
- [ ] Conflict of interest statement confirmed
- [ ] Data availability statement finalized with public archive link
- [ ] Code availability statement finalized with public archive link

## Figures and tables

- [ ] Figure files exported at submission quality
- [ ] Figure numbering matches manuscript callouts
- [ ] Table numbering matches manuscript callouts
- [ ] Table S1 / Table S2 actually prepared if cited
- [ ] Figure legends checked against final figures

## References

- [ ] DOI links verified
- [ ] Reference list checked against in-text citations
- [ ] Any missing classic genetic-code references added if desired
- [ ] Journal title abbreviations reviewed if required by production

## Submission metadata

- [ ] Cover letter completed: [cover_letter_JME.md](/storage/kiran-stuff/triplet-proof/cover_letter_JME.md)
- [ ] Suggested reviewers added
- [ ] Opposed reviewers added if needed
- [ ] Article type set to Original Article
- [ ] Short description / significance statement prepared if portal requests it

## Reproducibility package

- [ ] Canonical result files archived:
  - `results/phase1_demo.json`
  - `results/phase2_triplet.auto.json`
  - `results/phase2_doublet.raa10.json`
  - `results/phase2_quadruplet.auto.json`
  - `results/publication_controls.json`
  - `results/alphabet_control.json`
- [ ] Repository snapshot archived (Zenodo or equivalent)
- [ ] Commit or release tag recorded in manuscript notes

## Final technical checks

- [ ] `python -m pytest -q` passes
- [ ] No stale quantitative claims remain after final prose edits
- [ ] Manuscript numbers match:
  - `results/publication_controls.json`
  - `results/alphabet_control.json`
- [ ] README and manuscript are not in conflict on canonical outputs

## Remaining manual items in this repo

- [ ] Fill author/affiliation metadata in [manuscript_JME.md](/storage/kiran-stuff/triplet-proof/manuscript_JME.md)
- [ ] Fill reviewer suggestions in [cover_letter_JME.md](/storage/kiran-stuff/triplet-proof/cover_letter_JME.md)
- [ ] Insert public code/data archive links
- [ ] Prepare the cited supplementary tables

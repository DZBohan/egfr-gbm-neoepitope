# What is done and what is next

## Done

- Step 2, antigen processing and presentation at 9 and 13 residues. `docs/step2-9mer.md`,
  `docs/step2-13mer.md`, tables in `results/iedb_9mer/` and `results/iedb_13mer/`.
- Step 3, expression integration. `docs/step3-expression.md`, plots in `figures/`.
- Step 4, structural immunogenicity, both models. `docs/step4-immunogenicity.md`, table in
  `results/iedb_variant_comparison/`, plots in `figures/`.

- Step 5, benchmark validation against TCIA. `docs/step5-tcia-validation.md`, export in
  `results/tcia_gbm_egfr/`.

All five parts of the exercise are written up. What follows is optional.

## Open points

- ICERFIRE ran with its default `total_gene_tpm` of 6.071. EGFR is amplified and highly
  expressed in this cohort, so those percentiles are a floor. The web server has no field
  for expression, so correcting it needs the API or a local install.
- Step 5 compared our predictions against TCIA's predictions, not against observed spectra.
  The exercise calls the TCIA table mass spectrometry pull-down data and it is not. A real
  validation would need an immunopeptidomics dataset from HLA-A2-positive glioblastoma.
- `VVMGENNTL` is the one peptide both pipelines agree on, and it is also the one whose
  substitution sits at P1. Whether that is a presentable neo-epitope or a self-like surface
  is the open biological question, and it is an experiment, not a calculation.

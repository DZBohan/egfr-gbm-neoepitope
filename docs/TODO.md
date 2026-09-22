# What is done and what is next

## Done

- Step 2, antigen processing and presentation at 9 and 13 residues. `docs/step2-9mer.md`,
  `docs/step2-13mer.md`, tables in `results/iedb_9mer/` and `results/iedb_13mer/`.
- Step 3, expression integration. `docs/step3-expression.md`, plots in `figures/`.
- Step 4, structural immunogenicity, both models. `docs/step4-immunogenicity.md`, table in
  `results/iedb_variant_comparison/`, plots in `figures/`.

## Next: step 5, LC-MS/MS benchmark validation

The Cancer Immunome Atlas, glioblastoma, Neoantigens table, parent gene EGFR. From that
list:

1. Count how many of the step 4 peptides appear in observed pull-down data, true positives,
   and how many do not, false positives.
2. Compute precision, TP / (TP + FP), for the step 4 list.
3. Say whether that precision justifies testing these peptides experimentally.
4. Describe an in vivo approach to testing immunogenicity in glioblastoma, from current
   literature.

The step 4 candidate list to check against, best ICERFIRE percentile first: `TCPAVVMGE`
15.9, `AVVMGENNT` 26.1, `CPAVVMGEN` 29.1, `PAVVMGENN` 40.0, `VVMGENNTL` 42.3, `CVKTCPAVV`
50.8, `VKTCPAVVM` and `KTCPAVVMG` 52.5, `HCVKTCPAV` 56.0.

## Open points

- ICERFIRE ran with its default `total_gene_tpm` of 6.071. EGFR is amplified and highly
  expressed in this cohort, so those percentiles are a floor. The web server has no field
  for expression, so correcting it needs the API or a local install.
- Nothing from step 4 reaches a percentile usually treated as promising. Step 5 is where
  that gets checked against observed data rather than argued from the scores.

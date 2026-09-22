# Step 4, structural immunogenicity: what the two models said

Peptide Variant Comparison, Class I, HLA-A\*02:01. Input is `results/pairs_g598v.tsv`, the
nine 9-mer windows spanning residue 598 from step 2, wild type as peptide A and G598V as
peptide B on the same line. Three models were run: NetMHCpan 4.1 EL for context, Class I
pMHC Immunogenicity with the default mask of positions 1, 2 and the C terminus, and
Neo-Epitope Immunogenicity (ICERFIRE 1.0). Full output in
`results/iedb_variant_comparison/peptide_table.csv`.

The top ten presented peptides from step 2 are not in this input, and the reason is worth
stating: they are identical between the wild-type and mutant protein, so their pairs would
be two copies of the same sequence and every difference would be zero by construction. The
nine windows below are the only peptides where a pair can differ at all.

## Scores

| # | peptide A | peptide B | imm A | imm B | delta imm | EL %ile A | EL %ile B | ICERFIRE | ICERFIRE %ile |
|---|---|---|---|---|---|---|---|---|---|
| 1 | HCVKTCPAG | HCVKTCPAV | -0.20305 | -0.20305 | 0.00000 | 73 | 7.6 | 0.1550 | 55.97 |
| 2 | CVKTCPAGV | CVKTCPAVV | -0.04106 | -0.03674 | +0.00432 | 11 | 13 | 0.1619 | 50.75 |
| 3 | VKTCPAGVM | VKTCPAVVM | 0.03710 | 0.04334 | +0.00624 | 66 | 41 | 0.1596 | 52.45 |
| 4 | KTCPAGVMG | KTCPAVVMG | -0.02642 | -0.01946 | +0.00696 | 31 | 26 | 0.1596 | 52.45 |
| 5 | TCPAGVMGE | TCPAVVMGE | -0.02077 | -0.01357 | +0.00720 | 45 | 39 | 0.2265 | 15.86 |
| 6 | CPAGVMGEN | CPAVVMGEN | 0.00880 | 0.01624 | +0.00744 | 95 | 100 | 0.1952 | 29.14 |
| 7 | PAGVMGENN | PAVVMGENN | -0.00584 | -0.00344 | +0.00240 | 100 | 100 | 0.1771 | 40.00 |
| 8 | AGVMGENNT | AVVMGENNT | -0.04529 | -0.04529 | 0.00000 | 46 | 20 | 0.2011 | 26.11 |
| 9 | GVMGENNTL | VVMGENNTL | 0.08573 | 0.08573 | 0.00000 | 0.78 | 0.93 | 0.1737 | 42.28 |

Scatter plots, peptide B against peptide A, are `figures/step4_immunogenicity_score_scatter.jpg`
and `figures/step4_icerfire_el_rank_scatter.jpg`. In the first, the points sit on or beside
the diagonal. In the second, six pairs fall below it, meaning the mutant is better ranked,
two sit above it and one ties, with the largest move going from 80 to 10.1. That second plot
is ICERFIRE's internal EL rank, not its final percentile.

## Are the mutant peptides predicted to be more immunogenic?

By the pMHC immunogenicity model, no. Three pairs score exactly the same for wild type and
mutant, and the other six gain between 0.0024 and 0.0074, which is a rounding-level change
against scores that span 0.29 across this set. Nothing here is a change in predicted
immunogenicity.

The three exact zeros are the informative part, and they are predictable from position
alone. `HCVKTCPAG` to `HCVKTCPAV` puts the substitution on the C terminus, `AGVMGENNT` to
`AVVMGENNT` puts it at P2, and `GVMGENNTL` to `VVMGENNTL` puts it at P1. Those are the three
masked positions. The model masks them on the assumption that they anchor into the groove
and are not contacted by the TCR, so a substitution there cannot change its score. The zeros
therefore follow from the mask, and they are a statement about the model rather than a
measurement of what a receptor can reach: P1 in particular can be solvent exposed in real
HLA-A2 complexes.

ICERFIRE does give the mutants scores, and they are modest. The best is `TCPAVVMGE` at the
15.9th percentile, then `AVVMGENNT` at 26.1 and `CPAVVMGEN` at 29.1. The rest sit between
40 and 56. ICERFIRE publishes no validated cut-off for a promising neo-epitope, so these
are useful for ranking the nine against each other and not for declaring any of them
positive.

The two models also disagree about which peptide is the candidate. Step 2 singled out
`HCVKTCPAV`, whose EL percentile improves from 73 to 7.6, the largest presentation gain in
the set. ICERFIRE puts it at 55.97, near the bottom of these nine.

## Why the two models differ

**They see different positions.** The pMHC model masks P1, P2 and the C terminus and scores
only the outward-facing residues. ICERFIRE takes the wild-type and mutant peptide as a pair
and reports a difference for the three substitutions the first model cannot see. Its
published consensus model applies masking of its own, so the contrast is in what each model
masks and what else it uses, not in one having no mask at all. Either way, this is what
produces the three zero rows in one model and not the other.

**ICERFIRE includes presentation, the pMHC model does not.** The pMHC model scores the
peptide surface presented to the TCR and says nothing about whether the peptide is presented
at all. ICERFIRE combines a binding term with its other features, which is why a pair whose
EL rank moves from 80 to 10.1 is treated differently by the two.

**ICERFIRE reports similarity to self, and it tracks the disagreement here.** `HCVKTCPAV`
has an icore similarity score of 0.968, among the highest in the set, though not the maximum:
`VVMGENNTL` is higher at 0.971. Its presentation improves sharply while its ICERFIRE
percentile stays at 56. A neo-epitope closely resembling its own wild type is expected to
meet a repertoire already tolerised against it, which a binding-only view cannot express.
The similarity column is an output of the wrapper rather than a demonstrated term in the
published consensus model, so this is a consistent reading of the numbers, not a verified
mechanism.

**ICERFIRE scores the icore, not the submitted 9-mer.** Rows 3 and 4 return identical
ICERFIRE values because both collapse to the same icore, `KTCPAVVM`. Two different 9-mer
windows can therefore be indistinguishable to this model while the pMHC model treats them as
separate peptides.

**The expression term is at its default.** `total_gene_tpm` is 6.071 in every row, the value
used when none is supplied. ICERFIRE weights expression, and step 3 showed EGFR in this
cohort is amplified and highly expressed, so the real value for these tumours is far above
the default. The ICERFIRE scores here are therefore conservative, and a cohort-appropriate
TPM would raise them. The web server offers no field for it, so this cannot be corrected
from the interface used here, and the percentiles above should be read as a floor rather
than an estimate.

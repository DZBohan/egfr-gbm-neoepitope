# egfr-gbm-neoepitope

Working scripts for a proteogenomics course exercise on EGFR in glioblastoma: expression
integration from cBioPortal, antigen presentation prediction, and neo-epitope
immunogenicity.

The scripts are gene- and study-agnostic. EGFR in `gbm_cptac_2021` is the default because
that cohort's protein data is mass spectrometry rather than RPPA, and the exercise's final
question concerns mass spectrometry.

## What the analysis found

Write-ups are in `docs/`, one per part of the exercise, each naming the tool settings it
used and the file its numbers come from.

`docs/step2-9mer.md` and `docs/step2-13mer.md`, antigen processing and presentation. At both
peptide lengths the ten best-presented peptides are identical between wild type and G598V,
because none of them comes near residue 598. Of the windows that do span it, one moves
sharply at each length, and it is the one that places the substitution on the C terminus:
glycine gives the F pocket of HLA-A\*02:01 no anchor, valine does. The two runs also differ
from each other in a way the score columns explain: proteasome and TAP scores do not depend
on peptide length, while MHC binding does, because the groove holds nine to ten residues and
a 13-mer must bulge.

`docs/step3-expression.md`, expression integration. EGFR mRNA and protein are far higher in
mutated samples, and that contrast is confounded: every mutated sample with copy-number data
is also amplified, and amplified samples without a mutation sit at the same protein level.
The difference tracks copy number, not mutation status. The same file sets out why an
abundant mutant protein is still hard to confirm by mass spectrometry.

`docs/step4-immunogenicity.md`, structural immunogenicity. The pMHC immunogenicity model
reports no meaningful difference between wild-type and mutant peptides, and returns exactly
equal scores for the three windows that put the substitution on a masked anchor position.
ICERFIRE scores the same pairs differently, and the file works through the five reasons the
two models disagree.

`docs/step5-tcia-validation.md`, benchmark comparison. One of the nine peptides,
`VVMGENNTL`, is listed for four glioblastoma patients in The Cancer Immunome Atlas, giving a
precision of 1/9 on exact sequence identity, or 3/9 if a longer listed peptide containing
ours counts, which are answers to two different questions and both reported. The file also
sets out why a precision over nine overlapping windows on one substitution measures the
enumeration more than the biology, and why that atlas is a comparison against other
predictions rather than the mass spectrometry validation the exercise describes.

## Scripts

`scripts/01_expression.py` pulls mRNA, protein, copy number and mutations for one gene from
the cBioPortal API, compares mutated against non-mutated samples, and stratifies by
amplification.

`scripts/02_peptides.py` translates a wild-type/mutant CDS pair, reports the substitution it
actually finds, enumerates the k-mers spanning it, and writes the paired input file the
Peptide Variant Comparison server takes.

```bash
./scripts/01_expression.py --gene EGFR --entrez 1956 --study gbm_cptac_2021 --change G598V
./scripts/02_peptides.py --wt data/egfr_wt.fa --mut data/egfr_mut.fa --pairs results/pairs.tsv
```

## Why these scripts check what they check

Both exist because of errors made while working the same exercise on KRAS in colorectal
cancer, and both checks earn their place.

**The CDS was submitted to the binding predictor instead of the protein.** A, T, C and G are
all valid amino acid letters, so the tool accepted 567 nucleotides as a 567-residue protein
and returned 1,118 rows of nonsense without an error. `02_peptides.py` translates when the
input looks like DNA and prints residue numbering so the substitution can be confirmed
before anything is submitted.

**A file named G12D contained G13D.** KRAS carries glycines at both 12 and 13, so the
protein sequence looks correct either way. The script reports the substitution it finds
rather than the one the file name claims. On the original input it printed `G13D`, which is
what led to the correction.

## The amplification confounder

For EGFR in glioblastoma, mutated samples show much higher mRNA and protein than
non-mutated ones. Stratifying shows why: every mutated sample in this cohort is also
amplified, and amplified samples without a mutation sit at nearly the same protein level as
amplified samples with one. The difference tracks copy number, not mutation status.
`01_expression.py` prints that stratification alongside the raw contrast so the comparison
is not read at face value.

KRAS in colorectal cancer behaves in the opposite way and is a useful contrast: almost no
amplification, and no expression difference between mutated and non-mutated samples, which
is what an activating substitution that changes protein function rather than abundance
should look like.

## Anchor positions

The Class I pMHC immunogenicity model masks positions 1, 2 and the C terminus by default,
because those residues anchor into the MHC groove and are not contacted by the TCR. A
substitution landing there cannot change that model's score, while a neo-epitope model
taking the wild-type/mutant pair still registers it. `02_peptides.py` marks which windows
put the substitution in a masked position.

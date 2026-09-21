# egfr-gbm-neoepitope

Working scripts for a proteogenomics course exercise on EGFR in glioblastoma: expression
integration from cBioPortal, antigen presentation prediction, and neo-epitope
immunogenicity.

The scripts are gene- and study-agnostic. EGFR in `gbm_cptac_2021` is the default because
that cohort's protein data is mass spectrometry rather than RPPA, and the exercise's final
question concerns mass spectrometry.

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

# Step 3, expression integration: what the cohort shows

`gbm_cptac_2021`, EGFR, entrez 1956. That cohort is chosen because its protein data is mass
spectrometry rather than RPPA, and the question at the end of this step is about mass
spectrometry. Numbers below are from `./scripts/01_expression.py --gene EGFR --entrez 1956
--study gbm_cptac_2021 --change G598V`. The screenshot of the Plots tab, Protein vs mRNA
coloured by amino-acid change, sits alongside this file.

99 samples with both protein and mRNA, 96 with copy number, 17 mutated. The mutation
spectrum is G598V x6, R222C x2, A289D x2, and one each of G665D, H304Y, V774M, P772_H773dup,
D770_V774dup.

## Mutated samples look far higher, at both levels

| contrast | protein median | mRNA median |
|---|---|---|
| any mutation (n=17) | +2.055 | 8,019,678 |
| no mutation (n=82) | -0.098 | 527,698 |
| G598V (n=6) | +2.101 | 7,980,513 |
| rest (n=93) | +0.137 | 733,154 |

Protein and mRNA agree across the cohort, Pearson r = 0.817, so this is not a measurement
artefact at one level.

## Stratifying by copy number shows it is amplification, not mutation

47 of 96 samples are amplified, and every mutated sample with copy-number data is amplified,
15 of 15.

| group | n | protein median |
|---|---|---|
| amplified and mutated | 15 | +2.146 |
| amplified, no mutation | 32 | +1.835 |
| neither | 49 | -0.648 |

Amplified samples without a mutation sit at nearly the same protein level as amplified
samples with one, +1.835 against +2.146. The gap that matters is against the unamplified
group at -0.648. The contrast in the first table is confounded: mutation status in this
cohort is a near-perfect proxy for amplification, so the apparent effect of the mutation is
the effect of copy number. Nothing here says G598V raises EGFR expression.

## What this means for finding the mutant peptide by mass spectrometry

The favourable part is abundance. These tumours carry a heavily amplified, highly expressed
EGFR, which is close to the best case for detecting any EGFR peptide in a shotgun run.

The rest is harder, and abundance does not fix it.

**The mutant peptide has to be in the search database.** A search against a reference
proteome cannot report `TCPAVVMGENNTLVWK`, because that sequence is not in it. The spectrum
is either unassigned or misassigned to something else. Identifying it needs a
sample-specific database built from that tumour's own sequencing, which is why proteogenomic
studies build one per patient.

**Only one tryptic peptide covers the site.** Trypsin cuts after K and R, so residue 598
falls in `TCPAGVMGENNTLVWK`, residues 594 to 609, wild type, and `TCPAVVMGENNTLVWK` mutant.
There is no second peptide to fall back on if that one ionises poorly or is missed, so
coverage of this particular residue is all-or-nothing in a tryptic digest. A second protease
would give a different window over the same residue.

**Glycine to valine is a +42.047 Da shift.** That is close to trimethylation and to
acetylation at +42.011, and a search that allows those modifications can place them on the
wrong residue of the wild-type peptide and score it plausibly. The mutant and the modified
wild type need to be distinguished on fragment ions, not precursor mass.

**Both alleles are present.** Amplification means many copies, not exclusively mutant ones,
so the wild-type peptide is in the sample too, and usually in excess. The mutant is measured
against an abundant near-identical background, which is a harder detection problem than the
raw protein level suggests.

These points are why the binding predictions from step 2 are a starting point rather than an
answer. `HCVKTCPAV` scoring well on HLA-A\*02:01 says the peptide could be presented if it
is produced. Whether it is produced, and whether anyone can see it in a real run, is the
question this step raises and does not settle.

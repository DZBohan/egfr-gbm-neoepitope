# Step 5, benchmark validation against TCIA

The Cancer Immunome Atlas, disease GBM, Neoantigens table, parent gene EGFR, exported
filtered. The export is `results/tcia_gbm_egfr/neoantigens_gbm_egfr.tsv`: 111 rows, 78
distinct peptides, 35 patients, peptide lengths 8 to 11.

## What this table is, and what it is not

The exercise describes this as mass spectrometry pull-down data. It is not. The TCIA
neoantigen table is computed from TCGA exome and RNA sequencing: somatic variants are called
and candidate peptides are scored by binding prediction. No peptide in it was observed by
mass spectrometry, and the HLA-alleles column here is `NA` for every row.

The comparison below is therefore between our predictions and someone else's predictions,
made from a different cohort with a related but not identical pipeline. That is worth doing,
because agreement means two independent runs picked the same peptide out of the same protein,
but it is not experimental confirmation, and a precision computed this way is not a measured
precision.

## True positives, false positives, precision

The step 4 list is the nine mutant 9-mers spanning residue 598.

| step 4 peptide | in TCIA exactly | related TCIA entry |
|---|---|---|
| HCVKTCPAV | no | `GPHCVKTCPAV` contains it |
| CVKTCPAVV | no | none |
| VKTCPAVVM | no | none |
| KTCPAVVMG | no | none |
| TCPAVVMGE | no | none |
| CPAVVMGEN | no | none |
| PAVVMGENN | no | none |
| AVVMGENNT | no | `AVVMGENNTL` contains it |
| VVMGENNTL | yes, 4 patients | `AVVMGENNTL`, `VVMGENNTLV`, `VVMGENNTLVW` |

On exact sequence identity: TP = 1, FP = 8, precision = 1/9 = 0.111.

Counting a peptide as recovered when TCIA lists a longer peptide containing it: TP = 3,
FP = 6, precision = 3/9 = 0.333.

The first number is the honest one for a peptide-level claim, because a 9-mer and the 11-mer
containing it are different molecules that bind differently and are presented separately.
The second is the fair number for a residue-level claim, that the site was flagged. Both are
reported here rather than choosing the flattering one.

`VVMGENNTL` is the peptide both methods agree on, and it appears in four TCIA patients:
TCGA-06-0174, TCGA-12-0616, TCGA-19-2620 and TCGA-28-5213.

## Does this precision support experimental testing?

Read as a hit rate for the list as submitted, 0.111 does not. Eight of nine peptides would be
synthesised for nothing.

But the list was never a set of nine independent candidates. All nine are windows over one
substitution, they overlap heavily, and the point of enumerating them was to find which
register presents the site best, not to propose nine separate epitopes. Precision over an
overlapping window set mostly measures how many windows were enumerated: adding 13-mers to
the submission would have lowered it without changing anything biological.

The result that carries weight is convergence on a single peptide. `VVMGENNTL` leads the
nine mutant windows on processing total score, ranked 33 of 1,202 in step 2, it has the best
EL percentile of the nine at 0.93, and a separately built pipeline on a different cohort
listed it for four patients. That is worth testing. The other eight are not, and this comparison is a
reasonable way to have reached that conclusion.

Two cautions belong with it. First, `VVMGENNTL` puts the substitution at P1, an anchor
position, which is why the pMHC immunogenicity model in step 4 returned the same score for it
and for the wild type. That equality follows from the model's mask, not from a measurement of
what a T cell receptor can reach. Second, the substitution does not improve binding here. The
wild-type `GVMGENNTL` has the better EL percentile of the two, 0.78 against 0.93, and the
better predicted affinity, 1,509 nM against 1,808 nM. What the mutant leads on is the
processing total score. So the candidate is a peptide whose presentation is no better than
its wild-type counterpart and whose mutated residue points into the groove, and the
repertoire has plausibly been tolerised against the near-identical wild-type surface.

## How common is this mutation in the TCGA cohort

Six peptides in the export are compatible with the G598V sequence: `GPHCVKTCPAV`,
`KTCPAVVM`, `AVVMGENNTL`, `VVMGENNTL`, `VVMGENNTLV` and `VVMGENNTLVW`. They come from seven
patients, TCGA-06-0174,
TCGA-12-0616, TCGA-12-0619, TCGA-19-2620, TCGA-28-1747, TCGA-28-5213 and TCGA-32-1982, out of
the 35 patients in this EGFR table.

So this substitution recurs across patients rather than being private to one, which is what
makes a shared-epitope approach worth considering at all. The export carries no variant
coordinates or HLA assignments, so this is a count within a gene-filtered table, not a
prevalence estimate for TCGA glioblastoma. It also matches step 3, where G598V was 6 of 17
mutated samples in the CPTAC cohort.

## An in vivo approach to testing immunogenicity

The peptide is restricted to HLA-A\*02:01, which no mouse expresses, so the first requirement
is an HLA-A2 transgenic host, and the standard choice is the HHD line carrying a chimeric
HLA-A2.1/H-2Db molecule on a β2m-null background.

A workable sequence, from most to least direct:

1. Immunise HHD mice with `VVMGENNTL`, and in parallel with the wild-type `GVMGENNTL` as the
   control the whole question turns on. Read out by IFN-γ ELISpot and by tetramer staining
   against both peptides. The specific question is not whether the mutant is immunogenic but
   whether responses to it cross-react with wild type, because a response that does is a
   response to a self peptide.
2. Confirm the peptide is actually presented by the tumour, not only immunogenic when
   injected. Immunopeptidomics on an HLA-A2-positive, G598V-positive glioma line, or on the
   HHD tumour model, is the step that answers what step 3 raised and no prediction can.
3. Challenge with an orthotopic HLA-A2-positive glioma expressing EGFR G598V in vaccinated
   versus control mice. Intracranial placement matters here rather than being a formality:
   the brain's immune environment is the reason many systemically immunogenic glioma antigens
   have failed.

The clinical context to keep in view is that single-antigen EGFR vaccination in glioblastoma
has already been tried and failed. Rindopepimut targeted EGFRvIII and did not improve
survival in the phase 3 ACT IV trial, with antigen loss in recurrent tumours part of the
explanation. Current work has moved to personalised multi-epitope vaccination: a phase 1 DNA
vaccine trial used up to 40 neoantigens per patient and reported peripheral T cell expansion
in treated patients, and a phase Ib neoantigen-pulsed dendritic cell trial reported a median
progression-free survival of 16.2 months. The implication for a single peptide such as
`VVMGENNTL` is that it is a component of a multi-epitope design, not a candidate for
single-antigen vaccination on its own.

Sources for the clinical statements above:
[phase 1 multivalent neoantigen DNA vaccine](https://www.nature.com/articles/s43018-026-01163-w),
[phase Ib neoantigen-pulsed dendritic cells](https://pmc.ncbi.nlm.nih.gov/articles/PMC13469091/),
[review of glioblastoma immunotherapy trials](https://pmc.ncbi.nlm.nih.gov/articles/PMC12081429/).

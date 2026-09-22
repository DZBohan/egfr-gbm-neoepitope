# Step 2, 13-mers: what the run showed

Same submission as the 9-mer run except for peptide length: IEDB T Cell Prediction Class I,
HLA-A\*02:01, NetMHCpan 4.1 BA with MHC-I Processing (Basic Processing Predictions,
immunoproteasome), both proteins submitted together so `seq #` 1 is wild type and 2 is
G598V. 1,198 thirteen-mers per sequence.

## The top ten are identical between wild type and G598V

| # | peptide | start | proteasome | tap | mhc | total |
|---|---|---|---|---|---|---|
| 1 | FLSNMSMDFQNHL | 172 | 1.584 | 0.322 | -1.150 | 0.757 |
| 2 | RLLGICLTSTVQL | 776 | 1.598 | 0.536 | -2.109 | 0.025 |
| 3 | FLSLQRMFNNCEV | 48 | 1.082 | 0.091 | -1.559 | -0.386 |
| 4 | ALLAALCPASRAL | 13 | 1.419 | 0.517 | -2.454 | -0.519 |
| 5 | RMFNNCEVVLGNL | 53 | 1.450 | 0.603 | -2.697 | -0.644 |
| 6 | QLITQLMPFGCLL | 787 | 1.398 | 0.481 | -2.687 | -0.808 |
| 7 | SIATGMVGALLLL | 645 | 1.325 | 0.439 | -2.575 | -0.810 |
| 8 | TAGAALLALLAAL | 6 | 1.380 | 0.346 | -2.587 | -0.861 |
| 9 | VLIALNTVERIPL | 89 | 1.221 | 0.460 | -2.593 | -0.912 |
| 10 | QLMPFGCLLDYVR | 791 | 1.247 | 0.690 | -2.873 | -0.936 |

As at 9-mers, none of them comes near residue 598, so the substitution cannot alter them.
The lengthening changes which regions dominate the top of the list (172 and 776 are new,
813 and 649 drop out), but not the conclusion.

## The peptides that do span residue 598

| window | wild type | rank | total | G598V | rank | total | delta |
|---|---|---|---|---|---|---|---|
| 586-598 | IDGPHCVKTCPAG | 1120 | -4.622 | IDGPHCVKTCPAV | 510 | -3.466 | +1.156 |
| 587-599 | DGPHCVKTCPAGV | 581 | -3.584 | DGPHCVKTCPAVV | 484 | -3.410 | +0.174 |
| 588-600 | GPHCVKTCPAGVM | 573 | -3.574 | GPHCVKTCPAVVM | 454 | -3.331 | +0.244 |
| 589-601 | PHCVKTCPAGVMG | 979 | -4.372 | PHCVKTCPAVVMG | 975 | -4.361 | +0.010 |
| 590-602 | HCVKTCPAGVMGE | 1038 | -4.468 | HCVKTCPAVVMGE | 1036 | -4.459 | +0.009 |
| 591-603 | CVKTCPAGVMGEN | 822 | -4.074 | CVKTCPAVVMGEN | 841 | -4.115 | -0.041 |
| 592-604 | VKTCPAGVMGENN | 910 | -4.241 | VKTCPAVVMGENN | 908 | -4.233 | +0.008 |
| 593-605 | KTCPAGVMGENNT | 811 | -4.051 | KTCPAVVMGENNT | 790 | -3.998 | +0.053 |
| 594-606 | TCPAGVMGENNTL | 149 | -2.330 | TCPAVVMGENNTL | 144 | -2.318 | +0.012 |
| 595-607 | CPAGVMGENNTLV | 372 | -3.132 | CPAVVMGENNTLV | 371 | -3.135 | -0.002 |
| 596-608 | PAGVMGENNTLVW | 297 | -2.860 | PAVVMGENNTLVW | 255 | -2.738 | +0.122 |
| 597-609 | AGVMGENNTLVWK | 475 | -3.393 | AVVMGENNTLVWK | 398 | -3.187 | +0.206 |
| 598-610 | GVMGENNTLVWKY | 22 | -1.323 | VVMGENNTLVWKY | 19 | -1.216 | +0.107 |

Ranks are out of 1,198.

The same window wins, for the same reason. `IDGPHCVKTCPAG` to `IDGPHCVKTCPAV` moves 610
places, the only window that moves more than a few hundred, and it is the one window that
puts the substitution on the C terminus. Extending the peptide four residues at the N
terminus does not change what the F pocket sees, so the 9-mer result reappears at 13.

The other twelve windows move very little, and `GVMGENNTLVWKY` again ranks near the top for
both sequences because the substitution sits at P1.

## Reading these against the 9-mers

Rank positions are comparable between the two runs, absolute scores are not: a 13-mer is
scored against other 13-mers, and NetMHCpan handles 13-mers less well than 9-mers, which is
also the biological expectation. Most presented peptides on HLA-A\*02:01 are 9 to 10
residues long. The point of the 13-mer run is that the C-terminal anchor effect is not an
artefact of one peptide length.

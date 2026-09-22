# Step 2, 9-mers: what the run showed

IEDB T Cell Prediction Class I, HLA-A*02:01, NetMHCpan 4.1 BA with MHC-I Processing
(Basic Processing Predictions, immunoproteasome). Both the wild-type and G598V proteins
were submitted together, so `seq #` 1 is wild type and 2 is G598V. 1,202 nine-mers per
sequence.

## The top ten are identical between wild type and G598V

| # | peptide | start | proteasome | tap | mhc | total |
|---|---|---|---|---|---|---|
| 1 | YLLNWCVQI | 813 | 1.369 | 0.218 | -0.601 | 0.986 |
| 2 | ALLALLAAL | 10 | 1.380 | 0.563 | -1.024 | 0.919 |
| 3 | GMVGALLLL | 649 | 1.325 | 0.429 | -1.231 | 0.523 |
| 4 | QLMPFGCLL | 791 | 1.398 | 0.463 | -1.365 | 0.497 |
| 5 | SQYLLNWCV | 811 | 1.255 | 0.259 | -1.332 | 0.181 |
| 6 | RMFNNCEVV | 53 | 1.127 | 0.321 | -1.342 | 0.107 |
| 7 | LLLLLVVAL | 654 | 1.821 | 0.502 | -2.233 | 0.090 |
| 8 | YQDPHSTAV | 1125 | 1.018 | 0.073 | -1.300 | -0.209 |
| 9 | YVLIALNTV | 88 | 1.064 | 0.182 | -1.516 | -0.271 |
| 10 | ALLAALCPA | 13 | 0.739 | -0.131 | -0.926 | -0.318 |

Every one of these sits at residue 10, 13, 53, 88, 649, 654, 791, 811, 813 or 1125. None
comes within eighty residues of position 598, so the substitution cannot alter them. On
this allele the mutation does not change which peptides are best presented.

## The peptides that do span residue 598

| window | wild type | rank | total | G598V | rank | total | delta |
|---|---|---|---|---|---|---|---|
| 590-598 | HCVKTCPAG | 1028 | -4.412 | HCVKTCPAV | 280 | -2.824 | +1.588 |
| 591-599 | CVKTCPAGV | 189 | -2.474 | CVKTCPAVV | 210 | -2.543 | -0.069 |
| 592-600 | VKTCPAGVM | 439 | -3.285 | VKTCPAVVM | 351 | -3.006 | +0.279 |
| 593-601 | KTCPAGVMG | 757 | -3.928 | KTCPAVVMG | 718 | -3.838 | +0.090 |
| 594-602 | TCPAGVMGE | 988 | -4.351 | TCPAVVMGE | 974 | -4.318 | +0.033 |
| 595-603 | CPAGVMGEN | 924 | -4.234 | CPAVVMGEN | 951 | -4.286 | -0.052 |
| 596-604 | PAGVMGENN | 1097 | -4.546 | PAVVMGENN | 1041 | -4.431 | +0.115 |
| 597-605 | AGVMGENNT | 918 | -4.221 | AVVMGENNT | 752 | -3.910 | +0.311 |
| 598-606 | GVMGENNTL | 34 | -1.168 | VVMGENNTL | 33 | -1.142 | +0.026 |

Ranks are out of 1,202.

`HCVKTCPAG` to `HCVKTCPAV` moves 748 places. The substitution falls on the peptide's C
terminus, which anchors into the F pocket of HLA-A\*02:01, and that pocket prefers a
hydrophobic residue. Glycine gives it nothing; valine gives it an anchor. The eight other
windows place the substitution away from an anchor and move very little.

`GVMGENNTL` and `VVMGENNTL` are the only pair to reach the top fifty, and they rank
together because the substitution sits at P1 in both.

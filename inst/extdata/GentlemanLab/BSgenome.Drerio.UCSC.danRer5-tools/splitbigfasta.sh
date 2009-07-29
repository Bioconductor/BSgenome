#!/bin/bash
#
# Run with:
#   ./splitbigfasta.sh danRer5.fa
# Then test the result with:
#   cat chr[1-9].fa chr??.fa Zv7_NA.fa Zv7_scaffold.fa chrM.fa > merged.fa
#   cmp merged.fa danRer5.fa

set -e  # Exit immediately if a simple command exits with a non-zero status

# Note that in 'head -n NH $1 | tail -n NT' the NH and NT values
# were obtained from R with:
#   library(Biostrings)
#   fi <- fasta.info("danRer5.fa")
#   fi25 <- fi[1:25]
#   NT <- fi25 %/% 50 + 2
#   NH <- cumsum(NT)
head -n  1124095 $1 > chr1.fa
head -n  2211431 $1 | tail -n 1087336 > chr2.fa
head -n  3470057 $1 | tail -n 1258626 > chr3.fa
head -n  4322107 $1 | tail -n  852050 > chr4.fa
head -n  5729536 $1 | tail -n 1407429 > chr5.fa
head -n  6913551 $1 | tail -n 1184015 > chr6.fa
head -n  8318793 $1 | tail -n 1405242 > chr7.fa
head -n  9447929 $1 | tail -n 1129136 > chr8.fa
head -n 10477749 $1 | tail -n 1029820 > chr9.fa
head -n 11325342 $1 | tail -n  847593 > chr10.fa
head -n 12217671 $1 | tail -n  892329 > chr11.fa
head -n 13168147 $1 | tail -n  950476 > chr12.fa
head -n 14239096 $1 | tail -n 1070949 > chr13.fa
head -n 15369555 $1 | tail -n 1130459 > chr14.fa
head -n 16302145 $1 | tail -n  932590 > chr15.fa
head -n 17363560 $1 | tail -n 1061415 > chr16.fa
head -n 18409770 $1 | tail -n 1046210 > chr17.fa
head -n 19395399 $1 | tail -n  985629 > chr18.fa
head -n 20319025 $1 | tail -n  923626 > chr19.fa
head -n 21449600 $1 | tail -n 1130575 > chr20.fa
head -n 22370748 $1 | tail -n  921148 > chr21.fa
head -n 23150386 $1 | tail -n  779638 > chr22.fa
head -n 24078148 $1 | tail -n  927762 > chr23.fa
head -n 24884016 $1 | tail -n  805868 > chr24.fa
head -n 25541542 $1 | tail -n  657526 > chr25.fa
tail -n 333 $1 > chrM.fa

# NH and NT found after several try-and-fail attempts
head -n 27902547 $1 | tail -n 2361005 > Zv7_NA.fa
head -n 28818810 $1 | tail -n  916263 > Zv7_scaffold.fa


#!/bin/bash
#
# Run with:
#   ./split_human_b36_female.sh human_b36_female.fa
# Then test the result with:
#   cat chr[1-9].fa chr??.fa chrX.fa chrM.fa supercontigs.fa > merged.fa
#   cmp merged.fa human_b36_female.fa

set -e  # Exit immediately if a simple command exits with a non-zero status

# Note that in 'head -n NH $1 | tail -n NT' the NH and NT values
# were obtained from R with:
#   library(Biostrings)
#   fi <- fasta.info("human_b36_female.fa")
#   fi24 <- fi[1:24]
#   NT <- fi24 %/% 60L + 2L
#   NT[13] <- NT[13] - 1L
#   NH <- cumsum(NT)
head -n  4120830 $1 > chr1.fa
head -n  8170017 $1 | tail -n 4049187 > chr2.fa
head -n 11495049 $1 | tail -n 3325032 > chr3.fa
head -n 14682935 $1 | tail -n 3187886 > chr4.fa
head -n 17697234 $1 | tail -n 3014299 > chr5.fa
head -n 20545569 $1 | tail -n 2848335 > chr6.fa
head -n 23192594 $1 | tail -n 2647025 > chr7.fa
head -n 25630509 $1 | tail -n 2437915 > chr8.fa
head -n 27968398 $1 | tail -n 2337889 > chr9.fa
head -n 30224645 $1 | tail -n 2256247 > chr10.fa
head -n 32465520 $1 | tail -n 2240875 > chr11.fa
head -n 34671347 $1 | tail -n 2205827 > chr12.fa
head -n 36573731 $1 | tail -n 1902384 > chr13.fa
head -n 38346542 $1 | tail -n 1772811 > chr14.fa
head -n 40018859 $1 | tail -n 1672317 > chr15.fa
head -n 41499315 $1 | tail -n 1480456 > chr16.fa
head -n 42812229 $1 | tail -n 1312914 > chr17.fa
head -n 44080850 $1 | tail -n 1268621 > chr18.fa
head -n 45144379 $1 | tail -n 1063529 > chr19.fa
head -n 46184980 $1 | tail -n 1040601 > chr20.fa
head -n 46967387 $1 | tail -n  782407 > chr21.fa
head -n 47795579 $1 | tail -n  828192 > chr22.fa
head -n 50377476 $1 | tail -n 2581897 > chrX.fa
head -n 50377754 $1 | tail -n     278 > chrM.fa

tail -n 213994 $1 > supercontigs.fa


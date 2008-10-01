#!/bin/bash
#
# Run with:
#   ./splitbigfasta.sh assembly4_chromosomes.fa
# Then test the result with:
#   cat Group?.fa Group??.fa > merged.fa
#   cmp merged.fa assembly4_chromosomes.fa

set -e  # Exit immediately if a simple command exits with a non-zero status

tail -n       +1 $1 | head -n 427631 > Group1.fa
tail -n  +427632 $1 | head -n 229604 > Group2.fa
tail -n  +657236 $1 | head -n 194595 > Group3.fa
tail -n  +851831 $1 | head -n 175097 > Group4.fa
tail -n +1026928 $1 | head -n 207154 > Group5.fa
tail -n +1234082 $1 | head -n 253417 > Group6.fa
tail -n +1487499 $1 | head -n 183558 > Group7.fa
tail -n +1671057 $1 | head -n 188419 > Group8.fa
tail -n +1859476 $1 | head -n 158329 > Group9.fa
tail -n +2017805 $1 | head -n 180610 > Group10.fa
tail -n +2198415 $1 | head -n 207458 > Group11.fa
tail -n +2405873 $1 | head -n 161559 > Group12.fa
tail -n +2567432 $1 | head -n 146669 > Group13.fa
tail -n +2714101 $1 | head -n 142525 > Group14.fa
tail -n +2856626 $1 | head -n 145140 > Group15.fa
tail -n +3001766 $1 | head -n 101043 > Group16.fa
tail -n +3102809 $1 > GroupUn.fa


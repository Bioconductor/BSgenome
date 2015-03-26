#!/bin/bash

set -e  # Exit immediately if a simple command exits with a non-zero status

#for seed in BSgenome.*-seed; do
for seed in `ls BSgenome.*-seed | grep -v '\.masked-seed'`; do
	mv $seed ${seed}.old
	cat ${seed}.old | sed -r "s/^Version: 1\.4\.1$/Version: 1.4.2/" > $seed
	rm ${seed}.old
done


#!/bin/bash

set -e  # Exit immediately if a simple command exits with a non-zero status

for seed in BSgenome.*-seed; do
	mv $seed ${seed}.old
	cat ${seed}.old | sed -r "s/^Version: 1\.3\.16$/Version: 1.3.17/" > $seed
	rm ${seed}.old
done


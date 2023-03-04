#!/bin/bash

for i in 3 6 9 12; do
	sed "s/@/$i/g" Snakefile > snakeFile
	snakemake --snakefile snakeFile --cores $i --use-conda -r -p
done

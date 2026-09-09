#!/bin/bash

samples=(
    A1 B1 M1 N1
    A2 B2 M2 N2
    A3 B3 M3 N3
    A4 B4 M4 N4
    A5 B5 M5 N5
)

for sample in ${samples[@]}; do
	sbatch b2b_kallisto.slurm $sample
done

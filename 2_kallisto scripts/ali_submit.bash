#!/bin/bash

samples=(
    K1 K2 K5 K6
    L1 L2 L5 L6
    M1 M2 M5 M6
    N1 N2 N5 N6
    P1 P2 P5 P6
    Q1 Q2 Q5 Q6
    R1 R2 R5 R6
    T1 T2 T3 T4 T5 T6
    V1 V2 V3 V4 V5 V6
)

for sample in ${samples[@]}; do
	sbatch ali_kallisto.slurm $sample
done

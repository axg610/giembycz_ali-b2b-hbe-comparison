#!/bin/bash

for sample in /work/giembycz_lab/HBE_Form_redownload/*; do
	sbatch hbe_kallisto.slurm  "$(basename $sample)"
done

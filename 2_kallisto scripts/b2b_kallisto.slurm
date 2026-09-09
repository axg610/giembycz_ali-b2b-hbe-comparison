#!/bin/bash

#SBATCH --job-name=b2b_kallisto
#SBATCH --output=/work/giembycz_lab/ALI-B2B-HBE_comparison/kallisto/logs/_%j_ali_kallisto.out
#SBATCH --error=/work/giembycz_lab/ALI-B2B-HBE_comparison/kallisto/logs/_%j_ali_kallisto.err

#SBATCH --time=03:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1

module load kallisto/0.46.1

samp=$1 #e.g., A1
index=/work/giembycz_lab/ag_tools/kallisto_index_cdna_GRCh38_p13

start=$(date "+%Y-%m-%d %H:%M:%S")

cd /work/giembycz_lab/tamkeenko_analysis/fastq/${samp}

kallisto quant -i ${index} -o /work/giembycz_lab/ALI-B2B-HBE_comparison/kallisto/b2b_kallisto/${samp}/ \
     ./*L001_R1*.fastq.gz \
     ./*L001_R2*.fastq.gz \
     ./*L002_R1*.fastq.gz \
     ./*L002_R2*.fastq.gz \
     --threads 4 \
     -b 100

end=$(date "+%Y-%m-%d %H:%M:%S")

echo "==========================================================="
echo "Job ID: $SLURM_JOBID"
echo "Job name: $SLURM_JOBNAME"
echo "Sample: $samp"
echo "Start: $start"
echo "End: $end"


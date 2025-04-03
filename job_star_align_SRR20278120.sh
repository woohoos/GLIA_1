#!/bin/bash
#SBATCH --job-name=align_SRR20278120
#SBATCH --output=logs/align_SRR20278120.log
#SBATCH --error=logs/align_SRR20278120.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

# 1. Activate conda and your STAR environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate star_env

# 2. Run STAR for alignment
STAR --runThreadN 4 \
  --genomeDir ~/rna-seq-pipeline/data/reference/star_index \
  --readFilesIn results/trimmed/SRR20278120_1.trimmed.fastq.gz results/trimmed/SRR20278120_2.trimmed.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix results/mapping/SRR20278120_ \
  --outSAMtype BAM SortedByCoordinate

#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --output=logs/star_index.log
#SBATCH --error=logs/star_index.err
#SBATCH --time=02:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4

source ~/miniconda3/etc/profile.d/conda.sh
conda activate star_env

STAR --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir ~/rna-seq-pipeline/data/reference/star_index \
  --genomeFastaFiles ~/rna-seq-pipeline/data/reference/reference.fasta \
  --sjdbGTFfile ~/rna-seq-pipeline/data/reference/reference.gtf \
  --sjdbOverhang 99

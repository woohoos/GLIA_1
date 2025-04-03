<<<<<<< HEAD
# rna-seq-pipeline
=======
TODO:

get directory structure done

download files

start working on snakemake
>>>>>>> 79daf331c7fec23fc500efdf0ca63e471594ec9c

<<<<<<< HEAD
Your local README content
=======
Remote README content
>>>>>>> origin/main
ø

#RNA-seq Pipeline (GLIA_1)

This repository contains a reproducible RNA-seq analysis pipeline developed using [Snakemake](https://snakemake.readthedocs.io), [Conda](https://docs.conda.io), and SLURM HPC scheduler. It includes raw data processing, trimming, alignment, and read counting.

---

##Project Structure

rna-seq-pipeline/ ├── config/ # Configuration files ├── workflow/ # Snakemake rules ├── data/ │ ├── raw/ # Raw FASTQ files │ └── reference/ # Genome + GTF ├── results/ │ ├── trimmed/ # Cutadapt outputs │ ├── fastqc_trimmed/ # FastQC for trimmed files │ ├── mapping/ # STAR aligned BAMs │ └── counts/ # Gene counts (featureCounts) ├── logs/ # SLURM logs ├── job_star_index.sh # SLURM script for STAR genome index ├── job_star_align_SRR20278120.sh # SLURM script for STAR alignment └── README.md # This file


---

##Conda Environments

```bash
conda create -n cutadapt_env cutadapt
conda create -n fastqc_env fastqc
conda create -n star_env star
conda install -n star_env -c bioconda subread  # for featureCounts

Trimming (Cutadapt)

cutadapt -j 4 \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -o results/trimmed/SRR20278120_1.trimmed.fastq.gz \
  -p results/trimmed/SRR20278120_2.trimmed.fastq.gz \
  data/raw/SRR20278120_1.fastq.gz data/raw/SRR20278120_2.fastq.gz \
  > results/trimmed/cutadapt_report.txt

 Quality Check (FastQC)

fastqc results/trimmed/*.fastq.gz -o results/fastqc_trimmed
 Reference Genome

cd data/reference

# Genome FASTA
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/GRCm39.primary_assembly.genome.fa.gz
gunzip GRCm39.primary_assembly.genome.fa.gz
mv GRCm39.primary_assembly.genome.fa reference.fasta

# GTF Annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz
gunzip gencode.vM32.annotation.gtf.gz
mv gencode.vM32.annotation.gtf reference.gtf

STAR Genome Index (via SLURM)

sbatch job_star_index.sh
Or manually:

STAR --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir data/reference/star_index \
  --genomeFastaFiles data/reference/reference.fasta \
  --sjdbGTFfile data/reference/reference.gtf \
  --sjdbOverhang 99

STAR Alignment (via SLURM)

sbatch job_star_align_SRR20278120.sh
Or manually:

STAR --runThreadN 4 \
  --genomeDir data/reference/star_index \
  --readFilesIn results/trimmed/SRR20278120_1.trimmed.fastq.gz results/trimmed/SRR20278120_2.trimmed.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix results/mapping/SRR20278120_ \
  --outSAMtype BAM SortedByCoordinate

Gene Counts (featureCounts)

featureCounts -T 4 -p \
  -a data/reference/reference.gtf \
  -o results/counts/counts_SRR20278120.tsv \
  results/mapping/SRR20278120_Aligned.sortedByCoord.out.bam


Use Snakemake to re-run all steps:

snakemake --cores 8 --use-conda --rerun-incomplete

# RNA-seq Pipeline (GLIA_1)

This repository contains a reproducible RNA-seq analysis pipeline developed using [Snakemake](https://snakemake.readthedocs.io), [Conda](https://docs.conda.io), and SLURM HPC resources.

---

## Project Structure

rna-seq-pipeline/ ├── config/ # Configuration files ├── workflow/ # Snakemake rules ├── data/ │ ├── raw/ # Raw FASTQ files │ └── reference/ # Genome FASTA, GTF, STAR index ├── results/ │ ├── trimmed/ # Trimmed FASTQs │ ├── mapping/ # STAR BAM files │ ├── counts/ # featureCounts outputs │ ├── fastqc/ # Pre-trimming FastQC reports │ └── fastqc_trimmed/ # Post-trimming FastQC reports ├── logs/ # SLURM job logs └── Snakefile # Main Snakemake workflow


---

## Conda Environments

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

fastqc data/raw/*.fastq.gz -o results/fastqc
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
Or run manually:

STAR --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir data/reference/star_index \
  --genomeFastaFiles data/reference/reference.fasta \
  --sjdbGTFfile data/reference/reference.gtf \
  --sjdbOverhang 99
STAR Alignment (via SLURM)

sbatch job_star_align_SRR20278120.sh
Or run manually:

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
⚠️ Note: This dataset appears to be ATAC-seq rather than RNA-seq. Therefore, featureCounts may assign very few reads if used with RNA-seq annotations.

Reproducibility

Use Snakemake to re-run everything in one go:

snakemake --cores 8 --use-conda --rerun-incomplete

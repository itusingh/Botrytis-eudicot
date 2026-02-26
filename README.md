
Project
"A multi-plant transcriptomic atlas reveals conserved and lineage-specific defense architectures in response to Botrytis cinerea manuscript", full analysis code

RNA seq preprocessing and dual alignment pipeline for host and Botrytis cinerea
This repository contains SLURM-based scripts to perform quality control, adapter trimming, and alignment of paired end RNA seq data.

Raw data source
All raw sequencing data are publicly available in NCBI Sequence Read Archive.
BioProject ID: PRJNA1217477
Users can download raw fastq files directly from NCBI SRA.

Software requirements:
Operating system
Linux environment with SLURM scheduler and R environment

Required tools
FastQC
MultiQC
Trimmomatic
HISAT2 version 2.2.1
samtools
wget
conda or environment modules for tool management

Recommended
Install tools through conda or load via cluster modules to ensure version control.

Input files:
file_list.txt
A plain text file containing one SRA Run (per species) accession per line.

adapters.fa
Adapter sequences used by Trimmomatic.

Step 1. Download raw data and perform QC.sh

Step 2. Using the same raw data perform trimming

The script trim.sh performs the following steps for each Run ID listed in file_list.txt.

Download paired fastq files from SRA
Decompress fastq.gz files
Run Trimmomatic PE mode
Generate trimmed paired and unpaired fastq files

To run: sbatch trim.sh

Output
RunID_R1_trimmed_paired.fastq
RunID_R2_trimmed_paired.fastq
RunID_R1_trimmed_unpaired.fastq
RunID_R2_trimmed_unpaired.fastq

Step 3. Quality control

The script QC_pairedfastq.sh performs FastQC on trimmed paired reads and summarizes results using MultiQC.

To run

sbatch QC_pairedfastq.sh

Output

qc/fastqc_out directory containing FastQC reports
MultiQC summary report

Step 3. Alignment to host and Botrytis cinerea

The script alignment.sh performs dual alignment of trimmed paired reads.

Requirements before running

HISAT2 indices must be built for
Host reference genome
Botrytis cinerea reference genome

To build indices

hisat2-build host_genome.fa HostIndexPrefix
hisat2-build botrytis_genome.fa BotrytisIndexPrefix

Edit alignment.sh to specify correct index prefixes before running.

To run

sbatch alignment.sh

Output per sample

RunID_Host.bam
RunID_Bcin.bam
Alignment summary logs


Step 4: Host_normalization.R
This script cleans up the readcount data, performs TMM normalization, and calculates cpm for each organism within each sample. Needs to be adjusted per dataset.
Also make sure to do normalization of host and Bcin separately


Step 3a: model_means.R - Run Generalized Linear Mixed Model with Negative Binomial on host reads
This script inputs counts, sample IDs, and batch information for one host species and outputs for each gene:

anova results and variances for each model term
emmeans and standard errors
DEG data for the infected term, including log2FCs and p values
Run your sbatch command like this (you don't need to hard-code input files into the script)

sbatch sbatch_modelmeans.sh /path/to/counts.csv /path/to/sampleIDs.csv /path/to/batch.csv /path/to/output_dir/

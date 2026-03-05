
Project
Multi- eudicot plants transcriptomic atlas of Botrytis cinerea infection.

Goal: Compare host and Botrytis cinerea responses across ten plant species. Identify conserved and lineage specific defense programs. Quantify host and pathogen expression during infection.

Overview: This repository contains a SLURM based pipeline for preprocessing co transcriptome RNA seq data and performing dual alignment to host plant genomes and the Botrytis cinerea genome. The workflow includes raw data retrieval, quality control, adapter trimming, alignment, normalization, differential expression modeling, GO enrichment analysis, and orthology inference.

Raw data source
All raw sequencing data are publicly available in the NCBI Sequence Read Archive.
BioProject ID: PRJNA1217477
Users can download raw fastq files directly from NCBI SRA.

Software requirements:
Operating system: Linux environment with SLURM scheduler and R environment

Required tools
FastQC
MultiQC
Trimmomatic
HISAT2 version 2.2.1
samtools
wget

Tool installation: Tools can be installed using conda environments or loaded through cluster environment modules to maintain consistent software versions.

Input files:
1. file_list.txt: A plain text file containing one SRA Run (per species) accession per line.
2. adapters.fa: Adapter sequences used by Trimmomatic.
3. raw fasta file from ncbi SRA
   
Pipeline workflow
Step 0. Setup
Prepare the working environment, install required software, and prepare input files, including file_list.txt and adapters.fa.

Step 1. Download raw data and perform raw quality control
The script raw_QC.sh performs initial quality control. Run sbatch raw_QC.sh

Step 2. Adapter trimming and preprocessing
The script trim.sh processes each Run ID listed in file_list.txt and performs adapter and quality trimming using Trimmomatic in paired-end mode, and generate paired and unpaired trimmed reads. It also contains wget command through which you can download files directly from SRA. Update the path as per your path
To run: sbatch trim.sh

Output
RunID_R1_trimmed_paired.fastq
RunID_R2_trimmed_paired.fastq
RunID_R1_trimmed_unpaired.fastq
RunID_R2_trimmed_unpaired.fastq

Step 3. Quality control of trimmed reads
The script QC_pairedfastq.sh evaluates the quality of trimmed paired reads and summarizes results.
Run: sbatch QC_pairedfastq.sh

Outputs
qc/fastqc_out directory containing FastQC reports
MultiQC summary report

Step 4. Dual alignment to host and Botrytis cinerea genomes
The script alignment.sh aligns trimmed paired reads to both the host genome and the Botrytis cinerea genome. Also convert sam file to bam file using samtool

Requirements before running
HISAT2 indices must be generated for each reference genome.

Build indices
hisat2-build host_genome.fa HostIndexPrefix
hisat2-build botrytis_genome.fa BotrytisIndexPrefix

Edit alignment.sh to specify correct index prefixes before running.

Run: sbatch alignment.sh

Output per sample
RunID_Host.bam
RunID_Bcin.bam
Alignment summary logs

Step 5. Normalization of read counts
The script Host_normalization.R processes the read count data and performs normalization. Needs to be adjusted per dataset.

Functions
Clean read count tables
Perform TMM normalization
Calculate counts per million values
Host and Botrytis cinerea counts should be normalized separately.
Input
Gene level count matrices generated from aligned BAM files.

Step 6. Differential expression modeling
The script model_means.R performs statistical modeling of host gene expression using a generalized linear mixed model with a negative binomial distribution.

Inputs
counts.csv containing gene-level counts
sampleIDs.csv containing sample metadata
batch.csv containing batch information

Outputs per gene
ANOVA results and variance components for each model term
Estimated marginal means and standard errors
Differential expression results for the infection term, including log2 fold change and p values
Run your sbatch command like this (you don't need to hard-code input files into the script)
sbatch sbatch_modelmeans.sh /path/to/counts.csv /path/to/sampleIDs.csv /path/to/batch.csv /path/to/output_dir/

Step 7. GO enrichment analysis
Goal: Test which biological processes are overrepresented among infection-responsive genes.
Method
topGO enrichment using the Biological Process ontology.
Fisher exact test with Benjamini Hochberg correction.

Inputs
All expressed genes that were tested in the differential expression model.

Per species DEG table
A table with gene_id, log2FC, and FDR.

Gene to GO mapping
A gene2go mapping file for the same gene IDs used in the DEG table.

Outputs
Per species enrichment tables for up- and down-regulated sets.
Cross-species summary tables and heatmaps.

Step 8: Orthology inference using OrthoFinder
Goal: Identify orthologous gene groups across the ten host species. These orthogroups are used for cross-species comparisons of infection-responsive genes.

Input: Protein sequences for each species. Each species should have a single FASTA file containing primary transcript protein sequences.
module load orthofinder
orthofinder -f /group/kliebengrp/itusingh/orthofinder/primary_transcripts/ -t 16

Outputs
Orthogroups.tsv
Orthogroups.GeneCount.tsv
Single copy ortholog lists and inferred species tree.

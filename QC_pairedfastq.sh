#!/bin/bash
#SBATCH -D /group/kliebengrp/itusingh/eudicot/Capsicum_annum/Fel1_paired_fastq
#SBATCH -o /group/kliebengrp/itusingh/eudicot/Fel1_paired_fastq/QC_stdout-%j.txt
#SBATCH -e /group/kliebengrp/itusingh/eudicot/Fel1_paired_fastq/QC_stderr-%j.txt
#SBATCH -J alignment
#SBATCH -t 95:00:00
#SBATCH --mem 80G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=itusingh@ucdavis.edu


### Versions
# pip 23.0.1 from /usr/local/lib/python3.9/site-packages/pip (python 3.9)
# Python 2.7.16
# multiqc, version 1.12

### Set up environment
# conda create -n MultiQC
#eval "$(conda shell.bash hook)"
#conda activate MultiQC
# conda install MultiQC
# conda install fastqc

### Run QC
# make qc directories
mkdir -p qc/fastqc_out

# run fastqc
fastqc --threads 8 -o qc/fastqc_out *_paired.fastq 

# run multiqc
multiqc qc/fastqc_out

# deactivate conda environment
conda deactivate

### Proceed to alignment.sh working in the same directory and provide output in same directory

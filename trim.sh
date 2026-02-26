#!/bin/bash
#SBATCH -D /group/kliebengrp/itusingh/eudicot/chard
#SBATCH -o /group/kliebengrp/itusingh/chard_stdout-%j.txt
#SBATCH -e /group/kliebengrp/itusingh/chard_stderr-%j.txt
#SBATCH -J alignment
#SBATCH -t 17:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=itusingh@ucdavis.edu

# Load Conda environment (if needed)
module load conda
module activate trimmomatic
conda install -c bioconda trimmomatic

### Generate a list of the samples you want to run and save as file_list.txt
# Assign your file_list.txt to files variable
readarray -t files < file_list.txt

### Get raw fastq files
# Loop to download fastq files using your file list
for file in "${files[@]}"
do
  # Download R1 and R2 for the sample. Change path as needed

  wget -nv "http://slimsdata.genomecenter.ucdavis.edu/Data/vbr8zf85y/Unaligned/Project_DKRS_EUD_1/${file}_R1.fastq.gz"
  wget -nv "http://slimsdata.genomecenter.ucdavis.edu/Data/vbr8zf85y/Unaligned/Project_DKRS_EUD_1/${file}_R2.fastq.gz"
done

# Unzip the files
gunzip *.fastq.gz


### Run Trimmomatic
for file in "${files[@]}"
do
  echo ' '
  echo 'Trimming' $file '...'
  echo ' '
  trimmomatic PE -threads 8 ${file}_R1.fastq ${file}_R2.fastq \
  ${file}_R1_trimmed_paired.fastq ${file}_R1_trimmed_unpaired.fastq \
  ${file}_R2_trimmed_paired.fastq ${file}_R2_trimmed_unpaired.fastq \
  ILLUMINACLIP:adapters.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
  
done


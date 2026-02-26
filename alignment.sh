#!/bin/bash
#SBATCH -D /group/kliebengrp/itusingh/lettuce/slurm
#SBATCH -o /group/kliebengrp/itusingh/lettuce/slurm_stdout-%j.txt
#SBATCH -e /group/kliebengrp/itusingh/lettuce/slurm_stderr-%j.txt
#SBATCH -J alignment
#SBATCH -t 90:00:00
#SBATCH --mem 50G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=itusingh@ucdavis.edu

###Versions
#hisat2 version 2.2.1

###Set up environment
#conda create -n hisat2
eval "$(conda shell.bash hook)"
conda activate hisat2
#conda install hisat2
#conda install samtools

####### First, generate a list of the samples you want to run and save as file_list.txt
#assign your file_list.txt to files variable
readarray -t files < file_list.txt

####### LOOP to align each sample 2 fastqs to host and Bcin genome, convert to .bam
####### LOOP to align each sample 2 fastqs to host and Bcin genome, convert to .bam
for file in "${files[@]}"
do
  # Perform Host alignment
  echo ' '
  echo 'Aligning' $file 'to Host genome...'
  echo ' '
  hisat2  \
     -p 8 \
      -t \
      --summary-file 'LsHost_log.txt' \
      -I 5 \
      -5 10 \
      -3 5 \
      -x Lettuce_Salinas_V11 \
      -1 ${file}_R1_trimmed_paired.fastq \
      -2 ${file}_R2_trimmed_paired.fastq \
      --score-min L,0,-0.85 \
      -S ${file}_Host.sam \
     
  # Convert Host SAM to BAM
  echo ' '
  echo 'Converting .sam to .bam...'
  echo ' '
  samtools view -b ${file}_Host.sam > ${file}_Host.bam
  

  # Perform Bcin alignment
  echo ' '
  echo 'Aligning' $file 'to Botrytis_cinerea genome...'
  echo ' '
  hisat2 \
    -p 8 \
    -t \
    --summary-file 'Bcin_log.txt' \
    -I 5 \
    -5 10 \
    -3 5 \
    -x Botrytis_cinerea_index \ 
    -1 ${file}_R1_trimmed_paired.fastq \
    -2 ${file}_R2_trimmed_paired.fastq \
    --score-min L,0,-0.85 \
    -S ${file}_Bcin.sam \

  # Convert Bcin SAM to BAM
  echo ' '
  echo 'Converting .sam to .bam...'
  echo ' '
  samtools view -b ${file}_Bcin.sam > ${file}_Bcin.bam
  
done


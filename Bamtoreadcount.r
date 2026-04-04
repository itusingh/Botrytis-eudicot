##Loop for all .bamfile
## Install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")
install.packages("tidyverse")
BiocManager::install("rtracklayer")
install.packages("purrr")

## Load packages
library(tidyverse)
library(rtracklayer)
library(Rsubread)
library(tools)
library(purrr)

## Prepare references
## Convert gff to gtf (HanXRQr2.0-SUNRISE)
gff <- import.gff("GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.gff")
## Create a column "gene_id" that contains the gene name for every entry
gff$gene_id <- ifelse(is.na(gff$ID), gff$Parent, gff$ID)
## Export as gtf
export(gff, "GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.gtf", format = "gtf")

## Get read counts for the BAM files
bam_files <- list.files(pattern = "*_Host.bam")
for (bam_file in bam_files) {
  count_matrix <- featureCounts(file = bam_file,
                                annot.ext = "genomic.gtf",
                                isGTFAnnotationFile = TRUE,
                                isPairedEnd = TRUE)
  ## Format the counts dataframe
  counts <- count_matrix$counts
  counts <- as.data.frame(counts)
  counts <- rownames_to_column(counts, var = "transcript")
  spl_id <- sub(".bam$", "", bam_file)
  colnames(counts) <- c("transcript", paste0(spl_id))
  new_object_name <- paste0("counts_", spl_id)
  assign(new_object_name, counts)
  ## Save the counts dataframe to a CSV file
  write.csv(x = counts, file = paste0(spl_id, "_Host_readcounts.csv"))
}
## Get a list of all readcount CSV files
readcount_files <- list.files(pattern = "_Host_readcounts.csv")
## Create a list to store the readcount dataframes
readcount_list <- list()
## Loop through all readcount files and read them into dataframes
for (file in readcount_files) {
    readcount_df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
    readcount_list[[file]] <- readcount_df
}
## Merge all dataframes in the list by the "transcript" column
 merged_df <- Reduce(function(x, y) full_join(x, y, by = "transcript"), readcount_list)
 ## Save the merged dataframe as a CSV file
 write.csv(merged_df, file = "Fel2_Host_merged_readcounts.csv", row.names = FALSE)
 

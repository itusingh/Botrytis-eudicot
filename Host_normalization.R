library(tidyverse)
library(edgeR)

# ========== Step 1: Define File Paths ==========
# ========== Step 1: Define File Paths ==========
rep_files <- list(
  ZE1 = "D:/Postdoc_DJK/library/Cucumber/Bos1/Host_readcounts.csv",
  ZE2 = "D:/Postdoc_DJK/library/Cucumber/Bos2/Host_readcounts.csv",
  ZE3 = "D:/Postdoc_DJK/library/Cucumber/Bos3/Host_readcounts.csv"
)


# ========== Step 2: Function for Cleaning, Filtering, Normalizing ==========
normalize_bcin_cpmfilter <- function(file, rep_name) {
  df <- read.csv(file)
  
  # Clean transcript names and remove noncoding RNAs
  df$transcript <- sub("^transcript:", "", df$transcript)
  df <- df %>% filter(!grepl("^E", transcript))
  df <- df %>%
    mutate(gene = sub("\\.[0-9]+$", "", transcript)) %>%
    select(gene, everything(), -transcript)
  
  # Collapse transcript isoforms into gene-level counts
  df <- df %>%
    group_by(gene) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop")
  
  # Create DGEList
  dge <- DGEList(counts = df[,-1], genes = df$gene)
 
  #For normal i have 1 cpm in atleast 2 samples
  #for expressed I have 1 cpm in atleast 40 samples; used for DEGs and emmemeans
  # === CPM-based filtering: CPM ??? 1 in at least 40 samples ===
  keep <- rowSums(cpm(dge) >= 1) >= 40
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # === TMM normalization ===
  dge <- calcNormFactors(dge)
  norm_cpm <- cpm(dge, normalized.lib.sizes = TRUE)
  
  # Format normalized data
  norm_cpm_df <- as.data.frame(norm_cpm)
  norm_cpm_df <- norm_cpm_df %>%
    mutate(gene = dge$genes$gene) %>%
    relocate(gene)
  
  # Prefix sample columns
  colnames(norm_cpm_df)[-1] <- paste0(rep_name, "_", colnames(norm_cpm_df)[-1])
  
  return(norm_cpm_df)
}

# ========== Step 3: Normalize All Reps ==========
norm_list <- mapply(normalize_bcin_cpmfilter, rep_files, names(rep_files), SIMPLIFY = FALSE)

# ========== Step 4: Merge by Gene ==========
norm_merged <- reduce(norm_list, full_join, by = "gene")

# Replace NAs with 0 (optional, useful for heatmaps)
norm_merged[is.na(norm_merged)] <- 0
head(norm_merged)
# ========== Step 5: Write Output ==========
write.csv(norm_merged, "Cucumber_Host_Expressed_TMM_Normalized_CPM_Merged.csv", row.names = FALSE)




# ========== Step 5: CPM Density Plot ==========
ggplot(plot_data, aes(x = CPM, color = Replicate)) +
  geom_density() +
  scale_x_log10() +
  geom_vline(xintercept = c(0.5, 1, 2), linetype = "dashed", color = c("red", "blue", "darkgreen")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "CPM Density Distribution (Log10 Scale)",
    x = "CPM (log10)", y = "Density"
  )


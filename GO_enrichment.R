#!/usr/bin/env Rscript
# ============================================================
# Cross species GO enrichment with topGO and heatmap figure code
# Project: Multi eudicot transcriptomic atlas of Botrytis cinerea infection
#
# What this script does
# 1) Reads GO annotations from NCBI RefSeq GAF files and builds gene to GO mappings
# 2) Runs topGO enrichment per species for a provided gene universe and DEG list
# 3) Writes per species enrichment tables
# 4) Builds cross species binary heatmaps showing which GO terms are enriched
#
# What you must prepare before running
# A) Per species gene universe file: all genes tested in that species
#    Format: one gene ID per line
# B) Per species DEG file: significant genes for that species and contrast
#    Format: one gene ID per line
# C) Per species GAF file from NCBI RefSeq: gene_ontology.gaf.gz
#    IDs in GAF are usually NCBI GeneIDs
# D) If your DEG IDs are not NCBI GeneIDs, map them to GeneIDs first
#
# Output
# - results_topgo/<species>/<species>_<direction>_<ontology>_topgo.csv
# - figures/<direction>_<ontology>_top_terms_binary.pdf
# - figures/<direction>_<ontology>_curated_binary.pdf
#
# Recommended repo layout
# scripts/01_topgo_cross_species.R  (this file)
# inputs/<species>/genes_all.txt
# inputs/<species>/genes_deg_up.txt
# inputs/<species>/genes_deg_down.txt
# annotations/<species>/gene_ontology.gaf.gz
# results_topgo/  (created)
# figures/        (created)
# ============================================================

suppressPackageStartupMessages({
  library(topGO)
  library(tidyverse)
  library(readr)
  library(stringr)
  library(pheatmap)
})

# -----------------------------
# User configuration
# -----------------------------

# Your species folder names. These must match the folder names under inputs/ and annotations/
species_order <- c(
  "cowpea", "commonbean", "cucumber", "zucchini",
  "chard", "spinach", "lettuce", "sunflower", "pepper", "tomato"
)

# Choose which DEG set to analyze
# Use "Up" or "Down"
which_side <- "Up"

# topGO ontology
# Use "BP", "MF", or "CC"
ontology <- "BP"

# Enrichment significance cutoff for cross species heatmaps
fdr_threshold <- 0.05

# Heatmap selection parameters
min_species <- 7
top_n <- 80

# A curated list of GO terms for a focused panel heatmap
# Use the exact format: "GO:0000000 - term name"
curated_terms <- c(
  "GO:0006952 - defense response",
  "GO:0009620 - response to fungus",
  "GO:0010200 - response to chitin",
  "GO:0009607 - response to biotic stimulus",
  "GO:0050832 - defense response to fungus",
  "GO:0009753 - response to jasmonic acid",
  "GO:0009751 - response to salicylic acid",
  "GO:0009737 - response to abscisic acid",
  "GO:0009733 - response to auxin",
  "GO:0009873 - ethylene-activated signaling pathway",
  "GO:0009867 - jasmonic acid mediated signaling pathway",
  "GO:0009863 - salicylic acid mediated signaling pathway",
  "GO:0009755 - hormone-mediated signaling pathway"
)

# Input folders
inputs_dir <- "inputs"
annotations_dir <- "annotations"

# Output folders
results_dir <- "results_topgo"
fig_dir <- "figures"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Helpers
# -----------------------------

# Read a one-column text file with gene IDs
read_gene_list <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  x <- read_lines(path)
  x <- x[!is.na(x)]
  x <- str_trim(x)
  x <- x[x != ""]
  unique(x)
}

# Build gene2go list from an NCBI RefSeq GAF file (gene_ontology.gaf.gz)
# Returns a named list: names are GeneIDs, values are character vectors of GO IDs
read_gaf_gene2go <- function(gaf_gz) {
  if (!file.exists(gaf_gz)) stop("Missing GAF: ", gaf_gz)
  
  gaf <- read_tsv(
    gaf_gz,
    comment = "!",
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  # GAF 2.x key columns used here
  # V2  DB Object ID (commonly NCBI GeneID for RefSeq)
  # V5  GO ID
  # V12 DB Object Type
  names(gaf)[c(2, 5, 12)] <- c("DB_Object_ID", "GO_ID", "ObjectType")
  
  gaf2 <- gaf %>%
    transmute(
      GeneID = as.character(DB_Object_ID),
      GO = as.character(GO_ID),
      ObjectType = as.character(ObjectType)
    ) %>%
    filter(!is.na(GeneID), !is.na(GO)) %>%
    filter(GeneID != "", GO != "") %>%
    # keep gene level or related objects, depending on how RefSeq populated the file
    filter(ObjectType %in% c("gene", "protein", "transcript")) %>%
    distinct(GeneID, GO)
  
  gene2go_tbl <- gaf2 %>%
    group_by(GeneID) %>%
    summarise(GO = list(unique(GO)), .groups = "drop")
  
  setNames(gene2go_tbl$GO, gene2go_tbl$GeneID)
}

# Run topGO for one species
# all_genes and sig_genes must be in the same ID space as gene2go names
run_topgo_one <- function(all_genes, sig_genes, gene2go, ontology = "BP") {
  all_genes <- unique(all_genes)
  sig_genes <- intersect(unique(sig_genes), all_genes)
  
  geneList <- factor(as.integer(all_genes %in% sig_genes))
  names(geneList) <- all_genes
  
  GOdata <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.gene2GO,
    gene2GO = gene2go
  )
  
  # Use a conservative algorithm that accounts for GO graph structure
  res <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  
  tab <- GenTable(
    GOdata,
    Fisher = res,
    orderBy = "Fisher",
    topNodes = length(score(res))
  )
  
  tab %>%
    as_tibble() %>%
    mutate(
      Pvalue = suppressWarnings(as.numeric(Fisher)),
      FDR = p.adjust(Pvalue, method = "BH")
    ) %>%
    select(GO.ID, Term, Annotated, Significant, Expected, Pvalue, FDR)
}

# Read a topGO results CSV written by this script
read_topgo_result <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  
  species <- sub("_.*$", "", basename(path))
  
  df %>%
    mutate(
      Species = species,
      GO_Term = paste0(`GO.ID`, " - ", Term),
      Enriched = as.integer(!is.na(FDR) & FDR < fdr_threshold)
    ) %>%
    select(Species, `GO.ID`, Term, FDR, GO_Term, Enriched)
}

# Build cross species binary heatmaps from per species topGO outputs
make_binary_heatmaps <- function(result_files,
                                 out_prefix,
                                 species_order,
                                 curated_terms,
                                 min_species,
                                 top_n) {
  dat <- map_dfr(result_files, read_topgo_result)
  
  wide <- dat %>%
    group_by(GO_Term, Species) %>%
    summarise(Enriched = as.integer(any(Enriched == 1)), .groups = "drop") %>%
    pivot_wider(names_from = Species, values_from = Enriched, values_fill = 0)
  
  species_cols <- intersect(species_order, setdiff(names(wide), "GO_Term"))
  if (length(species_cols) == 0) stop("No species columns found in wide table")
  
  wide <- wide %>%
    mutate(count_enriched = rowSums(across(all_of(species_cols)), na.rm = TRUE))
  
  # Top terms heatmap
  sel_top <- wide %>%
    filter(count_enriched >= min_species) %>%
    arrange(desc(count_enriched), GO_Term) %>%
    slice_head(n = top_n)
  
  mat_top <- sel_top %>%
    select(GO_Term, all_of(species_cols)) %>%
    column_to_rownames("GO_Term") %>%
    as.matrix()
  
  pdf(paste0(out_prefix, "_top_terms_binary.pdf"), width = 12, height = 9)
  pheatmap(
    mat_top,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = paste0("Top GO terms enriched in at least ", min_species, " species (FDR<", fdr_threshold, ")"),
    fontsize_row = 6.5,
    fontsize_col = 9,
    border_color = "grey40",
    legend_breaks = c(0, 1),
    legend_labels = c("Not enriched", "Enriched")
  )
  dev.off()
  
  # Curated heatmap in fixed order
  sel_cur <- wide %>%
    filter(GO_Term %in% curated_terms) %>%
    mutate(GO_Term = factor(GO_Term, levels = curated_terms)) %>%
    arrange(GO_Term)
  
  mat_cur <- sel_cur %>%
    select(GO_Term, all_of(species_cols)) %>%
    column_to_rownames("GO_Term") %>%
    as.matrix()
  
  pdf(paste0(out_prefix, "_curated_binary.pdf"), width = 12, height = 9)
  pheatmap(
    mat_cur,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = paste0("Curated GO terms (FDR<", fdr_threshold, ")"),
    fontsize_row = 6.5,
    fontsize_col = 9,
    border_color = "grey40",
    legend_breaks = c(0, 1),
    legend_labels = c("Not enriched", "Enriched")
  )
  dev.off()
  
  invisible(wide)
}

# -----------------------------
# Main workflow
# -----------------------------

message("Running topGO across species")
message("Direction: ", which_side, ", Ontology: ", ontology)

result_files <- c()

for (sp in species_order) {
  message("Species: ", sp)
  
  sp_inputs <- file.path(inputs_dir, sp)
  sp_ann <- file.path(annotations_dir, sp)
  
  all_file <- file.path(sp_inputs, "genes_all.txt")
  deg_file <- if (which_side == "Up") {
    file.path(sp_inputs, "genes_deg_up.txt")
  } else {
    file.path(sp_inputs, "genes_deg_down.txt")
  }
  gaf_file <- file.path(sp_ann, "gene_ontology.gaf.gz")
  
  all_genes <- read_gene_list(all_file)
  sig_genes <- read_gene_list(deg_file)
  gene2go <- read_gaf_gene2go(gaf_file)
  
  # Keep only genes that exist in gene2go to avoid empty annotation issues
  all_genes2 <- intersect(all_genes, names(gene2go))
  sig_genes2 <- intersect(sig_genes, all_genes2)
  
  if (length(all_genes2) < 100) {
    stop("Too few annotated genes after intersect for ", sp, ". Check ID mapping.")
  }
  
  tab <- run_topgo_one(all_genes2, sig_genes2, gene2go, ontology = ontology)
  
  sp_out_dir <- file.path(results_dir, sp)
  dir.create(sp_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_csv <- file.path(
    sp_out_dir,
    paste0(sp, "_", which_side, "_", ontology, "_topgo.csv")
  )
  
  write_csv(tab, out_csv)
  result_files <- c(result_files, out_csv)
  
  message("Wrote: ", out_csv)
}

message("Building cross species heatmaps")
out_prefix <- file.path(fig_dir, paste0(which_side, "_", ontology))

make_binary_heatmaps(
  result_files = result_files,
  out_prefix = out_prefix,
  species_order = species_order,
  curated_terms = curated_terms,
  min_species = min_species,
  top_n = top_n
)

message("Done")


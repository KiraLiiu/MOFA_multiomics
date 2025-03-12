#############################################################################
# RNA-seq Analysis Pipeline using DESeq2
# 
# This script performs differential expression analysis using DESeq2 with support
# for multiple group comparisons and visualization.
#
# Author: Kira Liu
# Date: 2025-03-10
#############################################################################

### Step 1: Load Required Packages ###########################################
library(tidyverse)   # For data manipulation and visualization
library(DESeq2)      # For differential expression analysis
library(data.table)  # For efficient data reading (fread)
library(org.Mm.eg.db) # Mouse annotation database 

# Optional packages
# library(org.Hs.eg.db)  # Human annotation database
#library(openxlsx)    # For Excel file handling (if needed)
# library(AnnotationDbi) # Annotation database interface

# Load custom functions for RNA-seq analysis
source("./scripts/functions/RNAseq_DESeq2_Pipline_Multi_Groups.R")

# Set output directory
output_dir <- "./data/processed/RNAseq_DESeq2_results"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

### Step 2: Load and Prepare Data ##########################################
# Read count matrix
count_matrix <- fread("./data/raw/rsem_genes.txt", header = FALSE)

# format the count matrix
count_matrix <- data.frame(count_matrix)
rownames(count_matrix) <- count_matrix[,1]
count_matrix <- count_matrix[,-1]
colnames(count_matrix) <- count_matrix[1,]
count_matrix <- count_matrix[-1,]

# Format the count matrix
# - Set first column as rownames
# - Remove the first column after setting as rownames
sample <- colnames(count_matrix)

### Step 3: Create Sample Information #####################################
# Create sample information data frame
sample_info <- read.csv("./data/metadata/sample_info_filtered.csv")

# Check sample counts before running analysis
unique(sample_info$Group)
table(sample_info$Group)

### Step 4: Define Analysis Parameters ###################################
# Define groups and comparisons
groups <- c("R1AA", "R1CC", "R1MD", "WD2AA", "WD2CC", "WD2MD")

# sample_info <- sample_info %>%
#   filter(Group %in% groups)

# Define comparisons
# Format: list of c("Experimental", "Control") pairs
comparisons <- list(
  c("R1CC", "R1AA"),
  c("R1MD", "R1AA"),
  c("WD2CC", "WD2AA"),
  c("WD2MD", "WD2AA"),
  c("WD2AA", "R1AA"),
  c("WD2CC", "R1CC"),
  c("WD2MD", "R1MD")
  )

### Step 5: Run Analysis ################################################
# Option 1: Full analysis with annotation and duplicate handling
results_full <- analyze_rnaseq_multi(
  count_matrix = count_matrix,
  sample_info = sample_info,
  groups = groups,
  comparisons = comparisons,
  organism = "org.Mm.eg.db",     # Mouse database
  keytype = "ENSEMBL",           # ID type
  # annotation_cols = c("SYMBOL", "ENTREZID", "GENENAME"), # Get additional annotations "SYMBOL" "ENSEMBL' "ENTREZID' "REFSEQ' "UNIGENE"
  remove_duplicates = TRUE,      # Remove duplicates
  remove_na_symbol = TRUE,       # Remove NA symbols
  duplicate_agg = mean,          # Aggregate duplicates using mean or NULL
  report_stats = TRUE,           # Show statistics
  save_tables = TRUE,            # Save results to files
  output_dir = paste0(output_dir, "/Results_with_annotation")        # Output directory
)

# Option 2: Simple analysis without annotation
results_simple <- analyze_rnaseq_simple(
  count_matrix = count_matrix,  # Count matrix
  sample_info = sample_info,    # Sample information
  groups = groups,              # List of groups
  comparisons = comparisons,    # List of comparisons
  min_count = 10,               # Minimum count for filtering
  save_tables = TRUE,           # Save results to files
  output_dir = paste0(output_dir, "/Results_without_annotation")       # Output directory
)

# Select top 5000 variable genes for visualization
vst_counts <- as.matrix(vst_counts)
gene_vars <- rowVars(vst_counts)
top_var_genes <- order(gene_vars, decreasing = TRUE)[1:min(5000, length(gene_vars))]
vst_counts_top <- vst_counts[top_var_genes,]
write.csv(vst_counts_top, file = "./data/processed/RNAseq_DESeq2_results/transcriptomics_normalized_vst_top5000.csv")
transcriptomics_normalized <- vst_counts_top

# # Remove objects from memory
# rm(list = setdiff(ls(), c("transcriptomics_normalized")))
# gc()

### Step 6: Access Results #############################################
# The analysis creates several output files in your output directory:
# 1. DESeq2 results for each comparison (CSV files)
# 2. Normalized count matrices (VST and size-factor normalized)
# 3. Quality control plots in the 'plots' subdirectory:
#    - MA plots
#    - PCA plots (global and pairwise)
#    - Sample distance heatmaps
#    - Correlation heatmaps

# You can also access results programmatically:

# ## From annotated analysis (results_annotated):
# # Access DESeq2 results with annotation
# results_annotated <- results_full$annotated_results
# # Access condition table
# condition_table <- results_full$condition_table
# # Access DESeq2 object for custom analysis
# deseq_obj <- results_full$deseq_results$dds
# # Access VST object for custom visualization
# vst_obj <- results_full$deseq_results$vst
# 
# ## From simple analysis (results_simple):
# # Access raw DESeq2 results
# raw_results <- results_simple$raw_results
# # Access normalized count data
# vst_counts <- results_simple$vst_counts
# size_norm_counts <- results_simple$size_normalized_counts
# # Access DESeq2 object
# deseq_obj_simple <- results_simple$deseq_object
# # Access VST object
# vst_obj_simple <- results_simple$vsd_object
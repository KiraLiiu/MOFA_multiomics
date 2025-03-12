# PTM Data Processing Execution
# ============================
# This script executes the PTM data processing pipeline

# Load necessary libraries
# -----------------------
library(openxlsx)       # For Excel file handling
library(tidyverse)      # For data manipulation
library(stringr)        # For string manipulation
library(compositions)   # For CLR transformation

# Source the functions file
source("./scripts/functions/PTM_proteomics_cleaning_pipeline.R")

# Define input and output parameters
file_path <- "./data/raw/yoyo-ptm-stoichdata.xlsx"
sheet_number <- 2
output_dir <- "./data/processed/PTM_proteomics_normalized"

# Make sure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# 1. Load data
ptm_matrix <- read.xlsx(xlsxFile = file_path, sheet = sheet_number,
                        colNames = FALSE, rowNames = TRUE)

# 2. Process sample names
colnames(ptm_matrix) <- str_extract(ptm_matrix[1,], "\\d+[A-Z]")
colnames(ptm_matrix)
ptm_matrix <- ptm_matrix[-1,]

# Function to add leading zeros for single-digit numbers
format_sample_id <- function(id) {
  gsub("^([0-9])([A-Z])", "0\\1\\2", id)  # If a number is a single digit, add a leading zero
}

colnames(ptm_matrix) <- sapply(colnames(ptm_matrix), format_sample_id)
# Check the updated sample IDs
colnames(ptm_matrix)

# 3. Convert to numeric
ptm_matrix <- convert_to_numeric(ptm_matrix)

# 4. Extract histone families
histone_families <- extract_histone_families(ptm_matrix)

# 5. Apply CLR transformation
ptm_clr_transformed <- apply_clr_transformation(ptm_matrix, histone_families)

# 6. Save the transformed data
write.csv(ptm_clr_transformed, paste0(output_dir, "/ptm_clr_transformed.csv"), row.names = TRUE)

message("PTM data processing complete. Results saved to ", output_dir)

# 7 Boxplots showing distribution by histone family

plots <- list()

message("Creating boxplots...")
boxplot <- create_simple_boxplot(
  ptm_matrix, 
  ptm_clr_transformed,
  n_families = 8
)
plots[["boxplot"]] <- boxplot

# 8. Save all visualizations to a single PDF file
message("Saving all plots to PDF...", output_dir)
save_plots_to_pdf(plots, paste0(output_dir, "/ptm_transformation_visualization.pdf"))

ptm_proteomics_normalized <- ptm_clr_transformed

# # 9. Clean up intermediate objects
# rm(list = setdiff(ls(), "ptm_clr_transformed"))  # Remove intermediate objects
# gc()  # Run garbage collection to free memory

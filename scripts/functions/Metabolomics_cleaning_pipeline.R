#' Metabolomics Data Processing Pipeline (Core Version)
#' A streamlined pipeline for processing metabolomics data without visualizations
#' 
#' Author: Kira Liu
#' Date: 2025-03-07

# =========================================================================
# 1. Load Required Packages
# =========================================================================

required_packages <- c("tidyverse", "limma")
for(pkg in required_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# =========================================================================
# 2. Data Preprocessing Functions
# =========================================================================

#' Create design table from sample information
#' @param sample_info Sample information dataframe
#' @param id_col Column name containing sample IDs
#' @param group_col Column name containing group information
#' @param reference_group Optional reference group level
#' @return Dataframe with sample IDs and design information
create_design_table <- function(sample_info, id_col, group_col, reference_group = NULL) {
  # Check if specified columns exist in the dataset
  if(!(id_col %in% colnames(sample_info))) 
    stop(paste("ID Column", id_col, "not found in sample_info"))
  if(!(group_col %in% colnames(sample_info))) 
    stop(paste("Group Column", group_col, "not found in sample_info"))
  
  # Create design table
  design_df <- sample_info %>%
    select(ID = all_of(id_col), Group = all_of(group_col))
  
  # Handle factor levels for proper contrasts
  if(!is.null(reference_group)) {
    if(!(reference_group %in% design_df$Group)) {
      stop(paste("Reference group", reference_group, "not found in data"))
    }
    # Make sure reference group is the first level
    design_df$Group <- relevel(factor(design_df$Group), ref = reference_group)
  } else {
    design_df$Group <- factor(design_df$Group)
  }
  
  return(design_df)
}

#' Filter out metabolites with too many missing values
#' @param data_matrix Data matrix
#' @param max_missing_pct Maximum allowed percentage of missing values (0-100)
#' @return Filtered data matrix
filter_missing_values <- function(data_matrix, max_missing_pct = 50) {
  # Convert zeros to NA
  data_matrix[data_matrix == 0] <- NA
  
  # Calculate missing percentage for each feature
  missing_counts <- rowSums(is.na(data_matrix))
  missing_pct <- (missing_counts / ncol(data_matrix)) * 100
  
  # Keep features with acceptable missing value percentage
  keep_features <- missing_pct <= max_missing_pct
  filtered_matrix <- data_matrix[keep_features, , drop = FALSE]
  
  # Add attributes for reporting
  attr(filtered_matrix, "n_removed") <- sum(!keep_features)
  attr(filtered_matrix, "missing_pct") <- missing_pct
  
  return(filtered_matrix)
}

#' Impute missing values with row means
#' @param data_matrix Data matrix with missing values
#' @return Imputed data matrix
impute_missing_values <- function(data_matrix) {
  # Convert to matrix if needed
  data_matrix <- as.matrix(data_matrix)
  
  # Simple row-mean imputation
  for(i in 1:nrow(data_matrix)) {
    missing <- is.na(data_matrix[i, ])
    if(any(missing)) {
      row_mean <- mean(data_matrix[i, !missing], na.rm = TRUE)
      data_matrix[i, missing] <- row_mean
    }
  }
  
  return(data_matrix)
}


#' Filter low abundance metabolites
#' @param data_matrix Data matrix
#' @param min_abundance Minimum mean abundance threshold
#' @return Filtered data matrix
filter_low_abundance <- function(data_matrix, min_abundance = 1) {
  # Calculate mean abundance for each metabolite
  mean_abundance <- rowMeans(data_matrix, na.rm = TRUE)
  
  # Keep only metabolites with mean above threshold
  keep_features <- mean_abundance > min_abundance
  filtered_matrix <- data_matrix[keep_features, , drop = FALSE]
  
  # Add attributes for reporting
  attr(filtered_matrix, "n_low_abundance_removed") <- sum(!keep_features)
  
  return(filtered_matrix)
}

#' Normalize data using improved approach
#' @param data_matrix Data matrix
#' @param method Normalization method
#' @param output_dir Output directory for saving plots
#' @return Normalized data matrix



normalize_data <- function(data_matrix, 
                           method = "quantile", 
                           output_dir = output_dir
                           ) {
  require(limma)
  
  # Total Ion Count (TIC)
  tic_factors <- colSums(data_matrix)
  # Apply TIC normalization
  tic_normalized <- sweep(data_matrix, 2, tic_factors, "/") * 1e6
  # use limma "normalizeBetweenArrays" and method to remove batch effects
  normalized_matrix <- normalizeBetweenArrays(tic_normalized, method = method)
  
  # Save the 2 boxplots into 1 figure as a PDF file
  pdf(file = paste0(output_dir, "metabolomics_normalization_boxplot_comparison.pdf"), width=15, height=5)
  par(mfrow=c(1,3))  
  boxplot(data_matrix, outline=FALSE, notch=FALSE, las=2, main="Original Data")
  boxplot(tic_normalized, outline=FALSE, notch=FALSE, las=2, main="TIC Normalization")
  boxplot(normalized_matrix, outline=FALSE, notch=FALSE, las=2, main="Batch Normalization")
  dev.off()  # Close PDF file
  
  return(normalized_matrix)
}


# =========================================================================
# 3. Main Pipeline

# Define the main processing pipeline
process_metabolomics_data <- function(met_matrix, sample_info, output_dir, 
                                      max_missing_pct = 50, min_abundance = 10,
                                      normalization_method = "quantile") {
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Step 1: Creating Design Table\n")
  design_df <- create_design_table(sample_info = sample_info,
                                   id_col = "ID",
                                   group_col = "Group",
                                   reference_group = NULL)
  
  cat("Step 2: Reordering Columns in Data Matrix\n")
  met_matrix <- met_matrix[, design_df$ID, drop = FALSE]  # Ensure correct column order
  
  cat("Step 3: Filtering Missing Values (Threshold:", max_missing_pct, "%)\n")
  filtered_matrix <- filter_missing_values(data_matrix = met_matrix,
                                           max_missing_pct = max_missing_pct)
  
  cat("Step 4: Imputing Missing Values\n")
  filtered_matrix <- impute_missing_values(data_matrix = filtered_matrix)
  
  cat("Step 5: Filtering Low Abundance Features (Min:", min_abundance, ")\n")
  filtered_matrix <- filter_low_abundance(data_matrix = filtered_matrix,
                                          min_abundance = min_abundance)
  
  cat("Step 6: Normalizing Data (Method:", normalization_method, ")\n")
  normalized_matrix <- normalize_data(data_matrix = filtered_matrix,
                                      method = normalization_method,
                                      output_dir = output_dir)
  
  # Save intermediate and final outputs
  write.csv(filtered_matrix, file = file.path(output_dir, "filtered_matrix.csv"))
  write.csv(normalized_matrix, file = file.path(output_dir, "normalized_matrix.csv"))
  
  cat(paste0("Process completed! Results saved in: ", output_dir, "\n"))
  
  return(normalized_matrix)
}


#' Metabolomics Data Processing Pipeline (Enhanced Version)
#' A streamlined pipeline for processing metabolomics data with flexible normalization
#' 
#' Author: Kira Liu
#' Date: 2025-03-12 (Updated)

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
    dplyr::select(ID = all_of(id_col), Group = all_of(group_col))
  
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

# =========================================================================
# 3. Enhanced Normalization Functions
# =========================================================================

#' Apply Total Ion Count (TIC) normalization
#' @param data_matrix Data matrix
#' @param scale_factor Scale factor to apply after normalization (default: 1e6)
#' @return TIC-normalized data matrix
apply_tic_normalization <- function(data_matrix, scale_factor = 1e6) {
  # Calculate TIC for each sample
  tic_factors <- colSums(data_matrix, na.rm = TRUE)
  
  # Apply TIC normalization: divide each sample by its TIC and multiply by scale factor
  tic_normalized <- sweep(data_matrix, 2, tic_factors, "/") * scale_factor
  
  return(tic_normalized)
}

#' Apply generalized logarithm (glog) transformation without MSnbase dependency
#' @param data_matrix Data matrix
#' @param lambda Transformation parameter
#' @return glog-transformed data matrix
apply_glog_transformation <- function(data_matrix, lambda = 1e-10) {
  # Direct implementation of the generalized logarithm
  # glog(x) = log((x + sqrt(x^2 + lambda))/2)
  
  # Ensure data_matrix is a matrix
  data_matrix <- as.matrix(data_matrix)
  
  # Apply the glog transformation to each element
  glog_transform <- function(x, lambda) {
    log2((x + sqrt(x^2 + lambda))/2)
  }
  
  # Create a new matrix with the same dimensions
  result <- matrix(NA, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  rownames(result) <- rownames(data_matrix)
  colnames(result) <- colnames(data_matrix)
  
  # Apply transformation to each element
  for (i in 1:nrow(data_matrix)) {
    for (j in 1:ncol(data_matrix)) {
      result[i, j] <- glog_transform(data_matrix[i, j], lambda)
    }
  }
  
  return(result)
}

# Update the normalize_data_flexible function to use the standalone glog
normalize_data_flexible <- function(data_matrix, 
                                    method = c("tic_limma", "tic", "glog", "limma", "log2"),
                                    limma_method = "quantile",
                                    lambda = 1e-10,
                                    log_transform = FALSE,
                                    output_dir = NULL) {
  
  method <- match.arg(method)
  require(limma)
  
  # Original data
  original_data <- data_matrix
  
  # Store each step for plotting
  result_steps <- list(
    Original = original_data
  )
  
  # Apply log2 transformation if requested
  if(log_transform && method %in% c("tic", "tic_limma", "limma")) {
    data_matrix <- log2(data_matrix + 1)
    result_steps$Log2 <- data_matrix
  }
  
  # Apply selected normalization method
  if(method == "tic_limma") {
    # TIC normalization followed by limma
    tic_normalized <- apply_tic_normalization(data_matrix)
    result_steps$TIC <- tic_normalized
    
    normalized <- normalizeBetweenArrays(tic_normalized, method = limma_method)
    result_steps$Final <- normalized
    
  } else if(method == "tic") {
    # TIC normalization only
    normalized <- apply_tic_normalization(data_matrix)
    result_steps$Final <- normalized
    
  } else if(method == "glog") {
    # Generalized logarithm transformation using standalone function
    message("Applying standalone glog transformation with lambda = ", lambda)
    normalized <- apply_glog_transformation(data_matrix, lambda = lambda)
    result_steps$Final <- normalized
    
  } else if(method == "limma") {
    # limma normalization only
    normalized <- normalizeBetweenArrays(data_matrix, method = limma_method)
    result_steps$Final <- normalized
    
  } else if(method == "log2") {
    # Just log2 transformation
    normalized <- log2(data_matrix + 1)
    result_steps$Final <- normalized
  }
  
  # Create boxplots if output directory is provided
  if(!is.null(output_dir)) {
    create_normalization_boxplots(result_steps, method, output_dir)
  }
  
  return(normalized)
}


#' Create boxplots showing normalization steps
#' @param result_steps List of data matrices from each normalization step
#' @param method Normalization method used
#' @param output_dir Output directory
create_normalization_boxplots <- function(result_steps, method, output_dir) {
  # Create PDF file
  pdf_file <- file.path(output_dir, paste0("normalization_boxplots_", method, ".pdf"))
  pdf(file = pdf_file, width = 12, height = 7)
  
  # Set layout based on number of steps
  n_steps <- length(result_steps)
  if(n_steps <= 3) {
    par(mfrow = c(1, n_steps))
  } else {
    par(mfrow = c(2, ceiling(n_steps/2)))
  }
  
  # Create boxplot for each step
  for(i in seq_along(result_steps)) {
    step_name <- names(result_steps)[i]
    boxplot(result_steps[[i]], 
            outline = FALSE, 
            notch = FALSE, 
            las = 2, 
            main = paste0(step_name, " Data"),
            cex.axis = 0.7)
  }
  
  # Close PDF file
  dev.off()
  
  # Print message
  message("Normalization boxplots saved to: ", pdf_file)
}

# =========================================================================
# 4. Main Pipeline
# =========================================================================

#' Define the main processing pipeline
#' @param met_matrix Metabolomics data matrix
#' @param sample_info Sample information
#' @param output_dir Output directory
#' @param max_missing_pct Maximum percentage of missing values allowed
#' @param min_abundance Minimum abundance threshold
#' @param normalization_method Method for normalization: "tic_limma", "tic", "glog", "limma", or "log2"
#' @param limma_method Method for limma normalization
#' @param lambda Parameter for glog transformation
#' @param log_transform Apply log2 transformation before normalization
#' @return Normalized data matrix
process_metabolomics_data <- function(met_matrix, 
                                      sample_info, 
                                      output_dir, 
                                      max_missing_pct = 50, 
                                      min_abundance = 10,
                                      normalization_method = c("tic_limma", "tic", "glog", "limma", "log2"),
                                      limma_method = "quantile",
                                      lambda = 1e-10,
                                      log_transform = FALSE) {
  
  normalization_method <- match.arg(normalization_method)
  
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
  normalized_matrix <- normalize_data_flexible(
    data_matrix = filtered_matrix,
    method = normalization_method,
    limma_method = limma_method,
    lambda = lambda,
    log_transform = log_transform,
    output_dir = output_dir
  )
  
  # Save intermediate and final outputs
  write.csv(filtered_matrix, file = file.path(output_dir, "filtered_matrix.csv"))
  write.csv(normalized_matrix, file = file.path(output_dir, "normalized_matrix.csv"))
  
  cat(paste0("Process completed! Results saved in: ", output_dir, "\n"))
  
  return(normalized_matrix)
}

# =========================================================================
# 5. Example Usage
# =========================================================================

# Example usage:
# 
# # 1. TIC followed by limma normalization (default)
# normalized_data_tic_limma <- process_metabolomics_data(
#   met_matrix = met_matrix,
#   sample_info = sample_info,
#   output_dir = "results_tic_limma",
#   normalization_method = "tic_limma"
# )
#
# # 2. Glog transformation
# normalized_data_glog <- process_metabolomics_data(
#   met_matrix = met_matrix,
#   sample_info = sample_info,
#   output_dir = "results_glog",
#   normalization_method = "glog",
#   lambda = 1e-10
# )
#
# # 3. Simple TIC normalization
# normalized_data_tic <- process_metabolomics_data(
#   met_matrix = met_matrix,
#   sample_info = sample_info,
#   output_dir = "results_tic",
#   normalization_method = "tic"
# )
#
# # 4. Limma normalization only
# normalized_data_limma <- process_metabolomics_data(
#   met_matrix = met_matrix,
#   sample_info = sample_info,
#   output_dir = "results_limma",
#   normalization_method = "limma",
#   limma_method = "quantile",
#   log_transform = TRUE
# )

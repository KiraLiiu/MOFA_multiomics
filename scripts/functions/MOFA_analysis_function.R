
# Function to reorder and fill missing columns with NA
reorder_and_fill_na <- function(df, desired_order) {
  # Find missing columns
  missing_cols <- setdiff(desired_order, colnames(df))
  
  # Create a dataframe of NAs for missing columns
  na_df <- matrix(NA, nrow = nrow(df), ncol = length(missing_cols)) %>%
    as.data.frame()
  colnames(na_df) <- missing_cols
  
  # Merge original dataframe with missing columns
  df <- cbind(df, na_df)
  
  # Reorder columns to match desired order
  df <- df[, desired_order, drop = FALSE]
  
  return(df)
}


# Function for selecting high variance features
select_high_variance_features <- function(data_matrix, variance_threshold = 0.8) {
  # Compute variance for each feature (rows = features, columns = samples)
  feature_variance <- apply(data_matrix, 1, var, na.rm = TRUE)
  
  # Sort variances in descending order
  sorted_var <- sort(feature_variance, decreasing = TRUE)
  
  # Compute cumulative variance
  cumulative_var <- cumsum(sorted_var) / sum(sorted_var)
  
  # Find number of features explaining the given variance threshold
  num_features <- which(cumulative_var >= variance_threshold)[1]
  
  # Select the top features
  top_features <- names(sorted_var)[1:num_features]
  
  # Subset the original data matrix
  selected_features <- data_matrix[top_features, , drop = FALSE]
  
  # Return the filtered matrix
  return(as.data.frame(selected_features))
}



# function to extract top features for each factor
extract_top_features <- function(data, factor_column, value_column, factors = 1:9, output_dir = "./results/Top_Features", top_n = 1000) {
  
  # Create results directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Loop through each factor
  for (factor_id in factors) {
    factor_name <- paste0("Factor", factor_id)
    
    # Subset data for the current factor
    factor_weights <- subset(data, data[[factor_column]] == factor_name)
    
    if (nrow(factor_weights) == 0) {
      message(paste("No data found for", factor_name))
      next  # Skip if no data for this factor
    }
    
    # Save all features for the factor
    all_features_path <- file.path(output_dir, paste0(factor_name, "_all_features.xlsx"))
    write.xlsx(factor_weights, all_features_path)
    
    # Get top features based on absolute weight values
    factor_top_features <- factor_weights[order(-abs(factor_weights[[value_column]])), ][1:min(top_n, nrow(factor_weights)), ]
    
    # Save top features
    top_features_path <- file.path(output_dir, paste0(factor_name, "_top_features.xlsx"))
    write.xlsx(factor_top_features, top_features_path)
    
    # Store in list
    results_list[[factor_name]] <- list(
      all_features = factor_weights,
      top_features = factor_top_features
    )
    
    # Print progress
    message(paste("Saved results for", factor_name, "- All:", all_features_path, "Top:", top_features_path))
  }
  
  message("Feature extraction completed for all selected factors.")
  
  # Return structured list of all extracted features
  return(results_list)
}

# Example Usage:
# - `data`: Your MOFA weight matrix (e.g., `weights_df`)
# - `factor_column`: The column containing factor names (e.g., `"factor"`)
# - `value_column`: The column containing feature weight values (e.g., `"value"`)
# - `factors`: The factors to analyze (default: `1:9`)
# - `output_dir`: Where to save results
# - `top_n`: Number of top features to extract per factor (default: `1000`)

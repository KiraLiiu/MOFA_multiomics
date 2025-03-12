# Simplified PTM Data Processing Functions
# =======================================

#' Convert all values to numeric
#'
#' @param ptm_data Data frame with PTM data
#' @return Data frame with numeric values
convert_to_numeric <- function(ptm_data) {
  # Convert all values to numeric
  ptm_data[] <- lapply(ptm_data, function(x) as.numeric(as.character(x)))
  return(ptm_data)
}


#' Extract histone family information from row names
#'
#' @param ptm_data Data frame with row names containing histone info
#' @return Vector of histone family identifiers
extract_histone_families <- function(ptm_data) {
  # Extract histone family (assumes format like "H3_K4_HIST_me3")
  histone_families <- str_extract(rownames(ptm_data), "^[^_]+_[^_]+_[^_]+")
  return(histone_families)
}


#' Apply CLR transformation within each histone family
#'
#' @param ptm_data Data frame with PTM data
#' @param histone_families Vector of histone family identifiers
#' @return Data frame with CLR-transformed values
apply_clr_transformation <- function(ptm_data, histone_families) {
  # Remove samples with all missing values
  samples_with_missing <- colSums(is.na(ptm_data)) == nrow(ptm_data)
  missing_samples <- ptm_data[, samples_with_missing, drop = FALSE]
  ptm_data_filtered <- ptm_data[, !samples_with_missing, drop = FALSE]
  
  # Initialize matrix for transformed values
  clr_transformed <- matrix(NA, nrow = nrow(ptm_data_filtered), ncol = ncol(ptm_data_filtered))
  rownames(clr_transformed) <- rownames(ptm_data_filtered)
  colnames(clr_transformed) <- colnames(ptm_data_filtered)
  
  # Apply CLR transformation to each histone family separately
  unique_families <- unique(histone_families)
  
  for (family in unique_families) {
    family_indices <- which(histone_families == family)
    
    if (length(family_indices) > 1) {  # Apply CLR only if there are multiple PTMs in the family
      family_data <- ptm_data_filtered[family_indices, , drop = FALSE]  # Subset data for this family
      
      # Apply CLR transformation by column (sample)
      clr_transformed[family_indices, ] <- apply(family_data, 2, clr)
    } else {
      # For single PTMs, just use original values
      clr_transformed[family_indices, ] <- as.matrix(ptm_data_filtered[family_indices, , drop = FALSE])
    }
  }
  
  # Convert to dataframe
  clr_transformed <- as.data.frame(clr_transformed)
  
  # Add back missing samples
  if (ncol(missing_samples) > 0) {
    clr_transformed <- cbind(clr_transformed, missing_samples)
  }
  
  return(clr_transformed)
}


# Direct approach to create boxplot
create_simple_boxplot <- function(before_data, after_data, n_families = 8) {
  # Extract histone family directly from rownames
  histone_families <- str_extract(rownames(before_data), "^[^_]+_[^_]+_[^_]+")
  
  # Prepare long format data
  before_long <- before_data %>%
    as.data.frame() %>%
    rownames_to_column("PTM_ID") %>%
    mutate(histone_family = histone_families) %>%
    pivot_longer(cols = -c(PTM_ID, histone_family), 
                 names_to = "Sample", values_to = "Value") %>%
    mutate(Transformation = "Before CLR")
  
  after_long <- after_data %>%
    as.data.frame() %>%
    rownames_to_column("PTM_ID") %>%
    mutate(histone_family = histone_families) %>%
    pivot_longer(cols = -c(PTM_ID, histone_family), 
                 names_to = "Sample", values_to = "Value") %>%
    mutate(Transformation = "After CLR")
  
  # Combine and filter
  combined_data <- bind_rows(before_long, after_long) %>%
    filter(!is.na(Value), !is.na(histone_family))
  
  # Get top families
  family_counts <- table(combined_data$histone_family)
  top_families <- names(sort(family_counts, decreasing = TRUE)[1:min(n_families, length(family_counts))])
  
  # Plot
  ggplot(combined_data %>% filter(histone_family %in% top_families), 
         aes(x = histone_family, y = Value, fill = Transformation)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Before vs After CLR Transformation by Histone Family")
}



#' Save visualization plots to a PDF file
#'
#' @param plots List of ggplot objects
#' @param output_file Path to output PDF file
#' @return NULL
save_plots_to_pdf <- function(plots, output_file) {
  pdf(output_file, width = 12, height = 10)
  
  for (plot in plots) {
    if (inherits(plot, "ggplot")) {
      print(plot)
    }
  }
  
  dev.off()
  
  message(paste("Plots saved to", output_file))
}

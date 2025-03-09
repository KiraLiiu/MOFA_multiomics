# load custom functions
source("./scripts/functions/Metabolomics_cleaning_pipeline.R")

# Read in the data
met_matrix <- read.csv("./data/raw/final_mets.csv")

# pre-processing
rownames(met_matrix) <- met_matrix$Compound
met_matrix <- met_matrix[,-c(1:2)]
met_matrix <- t(met_matrix)

# Standardize sample IDs in data_matrix (Add leading zeros where needed) ？？？？
# Function to add leading zeros for single-digit numbers
format_sample_id <- function(id) {
  gsub("^([0-9])([A-Z])", "0\\1\\2", id)  # If a number is a single digit, add a leading zero
}

# Apply function to format sample IDs in design_table and data_matrix
colnames(met_matrix) <- sapply(colnames(met_matrix), format_sample_id)
# Check the updated sample IDs
colnames(met_matrix)

# Read sample information
sample_info <- read.csv("./data/metadata/sample_info_filtered.csv")
# Check unique group levels
unique(sample_info$Group)


# Filter sample information to include only relevant groups
sample_info <- sample_info %>%
  filter(Group %in% c("R1AA", "R1CC", "R1MD", "WD2AA", "WD2CC", "WD2MD"))


# Define the output directory
output_dir <- "data/processed/metabolomics_normalized/"

# Run the metabolomics data processing pipeline
normalized_metabolomics <- process_metabolomics_data(met_matrix = met_matrix, 
                                                     sample_info = sample_info,
                                                     output_dir = output_dir,
                                                     max_missing_pct = 50,
                                                     min_abundance = 10,
                                                     normalization_method = "quantile")

# remove intermediate objects and run garbage collection 
rm(list = setdiff(ls(), "normalized_metabolomics"))
gc()  # Run garbage collection to free memory



# # Output directory
# output_dir <- "data/processed/"
# 
# # Create a design table
# design_df <- create_design_table(sample_info = sample_info,
#                     id_col = "ID",
#                     group_col = "Group",
#                     reference_group = NULL)
# 
# # Reorder columns in data_matrix to match the order of IDs in design_table
# met_matrix <- met_matrix[,design_df$ID]
# 
# # Filter out metabolites with too many missing values
# filtered_matrix<- filter_missing_values(data_matrix = met_matrix,
#                                         max_missing_pct = 50 # percentage of maximum missing value
#                                         )
# 
# # Impute missing values
# filtered_matrix <-
#   impute_missing_values(data_matrix = filtered_matrix)
# 
# # Filter out low abundance features
# filtered_matrix <- filter_low_abundance(data_matrix = filtered_matrix,
#                                         min_abundance = 10)
# 
# # Normalize data
# normalized_matrix <- normalize_data(data_matrix = filtered_matrix,
#                                     method = "quantile",
#                                     output_dir = output_dir)
# 
# # Save the normalized matrix
# write.csv(normalized_matrix, file = paste0(output_dir, "normalized_matrix.csv"))
# 

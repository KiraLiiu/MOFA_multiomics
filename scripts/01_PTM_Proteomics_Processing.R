### load packages
library(openxlsx)
library(tidyverse)
library(DESeq2)
library(stringr)
library(data.table)

# ptm_matrix <- read.xlsx("./data/PTM/20250131_yoyo-liver-ptm_replicate-log2-vs-timematched-control.xlsx", colNames = FALSE, rowNames = TRUE)

ptm_matrix <- read.xlsx("./data/PTM/yoyo-ptm-stoichdata.xlsx", sheet = 2, colNames = FALSE, rowNames = TRUE)

name_extracted <- str_extract(ptm_matrix[1,], "\\d+[A-Z]")
length(name_extracted)
head(name_extracted)

name_formatted <- ifelse(!is.na(name_extracted),
                    sprintf("%02d%s", as.numeric(gsub("[A-Z]", "", name_extracted)), gsub("[0-9]", "", name_extracted)),
                    NA)

length(name_formatted)
head(name_formatted)

name_formatted

colnames(ptm_matrix) <- name_formatted
head(ptm_matrix)

# desired_order <- c("02N", "02R","05N", "05R", "11N", "11R",  # R1CC
#                    "07N", "07R", "09N", "09R", "12N", "12R", "13N", "13R",   # R1MD
#                    "21N", "21R", "34N", "34R", "43N", "43R",  # WD2CC
#                    "23N", "23R", "26N", "26R", "27N", "27R",  # WD2MD
#                    "01N", "01R", "04R", "06N", "06R", "08N", "08R",  # R1AA
#                    "10N", "14N", "15N", "15R", "16N", "16R"  # WD2AA
#                    )


desired_order <- c("01N", "01R", "04R", "06N", "06R", "08N", "08R",  # R1AA
                   "02N", "02R","05N", "05R", "11N", "11R",  # R1CC
                   "07N", "07R", "09N", "09R", "12N", "12R", "13N", "13R",   # R1MD
                   "10N", "14N", "15N", "15R", "16N", "16R",  # WD2AA
                   "21N", "21R", "34N", "34R", "43N", "43R",  # WD2CC
                   "23N", "23R", "26N", "26R", "27N", "27R"  # WD2MD
)



# # 5️⃣ 只保留 `desired_order` 里存在的列
# valid_order <- intersect(desired_order, name_formatted)
# valid_order
# 

# 5️⃣ 让 `tpm_matrix` 包含 `desired_order` 的所有列
missing_cols <- setdiff(desired_order, colnames(ptm_matrix))  # 找到缺失的列


# 6️⃣ 对缺失列填充 NA，并合并回 `tpm_matrix`
ptm_matrix[, missing_cols] <- NA


ptm_matrix <- ptm_matrix[-1,]

ptm_matrix <- ptm_matrix[, desired_order]

ptm_matrix[] <- lapply(ptm_matrix, function(x) as.numeric(as.character(x)))

head(ptm_matrix)

ptm_matrix[6,1] + ptm_matrix[7,1]


# Load necessary library
library(compositions)  # For CLR transformation


# Step 1: Identify and separate samples with all missing values
samples_with_missing <- colSums(is.na(ptm_matrix)) == nrow(ptm_matrix)  # TRUE if all values in a sample are NA
missing_samples <- ptm_matrix[, samples_with_missing, drop = FALSE]  # Store missing samples
ptm_matrix <- ptm_matrix[, !samples_with_missing]  # Keep only complete-case samples

# Extract histone family names
histone_families <- str_extract(rownames(ptm_matrix), "^[^_]+_[^_]+_[^_]+")

# Initialize a transformed data matrix
clr_transformed <- matrix(NA, nrow = nrow(ptm_matrix), ncol = ncol(ptm_matrix))
rownames(clr_transformed) <- rownames(ptm_matrix)
colnames(clr_transformed) <- colnames(ptm_matrix)

# Step 2: Apply CLR transformation within each histone family
unique_families <- unique(histone_families)

for (family in unique_families) {
  family_indices <- which(histone_families == family)
  
  if (length(family_indices) > 1) {  # Apply CLR only if there are multiple PTMs in the family
    family_data <- ptm_matrix[family_indices, , drop = FALSE]  # Subset data for this histone family
    
    # Apply CLR transformation
    clr_transformed[family_indices, ] <- apply(family_data, 2, clr)
  }
}

# Convert to dataframe
clr_transformed <- as.data.frame(clr_transformed)

# Step 3: Add back the missing samples (filled with NA)
clr_transformed <- cbind(clr_transformed, missing_samples)

# Ensure the original sample order is maintained
clr_transformed <- clr_transformed[, desired_order]

# Print transformed data
print(clr_transformed)

write.csv(clr_transformed, "./data/PTM/ptm_matrix_CLR_normalized.csv")

rm(list = ls())
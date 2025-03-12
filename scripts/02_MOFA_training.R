# ====
# This scripts is for MOFA Analysis
# ====

# Load the required libraries
library(MOFA2) # Load the MOFA2 package
library(MOFAdata) # Load the MOFAdata package
library(data.table) # Load the data.table package for efficient data manipulation
library(ggpubr) # involve ggline function
library(ggplot2) # involve ggplot function
library(tidyverse) # Load the tidyverse package for data manipulation

# # Clean the environment without removing the omics datasets
# rm(list = setdiff(ls(), c("metabolomics_normalized", "ptm_proteomics_normalized", "transcriptomics_normalized")))
# gc()  # Run garbage collection to free memory

# Load the normalized omics datasets
metabolomics_normalized <- read.csv("data/processed/metabolomics_normalized/metabolomics_normalized.csv", check.names = FALSE, row.names = 1)
ptm_proteomics_normalized <- read.csv("data/processed/ptm_proteomics_normalized/ptm_clr_transformed.csv", check.names = FALSE, row.names = 1)
transcriptomics_normalized <- read.csv("data/processed/RNAseq_DESeq2_results/Results_with_annotation/vst_normalized_counts.csv", check.names = FALSE, row.names = 1)

# Check the dimensions of the datasets
dim(metabolomics_normalized)
dim(ptm_proteomics_normalized)
dim(transcriptomics_normalized)

# Check the first few rows of the datasets
# desired_order <- c("02N", "02R","05N", "05R", "11N", "11R",  # R1CC
#                    "07N", "07R", "09N", "09R", "12N", "12R", "13N", "13R",   # R1MD
#                    "21N", "21R", "34N", "34R", "43N", "43R",  # WD2CC
#                    "23N", "23R", "26N", "26R", "27N", "27R",  # WD2MD
#                    "01N", "01R", "04R", "06N", "06R", "08N", "08R",  # R1AA
#                    "10N", "14N", "15N", "15R", "16N", "16R"  # WD2AA
# )

# seperate them between age
desired_order <- c("01N", "01R", "04R", "06N", "06R", "08N", "08R",  # R1AA
                   "02N", "02R","05N", "05R", "11N", "11R",  # R1CC
                   "07N", "07R", "09N", "09R", "12N", "12R", "13N", "13R",   # R1MD
                   "10N", "14N", "15N", "15R", "16N", "16R",  # WD2AA
                   "21N", "21R", "34N", "34R", "43N", "43R",  # WD2CC
                   "23N", "23R", "26N", "26R", "27N", "27R"  # WD2MD
)

# matched the order
metabolomics_normalized <- reorder_and_fill_na(metabolomics_normalized, desired_order)
ptm_proteomics_normalized <- reorder_and_fill_na(ptm_proteomics_normalized, desired_order)
transcriptomics_normalized <- reorder_and_fill_na(transcriptomics_normalized, desired_order)

# Apply function to RNA-seq (transcriptomics)
rna_top_features <- select_high_variance_features(transcriptomics_normalized, variance_threshold = 0.8)

# Apply function to metabolomics data
metabolomics_top_features <- select_high_variance_features(metabolomics_normalized, variance_threshold = 0.8)

# Apply function to PTM proteomics data
ptm_top_features <- select_high_variance_features(ptm_proteomics_normalized, variance_threshold = 0.8)

# Combine all selected features into a list
selected_features_list <- list(
  RNAseq = rna_top_features,
  Metabolomics = metabolomics_top_features,
  PTM = ptm_top_features
)

# # Load the MOFA data with original data
# mofa_data <- list(
#   rna = as.matrix(transcriptomics_normalized),
#   ptm = as.matrix(ptm_proteomics_normalized),
#   metabolite = as.matrix(metabolomics_normalized)
# )

# Load the MOFA data with selected features
mofa_data <- list(
  rna = as.matrix(rna_top_features),
  ptm = as.matrix(ptm_top_features),
  metabolite = as.matrix(metabolomics_top_features)
)

#######################
# Create MOFA object ##
#######################
MOFAobject <- create_mofa(mofa_data)

# Visualise data structure
plot_data_overview(MOFAobject)


# # Multi-Group
# N = ncol(data[[1]])
# groups = c(rep("Yo-yo Diet",26), rep("Western Diet", 13))
# groups = c(rep("During",21), rep("After", 18))
# MOFAobject <- create_mofa(MOFA_yoyo_data, groups=groups)


####################
## Define options ##
####################

# Data options
# - scale_views: if views have very different ranges/variances, it is good practice to scale each view to unit variance (default is FALSE)
data_opts <- get_default_data_options(MOFAobject)

data_opts

# Model options
# - likelihoods: likelihood per view (options are "gaussian","poisson","bernoulli"). "gaussian" is used by default
# - num_factors: number of factors. By default K=10
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 9
model_opts

# Training options
# - maxiter: number of iterations
# - convergence_mode: "fast", "medium", "slow". For exploration, the fast mode is good enough.
# - drop_factor_threshold: minimum variance explained criteria to drop factors while training. Default is -1 (no dropping of factors)
# - gpu_mode: use GPU mode? This needs cupy installed and a functional GPU, see https://biofam.github.io/MOFA2/gpu_training.html
# - seed: random seed
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42

train_opts

#########################
## Prepare MOFA object ##
#########################

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)

##############################
## Train and Save the model ##
##############################

# Train the model
MOFAobject <- run_mofa(MOFAobject, outfile="./results/MOFA_multi-omics.hdf5")

# Save the model
saveRDS(MOFAobject,"./results/MOFA_multi-omics.rds")


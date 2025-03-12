MOFAobject <- load_model("./results/MOFA_multi-omics.hdf5")


# # need to configure the environment and add mofapy2
# reticulate::py_config()
# reticulate::py_install("mofapy2", envname = "r-reticulate")
# reticulate::use_virtualenv("r-reticulate", required = TRUE)

plot_factor_cor(MOFAobject)

plot_variance_explained(MOFAobject, max_r2=25)

plot_variance_explained(MOFAobject, plot_total = T)[[2]]

# plot_variance_explained_per_feature(MOFAobject, view = "metabolites", features = 5)
# calculate_variance_explained(MOFAobject)
# calculate_variance_explained_per_sample(MOFAobject)


# MOFAobject@cache$variance_explained$r2_per_factor$After

r2 <- MOFAobject@cache$variance_explained$r2_per_factor[[1]]


r2.dt <- r2 %>%
  as.data.table %>% .[,factor:=as.factor(1:MOFAobject@dimensions$K)] %>%
  melt(id.vars=c("factor"), variable.name="view", value.name = "r2") %>%
  .[,cum_r2:=cumsum(r2), by="view"]

ggline(r2.dt, x="factor", y="cum_r2", color="view") +
  labs(x="Factor number", y="Cumulative variance explained (%)") +
  theme(
    legend.title = element_blank(), 
    legend.position = "top",
    axis.text = element_text(size=rel(0.8))
  )



# load metadata for the data frame (like: phenotype or health records) 
library(openxlsx)
MOFA_yoyo_metadata <- read.xlsx("./data/raw/data_yoyo.xlsx", sheet = 6)
samples_metadata(MOFAobject) <- MOFA_yoyo_metadata


# load required packages
# install.packages("psych")
library(psych)

colnames(MOFA_yoyo_metadata)

correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("Diet", 
                                                 "Intervention",
                                                 "Timepoint"), 
                                  factors = "all",
                                  groups = "all",
                                  abs = FALSE,
                                  plot="r",
                                  alpha = 0.05,
                                  transpose = F,
                                  col = colorRampPalette(c("blue3", "white", "red3"))(200),
                                  # color
                                  
)

correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("R1AA", "R1CC", "R1MD",
                                                 "WD2AA", "WD2CC", "WD2MD"), 
                                  plot="r",
                                  transpose = F,
                                  col = colorRampPalette(c("blue3", "white", "red3"))(200),  # color
)


correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("Diet","Intervention","Timepoint",
                                                 "R1AA", "WD2AA", 
                                                 "R1CC", "R1MD",
                                                 "WD2CC", "WD2MD"), 
                                  factors = "all",
                                  plot="r",
                                  alpha = 0.05,
                                  transpose = F,
                                  col = colorRampPalette(c("blue3", "white", "red3"))(200),  # color
)



correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("Diet","Intervention","Timepoint"), 
                                  plot="log_pval",
                                  abs = T
)


correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("R1AA", "R1CC", "R1MD",
                                                 "WD2AA", "WD2CC", "WD2MD"), 
                                  plot="log_pval"
)

# ggsave("./results/MOFA_yoyo/Factors_correlated_with_covariates.svg")
# ggsave("./results/MOFA_yoyo/Factors_correlated_with_covariates.pdf")


### Factor Plots
plot_factor(MOFAobject, 
            factors = 5, 
            color_by = "Factor5"
)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Diet"
)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Timepoint"
)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Timepoint*Diet"
)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Intervention*Timepoint"
)


### Weights Plot

plot_weights(MOFAobject,
             view = "metabolite",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
             view = "ptm",
             factor = 2,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "rna",
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)


plot_weights(MOFAobject,
             view = "metabolite",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
             view = "metabolite",
             factor = 1,
             nfeatures = 20,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)



plot_factor(MOFAobject, 
            factors = 1, 
            group_by = "Timepoint*Diet",
            color_by = "Cysteine",
            add_violin = TRUE,
            dodge = TRUE
)


### Extract the data from model

# "factors" is a list of matrices, one matrix per group (samples × factors)
factors <- get_factors(MOFAobject, factors = "all")

# Check dimensions (should be number of samples × number of factors)
lapply(factors, dim)

weights <- get_weights(MOFAobject, views = "all", factors = "all")

# Check dimensions (should be number of features × number of factors)
lapply(weights, dim)

weights$rna

weights_df <- get_weights(MOFAobject, as.data.frame = TRUE)
head(weights_df)

data <- get_data(MOFAobject)

# Check dimensions for each omics view
lapply(data, function(x) lapply(x, dim))[[1]]

head(data$rna)

data_df <- get_data(MOFAobject, as.data.frame = TRUE)
head(data_df)
write.xlsx(data_df, "./results/MOFA_multi_omics_results/MOFA_data.xlsx")


# Function to extract top features for each factor

# Load the function from the script
source("./scripts/functions/MOFA_analysis_function.R")

# load required packages
library(openxlsx)

# Save the results to a directory (create if it doesn't exist)
dir.create("./results/MOFA_multiomics_results", showWarnings = FALSE)
output_dir <- "./results/MOFA_multiomics_results"
# Run function for Factors 1-9

mofa_results <- extract_top_features(weights_df, 
                                     factor_column = "factor", value_column = "value", 
                                     factors = 1:9, 
                                     output_dir = output_dir, 
                                     top_n = 1000)

# 
# 
# ### GSEA Gene Set Enrichment Analysis ###
# library(data.table)
# library(purrr)
# library(ggplot2)
# library(cowplot)
# # library(MOFAdata)
# # library(MOFA2)
# 
# data("reactomeGS")
# head(rownames(reactomeGS), n=3)
# 
# # C2: curated gene sets from online pathway databases, publications in PubMed, and knowledge of domain experts.
# data("MSigDB_v6.0_C2_human") 
# 
# # C5: extracted from the Gene Ontology data.base
# data("MSigDB_v6.0_C5_human") 
# 
# head(rownames(MSigDB_v6.0_C2_human), n=3)
# 
# 
# # C2: curated gene sets from online pathway databases, publications in PubMed, and knowledge of domain experts.
# data("MSigDB_v6.0_C2_mouse") 
# 
# # C5: extracted from the Gene Ontology data.base
# data("MSigDB_v6.0_C5_mouse") 
# 
# head(rownames(MSigDB_v6.0_C2_mouse), n=3)
# 
# features_names(MOFAobject)[["rna"]] <- toupper(features_names(MOFAobject)[["rna"]])
# head(features_names(MOFAobject)[["rna"]])
# 
# 
# enrichment.parametric <- run_enrichment(MOFAobject,
#                                         view = "rna", factors = 1:6,
#                                         feature.sets = MSigDB_v6.0_C5_mouse,
#                                         sign = "all",
#                                         statistical.test = "parametric"
# )
# 
# # 
# # run_enrichment(MOFAobject, 
# #                feature.sets = reactomeGS, 
# #                view = "mRNA",
# #                sign = "positive"
# # )
# 
# 
# names(enrichment.parametric)


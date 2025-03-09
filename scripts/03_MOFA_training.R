library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(ggpubr) # involve ggline


#######################
# Create MOFA object ##
#######################

MOFAobject <- create_mofa(MOFA_yoyo_data)

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

MOFAobject <- run_mofa(MOFAobject, outfile="./results/MOFA_yoyo_RNA+PTM++Metabolite.hdf5")
saveRDS(MOFAobject,"./results/MOFA_yoyo_RNA+PTM+Metabolites_splited.rds")

# # need to configure the environment and add mofapy2
# reticulate::py_config()
# reticulate::py_install("mofapy2", envname = "r-reticulate")
# reticulate::use_virtualenv("r-reticulate", required = TRUE)

plot_factor_cor(MOFAobject)

plot_variance_explained(MOFAobject, max_r2=15)

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
MOFA_yoyo_metadata <- read.xlsx("./data/data_yoyo.xlsx", sheet = 6)
samples_metadata(MOFAobject) <- MOFA_yoyo_metadata


category.colors <- c(
  "Healthy, no antibiotics" = "#66C2A5", 
  "Healthy, antibiotics" = "#8DA0CB",
  "Sepsis" = "#E78AC3",
  "Non septic ICU" = "#FC8D62"
)

p <- plot_factors(MOFAobject, 
                  factors = c(1,3), 
                  color_by = "Timepoint", 
                  dot_size = 4
                  ) + scale_fill_manual(values=category.colors)

p + 
  # geom_density_2d(aes_string(color="color_by")) +
  stat_ellipse(aes(color=color_by), geom = "polygon", alpha=0.25) +
  scale_color_manual(values=category.colors)



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
             factor = 5,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
             view = "ptm",
             factor = 5,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "rna",
                 factor = 5,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)


plot_weights(MOFAobject,
             view = "metabolites",
             factor = 2,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
             view = "metabolites",
             factor = 2,
             nfeatures = 10,     # Top number of features to highlight
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


# Get feature weights for Factor 3
factor5_weights <- subset(weights_df, factor == "Factor5")
unique(factor5_weights$view)

write.xlsx(factor5_weights, "./results/MOFA_multi_omics_results/Factor5_all_features.xlsx")

# Get top 10 features (genes/metabolites) contributing most to Factor 3
factor5_top_features <- factor5_weights[order(-abs(factor5_weights$value)), ][1:1000, ]
print(factor5_top_features)

write.xlsx(factor5_top_features, "./results/MOFA_multi_omics_results/Factor5_top_features.xlsx")




# Get feature weights for Factor 1
factor1_weights <- subset(weights_df, factor == "Factor1")
unique(factor6_weights$view)

write.xlsx(factor1_weights, "./results/MOFA_multi_omics_results/Factor1_all_features.xlsx")

# Get top 10 features (genes/metabolites) contributing most to Factor 3
factor1_top_features <- factor1_weights[order(-abs(factor1_weights$value)), ][1:1000, ]
print(factor1_top_features)

write.xlsx(factor1_top_features, "./results/MOFA_multi_omics_results/Factor1_top_features.xlsx")



# Get feature weights for Factor 1
factor6_weights <- subset(weights_df, factor == "Factor6")
unique(factor6_weights$view)

write.xlsx(factor6_weights, "./results/MOFA_multi_omics_results/Factor6_all_features.xlsx")

# Get top 10 features (genes/metabolites) contributing most to Factor 6
factor6_top_features <- factor6_weights[order(-abs(factor6_weights$value)), ][1:1000, ]
print(factor6_top_features)

write.xlsx(factor6_top_features, "./results/MOFA_multi_omics_results/Factor6_top_features.xlsx")



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


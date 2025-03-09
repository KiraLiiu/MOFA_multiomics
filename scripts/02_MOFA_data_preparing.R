library(MOFAdata)
library(data.table)


### Check example data
# mofa_sample_data <- fread("./data/mofa_example_data.txt")
# CLL_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")
# utils::data("CLL_data")  
# 
# head(CLL_data$Drugs)
# view(CLL_data$mRNA)
# 
# unique(mofa_sample_data$group)
# unique(mofa_sample_data$feature)
# unique(mofa_sample_data$view)
# 
# view(results_R1CC_WD1AA$normalized_data)

# read data
met_normalized <- fread("./data/metabolomics_results/data/normalized_data.csv")
# rna_normalized <- fread("./data/RNAseq_DESeq2_normalized_with_groups/RNAseq_normalized_counts_matrix.csv")
rna_normalized <- fread("./RNAseq_DESeq2_normalized_with_groups/RNAseq_HVFs_5000.csv")
ptm_normalized <- fread("./data/PTM/ptm_matrix_CLR_normalized.csv")

# 1️⃣ 先转换为 `data.frame` 以便正确处理行名
met_normalized <- as.data.frame(met_normalized)
rna_normalized <- as.data.frame(rna_normalized)
ptm_normalized <- as.data.frame(ptm_normalized)

# 2️⃣ 设置第一列为行名（保留原数据）
rownames(met_normalized) <- met_normalized[[1]]
rownames(rna_normalized) <- rna_normalized[[1]]
rownames(ptm_normalized) <- ptm_normalized[[1]]

# 3️⃣ 删除第一列，但保持行名
met_normalized <- met_normalized[, -1, drop = FALSE]
rna_normalized <- rna_normalized[, -1, drop = FALSE]
ptm_normalized <- ptm_normalized[, -1, drop = FALSE]


# 4️⃣ 只保留 RNA 里存在的样本（列名）
# common_samples <- intersect(colnames(met_normalized), colnames(rna_normalized))

# # 定义目标列顺序(如果需要)
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

# data.table 使用
# met_normalized <- met_normalized[, ..desired_order]
# rna_normalized <- rna_normalized[, ..desired_order]

# 对于 data.frame 使用
met_normalized <- met_normalized[, desired_order]
rna_normalized <- rna_normalized[, desired_order]
ptm_normalized <- ptm_normalized[, desired_order]

# 5️⃣ 确保 `Met` 和 `RNA` 样本一致
print(dim(met_normalized))  # 检查 Met 样本数
print(dim(rna_normalized))  # 检查 RNA 样本数
print(dim(ptm_normalized))
# 
# met_normalized <- as.matrix(met_normalized)
# rna_normalized <- as.matrix(rna_normalized)

### Prepare the MOFA dataset
MOFA_yoyo_data <- list(
  rna = as.matrix(rna_normalized),
  ptm = as.matrix(ptm_normalized),
  metabolite = as.matrix(met_normalized)
)

MOFA_yoyo_data <- list(
  mRNA = as.matrix(mrna_normalized),
  ncRNA = as.matrix(ncrna_normalized),
  ptm = as.matrix(ptm_matrix),
  metabolite = as.matrix(met_normalized)
)

MOFA_yoyo_data <- list(
  mRNA = as.matrix(mrna_normalized),
  ptm = as.matrix(ptm_matrix),
  metabolite = as.matrix(met_normalized)
)

MOFA_yoyo_data <- list(
  ncRNA = as.matrix(ncrna_normalized),
  ptm = as.matrix(ptm_matrix),
  metabolite = as.matrix(met_normalized)
)




lapply(MOFA_yoyo_data,dim)


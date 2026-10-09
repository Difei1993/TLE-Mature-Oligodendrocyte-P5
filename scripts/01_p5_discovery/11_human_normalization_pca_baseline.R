## ============================================================
## 11_human_normalization_pca_baseline.R
##
## 目的：
## 1. 读取QC和双细胞过滤后的58,341个细胞核
## 2. 删除仅在极少数细胞中表达的基因
## 3. LogNormalize标准化
## 4. 选择3,000个高变基因
## 5. ScaleData
## 6. 运行未整合的50-PC PCA
##
## 本阶段不做Harmony、scVI、UMAP或聚类
## ============================================================

options(stringsAsFactors = FALSE)
set.seed(20260724)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

## ------------------------------------------------------------
## 1. 读取QC后的singlet对象
## ------------------------------------------------------------

input_file <- file.path(
  "02_objects",
  "human_qc_singlets.rds"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "找不到输入对象：",
      input_file
    )
  )
}

cat("正在读取QC后singlet对象……\n")

human <- readRDS(input_file)

DefaultAssay(human) <- "RNA"

cat("输入细胞数：", ncol(human), "\n")
cat("输入基因数：", nrow(human), "\n")

if (ncol(human) != 58341) {
  warning(
    paste0(
      "预期58,341个细胞核，实际为",
      ncol(human),
      "。"
    )
  )
}

## ------------------------------------------------------------
## 2. 合并Seurat v5中的多个counts layers
##
## 8个样本合并后通常存在多个counts layer。
## JoinLayers只在当前新对象中合并，不修改此前保存的对象。
## ------------------------------------------------------------

cat("\n正在合并RNA assay的counts layers……\n")

human <- JoinLayers(
  object = human,
  assay = "RNA"
)

counts_matrix <- LayerData(
  object = human,
  assay = "RNA",
  layer = "counts"
)

if (!inherits(counts_matrix, "sparseMatrix")) {
  counts_matrix <- as(
    counts_matrix,
    "dgCMatrix"
  )
}

## 检查原始计数
if (any(counts_matrix@x < 0)) {
  stop("counts矩阵中存在负数。")
}

if (
  !all(
    counts_matrix@x ==
    round(counts_matrix@x)
  )
) {
  stop("counts矩阵不是整数计数。")
}

## ------------------------------------------------------------
## 3. 基因层面的基础过滤
##
## 保留至少在3个细胞核中检测到的基因。
## 这不是细胞QC，不会删除任何细胞。
## ------------------------------------------------------------

cat("\n正在计算每个基因被多少细胞检测到……\n")

n_cells_detected <- Matrix::rowSums(
  counts_matrix > 0
)

genes_keep <- names(
  n_cells_detected
)[
  n_cells_detected >= 3
]

gene_filter_summary <- data.frame(
  n_genes_before = nrow(human),
  n_genes_after = length(genes_keep),
  n_genes_removed =
    nrow(human) - length(genes_keep),
  minimum_cells_required = 3,
  stringsAsFactors = FALSE
)

cat("\n========== 基因过滤 ==========\n")
print(
  gene_filter_summary,
  row.names = FALSE
)

write.csv(
  gene_filter_summary,
  "04_tables/human_gene_filter_summary.csv",
  row.names = FALSE
)

gene_detection_table <- data.frame(
  gene = names(n_cells_detected),
  n_cells_detected =
    as.numeric(n_cells_detected),
  keep = names(n_cells_detected) %in%
    genes_keep,
  stringsAsFactors = FALSE
)

write.csv(
  gene_detection_table,
  "04_tables/human_gene_detection_counts.csv",
  row.names = FALSE
)

human <- subset(
  human,
  features = genes_keep
)

rm(
  counts_matrix,
  n_cells_detected
)

invisible(gc())

cat("基因过滤后基因数：", nrow(human), "\n")
cat("细胞数保持为：", ncol(human), "\n")

## ------------------------------------------------------------
## 4. LogNormalize标准化
##
## 每个细胞按总UMI校正到10,000，
## 再进行log1p转换。
## 不回归condition、donor或library。
## ------------------------------------------------------------

cat("\n正在执行NormalizeData……\n")

human <- NormalizeData(
  object = human,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

## ------------------------------------------------------------
## 5. 寻找高变基因
## ------------------------------------------------------------

cat("\n正在选择3,000个高变基因……\n")

human <- FindVariableFeatures(
  object = human,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = TRUE
)

variable_genes <- VariableFeatures(human)

if (length(variable_genes) != 3000) {
  warning(
    paste0(
      "预期3,000个高变基因，实际为",
      length(variable_genes),
      "。"
    )
  )
}

variable_gene_table <- data.frame(
  rank = seq_along(variable_genes),
  gene = variable_genes,
  is_mitochondrial = grepl(
    "^MT-",
    variable_genes
  ),
  is_ribosomal = grepl(
    "^(RPS|RPL)",
    variable_genes
  ),
  stringsAsFactors = FALSE
)

write.csv(
  variable_gene_table,
  "04_tables/human_variable_genes_3000.csv",
  row.names = FALSE
)

cat("\n前30个高变基因：\n")
print(
  head(
    variable_gene_table,
    30
  ),
  row.names = FALSE
)

## ------------------------------------------------------------
## 6. 高变基因图
## ------------------------------------------------------------

top20_variable <- head(
  variable_genes,
  20
)

p_variable <- VariableFeaturePlot(
  human
)

p_variable_labeled <- LabelPoints(
  plot = p_variable,
  points = top20_variable,
  repel = TRUE
)

ggsave(
  "03_figures/human_variable_features_top20.png",
  p_variable_labeled,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  "03_figures/human_variable_features_top20.pdf",
  p_variable_labeled,
  width = 10,
  height = 7
)

## ------------------------------------------------------------
## 7. ScaleData
##
## 仅对3,000个高变基因进行中心化和标准化。
## 暂不回归nCount、percent.mt、donor或condition。
## ------------------------------------------------------------

cat("\n正在执行ScaleData……\n")

human <- ScaleData(
  object = human,
  features = variable_genes,
  verbose = TRUE
)

## ------------------------------------------------------------
## 8. PCA
## ------------------------------------------------------------

cat("\n正在运行50个PC的PCA……\n")

human <- RunPCA(
  object = human,
  features = variable_genes,
  npcs = 50,
  seed.use = 20260724,
  verbose = TRUE
)

## ------------------------------------------------------------
## 9. 保存PCA标准差和解释度
## ------------------------------------------------------------

pca_stdev <- human[["pca"]]@stdev

pca_variance <- pca_stdev^2

pca_variance_percent <-
  100 * pca_variance /
  sum(pca_variance)

pca_cumulative_percent <-
  cumsum(pca_variance_percent)

pca_summary <- data.frame(
  PC = seq_along(pca_stdev),
  standard_deviation = pca_stdev,
  variance_percent =
    pca_variance_percent,
  cumulative_variance_percent =
    pca_cumulative_percent,
  stringsAsFactors = FALSE
)

write.csv(
  pca_summary,
  "04_tables/human_pca_variance_summary.csv",
  row.names = FALSE
)

cat("\n========== 前20个PC ==========\n")
print(
  head(
    pca_summary,
    20
  ),
  row.names = FALSE
)

## ------------------------------------------------------------
## 10. Elbow plot
## ------------------------------------------------------------

p_elbow <- ElbowPlot(
  human,
  ndims = 50
) +
  labs(
    title = "PCA elbow plot: unintegrated human snRNA-seq"
  )

ggsave(
  "03_figures/human_pca_elbow_50pcs.png",
  p_elbow,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  "03_figures/human_pca_elbow_50pcs.pdf",
  p_elbow,
  width = 9,
  height = 6
)

## ------------------------------------------------------------
## 11. PCA按样本和condition着色
## ------------------------------------------------------------

p_pca_sample <- DimPlot(
  object = human,
  reduction = "pca",
  dims = c(1, 2),
  group.by = "sample_id",
  raster = TRUE
) +
  labs(
    title = "Unintegrated PCA by sample"
  )

p_pca_condition <- DimPlot(
  object = human,
  reduction = "pca",
  dims = c(1, 2),
  group.by = "condition",
  raster = TRUE
) +
  labs(
    title = "Unintegrated PCA by condition"
  )

p_pca_combined <- (
  p_pca_sample |
    p_pca_condition
)

ggsave(
  "03_figures/human_pca_by_sample_and_condition.png",
  p_pca_combined,
  width = 15,
  height = 7,
  dpi = 300
)

ggsave(
  "03_figures/human_pca_by_sample_and_condition.pdf",
  p_pca_combined,
  width = 15,
  height = 7
)

## ------------------------------------------------------------
## 12. PCA与技术变量关系
## ------------------------------------------------------------

p_pca_mt <- FeaturePlot(
  object = human,
  features = "percent.mt",
  reduction = "pca",
  dims = c(1, 2),
  raster = TRUE
) +
  labs(
    title = "PCA colored by percent.mt"
  )

p_pca_count <- FeaturePlot(
  object = human,
  features = "nCount_RNA",
  reduction = "pca",
  dims = c(1, 2),
  raster = TRUE
) +
  labs(
    title = "PCA colored by nCount_RNA"
  )

p_pca_feature <- FeaturePlot(
  object = human,
  features = "nFeature_RNA",
  reduction = "pca",
  dims = c(1, 2),
  raster = TRUE
) +
  labs(
    title = "PCA colored by nFeature_RNA"
  )

p_pca_technical <- (
  p_pca_mt |
    p_pca_count |
    p_pca_feature
)

ggsave(
  "03_figures/human_pca_technical_variables.png",
  p_pca_technical,
  width = 18,
  height = 6,
  dpi = 300
)

ggsave(
  "03_figures/human_pca_technical_variables.pdf",
  p_pca_technical,
  width = 18,
  height = 6
)

## ------------------------------------------------------------
## 13. PC loading
## ------------------------------------------------------------

p_loadings <- VizDimLoadings(
  object = human,
  dims = 1:5,
  reduction = "pca",
  nfeatures = 20,
  balanced = TRUE
)

ggsave(
  "03_figures/human_pca_loadings_pc1_to_pc5.png",
  p_loadings,
  width = 16,
  height = 10,
  dpi = 300
)

ggsave(
  "03_figures/human_pca_loadings_pc1_to_pc5.pdf",
  p_loadings,
  width = 16,
  height = 10
)

## ------------------------------------------------------------
## 14. 保存对象
## ------------------------------------------------------------

output_file <- file.path(
  "02_objects",
  "human_log_normalized_pca_baseline.rds"
)

saveRDS(
  human,
  file = output_file,
  compress = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  "05_logs/11_human_normalization_pca_baseline_sessionInfo.txt"
)

cat("\n========================================\n")
cat("标准化和未整合PCA完成。\n")
cat("细胞数：", ncol(human), "\n")
cat("基因数：", nrow(human), "\n")
cat("高变基因数：", length(variable_genes), "\n")
cat("PCA维数：", ncol(Embeddings(human, "pca")), "\n")
cat("对象保存于：\n")
cat(output_file, "\n")
cat("\n重要图片：\n")
cat("03_figures/human_variable_features_top20.png\n")
cat("03_figures/human_pca_elbow_50pcs.png\n")
cat("03_figures/human_pca_by_sample_and_condition.png\n")
cat("03_figures/human_pca_technical_variables.png\n")
cat("03_figures/human_pca_loadings_pc1_to_pc5.png\n")
cat("========================================\n")
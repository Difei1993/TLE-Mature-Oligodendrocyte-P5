## ============================================================
## 13_human_harmony_umap_clustering.R
##
## 目的：
## 1. 使用library_id对PCA空间进行Harmony校正
## 2. condition不参与批次校正
## 3. 在Harmony空间运行UMAP
## 4. 测试resolution 0.2、0.4、0.6
## 5. 比较未整合与Harmony结果的样本主导程度
##
## Harmony仅用于聚类和注释敏感性分析。
## 正式差异表达不得使用Harmony表达值。
## ============================================================

options(stringsAsFactors = FALSE)
set.seed(20260724)

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

## ------------------------------------------------------------
## 1. 读取未整合UMAP和聚类对象
## ------------------------------------------------------------

input_file <- file.path(
  "02_objects",
  "human_unintegrated_umap_multiresolution.rds"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "找不到输入对象：",
      input_file
    )
  )
}

cat("正在读取未整合对象……\n")

human <- readRDS(input_file)

DefaultAssay(human) <- "RNA"

required_metadata <- c(
  "library_id",
  "donor_id",
  "sample_id",
  "condition"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(human@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    paste(
      "缺少metadata字段：",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

cat("细胞数：", ncol(human), "\n")
cat("基因数：", nrow(human), "\n")

cat("\nLibrary数量：\n")
print(
  table(human$library_id)
)

cat("\nCondition数量：\n")
print(
  table(human$condition)
)

## 明确检查condition没有被当作batch
if (identical(
  as.character(human$library_id),
  as.character(human$condition)
)) {
  stop("library_id与condition完全相同，请检查metadata。")
}

## ------------------------------------------------------------
## 2. 运行Harmony
##
## 只校正library_id。
## 不校正condition、donor_id、nCount或nFeature。
## ------------------------------------------------------------

cat("\n正在运行Harmony……\n")

human <- RunHarmony(
  object = human,
  group.by.vars = "library_id",
  reduction.use = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  assay.use = "RNA",
  project.dim = FALSE,
  plot_convergence = TRUE,
  verbose = TRUE
)

if (!"harmony" %in% Reductions(human)) {
  stop("Harmony reduction未成功生成。")
}

cat(
  "Harmony维数：",
  ncol(Embeddings(human, "harmony")),
  "\n"
)

## ------------------------------------------------------------
## 3. Harmony邻居图
## ------------------------------------------------------------

cat("\n正在构建Harmony邻居图……\n")

human <- FindNeighbors(
  object = human,
  reduction = "harmony",
  dims = 1:30,
  k.param = 20,
  graph.name = c(
    "harmony30_nn",
    "harmony30_snn"
  ),
  verbose = TRUE
)

## ------------------------------------------------------------
## 4. Harmony UMAP
## ------------------------------------------------------------

cat("\n正在运行Harmony UMAP……\n")

human <- RunUMAP(
  object = human,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3,
  spread = 1,
  seed.use = 20260724,
  reduction.name = "umap.harmony",
  reduction.key = "UMAPHARM_",
  verbose = TRUE
)

## ------------------------------------------------------------
## 5. 测试三个resolution
## ------------------------------------------------------------

resolutions <- c(
  0.2,
  0.4,
  0.6
)

cluster_columns <- character(
  length(resolutions)
)

for (i in seq_along(resolutions)) {
  
  resolution_value <- resolutions[i]
  
  cluster_column <- paste0(
    "cluster_harmony30_r",
    gsub(
      "\\.",
      "p",
      as.character(resolution_value)
    )
  )
  
  cat(
    "\n运行Harmony resolution = ",
    resolution_value,
    "\n",
    sep = ""
  )
  
  human <- FindClusters(
    object = human,
    graph.name = "harmony30_snn",
    resolution = resolution_value,
    algorithm = 1,
    random.seed = 20260724,
    cluster.name = cluster_column,
    verbose = TRUE
  )
  
  cluster_columns[i] <- cluster_column
}

## ------------------------------------------------------------
## 6. 计算各resolution cluster构成
## ------------------------------------------------------------

resolution_summary_list <- list()
cluster_composition_list <- list()

for (i in seq_along(resolutions)) {
  
  resolution_value <- resolutions[i]
  cluster_column <- cluster_columns[i]
  
  cluster_values <- as.character(
    human@meta.data[[cluster_column]]
  )
  
  cluster_sizes <- table(
    cluster_values
  )
  
  resolution_summary_list[[i]] <- data.frame(
    resolution = resolution_value,
    cluster_column = cluster_column,
    n_clusters = length(cluster_sizes),
    smallest_cluster_cells = min(cluster_sizes),
    median_cluster_cells = median(cluster_sizes),
    largest_cluster_cells = max(cluster_sizes),
    stringsAsFactors = FALSE
  )
  
  temporary_data <- human@meta.data %>%
    tibble::rownames_to_column(
      "cell_id"
    ) %>%
    mutate(
      cluster = cluster_values
    )
  
  cluster_composition <- temporary_data %>%
    group_by(
      cluster
    ) %>%
    summarise(
      n_cells = n(),
      
      n_donors_present =
        n_distinct(donor_id),
      
      n_samples_present =
        n_distinct(sample_id),
      
      max_donor_fraction =
        max(table(donor_id)) / n(),
      
      max_sample_fraction =
        max(table(sample_id)) / n(),
      
      normal_fraction =
        mean(condition == "Normal"),
      
      tle_fraction =
        mean(condition == "TLE"),
      
      .groups = "drop"
    ) %>%
    mutate(
      resolution = resolution_value,
      cluster_column = cluster_column
    )
  
  cluster_composition_list[[i]] <-
    cluster_composition
}

harmony_resolution_summary <- bind_rows(
  resolution_summary_list
)

harmony_cluster_composition <- bind_rows(
  cluster_composition_list
)

write.csv(
  harmony_resolution_summary,
  "04_tables/human_harmony_resolution_summary.csv",
  row.names = FALSE
)

write.csv(
  harmony_cluster_composition,
  "04_tables/human_harmony_cluster_composition.csv",
  row.names = FALSE
)

## ------------------------------------------------------------
## 7. 聚类质量摘要
## ------------------------------------------------------------

harmony_quality_summary <- harmony_cluster_composition %>%
  group_by(
    resolution
  ) %>%
  summarise(
    n_clusters = n(),
    
    clusters_all_8_donors =
      sum(n_donors_present == 8),
    
    clusters_at_least_6_donors =
      sum(n_donors_present >= 6),
    
    clusters_fewer_than_4_donors =
      sum(n_donors_present < 4),
    
    clusters_sample_dominant_50 =
      sum(max_sample_fraction >= 0.50),
    
    clusters_sample_dominant_60 =
      sum(max_sample_fraction >= 0.60),
    
    clusters_sample_dominant_80 =
      sum(max_sample_fraction >= 0.80),
    
    clusters_condition_dominant_90 =
      sum(
        normal_fraction >= 0.90 |
          tle_fraction >= 0.90
      ),
    
    clusters_condition_dominant_95 =
      sum(
        normal_fraction >= 0.95 |
          tle_fraction >= 0.95
      ),
    
    clusters_under_300_cells =
      sum(n_cells < 300),
    
    clusters_under_500_cells =
      sum(n_cells < 500),
    
    .groups = "drop"
  )

write.csv(
  harmony_quality_summary,
  "04_tables/human_harmony_cluster_quality_summary.csv",
  row.names = FALSE
)

cat("\n========== Harmony分辨率 ==========\n")

print(
  harmony_resolution_summary,
  row.names = FALSE
)

cat("\n========== Harmony cluster质量 ==========\n")

print(
  harmony_quality_summary,
  row.names = FALSE
)

## ------------------------------------------------------------
## 8. UMAP按样本和condition
## ------------------------------------------------------------

p_sample <- DimPlot(
  object = human,
  reduction = "umap.harmony",
  group.by = "sample_id",
  shuffle = TRUE,
  seed = 20260724,
  raster = TRUE
) +
  labs(
    title = "Harmony UMAP by sample"
  )

p_condition <- DimPlot(
  object = human,
  reduction = "umap.harmony",
  group.by = "condition",
  shuffle = TRUE,
  seed = 20260724,
  raster = TRUE
) +
  labs(
    title = "Harmony UMAP by condition"
  )

p_combined <- (
  p_sample |
    p_condition
)

ggsave(
  "03_figures/human_harmony_umap_by_sample_condition.png",
  p_combined,
  width = 16,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  "03_figures/human_harmony_umap_by_sample_condition.pdf",
  p_combined,
  width = 16,
  height = 7
)

## ------------------------------------------------------------
## 9. Harmony UMAP按样本拆分
## ------------------------------------------------------------

p_split <- DimPlot(
  object = human,
  reduction = "umap.harmony",
  group.by = "condition",
  split.by = "sample_id",
  ncol = 4,
  raster = TRUE
) +
  labs(
    title = "Harmony UMAP split by sample"
  )

ggsave(
  "03_figures/human_harmony_umap_split_by_sample.png",
  p_split,
  width = 18,
  height = 9,
  dpi = 300,
  bg = "white"
)

## ------------------------------------------------------------
## 10. 多resolution聚类图
## ------------------------------------------------------------

cluster_plots <- lapply(
  seq_along(cluster_columns),
  function(i) {
    
    DimPlot(
      object = human,
      reduction = "umap.harmony",
      group.by = cluster_columns[i],
      label = TRUE,
      repel = TRUE,
      raster = TRUE
    ) +
      NoLegend() +
      labs(
        title = paste0(
          "Harmony resolution ",
          resolutions[i],
          " (",
          harmony_resolution_summary$n_clusters[i],
          " clusters)"
        )
      )
  }
)

p_cluster_grid <- wrap_plots(
  cluster_plots,
  ncol = 2
)

ggsave(
  "03_figures/human_harmony_umap_resolution_grid.png",
  p_cluster_grid,
  width = 15,
  height = 12,
  dpi = 300,
  bg = "white"
)

ggsave(
  "03_figures/human_harmony_umap_resolution_grid.pdf",
  p_cluster_grid,
  width = 15,
  height = 12
)

## ------------------------------------------------------------
## 11. 保存对象
## ------------------------------------------------------------

output_file <- file.path(
  "02_objects",
  "human_harmony_umap_multiresolution.rds"
)

saveRDS(
  human,
  output_file,
  compress = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  "05_logs/13_human_harmony_umap_clustering_sessionInfo.txt"
)

cat("\n========================================\n")
cat("Harmony整合敏感性分析完成。\n")
cat("condition没有作为batch进行校正。\n")
cat("Harmony只用于聚类和注释。\n")
cat("对象保存于：\n")
cat(output_file, "\n")
cat("\n主要图片：\n")
cat("03_figures/human_harmony_umap_by_sample_condition.png\n")
cat("03_figures/human_harmony_umap_split_by_sample.png\n")
cat("03_figures/human_harmony_umap_resolution_grid.png\n")
cat("========================================\n")

View(harmony_resolution_summary)
View(harmony_quality_summary)
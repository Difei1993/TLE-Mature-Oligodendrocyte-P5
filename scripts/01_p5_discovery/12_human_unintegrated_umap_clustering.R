## ============================================================
## 12_human_unintegrated_umap_clustering.R
##
## 目的：
## 1. 使用未整合PCA前30个PC构建邻居图
## 2. 运行未整合UMAP
## 3. 测试多个Leiden/Louvain聚类分辨率
## 4. 检查cluster数量及供者/条件构成
##
## 本阶段不进行细胞注释
## ============================================================

options(stringsAsFactors = FALSE)
set.seed(20260724)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

## ------------------------------------------------------------
## 1. 读取标准化和PCA对象
## ------------------------------------------------------------

input_file <- file.path(
  "02_objects",
  "human_log_normalized_pca_baseline.rds"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "找不到输入对象：",
      input_file
    )
  )
}

cat("正在读取PCA对象……\n")

human <- readRDS(input_file)

DefaultAssay(human) <- "RNA"

cat("细胞数：", ncol(human), "\n")
cat("基因数：", nrow(human), "\n")
cat(
  "已有PCA维数：",
  ncol(Embeddings(human, "pca")),
  "\n"
)

if (
  ncol(Embeddings(human, "pca")) < 30
) {
  stop("PCA维数不足30，请检查上一步。")
}

## ------------------------------------------------------------
## 2. 使用前30个PC构建邻居图
## ------------------------------------------------------------

dims_use <- 1:30

cat("\n正在使用前30个PC构建邻居图……\n")

human <- FindNeighbors(
  object = human,
  reduction = "pca",
  dims = dims_use,
  k.param = 20,
  graph.name = c(
    "pca30_nn",
    "pca30_snn"
  ),
  verbose = TRUE
)

## ------------------------------------------------------------
## 3. 运行未整合UMAP
## ------------------------------------------------------------

cat("\n正在运行未整合UMAP……\n")

human <- RunUMAP(
  object = human,
  reduction = "pca",
  dims = dims_use,
  n.neighbors = 30,
  min.dist = 0.3,
  spread = 1,
  seed.use = 20260724,
  reduction.name = "umap.unintegrated",
  reduction.key = "UMAPUNINT_",
  verbose = TRUE
)

## ------------------------------------------------------------
## 4. 测试多个聚类分辨率
## ------------------------------------------------------------

resolutions <- c(
  0.2,
  0.4,
  0.6,
  0.8,
  1.0
)

cluster_columns <- character(
  length(resolutions)
)

cat("\n正在测试多个聚类分辨率……\n")

for (i in seq_along(resolutions)) {
  
  resolution_value <- resolutions[i]
  
  cluster_column <- paste0(
    "cluster_pca30_r",
    gsub(
      "\\.",
      "p",
      as.character(resolution_value)
    )
  )
  
  cat(
    "运行resolution = ",
    resolution_value,
    "\n",
    sep = ""
  )
  
  human <- FindClusters(
    object = human,
    graph.name = "pca30_snn",
    resolution = resolution_value,
    algorithm = 1,
    random.seed = 20260724,
    cluster.name = cluster_column,
    verbose = TRUE
  )
  
  cluster_columns[i] <- cluster_column
}

## ------------------------------------------------------------
## 5. 生成各resolution的cluster数量汇总
## ------------------------------------------------------------

resolution_summary_list <- vector(
  mode = "list",
  length = length(resolutions)
)

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

resolution_summary <- bind_rows(
  resolution_summary_list
)

cluster_composition_all <- bind_rows(
  cluster_composition_list
)

write.csv(
  resolution_summary,
  "04_tables/human_unintegrated_resolution_summary.csv",
  row.names = FALSE
)

write.csv(
  cluster_composition_all,
  "04_tables/human_unintegrated_cluster_composition_all_resolutions.csv",
  row.names = FALSE
)

cat("\n========== 分辨率汇总 ==========\n")

print(
  resolution_summary,
  row.names = FALSE
)

## ------------------------------------------------------------
## 6. 绘制按样本和condition着色的UMAP
## ------------------------------------------------------------

p_umap_sample <- DimPlot(
  object = human,
  reduction = "umap.unintegrated",
  group.by = "sample_id",
  raster = TRUE,
  shuffle = TRUE,
  seed = 20260724
) +
  labs(
    title = "Unintegrated UMAP by sample"
  )

p_umap_condition <- DimPlot(
  object = human,
  reduction = "umap.unintegrated",
  group.by = "condition",
  raster = TRUE,
  shuffle = TRUE,
  seed = 20260724
) +
  labs(
    title = "Unintegrated UMAP by condition"
  )

p_umap_sample_condition <- (
  p_umap_sample |
    p_umap_condition
)

ggsave(
  "03_figures/human_unintegrated_umap_by_sample_condition.png",
  p_umap_sample_condition,
  width = 16,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  "03_figures/human_unintegrated_umap_by_sample_condition.pdf",
  p_umap_sample_condition,
  width = 16,
  height = 7
)

## ------------------------------------------------------------
## 7. 绘制各resolution聚类图
## ------------------------------------------------------------

cluster_plots <- lapply(
  seq_along(cluster_columns),
  function(i) {
    
    resolution_value <- resolutions[i]
    cluster_column <- cluster_columns[i]
    
    DimPlot(
      object = human,
      reduction = "umap.unintegrated",
      group.by = cluster_column,
      label = TRUE,
      repel = TRUE,
      raster = TRUE
    ) +
      NoLegend() +
      labs(
        title = paste0(
          "Resolution ",
          resolution_value,
          " (",
          resolution_summary$n_clusters[i],
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
  "03_figures/human_unintegrated_umap_resolution_grid.png",
  p_cluster_grid,
  width = 15,
  height = 18,
  dpi = 300,
  bg = "white"
)

ggsave(
  "03_figures/human_unintegrated_umap_resolution_grid.pdf",
  p_cluster_grid,
  width = 15,
  height = 18
)

## ------------------------------------------------------------
## 8. 按样本拆分UMAP
## ------------------------------------------------------------

p_umap_split_sample <- DimPlot(
  object = human,
  reduction = "umap.unintegrated",
  group.by = "condition",
  split.by = "sample_id",
  ncol = 4,
  raster = TRUE
) +
  labs(
    title = "Unintegrated UMAP split by sample"
  )

ggsave(
  "03_figures/human_unintegrated_umap_split_by_sample.png",
  p_umap_split_sample,
  width = 18,
  height = 9,
  dpi = 300,
  bg = "white"
)

## ------------------------------------------------------------
## 9. 保存对象
## ------------------------------------------------------------

output_file <- file.path(
  "02_objects",
  "human_unintegrated_umap_multiresolution.rds"
)

saveRDS(
  human,
  file = output_file,
  compress = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  "05_logs/12_human_unintegrated_umap_clustering_sessionInfo.txt"
)

cat("\n========================================\n")
cat("未整合UMAP和多分辨率聚类完成。\n")
cat("使用PC：1–30\n")
cat("细胞数：", ncol(human), "\n")
cat("对象保存于：\n")
cat(output_file, "\n")
cat("\n主要图片：\n")
cat(
  "03_figures/human_unintegrated_umap_by_sample_condition.png\n"
)
cat(
  "03_figures/human_unintegrated_umap_resolution_grid.png\n"
)
cat(
  "03_figures/human_unintegrated_umap_split_by_sample.png\n"
)
cat("========================================\n")

View(resolution_summary)
View(cluster_composition_all)
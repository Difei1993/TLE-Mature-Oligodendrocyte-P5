## ============================================================
## 07_human_doublet_detection.R
##
## 目的：
## 1. 对8个样本分别运行scDblFinder
## 2. 将双细胞评分加入Seurat对象
## 3. 本阶段只标记双细胞，不删除
## ============================================================

options(stringsAsFactors = FALSE)
set.seed(20260724)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
  library(dplyr)
  library(ggplot2)
})

## ------------------------------------------------------------
## 1. 文件和目录
## ------------------------------------------------------------

input_file <- file.path(
  "02_objects",
  "human_stage1_qc_filtered.rds"
)

output_object <- file.path(
  "02_objects",
  "human_stage1_with_doublet_scores.rds"
)

per_sample_dir <- file.path(
  "04_tables",
  "doublet_results_by_sample"
)

dir.create(
  per_sample_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "未找到输入对象：",
      input_file
    )
  )
}

## ------------------------------------------------------------
## 2. 读取第一阶段QC对象
## ------------------------------------------------------------

cat("正在读取第一阶段QC对象……\n")

human_stage1 <- readRDS(
  input_file
)

cat(
  "细胞数：",
  ncol(human_stage1),
  "\n"
)

cat(
  "基因数：",
  nrow(human_stage1),
  "\n"
)

if (
  !"sample_id" %in%
  colnames(human_stage1@meta.data)
) {
  stop("Seurat对象中缺少sample_id字段。")
}

sample_ids <- sort(
  unique(
    as.character(
      human_stage1$sample_id
    )
  )
)

cat("\n检测到样本数：", length(sample_ids), "\n")
print(sample_ids)

if (length(sample_ids) != 8) {
  stop(
    paste0(
      "预期8个样本，实际检测到",
      length(sample_ids),
      "个。"
    )
  )
}

## ------------------------------------------------------------
## 3. 逐样本运行scDblFinder
## ------------------------------------------------------------

all_doublet_results <- vector(
  mode = "list",
  length = length(sample_ids)
)

names(all_doublet_results) <- sample_ids

doublet_summary_list <- vector(
  mode = "list",
  length = length(sample_ids)
)

for (i in seq_along(sample_ids)) {
  
  sample_id_value <- sample_ids[i]
  
  cat(
    "\n========================================\n",
    "正在检测样本：",
    sample_id_value,
    "（",
    i,
    "/",
    length(sample_ids),
    "）\n",
    sep = ""
  )
  
  ## 提取单个样本
  sample_object <- subset(
    human_stage1,
    subset = sample_id == sample_id_value
  )
  
  cat(
    "样本细胞数：",
    ncol(sample_object),
    "\n"
  )
  
  ## Seurat v5合并后的对象可能包含多个counts layer
  ## 对当前单个样本重新合并layer
  sample_object <- JoinLayers(
    sample_object,
    assay = "RNA"
  )
  
  ## 提取原始整数计数
  counts_matrix <- LayerData(
    object = sample_object,
    assay = "RNA",
    layer = "counts"
  )
  
  if (!inherits(counts_matrix, "sparseMatrix")) {
    counts_matrix <- as(
      counts_matrix,
      "dgCMatrix"
    )
  }
  
  ## 基本检查
  if (any(counts_matrix@x < 0)) {
    stop(
      paste0(
        sample_id_value,
        "的counts矩阵存在负值。"
      )
    )
  }
  
  if (
    !all(
      counts_matrix@x ==
      round(counts_matrix@x)
    )
  ) {
    stop(
      paste0(
        sample_id_value,
        "的counts矩阵不是整数计数。"
      )
    )
  }
  
  ## 构建SingleCellExperiment
  sce <- SingleCellExperiment(
    assays = list(
      counts = counts_matrix
    )
  )
  
  colnames(sce) <- colnames(sample_object)
  
  colData(sce)$sample_id <- sample_id_value
  
  ## 每个样本使用不同但固定的随机种子
  set.seed(
    20260724 + i
  )
  
  ## clusters = FALSE：
  ## 使用随机人工双细胞模式，
  ## 避免依赖尚未完成的细胞聚类和注释
  sce <- scDblFinder(
    sce,
    clusters = FALSE,
    BPPARAM = SerialParam()
  )
  
  ## 提取所有scDblFinder字段
  doublet_columns <- grep(
    "^scDblFinder",
    colnames(colData(sce)),
    value = TRUE
  )
  
  sample_result <- as.data.frame(
    colData(sce)[
      ,
      doublet_columns,
      drop = FALSE
    ]
  )
  
  sample_result$cell_id <- rownames(
    sample_result
  )
  
  sample_result$sample_id <- sample_id_value
  
  ## 将cell_id放到第一列
  sample_result <- sample_result[
    ,
    c(
      "cell_id",
      "sample_id",
      setdiff(
        colnames(sample_result),
        c(
          "cell_id",
          "sample_id"
        )
      )
    )
  ]
  
  ## 保存单个样本结果
  write.csv(
    sample_result,
    file = file.path(
      per_sample_dir,
      paste0(
        sample_id_value,
        "_scDblFinder_results.csv"
      )
    ),
    row.names = FALSE
  )
  
  all_doublet_results[[i]] <- sample_result
  
  n_doublets <- sum(
    sample_result$scDblFinder.class ==
      "doublet"
  )
  
  n_singlets <- sum(
    sample_result$scDblFinder.class ==
      "singlet"
  )
  
  doublet_summary_list[[i]] <- data.frame(
    sample_id = sample_id_value,
    n_cells = nrow(sample_result),
    n_singlets = n_singlets,
    n_doublets = n_doublets,
    doublet_rate_percent =
      100 * n_doublets / nrow(sample_result),
    median_doublet_score = median(
      sample_result$scDblFinder.score,
      na.rm = TRUE
    ),
    max_doublet_score = max(
      sample_result$scDblFinder.score,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
  
  cat(
    "识别双细胞：",
    n_doublets,
    " / ",
    nrow(sample_result),
    "（",
    round(
      100 * n_doublets /
        nrow(sample_result),
      2
    ),
    "%）\n",
    sep = ""
  )
  
  rm(
    sample_object,
    counts_matrix,
    sce,
    sample_result
  )
  
  invisible(gc())
}

## ------------------------------------------------------------
## 4. 合并8个样本结果
## ------------------------------------------------------------

doublet_results <- bind_rows(
  all_doublet_results
)

doublet_summary <- bind_rows(
  doublet_summary_list
)

write.csv(
  doublet_results,
  "04_tables/human_scDblFinder_results_all_cells.csv",
  row.names = FALSE
)

write.csv(
  doublet_summary,
  "04_tables/human_scDblFinder_summary_by_sample.csv",
  row.names = FALSE
)

cat("\n========== 双细胞检测汇总 ==========\n")
print(
  doublet_summary,
  row.names = FALSE
)

## ------------------------------------------------------------
## 5. 将结果加入原Seurat对象
## ------------------------------------------------------------

match_index <- match(
  colnames(human_stage1),
  doublet_results$cell_id
)

if (anyNA(match_index)) {
  stop(
    paste0(
      "有",
      sum(is.na(match_index)),
      "个Seurat细胞没有双细胞结果。"
    )
  )
}

ordered_results <- doublet_results[
  match_index,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    ordered_results$cell_id,
    colnames(human_stage1)
  )
)

doublet_metadata <- ordered_results[
  ,
  c(
    "scDblFinder.score",
    "scDblFinder.class"
  ),
  drop = FALSE
]

rownames(doublet_metadata) <- ordered_results$cell_id

human_stage1 <- AddMetaData(
  object = human_stage1,
  metadata = doublet_metadata
)

## ------------------------------------------------------------
## 6. 比较双细胞与之前的高复杂度标记
## ------------------------------------------------------------

comparison_table <- human_stage1@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  group_by(
    sample_id,
    scDblFinder.class
  ) %>%
  summarise(
    n_cells = n(),
    
    median_nCount_RNA =
      median(nCount_RNA),
    
    median_nFeature_RNA =
      median(nFeature_RNA),
    
    percent_high_feature =
      100 * mean(high_feature_flag),
    
    percent_high_count =
      100 * mean(high_count_flag),
    
    .groups = "drop"
  )

write.csv(
  comparison_table,
  "04_tables/human_doublet_vs_high_complexity_flags.csv",
  row.names = FALSE
)

## ------------------------------------------------------------
## 7. 绘图
## ------------------------------------------------------------

plot_data <- human_stage1@meta.data %>%
  tibble::rownames_to_column("cell_id")

p_rate <- ggplot(
  doublet_summary,
  aes(
    x = sample_id,
    y = doublet_rate_percent,
    fill = grepl(
      "TLE",
      sample_id
    )
  )
) +
  geom_col() +
  geom_text(
    aes(
      label = sprintf(
        "%.1f%%",
        doublet_rate_percent
      )
    ),
    vjust = -0.25,
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = "#F8766D",
      "TRUE" = "#00BFC4"
    ),
    labels = c(
      "FALSE" = "Normal",
      "TRUE" = "TLE"
    )
  ) +
  labs(
    title = "scDblFinder doublet rate by sample",
    x = NULL,
    y = "Predicted doublets (%)",
    fill = "Condition"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  "03_figures/human_doublet_rate_by_sample.png",
  p_rate,
  width = 11,
  height = 6,
  dpi = 300
)

ggsave(
  "03_figures/human_doublet_rate_by_sample.pdf",
  p_rate,
  width = 11,
  height = 6
)

p_score <- ggplot(
  plot_data,
  aes(
    x = scDblFinder.score,
    fill = scDblFinder.class
  )
) +
  geom_histogram(
    bins = 80,
    position = "identity",
    alpha = 0.7
  ) +
  facet_wrap(
    ~ sample_id,
    ncol = 4,
    scales = "free_y"
  ) +
  labs(
    title = "scDblFinder score distributions",
    x = "scDblFinder score",
    y = "Number of nuclei",
    fill = "Classification"
  ) +
  theme_classic()

ggsave(
  "03_figures/human_doublet_score_distribution.png",
  p_score,
  width = 14,
  height = 8,
  dpi = 300
)

ggsave(
  "03_figures/human_doublet_score_distribution.pdf",
  p_score,
  width = 14,
  height = 8
)

## ------------------------------------------------------------
## 8. 保存带双细胞标记的对象
## 注意：目前尚未删除预测双细胞
## ------------------------------------------------------------

saveRDS(
  human_stage1,
  output_object,
  compress = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  "05_logs/07_human_doublet_detection_sessionInfo.txt"
)

cat("\n========================================\n")
cat("双细胞检测完成。\n")
cat("目前没有删除任何预测双细胞。\n")
cat("带评分对象保存于：\n")
cat(output_object, "\n")
cat("汇总表：\n")
cat("04_tables/human_scDblFinder_summary_by_sample.csv\n")
cat("主要图片：\n")
cat("03_figures/human_doublet_rate_by_sample.png\n")
cat("03_figures/human_doublet_score_distribution.png\n")
cat("========================================\n")

View(doublet_summary)
## ============================================================
## 02_human_qc_overview.R
## 目的：查看8个人类样本的原始QC分布
## 本阶段只绘图和汇总，不过滤任何细胞
## ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

## 1. 读取未过滤对象
input_file <- file.path(
  "02_objects",
  "human_raw_all_cells.rds"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "没有找到对象：",
      input_file
    )
  )
}

cat("正在读取未过滤对象……\n")

human_raw <- readRDS(input_file)

cat("细胞数：", ncol(human_raw), "\n")
cat("基因数：", nrow(human_raw), "\n")

## 2. 确认QC指标存在
required_metadata <- c(
  "sample_id",
  "donor_id",
  "condition",
  "nCount_RNA",
  "nFeature_RNA"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(human_raw@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    paste(
      "缺少metadata字段：",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

## 如果此前没有计算，则重新计算
if (!"percent.mt" %in% colnames(human_raw@meta.data)) {
  human_raw[["percent.mt"]] <- PercentageFeatureSet(
    human_raw,
    pattern = "^MT-"
  )
}

if (!"percent.ribo" %in% colnames(human_raw@meta.data)) {
  human_raw[["percent.ribo"]] <- PercentageFeatureSet(
    human_raw,
    pattern = "^(RPS|RPL)"
  )
}

## 3. 提取每个细胞的QC数据
qc_data <- human_raw@meta.data %>%
  tibble::rownames_to_column(
    var = "cell_id"
  ) %>%
  select(
    cell_id,
    sample_id,
    donor_id,
    condition,
    nCount_RNA,
    nFeature_RNA,
    percent.mt,
    percent.ribo
  )

## 4. 每个样本的QC汇总
qc_summary <- qc_data %>%
  group_by(
    sample_id,
    donor_id,
    condition
  ) %>%
  summarise(
    n_cells = n(),
    
    nFeature_min = min(nFeature_RNA),
    nFeature_q01 = quantile(nFeature_RNA, 0.01),
    nFeature_q05 = quantile(nFeature_RNA, 0.05),
    nFeature_median = median(nFeature_RNA),
    nFeature_q95 = quantile(nFeature_RNA, 0.95),
    nFeature_q99 = quantile(nFeature_RNA, 0.99),
    nFeature_max = max(nFeature_RNA),
    
    nCount_min = min(nCount_RNA),
    nCount_q01 = quantile(nCount_RNA, 0.01),
    nCount_q05 = quantile(nCount_RNA, 0.05),
    nCount_median = median(nCount_RNA),
    nCount_q95 = quantile(nCount_RNA, 0.95),
    nCount_q99 = quantile(nCount_RNA, 0.99),
    nCount_max = max(nCount_RNA),
    
    percent_mt_median = median(percent.mt),
    percent_mt_q95 = quantile(percent.mt, 0.95),
    percent_mt_q99 = quantile(percent.mt, 0.99),
    percent_mt_max = max(percent.mt),
    
    percent_ribo_median = median(percent.ribo),
    percent_ribo_q95 = quantile(percent.ribo, 0.95),
    percent_ribo_q99 = quantile(percent.ribo, 0.99),
    
    .groups = "drop"
  )

write.csv(
  qc_summary,
  file = file.path(
    "04_tables",
    "human_raw_qc_summary_by_sample.csv"
  ),
  row.names = FALSE
)

write.csv(
  qc_data,
  file = file.path(
    "04_tables",
    "human_raw_qc_per_cell.csv"
  ),
  row.names = FALSE
)

cat("\n========== 每个样本QC汇总 ==========\n")
print(
  qc_summary,
  row.names = FALSE
)

## 5. nFeature_RNA小提琴图
p_feature <- ggplot(
  qc_data,
  aes(
    x = sample_id,
    y = nFeature_RNA
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  labs(
    title = "Number of detected genes per nucleus",
    x = NULL,
    y = "nFeature_RNA"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

## 6. nCount_RNA小提琴图
p_count <- ggplot(
  qc_data,
  aes(
    x = sample_id,
    y = nCount_RNA
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  labs(
    title = "Total UMI counts per nucleus",
    x = NULL,
    y = "nCount_RNA"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

## 7. 线粒体比例小提琴图
p_mt <- ggplot(
  qc_data,
  aes(
    x = sample_id,
    y = percent.mt
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  labs(
    title = "Mitochondrial transcript percentage",
    x = NULL,
    y = "percent.mt"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

## 8. 核糖体比例小提琴图
p_ribo <- ggplot(
  qc_data,
  aes(
    x = sample_id,
    y = percent.ribo
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  labs(
    title = "Ribosomal transcript percentage",
    x = NULL,
    y = "percent.ribo"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

## 9. 合并QC图
qc_violin_plot <- (
  p_feature /
    p_count /
    p_mt /
    p_ribo
)

ggsave(
  filename = file.path(
    "03_figures",
    "human_raw_qc_violin_by_sample.pdf"
  ),
  plot = qc_violin_plot,
  width = 12,
  height = 16
)

ggsave(
  filename = file.path(
    "03_figures",
    "human_raw_qc_violin_by_sample.png"
  ),
  plot = qc_violin_plot,
  width = 12,
  height = 16,
  dpi = 300
)

## 10. nCount与nFeature关系
p_scatter <- ggplot(
  qc_data,
  aes(
    x = nCount_RNA,
    y = nFeature_RNA
  )
) +
  geom_point(
    alpha = 0.18,
    size = 0.25
  ) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(
    ~ sample_id,
    ncol = 4
  ) +
  labs(
    title = "Relationship between UMI counts and detected genes",
    x = "nCount_RNA, log10 scale",
    y = "nFeature_RNA, log10 scale"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    "03_figures",
    "human_raw_nCount_vs_nFeature_by_sample.pdf"
  ),
  plot = p_scatter,
  width = 14,
  height = 8
)

ggsave(
  filename = file.path(
    "03_figures",
    "human_raw_nCount_vs_nFeature_by_sample.png"
  ),
  plot = p_scatter,
  width = 14,
  height = 8,
  dpi = 300
)

## 11. 保存含QC指标、但仍未过滤的对象
qc_metric_object <- file.path(
  "02_objects",
  "human_raw_with_qc_metrics.rds"
)

saveRDS(
  human_raw,
  file = qc_metric_object,
  compress = FALSE
)

## 12. 保存运行环境
writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    "05_logs",
    "02_human_qc_overview_sessionInfo.txt"
  )
)

cat("\n========================================\n")
cat("原始QC汇总和绘图完成。\n")
cat("尚未过滤任何细胞。\n")
cat(
  "QC汇总表：",
  "04_tables/human_raw_qc_summary_by_sample.csv\n"
)
cat(
  "QC小提琴图：",
  "03_figures/human_raw_qc_violin_by_sample.png\n"
)
cat(
  "散点图：",
  "03_figures/human_raw_nCount_vs_nFeature_by_sample.png\n"
)
cat(
  "对象：",
  qc_metric_object,
  "\n"
)
cat("========================================\n")

View(qc_summary)
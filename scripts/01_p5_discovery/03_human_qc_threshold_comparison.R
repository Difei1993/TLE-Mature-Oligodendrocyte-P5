## ============================================================
## 03_human_qc_threshold_comparison.R
##
## 目的：
## 1. 按样本计算自适应QC阈值
## 2. 比较不同QC方案会保留/删除多少细胞
## 3. 绘制候选阈值和保留比例
##
## 本阶段只做预览，不过滤细胞
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
## 1. 读取带QC指标但尚未过滤的对象
## ------------------------------------------------------------

input_file <- file.path(
  "02_objects",
  "human_raw_with_qc_metrics.rds"
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

required_columns <- c(
  "sample_id",
  "donor_id",
  "condition",
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "percent.ribo"
)

missing_columns <- setdiff(
  required_columns,
  colnames(human_raw@meta.data)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "缺少以下metadata字段：",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}

qc_data <- human_raw@meta.data

qc_data$cell_id <- rownames(qc_data)

qc_data <- qc_data[
  ,
  c(
    "cell_id",
    "sample_id",
    "donor_id",
    "condition",
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "percent.ribo"
  )
]

## ------------------------------------------------------------
## 2. 定义稳健阈值函数
##
## nCount和nFeature均为右偏分布，因此在log10尺度上计算MAD。
## ------------------------------------------------------------

robust_log_limits <- function(
    x,
    nmads = 3
) {
  
  x <- as.numeric(x)
  
  log_x <- log10(x + 1)
  
  median_log <- median(
    log_x,
    na.rm = TRUE
  )
  
  mad_log <- mad(
    log_x,
    center = median_log,
    constant = 1.4826,
    na.rm = TRUE
  )
  
  ## 如果MAD为0，使用IQR估计替代
  if (
    !is.finite(mad_log) ||
    mad_log <= 0
  ) {
    
    mad_log <- IQR(
      log_x,
      na.rm = TRUE
    ) / 1.349
  }
  
  ## 如果仍然无法估计，退回1%和99%分位数
  if (
    !is.finite(mad_log) ||
    mad_log <= 0
  ) {
    
    lower <- quantile(
      x,
      0.01,
      na.rm = TRUE
    )
    
    upper <- quantile(
      x,
      0.99,
      na.rm = TRUE
    )
    
  } else {
    
    lower <- 10^(
      median_log - nmads * mad_log
    ) - 1
    
    upper <- 10^(
      median_log + nmads * mad_log
    ) - 1
  }
  
  c(
    lower = as.numeric(lower),
    upper = as.numeric(upper)
  )
}

## 线粒体比例阈值：
## 使用样本内 median + n×MAD，
## 但最低不低于5%，最高不超过设定的绝对上限。
robust_mt_limit <- function(
    x,
    nmads = 3,
    absolute_cap = 15,
    minimum_limit = 5
) {
  
  x <- as.numeric(x)
  
  median_x <- median(
    x,
    na.rm = TRUE
  )
  
  mad_x <- mad(
    x,
    center = median_x,
    constant = 1.4826,
    na.rm = TRUE
  )
  
  candidate <- median_x + nmads * mad_x
  
  if (
    !is.finite(candidate) ||
    candidate <= median_x
  ) {
    
    candidate <- quantile(
      x,
      0.99,
      na.rm = TRUE
    )
  }
  
  candidate <- max(
    candidate,
    minimum_limit
  )
  
  candidate <- min(
    candidate,
    absolute_cap
  )
  
  as.numeric(candidate)
}

## ------------------------------------------------------------
## 3. 为每个样本分别计算阈值
## ------------------------------------------------------------

sample_split <- split(
  qc_data,
  qc_data$sample_id
)

threshold_list <- lapply(
  sample_split,
  function(sample_data) {
    
    sample_id_value <- unique(
      sample_data$sample_id
    )
    
    donor_id_value <- unique(
      sample_data$donor_id
    )
    
    condition_value <- unique(
      sample_data$condition
    )
    
    feature_3mad <- robust_log_limits(
      sample_data$nFeature_RNA,
      nmads = 3
    )
    
    count_3mad <- robust_log_limits(
      sample_data$nCount_RNA,
      nmads = 3
    )
    
    feature_3.5mad <- robust_log_limits(
      sample_data$nFeature_RNA,
      nmads = 3.5
    )
    
    count_3.5mad <- robust_log_limits(
      sample_data$nCount_RNA,
      nmads = 3.5
    )
    
    data.frame(
      sample_id = sample_id_value,
      donor_id = donor_id_value,
      condition = condition_value,
      n_cells = nrow(sample_data),
      
      ## 3 MAD阈值
      feature_low_3mad = max(
        200,
        round(feature_3mad["lower"])
      ),
      
      feature_high_3mad = round(
        feature_3mad["upper"]
      ),
      
      count_low_3mad = max(
        500,
        round(count_3mad["lower"])
      ),
      
      count_high_3mad = round(
        count_3mad["upper"]
      ),
      
      mt_limit_3mad_cap15 = robust_mt_limit(
        sample_data$percent.mt,
        nmads = 3,
        absolute_cap = 15,
        minimum_limit = 5
      ),
      
      mt_limit_3mad_cap10 = robust_mt_limit(
        sample_data$percent.mt,
        nmads = 3,
        absolute_cap = 10,
        minimum_limit = 5
      ),
      
      ## 3.5 MAD宽松阈值
      feature_low_3.5mad = max(
        200,
        round(feature_3.5mad["lower"])
      ),
      
      feature_high_3.5mad = round(
        feature_3.5mad["upper"]
      ),
      
      count_low_3.5mad = max(
        500,
        round(count_3.5mad["lower"])
      ),
      
      count_high_3.5mad = round(
        count_3.5mad["upper"]
      ),
      
      mt_limit_3.5mad_cap15 = robust_mt_limit(
        sample_data$percent.mt,
        nmads = 3.5,
        absolute_cap = 15,
        minimum_limit = 5
      ),
      
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
)

thresholds <- bind_rows(
  threshold_list
)

thresholds <- thresholds[
  order(
    thresholds$condition,
    thresholds$sample_id
  ),
]

write.csv(
  thresholds,
  file = file.path(
    "04_tables",
    "human_qc_candidate_thresholds_by_sample.csv"
  ),
  row.names = FALSE
)

cat("\n========== 每个样本的候选阈值 ==========\n")
print(
  thresholds,
  row.names = FALSE
)

## ------------------------------------------------------------
## 4. 将阈值加入每个细胞的QC数据
## ------------------------------------------------------------

qc_flagged <- qc_data %>%
  left_join(
    thresholds,
    by = c(
      "sample_id",
      "donor_id",
      "condition"
    )
  )

## ------------------------------------------------------------
## 5. 建立4套候选QC规则
##
## 规则A：3 MAD + mt最高15%，相对平衡
## 规则B：3 MAD + mt最高10%，更严格
## 规则C：3.5 MAD + mt最高15%，较宽松
## 规则D：传统固定阈值，仅作为对照，不推荐直接采用
## ------------------------------------------------------------

qc_flagged$pass_adaptive_3mad_mt15 <- with(
  qc_flagged,
  
  nFeature_RNA >= feature_low_3mad &
    nFeature_RNA <= feature_high_3mad &
    nCount_RNA >= count_low_3mad &
    nCount_RNA <= count_high_3mad &
    percent.mt <= mt_limit_3mad_cap15
)

qc_flagged$pass_adaptive_3mad_mt10 <- with(
  qc_flagged,
  
  nFeature_RNA >= feature_low_3mad &
    nFeature_RNA <= feature_high_3mad &
    nCount_RNA >= count_low_3mad &
    nCount_RNA <= count_high_3mad &
    percent.mt <= mt_limit_3mad_cap10
)

qc_flagged$pass_adaptive_3.5mad_mt15 <- with(
  qc_flagged,
  
  nFeature_RNA >= `feature_low_3.5mad` &
    nFeature_RNA <= `feature_high_3.5mad` &
    nCount_RNA >= `count_low_3.5mad` &
    nCount_RNA <= `count_high_3.5mad` &
    percent.mt <= `mt_limit_3.5mad_cap15`
)

## 传统固定规则只用于比较
qc_flagged$pass_fixed_200_6000_mt5 <- with(
  qc_flagged,
  
  nFeature_RNA >= 200 &
    nFeature_RNA <= 6000 &
    percent.mt <= 5
)

## ------------------------------------------------------------
## 6. 生成每种规则的细胞保留情况
## ------------------------------------------------------------

summarise_rule <- function(
    data,
    pass_column,
    rule_name
) {
  
  data %>%
    group_by(
      sample_id,
      donor_id,
      condition
    ) %>%
    summarise(
      n_total = n(),
      
      n_keep = sum(
        .data[[pass_column]],
        na.rm = TRUE
      ),
      
      n_remove = n_total - n_keep,
      
      percent_keep = 100 * n_keep / n_total,
      
      percent_remove = 100 * n_remove / n_total,
      
      .groups = "drop"
    ) %>%
    mutate(
      rule = rule_name
    )
}

retention_summary <- bind_rows(
  summarise_rule(
    qc_flagged,
    "pass_adaptive_3mad_mt15",
    "Adaptive 3 MAD + mt cap 15%"
  ),
  
  summarise_rule(
    qc_flagged,
    "pass_adaptive_3mad_mt10",
    "Adaptive 3 MAD + mt cap 10%"
  ),
  
  summarise_rule(
    qc_flagged,
    "pass_adaptive_3.5mad_mt15",
    "Adaptive 3.5 MAD + mt cap 15%"
  ),
  
  summarise_rule(
    qc_flagged,
    "pass_fixed_200_6000_mt5",
    "Fixed 200–6000 genes + mt 5%"
  )
)

retention_summary <- retention_summary[
  order(
    retention_summary$rule,
    retention_summary$condition,
    retention_summary$sample_id
  ),
]

write.csv(
  retention_summary,
  file = file.path(
    "04_tables",
    "human_qc_candidate_retention_summary.csv"
  ),
  row.names = FALSE
)

cat("\n========== 各QC规则保留情况 ==========\n")
print(
  retention_summary,
  row.names = FALSE
)

## ------------------------------------------------------------
## 7. 查看平衡方案中每类删除原因
##
## 注意：同一个细胞可能同时满足多个删除原因，
## 因此各原因数量相加可能大于总删除数。
## ------------------------------------------------------------

qc_flagged$fail_low_feature <- with(
  qc_flagged,
  nFeature_RNA < feature_low_3mad
)

qc_flagged$fail_high_feature <- with(
  qc_flagged,
  nFeature_RNA > feature_high_3mad
)

qc_flagged$fail_low_count <- with(
  qc_flagged,
  nCount_RNA < count_low_3mad
)

qc_flagged$fail_high_count <- with(
  qc_flagged,
  nCount_RNA > count_high_3mad
)

qc_flagged$fail_high_mt <- with(
  qc_flagged,
  percent.mt > mt_limit_3mad_cap15
)

removal_reason_summary <- qc_flagged %>%
  group_by(
    sample_id,
    donor_id,
    condition
  ) %>%
  summarise(
    n_total = n(),
    
    n_fail_low_feature = sum(
      fail_low_feature,
      na.rm = TRUE
    ),
    
    n_fail_high_feature = sum(
      fail_high_feature,
      na.rm = TRUE
    ),
    
    n_fail_low_count = sum(
      fail_low_count,
      na.rm = TRUE
    ),
    
    n_fail_high_count = sum(
      fail_high_count,
      na.rm = TRUE
    ),
    
    n_fail_high_mt = sum(
      fail_high_mt,
      na.rm = TRUE
    ),
    
    n_fail_any = sum(
      !pass_adaptive_3mad_mt15,
      na.rm = TRUE
    ),
    
    percent_fail_any = 100 * n_fail_any / n_total,
    
    .groups = "drop"
  )

write.csv(
  removal_reason_summary,
  file = file.path(
    "04_tables",
    "human_qc_balanced_removal_reasons.csv"
  ),
  row.names = FALSE
)

## 保存逐细胞QC标记
write.csv(
  qc_flagged[
    ,
    c(
      "cell_id",
      "sample_id",
      "donor_id",
      "condition",
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      "percent.ribo",
      "pass_adaptive_3mad_mt15",
      "pass_adaptive_3mad_mt10",
      "pass_adaptive_3.5mad_mt15",
      "pass_fixed_200_6000_mt5",
      "fail_low_feature",
      "fail_high_feature",
      "fail_low_count",
      "fail_high_count",
      "fail_high_mt"
    )
  ],
  file = file.path(
    "04_tables",
    "human_qc_candidate_flags_per_cell.csv"
  ),
  row.names = FALSE
)

## ------------------------------------------------------------
## 8. 绘制不同规则的保留比例
## ------------------------------------------------------------

retention_summary$sample_id <- factor(
  retention_summary$sample_id,
  levels = unique(thresholds$sample_id)
)

p_retention <- ggplot(
  retention_summary,
  aes(
    x = sample_id,
    y = percent_keep,
    fill = condition
  )
) +
  geom_col(
    width = 0.75
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f%%",
        percent_keep
      )
    ),
    vjust = -0.25,
    size = 3
  ) +
  facet_wrap(
    ~ rule,
    ncol = 1
  ) +
  coord_cartesian(
    ylim = c(0, 105)
  ) +
  labs(
    title = "Cell retention under candidate QC rules",
    x = NULL,
    y = "Cells retained (%)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "top"
  )

ggsave(
  filename = file.path(
    "03_figures",
    "human_qc_candidate_retention_by_sample.pdf"
  ),
  plot = p_retention,
  width = 12,
  height = 14
)

ggsave(
  filename = file.path(
    "03_figures",
    "human_qc_candidate_retention_by_sample.png"
  ),
  plot = p_retention,
  width = 12,
  height = 14,
  dpi = 300
)

## ------------------------------------------------------------
## 9. 在nCount–nFeature散点图上标注平衡方案阈值
## ------------------------------------------------------------

p_scatter_threshold <- ggplot(
  qc_flagged,
  aes(
    x = nCount_RNA,
    y = nFeature_RNA
  )
) +
  geom_point(
    aes(
      color = pass_adaptive_3mad_mt15
    ),
    alpha = 0.18,
    size = 0.25
  ) +
  geom_vline(
    data = thresholds,
    aes(
      xintercept = count_low_3mad
    ),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  geom_vline(
    data = thresholds,
    aes(
      xintercept = count_high_3mad
    ),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  geom_hline(
    data = thresholds,
    aes(
      yintercept = feature_low_3mad
    ),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  geom_hline(
    data = thresholds,
    aes(
      yintercept = feature_high_3mad
    ),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(
    ~ sample_id,
    ncol = 4,
    scales = "free"
  ) +
  labs(
    title = "Adaptive 3-MAD QC boundaries",
    subtitle = "Dashed lines indicate candidate nCount and nFeature limits",
    x = "nCount_RNA, log10 scale",
    y = "nFeature_RNA, log10 scale",
    color = "Pass candidate QC"
  ) +
  theme_classic() +
  theme(
    legend.position = "top"
  )

ggsave(
  filename = file.path(
    "03_figures",
    "human_qc_adaptive_threshold_scatter.pdf"
  ),
  plot = p_scatter_threshold,
  width = 14,
  height = 8
)

ggsave(
  filename = file.path(
    "03_figures",
    "human_qc_adaptive_threshold_scatter.png"
  ),
  plot = p_scatter_threshold,
  width = 14,
  height = 8,
  dpi = 300
)

## ------------------------------------------------------------
## 10. 在线粒体比例小提琴图上标注候选阈值
## ------------------------------------------------------------

p_mt_threshold <- ggplot(
  qc_flagged,
  aes(
    x = sample_id,
    y = percent.mt,
    fill = condition
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
  geom_point(
    data = thresholds,
    aes(
      x = sample_id,
      y = mt_limit_3mad_cap15
    ),
    inherit.aes = FALSE,
    size = 2
  ) +
  labs(
    title = "Candidate mitochondrial thresholds",
    subtitle = "Points indicate adaptive thresholds with an absolute cap of 15%",
    x = NULL,
    y = "percent.mt"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "top"
  )

ggsave(
  filename = file.path(
    "03_figures",
    "human_qc_candidate_mt_thresholds.pdf"
  ),
  plot = p_mt_threshold,
  width = 12,
  height = 7
)

ggsave(
  filename = file.path(
    "03_figures",
    "human_qc_candidate_mt_thresholds.png"
  ),
  plot = p_mt_threshold,
  width = 12,
  height = 7,
  dpi = 300
)

## ------------------------------------------------------------
## 11. 保存运行环境
## ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    "05_logs",
    "03_human_qc_threshold_comparison_sessionInfo.txt"
  )
)

cat("\n========================================\n")
cat("候选QC阈值比较完成。\n")
cat("没有删除或过滤任何细胞。\n")
cat("\n主要结果：\n")
cat(
  "04_tables/human_qc_candidate_thresholds_by_sample.csv\n"
)
cat(
  "04_tables/human_qc_candidate_retention_summary.csv\n"
)
cat(
  "04_tables/human_qc_balanced_removal_reasons.csv\n"
)
cat(
  "03_figures/human_qc_candidate_retention_by_sample.png\n"
)
cat(
  "03_figures/human_qc_adaptive_threshold_scatter.png\n"
)
cat(
  "03_figures/human_qc_candidate_mt_thresholds.png\n"
)
cat("========================================\n")

View(thresholds)
View(retention_summary)
## ============================================================
## 04_human_stage1_qc_preview.R
##
## 第一阶段QC预览：
## 1. 样本自适应的低nFeature阈值
## 2. nCount_RNA至少500
## 3. percent.mt不超过5%
## 4. 高nFeature和高nCount仅标记，不立即删除
## ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

## 1. 读取未过滤对象
human_raw <- readRDS(
  "02_objects/human_raw_with_qc_metrics.rds"
)

## 2. 读取之前计算的样本阈值
thresholds <- read.csv(
  "04_tables/human_qc_candidate_thresholds_by_sample.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

## 3. 提取细胞metadata
qc_data <- human_raw@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
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

## 4. 加入样本阈值
qc_preview <- qc_data %>%
  left_join(
    thresholds %>%
      select(
        sample_id,
        feature_low_3mad,
        feature_high_3mad,
        count_high_3mad
      ),
    by = "sample_id"
  )

## 5. 定义第一阶段QC
qc_preview <- qc_preview %>%
  mutate(
    ## 明显低复杂度
    fail_low_feature =
      nFeature_RNA < feature_low_3mad,
    
    fail_low_count =
      nCount_RNA < 500,
    
    ## 高线粒体比例
    fail_high_mt =
      percent.mt > 5,
    
    ## 第一阶段是否保留
    pass_stage1_qc =
      !fail_low_feature &
      !fail_low_count &
      !fail_high_mt,
    
    ## 高复杂度只做标记，暂不删除
    high_feature_flag =
      nFeature_RNA > feature_high_3mad,
    
    high_count_flag =
      nCount_RNA > count_high_3mad
  )

## 6. 每个样本汇总
stage1_summary <- qc_preview %>%
  group_by(
    sample_id,
    donor_id,
    condition
  ) %>%
  summarise(
    n_total = n(),
    
    n_keep_stage1 = sum(
      pass_stage1_qc
    ),
    
    percent_keep_stage1 =
      100 * n_keep_stage1 / n_total,
    
    n_fail_low_feature = sum(
      fail_low_feature
    ),
    
    n_fail_low_count = sum(
      fail_low_count
    ),
    
    n_fail_high_mt = sum(
      fail_high_mt
    ),
    
    n_high_feature_flag = sum(
      high_feature_flag
    ),
    
    n_high_count_flag = sum(
      high_count_flag
    ),
    
    .groups = "drop"
  )

## 7. 保存结果
write.csv(
  qc_preview,
  "04_tables/human_stage1_qc_flags_per_cell.csv",
  row.names = FALSE
)

write.csv(
  stage1_summary,
  "04_tables/human_stage1_qc_preview_summary.csv",
  row.names = FALSE
)

cat("========== 第一阶段QC预览 ==========\n")

print(
  stage1_summary,
  row.names = FALSE
)

cat("\n注意：本脚本没有删除任何细胞。\n")
cat("高nFeature和高nCount细胞仅被标记，等待双细胞检测。\n")

View(stage1_summary)
## ============================================================
## 05_apply_human_stage1_qc.R
## 正式应用第一阶段QC：
## - 删除低nFeature
## - 删除低nCount
## - 删除percent.mt > 5%
## - 高nFeature和高nCount只保留标记，不删除
## ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

## 1. 读取未过滤对象和逐细胞QC标记
human_raw <- readRDS(
  "02_objects/human_raw_with_qc_metrics.rds"
)

qc_flags <- read.csv(
  "04_tables/human_stage1_qc_flags_per_cell.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "cell_id",
  "pass_stage1_qc",
  "fail_low_feature",
  "fail_low_count",
  "fail_high_mt",
  "high_feature_flag",
  "high_count_flag"
)

missing_columns <- setdiff(
  required_columns,
  colnames(qc_flags)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "QC标记表缺少字段：",
      paste(missing_columns, collapse = ", ")
    )
  )
}

## 2. 检查细胞名称是否一致
missing_in_flags <- setdiff(
  colnames(human_raw),
  qc_flags$cell_id
)

extra_in_flags <- setdiff(
  qc_flags$cell_id,
  colnames(human_raw)
)

if (length(missing_in_flags) > 0) {
  stop(
    paste0(
      "有",
      length(missing_in_flags),
      "个Seurat细胞在QC表中找不到。"
    )
  )
}

if (length(extra_in_flags) > 0) {
  stop(
    paste0(
      "QC表中有",
      length(extra_in_flags),
      "个细胞不在Seurat对象中。"
    )
  )
}

## 3. 按Seurat对象中的细胞顺序排列QC标记
qc_flags <- qc_flags[
  match(
    colnames(human_raw),
    qc_flags$cell_id
  ),
]

stopifnot(
  identical(
    qc_flags$cell_id,
    colnames(human_raw)
  )
)

## 4. 将QC标记加入Seurat metadata
qc_metadata <- qc_flags[
  ,
  c(
    "pass_stage1_qc",
    "fail_low_feature",
    "fail_low_count",
    "fail_high_mt",
    "high_feature_flag",
    "high_count_flag"
  )
]

rownames(qc_metadata) <- qc_flags$cell_id

human_raw <- AddMetaData(
  object = human_raw,
  metadata = qc_metadata
)

## 5. 保存被删除细胞名单
removed_cells <- human_raw@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  filter(!pass_stage1_qc) %>%
  select(
    cell_id,
    sample_id,
    donor_id,
    condition,
    nCount_RNA,
    nFeature_RNA,
    percent.mt,
    fail_low_feature,
    fail_low_count,
    fail_high_mt,
    high_feature_flag,
    high_count_flag
  )

write.csv(
  removed_cells,
  "04_tables/human_stage1_qc_removed_cells.csv",
  row.names = FALSE
)

## 6. 正式过滤
human_stage1 <- subset(
  human_raw,
  subset = pass_stage1_qc
)

## 7. 生成过滤前后汇总
before_summary <- human_raw@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  group_by(
    sample_id,
    donor_id,
    condition
  ) %>%
  summarise(
    n_before = n(),
    .groups = "drop"
  )

after_summary <- human_stage1@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  group_by(
    sample_id,
    donor_id,
    condition
  ) %>%
  summarise(
    n_after = n(),
    .groups = "drop"
  )

stage1_summary <- before_summary %>%
  left_join(
    after_summary,
    by = c(
      "sample_id",
      "donor_id",
      "condition"
    )
  ) %>%
  mutate(
    n_removed = n_before - n_after,
    percent_retained = 100 * n_after / n_before
  )

write.csv(
  stage1_summary,
  "04_tables/human_stage1_qc_final_summary.csv",
  row.names = FALSE
)

## 8. 基本验收
cat("========== 第一阶段QC结果 ==========\n")
cat("过滤前细胞数：", ncol(human_raw), "\n")
cat("过滤后细胞数：", ncol(human_stage1), "\n")
cat("删除细胞数：", ncol(human_raw) - ncol(human_stage1), "\n")
cat(
  "总保留率：",
  round(
    100 * ncol(human_stage1) / ncol(human_raw),
    2
  ),
  "%\n"
)

print(
  stage1_summary,
  row.names = FALSE
)

if (ncol(human_stage1) != 64936) {
  warning(
    paste0(
      "预期保留64936个细胞，实际为",
      ncol(human_stage1),
      "。请检查。"
    )
  )
}

## 9. 保存对象
saveRDS(
  human_stage1,
  file = "02_objects/human_stage1_qc_filtered.rds",
  compress = FALSE
)

## 10. 保存sessionInfo
writeLines(
  capture.output(sessionInfo()),
  "05_logs/05_apply_human_stage1_qc_sessionInfo.txt"
)

cat("\n对象已保存：\n")
cat("02_objects/human_stage1_qc_filtered.rds\n")
cat("\n下一步：按每个样本进行双细胞检测。\n")

View(stage1_summary)
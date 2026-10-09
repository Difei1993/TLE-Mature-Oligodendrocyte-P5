## ============================================================
## 09_apply_human_doublet_filter.R
##
## 正式移除scDblFinder预测的doublet
## 保留singlet进入后续标准化、降维和聚类
## ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

## 1. 读取带双细胞评分的对象
input_file <- file.path(
  "02_objects",
  "human_stage1_with_doublet_scores.rds"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "找不到输入对象：",
      input_file
    )
  )
}

human_scored <- readRDS(input_file)

required_columns <- c(
  "sample_id",
  "donor_id",
  "condition",
  "scDblFinder.score",
  "scDblFinder.class"
)

missing_columns <- setdiff(
  required_columns,
  colnames(human_scored@meta.data)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "对象缺少字段：",
      paste(missing_columns, collapse = ", ")
    )
  )
}

## 2. 保存被移除双细胞的信息
doublet_cells <- human_scored@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  filter(scDblFinder.class == "doublet") %>%
  select(
    cell_id,
    sample_id,
    donor_id,
    condition,
    nCount_RNA,
    nFeature_RNA,
    percent.mt,
    scDblFinder.score,
    scDblFinder.class,
    high_feature_flag,
    high_count_flag
  )

write.csv(
  doublet_cells,
  "04_tables/human_removed_predicted_doublets.csv",
  row.names = FALSE
)

## 3. 获取singlet细胞名称
singlet_cells <- rownames(
  human_scored@meta.data
)[
  human_scored$scDblFinder.class == "singlet"
]

if (length(singlet_cells) != 58341) {
  warning(
    paste0(
      "预期singlet为58341个，实际为",
      length(singlet_cells),
      "个。"
    )
  )
}

## 4. 生成主分析singlet对象
human_singlets <- subset(
  human_scored,
  cells = singlet_cells
)

## 5. 按样本汇总过滤前后数量
before_summary <- human_scored@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  group_by(
    sample_id,
    donor_id,
    condition
  ) %>%
  summarise(
    n_before_doublet_filter = n(),
    .groups = "drop"
  )

after_summary <- human_singlets@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  group_by(
    sample_id,
    donor_id,
    condition
  ) %>%
  summarise(
    n_after_doublet_filter = n(),
    .groups = "drop"
  )

doublet_filter_summary <- before_summary %>%
  left_join(
    after_summary,
    by = c(
      "sample_id",
      "donor_id",
      "condition"
    )
  ) %>%
  mutate(
    n_removed_doublets =
      n_before_doublet_filter -
      n_after_doublet_filter,
    
    doublet_rate_percent =
      100 *
      n_removed_doublets /
      n_before_doublet_filter,
    
    percent_retained =
      100 *
      n_after_doublet_filter /
      n_before_doublet_filter
  )

write.csv(
  doublet_filter_summary,
  "04_tables/human_doublet_filter_final_summary.csv",
  row.names = FALSE
)

## 6. 按条件汇总，仅用于检查过滤影响
condition_summary <- doublet_filter_summary %>%
  group_by(condition) %>%
  summarise(
    n_before = sum(n_before_doublet_filter),
    n_after = sum(n_after_doublet_filter),
    n_removed = sum(n_removed_doublets),
    doublet_rate_percent =
      100 * n_removed / n_before,
    .groups = "drop"
  )

write.csv(
  condition_summary,
  "04_tables/human_doublet_filter_by_condition.csv",
  row.names = FALSE
)

## 7. 基本验收
cat("========== 双细胞过滤结果 ==========\n")
cat("过滤前细胞数：", ncol(human_scored), "\n")
cat("移除预测双细胞：", nrow(doublet_cells), "\n")
cat("过滤后singlet数：", ncol(human_singlets), "\n")
cat(
  "总体保留率：",
  round(
    100 *
      ncol(human_singlets) /
      ncol(human_scored),
    2
  ),
  "%\n\n"
)

cat("========== 每个样本 ==========\n")
print(
  doublet_filter_summary,
  row.names = FALSE
)

cat("\n========== 按组汇总 ==========\n")
print(
  condition_summary,
  row.names = FALSE
)

## 8. 保存主分析对象
output_file <- file.path(
  "02_objects",
  "human_qc_singlets.rds"
)

saveRDS(
  human_singlets,
  file = output_file,
  compress = FALSE
)

## 保留细胞名单
write.csv(
  data.frame(
    cell_id = colnames(human_singlets),
    stringsAsFactors = FALSE
  ),
  "04_tables/human_final_qc_singlet_cell_ids.csv",
  row.names = FALSE
)

## 9. 保存运行环境
writeLines(
  capture.output(sessionInfo()),
  "05_logs/09_apply_human_doublet_filter_sessionInfo.txt"
)

cat("\n对象已保存：\n")
cat(output_file, "\n")
cat("\n下一步：检查过滤后QC分布，再进行标准化。\n")

View(doublet_filter_summary)
View(condition_summary)
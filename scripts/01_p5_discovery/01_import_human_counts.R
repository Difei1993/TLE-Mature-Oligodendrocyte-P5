## ============================================================
## 01_import_human_counts.R
## 导入GSE190452的8个原始计数矩阵
## 创建未经过QC过滤的人类Seurat对象
## ============================================================

options(stringsAsFactors = FALSE)
set.seed(20260724)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

## 1. 数据和输出路径
raw_dir <- file.path(
  "00_raw_data",
  "GSE190452_human"
)

object_dir <- "02_objects"
table_dir <- "04_tables"
log_dir <- "05_logs"

stopifnot(dir.exists(raw_dir))

## 2. 寻找全部barcode文件，并由文件名识别样本
barcode_files <- list.files(
  path = raw_dir,
  pattern = "_barcodes\\.tsv\\.gz$",
  full.names = TRUE
)

sample_ids <- sub(
  pattern = "_barcodes\\.tsv\\.gz$",
  replacement = "",
  x = basename(barcode_files)
)

sample_ids <- sort(sample_ids)

cat("检测到的样本数量：", length(sample_ids), "\n")
print(sample_ids)

if (length(sample_ids) != 8) {
  stop(
    paste0(
      "预期应有8个样本，但目前检测到",
      length(sample_ids),
      "个。请检查文件是否完整。"
    )
  )
}

## 3. 建立样本信息表
sample_table <- data.frame(
  sample_id = sample_ids,
  gsm = sub(
    pattern = "_.*$",
    replacement = "",
    x = sample_ids
  ),
  condition = ifelse(
    grepl("_TLE", sample_ids, ignore.case = TRUE),
    "TLE",
    "Normal"
  ),
  stringsAsFactors = FALSE
)

## 目前每个文库对应一名供者
## 后续会根据GEO临床信息再次核实
sample_table$library_id <- sample_table$sample_id
sample_table$donor_id <- sample_table$sample_id

cat("\n========== 样本信息 ==========\n")
print(sample_table, row.names = FALSE)

write.csv(
  sample_table,
  file = file.path(
    table_dir,
    "human_sample_manifest.csv"
  ),
  row.names = FALSE
)

## 4. 逐样本读取计数矩阵并创建Seurat对象
sample_objects <- vector(
  mode = "list",
  length = nrow(sample_table)
)

names(sample_objects) <- sample_table$sample_id

import_summary <- vector(
  mode = "list",
  length = nrow(sample_table)
)

for (i in seq_len(nrow(sample_table))) {
  
  sample_id <- sample_table$sample_id[i]
  
  cat(
    "\n========================================\n",
    "正在导入样本：", sample_id, "\n",
    sep = ""
  )
  
  matrix_file <- file.path(
    raw_dir,
    paste0(sample_id, "_matrix.mtx.gz")
  )
  
  feature_file <- file.path(
    raw_dir,
    paste0(sample_id, "_features.tsv.gz")
  )
  
  barcode_file <- file.path(
    raw_dir,
    paste0(sample_id, "_barcodes.tsv.gz")
  )
  
  required_files <- c(
    matrix_file,
    feature_file,
    barcode_file
  )
  
  if (!all(file.exists(required_files))) {
    missing_files <- required_files[
      !file.exists(required_files)
    ]
    
    stop(
      paste(
        "样本缺少以下文件：",
        paste(missing_files, collapse = "\n")
      )
    )
  }
  
  ## ReadMtx直接读取压缩的MTX、feature和barcode文件
  counts_matrix <- ReadMtx(
    mtx = matrix_file,
    cells = barcode_file,
    features = feature_file,
    cell.column = 1,
    feature.column = 2,
    unique.features = TRUE
  )
  
  cat(
    "基因数：", nrow(counts_matrix),
    "；细胞数：", ncol(counts_matrix),
    "\n"
  )
  
  ## 导入阶段不进行任何细胞或基因过滤
  sample_object <- CreateSeuratObject(
    counts = counts_matrix,
    project = sample_id,
    assay = "RNA",
    min.cells = 0,
    min.features = 0
  )
  
  ## 添加样本级metadata
  sample_object$sample_id <- sample_id
  sample_object$library_id <- sample_table$library_id[i]
  sample_object$donor_id <- sample_table$donor_id[i]
  sample_object$gsm <- sample_table$gsm[i]
  sample_object$condition <- sample_table$condition[i]
  
  ## 计算QC指标，但暂不据此过滤
  sample_object[["percent.mt"]] <- PercentageFeatureSet(
    sample_object,
    pattern = "^MT-"
  )
  
  sample_object[["percent.ribo"]] <- PercentageFeatureSet(
    sample_object,
    pattern = "^(RPS|RPL)"
  )
  
  sample_objects[[sample_id]] <- sample_object
  
  import_summary[[i]] <- data.frame(
    sample_id = sample_id,
    gsm = sample_table$gsm[i],
    condition = sample_table$condition[i],
    n_genes_matrix = nrow(counts_matrix),
    n_cells_matrix = ncol(counts_matrix),
    n_genes_seurat = nrow(sample_object),
    n_cells_seurat = ncol(sample_object),
    median_nCount_RNA = median(
      sample_object$nCount_RNA
    ),
    median_nFeature_RNA = median(
      sample_object$nFeature_RNA
    ),
    median_percent_mt = median(
      sample_object$percent.mt
    ),
    stringsAsFactors = FALSE
  )
  
  rm(
    counts_matrix,
    sample_object
  )
  
  invisible(gc())
}

import_summary <- do.call(
  rbind,
  import_summary
)

write.csv(
  import_summary,
  file = file.path(
    table_dir,
    "human_import_summary_by_sample.csv"
  ),
  row.names = FALSE
)

cat("\n========== 每个样本导入结果 ==========\n")
print(import_summary, row.names = FALSE)

## 5. 合并8个样本
cat("\n正在合并8个Seurat对象……\n")

human_raw <- merge(
  x = sample_objects[[1]],
  y = sample_objects[-1],
  add.cell.ids = names(sample_objects),
  merge.data = FALSE
)

DefaultAssay(human_raw) <- "RNA"

## 6. 检查合并结果
cat("\n========== 合并后对象 ==========\n")
cat("总基因数：", nrow(human_raw), "\n")
cat("总细胞数：", ncol(human_raw), "\n")

cat("\n各组细胞数：\n")
print(
  table(human_raw$condition)
)

cat("\n各样本细胞数：\n")
print(
  table(human_raw$sample_id)
)

## 7. 检查细胞名称是否唯一
if (anyDuplicated(colnames(human_raw)) > 0) {
  stop("合并后出现重复细胞名称。")
}

## 8. 保存完全未过滤的对象
raw_object_file <- file.path(
  object_dir,
  "human_raw_all_cells.rds"
)

cat(
  "\n正在保存未过滤对象：\n",
  raw_object_file,
  "\n"
)

saveRDS(
  human_raw,
  file = raw_object_file,
  compress = FALSE
)

## 9. 保存合并后的样本统计
merged_summary <- as.data.frame(
  table(
    sample_id = human_raw$sample_id,
    condition = human_raw$condition
  ),
  stringsAsFactors = FALSE
)

merged_summary <- merged_summary[
  merged_summary$Freq > 0,
]

colnames(merged_summary)[
  colnames(merged_summary) == "Freq"
] <- "n_cells"

write.csv(
  merged_summary,
  file = file.path(
    table_dir,
    "human_merged_cell_counts.csv"
  ),
  row.names = FALSE
)

## 10. 保存sessionInfo
session_text <- capture.output(
  sessionInfo()
)

writeLines(
  session_text,
  con = file.path(
    log_dir,
    "01_import_human_counts_sessionInfo.txt"
  )
)

cat("\n========================================\n")
cat("人类原始计数导入完成。\n")
cat("总细胞数：", ncol(human_raw), "\n")
cat("总基因数：", nrow(human_raw), "\n")
cat("对象保存位置：", raw_object_file, "\n")
cat("本阶段没有执行QC过滤。\n")
cat("========================================\n")
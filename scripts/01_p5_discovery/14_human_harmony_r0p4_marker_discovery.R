## ============================================================
## 14_human_harmony_r0p4_marker_discovery.R
##
## 目的：
## 1. 使用Harmony resolution 0.4的19个cluster
## 2. 在RNA标准化表达层面寻找cluster marker
## 3. 绘制经典脑细胞类型marker DotPlot和FeaturePlot
## 4. 检查Harmony cluster与未整合cluster的对应关系
##
## 注意：
## - Harmony只定义低维空间和cluster
## - marker来自RNA assay，不使用Harmony值
## - marker的细胞级P值只用于注释，不用于TLE疾病推断
## ============================================================

options(stringsAsFactors = FALSE)
set.seed(20260724)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

## ------------------------------------------------------------
## 1. 读取Harmony对象
## ------------------------------------------------------------

input_file <- file.path(
  "02_objects",
  "human_harmony_umap_multiresolution.rds"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "找不到输入对象：",
      input_file
    )
  )
}

cat("正在读取Harmony对象……\n")

human <- readRDS(input_file)

DefaultAssay(human) <- "RNA"

cluster_column <- "cluster_harmony30_r0p4"

if (
  !cluster_column %in%
  colnames(human@meta.data)
) {
  stop(
    paste0(
      "对象中找不到聚类字段：",
      cluster_column
    )
  )
}

cat("细胞数：", ncol(human), "\n")
cat("基因数：", nrow(human), "\n")

## ------------------------------------------------------------
## 2. 建立暂定主cluster字段
## ------------------------------------------------------------

cluster_values <- as.character(
  human@meta.data[[cluster_column]]
)

cluster_levels <- as.character(
  sort(
    unique(
      as.integer(cluster_values)
    )
  )
)

human$cluster_harmony_r0p4_provisional <- factor(
  cluster_values,
  levels = cluster_levels
)

Idents(human) <- "cluster_harmony_r0p4_provisional"

cat("\n========== 暂定cluster大小 ==========\n")

cluster_sizes <- data.frame(
  cluster = names(table(Idents(human))),
  n_cells = as.integer(table(Idents(human))),
  stringsAsFactors = FALSE
)

print(
  cluster_sizes,
  row.names = FALSE
)

write.csv(
  cluster_sizes,
  "04_tables/human_harmony_r0p4_cluster_sizes.csv",
  row.names = FALSE
)

## ------------------------------------------------------------
## 3. 合并Seurat v5的RNA layers
## ------------------------------------------------------------

cat("\n正在合并RNA assay layers……\n")

human <- JoinLayers(
  object = human,
  assay = "RNA"
)

## 检查标准化data layer
rna_layers <- Layers(
  human[["RNA"]]
)

cat("RNA layers：\n")
print(rna_layers)

if (!"data" %in% rna_layers) {
  stop(
    "RNA assay中没有data layer，无法进行marker分析。"
  )
}

## ------------------------------------------------------------
## 4. cluster × donor与cluster × condition构成
##
## 这些是细胞数量描述，不是疾病统计检验。
## ------------------------------------------------------------

metadata_for_summary <- human@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(
    cluster =
      as.character(
        cluster_harmony_r0p4_provisional
      )
  )

cluster_donor_counts <- metadata_for_summary %>%
  count(
    cluster,
    donor_id,
    sample_id,
    condition,
    name = "n_cells"
  )

write.csv(
  cluster_donor_counts,
  "04_tables/human_harmony_r0p4_cluster_donor_counts.csv",
  row.names = FALSE
)

cluster_summary <- metadata_for_summary %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    
    n_donors =
      n_distinct(donor_id),
    
    n_normal_donors =
      n_distinct(
        donor_id[condition == "Normal"]
      ),
    
    n_tle_donors =
      n_distinct(
        donor_id[condition == "TLE"]
      ),
    
    max_donor_fraction =
      max(table(donor_id)) / n(),
    
    normal_cell_fraction =
      mean(condition == "Normal"),
    
    tle_cell_fraction =
      mean(condition == "TLE"),
    
    .groups = "drop"
  )

write.csv(
  cluster_summary,
  "04_tables/human_harmony_r0p4_cluster_summary.csv",
  row.names = FALSE
)

cat("\n========== cluster供者构成 ==========\n")

print(
  cluster_summary,
  row.names = FALSE
)

## ------------------------------------------------------------
## 5. 寻找全部阳性marker
##
## max.cells.per.ident = 5000用于加速初筛。
## 最终注释后可对关键cluster用全部细胞复核。
## ------------------------------------------------------------

marker_output_file <- file.path(
  "04_tables",
  "human_harmony_r0p4_all_positive_markers.csv"
)

if (file.exists(marker_output_file)) {
  
  cat(
    "\n已找到现有marker结果，直接读取，不重复计算。\n"
  )
  
  all_markers <- read.csv(
    marker_output_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
} else {
  
  cat(
    "\n开始寻找19个cluster的阳性marker……\n"
  )
  
  all_markers <- FindAllMarkers(
    object = human,
    assay = "RNA",
    slot = "data",
    only.pos = TRUE,
    test.use = "wilcox",
    min.pct = 0.20,
    logfc.threshold = 0.25,
    max.cells.per.ident = 5000,
    random.seed = 20260724,
    verbose = TRUE
  )
  
  write.csv(
    all_markers,
    marker_output_file,
    row.names = FALSE
  )
}

if (nrow(all_markers) == 0) {
  stop("没有获得任何cluster marker。")
}

cat(
  "marker结果总行数：",
  nrow(all_markers),
  "\n"
)

## ------------------------------------------------------------
## 6. 每个cluster的top20和top50 marker
## ------------------------------------------------------------

fc_column <- intersect(
  c(
    "avg_log2FC",
    "avg_logFC"
  ),
  colnames(all_markers)
)

if (length(fc_column) == 0) {
  stop(
    "marker结果中找不到avg_log2FC或avg_logFC。"
  )
}

fc_column <- fc_column[1]

top20_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = .data[[fc_column]],
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

top50_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = .data[[fc_column]],
    n = 50,
    with_ties = FALSE
  ) %>%
  ungroup()

write.csv(
  top20_markers,
  "04_tables/human_harmony_r0p4_top20_markers.csv",
  row.names = FALSE
)

write.csv(
  top50_markers,
  "04_tables/human_harmony_r0p4_top50_markers.csv",
  row.names = FALSE
)

cat("\n========== 每个cluster前5个marker ==========\n")

top5_print <- all_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = .data[[fc_column]],
    n = 5,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    cluster,
    gene,
    all_of(fc_column),
    pct.1,
    pct.2,
    p_val_adj
  )

print(
  top5_print,
  n = Inf
)

write.csv(
  top5_print,
  "04_tables/human_harmony_r0p4_top5_markers_for_review.csv",
  row.names = FALSE
)

## ------------------------------------------------------------
## 7. 经典脑细胞marker面板
## ------------------------------------------------------------

marker_panels <- list(
  
  Neuronal = c(
    ## 兴奋性神经元
    "SLC17A7",
    "SATB2",
    "CUX1",
    "CUX2",
    "TBR1",
    "BCL11B",
    
    ## 抑制性神经元
    "GAD1",
    "GAD2",
    "SLC6A1",
    "PVALB",
    "SST",
    "VIP",
    "RELN",
    "LAMP5"
  ),
  
  Glial = c(
    ## 星形胶质
    "AQP4",
    "ALDH1L1",
    "SLC1A2",
    "SLC1A3",
    "GFAP",
    "GJA1",
    
    ## 少突胶质
    "PLP1",
    "MBP",
    "MOG",
    "MAG",
    "MOBP",
    "SOX10",
    
    ## OPC
    "PDGFRA",
    "CSPG4",
    "VCAN",
    "OLIG1",
    "OLIG2"
  ),
  
  Myeloid_Vascular_Other = c(
    ## 小胶质/髓系
    "P2RY12",
    "TMEM119",
    "CX3CR1",
    "AIF1",
    "C1QA",
    "C1QB",
    "C1QC",
    "INPP5D",
    "TYROBP",
    "FCER1G",
    "CD68",
    "APOE",
    "LPL",
    
    ## 内皮
    "CLDN5",
    "PECAM1",
    "VWF",
    "EMCN",
    "KDR",
    "FLT1",
    
    ## 周细胞
    "PDGFRB",
    "RGS5",
    "MCAM",
    "ABCC9",
    "KCNJ8",
    
    ## 血管脑膜细胞
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    
    ## 室管膜
    "FOXJ1",
    "PIFO",
    "TPPP3",
    "CAPS",
    
    ## 增殖/应激/免疫浸润
    "MKI67",
    "TOP2A",
    "FOS",
    "JUN",
    "HSPA1A",
    "CD3D",
    "NKG7",
    "MS4A1"
  )
)

## 只保留对象中实际存在的基因
marker_panels_present <- lapply(
  marker_panels,
  function(x) {
    intersect(
      x,
      rownames(human)
    )
  }
)

marker_presence <- bind_rows(
  lapply(
    names(marker_panels),
    function(panel_name) {
      
      data.frame(
        panel = panel_name,
        gene = marker_panels[[panel_name]],
        present =
          marker_panels[[panel_name]] %in%
          rownames(human),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  marker_presence,
  "04_tables/human_annotation_marker_presence.csv",
  row.names = FALSE
)

## ------------------------------------------------------------
## 8. 绘制三个marker DotPlot
## ------------------------------------------------------------

make_dotplot <- function(
    features,
    title_text
) {
  
  DotPlot(
    object = human,
    features = features,
    assay = "RNA",
    group.by =
      "cluster_harmony_r0p4_provisional",
    dot.scale = 6
  ) +
    RotatedAxis() +
    labs(
      title = title_text,
      x = NULL,
      y = "Harmony r0.4 cluster"
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      )
    )
}

p_dot_neuronal <- make_dotplot(
  marker_panels_present$Neuronal,
  "Canonical neuronal markers"
)

p_dot_glial <- make_dotplot(
  marker_panels_present$Glial,
  "Canonical glial markers"
)

p_dot_other <- make_dotplot(
  marker_panels_present$Myeloid_Vascular_Other,
  "Myeloid, vascular and quality-control markers"
)

ggsave(
  "03_figures/human_harmony_r0p4_dotplot_neuronal.png",
  p_dot_neuronal,
  width = 14,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  "03_figures/human_harmony_r0p4_dotplot_glial.png",
  p_dot_glial,
  width = 13,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  "03_figures/human_harmony_r0p4_dotplot_myeloid_vascular_other.png",
  p_dot_other,
  width = 19,
  height = 8,
  dpi = 300,
  bg = "white"
)

## PDF版本
ggsave(
  "03_figures/human_harmony_r0p4_dotplot_neuronal.pdf",
  p_dot_neuronal,
  width = 14,
  height = 8
)

ggsave(
  "03_figures/human_harmony_r0p4_dotplot_glial.pdf",
  p_dot_glial,
  width = 13,
  height = 8
)

ggsave(
  "03_figures/human_harmony_r0p4_dotplot_myeloid_vascular_other.pdf",
  p_dot_other,
  width = 19,
  height = 8
)

## ------------------------------------------------------------
## 9. 核心marker FeaturePlot
## ------------------------------------------------------------

core_feature_genes <- c(
  "SLC17A7",
  "GAD1",
  "AQP4",
  "PLP1",
  "PDGFRA",
  "P2RY12",
  "CLDN5",
  "RGS5",
  "FOXJ1",
  "MKI67",
  "FOS",
  "CD3D"
)

core_feature_genes <- intersect(
  core_feature_genes,
  rownames(human)
)

p_core_features <- FeaturePlot(
  object = human,
  features = core_feature_genes,
  reduction = "umap.harmony",
  order = TRUE,
  raster = TRUE,
  ncol = 4
)

ggsave(
  "03_figures/human_harmony_r0p4_core_marker_featureplots.png",
  p_core_features,
  width = 16,
  height = 13,
  dpi = 300,
  bg = "white"
)

ggsave(
  "03_figures/human_harmony_r0p4_core_marker_featureplots.pdf",
  p_core_features,
  width = 16,
  height = 13
)

## ------------------------------------------------------------
## 10. Harmony r0.4与未整合r0.2对应关系
## ------------------------------------------------------------

unintegrated_column <- "cluster_pca30_r0p2"
harmony_column <- "cluster_harmony_r0p4_provisional"

if (
  unintegrated_column %in% colnames(human@meta.data) &&
  harmony_column %in% colnames(human@meta.data)
) {
  
  cluster_cross_tab <- table(
    unintegrated = human@meta.data[[unintegrated_column]],
    harmony = human@meta.data[[harmony_column]]
  )
  
  write.csv(
    as.data.frame.matrix(cluster_cross_tab),
    "04_tables/human_unintegrated_r0p2_vs_harmony_r0p4_counts.csv",
    row.names = TRUE
  )
  
  ## 每个未整合cluster内部，各Harmony cluster所占比例
  row_fraction <- prop.table(
    cluster_cross_tab,
    margin = 1
  )
  
  write.csv(
    as.data.frame.matrix(row_fraction),
    "04_tables/human_unintegrated_r0p2_vs_harmony_r0p4_row_fraction.csv",
    row.names = TRUE
  )
  
  cat(
    "\n已保存未整合r0.2与Harmony r0.4的对应表。\n"
  )
  
} else {
  
  missing_cluster_columns <- setdiff(
    c(unintegrated_column, harmony_column),
    colnames(human@meta.data)
  )
  
  warning(
    paste(
      "无法生成新旧cluster对应表，缺少字段：",
      paste(missing_cluster_columns, collapse = ", ")
    )
  )
}
## ------------------------------------------------------------
## 11. 保存marker准备对象
## ------------------------------------------------------------

output_file <- file.path(
  "02_objects",
  "human_harmony_r0p4_marker_ready.rds"
)

saveRDS(
  human,
  output_file,
  compress = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  "05_logs/14_human_harmony_r0p4_marker_discovery_sessionInfo.txt"
)

cat("\n========================================\n")
cat("Harmony r0.4 marker初筛完成。\n")
cat("暂定cluster数量：", length(levels(Idents(human))), "\n")
cat("marker总行数：", nrow(all_markers), "\n")
cat("对象保存于：\n")
cat(output_file, "\n")
cat("\n主要表格：\n")
cat("04_tables/human_harmony_r0p4_top5_markers_for_review.csv\n")
cat("04_tables/human_harmony_r0p4_top20_markers.csv\n")
cat("\n主要图片：\n")
cat("03_figures/human_harmony_r0p4_dotplot_neuronal.png\n")
cat("03_figures/human_harmony_r0p4_dotplot_glial.png\n")
cat(
  "03_figures/",
  "human_harmony_r0p4_dotplot_myeloid_vascular_other.png\n",
  sep = ""
)
cat("03_figures/human_harmony_r0p4_core_marker_featureplots.png\n")
cat("========================================\n")

View(top5_print)
View(cluster_summary)
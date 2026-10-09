## ============================================================
## 15_human_assign_level1_annotation.R
##
## 目的：
## 1. 为Harmony resolution 0.4的19个cluster写入Level 1注释
## 2. 保存注释证据和供者构成
## 3. 绘制Level 1 UMAP和marker DotPlot
##
## 本阶段只进行宽细胞大类注释
## 不进行细胞亚型、状态或衰老命名
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
## 1. 读取marker分析后的对象
## ------------------------------------------------------------

input_file <- file.path(
  "02_objects",
  "human_harmony_r0p4_marker_ready.rds"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "找不到输入对象：",
      input_file
    )
  )
}

human <- readRDS(input_file)

DefaultAssay(human) <- "RNA"

cluster_column <- "cluster_harmony_r0p4_provisional"

if (!cluster_column %in% colnames(human@meta.data)) {
  stop(
    paste0(
      "对象中没有聚类字段：",
      cluster_column
    )
  )
}

## ------------------------------------------------------------
## 2. 建立注释证据表
## ------------------------------------------------------------

annotation_table <- data.frame(
  cluster = as.character(0:18),
  
  annotation_level1 = c(
    "Oligodendrocytes",       # 0
    "Excitatory neurons",     # 1
    "Microglia/macrophage",   # 2
    "Astrocytes",             # 3
    "Excitatory neurons",     # 4
    "OPC",                    # 5
    "Inhibitory neurons",     # 6
    "Inhibitory neurons",     # 7
    "Excitatory neurons",     # 8
    "Excitatory neurons",     # 9
    "Inhibitory neurons",     # 10
    "Astrocytes",             # 11
    "Inhibitory neurons",     # 12
    "Astrocytes",             # 13
    "Inhibitory neurons",     # 14
    "Excitatory neurons",     # 15
    "Excitatory neurons",     # 16
    "Vascular cells",         # 17
    "T lymphocytes"           # 18
  ),
  
  confidence = c(
    "High",
    "High",
    "High",
    "High",
    "High",
    "High",
    "High",
    "High",
    "Moderate",
    "Moderate",
    "High",
    "High",
    "High",
    "Moderate",
    "High",
    "Moderate",
    "Moderate",
    "Moderate",
    "High"
  ),
  
  marker_evidence = c(
    "MOBP; MOG; MBP; OPALIN; CLDN11; MYRF",
    "SLC17A7; SATB2; CUX1; CUX2",
    "CX3CR1; CSF1R; INPP5D; ITGAM; MS4A7; DOCK8",
    "SLC1A2; SLC4A4; ADGRV1; ALDH1A1",
    "SLC17A7; NEFL; NEFM; NPTX1",
    "PDGFRA; VCAN; CSPG4; GPR17",
    "VIP; TAC3; DLX6-AS1; SYNPR",
    "LHX6; SST; GRIK1; ARX",
    "SLC17A7; SATB2; CUX1-associated neuronal profile",
    "SLC17A7-associated neuronal profile",
    "GAD2; KCNS3; TAC1; LHX6",
    "ALDH1L1; SLC1A2; GFAP; ATP1A2",
    "NDNF; RELN; CALB2; inhibitory neuronal profile",
    "GFAP; AQP4; SLC1A2; SLC1A3; GJA1",
    "GAD2; KIT; TACR1; NXPH2",
    "SLC17A7; SATB2; TBR1; CUX1",
    "SLC17A7; BCL11B; TBR1; SATB2",
    "CLDN5; VWF; PECAM1; RGS5; PDGFRB",
    "CD3E; CD2; IL7R; LCK; ZAP70; THEMIS"
  ),
  
  annotation_note = c(
    "",
    "",
    "Do not assign reactive or senescence state at Level 1",
    "",
    "",
    "",
    "",
    "",
    "Fine neuronal subtype remains unresolved",
    "Fine neuronal subtype remains unresolved",
    "",
    "",
    "",
    "Possible reactive-associated features; not formally named",
    "",
    "Fine neuronal subtype remains unresolved",
    "Fine neuronal subtype remains unresolved",
    "Endothelial and pericyte signals are combined; requires subclustering",
    "Small cluster; no fine immune subtype assigned"
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  annotation_table,
  "04_tables/human_harmony_r0p4_level1_annotation_evidence.csv",
  row.names = FALSE
)

## ------------------------------------------------------------
## 3. 将注释写入对象
## ------------------------------------------------------------

cluster_values <- as.character(
  human@meta.data[[cluster_column]]
)

annotation_index <- match(
  cluster_values,
  annotation_table$cluster
)

if (anyNA(annotation_index)) {
  missing_clusters <- unique(
    cluster_values[
      is.na(annotation_index)
    ]
  )
  
  stop(
    paste(
      "以下cluster没有注释：",
      paste(missing_clusters, collapse = ", ")
    )
  )
}

human$annotation_level1_manual <-
  annotation_table$annotation_level1[
    annotation_index
  ]

human$annotation_level1_confidence <-
  annotation_table$confidence[
    annotation_index
  ]

human$annotation_level1_note <-
  annotation_table$annotation_note[
    annotation_index
  ]

annotation_order <- c(
  "Excitatory neurons",
  "Inhibitory neurons",
  "Astrocytes",
  "Oligodendrocytes",
  "OPC",
  "Microglia/macrophage",
  "Vascular cells",
  "T lymphocytes"
)

human$annotation_level1_manual <- factor(
  human$annotation_level1_manual,
  levels = annotation_order
)

## ------------------------------------------------------------
## 4. 按细胞大类和供者汇总
## ------------------------------------------------------------

metadata <- human@meta.data %>%
  tibble::rownames_to_column("cell_id")

annotation_counts <- metadata %>%
  count(
    annotation_level1_manual,
    name = "n_cells"
  ) %>%
  mutate(
    percent_of_all_cells =
      100 * n_cells / sum(n_cells)
  )

annotation_donor_counts <- metadata %>%
  count(
    annotation_level1_manual,
    donor_id,
    sample_id,
    condition,
    name = "n_cells"
  )

annotation_donor_presence <- annotation_donor_counts %>%
  group_by(annotation_level1_manual) %>%
  summarise(
    n_cells = sum(n_cells),
    n_donors_present =
      n_distinct(donor_id[n_cells > 0]),
    n_normal_donors =
      n_distinct(
        donor_id[
          condition == "Normal" &
            n_cells > 0
        ]
      ),
    n_tle_donors =
      n_distinct(
        donor_id[
          condition == "TLE" &
            n_cells > 0
        ]
      ),
    maximum_donor_fraction =
      max(n_cells) / sum(n_cells),
    .groups = "drop"
  )

write.csv(
  annotation_counts,
  "04_tables/human_level1_cell_counts.csv",
  row.names = FALSE
)

write.csv(
  annotation_donor_counts,
  "04_tables/human_level1_donor_cell_counts.csv",
  row.names = FALSE
)

write.csv(
  annotation_donor_presence,
  "04_tables/human_level1_donor_presence_summary.csv",
  row.names = FALSE
)

cat("========== Level 1细胞数量 ==========\n")

print(
  annotation_counts,
  row.names = FALSE
)

cat("\n========== Level 1供者覆盖 ==========\n")

print(
  annotation_donor_presence,
  row.names = FALSE
)

## ------------------------------------------------------------
## 5. Level 1 UMAP
## ------------------------------------------------------------

p_annotation <- DimPlot(
  object = human,
  reduction = "umap.harmony",
  group.by = "annotation_level1_manual",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  labs(
    title = "Manual Level 1 annotation",
    subtitle = "Harmony space used for clustering; RNA markers used for annotation"
  ) +
  theme(
    legend.position = "right"
  )

ggsave(
  "03_figures/human_harmony_level1_annotation_umap.png",
  p_annotation,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  "03_figures/human_harmony_level1_annotation_umap.pdf",
  p_annotation,
  width = 12,
  height = 8
)

## ------------------------------------------------------------
## 6. Level 1 UMAP按condition拆分
## ------------------------------------------------------------

p_annotation_condition <- DimPlot(
  object = human,
  reduction = "umap.harmony",
  group.by = "annotation_level1_manual",
  split.by = "condition",
  ncol = 2,
  raster = TRUE
) +
  labs(
    title = "Level 1 annotation split by condition"
  )

ggsave(
  "03_figures/human_harmony_level1_annotation_split_condition.png",
  p_annotation_condition,
  width = 16,
  height = 8,
  dpi = 300,
  bg = "white"
)

## ------------------------------------------------------------
## 7. Level 1经典marker DotPlot
## ------------------------------------------------------------

level1_marker_panel <- c(
  "SLC17A7",
  "SATB2",
  "GAD1",
  "GAD2",
  "AQP4",
  "ALDH1L1",
  "PLP1",
  "MOBP",
  "PDGFRA",
  "VCAN",
  "P2RY12",
  "INPP5D",
  "CLDN5",
  "PECAM1",
  "RGS5",
  "PDGFRB",
  "CD3D",
  "CD3E"
)

level1_marker_panel <- intersect(
  level1_marker_panel,
  rownames(human)
)

p_level1_dotplot <- DotPlot(
  object = human,
  features = level1_marker_panel,
  assay = "RNA",
  group.by = "annotation_level1_manual",
  dot.scale = 7
) +
  RotatedAxis() +
  labs(
    title = "Canonical markers across Level 1 annotations",
    x = NULL,
    y = NULL
  )

ggsave(
  "03_figures/human_level1_canonical_marker_dotplot.png",
  p_level1_dotplot,
  width = 14,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  "03_figures/human_level1_canonical_marker_dotplot.pdf",
  p_level1_dotplot,
  width = 14,
  height = 6
)

## ------------------------------------------------------------
## 8. 保存带Level 1注释的对象
## ------------------------------------------------------------

output_file <- file.path(
  "02_objects",
  "human_harmony_r0p4_level1_annotated.rds"
)

saveRDS(
  human,
  output_file,
  compress = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  "05_logs/15_human_assign_level1_annotation_sessionInfo.txt"
)

cat("\n========================================\n")
cat("Level 1宽细胞类型注释完成。\n")
cat("注释类别数：", nlevels(human$annotation_level1_manual), "\n")
cat("对象保存于：\n")
cat(output_file, "\n")
cat("\n主要图片：\n")
cat("03_figures/human_harmony_level1_annotation_umap.png\n")
cat("03_figures/human_level1_canonical_marker_dotplot.png\n")
cat("========================================\n")

View(annotation_table)
View(annotation_counts)
View(annotation_donor_presence)


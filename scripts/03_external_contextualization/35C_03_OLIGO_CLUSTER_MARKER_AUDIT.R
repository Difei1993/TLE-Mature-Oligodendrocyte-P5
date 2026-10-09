###############################################################################
# 35C-03
# GSE140393 oligodendroglial cluster marker audit
#
# Purpose:
#   - Calculate cluster-level differential markers
#   - Generate canonical-marker DotPlot
#   - Audit donor composition
#   - Produce a manual annotation worksheet
#
# Important:
#   - No P5 score is calculated here
#   - No cell is removed automatically
#   - The previous automatic V2 stage labels are not used
###############################################################################

options(
  stringsAsFactors = FALSE,
  warn = 1
)

set.seed(20260731)

###############################################################################
# 1. Packages
###############################################################################

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tibble",
  "readr",
  "ggplot2",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

###############################################################################
# 2. Paths
###############################################################################

project_root <- Sys.getenv("P5_WORK_ROOT")

work_root <- file.path(
  project_root,
  "35C_P5_CROSS_REFERENCE"
)

input_object <- file.path(
  work_root,
  "03_GSE140393_QC_AND_OLIGO_ANNOTATION",
  "objects",
  "GSE140393_TLE_OLIGODENDROGLIA_ANNOTATED.rds"
)

output_root <- file.path(
  work_root,
  "05_GSE140393_OLIGO_MANUAL_ANNOTATION_AUDIT"
)

table_root <- file.path(output_root, "tables")
figure_root <- file.path(output_root, "figures")
log_root <- file.path(output_root, "logs")

for (directory in c(
  output_root,
  table_root,
  figure_root,
  log_root
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (!file.exists(input_object)) {
  stop(
    "Input object not found: ",
    input_object
  )
}

###############################################################################
# 3. Read object
###############################################################################

message("Reading oligodendroglial object.")

oligo <- readRDS(input_object)

DefaultAssay(oligo) <- "RNA"

if (
  packageVersion("SeuratObject") >=
  package_version("5.0.0")
) {
  oligo <- tryCatch(
    JoinLayers(oligo),
    error = function(error) {
      message(
        "JoinLayers skipped safely: ",
        conditionMessage(error)
      )
      oligo
    }
  )
}

required_metadata <- c(
  "seurat_clusters",
  "donor_id",
  "sample_id"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(oligo@meta.data)
)

if (length(missing_metadata) > 0L) {
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

cluster_values <- as.character(
  oligo$seurat_clusters
)

if (anyNA(cluster_values)) {
  stop("Some cells have missing cluster labels.")
}

cluster_numeric <- suppressWarnings(
  as.integer(cluster_values)
)

if (all(!is.na(cluster_numeric))) {
  cluster_levels <- as.character(
    sort(unique(cluster_numeric))
  )
} else {
  cluster_levels <- sort(
    unique(cluster_values)
  )
}

Idents(oligo) <- factor(
  cluster_values,
  levels = cluster_levels
)

message("Cells: ", ncol(oligo))
message("Clusters: ", paste(cluster_levels, collapse = ", "))

###############################################################################
# 4. Cluster and donor composition
###############################################################################

cell_metadata <- oligo@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  mutate(
    seurat_cluster = as.character(seurat_clusters),
    donor_id = as.character(donor_id),
    sample_id = as.character(sample_id)
  )

cluster_counts <- cell_metadata %>%
  count(
    seurat_cluster,
    name = "cell_count"
  ) %>%
  arrange(
    suppressWarnings(
      as.integer(seurat_cluster)
    )
  )

donor_cluster_counts <- cell_metadata %>%
  count(
    seurat_cluster,
    donor_id,
    name = "cell_count"
  ) %>%
  group_by(seurat_cluster) %>%
  mutate(
    cluster_total = sum(cell_count),
    donor_fraction = cell_count / cluster_total
  ) %>%
  ungroup() %>%
  arrange(
    suppressWarnings(
      as.integer(seurat_cluster)
    ),
    desc(cell_count)
  )

dominant_donor <- donor_cluster_counts %>%
  group_by(seurat_cluster) %>%
  slice_max(
    order_by = cell_count,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  transmute(
    seurat_cluster,
    dominant_donor = donor_id,
    dominant_donor_cells = cell_count,
    dominant_donor_fraction = donor_fraction
  )

write_csv(
  cluster_counts,
  file.path(
    table_root,
    "GSE140393_OLIGO_CLUSTER_CELL_COUNTS.csv"
  )
)

write_csv(
  donor_cluster_counts,
  file.path(
    table_root,
    "GSE140393_OLIGO_DONOR_CLUSTER_COUNTS.csv"
  )
)

###############################################################################
# 5. Differential marker analysis
###############################################################################

message("Calculating cluster markers.")

marker_table <- FindAllMarkers(
  object = oligo,
  assay = "RNA",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25,
  return.thresh = 1,
  verbose = FALSE
)

if (nrow(marker_table) == 0L) {
  stop("FindAllMarkers returned no markers.")
}

fold_change_column <- intersect(
  c(
    "avg_log2FC",
    "avg_logFC"
  ),
  colnames(marker_table)
)

if (length(fold_change_column) == 0L) {
  stop(
    "No avg_log2FC or avg_logFC column was found."
  )
}

fold_change_column <- fold_change_column[1L]

marker_table <- marker_table %>%
  mutate(
    cluster = as.character(cluster),
    gene = as.character(gene)
  ) %>%
  arrange(
    suppressWarnings(
      as.integer(cluster)
    ),
    desc(.data[[fold_change_column]])
  )

write_csv(
  marker_table,
  file.path(
    table_root,
    "GSE140393_OLIGO_ALL_CLUSTER_MARKERS.csv"
  )
)

top_marker_table <- marker_table %>%
  group_by(cluster) %>%
  slice_max(
    order_by = .data[[fold_change_column]],
    n = 25,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    suppressWarnings(
      as.integer(cluster)
    ),
    desc(.data[[fold_change_column]])
  )

write_csv(
  top_marker_table,
  file.path(
    table_root,
    "GSE140393_OLIGO_TOP25_MARKERS_BY_CLUSTER.csv"
  )
)

top_marker_summary <- top_marker_table %>%
  group_by(cluster) %>%
  summarise(
    top_15_markers = paste(
      head(gene, 15L),
      collapse = "|"
    ),
    .groups = "drop"
  ) %>%
  rename(
    seurat_cluster = cluster
  )

###############################################################################
# 6. Canonical marker panel
###############################################################################

marker_groups <- list(
  Oligodendroglial_identity = c(
    "SOX10",
    "OLIG1",
    "OLIG2"
  ),
  
  OPC = c(
    "PDGFRA",
    "CSPG4",
    "VCAN",
    "PTPRZ1"
  ),
  
  Transitional_oligodendrocyte = c(
    "BCAS1",
    "TCF7L2",
    "ENPP6",
    "GPR17"
  ),
  
  Mature_oligodendrocyte = c(
    "PLP1",
    "CNP",
    "MOBP",
    "MOG",
    "MBP",
    "MAG",
    "CLDN11",
    "ERMN"
  ),
  
  Astrocyte = c(
    "AQP4",
    "ALDH1L1",
    "SLC1A2",
    "SLC1A3",
    "GFAP"
  ),
  
  Microglia = c(
    "P2RY12",
    "TMEM119",
    "C1QA",
    "C1QB",
    "TYROBP",
    "AIF1"
  ),
  
  Neuron = c(
    "SNAP25",
    "RBFOX3",
    "SYT1",
    "SLC17A7",
    "GAD1",
    "GAD2"
  ),
  
  Endothelial_pericyte = c(
    "CLDN5",
    "PECAM1",
    "VWF",
    "KDR",
    "RGS5",
    "PDGFRB"
  ),
  
  P5_representative = c(
    "CRYAB",
    "APOD",
    "PTGDS",
    "SCD",
    "SELENOP",
    "FTH1",
    "FTL"
  )
)

marker_panel <- unique(
  unlist(marker_groups)
)

marker_panel_present <- intersect(
  marker_panel,
  rownames(oligo)
)

marker_panel_missing <- setdiff(
  marker_panel,
  rownames(oligo)
)

marker_registry <- tibble(
  gene = marker_panel,
  present = marker_panel %in% rownames(oligo)
)

write_csv(
  marker_registry,
  file.path(
    table_root,
    "GSE140393_OLIGO_CANONICAL_MARKER_REGISTRY.csv"
  )
)

if (length(marker_panel_present) == 0L) {
  stop("No canonical marker genes were detected.")
}

###############################################################################
# 7. Marker DotPlot
###############################################################################

marker_dotplot <- DotPlot(
  object = oligo,
  features = marker_panel_present,
  group.by = "seurat_clusters",
  dot.scale = 5,
  scale = TRUE
) +
  RotatedAxis() +
  labs(
    title = "GSE140393 oligodendroglial candidate clusters",
    subtitle = paste0(
      "Canonical lineage, maturation, contaminant and ",
      "P5-representative markers"
    ),
    x = NULL,
    y = "Seurat cluster"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1
    )
  )

ggsave(
  filename = file.path(
    figure_root,
    "GSE140393_OLIGO_CLUSTER_CANONICAL_MARKER_DOTPLOT.png"
  ),
  plot = marker_dotplot,
  width = 20,
  height = 9,
  dpi = 320
)

ggsave(
  filename = file.path(
    figure_root,
    "GSE140393_OLIGO_CLUSTER_CANONICAL_MARKER_DOTPLOT.pdf"
  ),
  plot = marker_dotplot,
  width = 20,
  height = 9
)

###############################################################################
# 8. UMAP plots
###############################################################################

cluster_umap <- DimPlot(
  oligo,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Oligodendroglial candidate clusters")

donor_umap <- DimPlot(
  oligo,
  reduction = "umap",
  group.by = "donor_id"
) +
  ggtitle("Donor composition")

combined_umap <- cluster_umap + donor_umap

ggsave(
  filename = file.path(
    figure_root,
    "GSE140393_OLIGO_CLUSTER_AND_DONOR_UMAP.png"
  ),
  plot = combined_umap,
  width = 14,
  height = 6,
  dpi = 320
)

ggsave(
  filename = file.path(
    figure_root,
    "GSE140393_OLIGO_CLUSTER_AND_DONOR_UMAP.pdf"
  ),
  plot = combined_umap,
  width = 14,
  height = 6
)

###############################################################################
# 9. Manual annotation worksheet
###############################################################################

annotation_worksheet <- cluster_counts %>%
  left_join(
    dominant_donor,
    by = "seurat_cluster"
  ) %>%
  left_join(
    top_marker_summary,
    by = "seurat_cluster"
  ) %>%
  mutate(
    manual_label = "",
    include_in_formal_p5_analysis = "",
    annotation_reason = "",
    reviewer_notes = ""
  ) %>%
  arrange(
    suppressWarnings(
      as.integer(seurat_cluster)
    )
  )

write_csv(
  annotation_worksheet,
  file.path(
    table_root,
    "GSE140393_OLIGO_MANUAL_ANNOTATION_WORKSHEET.csv"
  )
)

###############################################################################
# 10. Console summary
###############################################################################

console_summary <- annotation_worksheet %>%
  select(
    seurat_cluster,
    cell_count,
    dominant_donor,
    dominant_donor_fraction,
    top_15_markers
  )

message("")
message("===================================================")
message("Oligodendroglial cluster marker audit completed")
message("===================================================")
message("")
message("Cluster summary:")

print(
  as.data.frame(console_summary),
  row.names = FALSE
)

message("")
message("Missing canonical marker genes:")

if (length(marker_panel_missing) == 0L) {
  message("None")
} else {
  message(
    paste(
      marker_panel_missing,
      collapse = ", "
    )
  )
}

message("")
message("Outputs:")
message(
  file.path(
    table_root,
    "GSE140393_OLIGO_TOP25_MARKERS_BY_CLUSTER.csv"
  )
)
message(
  file.path(
    table_root,
    "GSE140393_OLIGO_MANUAL_ANNOTATION_WORKSHEET.csv"
  )
)
message(
  file.path(
    figure_root,
    "GSE140393_OLIGO_CLUSTER_CANONICAL_MARKER_DOTPLOT.png"
  )
)

###############################################################################
# 11. Session information
###############################################################################

capture.output(
  sessionInfo(),
  file = file.path(
    log_root,
    "sessionInfo.txt"
  )
)
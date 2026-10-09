###############################################################################
# 35C-04
# Apply manually reviewed oligodendroglial annotations to GSE140393
#
# No clustering, QC or P5 scoring is performed in this script.
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
  "06_GSE140393_FINAL_MANUAL_OLIGO_ANNOTATION"
)

table_root <- file.path(output_root, "tables")
figure_root <- file.path(output_root, "figures")
object_root <- file.path(output_root, "objects")
log_root <- file.path(output_root, "logs")

for (directory in c(
  output_root,
  table_root,
  figure_root,
  object_root,
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

message("Reading GSE140393 candidate oligodendroglial object.")

oligo_candidate <- readRDS(input_object)

DefaultAssay(oligo_candidate) <- "RNA"

if (
  packageVersion("SeuratObject") >=
  package_version("5.0.0")
) {
  oligo_candidate <- tryCatch(
    JoinLayers(oligo_candidate),
    error = function(error) {
      message(
        "JoinLayers skipped safely: ",
        conditionMessage(error)
      )
      
      oligo_candidate
    }
  )
}

if (
  !"seurat_clusters" %in%
  colnames(oligo_candidate@meta.data)
) {
  stop("seurat_clusters is absent from metadata.")
}

###############################################################################
# 4. Frozen manual annotation map
###############################################################################

annotation_map <- tibble::tribble(
  ~seurat_cluster, ~manual_label,
  ~broad_class, ~formal_p5_role, ~annotation_evidence,
  
  "0", "Mature_oligodendrocyte",
  "Oligodendroglia", "PRIMARY_MATURE",
  "NKX6-2|HHIP|C21orf91|ABCA8",
  
  "1", "OPC",
  "Oligodendroglia", "SECONDARY_OPC",
  "PDGFRA|PTPRZ1|CSPG4|VCAN|ASCL1",
  
  "2", "T_NK_cell",
  "Contaminant", "EXCLUDE",
  "KLRD1|TRDC|PRF1|CD3G|CCL5",
  
  "3", "Microglia_macrophage",
  "Contaminant", "EXCLUDE",
  "C1QA|C1QB|C1QC|TREM2|AIF1",
  
  "4", "Endothelial_cell",
  "Contaminant", "EXCLUDE",
  "CLDN5|CA4|ESAM|MFSD2A",
  
  "5", "Mural_fibroblast_cell",
  "Contaminant", "EXCLUDE",
  "ACTA2|TAGLN|CNN1|COL14A1",
  
  "6", "GPR17_transitional_oligodendrocyte",
  "Oligodendroglia", "EXPLORATORY_TRANSITIONAL",
  "GPR17|SLC1A1|TMEM108",
  
  "7", "Plasma_cell",
  "Contaminant", "EXCLUDE",
  "IGHG1|IGHG3|MZB1|JCHAIN",
  
  "8", "ENPP6_maturing_oligodendrocyte",
  "Oligodendroglia", "EXPLORATORY_MATURING",
  "ENPP6|HAPLN2|SERPINI1",
  
  "9", "Fibroblast",
  "Contaminant", "EXCLUDE",
  "COL1A1|COL3A1|LUM|FMOD|ASPN",
  
  "10", "Cycling_lineage_unresolved",
  "Unresolved", "EXCLUDE",
  "MKI67|TOP2A|CDK1|BIRC5",
  
  "11", "B_cell",
  "Contaminant", "EXCLUDE",
  "CD79A|CD79B|MS4A1|CD19",
  
  "12", "Cycling_lineage_unresolved",
  "Unresolved", "EXCLUDE",
  "FOXM1|UBE2C|CCNB1|TOP2A",
  
  "13", "Endothelial_cell",
  "Contaminant", "EXCLUDE",
  "MFSD2A|CLEC14A|RAMP2|ABCG2",
  
  "14", "Plasmacytoid_dendritic_cell",
  "Contaminant", "EXCLUDE",
  "LILRA4|SPIB|IRF4|GZ?",
  
  "15", "Microglia_vascular_mixed",
  "Contaminant", "EXCLUDE",
  "TREM2|C1QA|IFI27|SLC9A3R2"
)

annotation_map <- annotation_map %>%
  mutate(
    seurat_cluster = as.character(seurat_cluster)
  )

observed_clusters <- sort(
  unique(
    as.character(
      oligo_candidate$seurat_clusters
    )
  )
)

missing_annotations <- setdiff(
  observed_clusters,
  annotation_map$seurat_cluster
)

extra_annotations <- setdiff(
  annotation_map$seurat_cluster,
  observed_clusters
)

if (length(missing_annotations) > 0L) {
  stop(
    "The following observed clusters are not annotated: ",
    paste(missing_annotations, collapse = ", ")
  )
}

if (length(extra_annotations) > 0L) {
  stop(
    "The annotation map contains absent clusters: ",
    paste(extra_annotations, collapse = ", ")
  )
}

###############################################################################
# 5. Map annotations to cells
###############################################################################

cluster_index <- match(
  as.character(
    oligo_candidate$seurat_clusters
  ),
  annotation_map$seurat_cluster
)

oligo_candidate$manual_label_35c <-
  annotation_map$manual_label[cluster_index]

oligo_candidate$broad_class_35c <-
  annotation_map$broad_class[cluster_index]

oligo_candidate$formal_p5_role_35c <-
  annotation_map$formal_p5_role[cluster_index]

if (
  anyNA(oligo_candidate$manual_label_35c) ||
  anyNA(oligo_candidate$formal_p5_role_35c)
) {
  stop("Manual annotation mapping produced missing values.")
}

oligo_candidate$is_bona_fide_oligodendroglia_35c <-
  oligo_candidate$broad_class_35c ==
  "Oligodendroglia"

###############################################################################
# 6. Create clean subsets
###############################################################################

bona_fide_cells <- colnames(
  oligo_candidate
)[
  oligo_candidate$is_bona_fide_oligodendroglia_35c
]

mature_cells <- colnames(
  oligo_candidate
)[
  oligo_candidate$formal_p5_role_35c ==
    "PRIMARY_MATURE"
]

opc_cells <- colnames(
  oligo_candidate
)[
  oligo_candidate$formal_p5_role_35c ==
    "SECONDARY_OPC"
]

transitional_cells <- colnames(
  oligo_candidate
)[
  oligo_candidate$formal_p5_role_35c %in%
    c(
      "EXPLORATORY_TRANSITIONAL",
      "EXPLORATORY_MATURING"
    )
]

if (length(bona_fide_cells) == 0L) {
  stop("No bona fide oligodendroglial cells were selected.")
}

if (length(mature_cells) == 0L) {
  stop("No mature oligodendrocytes were selected.")
}

if (length(opc_cells) == 0L) {
  stop("No OPC cells were selected.")
}

oligo_clean <- subset(
  oligo_candidate,
  cells = bona_fide_cells
)

mature_oligo <- subset(
  oligo_candidate,
  cells = mature_cells
)

opc <- subset(
  oligo_candidate,
  cells = opc_cells
)

transitional_oligo <- subset(
  oligo_candidate,
  cells = transitional_cells
)

###############################################################################
# 7. Cell-count tables
###############################################################################

candidate_cluster_counts <- oligo_candidate@meta.data %>%
  count(
    seurat_clusters,
    manual_label_35c,
    broad_class_35c,
    formal_p5_role_35c,
    name = "cell_count"
  ) %>%
  mutate(
    seurat_clusters =
      as.character(seurat_clusters)
  ) %>%
  left_join(
    annotation_map %>%
      select(
        seurat_cluster,
        annotation_evidence
      ),
    by = c(
      "seurat_clusters" =
        "seurat_cluster"
    )
  ) %>%
  arrange(
    suppressWarnings(
      as.integer(seurat_clusters)
    )
  )

donor_label_counts <- oligo_candidate@meta.data %>%
  count(
    sample_id,
    donor_id,
    age,
    sex,
    age_group,
    original_paper_status,
    manual_label_35c,
    formal_p5_role_35c,
    name = "cell_count"
  ) %>%
  arrange(
    donor_id,
    formal_p5_role_35c
  )

donor_formal_counts <- oligo_candidate@meta.data %>%
  filter(
    formal_p5_role_35c %in%
      c(
        "PRIMARY_MATURE",
        "SECONDARY_OPC",
        "EXPLORATORY_TRANSITIONAL",
        "EXPLORATORY_MATURING"
      )
  ) %>%
  count(
    sample_id,
    donor_id,
    age,
    sex,
    age_group,
    original_paper_status,
    formal_p5_role_35c,
    name = "cell_count"
  ) %>%
  mutate(
    formal_donor_eligible =
      cell_count >= 30L
  ) %>%
  arrange(
    formal_p5_role_35c,
    donor_id
  )

write_csv(
  annotation_map,
  file.path(
    table_root,
    "GSE140393_MANUAL_CLUSTER_ANNOTATION_MAP.csv"
  )
)

write_csv(
  candidate_cluster_counts,
  file.path(
    table_root,
    "GSE140393_FINAL_CLUSTER_CELL_COUNTS.csv"
  )
)

write_csv(
  donor_label_counts,
  file.path(
    table_root,
    "GSE140393_DONOR_MANUAL_LABEL_COUNTS.csv"
  )
)

write_csv(
  donor_formal_counts,
  file.path(
    table_root,
    "GSE140393_DONOR_FORMAL_P5_COMPARTMENT_COUNTS.csv"
  )
)

###############################################################################
# 8. White-background diagnostic plots
###############################################################################

manual_umap <- DimPlot(
  object = oligo_candidate,
  reduction = "umap",
  group.by = "manual_label_35c",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Final manual annotation") +
  theme_classic()

clean_umap <- DimPlot(
  object = oligo_clean,
  reduction = "umap",
  group.by = "manual_label_35c",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Bona fide oligodendroglial cells") +
  theme_classic()

combined_umap <- manual_umap + clean_umap

ggsave(
  filename = file.path(
    figure_root,
    "GSE140393_FINAL_MANUAL_ANNOTATION_UMAP.png"
  ),
  plot = combined_umap,
  width = 17,
  height = 7,
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = file.path(
    figure_root,
    "GSE140393_FINAL_MANUAL_ANNOTATION_UMAP.pdf"
  ),
  plot = combined_umap,
  width = 17,
  height = 7,
  bg = "white"
)

diagnostic_markers <- c(
  "SOX10", "OLIG1", "OLIG2",
  "PDGFRA", "CSPG4", "VCAN", "PTPRZ1",
  "GPR17", "BCAS1", "TCF7L2", "ENPP6",
  "PLP1", "CNP", "MOBP", "MOG",
  "MBP", "MAG", "CLDN11", "ERMN",
  "CRYAB", "APOD", "PTGDS", "SCD",
  "SELENOP", "FTH1", "FTL"
)

diagnostic_markers <- intersect(
  diagnostic_markers,
  rownames(oligo_clean)
)

Idents(oligo_clean) <-
  oligo_clean$manual_label_35c

diagnostic_dotplot <- DotPlot(
  object = oligo_clean,
  features = diagnostic_markers,
  group.by = "manual_label_35c",
  dot.scale = 6,
  scale = TRUE
) +
  RotatedAxis() +
  ggtitle(
    "Final oligodendroglial-stage marker audit"
  ) +
  theme_classic() +
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
    "GSE140393_FINAL_OLIGO_MARKER_DOTPLOT.png"
  ),
  plot = diagnostic_dotplot,
  width = 15,
  height = 6,
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = file.path(
    figure_root,
    "GSE140393_FINAL_OLIGO_MARKER_DOTPLOT.pdf"
  ),
  plot = diagnostic_dotplot,
  width = 15,
  height = 6,
  bg = "white"
)

###############################################################################
# 9. Save objects
###############################################################################

saveRDS(
  oligo_candidate,
  file.path(
    object_root,
    "GSE140393_ALL_CANDIDATE_CELLS_MANUALLY_ANNOTATED.rds"
  ),
  compress = FALSE
)

saveRDS(
  oligo_clean,
  file.path(
    object_root,
    "GSE140393_BONA_FIDE_OLIGODENDROGLIA.rds"
  ),
  compress = FALSE
)

saveRDS(
  mature_oligo,
  file.path(
    object_root,
    "GSE140393_MATURE_OLIGODENDROCYTES.rds"
  ),
  compress = FALSE
)

saveRDS(
  opc,
  file.path(
    object_root,
    "GSE140393_OPC.rds"
  ),
  compress = FALSE
)

saveRDS(
  transitional_oligo,
  file.path(
    object_root,
    "GSE140393_TRANSITIONAL_AND_MATURING_OLIGO.rds"
  ),
  compress = FALSE
)

###############################################################################
# 10. Final checks
###############################################################################

expected_clean_count <- sum(
  candidate_cluster_counts$cell_count[
    candidate_cluster_counts$broad_class_35c ==
      "Oligodendroglia"
  ]
)

if (ncol(oligo_clean) != expected_clean_count) {
  stop(
    "Clean subset count does not match the annotation map."
  )
}

message("")
message("===================================================")
message("Final manual oligodendroglial annotation completed")
message("===================================================")
message("")
message("Candidate cells: ", ncol(oligo_candidate))
message("Bona fide oligodendroglial cells: ", ncol(oligo_clean))
message("Mature oligodendrocytes: ", ncol(mature_oligo))
message("OPC: ", ncol(opc))
message(
  "Transitional/maturing oligodendrocytes: ",
  ncol(transitional_oligo)
)

message("")
message("Donor-level formal compartment counts:")

print(
  as.data.frame(donor_formal_counts),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    log_root,
    "sessionInfo.txt"
  )
)
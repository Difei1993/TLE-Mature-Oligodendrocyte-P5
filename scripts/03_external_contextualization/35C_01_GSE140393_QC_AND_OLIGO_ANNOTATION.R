###############################################################################
# 35C-01
# GSE140393 five-TLE scRNA-seq:
# per-donor QC, clustering, major-lineage annotation,
# oligodendroglial extraction and maturation-stage annotation
###############################################################################

options(
  stringsAsFactors = FALSE,
  warn = 1,
  future.globals.maxSize = 12 * 1024^3
)

set.seed(20260731)

###############################################################################
# 1. Paths
###############################################################################

project_root <- Sys.getenv("P5_WORK_ROOT")

work_root <- file.path(
  project_root,
  "35C_P5_CROSS_REFERENCE"
)

input_root <- file.path(
  work_root,
  "data",
  "raw",
  "GSE140393"
)

sample_root <- file.path(
  input_root,
  "five_TLE_samples"
)

manifest_path <- file.path(
  input_root,
  "GSE140393_TLE_DONOR_MANIFEST.csv"
)

output_root <- file.path(
  work_root,
  "03_GSE140393_QC_AND_OLIGO_ANNOTATION"
)

figure_root <- file.path(output_root, "figures")
table_root <- file.path(output_root, "tables")
object_root <- file.path(output_root, "objects")
log_root <- file.path(output_root, "logs")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(table_root, recursive = TRUE, showWarnings = FALSE)
dir.create(object_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 2. Packages
###############################################################################

cran_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "data.table",
  "dplyr",
  "tibble",
  "readr",
  "ggplot2",
  "patchwork",
  "harmony",
  "future"
)

missing_packages <- cran_packages[
  !vapply(
    cran_packages,
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
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(future)
})

future::plan("sequential")

###############################################################################
# 3. Helper functions
###############################################################################
normalize_cluster_id <- function(x) {
  x <- as.character(x)
  
  sub(
    pattern = "^g(?=[0-9]+$)",
    replacement = "",
    x = x,
    perl = TRUE
  )
}
safe_mad <- function(x) {
  value <- stats::mad(
    x,
    center = stats::median(x, na.rm = TRUE),
    constant = 1,
    na.rm = TRUE
  )
  
  if (!is.finite(value) || value <= 0) {
    value <- stats::IQR(x, na.rm = TRUE) / 1.349
  }
  
  if (!is.finite(value) || value <= 0) {
    value <- stats::sd(x, na.rm = TRUE)
  }
  
  if (!is.finite(value) || value <= 0) {
    value <- 1
  }
  
  value
}

make_qc_limits <- function(metadata) {
  feature_median <- median(
    metadata$nFeature_RNA,
    na.rm = TRUE
  )
  
  count_median <- median(
    metadata$nCount_RNA,
    na.rm = TRUE
  )
  
  mt_median <- median(
    metadata$percent.mt,
    na.rm = TRUE
  )
  
  feature_mad <- safe_mad(metadata$nFeature_RNA)
  count_mad <- safe_mad(metadata$nCount_RNA)
  mt_mad <- safe_mad(metadata$percent.mt)
  
  feature_lower <- max(
    200,
    as.numeric(
      quantile(
        metadata$nFeature_RNA,
        0.01,
        na.rm = TRUE
      )
    ),
    feature_median - 4 * feature_mad
  )
  
  feature_upper <- min(
    as.numeric(
      quantile(
        metadata$nFeature_RNA,
        0.995,
        na.rm = TRUE
      )
    ),
    feature_median + 4 * feature_mad
  )
  
  count_upper <- min(
    as.numeric(
      quantile(
        metadata$nCount_RNA,
        0.995,
        na.rm = TRUE
      )
    ),
    count_median + 4 * count_mad
  )
  
  mt_upper <- min(
    25,
    max(
      10,
      mt_median + 4 * mt_mad
    )
  )
  
  data.frame(
    feature_lower = floor(feature_lower),
    feature_upper = ceiling(feature_upper),
    count_upper = ceiling(count_upper),
    mt_upper = mt_upper,
    stringsAsFactors = FALSE
  )
}

join_layers_if_needed <- function(object) {
  if (
    packageVersion("SeuratObject") >= package_version("5.0.0")
  ) {
    object <- tryCatch(
      JoinLayers(object),
      error = function(e) {
        message(
          "JoinLayers未执行：",
          conditionMessage(e)
        )
        object
      }
    )
  }
  
  object
}

existing_genes <- function(object, genes) {
  intersect(
    genes,
    rownames(object)
  )
}

mean_marker_expression <- function(
    average_matrix,
    marker_genes
) {
  genes <- intersect(
    marker_genes,
    rownames(average_matrix)
  )
  
  if (length(genes) == 0L) {
    return(
      rep(
        NA_real_,
        ncol(average_matrix)
      )
    )
  }
  
  colMeans(
    average_matrix[genes, , drop = FALSE],
    na.rm = TRUE
  )
}

save_plot <- function(
    plot,
    file_stub,
    width = 10,
    height = 7
) {
  ggsave(
    filename = file.path(
      figure_root,
      paste0(file_stub, ".png")
    ),
    plot = plot,
    width = width,
    height = height,
    dpi = 320
  )
  
  ggsave(
    filename = file.path(
      figure_root,
      paste0(file_stub, ".pdf")
    ),
    plot = plot,
    width = width,
    height = height
  )
}

###############################################################################
# 4. Read donor manifest
###############################################################################

if (!file.exists(manifest_path)) {
  stop(
    "Donor manifest not found: ",
    manifest_path
  )
}

manifest <- readr::read_csv(
  manifest_path,
  show_col_types = FALSE
) %>%
  mutate(
    sample_id = as.character(sample_id),
    donor_id = as.character(donor_id),
    age = as.numeric(age),
    sex = as.character(sex),
    age_group = as.character(age_group),
    tenx_version = as.character(tenx_version),
    original_paper_status =
      as.character(original_paper_status)
  )

expected_samples <- c(
  "GSM4160577",
  "GSM4160578",
  "GSM4160579",
  "GSM4160580",
  "GSM4160581"
)

if (!setequal(manifest$sample_id, expected_samples)) {
  stop(
    "Manifest sample set does not match the expected five samples."
  )
}

manifest <- manifest[
  match(expected_samples, manifest$sample_id),
  ,
  drop = FALSE
]

###############################################################################
# 5. Read each 10x matrix and perform per-donor QC
###############################################################################

sample_objects <- list()
qc_threshold_rows <- list()
qc_summary_rows <- list()

for (index in seq_len(nrow(manifest))) {
  sample_info <- manifest[index, ]
  
  sample_id <- sample_info$sample_id
  donor_id <- sample_info$donor_id
  
  sample_directory <- file.path(
    sample_root,
    sample_id
  )
  
  required_files <- file.path(
    sample_directory,
    c(
      "matrix.mtx.gz",
      "barcodes.tsv.gz",
      "features.tsv.gz"
    )
  )
  
  if (!all(file.exists(required_files))) {
    stop(
      "Incomplete 10x input for ",
      sample_id
    )
  }
  
  message("")
  message(
    "Reading ",
    sample_id,
    " / donor ",
    donor_id
  )
  
  counts <- Read10X(
    data.dir = sample_directory,
    gene.column = 2,
    unique.features = TRUE
  )
  
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1L]]
    }
  }
  
  if (!inherits(counts, "sparseMatrix")) {
    counts <- Matrix::Matrix(
      counts,
      sparse = TRUE
    )
  }
  
  object <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 3,
    min.features = 100
  )
  
  object <- RenameCells(
    object,
    add.cell.id = sample_id
  )
  
  object$sample_id <- sample_id
  object$donor_id <- donor_id
  object$age <- sample_info$age
  object$sex <- sample_info$sex
  object$age_group <- sample_info$age_group
  object$tenx_version <- sample_info$tenx_version
  object$original_paper_status <-
    sample_info$original_paper_status
  object$group <- "TLE"
  object$dataset <- "GSE140393"
  
  object[["percent.mt"]] <- PercentageFeatureSet(
    object,
    pattern = "^MT-"
  )
  
  limits <- make_qc_limits(
    object@meta.data
  )
  
  limits$sample_id <- sample_id
  limits$donor_id <- donor_id
  
  qc_threshold_rows[[sample_id]] <- limits
  
  before_cells <- ncol(object)
  
  keep_cells <- with(
    object@meta.data,
    nFeature_RNA >= limits$feature_lower &
      nFeature_RNA <= limits$feature_upper &
      nCount_RNA <= limits$count_upper &
      percent.mt <= limits$mt_upper
  )
  
  object$pass_qc_35c <- keep_cells
  
  object <- subset(
    object,
    subset = pass_qc_35c
  )
  
  after_cells <- ncol(object)
  
  qc_summary_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    donor_id = donor_id,
    cells_before_qc = before_cells,
    cells_after_qc = after_cells,
    retained_fraction = after_cells / before_cells,
    median_nFeature_RNA = median(
      object$nFeature_RNA
    ),
    median_nCount_RNA = median(
      object$nCount_RNA
    ),
    median_percent_mt = median(
      object$percent.mt
    ),
    stringsAsFactors = FALSE
  )
  
  sample_objects[[sample_id]] <- object
  
  rm(counts, object)
  invisible(gc())
}

qc_threshold_table <- bind_rows(
  qc_threshold_rows
) %>%
  select(
    sample_id,
    donor_id,
    everything()
  )

qc_summary_table <- bind_rows(
  qc_summary_rows
)

write_csv(
  qc_threshold_table,
  file.path(
    table_root,
    "GSE140393_PER_SAMPLE_QC_THRESHOLDS.csv"
  )
)

write_csv(
  qc_summary_table,
  file.path(
    table_root,
    "GSE140393_PER_SAMPLE_QC_SUMMARY.csv"
  )
)

###############################################################################
# 6. Merge five donors
###############################################################################

message("")
message("Merging five TLE donors.")

merged <- merge(
  x = sample_objects[[1L]],
  y = sample_objects[-1L],
  merge.data = FALSE,
  project = "GSE140393_TLE"
)

merged <- join_layers_if_needed(merged)

DefaultAssay(merged) <- "RNA"

###############################################################################
# 7. Standard preprocessing and donor-aware embedding
###############################################################################

merged <- NormalizeData(
  merged,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

merged <- FindVariableFeatures(
  merged,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

merged <- ScaleData(
  merged,
  vars.to.regress = c(
    "nCount_RNA",
    "percent.mt"
  ),
  verbose = FALSE
)

merged <- RunPCA(
  merged,
  npcs = 50,
  verbose = FALSE
)

merged <- RunHarmony(
  object = merged,
  group.by.vars = "donor_id",
  reduction.use = "pca",
  dims.use = 1:40,
  assay.use = "RNA",
  project.dim = FALSE,
  verbose = FALSE
)

merged <- RunUMAP(
  merged,
  reduction = "harmony",
  dims = 1:30,
  seed.use = 20260731,
  verbose = FALSE
)

merged <- FindNeighbors(
  merged,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

merged <- FindClusters(
  merged,
  resolution = 0.5,
  random.seed = 20260731,
  verbose = FALSE
)

###############################################################################
# 8. Major lineage marker audit
###############################################################################

major_marker_sets <- list(
  Oligodendroglia = c(
    "SOX10", "OLIG1", "OLIG2",
    "PLP1", "CNP", "MOBP",
    "MOG", "MBP", "MAG", "CLDN11"
  ),
  OPC = c(
    "PDGFRA", "CSPG4",
    "VCAN", "PTPRZ1"
  ),
  Astrocyte = c(
    "AQP4", "ALDH1L1",
    "SLC1A2", "SLC1A3",
    "GFAP", "GJA1"
  ),
  Microglia = c(
    "P2RY12", "TMEM119",
    "CX3CR1", "AIF1",
    "TYROBP", "C1QA"
  ),
  Excitatory_neuron = c(
    "SLC17A7", "CAMK2A",
    "SATB2", "CUX2"
  ),
  Inhibitory_neuron = c(
    "GAD1", "GAD2",
    "SLC6A1", "DLX1"
  ),
  Endothelial = c(
    "CLDN5", "PECAM1",
    "VWF", "FLT1", "KDR"
  ),
  Pericyte = c(
    "PDGFRB", "RGS5",
    "MCAM", "COL4A1"
  ),
  Ependymal = c(
    "FOXJ1", "PIFO",
    "TPPP3"
  )
)

major_markers_present <- unique(
  unlist(
    lapply(
      major_marker_sets,
      function(genes) {
        existing_genes(merged, genes)
      }
    )
  )
)

average_expression <- AverageExpression(
  merged,
  assays = "RNA",
  features = major_markers_present,
  group.by = "seurat_clusters",
  slot = "data",
  verbose = FALSE
)$RNA

cluster_score_table <- data.frame(
  cluster = normalize_cluster_id(
    colnames(average_expression)
  ),
  stringsAsFactors = FALSE
)

for (lineage in names(major_marker_sets)) {
  cluster_score_table[[lineage]] <-
    mean_marker_expression(
      average_expression,
      major_marker_sets[[lineage]]
    )
}

lineage_columns <- names(major_marker_sets)

cluster_score_table$preliminary_lineage <- apply(
  cluster_score_table[
    ,
    lineage_columns,
    drop = FALSE
  ],
  1,
  function(values) {
    if (all(is.na(values))) {
      return("Unresolved")
    }
    
    lineage_columns[
      which.max(values)
    ]
  }
)

cluster_score_table$cluster <- as.character(
  cluster_score_table$cluster
)

valid_global_clusters <- sort(
  unique(
    as.character(
      merged$seurat_clusters
    )
  )
)

unmatched_global_clusters <- setdiff(
  cluster_score_table$cluster,
  valid_global_clusters
)

if (length(unmatched_global_clusters) > 0L) {
  stop(
    "AverageExpression cluster labels do not match ",
    "Seurat cluster identities: ",
    paste(
      unmatched_global_clusters,
      collapse = ", "
    )
  )
}

cluster_to_lineage <- setNames(
  cluster_score_table$preliminary_lineage,
  cluster_score_table$cluster
)

merged$preliminary_lineage <- unname(
  cluster_to_lineage[
    as.character(
      merged$seurat_clusters
    )
  ]
)

if (anyNA(merged$preliminary_lineage)) {
  missing_clusters <- unique(
    as.character(
      merged$seurat_clusters[
        is.na(
          merged$preliminary_lineage
        )
      ]
    )
  )
  
  stop(
    "Some cells could not be assigned a preliminary lineage. ",
    "Affected clusters: ",
    paste(
      missing_clusters,
      collapse = ", "
    )
  )
}

write_csv(
  cluster_score_table,
  file.path(
    table_root,
    "GSE140393_MAJOR_LINEAGE_CLUSTER_SCORES.csv"
  )
)

###############################################################################
# 9. Major-lineage plots
###############################################################################

p_umap_cluster <- DimPlot(
  merged,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("GSE140393 TLE: clusters")

p_umap_donor <- DimPlot(
  merged,
  reduction = "umap",
  group.by = "donor_id"
) +
  ggtitle("GSE140393 TLE: donors")

p_umap_lineage <- DimPlot(
  merged,
  reduction = "umap",
  group.by = "preliminary_lineage",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Preliminary major-lineage annotation")

save_plot(
  p_umap_cluster +
    p_umap_donor +
    p_umap_lineage,
  "GSE140393_GLOBAL_UMAP",
  width = 18,
  height = 6
)

dotplot_genes <- unique(
  unlist(major_marker_sets)
)

dotplot_genes <- existing_genes(
  merged,
  dotplot_genes
)

p_dotplot <- DotPlot(
  merged,
  features = dotplot_genes,
  group.by = "seurat_clusters",
  dot.scale = 5
) +
  RotatedAxis() +
  ggtitle("Major-lineage marker audit")

save_plot(
  p_dotplot,
  "GSE140393_MAJOR_LINEAGE_DOTPLOT",
  width = 18,
  height = 8
)

###############################################################################
# 10. Extract oligodendroglial lineage
###############################################################################

oligodendroglial_clusters <- cluster_score_table %>%
  filter(
    preliminary_lineage %in%
      c(
        "Oligodendroglia",
        "OPC"
      )
  ) %>%
  pull(cluster) %>%
  as.character() %>%
  unique()

valid_global_clusters <- unique(
  as.character(
    merged$seurat_clusters
  )
)

invalid_oligodendroglial_clusters <- setdiff(
  oligodendroglial_clusters,
  valid_global_clusters
)

if (
  length(
    invalid_oligodendroglial_clusters
  ) > 0L
) {
  stop(
    "Invalid oligodendroglial cluster labels: ",
    paste(
      invalid_oligodendroglial_clusters,
      collapse = ", "
    )
  )
}

if (length(oligodendroglial_clusters) == 0L) {
  stop(
    "No oligodendroglial clusters were identified."
  )
}

oligodendroglial_cells <- rownames(
  merged@meta.data
)[
  as.character(
    merged$seurat_clusters
  ) %in% oligodendroglial_clusters
]

if (length(oligodendroglial_cells) == 0L) {
  stop(
    "Oligodendroglial cluster labels were identified, ",
    "but no matching cells were found."
  )
}

message(
  "Selected oligodendroglial clusters: ",
  paste(
    sort(
      oligodendroglial_clusters
    ),
    collapse = ", "
  )
)

message(
  "Selected oligodendroglial cells: ",
  length(
    oligodendroglial_cells
  )
)

oligo <- subset(
  merged,
  cells = oligodendroglial_cells
)

oligo <- join_layers_if_needed(oligo)

DefaultAssay(oligo) <- "RNA"

###############################################################################
# 11. Reprocess oligodendroglial subset
###############################################################################

oligo <- NormalizeData(
  oligo,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

oligo <- FindVariableFeatures(
  oligo,
  selection.method = "vst",
  nfeatures = min(
    2500,
    nrow(oligo)
  ),
  verbose = FALSE
)

oligo <- ScaleData(
  oligo,
  vars.to.regress = c(
    "nCount_RNA",
    "percent.mt"
  ),
  verbose = FALSE
)

oligo <- RunPCA(
  oligo,
  npcs = 40,
  verbose = FALSE
)

oligo <- RunHarmony(
  object = oligo,
  group.by.vars = "donor_id",
  reduction.use = "pca",
  dims.use = 1:30,
  assay.use = "RNA",
  project.dim = FALSE,
  verbose = FALSE
)

oligo <- RunUMAP(
  oligo,
  reduction = "harmony",
  dims = 1:25,
  seed.use = 20260731,
  verbose = FALSE
)

oligo <- FindNeighbors(
  oligo,
  reduction = "harmony",
  dims = 1:25,
  verbose = FALSE
)

oligo <- FindClusters(
  oligo,
  resolution = 0.4,
  random.seed = 20260731,
  verbose = FALSE
)

###############################################################################
# 12. Oligodendroglial maturation-stage annotation
###############################################################################

stage_marker_sets <- list(
  OPC = c(
    "PDGFRA", "CSPG4",
    "VCAN", "PTPRZ1"
  ),
  Transitional_oligodendrocyte = c(
    "BCAS1", "TCF7L2",
    "ENPP6", "GPR17",
    "SOX6"
  ),
  Mature_oligodendrocyte = c(
    "MOG", "MOBP",
    "MBP", "PLP1",
    "CNP", "MAG",
    "CLDN11", "ERMN"
  )
)

p5_representative_genes <- c(
  "CRYAB",
  "APOD",
  "PTGDS",
  "SCD",
  "SELENOP",
  "FTH1",
  "FTL",
  "PLP1",
  "CNP",
  "ERMN"
)

stage_features <- unique(
  c(
    unlist(stage_marker_sets),
    p5_representative_genes
  )
)

stage_features <- existing_genes(
  oligo,
  stage_features
)

oligo_average <- AverageExpression(
  oligo,
  assays = "RNA",
  features = stage_features,
  group.by = "seurat_clusters",
  slot = "data",
  verbose = FALSE
)$RNA

stage_score_table <- data.frame(
  oligo_cluster = normalize_cluster_id(
    colnames(oligo_average)
  ),
  stringsAsFactors = FALSE
)

for (stage in names(stage_marker_sets)) {
  stage_score_table[[stage]] <-
    mean_marker_expression(
      oligo_average,
      stage_marker_sets[[stage]]
    )
}

stage_columns <- names(stage_marker_sets)

stage_score_table$oligodendroglial_stage <- apply(
  stage_score_table[
    ,
    stage_columns,
    drop = FALSE
  ],
  1,
  function(values) {
    if (all(is.na(values))) {
      return("OLIGO_UNRESOLVED")
    }
    
    stage_columns[
      which.max(values)
    ]
  }
)

cluster_to_stage <- setNames(
  stage_score_table$oligodendroglial_stage,
  stage_score_table$oligo_cluster
)

valid_oligo_clusters <- sort(
  unique(
    as.character(
      oligo$seurat_clusters
    )
  )
)

unmatched_oligo_clusters <- setdiff(
  stage_score_table$oligo_cluster,
  valid_oligo_clusters
)

if (length(unmatched_oligo_clusters) > 0L) {
  stop(
    "Oligodendroglial AverageExpression labels ",
    "do not match Seurat cluster identities: ",
    paste(
      unmatched_oligo_clusters,
      collapse = ", "
    )
  )
}

cluster_to_stage <- setNames(
  stage_score_table$oligodendroglial_stage,
  stage_score_table$oligo_cluster
)

oligo$oligodendroglial_stage <- unname(
  cluster_to_stage[
    as.character(
      oligo$seurat_clusters
    )
  ]
)

if (anyNA(oligo$oligodendroglial_stage)) {
  missing_oligo_clusters <- unique(
    as.character(
      oligo$seurat_clusters[
        is.na(
          oligo$oligodendroglial_stage
        )
      ]
    )
  )
  
  stop(
    "Some oligodendroglial cells could not be assigned ",
    "a maturation stage. Affected clusters: ",
    paste(
      missing_oligo_clusters,
      collapse = ", "
    )
  )
}

write_csv(
  stage_score_table,
  file.path(
    table_root,
    "GSE140393_OLIGODENDROGLIAL_STAGE_CLUSTER_SCORES.csv"
  )
)

###############################################################################
# 13. Oligodendroglial plots
###############################################################################

p_oligo_cluster <- DimPlot(
  oligo,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Oligodendroglial clusters")

p_oligo_stage <- DimPlot(
  oligo,
  reduction = "umap",
  group.by = "oligodendroglial_stage",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Oligodendroglial maturation stages")

p_oligo_donor <- DimPlot(
  oligo,
  reduction = "umap",
  group.by = "donor_id"
) +
  ggtitle("Oligodendroglial cells by donor")

save_plot(
  p_oligo_cluster +
    p_oligo_stage +
    p_oligo_donor,
  "GSE140393_OLIGODENDROGLIAL_UMAP",
  width = 18,
  height = 6
)

p_oligo_dotplot <- DotPlot(
  oligo,
  features = stage_features,
  group.by = "oligodendroglial_stage",
  dot.scale = 6
) +
  RotatedAxis() +
  ggtitle(
    "Oligodendroglial stage and P5 representative markers"
  )

save_plot(
  p_oligo_dotplot,
  "GSE140393_OLIGODENDROGLIAL_STAGE_DOTPLOT",
  width = 16,
  height = 6
)

###############################################################################
# 14. Donor-level cell counts
###############################################################################

donor_stage_counts <- oligo@meta.data %>%
  count(
    sample_id,
    donor_id,
    age,
    sex,
    age_group,
    original_paper_status,
    oligodendroglial_stage,
    name = "cell_count"
  ) %>%
  arrange(
    donor_id,
    oligodendroglial_stage
  )

donor_total_counts <- oligo@meta.data %>%
  count(
    sample_id,
    donor_id,
    age,
    sex,
    age_group,
    original_paper_status,
    name = "oligodendroglial_cell_count"
  ) %>%
  arrange(donor_id)

write_csv(
  donor_stage_counts,
  file.path(
    table_root,
    "GSE140393_DONOR_STAGE_CELL_COUNTS.csv"
  )
)

write_csv(
  donor_total_counts,
  file.path(
    table_root,
    "GSE140393_DONOR_OLIGODENDROGLIAL_CELL_COUNTS.csv"
  )
)

###############################################################################
# 15. Save cell metadata and objects
###############################################################################

oligo_metadata <- oligo@meta.data %>%
  rownames_to_column("cell_barcode")

write_csv(
  oligo_metadata,
  file.path(
    table_root,
    "GSE140393_OLIGODENDROGLIAL_CELL_METADATA.csv"
  )
)

saveRDS(
  merged,
  file.path(
    object_root,
    "GSE140393_TLE_ALL_CELLS_CLUSTERED.rds"
  ),
  compress = FALSE
)

saveRDS(
  oligo,
  file.path(
    object_root,
    "GSE140393_TLE_OLIGODENDROGLIA_ANNOTATED.rds"
  ),
  compress = FALSE
)

###############################################################################
# 16. Final summary
###############################################################################

summary_lines <- c(
  "# GSE140393 TLE QC and oligodendroglial annotation",
  "",
  paste0(
    "- Total QC-passed cells: ",
    ncol(merged)
  ),
  paste0(
    "- Oligodendroglial cells: ",
    ncol(oligo)
  ),
  paste0(
    "- Oligodendroglial clusters: ",
    paste(
      oligodendroglial_clusters,
      collapse = ", "
    )
  ),
  "",
  "## Donor-level oligodendroglial counts",
  "",
  paste0(
    "- ",
    donor_total_counts$donor_id,
    ": ",
    donor_total_counts$oligodendroglial_cell_count
  ),
  "",
  "## Important boundary",
  "",
  paste0(
    "- GSM4160577 remains in the full analysis but is ",
    "pre-labelled as low complexity based on the original paper."
  ),
  paste0(
    "- Its exclusion will be evaluated later as a ",
    "prespecified sensitivity analysis."
  ),
  "- P5 was not scored in this script."
)

writeLines(
  summary_lines,
  file.path(
    output_root,
    "GSE140393_QC_AND_OLIGO_SUMMARY.md"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    log_root,
    "sessionInfo.txt"
  )
)

message("")
message("===================================================")
message("GSE140393 QC and oligodendroglial annotation done")
message("===================================================")
message("")
message("QC-passed cells: ", ncol(merged))
message("Oligodendroglial cells: ", ncol(oligo))
message("")
message("Donor-level oligodendroglial counts:")

print(
  donor_total_counts,
  n = Inf
)

message("")
message("Output directory:")
message(output_root)
###############################################################################
# 35C-07
# P5 cell-number downsampling, UMI-depth downsampling and
# SEA-AD nuclear-detectability sensitivity analysis
#
# Primary comparison:
#   Adult TLE mature oligodendrocytes:
#     14431, 19619, 20188
#   versus SEA-AD external MTG reference:
#     H18.30.002, H19.30.001, H19.30.002
#
# No clustering, annotation, gene reselection or cell-level inference.
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
  "Matrix",
  "dplyr",
  "readr",
  "ggplot2"
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
  library(Matrix)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

###############################################################################
# 2. Paths and parameters
###############################################################################

project_root <- Sys.getenv("P5_WORK_ROOT")

work_root <- file.path(
  project_root,
  "35C_P5_CROSS_REFERENCE"
)

tle_object_path <- file.path(
  work_root,
  "06_GSE140393_FINAL_MANUAL_OLIGO_ANNOTATION",
  "objects",
  "GSE140393_MATURE_OLIGODENDROCYTES.rds"
)

seaad_root <- file.path(
  work_root,
  "07_SEAAD_R_MATRIX_EXPORT",
  "Mature_oligodendrocyte"
)

p5_path <- file.path(
  work_root,
  "08_FROZEN_P5_DONOR_SCORING",
  "00_input_audit",
  "FROZEN_P5_TOP100_AND_WEIGHTS.csv"
)

output_root <- file.path(
  work_root,
  "10_DOWNSAMPLING_AND_NUCLEAR_DETECTABILITY"
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

adult_tle_ids <- c(
  "14431",
  "19619",
  "20188"
)

reference_ids <- c(
  "H18.30.002",
  "H19.30.001",
  "H19.30.002"
)

n_replicates <- 100L
minimum_p5_genes <- 20L

required_paths <- c(
  tle_object_path,
  p5_path,
  file.path(seaad_root, "matrix.mtx.gz"),
  file.path(seaad_root, "features.tsv.gz"),
  file.path(seaad_root, "barcodes.tsv.gz"),
  file.path(seaad_root, "cell_metadata.csv.gz")
)

missing_paths <- required_paths[
  !file.exists(required_paths)
]

if (length(missing_paths) > 0L) {
  stop(
    "Missing required inputs:\n",
    paste(missing_paths, collapse = "\n")
  )
}

###############################################################################
# 3. General helpers
###############################################################################

canonicalize_gene_symbol <- function(symbols) {
  symbols <- toupper(trimws(as.character(symbols)))
  symbols[is.na(symbols)] <- ""
  
  sub(
    pattern = "\\.[0-9]+$",
    replacement = "",
    x = symbols,
    perl = TRUE
  )
}

collapse_rows_by_symbol <- function(counts, symbols) {
  symbols <- canonicalize_gene_symbol(symbols)
  
  keep <- nzchar(symbols)
  counts <- counts[keep, , drop = FALSE]
  symbols <- symbols[keep]
  
  unique_symbols <- unique(symbols)
  
  if (length(unique_symbols) == length(symbols)) {
    rownames(counts) <- symbols
    return(counts)
  }
  
  aggregation_matrix <- Matrix::sparseMatrix(
    i = match(symbols, unique_symbols),
    j = seq_along(symbols),
    x = 1,
    dims = c(
      length(unique_symbols),
      length(symbols)
    )
  )
  
  collapsed <- aggregation_matrix %*% counts
  
  rownames(collapsed) <- unique_symbols
  colnames(collapsed) <- colnames(counts)
  
  methods::as(collapsed, "dgCMatrix")
}

get_seurat_counts <- function(object) {
  DefaultAssay(object) <- "RNA"
  
  if (
    packageVersion("SeuratObject") >=
    package_version("5.0.0")
  ) {
    object <- tryCatch(
      JoinLayers(object),
      error = function(error) object
    )
    
    counts <- LayerData(
      object = object,
      assay = "RNA",
      layer = "counts"
    )
  } else {
    counts <- GetAssayData(
      object = object,
      assay = "RNA",
      slot = "counts"
    )
  }
  
  if (!inherits(counts, "sparseMatrix")) {
    counts <- Matrix::Matrix(
      counts,
      sparse = TRUE
    )
  }
  
  methods::as(counts, "dgCMatrix")
}

read_gzip_lines <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  
  readLines(connection, warn = FALSE)
}

read_gzip_table <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  
  read.delim(
    connection,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

read_matrix_market_gzip <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  
  matrix <- Matrix::readMM(connection)
  
  methods::as(matrix, "dgCMatrix")
}

aggregate_counts_by_donor <- function(
    counts,
    donor_ids,
    donors
) {
  donor_ids <- as.character(donor_ids)
  
  missing_donors <- setdiff(
    donors,
    unique(donor_ids)
  )
  
  if (length(missing_donors) > 0L) {
    stop(
      "Missing donors: ",
      paste(missing_donors, collapse = ", ")
    )
  }
  
  result <- matrix(
    0,
    nrow = nrow(counts),
    ncol = length(donors),
    dimnames = list(
      rownames(counts),
      donors
    )
  )
  
  for (donor_index in seq_along(donors)) {
    donor <- donors[donor_index]
    cell_indices <- which(donor_ids == donor)
    
    result[, donor_index] <- Matrix::rowSums(
      counts[, cell_indices, drop = FALSE]
    )
  }
  
  Matrix::Matrix(
    result,
    sparse = TRUE
  )
}

log_cpm <- function(counts) {
  counts <- as.matrix(counts)
  library_sizes <- colSums(counts)
  
  if (
    any(!is.finite(library_sizes)) ||
    any(library_sizes <= 0)
  ) {
    stop("Invalid pseudobulk library size.")
  }
  
  cpm <- sweep(
    counts,
    MARGIN = 2,
    STATS = library_sizes,
    FUN = "/"
  ) * 1e6
  
  log1p(cpm)
}

score_rank_program <- function(
    pseudobulk,
    universe_genes,
    program_genes,
    program_weights
) {
  program_genes <- intersect(
    program_genes,
    universe_genes
  )
  
  if (length(program_genes) < minimum_p5_genes) {
    stop(
      "Too few program genes for scoring: ",
      length(program_genes)
    )
  }
  
  weights <- program_weights[
    match(
      program_genes,
      names(program_weights)
    )
  ]
  
  weights <- weights / sum(weights)
  
  expression_matrix <- log_cpm(
    pseudobulk[
      universe_genes,
      ,
      drop = FALSE
    ]
  )
  
  rank_matrix <- apply(
    expression_matrix,
    MARGIN = 2,
    FUN = function(values) {
      (
        rank(
          values,
          ties.method = "average"
        ) - 1
      ) / (
        length(values) - 1
      )
    }
  )
  
  rank_matrix <- as.matrix(rank_matrix)
  
  rownames(rank_matrix) <- universe_genes
  colnames(rank_matrix) <- colnames(expression_matrix)
  
  weighted_scores <- colSums(
    sweep(
      rank_matrix[
        program_genes,
        ,
        drop = FALSE
      ],
      MARGIN = 1,
      STATS = weights,
      FUN = "*"
    )
  )
  
  equal_weight_scores <- colMeans(
    rank_matrix[
      program_genes,
      ,
      drop = FALSE
    ]
  )
  
  list(
    weighted = weighted_scores,
    equal_weight = equal_weight_scores,
    mapped_genes = length(program_genes)
  )
}

hedges_g <- function(tle_values, reference_values) {
  tle_values <- as.numeric(tle_values)
  reference_values <- as.numeric(reference_values)
  
  n_tle <- length(tle_values)
  n_reference <- length(reference_values)
  
  pooled_variance <- (
    (n_tle - 1) * stats::var(tle_values) +
      (n_reference - 1) * stats::var(reference_values)
  ) / (
    n_tle + n_reference - 2
  )
  
  if (
    !is.finite(pooled_variance) ||
    pooled_variance <= 0
  ) {
    return(NA_real_)
  }
  
  cohen_d <- (
    mean(tle_values) -
      mean(reference_values)
  ) / sqrt(pooled_variance)
  
  correction <- 1 - (
    3 /
      (
        4 * (n_tle + n_reference) - 9
      )
  )
  
  correction * cohen_d
}

exact_permutation_p <- function(
    tle_values,
    reference_values
) {
  all_values <- c(
    as.numeric(tle_values),
    as.numeric(reference_values)
  )
  
  n_tle <- length(tle_values)
  
  observed <- (
    mean(tle_values) -
      mean(reference_values)
  )
  
  combinations <- combn(
    seq_along(all_values),
    n_tle
  )
  
  permuted <- apply(
    combinations,
    2,
    function(tle_indices) {
      mean(all_values[tle_indices]) -
        mean(all_values[-tle_indices])
    }
  )
  
  mean(
    abs(permuted) >=
      abs(observed) -
      1e-12
  )
}

summarize_contrast <- function(
    scores,
    tle_ids,
    reference_ids
) {
  tle_values <- scores[tle_ids]
  reference_values <- scores[reference_ids]
  
  mean_difference <- (
    mean(tle_values) -
      mean(reference_values)
  )
  
  data.frame(
    tle_mean = mean(tle_values),
    reference_mean = mean(reference_values),
    mean_difference = mean_difference,
    hedges_g = hedges_g(
      tle_values,
      reference_values
    ),
    exact_permutation_p = exact_permutation_p(
      tle_values,
      reference_values
    ),
    direction = ifelse(
      mean_difference > 0,
      "TLE_HIGHER",
      ifelse(
        mean_difference < 0,
        "TLE_LOWER",
        "NO_DIFFERENCE"
      )
    ),
    complete_separation = (
      min(tle_values) >
        max(reference_values)
    ),
    stringsAsFactors = FALSE
  )
}

###############################################################################
# 4. Read TLE mature oligodendrocytes
###############################################################################

message("Reading GSE140393 mature oligodendrocytes.")

tle_object <- readRDS(tle_object_path)
tle_counts <- get_seurat_counts(tle_object)

tle_metadata <- tle_object@meta.data[
  colnames(tle_counts),
  ,
  drop = FALSE
]

if (!"donor_id" %in% colnames(tle_metadata)) {
  stop("TLE object lacks donor_id.")
}

tle_counts <- collapse_rows_by_symbol(
  counts = tle_counts,
  symbols = rownames(tle_counts)
)

tle_donor_ids <- as.character(
  tle_metadata$donor_id
)

tle_keep <- tle_donor_ids %in%
  adult_tle_ids

tle_counts <- tle_counts[
  ,
  tle_keep,
  drop = FALSE
]

tle_donor_ids <- tle_donor_ids[
  tle_keep
]

###############################################################################
# 5. Read SEA-AD mature oligodendrocytes
###############################################################################

message("Reading SEA-AD mature oligodendrocytes.")

seaad_counts <- read_matrix_market_gzip(
  file.path(
    seaad_root,
    "matrix.mtx.gz"
  )
)

seaad_features <- read_gzip_table(
  file.path(
    seaad_root,
    "features.tsv.gz"
  )
)

seaad_barcodes <- read_gzip_lines(
  file.path(
    seaad_root,
    "barcodes.tsv.gz"
  )
)

seaad_metadata <- suppressMessages(
  readr::read_csv(
    file.path(
      seaad_root,
      "cell_metadata.csv.gz"
    ),
    show_col_types = FALSE
  )
)

if (nrow(seaad_features) != nrow(seaad_counts)) {
  stop("SEA-AD feature count mismatch.")
}

if (length(seaad_barcodes) != ncol(seaad_counts)) {
  stop("SEA-AD barcode count mismatch.")
}

seaad_gene_symbols <- if (
  ncol(seaad_features) >= 2L
) {
  seaad_features[[2L]]
} else {
  seaad_features[[1L]]
}

rownames(seaad_counts) <- canonicalize_gene_symbol(
  seaad_gene_symbols
)

colnames(seaad_counts) <- seaad_barcodes

metadata_index <- match(
  seaad_barcodes,
  seaad_metadata$cell_barcode
)

if (anyNA(metadata_index)) {
  stop("Some SEA-AD barcodes were not found in metadata.")
}

seaad_metadata <- seaad_metadata[
  metadata_index,
  ,
  drop = FALSE
]

seaad_donor_ids <- as.character(
  seaad_metadata$donor_id_35c
)

seaad_keep <- seaad_donor_ids %in%
  reference_ids

seaad_counts <- seaad_counts[
  ,
  seaad_keep,
  drop = FALSE
]

seaad_donor_ids <- seaad_donor_ids[
  seaad_keep
]

seaad_counts <- collapse_rows_by_symbol(
  counts = seaad_counts,
  symbols = rownames(seaad_counts)
)

###############################################################################
# 6. Read frozen P5
###############################################################################

p5_table <- suppressMessages(
  readr::read_csv(
    p5_path,
    show_col_types = FALSE
  )
)

required_p5_columns <- c(
  "gene",
  "normalized_weight"
)

if (
  !all(
    required_p5_columns %in%
    colnames(p5_table)
  )
) {
  stop(
    "Frozen P5 table lacks required columns."
  )
}

p5_table <- p5_table %>%
  transmute(
    gene = canonicalize_gene_symbol(gene),
    normalized_weight =
      as.numeric(normalized_weight)
  ) %>%
  filter(
    nzchar(gene),
    is.finite(normalized_weight),
    normalized_weight > 0
  ) %>%
  group_by(gene) %>%
  summarise(
    normalized_weight =
      max(normalized_weight),
    .groups = "drop"
  )

p5_weight_vector <- p5_table$normalized_weight
names(p5_weight_vector) <- p5_table$gene

###############################################################################
# 7. Common genes and full pseudobulk
###############################################################################

common_genes <- intersect(
  rownames(tle_counts),
  rownames(seaad_counts)
)

tle_counts <- tle_counts[
  common_genes,
  ,
  drop = FALSE
]

seaad_counts <- seaad_counts[
  common_genes,
  ,
  drop = FALSE
]

tle_full_pseudobulk <- aggregate_counts_by_donor(
  counts = tle_counts,
  donor_ids = tle_donor_ids,
  donors = adult_tle_ids
)

reference_full_pseudobulk <- aggregate_counts_by_donor(
  counts = seaad_counts,
  donor_ids = seaad_donor_ids,
  donors = reference_ids
)

full_pseudobulk <- cbind(
  tle_full_pseudobulk,
  reference_full_pseudobulk
)

stable_universe <- rownames(
  full_pseudobulk
)[
  Matrix::rowSums(
    full_pseudobulk > 0
  ) == ncol(full_pseudobulk)
]

if (length(stable_universe) < 3000L) {
  stop(
    "Too few genes detected in all six donor pseudobulks: ",
    length(stable_universe)
  )
}

stable_p5_genes <- intersect(
  names(p5_weight_vector),
  stable_universe
)

if (length(stable_p5_genes) < minimum_p5_genes) {
  stop(
    "Too few stable P5 genes: ",
    length(stable_p5_genes)
  )
}

###############################################################################
# 8. Full-data reference result
###############################################################################

full_score_result <- score_rank_program(
  pseudobulk = full_pseudobulk,
  universe_genes = stable_universe,
  program_genes = stable_p5_genes,
  program_weights = p5_weight_vector
)

full_weighted_result <- summarize_contrast(
  scores = full_score_result$weighted,
  tle_ids = adult_tle_ids,
  reference_ids = reference_ids
)

full_equal_result <- summarize_contrast(
  scores = full_score_result$equal_weight,
  tle_ids = adult_tle_ids,
  reference_ids = reference_ids
)

full_results <- bind_rows(
  cbind(
    scoring_method = "FROZEN_WEIGHTED",
    full_weighted_result
  ),
  cbind(
    scoring_method = "EQUAL_WEIGHT",
    full_equal_result
  )
)

full_results$mapped_p5_genes <-
  full_score_result$mapped_genes

###############################################################################
# 9. Cell-number downsampling
###############################################################################

cell_counts <- c(
  table(
    factor(
      tle_donor_ids,
      levels = adult_tle_ids
    )
  ),
  table(
    factor(
      seaad_donor_ids,
      levels = reference_ids
    )
  )
)

target_cell_count <- min(
  as.integer(cell_counts)
)

message(
  "Cell-number downsampling target: ",
  target_cell_count,
  " cells per donor."
)

sample_cells_by_donor <- function(
    counts,
    donor_ids,
    donors,
    target_cells
) {
  result <- matrix(
    0,
    nrow = nrow(counts),
    ncol = length(donors),
    dimnames = list(
      rownames(counts),
      donors
    )
  )
  
  for (donor_index in seq_along(donors)) {
    donor <- donors[donor_index]
    available_indices <- which(
      donor_ids == donor
    )
    
    if (length(available_indices) < target_cells) {
      stop(
        donor,
        " has fewer cells than the downsampling target."
      )
    }
    
    selected_indices <- sample(
      available_indices,
      size = target_cells,
      replace = FALSE
    )
    
    result[, donor_index] <- Matrix::rowSums(
      counts[
        ,
        selected_indices,
        drop = FALSE
      ]
    )
  }
  
  Matrix::Matrix(
    result,
    sparse = TRUE
  )
}

cell_downsampling_rows <- vector(
  mode = "list",
  length = n_replicates * 2L
)

row_index <- 1L

set.seed(20260731)

for (replicate_index in seq_len(n_replicates)) {
  tle_sampled <- sample_cells_by_donor(
    counts = tle_counts,
    donor_ids = tle_donor_ids,
    donors = adult_tle_ids,
    target_cells = target_cell_count
  )
  
  reference_sampled <- sample_cells_by_donor(
    counts = seaad_counts,
    donor_ids = seaad_donor_ids,
    donors = reference_ids,
    target_cells = target_cell_count
  )
  
  sampled_pseudobulk <- cbind(
    tle_sampled,
    reference_sampled
  )
  
  sampled_score <- score_rank_program(
    pseudobulk = sampled_pseudobulk,
    universe_genes = stable_universe,
    program_genes = stable_p5_genes,
    program_weights = p5_weight_vector
  )
  
  weighted_result <- summarize_contrast(
    scores = sampled_score$weighted,
    tle_ids = adult_tle_ids,
    reference_ids = reference_ids
  )
  
  equal_result <- summarize_contrast(
    scores = sampled_score$equal_weight,
    tle_ids = adult_tle_ids,
    reference_ids = reference_ids
  )
  
  cell_downsampling_rows[[row_index]] <- cbind(
    replicate = replicate_index,
    scoring_method = "FROZEN_WEIGHTED",
    target_cells = target_cell_count,
    weighted_result
  )
  
  row_index <- row_index + 1L
  
  cell_downsampling_rows[[row_index]] <- cbind(
    replicate = replicate_index,
    scoring_method = "EQUAL_WEIGHT",
    target_cells = target_cell_count,
    equal_result
  )
  
  row_index <- row_index + 1L
}

cell_downsampling_results <- bind_rows(
  cell_downsampling_rows
)

###############################################################################
# 10. UMI-depth downsampling
###############################################################################

full_library_sizes <- colSums(
  full_pseudobulk
)

target_umi_depth <- floor(
  min(full_library_sizes)
)

message(
  "UMI downsampling target: ",
  target_umi_depth,
  " UMIs per donor pseudobulk."
)

downsample_count_vector <- function(
    count_vector,
    target_depth
) {
  count_vector <- as.numeric(count_vector)
  total_count <- sum(count_vector)
  
  if (total_count < target_depth) {
    stop(
      "Pseudobulk library is smaller than target depth."
    )
  }
  
  probabilities <- count_vector / total_count
  
  as.numeric(
    stats::rmultinom(
      n = 1,
      size = target_depth,
      prob = probabilities
    )[, 1L]
  )
}

umi_downsampling_rows <- vector(
  mode = "list",
  length = n_replicates * 2L
)

row_index <- 1L

set.seed(20260801)

for (replicate_index in seq_len(n_replicates)) {
  sampled_matrix <- matrix(
    0,
    nrow = nrow(full_pseudobulk),
    ncol = ncol(full_pseudobulk),
    dimnames = dimnames(full_pseudobulk)
  )
  
  for (donor_index in seq_len(ncol(full_pseudobulk))) {
    sampled_matrix[, donor_index] <-
      downsample_count_vector(
        count_vector =
          full_pseudobulk[, donor_index],
        target_depth = target_umi_depth
      )
  }
  
  sampled_matrix <- Matrix::Matrix(
    sampled_matrix,
    sparse = TRUE
  )
  
  sampled_score <- score_rank_program(
    pseudobulk = sampled_matrix,
    universe_genes = stable_universe,
    program_genes = stable_p5_genes,
    program_weights = p5_weight_vector
  )
  
  weighted_result <- summarize_contrast(
    scores = sampled_score$weighted,
    tle_ids = adult_tle_ids,
    reference_ids = reference_ids
  )
  
  equal_result <- summarize_contrast(
    scores = sampled_score$equal_weight,
    tle_ids = adult_tle_ids,
    reference_ids = reference_ids
  )
  
  umi_downsampling_rows[[row_index]] <- cbind(
    replicate = replicate_index,
    scoring_method = "FROZEN_WEIGHTED",
    target_umi_depth = target_umi_depth,
    weighted_result
  )
  
  row_index <- row_index + 1L
  
  umi_downsampling_rows[[row_index]] <- cbind(
    replicate = replicate_index,
    scoring_method = "EQUAL_WEIGHT",
    target_umi_depth = target_umi_depth,
    equal_result
  )
  
  row_index <- row_index + 1L
}

umi_downsampling_results <- bind_rows(
  umi_downsampling_rows
)

###############################################################################
# 11. SEA-AD nuclear-detectability audit
###############################################################################

p5_common_genes <- intersect(
  names(p5_weight_vector),
  rownames(seaad_counts)
)

reference_detection_rows <- list()

for (donor in reference_ids) {
  donor_indices <- which(
    seaad_donor_ids == donor
  )
  
  donor_detection <- Matrix::rowMeans(
    seaad_counts[
      p5_common_genes,
      donor_indices,
      drop = FALSE
    ] > 0
  )
  
  reference_detection_rows[[donor]] <- data.frame(
    gene = p5_common_genes,
    donor_id = donor,
    detection_fraction = as.numeric(
      donor_detection
    ),
    stringsAsFactors = FALSE
  )
}

reference_detection <- bind_rows(
  reference_detection_rows
)

reference_detection_wide <- reference_detection %>%
  group_by(gene) %>%
  summarise(
    minimum_detection =
      min(detection_fraction),
    median_detection =
      median(detection_fraction),
    donors_detected_1pct =
      sum(detection_fraction >= 0.01),
    donors_detected_5pct =
      sum(detection_fraction >= 0.05),
    donors_detected_10pct =
      sum(detection_fraction >= 0.10),
    .groups = "drop"
  )

reference_gene_sets <- list(
  REF_DETECT_1PCT_ALL3 =
    reference_detection_wide$gene[
      reference_detection_wide$
        donors_detected_1pct == 3L
    ],
  
  REF_DETECT_5PCT_ALL3 =
    reference_detection_wide$gene[
      reference_detection_wide$
        donors_detected_5pct == 3L
    ],
  
  REF_DETECT_10PCT_ALL3 =
    reference_detection_wide$gene[
      reference_detection_wide$
        donors_detected_10pct == 3L
    ],
  
  REF_DETECT_5PCT_AT_LEAST2 =
    reference_detection_wide$gene[
      reference_detection_wide$
        donors_detected_5pct >= 2L
    ]
)

detectability_result_rows <- list()

for (set_name in names(reference_gene_sets)) {
  selected_genes <- intersect(
    reference_gene_sets[[set_name]],
    stable_p5_genes
  )
  
  if (length(selected_genes) < minimum_p5_genes) {
    detectability_result_rows[[set_name]] <-
      data.frame(
        gene_set = set_name,
        mapped_genes = length(selected_genes),
        scoring_status = "INSUFFICIENT_GENES",
        tle_mean = NA_real_,
        reference_mean = NA_real_,
        mean_difference = NA_real_,
        hedges_g = NA_real_,
        exact_permutation_p = NA_real_,
        direction = NA_character_,
        complete_separation = NA,
        stringsAsFactors = FALSE
      )
    
    next
  }
  
  score_result <- score_rank_program(
    pseudobulk = full_pseudobulk,
    universe_genes = stable_universe,
    program_genes = selected_genes,
    program_weights = p5_weight_vector
  )
  
  contrast_result <- summarize_contrast(
    scores = score_result$weighted,
    tle_ids = adult_tle_ids,
    reference_ids = reference_ids
  )
  
  detectability_result_rows[[set_name]] <-
    cbind(
      gene_set = set_name,
      mapped_genes = length(selected_genes),
      scoring_status = "SCORED",
      contrast_result
    )
}

detectability_results <- bind_rows(
  detectability_result_rows
)

###############################################################################
# 12. Predefined gene-class sensitivity
###############################################################################

ribosomal_pattern <- paste0(
  "^(RPL|RPS|MRPL|MRPS)"
)

oxidative_pattern <- paste0(
  "^(ATP5|COX|NDUF|UQCR|SDH|CYC1$)"
)

gene_class_sets <- list(
  FULL_STABLE_P5 = stable_p5_genes,
  
  DROP_RIBOSOMAL = stable_p5_genes[
    !grepl(
      ribosomal_pattern,
      stable_p5_genes
    )
  ],
  
  DROP_OXPHOS = stable_p5_genes[
    !grepl(
      oxidative_pattern,
      stable_p5_genes
    )
  ],
  
  DROP_RIBOSOMAL_AND_OXPHOS =
    stable_p5_genes[
      !grepl(
        paste0(
          ribosomal_pattern,
          "|",
          oxidative_pattern
        ),
        stable_p5_genes
      )
    ]
)

gene_class_result_rows <- list()

for (set_name in names(gene_class_sets)) {
  selected_genes <- gene_class_sets[[set_name]]
  
  score_result <- score_rank_program(
    pseudobulk = full_pseudobulk,
    universe_genes = stable_universe,
    program_genes = selected_genes,
    program_weights = p5_weight_vector
  )
  
  contrast_result <- summarize_contrast(
    scores = score_result$weighted,
    tle_ids = adult_tle_ids,
    reference_ids = reference_ids
  )
  
  gene_class_result_rows[[set_name]] <-
    cbind(
      gene_set = set_name,
      mapped_genes = length(selected_genes),
      contrast_result
    )
}

gene_class_results <- bind_rows(
  gene_class_result_rows
)

###############################################################################
# 13. Summaries
###############################################################################

downsampling_summary <- bind_rows(
  cell_downsampling_results %>%
    group_by(scoring_method) %>%
    summarise(
      analysis = "CELL_NUMBER_DOWNSAMPLING",
      replicates = n(),
      target_value = first(target_cells),
      median_difference =
        median(mean_difference),
      minimum_difference =
        min(mean_difference),
      maximum_difference =
        max(mean_difference),
      proportion_TLE_higher =
        mean(mean_difference > 0),
      proportion_complete_separation =
        mean(complete_separation),
      median_hedges_g =
        median(hedges_g, na.rm = TRUE),
      .groups = "drop"
    ),
  
  umi_downsampling_results %>%
    group_by(scoring_method) %>%
    summarise(
      analysis = "UMI_DEPTH_DOWNSAMPLING",
      replicates = n(),
      target_value = first(target_umi_depth),
      median_difference =
        median(mean_difference),
      minimum_difference =
        min(mean_difference),
      maximum_difference =
        max(mean_difference),
      proportion_TLE_higher =
        mean(mean_difference > 0),
      proportion_complete_separation =
        mean(complete_separation),
      median_hedges_g =
        median(hedges_g, na.rm = TRUE),
      .groups = "drop"
    )
)

###############################################################################
# 14. Write outputs
###############################################################################

write_csv(
  full_results,
  file.path(
    table_root,
    "FULL_DATA_STABLE_GENE_RESULTS.csv"
  )
)

write_csv(
  cell_downsampling_results,
  file.path(
    table_root,
    "CELL_NUMBER_DOWNSAMPLING_RESULTS.csv"
  )
)

write_csv(
  umi_downsampling_results,
  file.path(
    table_root,
    "UMI_DEPTH_DOWNSAMPLING_RESULTS.csv"
  )
)

write_csv(
  downsampling_summary,
  file.path(
    table_root,
    "DOWNSAMPLING_SUMMARY.csv"
  )
)

write_csv(
  reference_detection,
  file.path(
    table_root,
    "SEAAD_P5_GENE_DETECTION_BY_DONOR.csv"
  )
)

write_csv(
  reference_detection_wide,
  file.path(
    table_root,
    "SEAAD_P5_REFERENCE_DETECTABILITY_SUMMARY.csv"
  )
)

write_csv(
  detectability_results,
  file.path(
    table_root,
    "REFERENCE_DETECTABLE_P5_RESULTS.csv"
  )
)

write_csv(
  gene_class_results,
  file.path(
    table_root,
    "P5_GENE_CLASS_SENSITIVITY_RESULTS.csv"
  )
)

###############################################################################
# 15. Figures
###############################################################################

cell_plot <- ggplot(
  cell_downsampling_results,
  aes(
    x = mean_difference,
    fill = scoring_method
  )
) +
  geom_histogram(
    bins = 30,
    alpha = 0.55,
    position = "identity"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  labs(
    title = "P5 cell-number downsampling",
    subtitle = paste0(
      target_cell_count,
      " mature oligodendrocytes per donor"
    ),
    x = "TLE minus external-reference score",
    y = "Replicates",
    fill = "Scoring"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    figure_root,
    "P5_CELL_NUMBER_DOWNSAMPLING.png"
  ),
  plot = cell_plot,
  width = 8,
  height = 5.5,
  dpi = 320,
  bg = "white"
)

umi_plot <- ggplot(
  umi_downsampling_results,
  aes(
    x = mean_difference,
    fill = scoring_method
  )
) +
  geom_histogram(
    bins = 30,
    alpha = 0.55,
    position = "identity"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  labs(
    title = "P5 pseudobulk UMI-depth downsampling",
    subtitle = paste0(
      target_umi_depth,
      " UMIs per donor"
    ),
    x = "TLE minus external-reference score",
    y = "Replicates",
    fill = "Scoring"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    figure_root,
    "P5_UMI_DEPTH_DOWNSAMPLING.png"
  ),
  plot = umi_plot,
  width = 8,
  height = 5.5,
  dpi = 320,
  bg = "white"
)

###############################################################################
# 16. Final interpretation summary
###############################################################################

weighted_cell_summary <- downsampling_summary %>%
  filter(
    analysis ==
      "CELL_NUMBER_DOWNSAMPLING",
    scoring_method ==
      "FROZEN_WEIGHTED"
  )

weighted_umi_summary <- downsampling_summary %>%
  filter(
    analysis ==
      "UMI_DEPTH_DOWNSAMPLING",
    scoring_method ==
      "FROZEN_WEIGHTED"
  )

summary_lines <- c(
  "# 35C-07 Downsampling and nuclear detectability",
  "",
  paste0(
    "- Adult TLE donors: ",
    paste(adult_tle_ids, collapse = "|")
  ),
  paste0(
    "- SEA-AD reference donors: ",
    paste(reference_ids, collapse = "|")
  ),
  paste0(
    "- Common genes: ",
    length(common_genes)
  ),
  paste0(
    "- Genes detected in all six pseudobulks: ",
    length(stable_universe)
  ),
  paste0(
    "- Stable mapped P5 genes: ",
    length(stable_p5_genes)
  ),
  "",
  "## Cell-number downsampling",
  "",
  paste0(
    "- Target cells per donor: ",
    target_cell_count
  ),
  paste0(
    "- Replicates with TLE higher: ",
    weighted_cell_summary$
      proportion_TLE_higher
  ),
  paste0(
    "- Replicates with complete separation: ",
    weighted_cell_summary$
      proportion_complete_separation
  ),
  "",
  "## UMI-depth downsampling",
  "",
  paste0(
    "- Target UMIs per donor: ",
    target_umi_depth
  ),
  paste0(
    "- Replicates with TLE higher: ",
    weighted_umi_summary$
      proportion_TLE_higher
  ),
  paste0(
    "- Replicates with complete separation: ",
    weighted_umi_summary$
      proportion_complete_separation
  ),
  "",
  "## Interpretation boundary",
  "",
  paste0(
    "- These analyses test cell-number, library-depth and ",
    "reference nuclear-detectability sensitivity."
  ),
  paste0(
    "- They do not remove complete confounding between ",
    "disease, dataset and whole-cell versus nuclear modality."
  )
)

writeLines(
  summary_lines,
  file.path(
    output_root,
    "35C_07_DOWNSAMPLING_SUMMARY.md"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    log_root,
    "sessionInfo.txt"
  )
)

###############################################################################
# 17. Console output
###############################################################################

message("")
message("===================================================")
message("Downsampling and nuclear detectability completed")
message("===================================================")
message("")
message("Target cells per donor: ", target_cell_count)
message("Target UMI depth: ", target_umi_depth)
message("Stable universe genes: ", length(stable_universe))
message("Stable mapped P5 genes: ", length(stable_p5_genes))
message("")
message("Downsampling summary:")

print(
  as.data.frame(downsampling_summary),
  row.names = FALSE
)

message("")
message("Reference-detectable P5 results:")

print(
  as.data.frame(detectability_results),
  row.names = FALSE
)

message("")
message("Gene-class sensitivity:")

print(
  as.data.frame(gene_class_results),
  row.names = FALSE
)

message("")
message("Output directory:")
message(output_root)
###############################################################################
# 35C-05
# Frozen P5 donor-level scoring:
# GSE140393 TLE scRNA-seq versus SEA-AD external MTG snRNA-seq reference
#
# Primary comparison:
#   Adult TLE mature oligodendrocytes:
#     14431, 19619, 20188
#   versus SEA-AD mature oligodendrocytes:
#     H18.30.002, H19.30.001, H19.30.002
#
# Evidence class:
#   EXTERNAL_CROSS_STUDY_CASE_REFERENCE_COMPARISON
#
# Important:
#   - Donor is the statistical unit.
#   - Frozen P5 genes and weights are not reselected or reweighted.
#   - scRNA-seq and snRNA-seq matrices are processed separately.
#   - No cell-level pseudoreplication is used.
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
  "tibble",
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
  library(tibble)
  library(readr)
  library(ggplot2)
})

###############################################################################
# 2. Paths
###############################################################################

project_root <- Sys.getenv("P5_WORK_ROOT")

work_root <- file.path(
  project_root,
  "35C_P5_CROSS_REFERENCE"
)

frozen_root <- file.path(
  project_root,
  "99_frozen_releases"
)

tle_object_root <- file.path(
  work_root,
  "06_GSE140393_FINAL_MANUAL_OLIGO_ANNOTATION",
  "objects"
)

tle_mature_path <- file.path(
  tle_object_root,
  "GSE140393_MATURE_OLIGODENDROCYTES.rds"
)

tle_opc_path <- file.path(
  tle_object_root,
  "GSE140393_OPC.rds"
)

seaad_export_root <- file.path(
  work_root,
  "07_SEAAD_R_MATRIX_EXPORT"
)

seaad_mature_root <- file.path(
  seaad_export_root,
  "Mature_oligodendrocyte"
)

seaad_opc_root <- file.path(
  seaad_export_root,
  "OPC"
)

output_root <- file.path(
  work_root,
  "08_FROZEN_P5_DONOR_SCORING"
)

input_audit_root <- file.path(
  output_root,
  "00_input_audit"
)

score_root <- file.path(
  output_root,
  "01_donor_scores"
)

statistics_root <- file.path(
  output_root,
  "02_statistics"
)

figure_root <- file.path(
  output_root,
  "03_figures"
)

log_root <- file.path(
  output_root,
  "04_logs"
)

temporary_input_root <- file.path(
  output_root,
  "temporary_frozen_input"
)

for (directory in c(
  output_root,
  input_audit_root,
  score_root,
  statistics_root,
  figure_root,
  log_root,
  temporary_input_root
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

required_paths <- c(
  tle_mature_path,
  tle_opc_path,
  file.path(
    seaad_mature_root,
    "matrix.mtx.gz"
  ),
  file.path(
    seaad_mature_root,
    "features.tsv.gz"
  ),
  file.path(
    seaad_mature_root,
    "barcodes.tsv.gz"
  ),
  file.path(
    seaad_mature_root,
    "cell_metadata.csv.gz"
  ),
  file.path(
    seaad_opc_root,
    "matrix.mtx.gz"
  ),
  file.path(
    seaad_opc_root,
    "features.tsv.gz"
  ),
  file.path(
    seaad_opc_root,
    "barcodes.tsv.gz"
  ),
  file.path(
    seaad_opc_root,
    "cell_metadata.csv.gz"
  )
)

missing_paths <- required_paths[
  !file.exists(required_paths)
]

if (length(missing_paths) > 0L) {
  stop(
    "Missing required input files:\n",
    paste(
      missing_paths,
      collapse = "\n"
    )
  )
}

###############################################################################
# 3. General helper functions
###############################################################################

canonicalize_gene_symbol <- function(symbols) {
  symbols <- toupper(
    trimws(
      as.character(symbols)
    )
  )
  
  symbols[is.na(symbols)] <- ""
  
  # Read10X(unique.features = TRUE) may append .1, .2, ...
  symbols <- sub(
    pattern = "\\.[0-9]+$",
    replacement = "",
    x = symbols,
    perl = TRUE
  )
  
  symbols
}

collapse_rows_by_symbol <- function(
    counts,
    symbols
) {
  symbols <- canonicalize_gene_symbol(
    symbols
  )
  
  keep <- nzchar(symbols)
  
  counts <- counts[
    keep,
    ,
    drop = FALSE
  ]
  
  symbols <- symbols[keep]
  
  unique_symbols <- unique(symbols)
  
  if (length(unique_symbols) == length(symbols)) {
    rownames(counts) <- symbols
    
    return(counts)
  }
  
  symbol_index <- match(
    symbols,
    unique_symbols
  )
  
  aggregation_matrix <- Matrix::sparseMatrix(
    i = symbol_index,
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
  
  collapsed
}

get_seurat_counts <- function(object) {
  DefaultAssay(object) <- "RNA"
  
  if (
    packageVersion("SeuratObject") >=
    package_version("5.0.0")
  ) {
    object <- tryCatch(
      JoinLayers(object),
      error = function(error) {
        message(
          "JoinLayers skipped safely: ",
          conditionMessage(error)
        )
        
        object
      }
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
  
  counts
}

read_gzip_lines <- function(path) {
  connection <- gzfile(
    path,
    open = "rt"
  )
  
  on.exit(
    close(connection),
    add = TRUE
  )
  
  readLines(
    connection,
    warn = FALSE
  )
}

read_gzip_table <- function(
    path,
    header = FALSE,
    sep = "\t"
) {
  connection <- gzfile(
    path,
    open = "rt"
  )
  
  on.exit(
    close(connection),
    add = TRUE
  )
  
  read.delim(
    connection,
    header = header,
    sep = sep,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

read_matrix_market_gzip <- function(path) {
  connection <- gzfile(
    path,
    open = "rt"
  )
  
  on.exit(
    close(connection),
    add = TRUE
  )
  
  matrix <- Matrix::readMM(
    connection
  )
  
  methods::as(
    matrix,
    "CsparseMatrix"
  )
}

aggregate_counts_by_donor <- function(
    counts,
    donor_ids
) {
  donor_ids <- as.character(
    donor_ids
  )
  
  if (length(donor_ids) != ncol(counts)) {
    stop(
      "Donor vector length does not match matrix columns."
    )
  }
  
  if (anyNA(donor_ids) || any(!nzchar(donor_ids))) {
    stop(
      "Missing donor IDs were detected."
    )
  }
  
  donors <- sort(
    unique(donor_ids)
  )
  
  donor_factor <- factor(
    donor_ids,
    levels = donors
  )
  
  design <- Matrix::sparseMatrix(
    i = seq_along(donor_factor),
    j = as.integer(donor_factor),
    x = 1,
    dims = c(
      length(donor_factor),
      length(donors)
    ),
    dimnames = list(
      colnames(counts),
      donors
    )
  )
  
  pseudobulk <- counts %*% design
  
  colnames(pseudobulk) <- donors
  
  donor_cell_counts <- table(
    factor(
      donor_ids,
      levels = donors
    )
  )
  
  list(
    counts = pseudobulk,
    donors = donors,
    cell_counts = setNames(
      as.integer(donor_cell_counts),
      donors
    )
  )
}

log_cpm <- function(pseudobulk_counts) {
  dense_counts <- as.matrix(
    pseudobulk_counts
  )
  
  library_sizes <- colSums(
    dense_counts
  )
  
  if (any(!is.finite(library_sizes))) {
    stop(
      "Non-finite pseudobulk library size detected."
    )
  }
  
  if (any(library_sizes <= 0)) {
    stop(
      "Zero pseudobulk library size detected."
    )
  }
  
  cpm <- sweep(
    dense_counts,
    MARGIN = 2,
    STATS = library_sizes,
    FUN = "/"
  ) * 1e6
  
  log1p(cpm)
}

###############################################################################
# 4. Locate and read frozen P5
###############################################################################

first_existing_column <- function(
    column_names,
    candidates
) {
  result <- intersect(
    candidates,
    column_names
  )
  
  if (length(result) == 0L) {
    return(NULL)
  }
  
  result[1L]
}

read_candidate_p5_file <- function(path) {
  data <- suppressMessages(
    readr::read_csv(
      path,
      show_col_types = FALSE
    )
  )
  
  gene_column <- first_existing_column(
    colnames(data),
    c(
      "gene",
      "gene_symbol",
      "symbol",
      "human_gene",
      "feature"
    )
  )
  
  weight_column <- first_existing_column(
    colnames(data),
    c(
      "meta_program_score",
      "weight",
      "program_weight",
      "score"
    )
  )
  
  program_column <- first_existing_column(
    colnames(data),
    c(
      "program",
      "meta_program",
      "program_name"
    )
  )
  
  if (
    is.null(gene_column) ||
    is.null(weight_column)
  ) {
    return(NULL)
  }
  
  if (!is.null(program_column)) {
    program_values <- as.character(
      data[[program_column]]
    )
    
    p5_rows <- toupper(
      trimws(program_values)
    ) == "OLIGODENDROCYTES__P5"
    
    if (any(p5_rows)) {
      data <- data[
        p5_rows,
        ,
        drop = FALSE
      ]
    }
  }
  
  result <- data.frame(
    gene = canonicalize_gene_symbol(
      data[[gene_column]]
    ),
    weight = suppressWarnings(
      as.numeric(
        data[[weight_column]]
      )
    ),
    stringsAsFactors = FALSE
  ) %>%
    filter(
      nzchar(gene),
      is.finite(weight),
      weight >= 0
    ) %>%
    group_by(gene) %>%
    summarise(
      weight = max(weight),
      .groups = "drop"
    )
  
  if (nrow(result) < 50L) {
    return(NULL)
  }
  
  list(
    path = path,
    data = result
  )
}

locate_frozen_p5 <- function() {
  csv_candidates <- character(0)
  
  if (dir.exists(frozen_root)) {
    csv_candidates <- list.files(
      frozen_root,
      pattern = "34A5C_FROZEN_HUMAN_TOP100\\.csv$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  }
  
  if (length(csv_candidates) == 0L) {
    csv_candidates <- list.files(
      project_root,
      pattern = "34A5C_FROZEN_HUMAN_TOP100\\.csv$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  }
  
  if (length(csv_candidates) > 0L) {
    csv_candidates <- csv_candidates[
      order(
        file.info(csv_candidates)$mtime,
        decreasing = TRUE
      )
    ]
    
    for (candidate in csv_candidates) {
      result <- tryCatch(
        read_candidate_p5_file(candidate),
        error = function(error) {
          NULL
        }
      )
      
      if (!is.null(result)) {
        return(result)
      }
    }
  }
  
  zip_candidates <- character(0)
  
  if (dir.exists(frozen_root)) {
    zip_candidates <- list.files(
      frozen_root,
      pattern = "34A5C.*\\.zip$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  }
  
  if (length(zip_candidates) == 0L) {
    stop(
      "Frozen P5 Top100 CSV or 34A5C ZIP was not found."
    )
  }
  
  zip_candidates <- zip_candidates[
    order(
      file.info(zip_candidates)$mtime,
      decreasing = TRUE
    )
  ]
  
  expected_basename <- tolower(
    "34A5C_FROZEN_HUMAN_TOP100.csv"
  )
  
  for (zip_path in zip_candidates) {
    zip_listing <- tryCatch(
      utils::unzip(
        zip_path,
        list = TRUE
      ),
      error = function(error) {
        NULL
      }
    )
    
    if (is.null(zip_listing)) {
      next
    }
    
    matching_entries <- zip_listing$Name[
      tolower(
        basename(zip_listing$Name)
      ) == expected_basename
    ]
    
    if (length(matching_entries) == 0L) {
      next
    }
    
    entry <- matching_entries[1L]
    
    utils::unzip(
      zipfile = zip_path,
      files = entry,
      exdir = temporary_input_root,
      overwrite = TRUE
    )
    
    extracted_path <- file.path(
      temporary_input_root,
      entry
    )
    
    result <- tryCatch(
      read_candidate_p5_file(
        extracted_path
      ),
      error = function(error) {
        NULL
      }
    )
    
    if (!is.null(result)) {
      result$source_zip <- zip_path
      
      return(result)
    }
  }
  
  stop(
    "A valid frozen P5 Top100 table could not be identified."
  )
}

p5_result <- locate_frozen_p5()

p5_table <- p5_result$data

if (nrow(p5_table) != 100L) {
  warning(
    "Frozen P5 table contains ",
    nrow(p5_table),
    " unique genes rather than exactly 100."
  )
}

if (sum(p5_table$weight) <= 0) {
  stop(
    "Frozen P5 weights sum to zero."
  )
}

p5_table <- p5_table %>%
  mutate(
    normalized_weight =
      weight / sum(weight)
  ) %>%
  arrange(
    desc(normalized_weight)
  )

write_csv(
  p5_table,
  file.path(
    input_audit_root,
    "FROZEN_P5_TOP100_AND_WEIGHTS.csv"
  )
)

writeLines(
  c(
    paste0(
      "P5 source file: ",
      p5_result$path
    ),
    paste0(
      "P5 unique genes: ",
      nrow(p5_table)
    )
  ),
  file.path(
    input_audit_root,
    "FROZEN_P5_SOURCE.txt"
  )
)

message("")
message("Frozen P5 source:")
message(p5_result$path)
message("Frozen P5 genes: ", nrow(p5_table))

###############################################################################
# 5. Read GSE140393 TLE compartments
###############################################################################

read_tle_compartment <- function(
    object_path,
    compartment_name
) {
  message("")
  message(
    "Reading GSE140393 ",
    compartment_name,
    " object."
  )
  
  object <- readRDS(
    object_path
  )
  
  counts <- get_seurat_counts(
    object
  )
  
  metadata <- object@meta.data[
    colnames(counts),
    ,
    drop = FALSE
  ]
  
  if (!"donor_id" %in% colnames(metadata)) {
    stop(
      compartment_name,
      " object lacks donor_id."
    )
  }
  
  counts <- collapse_rows_by_symbol(
    counts = counts,
    symbols = rownames(counts)
  )
  
  list(
    counts = counts,
    donor_ids = as.character(
      metadata$donor_id
    ),
    metadata = metadata,
    compartment = compartment_name,
    dataset = "GSE140393_TLE"
  )
}

tle_mature <- read_tle_compartment(
  object_path = tle_mature_path,
  compartment_name = "Mature_oligodendrocyte"
)

tle_opc <- read_tle_compartment(
  object_path = tle_opc_path,
  compartment_name = "OPC"
)

###############################################################################
# 6. Read SEA-AD exported matrices
###############################################################################

read_seaad_compartment <- function(
    directory,
    compartment_name
) {
  message("")
  message(
    "Reading SEA-AD ",
    compartment_name,
    " export."
  )
  
  matrix_path <- file.path(
    directory,
    "matrix.mtx.gz"
  )
  
  feature_path <- file.path(
    directory,
    "features.tsv.gz"
  )
  
  barcode_path <- file.path(
    directory,
    "barcodes.tsv.gz"
  )
  
  metadata_path <- file.path(
    directory,
    "cell_metadata.csv.gz"
  )
  
  counts <- read_matrix_market_gzip(
    matrix_path
  )
  
  features <- read_gzip_table(
    feature_path,
    header = FALSE,
    sep = "\t"
  )
  
  barcodes <- read_gzip_lines(
    barcode_path
  )
  
  metadata <- suppressMessages(
    readr::read_csv(
      metadata_path,
      show_col_types = FALSE
    )
  )
  
  if (nrow(features) != nrow(counts)) {
    stop(
      compartment_name,
      ": feature count does not match matrix rows."
    )
  }
  
  if (length(barcodes) != ncol(counts)) {
    stop(
      compartment_name,
      ": barcode count does not match matrix columns."
    )
  }
  
  if (ncol(features) >= 2L) {
    gene_symbols <- features[[2L]]
  } else {
    gene_symbols <- features[[1L]]
  }
  
  rownames(counts) <- canonicalize_gene_symbol(
    gene_symbols
  )
  
  colnames(counts) <- barcodes
  
  if (
    !"cell_barcode" %in%
    colnames(metadata)
  ) {
    stop(
      compartment_name,
      ": cell_metadata lacks cell_barcode."
    )
  }
  
  if (
    !"donor_id_35c" %in%
    colnames(metadata)
  ) {
    stop(
      compartment_name,
      ": cell_metadata lacks donor_id_35c."
    )
  }
  
  metadata_index <- match(
    barcodes,
    metadata$cell_barcode
  )
  
  if (anyNA(metadata_index)) {
    stop(
      compartment_name,
      ": some barcodes were not found in metadata."
    )
  }
  
  metadata <- metadata[
    metadata_index,
    ,
    drop = FALSE
  ]
  
  counts <- collapse_rows_by_symbol(
    counts = counts,
    symbols = rownames(counts)
  )
  
  list(
    counts = counts,
    donor_ids = as.character(
      metadata$donor_id_35c
    ),
    metadata = metadata,
    compartment = compartment_name,
    dataset = "SEAAD_EXTERNAL_REFERENCE"
  )
}

seaad_mature <- read_seaad_compartment(
  directory = seaad_mature_root,
  compartment_name = "Mature_oligodendrocyte"
)

seaad_opc <- read_seaad_compartment(
  directory = seaad_opc_root,
  compartment_name = "OPC"
)

###############################################################################
# 7. Score one compartment
###############################################################################

score_compartment <- function(
    tle_input,
    reference_input,
    p5_weights,
    compartment_name,
    minimum_cells = 30L,
    background_replicates = 500L
) {
  tle_aggregate <- aggregate_counts_by_donor(
    counts = tle_input$counts,
    donor_ids = tle_input$donor_ids
  )
  
  reference_aggregate <- aggregate_counts_by_donor(
    counts = reference_input$counts,
    donor_ids = reference_input$donor_ids
  )
  
  tle_eligible <- names(
    tle_aggregate$cell_counts
  )[
    tle_aggregate$cell_counts >= minimum_cells
  ]
  
  reference_eligible <- names(
    reference_aggregate$cell_counts
  )[
    reference_aggregate$cell_counts >= minimum_cells
  ]
  
  if (length(tle_eligible) < 2L) {
    stop(
      compartment_name,
      ": fewer than two eligible TLE donors."
    )
  }
  
  if (length(reference_eligible) < 2L) {
    stop(
      compartment_name,
      ": fewer than two eligible reference donors."
    )
  }
  
  common_genes <- intersect(
    rownames(tle_aggregate$counts),
    rownames(reference_aggregate$counts)
  )
  
  if (length(common_genes) < 5000L) {
    stop(
      compartment_name,
      ": too few common genes between datasets."
    )
  }
  
  tle_pseudobulk <- tle_aggregate$counts[
    common_genes,
    tle_eligible,
    drop = FALSE
  ]
  
  reference_pseudobulk <- reference_aggregate$counts[
    common_genes,
    reference_eligible,
    drop = FALSE
  ]
  
  expressed_in_tle <- Matrix::rowSums(
    tle_pseudobulk
  ) > 0
  
  expressed_in_reference <- Matrix::rowSums(
    reference_pseudobulk
  ) > 0
  
  universe_genes <- common_genes[
    expressed_in_tle &
      expressed_in_reference
  ]
  
  if (length(universe_genes) < 3000L) {
    stop(
      compartment_name,
      ": too few jointly expressed common genes."
    )
  }
  
  tle_logcpm <- log_cpm(
    tle_pseudobulk[
      universe_genes,
      ,
      drop = FALSE
    ]
  )
  
  reference_logcpm <- log_cpm(
    reference_pseudobulk[
      universe_genes,
      ,
      drop = FALSE
    ]
  )
  
  combined_logcpm <- cbind(
    tle_logcpm,
    reference_logcpm
  )
  
  p5_genes <- intersect(
    p5_weights$gene,
    universe_genes
  )
  
  if (length(p5_genes) < 50L) {
    stop(
      compartment_name,
      ": fewer than 50 frozen P5 genes are mappable."
    )
  }
  
  p5_weight_vector <- p5_weights$normalized_weight[
    match(
      p5_genes,
      p5_weights$gene
    )
  ]
  
  p5_weight_vector <- (
    p5_weight_vector /
      sum(p5_weight_vector)
  )
  
  rank_percentiles <- apply(
    combined_logcpm,
    MARGIN = 2,
    FUN = function(values) {
      if (length(values) <= 1L) {
        return(
          rep(
            NA_real_,
            length(values)
          )
        )
      }
      
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
  
  rank_percentiles <- as.matrix(
    rank_percentiles
  )
  
  rownames(rank_percentiles) <-
    rownames(combined_logcpm)
  
  colnames(rank_percentiles) <-
    colnames(combined_logcpm)
  
  raw_p5_rank_score <- colSums(
    sweep(
      rank_percentiles[
        p5_genes,
        ,
        drop = FALSE
      ],
      MARGIN = 1,
      STATS = p5_weight_vector,
      FUN = "*"
    )
  )
  
  mean_expression <- rowMeans(
    combined_logcpm,
    na.rm = TRUE
  )
  
  expression_rank <- rank(
    mean_expression,
    ties.method = "average"
  )
  
  expression_bin <- ceiling(
    expression_rank /
      length(expression_rank) *
      10
  )
  
  expression_bin <- pmin(
    10L,
    pmax(
      1L,
      expression_bin
    )
  )
  
  names(expression_bin) <-
    rownames(combined_logcpm)
  
  p5_bins <- expression_bin[
    p5_genes
  ]
  
  non_p5_genes <- setdiff(
    universe_genes,
    p5_genes
  )
  
  background_scores <- matrix(
    NA_real_,
    nrow = ncol(rank_percentiles),
    ncol = background_replicates,
    dimnames = list(
      colnames(rank_percentiles),
      paste0(
        "background_",
        seq_len(background_replicates)
      )
    )
  )
  
  set.seed(20260731)
  
  for (
    replicate_index in
    seq_len(background_replicates)
  ) {
    sampled_background <- character(
      length(p5_genes)
    )
    
    for (
      bin_value in
      sort(unique(p5_bins))
    ) {
      p5_indices <- which(
        p5_bins == bin_value
      )
      
      candidate_pool <- non_p5_genes[
        expression_bin[
          non_p5_genes
        ] == bin_value
      ]
      
      if (length(candidate_pool) == 0L) {
        candidate_pool <- non_p5_genes
      }
      
      sampled_background[
        p5_indices
      ] <- sample(
        candidate_pool,
        size = length(p5_indices),
        replace = (
          length(candidate_pool) <
            length(p5_indices)
        )
      )
    }
    
    background_scores[
      ,
      replicate_index
    ] <- colSums(
      sweep(
        rank_percentiles[
          sampled_background,
          ,
          drop = FALSE
        ],
        MARGIN = 1,
        STATS = p5_weight_vector,
        FUN = "*"
      )
    )
  }
  
  background_mean <- rowMeans(
    background_scores,
    na.rm = TRUE
  )
  
  background_sd <- apply(
    background_scores,
    MARGIN = 1,
    FUN = stats::sd,
    na.rm = TRUE
  )
  
  background_sd[
    !is.finite(background_sd) |
      background_sd <= 0
  ] <- NA_real_
  
  calibrated_rank_score <- (
    raw_p5_rank_score -
      background_mean
  )
  
  p5_nes <- calibrated_rank_score /
    background_sd
  
  all_donors <- colnames(
    combined_logcpm
  )
  
  tle_cell_counts <- tle_aggregate$cell_counts
  reference_cell_counts <-
    reference_aggregate$cell_counts
  
  donor_group <- ifelse(
    all_donors %in% tle_eligible,
    "TLE",
    "External_reference"
  )
  
  donor_dataset <- ifelse(
    donor_group == "TLE",
    "GSE140393",
    "SEA_AD"
  )
  
  donor_cell_count <- ifelse(
    donor_group == "TLE",
    tle_cell_counts[
      all_donors
    ],
    reference_cell_counts[
      all_donors
    ]
  )
  
  donor_scores <- data.frame(
    compartment = compartment_name,
    donor_id = all_donors,
    group = donor_group,
    dataset = donor_dataset,
    cell_count = as.integer(
      donor_cell_count
    ),
    common_gene_universe =
      length(universe_genes),
    mapped_p5_genes =
      length(p5_genes),
    raw_weighted_rank_score =
      as.numeric(
        raw_p5_rank_score[
          all_donors
        ]
      ),
    matched_background_mean =
      as.numeric(
        background_mean[
          all_donors
        ]
      ),
    matched_background_sd =
      as.numeric(
        background_sd[
          all_donors
        ]
      ),
    calibrated_rank_score =
      as.numeric(
        calibrated_rank_score[
          all_donors
        ]
      ),
    p5_nes =
      as.numeric(
        p5_nes[
          all_donors
        ]
      ),
    stringsAsFactors = FALSE
  )
  
  p5_mapping <- p5_weights %>%
    mutate(
      compartment = compartment_name,
      present_in_tle =
        gene %in%
        rownames(tle_aggregate$counts),
      present_in_reference =
        gene %in%
        rownames(reference_aggregate$counts),
      present_in_common_universe =
        gene %in%
        universe_genes,
      used_for_scoring =
        gene %in%
        p5_genes
    )
  
  eligibility <- bind_rows(
    data.frame(
      compartment = compartment_name,
      donor_id =
        names(tle_aggregate$cell_counts),
      group = "TLE",
      cell_count =
        as.integer(
          tle_aggregate$cell_counts
        ),
      eligible =
        tle_aggregate$cell_counts >=
        minimum_cells,
      stringsAsFactors = FALSE
    ),
    data.frame(
      compartment = compartment_name,
      donor_id =
        names(
          reference_aggregate$cell_counts
        ),
      group = "External_reference",
      cell_count =
        as.integer(
          reference_aggregate$cell_counts
        ),
      eligible =
        reference_aggregate$cell_counts >=
        minimum_cells,
      stringsAsFactors = FALSE
    )
  )
  
  list(
    scores = donor_scores,
    mapping = p5_mapping,
    eligibility = eligibility,
    tle_eligible = tle_eligible,
    reference_eligible =
      reference_eligible,
    universe_gene_count =
      length(universe_genes),
    mapped_p5_gene_count =
      length(p5_genes)
  )
}

mature_result <- score_compartment(
  tle_input = tle_mature,
  reference_input = seaad_mature,
  p5_weights = p5_table,
  compartment_name =
    "Mature_oligodendrocyte"
)

opc_result <- score_compartment(
  tle_input = tle_opc,
  reference_input = seaad_opc,
  p5_weights = p5_table,
  compartment_name = "OPC"
)

all_donor_scores <- bind_rows(
  mature_result$scores,
  opc_result$scores
)

all_p5_mapping <- bind_rows(
  mature_result$mapping,
  opc_result$mapping
)

all_eligibility <- bind_rows(
  mature_result$eligibility,
  opc_result$eligibility
)

write_csv(
  all_donor_scores,
  file.path(
    score_root,
    "P5_DONOR_LEVEL_SCORES.csv"
  )
)

write_csv(
  all_p5_mapping,
  file.path(
    input_audit_root,
    "P5_GENE_MAPPING_BY_COMPARTMENT.csv"
  )
)

write_csv(
  all_eligibility,
  file.path(
    input_audit_root,
    "DONOR_CELL_COUNT_ELIGIBILITY.csv"
  )
)

###############################################################################
# 8. Statistical helpers
###############################################################################

exact_permutation_test <- function(
    tle_values,
    reference_values
) {
  tle_values <- as.numeric(
    tle_values
  )
  
  reference_values <- as.numeric(
    reference_values
  )
  
  all_values <- c(
    tle_values,
    reference_values
  )
  
  n_tle <- length(
    tle_values
  )
  
  observed_difference <- (
    mean(tle_values) -
      mean(reference_values)
  )
  
  combinations <- combn(
    seq_along(all_values),
    n_tle
  )
  
  permutation_differences <- apply(
    combinations,
    MARGIN = 2,
    FUN = function(tle_index) {
      mean(
        all_values[tle_index]
      ) -
        mean(
          all_values[
            -tle_index
          ]
        )
    }
  )
  
  p_value <- mean(
    abs(permutation_differences) >=
      abs(observed_difference) -
      1e-12
  )
  
  list(
    observed_difference =
      observed_difference,
    p_value = p_value,
    permutations =
      length(
        permutation_differences
      )
  )
}

hedges_g <- function(
    tle_values,
    reference_values
) {
  tle_values <- as.numeric(
    tle_values
  )
  
  reference_values <- as.numeric(
    reference_values
  )
  
  n_tle <- length(tle_values)
  n_reference <-
    length(reference_values)
  
  if (
    n_tle < 2L ||
    n_reference < 2L
  ) {
    return(NA_real_)
  }
  
  pooled_variance <- (
    (
      (n_tle - 1) *
        stats::var(tle_values)
    ) +
      (
        (n_reference - 1) *
          stats::var(
            reference_values
          )
      )
  ) / (
    n_tle +
      n_reference -
      2
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
        4 *
          (
            n_tle +
              n_reference
          ) -
          9
      )
  )
  
  correction * cohen_d
}

run_contrast <- function(
    scores,
    comparison_name,
    tle_ids,
    reference_ids
) {
  available_ids <- scores$donor_id
  
  missing_tle <- setdiff(
    tle_ids,
    available_ids
  )
  
  missing_reference <- setdiff(
    reference_ids,
    available_ids
  )
  
  if (
    length(missing_tle) > 0L ||
    length(missing_reference) > 0L
  ) {
    stop(
      comparison_name,
      ": missing requested donor IDs."
    )
  }
  
  tle_values <- scores$p5_nes[
    match(
      tle_ids,
      scores$donor_id
    )
  ]
  
  reference_values <- scores$p5_nes[
    match(
      reference_ids,
      scores$donor_id
    )
  ]
  
  if (
    any(!is.finite(tle_values)) ||
    any(!is.finite(reference_values))
  ) {
    stop(
      comparison_name,
      ": non-finite P5 NES detected."
    )
  }
  
  permutation_result <-
    exact_permutation_test(
      tle_values = tle_values,
      reference_values =
        reference_values
    )
  
  observed_sign <- sign(
    permutation_result$
      observed_difference
  )
  
  all_ids <- c(
    tle_ids,
    reference_ids
  )
  
  lodo_rows <- vector(
    mode = "list",
    length = length(all_ids)
  )
  
  for (index in seq_along(all_ids)) {
    removed_id <- all_ids[index]
    
    reduced_tle <- setdiff(
      tle_ids,
      removed_id
    )
    
    reduced_reference <- setdiff(
      reference_ids,
      removed_id
    )
    
    reduced_tle_values <-
      scores$p5_nes[
        match(
          reduced_tle,
          scores$donor_id
        )
      ]
    
    reduced_reference_values <-
      scores$p5_nes[
        match(
          reduced_reference,
          scores$donor_id
        )
      ]
    
    reduced_difference <- (
      mean(reduced_tle_values) -
        mean(
          reduced_reference_values
        )
    )
    
    lodo_rows[[index]] <- data.frame(
      comparison =
        comparison_name,
      removed_donor =
        removed_id,
      removed_group =
        ifelse(
          removed_id %in%
            tle_ids,
          "TLE",
          "External_reference"
        ),
      remaining_tle_n =
        length(reduced_tle),
      remaining_reference_n =
        length(
          reduced_reference
        ),
      mean_difference =
        reduced_difference,
      direction_preserved =
        sign(
          reduced_difference
        ) == observed_sign,
      stringsAsFactors = FALSE
    )
  }
  
  lodo_table <- bind_rows(
    lodo_rows
  )
  
  statistics <- data.frame(
    comparison =
      comparison_name,
    compartment =
      unique(
        scores$compartment
      ),
    endpoint =
      "P5_matched_background_NES",
    tle_donors =
      paste(
        tle_ids,
        collapse = "|"
      ),
    reference_donors =
      paste(
        reference_ids,
        collapse = "|"
      ),
    tle_n =
      length(tle_values),
    reference_n =
      length(reference_values),
    tle_mean =
      mean(tle_values),
    reference_mean =
      mean(reference_values),
    mean_difference =
      mean(tle_values) -
      mean(reference_values),
    tle_median =
      median(tle_values),
    reference_median =
      median(reference_values),
    median_difference =
      median(tle_values) -
      median(reference_values),
    hedges_g =
      hedges_g(
        tle_values,
        reference_values
      ),
    exact_permutation_p =
      permutation_result$p_value,
    permutation_count =
      permutation_result$
      permutations,
    lodo_direction_stability =
      mean(
        lodo_table$
          direction_preserved
      ),
    direction =
      ifelse(
        permutation_result$
          observed_difference > 0,
        "TLE_HIGHER",
        ifelse(
          permutation_result$
            observed_difference < 0,
          "TLE_LOWER",
          "NO_DIFFERENCE"
        )
      ),
    evidence_type =
      paste0(
        "EXTERNAL_CROSS_STUDY_",
        "CASE_REFERENCE_COMPARISON"
      ),
    stringsAsFactors = FALSE
  )
  
  list(
    statistics = statistics,
    lodo = lodo_table
  )
}

###############################################################################
# 9. Prespecified comparisons
###############################################################################

mature_scores <- all_donor_scores %>%
  filter(
    compartment ==
      "Mature_oligodendrocyte"
  )

opc_scores <- all_donor_scores %>%
  filter(
    compartment == "OPC"
  )

reference_mature_ids <- mature_scores %>%
  filter(
    group ==
      "External_reference"
  ) %>%
  pull(donor_id)

reference_opc_ids <- opc_scores %>%
  filter(
    group ==
      "External_reference"
  ) %>%
  pull(donor_id)

adult_tle_ids <- c(
  "14431",
  "19619",
  "20188"
)

all_tle_ids <- c(
  "12814",
  "13059",
  "14431",
  "19619",
  "20188"
)

paper_included_tle_ids <- c(
  "13059",
  "14431",
  "19619",
  "20188"
)

comparison_results <- list()

comparison_results[["MATURE_ADULT_3_VS_3"]] <- run_contrast(
  scores = mature_scores,
  comparison_name =
    "MATURE_ADULT_3_VS_3",
  tle_ids = adult_tle_ids,
  reference_ids =
    reference_mature_ids
)

comparison_results[["MATURE_FULL_5_VS_3"]] <- run_contrast(
  scores = mature_scores,
  comparison_name =
    "MATURE_FULL_5_VS_3",
  tle_ids = all_tle_ids,
  reference_ids =
    reference_mature_ids
)

comparison_results[["MATURE_PAPER_INCLUDED_4_VS_3"]] <- run_contrast(
  scores = mature_scores,
  comparison_name =
    "MATURE_PAPER_INCLUDED_4_VS_3",
  tle_ids =
    paper_included_tle_ids,
  reference_ids =
    reference_mature_ids
)

eligible_opc_tle_ids <- opc_scores %>%
  filter(
    group == "TLE"
  ) %>%
  pull(donor_id)

if (
  length(eligible_opc_tle_ids) >= 3L &&
  length(reference_opc_ids) >= 3L
) {
  comparison_results[["OPC_ELIGIBLE_TLE_VS_3"]] <- run_contrast(
    scores = opc_scores,
    comparison_name =
      "OPC_ELIGIBLE_TLE_VS_3",
    tle_ids =
      eligible_opc_tle_ids,
    reference_ids =
      reference_opc_ids
  )
}

statistics_table <- bind_rows(
  lapply(
    comparison_results,
    function(result) {
      result$statistics
    }
  )
)

statistics_table$p_adjust_BH <- p.adjust(
  statistics_table$
    exact_permutation_p,
  method = "BH"
)

lodo_table <- bind_rows(
  lapply(
    comparison_results,
    function(result) {
      result$lodo
    }
  )
)

write_csv(
  statistics_table,
  file.path(
    statistics_root,
    "P5_DONOR_LEVEL_COMPARISON_RESULTS.csv"
  )
)

write_csv(
  lodo_table,
  file.path(
    statistics_root,
    "P5_LEAVE_ONE_DONOR_OUT_RESULTS.csv"
  )
)

###############################################################################
# 10. Figures
###############################################################################

plot_donor_scores <- function(
    score_table,
    plot_title,
    output_stub
) {
  plot_data <- score_table %>%
    mutate(
      group = factor(
        group,
        levels = c(
          "External_reference",
          "TLE"
        )
      )
    )
  
  plot <- ggplot(
    plot_data,
    aes(
      x = group,
      y = p5_nes,
      color = group
    )
  ) +
    geom_boxplot(
      width = 0.45,
      outlier.shape = NA,
      alpha = 0.15
    ) +
    geom_jitter(
      width = 0.08,
      height = 0,
      size = 3.2
    ) +
    geom_text(
      aes(
        label = donor_id
      ),
      position = position_jitter(
        width = 0.08,
        height = 0,
        seed = 20260731
      ),
      vjust = -0.9,
      size = 3.2,
      show.legend = FALSE
    ) +
    labs(
      title = plot_title,
      subtitle = paste0(
        "Matched-background frozen P5 rank NES; ",
        "donor is the statistical unit"
      ),
      x = NULL,
      y = "Frozen P5 matched-background NES"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      legend.position = "none"
    )
  
  ggsave(
    filename = file.path(
      figure_root,
      paste0(
        output_stub,
        ".png"
      )
    ),
    plot = plot,
    width = 7,
    height = 6,
    dpi = 320,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      figure_root,
      paste0(
        output_stub,
        ".pdf"
      )
    ),
    plot = plot,
    width = 7,
    height = 6,
    bg = "white"
  )
}

plot_donor_scores(
  score_table = mature_scores,
  plot_title =
    "Mature oligodendrocyte frozen P5 scores",
  output_stub =
    "MATURE_OLIGODENDROCYTE_P5_DONOR_SCORES"
)

plot_donor_scores(
  score_table = opc_scores,
  plot_title =
    "OPC frozen P5 scores",
  output_stub =
    "OPC_P5_DONOR_SCORES"
)

###############################################################################
# 11. Summary
###############################################################################

primary_result <- statistics_table %>%
  filter(
    comparison ==
      "MATURE_ADULT_3_VS_3"
  )

if (nrow(primary_result) != 1L) {
  stop(
    "The primary comparison result was not uniquely identified."
  )
}

summary_lines <- c(
  "# 35C Frozen P5 donor-level scoring",
  "",
  paste0(
    "- Evidence type: ",
    "EXTERNAL_CROSS_STUDY_CASE_REFERENCE_COMPARISON"
  ),
  paste0(
    "- Frozen P5 source: ",
    p5_result$path
  ),
  paste0(
    "- Frozen P5 genes: ",
    nrow(p5_table)
  ),
  paste0(
    "- Mature common P5 genes: ",
    mature_result$mapped_p5_gene_count
  ),
  paste0(
    "- OPC common P5 genes: ",
    opc_result$mapped_p5_gene_count
  ),
  "",
  "## Primary comparison",
  "",
  paste0(
    "- Comparison: ",
    primary_result$comparison
  ),
  paste0(
    "- TLE donors: ",
    primary_result$tle_donors
  ),
  paste0(
    "- External reference donors: ",
    primary_result$reference_donors
  ),
  paste0(
    "- TLE mean P5 NES: ",
    signif(
      primary_result$tle_mean,
      6
    )
  ),
  paste0(
    "- Reference mean P5 NES: ",
    signif(
      primary_result$reference_mean,
      6
    )
  ),
  paste0(
    "- Mean difference: ",
    signif(
      primary_result$mean_difference,
      6
    )
  ),
  paste0(
    "- Hedges g: ",
    signif(
      primary_result$hedges_g,
      6
    )
  ),
  paste0(
    "- Exact permutation P: ",
    signif(
      primary_result$exact_permutation_p,
      6
    )
  ),
  paste0(
    "- LODO direction stability: ",
    signif(
      primary_result$lodo_direction_stability,
      6
    )
  ),
  paste0(
    "- Direction: ",
    primary_result$direction
  ),
  "",
  "## Interpretation boundary",
  "",
  paste0(
    "- GSE140393 is whole-cell scRNA-seq from surgical TLE tissue."
  ),
  paste0(
    "- SEA-AD is postmortem snRNA-seq from an external MTG reference."
  ),
  paste0(
    "- Disease status is fully confounded with dataset and modality."
  ),
  paste0(
    "- Results must not be described as an unconfounded ",
    "independent case-control replication."
  ),
  paste0(
    "- The analysis tests external-reference deviation only."
  ),
  paste0(
    "- Cell-level P5-high fractions were not calculated in this script."
  )
)

writeLines(
  summary_lines,
  file.path(
    output_root,
    "35C_FROZEN_P5_DONOR_SCORING_SUMMARY.md"
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
# 12. Console output
###############################################################################

message("")
message("===================================================")
message("Frozen P5 donor-level scoring completed")
message("===================================================")
message("")
message(
  "Mature mappable P5 genes: ",
  mature_result$mapped_p5_gene_count
)
message(
  "OPC mappable P5 genes: ",
  opc_result$mapped_p5_gene_count
)
message("")
message("Donor scores:")

print(
  as.data.frame(
    all_donor_scores %>%
      select(
        compartment,
        donor_id,
        group,
        cell_count,
        mapped_p5_genes,
        raw_weighted_rank_score,
        calibrated_rank_score,
        p5_nes
      )
  ),
  row.names = FALSE
)

message("")
message("Comparison results:")

print(
  as.data.frame(
    statistics_table %>%
      select(
        comparison,
        tle_n,
        reference_n,
        mean_difference,
        hedges_g,
        exact_permutation_p,
        p_adjust_BH,
        lodo_direction_stability,
        direction
      )
  ),
  row.names = FALSE
)

message("")
message("Output directory:")
message(output_root)

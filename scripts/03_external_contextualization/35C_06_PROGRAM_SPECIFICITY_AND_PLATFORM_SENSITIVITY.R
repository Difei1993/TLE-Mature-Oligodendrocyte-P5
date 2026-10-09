###############################################################################
# 35C-06
# Frozen P5 program-specificity and cross-platform sensitivity analysis
#
# Primary question:
# Is the GSE140393-versus-SEA-AD deviation specific to P5, or do most frozen
# oligodendrocyte programs shift in the same direction?
#
# No clustering, annotation, gene reselection or cell-level testing is done.
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
  "09_PROGRAM_SPECIFICITY_AND_PLATFORM_SENSITIVITY"
)

audit_root <- file.path(output_root, "00_input_audit")
score_root <- file.path(output_root, "01_program_scores")
result_root <- file.path(output_root, "02_results")
figure_root <- file.path(output_root, "03_figures")
log_root <- file.path(output_root, "04_logs")

for (directory in c(
  output_root,
  audit_root,
  score_root,
  result_root,
  figure_root,
  log_root
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
  file.path(seaad_mature_root, "matrix.mtx.gz"),
  file.path(seaad_mature_root, "features.tsv.gz"),
  file.path(seaad_mature_root, "barcodes.tsv.gz"),
  file.path(seaad_mature_root, "cell_metadata.csv.gz"),
  file.path(seaad_opc_root, "matrix.mtx.gz"),
  file.path(seaad_opc_root, "features.tsv.gz"),
  file.path(seaad_opc_root, "barcodes.tsv.gz"),
  file.path(seaad_opc_root, "cell_metadata.csv.gz")
)

missing_paths <- required_paths[
  !file.exists(required_paths)
]

if (length(missing_paths) > 0L) {
  stop(
    "Missing required input files:\n",
    paste(missing_paths, collapse = "\n")
  )
}

###############################################################################
# 3. General helpers
###############################################################################

canonicalize_gene_symbol <- function(symbols) {
  symbols <- toupper(
    trimws(
      as.character(symbols)
    )
  )
  
  symbols[is.na(symbols)] <- ""
  
  symbols <- sub(
    pattern = "\\.[0-9]+$",
    replacement = "",
    x = symbols,
    perl = TRUE
  )
  
  symbols
}

canonicalize_program_name <- function(values) {
  values <- as.character(values)
  
  output <- rep(
    NA_character_,
    length(values)
  )
  
  for (index in seq_along(values)) {
    value <- values[index]
    
    match_result <- regexec(
      pattern = "(?i)oligodendrocytes?.*?p[_ -]*([0-9]+)",
      text = value,
      perl = TRUE
    )
    
    matches <- regmatches(
      value,
      match_result
    )[[1]]
    
    if (length(matches) >= 2L) {
      output[index] <- paste0(
        "Oligodendrocytes__P",
        matches[2]
      )
    }
  }
  
  output
}

first_existing_column <- function(
    column_names,
    candidates
) {
  hits <- intersect(
    candidates,
    column_names
  )
  
  if (length(hits) == 0L) {
    return(NULL)
  }
  
  hits[1L]
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
    check.names = FALSE,
    stringsAsFactors = FALSE
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
  donor_ids <- as.character(donor_ids)
  
  if (length(donor_ids) != ncol(counts)) {
    stop(
      "Donor vector length does not match matrix columns."
    )
  }
  
  if (anyNA(donor_ids) || any(!nzchar(donor_ids))) {
    stop("Missing donor IDs were detected.")
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
    )
  )
  
  pseudobulk <- counts %*% design
  colnames(pseudobulk) <- donors
  
  cell_counts <- table(
    factor(
      donor_ids,
      levels = donors
    )
  )
  
  list(
    counts = pseudobulk,
    cell_counts = setNames(
      as.integer(cell_counts),
      donors
    )
  )
}

log_cpm <- function(counts) {
  dense_counts <- as.matrix(counts)
  
  library_sizes <- colSums(
    dense_counts
  )
  
  if (
    any(!is.finite(library_sizes)) ||
    any(library_sizes <= 0)
  ) {
    stop("Invalid pseudobulk library size.")
  }
  
  cpm <- sweep(
    dense_counts,
    MARGIN = 2,
    STATS = library_sizes,
    FUN = "/"
  ) * 1e6
  
  log1p(cpm)
}

rank_percentile_matrix <- function(expression_matrix) {
  ranked <- apply(
    expression_matrix,
    MARGIN = 2,
    FUN = function(values) {
      if (length(values) <= 1L) {
        return(rep(NA_real_, length(values)))
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
  
  ranked <- as.matrix(ranked)
  
  rownames(ranked) <- rownames(expression_matrix)
  colnames(ranked) <- colnames(expression_matrix)
  
  ranked
}

exact_permutation_test <- function(
    tle_values,
    reference_values
) {
  tle_values <- as.numeric(tle_values)
  reference_values <- as.numeric(reference_values)
  
  all_values <- c(
    tle_values,
    reference_values
  )
  
  n_tle <- length(tle_values)
  
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
    FUN = function(tle_indices) {
      mean(all_values[tle_indices]) -
        mean(all_values[-tle_indices])
    }
  )
  
  p_value <- mean(
    abs(permutation_differences) >=
      abs(observed_difference) -
      1e-12
  )
  
  list(
    difference = observed_difference,
    p_value = p_value,
    permutation_count = length(permutation_differences)
  )
}

hedges_g <- function(
    tle_values,
    reference_values
) {
  tle_values <- as.numeric(tle_values)
  reference_values <- as.numeric(reference_values)
  
  n_tle <- length(tle_values)
  n_reference <- length(reference_values)
  
  if (n_tle < 2L || n_reference < 2L) {
    return(NA_real_)
  }
  
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

lodo_direction_stability <- function(
    scores,
    tle_ids,
    reference_ids
) {
  original_difference <- (
    mean(scores[tle_ids]) -
      mean(scores[reference_ids])
  )
  
  original_direction <- sign(
    original_difference
  )
  
  all_ids <- c(
    tle_ids,
    reference_ids
  )
  
  direction_preserved <- logical(
    length(all_ids)
  )
  
  for (index in seq_along(all_ids)) {
    removed_id <- all_ids[index]
    
    retained_tle <- setdiff(
      tle_ids,
      removed_id
    )
    
    retained_reference <- setdiff(
      reference_ids,
      removed_id
    )
    
    reduced_difference <- (
      mean(scores[retained_tle]) -
        mean(scores[retained_reference])
    )
    
    direction_preserved[index] <- (
      sign(reduced_difference) ==
        original_direction
    )
  }
  
  mean(direction_preserved)
}

###############################################################################
# 4. Read TLE and SEA-AD compartments
###############################################################################

read_tle_compartment <- function(
    object_path,
    compartment
) {
  object <- readRDS(object_path)
  
  counts <- get_seurat_counts(
    object
  )
  
  metadata <- object@meta.data[
    colnames(counts),
    ,
    drop = FALSE
  ]
  
  if (!"donor_id" %in% colnames(metadata)) {
    stop(compartment, " object lacks donor_id.")
  }
  
  counts <- collapse_rows_by_symbol(
    counts,
    rownames(counts)
  )
  
  list(
    counts = counts,
    donor_ids = as.character(metadata$donor_id),
    compartment = compartment
  )
}

read_seaad_compartment <- function(
    directory,
    compartment
) {
  counts <- read_matrix_market_gzip(
    file.path(directory, "matrix.mtx.gz")
  )
  
  features <- read_gzip_table(
    file.path(directory, "features.tsv.gz"),
    header = FALSE
  )
  
  barcodes <- read_gzip_lines(
    file.path(directory, "barcodes.tsv.gz")
  )
  
  metadata <- suppressMessages(
    readr::read_csv(
      file.path(directory, "cell_metadata.csv.gz"),
      show_col_types = FALSE
    )
  )
  
  if (nrow(features) != nrow(counts)) {
    stop(compartment, ": feature count mismatch.")
  }
  
  if (length(barcodes) != ncol(counts)) {
    stop(compartment, ": barcode count mismatch.")
  }
  
  gene_symbols <- if (ncol(features) >= 2L) {
    features[[2L]]
  } else {
    features[[1L]]
  }
  
  rownames(counts) <- canonicalize_gene_symbol(
    gene_symbols
  )
  
  colnames(counts) <- barcodes
  
  if (
    !"cell_barcode" %in% colnames(metadata) ||
    !"donor_id_35c" %in% colnames(metadata)
  ) {
    stop(compartment, ": required metadata columns are missing.")
  }
  
  metadata_index <- match(
    barcodes,
    metadata$cell_barcode
  )
  
  if (anyNA(metadata_index)) {
    stop(compartment, ": unmatched cell barcodes.")
  }
  
  metadata <- metadata[
    metadata_index,
    ,
    drop = FALSE
  ]
  
  counts <- collapse_rows_by_symbol(
    counts,
    rownames(counts)
  )
  
  list(
    counts = counts,
    donor_ids = as.character(metadata$donor_id_35c),
    compartment = compartment
  )
}

tle_mature <- read_tle_compartment(
  tle_mature_path,
  "Mature_oligodendrocyte"
)

tle_opc <- read_tle_compartment(
  tle_opc_path,
  "OPC"
)

seaad_mature <- read_seaad_compartment(
  seaad_mature_root,
  "Mature_oligodendrocyte"
)

seaad_opc <- read_seaad_compartment(
  seaad_opc_root,
  "OPC"
)

###############################################################################
# 5. Locate the official frozen P5
###############################################################################

read_program_csv <- function(path) {
  header <- tryCatch(
    names(
      suppressMessages(
        readr::read_csv(
          path,
          n_max = 0,
          show_col_types = FALSE
        )
      )
    ),
    error = function(error) {
      character(0)
    }
  )
  
  if (length(header) == 0L) {
    return(NULL)
  }
  
  gene_column <- first_existing_column(
    header,
    c(
      "gene",
      "gene_symbol",
      "symbol",
      "human_gene",
      "feature"
    )
  )
  
  weight_column <- first_existing_column(
    header,
    c(
      "meta_program_score",
      "weight",
      "program_weight",
      "gene_weight",
      "loading",
      "score"
    )
  )
  
  program_column <- first_existing_column(
    header,
    c(
      "program",
      "program_id",
      "program_name",
      "meta_program",
      "meta_program_name"
    )
  )
  
  if (
    is.null(gene_column) ||
    is.null(weight_column)
  ) {
    return(NULL)
  }
  
  data <- tryCatch(
    suppressMessages(
      readr::read_csv(
        path,
        show_col_types = FALSE
      )
    ),
    error = function(error) {
      NULL
    }
  )
  
  if (is.null(data)) {
    return(NULL)
  }
  
  inferred_program <- canonicalize_program_name(
    basename(path)
  )[1L]
  
  if (!is.null(program_column)) {
    programs <- canonicalize_program_name(
      data[[program_column]]
    )
  } else {
    programs <- rep(
      inferred_program,
      nrow(data)
    )
  }
  
  result <- data.frame(
    program = programs,
    gene = canonicalize_gene_symbol(
      data[[gene_column]]
    ),
    weight = suppressWarnings(
      as.numeric(data[[weight_column]])
    ),
    source_file = normalizePath(
      path,
      winslash = "/",
      mustWork = FALSE
    ),
    stringsAsFactors = FALSE
  ) %>%
    filter(
      !is.na(program),
      nzchar(gene),
      is.finite(weight),
      weight > 0
    )
  
  if (nrow(result) == 0L) {
    return(NULL)
  }
  
  result
}

official_p5_files <- list.files(
  frozen_root,
  pattern = "34A5C_FROZEN_HUMAN_TOP100\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(official_p5_files) == 0L) {
  stop(
    "Official 34A5C frozen P5 Top100 CSV was not found."
  )
}

official_p5_files <- official_p5_files[
  order(
    file.info(official_p5_files)$mtime,
    decreasing = TRUE
  )
]

official_p5_raw <- read_program_csv(
  official_p5_files[1L]
)

if (is.null(official_p5_raw)) {
  stop("The official P5 file could not be parsed.")
}

official_p5 <- official_p5_raw %>%
  transmute(
    program = "Oligodendrocytes__P5",
    gene,
    weight,
    source_file
  ) %>%
  group_by(program, gene, source_file) %>%
  summarise(
    weight = max(weight),
    .groups = "drop"
  ) %>%
  arrange(desc(weight)) %>%
  slice_head(n = 100L)

if (nrow(official_p5) != 100L) {
  stop(
    "Official frozen P5 does not contain exactly 100 unique genes."
  )
}

###############################################################################
# 6. Discover other frozen oligodendrocyte programs
###############################################################################

candidate_csv_files <- list.files(
  frozen_root,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

candidate_csv_files <- candidate_csv_files[
  grepl(
    pattern = paste0(
      "program|meta.program|cnmf|nmf|top100|",
      "gene.*weight|loading"
    ),
    x = basename(candidate_csv_files),
    ignore.case = TRUE
  )
]

candidate_csv_files <- candidate_csv_files[
  file.info(candidate_csv_files)$size <
    100 * 1024^2
]

candidate_file_audit <- data.frame(
  file = candidate_csv_files,
  size_bytes = file.info(candidate_csv_files)$size,
  parsed = FALSE,
  oligodendrocyte_program_rows = 0L,
  stringsAsFactors = FALSE
)

program_rows <- list()

for (file_index in seq_along(candidate_csv_files)) {
  path <- candidate_csv_files[file_index]
  
  parsed <- read_program_csv(path)
  
  if (is.null(parsed)) {
    next
  }
  
  parsed <- parsed %>%
    filter(
      grepl(
        "^Oligodendrocytes__P[0-9]+$",
        program
      )
    )
  
  candidate_file_audit$parsed[file_index] <- TRUE
  candidate_file_audit$
    oligodendrocyte_program_rows[file_index] <- nrow(parsed)
  
  if (nrow(parsed) > 0L) {
    program_rows[[length(program_rows) + 1L]] <- parsed
  }
}

write_csv(
  candidate_file_audit,
  file.path(
    audit_root,
    "OLIGODENDROCYTE_PROGRAM_CANDIDATE_FILES.csv"
  )
)

if (length(program_rows) == 0L) {
  stop(
    paste0(
      "No additional frozen oligodendrocyte program table was found.\n",
      "Inspect: ",
      file.path(
        audit_root,
        "OLIGODENDROCYTE_PROGRAM_CANDIDATE_FILES.csv"
      )
    )
  )
}

discovered_program_rows <- bind_rows(
  program_rows
)

source_summary <- discovered_program_rows %>%
  group_by(program, source_file) %>%
  summarise(
    gene_count = n_distinct(gene),
    .groups = "drop"
  )

selected_sources <- source_summary %>%
  group_by(program) %>%
  arrange(desc(gene_count), source_file) %>%
  slice_head(n = 1L) %>%
  ungroup()

all_programs <- discovered_program_rows %>%
  inner_join(
    selected_sources %>%
      select(program, source_file),
    by = c("program", "source_file")
  ) %>%
  group_by(program, gene) %>%
  summarise(
    weight = max(weight),
    source_file = first(source_file),
    .groups = "drop"
  ) %>%
  group_by(program) %>%
  arrange(desc(weight), .by_group = TRUE) %>%
  slice_head(n = 100L) %>%
  ungroup()

# The official frozen P5 always overrides any other P5 copy.
all_programs <- all_programs %>%
  filter(program != "Oligodendrocytes__P5") %>%
  bind_rows(official_p5) %>%
  arrange(
    as.integer(
      sub(
        "^Oligodendrocytes__P",
        "",
        program
      )
    ),
    desc(weight)
  )

program_registry <- all_programs %>%
  group_by(program, source_file) %>%
  summarise(
    gene_count = n_distinct(gene),
    weight_sum = sum(weight),
    .groups = "drop"
  )

write_csv(
  program_registry,
  file.path(
    audit_root,
    "FROZEN_OLIGODENDROCYTE_PROGRAM_REGISTRY.csv"
  )
)

write_csv(
  all_programs,
  file.path(
    audit_root,
    "FROZEN_OLIGODENDROCYTE_PROGRAM_TOP_GENES.csv"
  )
)

program_names <- unique(
  all_programs$program
)

if (
  !"Oligodendrocytes__P5" %in%
  program_names
) {
  stop("P5 was lost during program discovery.")
}

if (length(program_names) < 2L) {
  stop(
    "Only P5 was found; program-specificity analysis cannot proceed."
  )
}

message("")
message(
  "Frozen oligodendrocyte programs found: ",
  paste(program_names, collapse = ", ")
)

###############################################################################
# 7. Prepare one compartment
###############################################################################

prepare_compartment <- function(
    tle_input,
    reference_input,
    minimum_cells = 30L
) {
  tle_aggregate <- aggregate_counts_by_donor(
    tle_input$counts,
    tle_input$donor_ids
  )
  
  reference_aggregate <- aggregate_counts_by_donor(
    reference_input$counts,
    reference_input$donor_ids
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
  
  common_genes <- intersect(
    rownames(tle_aggregate$counts),
    rownames(reference_aggregate$counts)
  )
  
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
  
  jointly_expressed <- (
    Matrix::rowSums(tle_pseudobulk) > 0 &
      Matrix::rowSums(reference_pseudobulk) > 0
  )
  
  universe_genes <- common_genes[
    jointly_expressed
  ]
  
  if (length(universe_genes) < 3000L) {
    stop("Too few jointly expressed common genes.")
  }
  
  combined_expression <- cbind(
    log_cpm(
      tle_pseudobulk[
        universe_genes,
        ,
        drop = FALSE
      ]
    ),
    log_cpm(
      reference_pseudobulk[
        universe_genes,
        ,
        drop = FALSE
      ]
    )
  )
  
  rank_matrix <- rank_percentile_matrix(
    combined_expression
  )
  
  mean_expression <- rowMeans(
    combined_expression
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
    pmax(1L, expression_bin)
  )
  
  names(expression_bin) <- universe_genes
  
  donor_groups <- ifelse(
    colnames(combined_expression) %in%
      tle_eligible,
    "TLE",
    "External_reference"
  )
  
  names(donor_groups) <- colnames(combined_expression)
  
  list(
    rank_matrix = rank_matrix,
    expression_bin = expression_bin,
    universe_genes = universe_genes,
    donor_groups = donor_groups,
    tle_eligible = tle_eligible,
    reference_eligible = reference_eligible,
    tle_cell_counts = tle_aggregate$cell_counts,
    reference_cell_counts = reference_aggregate$cell_counts
  )
}

mature_prepared <- prepare_compartment(
  tle_mature,
  seaad_mature
)

opc_prepared <- prepare_compartment(
  tle_opc,
  seaad_opc
)

###############################################################################
# 8. Score one frozen program
###############################################################################

score_one_program <- function(
    prepared,
    program_table,
    background_replicates = 200L
) {
  mapped_genes <- intersect(
    program_table$gene,
    prepared$universe_genes
  )
  
  if (length(mapped_genes) < 20L) {
    return(NULL)
  }
  
  weights <- program_table$weight[
    match(
      mapped_genes,
      program_table$gene
    )
  ]
  
  weights <- weights / sum(weights)
  
  raw_weighted_score <- colSums(
    sweep(
      prepared$rank_matrix[
        mapped_genes,
        ,
        drop = FALSE
      ],
      MARGIN = 1,
      STATS = weights,
      FUN = "*"
    )
  )
  
  unweighted_score <- colMeans(
    prepared$rank_matrix[
      mapped_genes,
      ,
      drop = FALSE
    ]
  ) - 0.5
  
  program_bins <- prepared$expression_bin[
    mapped_genes
  ]
  
  non_program_genes <- setdiff(
    prepared$universe_genes,
    mapped_genes
  )
  
  background_scores <- matrix(
    NA_real_,
    nrow = ncol(prepared$rank_matrix),
    ncol = background_replicates,
    dimnames = list(
      colnames(prepared$rank_matrix),
      paste0("background_", seq_len(background_replicates))
    )
  )
  
  for (replicate_index in seq_len(background_replicates)) {
    sampled_genes <- character(
      length(mapped_genes)
    )
    
    for (bin_value in sort(unique(program_bins))) {
      program_indices <- which(
        program_bins == bin_value
      )
      
      candidate_pool <- non_program_genes[
        prepared$expression_bin[
          non_program_genes
        ] == bin_value
      ]
      
      if (length(candidate_pool) == 0L) {
        candidate_pool <- non_program_genes
      }
      
      sampled_genes[program_indices] <- sample(
        candidate_pool,
        size = length(program_indices),
        replace = (
          length(candidate_pool) <
            length(program_indices)
        )
      )
    }
    
    background_scores[, replicate_index] <- colSums(
      sweep(
        prepared$rank_matrix[
          sampled_genes,
          ,
          drop = FALSE
        ],
        MARGIN = 1,
        STATS = weights,
        FUN = "*"
      )
    )
  }
  
  background_mean <- rowMeans(
    background_scores
  )
  
  background_sd <- apply(
    background_scores,
    MARGIN = 1,
    FUN = stats::sd
  )
  
  background_sd[
    !is.finite(background_sd) |
      background_sd <= 0
  ] <- NA_real_
  
  nes <- (
    raw_weighted_score -
      background_mean
  ) / background_sd
  
  data.frame(
    donor_id = names(raw_weighted_score),
    group = unname(
      prepared$donor_groups[
        names(raw_weighted_score)
      ]
    ),
    mapped_genes = length(mapped_genes),
    raw_weighted_rank = as.numeric(raw_weighted_score),
    unweighted_rank_centered = as.numeric(unweighted_score),
    matched_background_nes = as.numeric(nes),
    stringsAsFactors = FALSE
  )
}

###############################################################################
# 9. Score all programs in mature oligodendrocytes and OPC
###############################################################################

score_all_programs <- function(
    prepared,
    program_definitions,
    compartment
) {
  output <- list()
  
  programs <- unique(
    program_definitions$program
  )
  
  for (program_name in programs) {
    message(
      "Scoring ",
      compartment,
      " / ",
      program_name
    )
    
    program_table <- program_definitions %>%
      filter(program == program_name)
    
    scored <- score_one_program(
      prepared = prepared,
      program_table = program_table
    )
    
    if (is.null(scored)) {
      next
    }
    
    scored$program <- program_name
    scored$compartment <- compartment
    
    output[[length(output) + 1L]] <- scored
  }
  
  bind_rows(output)
}

mature_program_scores <- score_all_programs(
  prepared = mature_prepared,
  program_definitions = all_programs,
  compartment = "Mature_oligodendrocyte"
)

opc_program_scores <- score_all_programs(
  prepared = opc_prepared,
  program_definitions = all_programs,
  compartment = "OPC"
)

all_program_scores <- bind_rows(
  mature_program_scores,
  opc_program_scores
)

write_csv(
  all_program_scores,
  file.path(
    score_root,
    "ALL_OLIGODENDROCYTE_PROGRAM_DONOR_SCORES.csv"
  )
)

###############################################################################
# 10. Program-level comparisons
###############################################################################

run_program_contrast <- function(
    score_table,
    tle_ids,
    reference_ids,
    endpoint
) {
  score_vector <- score_table[[endpoint]]
  names(score_vector) <- score_table$donor_id
  
  missing_ids <- setdiff(
    c(tle_ids, reference_ids),
    names(score_vector)
  )
  
  if (length(missing_ids) > 0L) {
    stop(
      "Missing donors: ",
      paste(missing_ids, collapse = ", ")
    )
  }
  
  tle_values <- score_vector[tle_ids]
  reference_values <- score_vector[reference_ids]
  
  if (
    any(!is.finite(tle_values)) ||
    any(!is.finite(reference_values))
  ) {
    return(NULL)
  }
  
  permutation <- exact_permutation_test(
    tle_values,
    reference_values
  )
  
  data.frame(
    endpoint = endpoint,
    tle_n = length(tle_values),
    reference_n = length(reference_values),
    tle_mean = mean(tle_values),
    reference_mean = mean(reference_values),
    mean_difference = (
      mean(tle_values) -
        mean(reference_values)
    ),
    hedges_g = hedges_g(
      tle_values,
      reference_values
    ),
    exact_permutation_p = permutation$p_value,
    permutation_count = permutation$permutation_count,
    lodo_direction_stability =
      lodo_direction_stability(
        scores = score_vector,
        tle_ids = tle_ids,
        reference_ids = reference_ids
      ),
    direction = ifelse(
      permutation$difference > 0,
      "TLE_HIGHER",
      ifelse(
        permutation$difference < 0,
        "TLE_LOWER",
        "NO_DIFFERENCE"
      )
    ),
    stringsAsFactors = FALSE
  )
}

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

mature_reference_ids <- mature_prepared$reference_eligible
opc_reference_ids <- opc_prepared$reference_eligible
opc_tle_ids <- opc_prepared$tle_eligible

comparison_definitions <- list(
  list(
    name = "MATURE_ADULT_3_VS_3",
    compartment = "Mature_oligodendrocyte",
    tle_ids = adult_tle_ids,
    reference_ids = mature_reference_ids
  ),
  list(
    name = "MATURE_FULL_5_VS_3",
    compartment = "Mature_oligodendrocyte",
    tle_ids = all_tle_ids,
    reference_ids = mature_reference_ids
  ),
  list(
    name = "MATURE_PAPER_INCLUDED_4_VS_3",
    compartment = "Mature_oligodendrocyte",
    tle_ids = paper_included_tle_ids,
    reference_ids = mature_reference_ids
  ),
  list(
    name = "OPC_ELIGIBLE_TLE_VS_3",
    compartment = "OPC",
    tle_ids = opc_tle_ids,
    reference_ids = opc_reference_ids
  )
)

endpoints <- c(
  "raw_weighted_rank",
  "unweighted_rank_centered",
  "matched_background_nes"
)

comparison_rows <- list()

for (comparison in comparison_definitions) {
  compartment_scores <- all_program_scores %>%
    filter(
      compartment == comparison$compartment
    )
  
  for (program_name in unique(compartment_scores$program)) {
    program_scores <- compartment_scores %>%
      filter(program == program_name)
    
    for (endpoint in endpoints) {
      result <- run_program_contrast(
        score_table = program_scores,
        tle_ids = comparison$tle_ids,
        reference_ids = comparison$reference_ids,
        endpoint = endpoint
      )
      
      if (is.null(result)) {
        next
      }
      
      result$comparison <- comparison$name
      result$compartment <- comparison$compartment
      result$program <- program_name
      
      comparison_rows[[length(comparison_rows) + 1L]] <- result
    }
  }
}

program_comparison_results <- bind_rows(
  comparison_rows
) %>%
  group_by(
    comparison,
    compartment,
    endpoint
  ) %>%
  mutate(
    p_adjust_BH = p.adjust(
      exact_permutation_p,
      method = "BH"
    )
  ) %>%
  ungroup()

write_csv(
  program_comparison_results,
  file.path(
    result_root,
    "ALL_PROGRAM_COMPARISON_RESULTS.csv"
  )
)

###############################################################################
# 11. Mature adult 3-vs-3 program specificity summary
###############################################################################

primary_program_results <- program_comparison_results %>%
  filter(
    comparison == "MATURE_ADULT_3_VS_3",
    endpoint == "raw_weighted_rank"
  ) %>%
  arrange(desc(mean_difference)) %>%
  mutate(
    effect_rank = row_number(),
    effect_percentile = ifelse(
      n() > 1L,
      1 - (effect_rank - 1) / (n() - 1),
      1
    )
  )

p5_primary_row <- primary_program_results %>%
  filter(program == "Oligodendrocytes__P5")

if (nrow(p5_primary_row) != 1L) {
  stop(
    "P5 was not uniquely identified in the primary program comparison."
  )
}

p5_direction <- p5_primary_row$direction

same_direction_fraction <- mean(
  primary_program_results$direction ==
    p5_direction
)

specificity_summary <- data.frame(
  program_count = nrow(primary_program_results),
  p5_effect_rank = p5_primary_row$effect_rank,
  p5_effect_percentile = p5_primary_row$effect_percentile,
  p5_mean_difference = p5_primary_row$mean_difference,
  p5_hedges_g = p5_primary_row$hedges_g,
  p5_exact_permutation_p =
    p5_primary_row$exact_permutation_p,
  p5_lodo_direction_stability =
    p5_primary_row$lodo_direction_stability,
  p5_direction = p5_direction,
  fraction_of_programs_same_direction =
    same_direction_fraction,
  global_same_direction_warning =
    same_direction_fraction >= 0.75,
  p5_in_top_quartile =
    p5_primary_row$effect_rank <=
    ceiling(nrow(primary_program_results) / 4),
  stringsAsFactors = FALSE
)

write_csv(
  primary_program_results,
  file.path(
    result_root,
    "MATURE_ADULT_3_VS_3_PROGRAM_SPECIFICITY_RANKING.csv"
  )
)

write_csv(
  specificity_summary,
  file.path(
    result_root,
    "P5_PROGRAM_SPECIFICITY_SUMMARY.csv"
  )
)

###############################################################################
# 12. P5 scoring-method and gene-set sensitivity
###############################################################################

major_myelin_genes <- c(
  "PLP1",
  "CNP",
  "ERMN",
  "MBP",
  "MOG",
  "MOBP",
  "MAG",
  "CLDN11"
)

p5_ranked <- official_p5 %>%
  arrange(desc(weight))

p5_top10_genes <- head(
  p5_ranked$gene,
  10L
)

p5_sensitivity_programs <- bind_rows(
  official_p5 %>%
    transmute(
      program = "P5_FROZEN_WEIGHTED",
      gene,
      weight
    ),
  
  official_p5 %>%
    transmute(
      program = "P5_EQUAL_WEIGHT",
      gene,
      weight = 1
    ),
  
  official_p5 %>%
    filter(
      !gene %in% p5_top10_genes
    ) %>%
    transmute(
      program = "P5_DROP_TOP10_WEIGHTED",
      gene,
      weight
    ),
  
  official_p5 %>%
    filter(
      !gene %in% major_myelin_genes
    ) %>%
    transmute(
      program = "P5_DROP_MAJOR_MYELIN",
      gene,
      weight
    )
)

p5_sensitivity_scores <- score_all_programs(
  prepared = mature_prepared,
  program_definitions = p5_sensitivity_programs,
  compartment = "Mature_oligodendrocyte"
)

p5_sensitivity_rows <- list()

for (
  sensitivity_name in
  unique(p5_sensitivity_scores$program)
) {
  sensitivity_scores <- p5_sensitivity_scores %>%
    filter(program == sensitivity_name)
  
  for (endpoint in endpoints) {
    result <- run_program_contrast(
      score_table = sensitivity_scores,
      tle_ids = adult_tle_ids,
      reference_ids = mature_reference_ids,
      endpoint = endpoint
    )
    
    if (is.null(result)) {
      next
    }
    
    result$sensitivity_definition <- sensitivity_name
    
    p5_sensitivity_rows[[length(p5_sensitivity_rows) + 1L]] <- result
  }
}

p5_sensitivity_results <- bind_rows(
  p5_sensitivity_rows
)

write_csv(
  p5_sensitivity_scores,
  file.path(
    score_root,
    "P5_METHOD_AND_GENESET_SENSITIVITY_DONOR_SCORES.csv"
  )
)

write_csv(
  p5_sensitivity_results,
  file.path(
    result_root,
    "P5_METHOD_AND_GENESET_SENSITIVITY_RESULTS.csv"
  )
)

###############################################################################
# 13. Per-gene contribution to the mature adult 3-vs-3 difference
###############################################################################

p5_mapped_genes <- intersect(
  official_p5$gene,
  mature_prepared$universe_genes
)

p5_weights <- official_p5$weight[
  match(
    p5_mapped_genes,
    official_p5$gene
  )
]

p5_weights <- p5_weights / sum(p5_weights)

adult_tle_mean_ranks <- rowMeans(
  mature_prepared$rank_matrix[
    p5_mapped_genes,
    adult_tle_ids,
    drop = FALSE
  ]
)

reference_mean_ranks <- rowMeans(
  mature_prepared$rank_matrix[
    p5_mapped_genes,
    mature_reference_ids,
    drop = FALSE
  ]
)

gene_contribution_table <- data.frame(
  gene = p5_mapped_genes,
  normalized_weight = p5_weights,
  tle_mean_rank_percentile = adult_tle_mean_ranks,
  reference_mean_rank_percentile =
    reference_mean_ranks,
  rank_difference = (
    adult_tle_mean_ranks -
      reference_mean_ranks
  ),
  weighted_contribution = (
    adult_tle_mean_ranks -
      reference_mean_ranks
  ) * p5_weights,
  major_myelin_gene =
    p5_mapped_genes %in%
    major_myelin_genes,
  top10_weight_gene =
    p5_mapped_genes %in%
    p5_top10_genes,
  stringsAsFactors = FALSE
) %>%
  arrange(
    desc(abs(weighted_contribution))
  )

write_csv(
  gene_contribution_table,
  file.path(
    result_root,
    "P5_PER_GENE_WEIGHTED_CONTRIBUTIONS.csv"
  )
)

myelin_contribution_fraction <- (
  sum(
    gene_contribution_table$weighted_contribution[
      gene_contribution_table$major_myelin_gene
    ]
  ) /
    sum(
      gene_contribution_table$weighted_contribution
    )
)

top10_contribution_fraction <- (
  sum(
    gene_contribution_table$weighted_contribution[
      gene_contribution_table$top10_weight_gene
    ]
  ) /
    sum(
      gene_contribution_table$weighted_contribution
    )
)

contribution_summary <- data.frame(
  mapped_p5_genes = nrow(gene_contribution_table),
  net_weighted_difference =
    sum(gene_contribution_table$weighted_contribution),
  major_myelin_contribution_fraction =
    myelin_contribution_fraction,
  top10_weight_contribution_fraction =
    top10_contribution_fraction,
  stringsAsFactors = FALSE
)

write_csv(
  contribution_summary,
  file.path(
    result_root,
    "P5_GENE_CONTRIBUTION_SUMMARY.csv"
  )
)

###############################################################################
# 14. Figures
###############################################################################

program_effect_plot <- ggplot(
  primary_program_results,
  aes(
    x = reorder(
      program,
      mean_difference
    ),
    y = mean_difference
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  coord_flip() +
  labs(
    title = paste0(
      "Frozen oligodendrocyte program deviations: ",
      "adult TLE versus SEA-AD reference"
    ),
    subtitle = paste0(
      "Raw weighted within-donor rank score; ",
      "donor is the statistical unit"
    ),
    x = NULL,
    y = "Mean difference: TLE minus reference"
  ) +
  theme_classic(
    base_size = 11
  )

ggsave(
  filename = file.path(
    figure_root,
    "ALL_PROGRAM_EFFECT_RANKING.png"
  ),
  plot = program_effect_plot,
  width = 9,
  height = 6,
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = file.path(
    figure_root,
    "ALL_PROGRAM_EFFECT_RANKING.pdf"
  ),
  plot = program_effect_plot,
  width = 9,
  height = 6,
  bg = "white"
)

sensitivity_plot_data <- p5_sensitivity_results %>%
  filter(endpoint == "raw_weighted_rank")

sensitivity_plot <- ggplot(
  sensitivity_plot_data,
  aes(
    x = reorder(
      sensitivity_definition,
      mean_difference
    ),
    y = mean_difference
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  coord_flip() +
  labs(
    title = "P5 scoring and gene-set sensitivity",
    subtitle = "Adult TLE 3 versus SEA-AD reference 3",
    x = NULL,
    y = "Mean difference: TLE minus reference"
  ) +
  theme_classic(
    base_size = 11
  )

ggsave(
  filename = file.path(
    figure_root,
    "P5_METHOD_GENESET_SENSITIVITY.png"
  ),
  plot = sensitivity_plot,
  width = 9,
  height = 5.5,
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = file.path(
    figure_root,
    "P5_METHOD_GENESET_SENSITIVITY.pdf"
  ),
  plot = sensitivity_plot,
  width = 9,
  height = 5.5,
  bg = "white"
)

###############################################################################
# 15. Summary
###############################################################################

summary_lines <- c(
  "# 35C-06 Program specificity and platform sensitivity",
  "",
  paste0(
    "- Frozen oligodendrocyte programs analysed: ",
    nrow(primary_program_results)
  ),
  paste0(
    "- P5 effect rank: ",
    specificity_summary$p5_effect_rank,
    " of ",
    specificity_summary$program_count
  ),
  paste0(
    "- P5 effect percentile: ",
    signif(
      specificity_summary$p5_effect_percentile,
      5
    )
  ),
  paste0(
    "- P5 mean rank-score difference: ",
    signif(
      specificity_summary$p5_mean_difference,
      6
    )
  ),
  paste0(
    "- Fraction of programs with the P5 direction: ",
    signif(
      specificity_summary$
        fraction_of_programs_same_direction,
      5
    )
  ),
  paste0(
    "- Global same-direction warning: ",
    specificity_summary$
      global_same_direction_warning
  ),
  paste0(
    "- P5 in top program-effect quartile: ",
    specificity_summary$
      p5_in_top_quartile
  ),
  "",
  "## Gene-contribution sensitivity",
  "",
  paste0(
    "- Major-myelin-gene contribution fraction: ",
    signif(
      myelin_contribution_fraction,
      6
    )
  ),
  paste0(
    "- Top-10-weight-gene contribution fraction: ",
    signif(
      top10_contribution_fraction,
      6
    )
  ),
  "",
  "## Interpretation",
  "",
  paste0(
    "- This analysis evaluates program specificity and ",
    "cross-platform sensitivity."
  ),
  paste0(
    "- It does not remove the complete confounding between ",
    "dataset, modality and disease status."
  ),
  paste0(
    "- Cell-level P5-high fractions remain deferred."
  )
)

writeLines(
  summary_lines,
  file.path(
    output_root,
    "35C_06_PROGRAM_SPECIFICITY_SUMMARY.md"
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
# 16. Console output
###############################################################################

message("")
message("===================================================")
message("Program specificity and sensitivity analysis done")
message("===================================================")
message("")
message(
  "Programs analysed: ",
  nrow(primary_program_results)
)
message(
  "P5 effect rank: ",
  specificity_summary$p5_effect_rank,
  " / ",
  specificity_summary$program_count
)
message(
  "P5 effect percentile: ",
  signif(
    specificity_summary$p5_effect_percentile,
    5
  )
)
message(
  "Programs with same direction as P5: ",
  signif(
    specificity_summary$
      fraction_of_programs_same_direction,
    5
  )
)
message(
  "Global same-direction warning: ",
  specificity_summary$
    global_same_direction_warning
)
message("")
message("P5 sensitivity results:")

print(
  as.data.frame(
    p5_sensitivity_results %>%
      filter(
        endpoint == "raw_weighted_rank"
      ) %>%
      select(
        sensitivity_definition,
        mean_difference,
        hedges_g,
        exact_permutation_p,
        lodo_direction_stability,
        direction
      )
  ),
  row.names = FALSE
)

message("")
message("Top P5 gene contributions:")

print(
  head(
    as.data.frame(
      gene_contribution_table %>%
        select(
          gene,
          normalized_weight,
          rank_difference,
          weighted_contribution,
          major_myelin_gene,
          top10_weight_gene
        )
    ),
    20L
  ),
  row.names = FALSE
)

message("")
message("Output directory:")
message(output_root)

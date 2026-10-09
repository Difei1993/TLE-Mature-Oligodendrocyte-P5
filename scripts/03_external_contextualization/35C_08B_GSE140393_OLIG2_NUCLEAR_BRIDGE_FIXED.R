###############################################################################
# 35C-08B
# GSE140393 OLIG2-positive sorted-nuclei P5 bridge analysis
# Fixed version using exact internal matrix column names
###############################################################################

options(
  stringsAsFactors = FALSE,
  warn = 1
)

set.seed(20260801)

###############################################################################
# 1. Packages
###############################################################################

required_packages <- c(
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
  stop(
    "缺少R包：",
    paste(missing_packages, collapse = ", "),
    "\n请先安装后重新运行。",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
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

matrix_path <- file.path(
  work_root,
  "data",
  "raw",
  "GSE140393",
  "GSE140393_rldnormalized_nuclei.txt.gz"
)

p5_path <- file.path(
  work_root,
  "08_FROZEN_P5_DONOR_SCORING",
  "00_input_audit",
  "FROZEN_P5_TOP100_AND_WEIGHTS.csv"
)

output_root <- file.path(
  work_root,
  "11_GSE140393_OLIG2_NUCLEAR_BRIDGE_FIXED"
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

if (!file.exists(matrix_path)) {
  stop(
    "找不到核转录组矩阵：",
    matrix_path,
    call. = FALSE
  )
}

if (!file.exists(p5_path)) {
  stop(
    "找不到冻结P5文件：",
    p5_path,
    call. = FALSE
  )
}

###############################################################################
# 3. Fixed OLIG2 sample-to-column mapping
###############################################################################

sample_metadata <- data.frame(
  gsm = c(
    "GSM4160557",
    "GSM4160560",
    "GSM4160563",
    "GSM4160566",
    "GSM4160569",
    "GSM4160572"
  ),
  
  donor_id = c(
    "10662",
    "10355",
    "10997",
    "12726",
    "12321",
    "10308"
  ),
  
  group = c(
    "Control",
    "Control",
    "Control",
    "TLE",
    "TLE",
    "TLE"
  ),
  
  matrix_column = c(
    "Sample_A1-E1-10308-TCCGGAGA-ATAGAGGC",
    "Sample_A2-E2-14431-TCCGGAGA-ATAGAGGC",
    "Sample_A3-E4-TS-TCCGGAGA-ATAGAGGC",
    
    "Sample_A1-E1-10308-ATTCAGAA-ATAGAGGC",
    "Sample_A2-E2-14431-ATTCAGAA-ATAGAGGC",
    "Sample_A3-E4-TS-ATTCAGAA-ATAGAGGC"
  ),
  
  lineage = "OLIG2_positive",
  
  stringsAsFactors = FALSE
)

if (
  sum(sample_metadata$group == "Control") != 3L ||
  sum(sample_metadata$group == "TLE") != 3L
) {
  stop(
    "样本映射表不是3例Control和3例TLE。",
    call. = FALSE
  )
}

if (
  anyDuplicated(sample_metadata$matrix_column) ||
  anyDuplicated(sample_metadata$donor_id)
) {
  stop(
    "样本映射表包含重复列或重复供者。",
    call. = FALSE
  )
}

write_csv(
  sample_metadata,
  file.path(
    table_root,
    "GSE140393_FIXED_OLIG2_SAMPLE_MAPPING.csv"
  )
)

###############################################################################
# 4. Helper functions
###############################################################################

canonicalize_gene <- function(x) {
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  x[is.na(x)] <- ""
  
  sub(
    pattern = "\\.[0-9]+$",
    replacement = "",
    x = x,
    perl = TRUE
  )
}

collapse_duplicate_genes <- function(expression_matrix) {
  gene_symbols <- rownames(expression_matrix)
  
  if (!anyDuplicated(gene_symbols)) {
    return(expression_matrix)
  }
  
  duplicate_counts <- table(gene_symbols)
  
  summed_matrix <- rowsum(
    expression_matrix,
    group = gene_symbols,
    reorder = FALSE
  )
  
  sweep(
    summed_matrix,
    MARGIN = 1,
    STATS = as.numeric(
      duplicate_counts[
        rownames(summed_matrix)
      ]
    ),
    FUN = "/"
  )
}

exact_permutation_test <- function(
    tle_values,
    control_values
) {
  tle_values <- as.numeric(tle_values)
  control_values <- as.numeric(control_values)
  
  all_values <- c(
    tle_values,
    control_values
  )
  
  n_tle <- length(tle_values)
  
  observed_difference <- (
    mean(tle_values) -
      mean(control_values)
  )
  
  label_combinations <- combn(
    seq_along(all_values),
    n_tle
  )
  
  permutation_differences <- apply(
    label_combinations,
    2,
    function(tle_indices) {
      mean(all_values[tle_indices]) -
        mean(all_values[-tle_indices])
    }
  )
  
  list(
    difference = observed_difference,
    
    p_value = mean(
      abs(permutation_differences) >=
        abs(observed_difference) -
        1e-12
    ),
    
    permutation_count =
      length(permutation_differences)
  )
}

hedges_g <- function(
    tle_values,
    control_values
) {
  tle_values <- as.numeric(tle_values)
  control_values <- as.numeric(control_values)
  
  n_tle <- length(tle_values)
  n_control <- length(control_values)
  
  pooled_variance <- (
    (n_tle - 1L) * stats::var(tle_values) +
      (n_control - 1L) * stats::var(control_values)
  ) / (
    n_tle + n_control - 2L
  )
  
  if (
    !is.finite(pooled_variance) ||
    pooled_variance <= 0
  ) {
    return(NA_real_)
  }
  
  cohen_d <- (
    mean(tle_values) -
      mean(control_values)
  ) / sqrt(pooled_variance)
  
  correction <- 1 - (
    3 /
      (
        4 * (n_tle + n_control) -
          9
      )
  )
  
  correction * cohen_d
}

run_lodo <- function(
    score_vector,
    tle_ids,
    control_ids,
    endpoint
) {
  score_names <- names(score_vector)
  
  if (
    is.null(score_names) ||
    any(!nzchar(score_names))
  ) {
    stop(
      endpoint,
      "的评分向量缺少供者名称。",
      call. = FALSE
    )
  }
  
  score_vector <- as.numeric(score_vector)
  names(score_vector) <- score_names
  
  original_difference <- (
    mean(score_vector[tle_ids]) -
      mean(score_vector[control_ids])
  )
  
  original_direction <- sign(
    original_difference
  )
  
  all_ids <- c(
    tle_ids,
    control_ids
  )
  
  output <- vector(
    "list",
    length(all_ids)
  )
  
  for (index in seq_along(all_ids)) {
    removed_donor <- all_ids[index]
    
    remaining_tle <- setdiff(
      tle_ids,
      removed_donor
    )
    
    remaining_control <- setdiff(
      control_ids,
      removed_donor
    )
    
    reduced_difference <- (
      mean(score_vector[remaining_tle]) -
        mean(score_vector[remaining_control])
    )
    
    output[[index]] <- data.frame(
      endpoint = endpoint,
      removed_donor = removed_donor,
      
      removed_group = ifelse(
        removed_donor %in% tle_ids,
        "TLE",
        "Control"
      ),
      
      remaining_tle_n =
        length(remaining_tle),
      
      remaining_control_n =
        length(remaining_control),
      
      mean_difference =
        reduced_difference,
      
      direction_preserved = (
        sign(reduced_difference) ==
          original_direction
      ),
      
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(output)
}

analyse_endpoint <- function(
    endpoint,
    score_vector,
    tle_ids,
    control_ids
) {
  score_names <- names(score_vector)
  
  if (
    is.null(score_names) ||
    any(!nzchar(score_names))
  ) {
    stop(
      endpoint,
      "的评分向量缺少供者名称。",
      call. = FALSE
    )
  }
  
  score_vector <- as.numeric(score_vector)
  names(score_vector) <- score_names
  
  requested_ids <- c(
    tle_ids,
    control_ids
  )
  
  missing_ids <- setdiff(
    requested_ids,
    names(score_vector)
  )
  
  if (length(missing_ids) > 0L) {
    stop(
      endpoint,
      "缺少供者：",
      paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }
  
  tle_values <- score_vector[tle_ids]
  control_values <- score_vector[control_ids]
  
  if (
    any(!is.finite(tle_values)) ||
    any(!is.finite(control_values))
  ) {
    stop(
      endpoint,
      "包含非有限评分。",
      call. = FALSE
    )
  }
  
  permutation <- exact_permutation_test(
    tle_values,
    control_values
  )
  
  lodo_table <- run_lodo(
    score_vector = score_vector,
    tle_ids = tle_ids,
    control_ids = control_ids,
    endpoint = endpoint
  )
  
  mean_difference <- (
    mean(tle_values) -
      mean(control_values)
  )
  
  result <- data.frame(
    endpoint = endpoint,
    
    tle_n = length(tle_values),
    control_n = length(control_values),
    
    tle_mean = mean(tle_values),
    control_mean = mean(control_values),
    
    mean_difference =
      mean_difference,
    
    tle_median =
      median(tle_values),
    
    control_median =
      median(control_values),
    
    median_difference = (
      median(tle_values) -
        median(control_values)
    ),
    
    hedges_g =
      hedges_g(
        tle_values,
        control_values
      ),
    
    exact_permutation_p =
      permutation$p_value,
    
    permutation_count =
      permutation$permutation_count,
    
    lodo_direction_stability =
      mean(
        lodo_table$direction_preserved
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
        max(control_values)
    ) || (
      max(tle_values) <
        min(control_values)
    ),
    
    stringsAsFactors = FALSE
  )
  
  list(
    result = result,
    lodo = lodo_table
  )
}

###############################################################################
# 5. Read nuclear matrix
###############################################################################

message(
  "Reading GSE140393 RLD-normalized nuclear matrix."
)

connection <- gzfile(
  matrix_path,
  open = "rt"
)

input_table <- tryCatch(
  {
    read.delim(
      connection,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  },
  finally = {
    close(connection)
  }
)

message(
  "Input dimensions: ",
  nrow(input_table),
  " rows x ",
  ncol(input_table),
  " columns"
)

if (nrow(input_table) < 1000L) {
  stop(
    "输入矩阵基因行数异常。",
    call. = FALSE
  )
}

###############################################################################
# 6. Validate gene and sample columns
###############################################################################

gene_column <- "hgnc_symbol"

if (!gene_column %in% colnames(input_table)) {
  stop(
    "矩阵中没有hgnc_symbol列。\n实际列名：",
    paste(
      colnames(input_table),
      collapse = " | "
    ),
    call. = FALSE
  )
}

missing_sample_columns <- setdiff(
  sample_metadata$matrix_column,
  colnames(input_table)
)

if (length(missing_sample_columns) > 0L) {
  stop(
    "以下预设OLIG2列不存在：\n",
    paste(
      missing_sample_columns,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

genes <- canonicalize_gene(
  input_table[[gene_column]]
)

sample_data <- input_table[
  ,
  sample_metadata$matrix_column,
  drop = FALSE
]

numeric_columns <- lapply(
  sample_data,
  function(column) {
    suppressWarnings(
      as.numeric(
        as.character(column)
      )
    )
  }
)

nonfinite_counts <- vapply(
  numeric_columns,
  function(column) {
    sum(!is.finite(column))
  },
  integer(1)
)

conversion_audit <- data.frame(
  gsm = sample_metadata$gsm,
  donor_id = sample_metadata$donor_id,
  group = sample_metadata$group,
  matrix_column =
    sample_metadata$matrix_column,
  row_count =
    nrow(input_table),
  nonfinite_count =
    as.integer(nonfinite_counts),
  finite_fraction = (
    nrow(input_table) -
      nonfinite_counts
  ) / nrow(input_table),
  stringsAsFactors = FALSE
)

write_csv(
  conversion_audit,
  file.path(
    table_root,
    "GSE140393_FIXED_OLIG2_NUMERIC_AUDIT.csv"
  )
)

###############################################################################
# Handle non-finite RLD rows conservatively
###############################################################################

olig2_expression <- do.call(
  cbind,
  numeric_columns
)

storage.mode(olig2_expression) <-
  "double"

rownames(olig2_expression) <-
  genes

colnames(olig2_expression) <-
  sample_metadata$donor_id

finite_matrix <- is.finite(
  olig2_expression
)

finite_sample_count <- rowSums(
  finite_matrix
)

sample_count <- ncol(
  olig2_expression
)

all_nonfinite_rows <- (
  finite_sample_count == 0L
)

partially_nonfinite_rows <- (
  finite_sample_count > 0L &
    finite_sample_count < sample_count
)

nonfinite_row_audit <- data.frame(
  original_row_index =
    seq_len(nrow(olig2_expression)),
  
  gene_symbol =
    rownames(olig2_expression),
  
  finite_sample_count =
    finite_sample_count,
  
  nonfinite_sample_count = (
    sample_count -
      finite_sample_count
  ),
  
  nonfinite_donors = apply(
    !finite_matrix,
    MARGIN = 1,
    FUN = function(values) {
      paste(
        colnames(olig2_expression)[values],
        collapse = "|"
      )
    }
  ),
  
  classification = ifelse(
    all_nonfinite_rows,
    "ALL_SIX_NONFINITE",
    ifelse(
      partially_nonfinite_rows,
      "PARTIALLY_NONFINITE",
      "ALL_FINITE"
    )
  ),
  
  stringsAsFactors = FALSE
)

nonfinite_row_audit <- nonfinite_row_audit[
  nonfinite_row_audit$classification !=
    "ALL_FINITE",
  ,
  drop = FALSE
]

write_csv(
  nonfinite_row_audit,
  file.path(
    table_root,
    "GSE140393_FIXED_OLIG2_NONFINITE_ROW_AUDIT.csv"
  )
)

# Partial missingness cannot be repaired without making assumptions.
if (any(partially_nonfinite_rows)) {
  partial_audit <- nonfinite_row_audit[
    nonfinite_row_audit$classification ==
      "PARTIALLY_NONFINITE",
    ,
    drop = FALSE
  ]
  
  stop(
    paste0(
      "发现",
      nrow(partial_audit),
      "个仅在部分OLIG2样本中为NA/Inf的基因行。",
      "\n不能直接删除、填0或插补。",
      "\n请检查：",
      file.path(
        table_root,
        "GSE140393_FIXED_OLIG2_NONFINITE_ROW_AUDIT.csv"
      )
    ),
    call. = FALSE
  )
}

# Rows that are non-finite in all six samples contain no usable information.
all_nonfinite_count <- sum(
  all_nonfinite_rows
)

p5_audit_table <- suppressMessages(
  read_csv(
    p5_path,
    show_col_types = FALSE
  )
)

if ("gene" %in% colnames(p5_audit_table)) {
  frozen_p5_genes_for_audit <- unique(
    canonicalize_gene(
      p5_audit_table$gene
    )
  )
} else {
  frozen_p5_genes_for_audit <-
    character(0)
}

dropped_gene_symbols <- unique(
  rownames(olig2_expression)[
    all_nonfinite_rows
  ]
)

dropped_gene_symbols <- dropped_gene_symbols[
  nzchar(dropped_gene_symbols)
]

dropped_p5_genes <- intersect(
  dropped_gene_symbols,
  frozen_p5_genes_for_audit
)

nonfinite_summary <- data.frame(
  total_input_rows =
    nrow(olig2_expression),
  
  all_six_nonfinite_rows =
    all_nonfinite_count,
  
  partially_nonfinite_rows =
    sum(partially_nonfinite_rows),
  
  dropped_unique_gene_symbols =
    length(dropped_gene_symbols),
  
  dropped_frozen_p5_gene_count =
    length(dropped_p5_genes),
  
  dropped_frozen_p5_genes = paste(
    dropped_p5_genes,
    collapse = "|"
  ),
  
  stringsAsFactors = FALSE
)

write_csv(
  nonfinite_summary,
  file.path(
    table_root,
    "GSE140393_FIXED_OLIG2_NONFINITE_SUMMARY.csv"
  )
)

if (all_nonfinite_count > 0L) {
  message(
    "Removing ",
    all_nonfinite_count,
    " rows that are non-finite in all six OLIG2 samples."
  )
  
  message(
    "Frozen P5 genes among removed rows: ",
    ifelse(
      length(dropped_p5_genes) == 0L,
      "none",
      paste(
        dropped_p5_genes,
        collapse = ", "
      )
    )
  )
  
  olig2_expression <- olig2_expression[
    !all_nonfinite_rows,
    ,
    drop = FALSE
  ]
}

if (any(!is.finite(olig2_expression))) {
  stop(
    "删除全样本非有限行后，表达矩阵仍包含NA或Inf。",
    call. = FALSE
  )
}

valid_gene_rows <- (
  nzchar(rownames(olig2_expression)) &
    apply(
      olig2_expression,
      1,
      function(values) {
        all(is.finite(values))
      }
    )
)

olig2_expression <- olig2_expression[
  valid_gene_rows,
  ,
  drop = FALSE
]

olig2_expression <- collapse_duplicate_genes(
  olig2_expression
)

if (
  !identical(
    colnames(olig2_expression),
    sample_metadata$donor_id
  )
) {
  stop(
    "供者顺序与表达矩阵不一致。",
    call. = FALSE
  )
}

message(
  "Selected OLIG2 donors: ",
  paste(
    colnames(olig2_expression),
    collapse = ", "
  )
)

message(
  "Genes after symbol processing: ",
  nrow(olig2_expression)
)

write_csv(
  data.frame(
    gene_column = gene_column,
    expression_genes =
      nrow(olig2_expression),
    selected_samples =
      ncol(olig2_expression),
    stringsAsFactors = FALSE
  ),
  file.path(
    table_root,
    "GSE140393_FIXED_NUCLEAR_MATRIX_STRUCTURE.csv"
  )
)

###############################################################################
# 7. Read frozen P5
###############################################################################

p5_table <- suppressMessages(
  read_csv(
    p5_path,
    show_col_types = FALSE
  )
)

required_p5_columns <- c(
  "gene",
  "normalized_weight"
)

missing_p5_columns <- setdiff(
  required_p5_columns,
  colnames(p5_table)
)

if (length(missing_p5_columns) > 0L) {
  stop(
    "冻结P5文件缺少列：",
    paste(
      missing_p5_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

p5_table <- p5_table %>%
  transmute(
    gene = canonicalize_gene(gene),
    
    normalized_weight =
      suppressWarnings(
        as.numeric(
          normalized_weight
        )
      )
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

mapped_p5_genes <- intersect(
  p5_table$gene,
  rownames(olig2_expression)
)

if (length(mapped_p5_genes) < 50L) {
  stop(
    "可映射P5基因少于50个：",
    length(mapped_p5_genes),
    call. = FALSE
  )
}

p5_weights <- p5_table$normalized_weight[
  match(
    mapped_p5_genes,
    p5_table$gene
  )
]

p5_weights <- p5_weights /
  sum(p5_weights)

names(p5_weights) <-
  mapped_p5_genes

write_csv(
  p5_table %>%
    mutate(
      present_in_nuclear_matrix = (
        gene %in%
          rownames(olig2_expression)
      ),
      
      used_for_scoring = (
        gene %in%
          mapped_p5_genes
      )
    ),
  file.path(
    table_root,
    "GSE140393_FIXED_OLIG2_P5_GENE_MAPPING.csv"
  )
)

###############################################################################
# 8. Within-sample rank scoring
###############################################################################

rank_matrix <- apply(
  olig2_expression,
  2,
  function(values) {
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

rank_matrix <- as.matrix(
  rank_matrix
)

rownames(rank_matrix) <-
  rownames(olig2_expression)

colnames(rank_matrix) <-
  colnames(olig2_expression)

raw_weighted_rank <- colSums(
  sweep(
    rank_matrix[
      mapped_p5_genes,
      ,
      drop = FALSE
    ],
    MARGIN = 1,
    STATS = p5_weights,
    FUN = "*"
  )
)

equal_weight_rank <- colMeans(
  rank_matrix[
    mapped_p5_genes,
    ,
    drop = FALSE
  ]
)

###############################################################################
# 9. Expression-matched background calibration
###############################################################################

mean_expression <- rowMeans(
  olig2_expression,
  na.rm = TRUE
)

expression_order <- rank(
  mean_expression,
  ties.method = "average"
)

expression_bin <- ceiling(
  expression_order /
    length(expression_order) *
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
  rownames(olig2_expression)

p5_bins <- expression_bin[
  mapped_p5_genes
]

non_p5_genes <- setdiff(
  rownames(olig2_expression),
  mapped_p5_genes
)

background_replicates <- 500L

background_scores <- matrix(
  NA_real_,
  nrow = ncol(rank_matrix),
  ncol = background_replicates,
  dimnames = list(
    colnames(rank_matrix),
    paste0(
      "background_",
      seq_len(background_replicates)
    )
  )
)

set.seed(20260801)

for (
  replicate_index in
  seq_len(background_replicates)
) {
  sampled_background <- character(
    length(mapped_p5_genes)
  )
  
  for (
    bin_value in
    sort(unique(p5_bins))
  ) {
    p5_positions <- which(
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
      p5_positions
    ] <- sample(
      candidate_pool,
      size = length(p5_positions),
      replace = (
        length(candidate_pool) <
          length(p5_positions)
      )
    )
  }
  
  background_scores[
    ,
    replicate_index
  ] <- colSums(
    sweep(
      rank_matrix[
        sampled_background,
        ,
        drop = FALSE
      ],
      MARGIN = 1,
      STATS = p5_weights,
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
  1,
  stats::sd,
  na.rm = TRUE
)

background_sd[
  !is.finite(background_sd) |
    background_sd <= 0
] <- NA_real_

calibrated_rank_score <- (
  raw_weighted_rank -
    background_mean
)

p5_nes <- calibrated_rank_score /
  background_sd

if (any(!is.finite(p5_nes))) {
  stop(
    "背景校准后出现非有限P5 NES。",
    call. = FALSE
  )
}

###############################################################################
# 10. Donor scores
###############################################################################

score_table <- sample_metadata %>%
  mutate(
    mapped_p5_genes =
      length(mapped_p5_genes),
    
    raw_weighted_rank =
      raw_weighted_rank[donor_id],
    
    equal_weight_rank =
      equal_weight_rank[donor_id],
    
    matched_background_mean =
      background_mean[donor_id],
    
    matched_background_sd =
      background_sd[donor_id],
    
    calibrated_rank_score =
      calibrated_rank_score[donor_id],
    
    p5_nes =
      p5_nes[donor_id]
  )

write_csv(
  score_table,
  file.path(
    table_root,
    "GSE140393_FIXED_OLIG2_P5_DONOR_SCORES.csv"
  )
)

###############################################################################
# 11. Donor-level statistics
###############################################################################

tle_ids <- sample_metadata$donor_id[
  sample_metadata$group == "TLE"
]

control_ids <- sample_metadata$donor_id[
  sample_metadata$group == "Control"
]

raw_analysis <- analyse_endpoint(
  endpoint =
    "FROZEN_WEIGHTED_RANK",
  score_vector =
    raw_weighted_rank,
  tle_ids =
    tle_ids,
  control_ids =
    control_ids
)

equal_analysis <- analyse_endpoint(
  endpoint =
    "EQUAL_WEIGHT_RANK",
  score_vector =
    equal_weight_rank,
  tle_ids =
    tle_ids,
  control_ids =
    control_ids
)

nes_analysis <- analyse_endpoint(
  endpoint =
    "MATCHED_BACKGROUND_P5_NES",
  score_vector =
    p5_nes,
  tle_ids =
    tle_ids,
  control_ids =
    control_ids
)

result_table <- bind_rows(
  raw_analysis$result,
  equal_analysis$result,
  nes_analysis$result
)

result_table$p_adjust_BH <- p.adjust(
  result_table$exact_permutation_p,
  method = "BH"
)

lodo_table <- bind_rows(
  raw_analysis$lodo,
  equal_analysis$lodo,
  nes_analysis$lodo
)

write_csv(
  result_table,
  file.path(
    table_root,
    "GSE140393_FIXED_OLIG2_P5_COMPARISON_RESULTS.csv"
  )
)

write_csv(
  lodo_table,
  file.path(
    table_root,
    "GSE140393_FIXED_OLIG2_P5_LODO_RESULTS.csv"
  )
)

###############################################################################
# 12. Per-gene contribution
###############################################################################

tle_mean_rank <- rowMeans(
  rank_matrix[
    mapped_p5_genes,
    tle_ids,
    drop = FALSE
  ]
)

control_mean_rank <- rowMeans(
  rank_matrix[
    mapped_p5_genes,
    control_ids,
    drop = FALSE
  ]
)

gene_contributions <- data.frame(
  gene =
    mapped_p5_genes,
  
  normalized_weight =
    p5_weights[
      mapped_p5_genes
    ],
  
  tle_mean_rank =
    tle_mean_rank,
  
  control_mean_rank =
    control_mean_rank,
  
  rank_difference = (
    tle_mean_rank -
      control_mean_rank
  ),
  
  weighted_contribution = (
    tle_mean_rank -
      control_mean_rank
  ) *
    p5_weights[
      mapped_p5_genes
    ],
  
  stringsAsFactors = FALSE
) %>%
  arrange(
    desc(
      abs(weighted_contribution)
    )
  )

write_csv(
  gene_contributions,
  file.path(
    table_root,
    "GSE140393_FIXED_OLIG2_P5_GENE_CONTRIBUTIONS.csv"
  )
)

###############################################################################
# 13. Figure
###############################################################################

plot_data <- score_table %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "Control",
        "TLE"
      )
    )
  )

score_plot <- ggplot(
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
    width = 0.06,
    height = 0,
    size = 3.2
  ) +
  geom_text(
    aes(
      label = donor_id
    ),
    position = position_jitter(
      width = 0.06,
      height = 0,
      seed = 20260801
    ),
    vjust = -0.85,
    size = 3.2,
    show.legend = FALSE
  ) +
  labs(
    title = paste0(
      "Frozen P5 in GSE140393 ",
      "OLIG2-positive nuclei"
    ),
    subtitle = paste0(
      "Three TLE versus three temporal-lobe controls; ",
      "donor is the statistical unit"
    ),
    x = NULL,
    y = "Matched-background frozen P5 NES"
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
    "GSE140393_FIXED_OLIG2_P5_DONOR_SCORES.png"
  ),
  plot = score_plot,
  width = 7,
  height = 6,
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = file.path(
    figure_root,
    "GSE140393_FIXED_OLIG2_P5_DONOR_SCORES.pdf"
  ),
  plot = score_plot,
  width = 7,
  height = 6,
  bg = "white"
)

###############################################################################
# 14. Summary and session information
###############################################################################

primary_result <- result_table %>%
  filter(
    endpoint ==
      "MATCHED_BACKGROUND_P5_NES"
  )

if (nrow(primary_result) != 1L) {
  stop(
    "无法唯一识别主要P5 NES结果。",
    call. = FALSE
  )
}

summary_lines <- c(
  "# GSE140393 OLIG2-positive nuclear bridge — fixed mapping",
  "",
  paste0(
    "- Matrix: ",
    matrix_path
  ),
  paste0(
    "- Gene column: ",
    gene_column
  ),
  paste0(
    "- Expression genes: ",
    nrow(olig2_expression)
  ),
  paste0(
    "- Mapped frozen P5 genes: ",
    length(mapped_p5_genes)
  ),
  paste0(
    "- TLE donors: ",
    paste(
      tle_ids,
      collapse = "|"
    )
  ),
  paste0(
    "- Control donors: ",
    paste(
      control_ids,
      collapse = "|"
    )
  ),
  "",
  "## Primary matched-background result",
  "",
  paste0(
    "- TLE mean P5 NES: ",
    signif(
      primary_result$tle_mean,
      6
    )
  ),
  paste0(
    "- Control mean P5 NES: ",
    signif(
      primary_result$control_mean,
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
  paste0(
    "- Complete separation: ",
    primary_result$complete_separation
  ),
  "",
  "## Interpretation boundary",
  "",
  paste0(
    "- This is a donor-level comparison within one ",
    "processed nuclear RNA-seq matrix."
  ),
  paste0(
    "- OLIG2-positive nuclei are OPC-enriched and do not ",
    "represent purified mature oligodendrocytes."
  ),
  paste0(
    "- TLE tissue was surgical whereas controls were ",
    "postmortem; tissue-source confounding remains."
  )
)

writeLines(
  summary_lines,
  file.path(
    output_root,
    "GSE140393_FIXED_OLIG2_NUCLEAR_BRIDGE_SUMMARY.md"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    log_root,
    "sessionInfo.txt"
  )
)

###############################################################################
# 15. Console output
###############################################################################

message("")
message("======================================================")
message("GSE140393 fixed OLIG2 nuclear bridge completed")
message("======================================================")
message("")
message("Expression genes: ", nrow(olig2_expression))
message("Mapped P5 genes: ", length(mapped_p5_genes))
message("")
message("Fixed sample mapping:")

print(
  as.data.frame(
    sample_metadata %>%
      select(
        gsm,
        donor_id,
        group,
        matrix_column
      )
  ),
  row.names = FALSE
)

message("")
message("Donor scores:")

print(
  as.data.frame(
    score_table %>%
      select(
        donor_id,
        group,
        raw_weighted_rank,
        equal_weight_rank,
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
    result_table %>%
      select(
        endpoint,
        tle_n,
        control_n,
        mean_difference,
        hedges_g,
        exact_permutation_p,
        p_adjust_BH,
        lodo_direction_stability,
        direction,
        complete_separation
      )
  ),
  row.names = FALSE
)

message("")
message("Output directory:")
message(output_root)
# P5_RESTART5_GSE295237_SPATIAL_QUALIFICATION_FIX2.R
# Pure-bioinformatics manuscript: Oligodendrocytes__P5
#
# Goal
# ----
# Spatial/reference qualification of the frozen PRIMARY mechanism module:
#   APOE -> LRP1
# in an independent human temporal-neocortex epilepsy spatial dataset.
#
# Dataset:
#   GSE295237 (10x Visium)
#   Spatial samples: Control_2, Control_3, FCDI_1, FCDI_2
#
# Important boundary:
#   GSE295237 is FCD type I-related epilepsy, NOT TLE.
#   Therefore this is cross-epilepsy temporal-neocortex spatial qualification,
#   not independent TLE replication.
#
# Frozen upstream:
#   - P5 Top100/loading from RESTART2
#   - mechanism axis freeze from RESTART4B
#   - lead module: APOE -> LRP1
#
# Main spatial question:
#   In oligodendrocyte-rich Visium spots, is the OL-abundance-adjusted P5 signal
#   spatially compatible with:
#     (a) OL/LRP1 receiver context, and
#     (b) neighboring APOE-rich astrocytic or myeloid context?
#
# Spot-level correlations are descriptive. Sample/donor is the biological unit.
# With 2 FCDI and 2 controls, no formal disease-group p-value is authorized.

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(20261006)

required_pkgs <- c("data.table", "Matrix")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs)) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nUse the established R environment; do not update packages only for RESTART5."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

ROOT <- Sys.getenv("P5_WORK_ROOT")

PROJECT <- file.path(
  ROOT,
  "P5_MATURE_OLIGODENDROCYTE_PURE_BIOINFO_RESTART"
)

R2 <- file.path(PROJECT, "RESTART2_FUNCTIONAL_ARCHITECTURE")
OUT <- file.path(PROJECT, "RESTART5_GSE295237_SPATIAL_QUALIFICATION")
RAW <- file.path(OUT, "00_RAW_GSE295237")
EXTRACT <- file.path(RAW, "extracted")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(RAW, recursive = TRUE, showWarnings = FALSE)
dir.create(EXTRACT, recursive = TRUE, showWarnings = FALSE)

P5_FILE <- file.path(R2, "P5_RESTART2_TOP100_WITH_LOADINGS.tsv")

if (!file.exists(P5_FILE)) {
  stop("Frozen P5 Top100/loading file not found: ", P5_FILE)
}

# -------------------------------------------------------------------------
# 1. Frozen signatures / primary axis
# -------------------------------------------------------------------------

p5 <- fread(P5_FILE)

if (!all(c("gene", "meta_program_score") %in% names(p5))) {
  stop("P5 Top100/loading table lacks gene/meta_program_score.")
}

p5[, gene := toupper(trimws(as.character(gene)))]
p5[, meta_program_score := as.numeric(meta_program_score)]
p5 <- p5[
  !is.na(gene) & nzchar(gene) &
    is.finite(meta_program_score)
]

if (uniqueN(p5$gene) != 100L) {
  stop("Expected exactly 100 frozen P5 genes.")
}

p5[, weight := pmax(meta_program_score, 0)]
if (sum(p5$weight) <= 0) {
  p5[, weight := 1]
}
p5[, weight := weight / sum(weight)]

OL_MARKERS <- c(
  "PLP1","MBP","MOG","MOBP","MAG","CNP",
  "CLDN11","ERMN","SOX10","OLIG1","OLIG2"
)

ASTRO_MARKERS <- c(
  "AQP4","GFAP","SLC1A2","SLC1A3","ALDH1L1","FGFR3","GLUL"
)

MYELOID_MARKERS <- c(
  "P2RY12","CX3CR1","CSF1R","C1QA","C1QB","TYROBP","AIF1"
)

PRIMARY_LIGAND <- "APOE"
PRIMARY_RECEPTOR <- "LRP1"

# -------------------------------------------------------------------------
# 2. Download / extract public processed Space Ranger files
# -------------------------------------------------------------------------

TAR_FILE <- file.path(RAW, "GSE295237_RAW.tar")

# Large GEO supplement: R's default network timeout is often only ~60 s.
# GSE295237_RAW.tar is ~204.8 MB, so use a long timeout and bounded retries.
options(
  timeout = max(
    getOption("timeout"),
    7200
  )
)

EXPECTED_MIN_TAR_BYTES <- 190e6

if (
  !file.exists(TAR_FILE) ||
  !is.finite(file.info(TAR_FILE)$size) ||
  file.info(TAR_FILE)$size < EXPECTED_MIN_TAR_BYTES
) {

  message(
    "[RESTART5] 1/5 Downloading GSE295237_RAW.tar (~195 MB). ",
    "A slow connection may take several minutes."
  )

  urls <- c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295237/suppl/GSE295237_RAW.tar",
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE295237&format=file&file=GSE295237_RAW.tar"
  )

  success <- FALSE
  tmp <- paste0(TAR_FILE, ".part")

  # Do not mistake an old incomplete .part file for a complete archive.
  if (
    file.exists(tmp) &&
    is.finite(file.info(tmp)$size) &&
    file.info(tmp)$size < EXPECTED_MIN_TAR_BYTES
  ) {
    unlink(tmp)
  }

  for (attempt in seq_len(3L)) {

    u <- urls[
      ((attempt - 1L) %% length(urls)) + 1L
    ]

    message(
      "[RESTART5]   download attempt ",
      attempt,
      "/3"
    )

    if (file.exists(tmp)) unlink(tmp)

    status <- tryCatch(
      suppressWarnings(
        utils::download.file(
          url = u,
          destfile = tmp,
          mode = "wb",
          method = "libcurl",
          quiet = FALSE
        )
      ),
      error = function(e) {
        message(
          "[RESTART5]   download error: ",
          conditionMessage(e)
        )
        NA_integer_
      }
    )

    downloaded_bytes <- if (
      file.exists(tmp) &&
      is.finite(file.info(tmp)$size)
    ) {
      as.numeric(file.info(tmp)$size)
    } else {
      0
    }

    message(
      "[RESTART5]   downloaded ",
      round(downloaded_bytes / 1024^2, 1),
      " MB"
    )

    if (
      !is.na(status) &&
      identical(as.integer(status), 0L) &&
      downloaded_bytes >= EXPECTED_MIN_TAR_BYTES
    ) {

      if (file.exists(TAR_FILE)) unlink(TAR_FILE)

      if (!file.rename(tmp, TAR_FILE)) {
        stop(
          "Download completed but the temporary archive could not be renamed."
        )
      }

      success <- TRUE
      break
    }

    if (attempt < 3L) {
      message("[RESTART5]   retrying in 5 seconds...")
      Sys.sleep(5)
    }
  }

  if (!success) {
    stop(
      "GSE295237_RAW.tar did not complete after 3 long-timeout attempts.\n",
      "This is a network/download problem, not an analysis failure.\n",
      "If needed, download GSE295237_RAW.tar manually and place it at:\n",
      TAR_FILE,
      "\nThen rerun the same script; it will skip the download."
    )
  }
}

if (length(list.files(EXTRACT, recursive = TRUE, full.names = TRUE)) < 10L) {
  message("[RESTART5] 2/5 Extracting GSE295237 processed spatial files...")
  utils::untar(TAR_FILE, exdir = EXTRACT)
}

all_files <- list.files(
  EXTRACT,
  recursive = TRUE,
  full.names = TRUE
)

# -------------------------------------------------------------------------
# 3. Sample registry / robust Space Ranger parser
# -------------------------------------------------------------------------

SAMPLES <- data.table(
  gsm = c(
    "GSM8944108",
    "GSM8944109",
    "GSM8944110",
    "GSM8944111"
  ),
  sample = c(
    "Control_2",
    "Control_3",
    "FCDI_1",
    "FCDI_2"
  ),
  group = c(
    "Control",
    "Control",
    "FCDI",
    "FCDI"
  )
)

find_one_sample_file <- function(gsm, pattern) {

  h <- all_files[
    grepl(gsm, basename(all_files), fixed = TRUE) &
      grepl(pattern, basename(all_files), ignore.case = TRUE)
  ]

  if (length(h) != 1L) {
    stop(
      "Could not uniquely identify file for ",
      gsm, " pattern=", pattern,
      "\nHits:\n", paste(h, collapse = "\n")
    )
  }

  h
}

read_text_auto <- function(path, header = "auto") {

  if (grepl("\\.gz$", path, ignore.case = TRUE)) {

    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)

    txt <- readLines(
      con,
      warn = FALSE
    )

    data.table::fread(
      text = paste(
        txt,
        collapse = "\n"
      ),
      header = header
    )

  } else {

    data.table::fread(
      path,
      header = header
    )
  }
}

read_mtx_auto <- function(path) {

  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    Matrix::readMM(con)
  } else {
    Matrix::readMM(path)
  }
}

read_positions <- function(path) {

  z <- tryCatch(
    read_text_auto(path, header = TRUE),
    error = function(e) NULL
  )

  canonical <- c(
    "barcode","in_tissue","array_row","array_col",
    "pxl_row_in_fullres","pxl_col_in_fullres"
  )

  if (
    is.null(z) ||
    !all(c("barcode","in_tissue","array_row","array_col") %in% names(z))
  ) {

    z <- read_text_auto(path, header = FALSE)

    if (ncol(z) < 6L) {
      stop("Unrecognized tissue_positions format: ", path)
    }

    setnames(
      z,
      names(z)[1:6],
      canonical
    )
  }

  z[, barcode := as.character(barcode)]
  z[, in_tissue := as.integer(in_tissue)]
  z[, array_row := as.integer(array_row)]
  z[, array_col := as.integer(array_col)]

  z
}

collapse_duplicate_genes <- function(counts, genes) {

  genes <- toupper(trimws(as.character(genes)))

  keep <- !is.na(genes) & nzchar(genes)
  counts <- counts[keep, , drop = FALSE]
  genes <- genes[keep]

  if (!anyDuplicated(genes)) {
    rownames(counts) <- genes
    return(counts)
  }

  lev <- unique(genes)
  map <- match(genes, lev)

  M <- sparseMatrix(
    i = map,
    j = seq_along(map),
    x = 1,
    dims = c(length(lev), length(map))
  )

  out <- M %*% counts
  out <- as(out, "dgCMatrix")
  rownames(out) <- lev
  out
}

zscore_rows <- function(mat) {

  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 1, sd, na.rm = TRUE)

  good <- is.finite(sdv) & sdv > 0

  z <- matrix(
    0,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )

  z[good, ] <- (
    mat[good, , drop = FALSE] -
      mu[good]
  ) / sdv[good]

  z
}

score_mean <- function(z, genes, min_genes = 3L) {

  g <- intersect(genes, rownames(z))

  if (length(g) < min_genes) {
    return(rep(NA_real_, ncol(z)))
  }

  colMeans(
    z[g, , drop = FALSE],
    na.rm = TRUE
  )
}

neighbor_mean_visium <- function(v, row, col) {

  offsets <- rbind(
    c(0,  2),
    c(0, -2),
    c(1,  1),
    c(1, -1),
    c(-1, 1),
    c(-1,-1)
  )

  key <- paste(row, col, sep = "_")
  lookup <- setNames(seq_along(key), key)

  acc <- rep(0, length(v))
  n <- integer(length(v))

  for (k in seq_len(nrow(offsets))) {

    target <- paste(
      row + offsets[k, 1],
      col + offsets[k, 2],
      sep = "_"
    )

    idx <- unname(lookup[target])
    ok <- !is.na(idx)

    acc[ok] <- acc[ok] + v[idx[ok]]
    n[ok] <- n[ok] + 1L
  }

  out <- rep(NA_real_, length(v))
  ok <- n > 0
  out[ok] <- acc[ok] / n[ok]
  out
}

safe_spearman <- function(x, y) {

  ok <- is.finite(x) & is.finite(y)

  if (sum(ok) < 30L) return(NA_real_)

  suppressWarnings(
    cor(
      x[ok],
      y[ok],
      method = "spearman"
    )
  )
}

# -------------------------------------------------------------------------
# 4. Score each spatial sample
# -------------------------------------------------------------------------

message("[RESTART5] 3/5 Scoring P5 and APOE-LRP1 spatial context...")

sample_stats <- list()
spot_outputs <- list()

for (ii in seq_len(nrow(SAMPLES))) {

  sm <- SAMPLES[ii]
  gsm <- sm$gsm

  message("[RESTART5]   ", gsm, " ", sm$sample)

  mtx_file <- find_one_sample_file(
    gsm,
    "matrix\\.mtx(\\.gz)?$"
  )

  features_file <- find_one_sample_file(
    gsm,
    "features\\.tsv(\\.gz)?$"
  )

  barcodes_file <- find_one_sample_file(
    gsm,
    "barcodes\\.tsv(\\.gz)?$"
  )

  positions_file <- find_one_sample_file(
    gsm,
    "tissue_positions.*\\.csv(\\.gz)?$"
  )

  counts <- read_mtx_auto(mtx_file)
  counts <- as(counts, "dgCMatrix")

  feat <- read_text_auto(
    features_file,
    header = FALSE
  )

  bc <- read_text_auto(
    barcodes_file,
    header = FALSE
  )[[1]]

  if (nrow(feat) != nrow(counts) && nrow(feat) == ncol(counts)) {
    counts <- t(counts)
    counts <- as(counts, "dgCMatrix")
  }

  if (nrow(feat) != nrow(counts) || length(bc) != ncol(counts)) {
    stop("10x dimensions do not match for ", gsm)
  }

  gene_symbol <- as.character(
    feat[[min(2L, ncol(feat))]]
  )

  colnames(counts) <- as.character(bc)
  counts <- collapse_duplicate_genes(
    counts,
    gene_symbol
  )

  pos <- read_positions(positions_file)
  pos <- pos[
    in_tissue == 1L
  ]

  common <- intersect(
    colnames(counts),
    pos$barcode
  )

  if (length(common) < 500L) {
    stop(
      "Too few in-tissue spots for ",
      gsm, ": ", length(common)
    )
  }

  counts <- counts[
    ,
    common,
    drop = FALSE
  ]

  pos <- pos[
    match(common, barcode)
  ]

  lib <- Matrix::colSums(counts)

  keep_spot <- is.finite(lib) & lib > 0
  counts <- counts[, keep_spot, drop = FALSE]
  lib <- lib[keep_spot]
  pos <- pos[keep_spot]

  wanted <- unique(c(
    p5$gene,
    OL_MARKERS,
    ASTRO_MARKERS,
    MYELOID_MARKERS,
    PRIMARY_LIGAND,
    PRIMARY_RECEPTOR
  ))

  present <- intersect(
    wanted,
    rownames(counts)
  )

  small_counts <- as.matrix(
    counts[
      present,
      ,
      drop = FALSE
    ]
  )

  logcpm <- log1p(
    sweep(
      small_counts,
      2,
      lib,
      "/"
    ) * 1e6
  )

  z <- zscore_rows(logcpm)

  # Frozen weighted P5 score.
  p5_present <- intersect(
    p5$gene,
    rownames(z)
  )

  if (length(p5_present) < 60L) {
    stop(
      "P5 coverage <60/100 for ",
      gsm, ": ", length(p5_present)
    )
  }

  pw <- p5[
    match(p5_present, gene),
    weight
  ]
  pw <- pw / sum(pw)

  p5_score <- as.numeric(
    crossprod(
      pw,
      z[p5_present, , drop = FALSE]
    )
  )

  ol_score <- score_mean(
    z,
    OL_MARKERS,
    min_genes = 5L
  )

  astro_score <- score_mean(
    z,
    ASTRO_MARKERS,
    min_genes = 4L
  )

  myeloid_score <- score_mean(
    z,
    MYELOID_MARKERS,
    min_genes = 4L
  )

  if (
    any(!is.finite(ol_score)) ||
    any(!is.finite(astro_score)) ||
    any(!is.finite(myeloid_score))
  ) {
    stop("Lineage marker score failed for ", gsm)
  }

  # Remove OL-abundance component from P5 before spatial state analysis.
  p5_resid <- residuals(
    lm(
      p5_score ~ ol_score + log10(lib)
    )
  )

  apoe_z <- if (PRIMARY_LIGAND %in% rownames(z)) {
    as.numeric(z[PRIMARY_LIGAND, ])
  } else {
    rep(NA_real_, ncol(z))
  }

  lrp1_z <- if (PRIMARY_RECEPTOR %in% rownames(z)) {
    as.numeric(z[PRIMARY_RECEPTOR, ])
  } else {
    rep(NA_real_, ncol(z))
  }

  if (
    all(!is.finite(apoe_z)) ||
    all(!is.finite(lrp1_z))
  ) {
    stop("APOE or LRP1 is not detected in ", gsm)
  }

  astro_apoe_sender <- rowMeans(
    cbind(apoe_z, astro_score),
    na.rm = TRUE
  )

  myeloid_apoe_sender <- rowMeans(
    cbind(apoe_z, myeloid_score),
    na.rm = TRUE
  )

  ol_lrp1_receiver <- rowMeans(
    cbind(lrp1_z, ol_score),
    na.rm = TRUE
  )

  n_astro_apoe <- neighbor_mean_visium(
    astro_apoe_sender,
    pos$array_row,
    pos$array_col
  )

  n_myeloid_apoe <- neighbor_mean_visium(
    myeloid_apoe_sender,
    pos$array_row,
    pos$array_col
  )

  ol_cut <- quantile(
    ol_score,
    probs = 0.75,
    na.rm = TRUE,
    names = FALSE
  )

  ol_rich <- ol_score >= ol_cut

  if (sum(ol_rich) < 100L) {
    stop("Too few OL-rich spots for ", gsm)
  }

  rho_receiver <- safe_spearman(
    p5_resid[ol_rich],
    ol_lrp1_receiver[ol_rich]
  )

  rho_astro_neighbor <- safe_spearman(
    p5_resid[ol_rich],
    n_astro_apoe[ol_rich]
  )

  rho_myeloid_neighbor <- safe_spearman(
    p5_resid[ol_rich],
    n_myeloid_apoe[ol_rich]
  )

  sample_stats[[ii]] <- data.table(
    gsm = gsm,
    sample = sm$sample,
    group = sm$group,
    n_spots = ncol(counts),
    P5_coverage = length(p5_present),
    n_OL_markers = length(intersect(OL_MARKERS, rownames(z))),
    n_astro_markers = length(intersect(ASTRO_MARKERS, rownames(z))),
    n_myeloid_markers = length(intersect(MYELOID_MARKERS, rownames(z))),
    n_OL_rich_spots = sum(ol_rich),
    rho_P5resid_vs_OL_LRP1_receiver = rho_receiver,
    rho_P5resid_vs_neighbor_Astro_APOE = rho_astro_neighbor,
    rho_P5resid_vs_neighbor_Myeloid_APOE = rho_myeloid_neighbor,
    max_neighbor_APOE_rho = max(
      c(rho_astro_neighbor, rho_myeloid_neighbor),
      na.rm = TRUE
    ),
    dominant_neighbor_sender = ifelse(
      rho_astro_neighbor >= rho_myeloid_neighbor,
      "Astrocyte_APOE",
      "Myeloid_APOE"
    )
  )

  spot_outputs[[ii]] <- data.table(
    gsm = gsm,
    sample = sm$sample,
    group = sm$group,
    barcode = colnames(counts),
    array_row = pos$array_row,
    array_col = pos$array_col,
    pxl_row_in_fullres = pos$pxl_row_in_fullres,
    pxl_col_in_fullres = pos$pxl_col_in_fullres,
    library_size = lib,
    P5_score = p5_score,
    OL_identity = ol_score,
    P5_residual = p5_resid,
    APOE_z = apoe_z,
    LRP1_z = lrp1_z,
    Astro_identity = astro_score,
    Myeloid_identity = myeloid_score,
    Astro_APOE_sender = astro_apoe_sender,
    Myeloid_APOE_sender = myeloid_apoe_sender,
    neighbor_Astro_APOE = n_astro_apoe,
    neighbor_Myeloid_APOE = n_myeloid_apoe,
    OL_LRP1_receiver = ol_lrp1_receiver,
    OL_rich = ol_rich
  )
}

stats <- rbindlist(sample_stats)
spots <- rbindlist(spot_outputs)

fwrite(
  stats,
  file.path(
    OUT,
    "P5_RESTART5_SAMPLE_SPATIAL_STATISTICS.tsv"
  ),
  sep = "\t"
)

fwrite(
  spots,
  file.path(
    OUT,
    "P5_RESTART5_SPOT_SCORES.tsv"
  ),
  sep = "\t"
)

# -------------------------------------------------------------------------
# 5. Sample-level synthesis + simple spatial figures
# -------------------------------------------------------------------------

message("[RESTART5] 4/5 Summarizing donor-level spatial direction...")

summary <- stats[
  ,
  .(
    n_samples = .N,
    receiver_positive = sum(
      rho_P5resid_vs_OL_LRP1_receiver > 0,
      na.rm = TRUE
    ),
    neighbor_positive = sum(
      max_neighbor_APOE_rho > 0,
      na.rm = TRUE
    ),
    median_receiver_rho = median(
      rho_P5resid_vs_OL_LRP1_receiver,
      na.rm = TRUE
    ),
    median_neighbor_rho = median(
      max_neighbor_APOE_rho,
      na.rm = TRUE
    )
  ),
  by = group
]

fwrite(
  summary,
  file.path(
    OUT,
    "P5_RESTART5_GROUP_DIRECTION_SUMMARY.tsv"
  ),
  sep = "\t"
)

# Structural support across all four donors.
all_receiver_positive <- sum(
  stats$rho_P5resid_vs_OL_LRP1_receiver > 0,
  na.rm = TRUE
)

all_neighbor_positive <- sum(
  stats$max_neighbor_APOE_rho > 0,
  na.rm = TRUE
)

FINAL_CLASS <- if (
  all_receiver_positive >= 3L &&
  all_neighbor_positive >= 3L
) {
  "SPATIAL_REFERENCE_SUPPORT_FOR_APOE_LRP1_P5_COMPATIBILITY"
} else if (
  all_receiver_positive >= 3L ||
  all_neighbor_positive >= 3L
) {
  "PARTIAL_SPATIAL_REFERENCE_SUPPORT"
} else {
  "NO_CONSISTENT_SPATIAL_REFERENCE_SUPPORT"
}

# Compact base-R maps.
FIG_DIR <- file.path(OUT, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

plot_spatial <- function(z, value, title, file) {

  v <- z[[value]]

  png(
    file,
    width = 1800,
    height = 1500,
    res = 220
  )

  pal <- grDevices::hcl.colors(
    100,
    "Viridis"
  )

  br <- quantile(
    v[is.finite(v)],
    probs = seq(0, 1, length.out = 101),
    na.rm = TRUE
  )

  br <- unique(br)

  if (length(br) < 3L) {
    cols <- rep(pal[50], length(v))
  } else {
    idx <- cut(
      v,
      breaks = br,
      include.lowest = TRUE,
      labels = FALSE
    )
    idx <- pmin(
      pmax(idx, 1),
      length(pal)
    )
    cols <- pal[idx]
  }

  plot(
    z$pxl_col_in_fullres,
    -z$pxl_row_in_fullres,
    pch = 16,
    cex = 0.45,
    col = cols,
    asp = 1,
    xlab = "",
    ylab = "",
    axes = FALSE,
    main = title
  )

  dev.off()
}

# data.table-safe loop
for (g in unique(spots$gsm)) {

  z <- spots[spots$gsm == g]
  sm <- unique(z$sample)
  gp <- unique(z$group)

  plot_spatial(
    z,
    "P5_residual",
    paste0(
      g, " ", sm, " [", gp, "]\n",
      "OL-adjusted P5 residual"
    ),
    file.path(
      FIG_DIR,
      paste0(g, "_P5_residual.png")
    )
  )

  plot_spatial(
    z,
    "OL_LRP1_receiver",
    paste0(
      g, " ", sm, " [", gp, "]\n",
      "OL/LRP1 receiver context"
    ),
    file.path(
      FIG_DIR,
      paste0(g, "_LRP1_receiver.png")
    )
  )
}

message("[RESTART5] 5/5 Freezing interpretation...")

decision <- data.table(
  item = c(
    "dataset",
    "spatial_disease_context",
    "primary_axis",
    "n_samples",
    "n_FCDI",
    "n_Control",
    "receiver_positive_samples",
    "neighbor_APOE_positive_samples",
    "FINAL_CLASS",
    "CLAIM_BOUNDARY"
  ),
  value = c(
    "GSE295237",
    "FCDI_related_epilepsy_temporal_neocortex_NOT_TLE",
    "APOE_LRP1",
    nrow(stats),
    sum(stats$group == "FCDI"),
    sum(stats$group == "Control"),
    all_receiver_positive,
    all_neighbor_positive,
    FINAL_CLASS,
    "cross_epilepsy_spatial_reference_only_not_TLE_replication"
  )
)

fwrite(
  decision,
  file.path(
    OUT,
    "P5_RESTART5_DECISION.tsv"
  ),
  sep = "\t"
)

writeLines(
  c(
    "P5 MATURE OLIGODENDROCYTE PURE-BIOINFORMATICS MANUSCRIPT",
    "RESTART-5 | GSE295237 SPATIAL QUALIFICATION",
    "============================================",
    "",
    "DATASET",
    "-------",
    "- GSE295237 10x Visium temporal neocortex.",
    "- Spatial samples: Control_2, Control_3, FCDI_1, FCDI_2.",
    "- FCDI is epilepsy-related cortical pathology but is NOT TLE.",
    "- Therefore this stage is spatial/reference qualification, not TLE replication.",
    "",
    "FROZEN PRIMARY MECHANISM",
    "------------------------",
    "- APOE -> LRP1",
    "- Sender context: Astrocytes and Microglia/macrophage.",
    "- Receiver context: oligodendrocyte-rich spots.",
    "- P5 is adjusted for oligodendrocyte identity and library size before spatial state analysis.",
    "",
    "SAMPLE RESULTS",
    "--------------",
    capture.output(print(stats)),
    "",
    "GROUP DIRECTION SUMMARY",
    "-----------------------",
    capture.output(print(summary)),
    "",
    paste0("FINAL_CLASS=", FINAL_CLASS),
    "",
    "INTERPRETATION BOUNDARY",
    "-----------------------",
    "- Visium spots are multicellular and do not prove direct cell-cell signaling.",
    "- APOE source assignment is approximated by neighboring astrocyte/myeloid marker context.",
    "- Spot-level correlations are descriptive; sample/donor is the biological unit.",
    "- n=2 FCDI and n=2 controls is insufficient for a formal disease-group test.",
    "- Positive results support anatomical compatibility of the frozen APOE-LRP1-P5 model only.",
    "- No secretion, receptor activation, physical binding or causal claim is authorized."
  ),
  file.path(
    OUT,
    "P5_RESTART5_REPORT.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUT,
    "P5_RESTART5_SESSIONINFO.txt"
  )
)

message("P5 RESTART-5 completed.")
message("Output: ", OUT)

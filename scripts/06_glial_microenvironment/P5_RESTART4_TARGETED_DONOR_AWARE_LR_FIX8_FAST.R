# P5_RESTART4_TARGETED_DONOR_AWARE_LR_FIX8_FAST.R
# Uses the exact primary cNMF H5AD inputs rather than the blinded CP1R object.
# This avoids invalid direct matching between CP1R B01-B08 and the donor_id
# carried by the frozen cNMF/P5 donor-score tables.

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(20261005)

required_pkgs <- c("data.table", "edgeR")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs)) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

ROOT <- Sys.getenv("P5_WORK_ROOT")
PROJECT <- file.path(ROOT, "P5_MATURE_OLIGODENDROCYTE_PURE_BIOINFO_RESTART")
R2 <- file.path(PROJECT, "RESTART2_FUNCTIONAL_ARCHITECTURE")
R3 <- file.path(PROJECT, "RESTART3_REGULATOR_MICROENVIRONMENT_PRIORITIZATION")
OUT <- file.path(PROJECT, "RESTART4_TARGETED_DONOR_AWARE_LR")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

P5 <- "Oligodendrocytes__P5"

SENDER_SHORTLIST <- file.path(R3, "P5_RESTART3_SENDER_LINEAGE_SHORTLIST.tsv")
DONOR_PROGRAM_LONG <- file.path(R3, "P5_RESTART3_FIX1_DONOR_PROGRAM_SCORE_LONG.tsv")
P5_TOP100_FILE <- file.path(R2, "P5_RESTART2_TOP100_WITH_LOADINGS.tsv")

H5ROOT <- file.path(
  ROOT, "04_tables", "32B_CNMf_WORK_20260728_215801", "_h5ad_inputs"
)

H5 <- c(
  Astrocytes = file.path(H5ROOT, "Astrocytes_primary.h5ad"),
  Microglia_macrophage = file.path(H5ROOT, "Microglia_macrophage_primary.h5ad"),
  Oligodendrocytes = file.path(H5ROOT, "Oligodendrocytes_primary.h5ad")
)

req <- c(SENDER_SHORTLIST, DONOR_PROGRAM_LONG, P5_TOP100_FILE, H5)
miss <- req[!file.exists(req)]
if (length(miss)) stop("Missing frozen input(s):\n", paste(miss, collapse = "\n"))

MIN_CELLS <- 20L
MIN_CPM_AVAILABLE <- 1
MIN_DONORS_AVAILABLE <- 6L
MIN_CORR_ABS <- 0.30
MIN_LODO <- 0.75
NICHENET_TOP_TARGETS <- 250L
NICHENET_FDR_SUPPORT <- 0.10
MIN_P5_TARGET_OVERLAP <- 2L

normalize_condition <- function(z) {
  y <- tolower(trimws(as.character(z)))
  out <- rep(NA_character_, length(y))
  out[grepl("tle|epilep|case", y)] <- "TLE"
  out[grepl("control|normal|ctrl|non.?epilep", y)] <- "Control"
  out
}

# -------------------------------------------------------------------------
# 1. Frozen sender scope + P5 donor vector
# -------------------------------------------------------------------------

sender_shortlist <- fread(SENDER_SHORTLIST)
sender_shortlist[, next_LR_scope := as.logical(next_LR_scope)]
senders <- sender_shortlist[next_LR_scope == TRUE, partner_lineage]

if (!setequal(senders, c("Astrocytes", "Microglia_macrophage"))) {
  stop("Unexpected frozen sender scope: ", paste(senders, collapse = ", "))
}

program_long <- fread(DONOR_PROGRAM_LONG)
if (!all(c("donor","program","score","condition") %in% names(program_long))) {
  stop("RESTART3 donor-program table schema changed.")
}

p5_score <- program_long[
  program == P5,
  .(
    donor = trimws(as.character(donor)),
    P5_score = as.numeric(score),
    condition = trimws(as.character(condition))
  )
]

if (nrow(p5_score) != 8L || uniqueN(p5_score$donor) != 8L) {
  stop("Expected 8 unique P5 donors; found ", nrow(p5_score))
}

p5_score[, condition_norm := normalize_condition(condition)]

if (
  p5_score[condition_norm == "TLE", .N] != 4L ||
  p5_score[condition_norm == "Control", .N] != 4L
) {
  stop("Expected 4 TLE + 4 Control in frozen P5 scores.")
}

p5_top <- fread(P5_TOP100_FILE)
if (!"gene" %in% names(p5_top)) stop("P5 Top100 file lacks gene column.")

P5_TOP100 <- unique(toupper(trimws(as.character(p5_top$gene))))
P5_TOP100 <- P5_TOP100[nzchar(P5_TOP100)]
if (length(P5_TOP100) != 100L) stop("Frozen P5 Top100 != 100 genes.")

# -------------------------------------------------------------------------
# 2. Use established Python environment to export donor pseudobulk counts
message("[RESTART4] 1/6 Exporting frozen H5AD donor pseudobulk counts...")
# -------------------------------------------------------------------------

py_candidates <- unique(c(
  file.path(ROOT, "04_tables", "python_envs", "scvi312", "python.exe"),
  file.path(ROOT, "04_tables", "python_envs", "cnmf310", "python.exe"),
  Sys.which("python"),
  Sys.which("python3")
))
py_candidates <- py_candidates[nzchar(py_candidates) & file.exists(py_candidates)]

if (!length(py_candidates)) stop("No Python executable found.")

py_ok <- vapply(
  py_candidates,
  function(py) {
    z <- tryCatch(
      suppressWarnings(
        system2(
        py,
        c("-c", shQuote("import anndata,numpy,pandas,scipy; print('PY_OK')")),
        stdout = TRUE,
        stderr = TRUE
        )
      ),
      error = function(e) character()
    )
    any(grepl("PY_OK", z, fixed = TRUE))
  },
  logical(1)
)

if (!any(py_ok)) {
  stop("No established Python can import anndata/numpy/pandas/scipy.")
}
PY <- py_candidates[which(py_ok)[1]]

EXPORT <- file.path(OUT, "00_H5AD_PSEUDOBULK_EXPORT")
dir.create(EXPORT, recursive = TRUE, showWarnings = FALSE)

PYFILE <- file.path(EXPORT, "export_pseudobulk.py")

py_code <- c(
  "from pathlib import Path",
  "import sys",
  "import anndata as ad",
  "import numpy as np",
  "import pandas as pd",
  "from scipy import sparse",
  "",
  "astro, micro, oligo, outdir = sys.argv[1:5]",
  "files = {",
  "    \"Astrocytes\": Path(astro),",
  "    \"Microglia_macrophage\": Path(micro),",
  "    \"Oligodendrocytes\": Path(oligo),",
  "}",
  "out = Path(outdir)",
  "out.mkdir(parents=True, exist_ok=True)",
  "",
  "def vals(mat):",
  "    if sparse.issparse(mat):",
  "        x = np.asarray(mat.data, dtype=float)",
  "    else:",
  "        a = np.asarray(mat)",
  "        step = max(1, a.size // 200000)",
  "        x = np.asarray(a.ravel()[::step], dtype=float)",
  "    return x[np.isfinite(x)]",
  "",
  "def count_like(mat):",
  "    x = vals(mat)",
  "    if x.size == 0:",
  "        return True",
  "    return (x.min() >= -1e-8) and (np.max(np.abs(x - np.rint(x))) < 1e-6)",
  "",
  "def choose_counts(a):",
  "    for nm in (\"counts\", \"raw_counts\", \"count\", \"raw_count\"):",
  "        if nm in a.layers and count_like(a.layers[nm]):",
  "            return a.layers[nm], \"layer:\" + nm",
  "    if count_like(a.X):",
  "        return a.X, \"X\"",
  "    raise RuntimeError(\"No count-like matrix found.\")",
  "",
  "audit = []",
  "",
  "for lineage, path in files.items():",
  "    a = ad.read_h5ad(path)",
  "    if \"donor_id\" not in a.obs.columns:",
  "        raise RuntimeError(f\"{lineage}: missing obs column donor_id\")",
  "",
  "    donors = a.obs[\"donor_id\"].astype(str)",
  "    donor_values = sorted(donors.unique().tolist())",
  "    if len(donor_values) != 8:",
  "        raise RuntimeError(f\"{lineage}: expected 8 donors, got {len(donor_values)}\")",
  "",
  "    mat, source = choose_counts(a)",
  "    mat = sparse.csr_matrix(mat)",
  "",
  "    genes = np.asarray(a.var_names.astype(str))",
  "    if pd.Index(genes).duplicated().any():",
  "        raise RuntimeError(f\"{lineage}: duplicated feature names\")",
  "",
  "    sums = []",
  "",
  "    for d in donor_values:",
  "        mask = np.asarray(donors == d)",
  "        v = np.asarray(mat[mask, :].sum(axis=0)).ravel()",
  "        sums.append(v)",
  "        audit.append({",
  "            \"lineage\": lineage,",
  "            \"donor\": d,",
  "            \"n_cells\": int(mask.sum()),",
  "            \"n_genes\": int(a.n_vars),",
  "            \"count_source\": source,",
  "            \"total_counts\": float(v.sum()),",
  "        })",
  "",
  "    pb = np.vstack(sums).T",
  "    df = pd.DataFrame(pb, columns=donor_values)",
  "    df.insert(0, \"gene\", genes)",
  "    df.to_csv(out / f\"{lineage}__DONOR_PSEUDOBULK_COUNTS.tsv\", sep=\"\\t\", index=False)",
  "",
  "pd.DataFrame(audit).to_csv(out / \"H5AD_DONOR_LINEAGE_AUDIT.tsv\", sep=\"\\t\", index=False)",
  "print(\"EXPORT_OK\")"
)

writeLines(py_code, PYFILE)

pylog <- system2(
  PY,
  c(
    shQuote(PYFILE),
    shQuote(H5[["Astrocytes"]]),
    shQuote(H5[["Microglia_macrophage"]]),
    shQuote(H5[["Oligodendrocytes"]]),
    shQuote(EXPORT)
  ),
  stdout = TRUE,
  stderr = TRUE
)

writeLines(pylog, file.path(OUT, "P5_RESTART4_FIX8_PYTHON_EXPORT_LOG.txt"))

if (!any(grepl("EXPORT_OK", pylog, fixed = TRUE))) {
  stop(
    "H5AD export failed:\n",
    paste(tail(pylog, 30), collapse = "\n"),
    "\nSee P5_RESTART4_FIX8_PYTHON_EXPORT_LOG.txt"
  )
}

audit <- fread(file.path(EXPORT, "H5AD_DONOR_LINEAGE_AUDIT.tsv"))

fwrite(
  audit,
  file.path(OUT, "P5_RESTART4_FIX8_DONOR_LINEAGE_AUDIT.tsv"),
  sep = "\t"
)

# -------------------------------------------------------------------------
# 3. Exact donor/condition alignment
message("[RESTART4] 2/6 Checking donor/condition alignment...")
# -------------------------------------------------------------------------

h5don <- sort(unique(audit$donor))
p5don <- sort(unique(p5_score$donor))

fwrite(
  data.table(
    source = c("PRIMARY_H5AD", "FROZEN_P5"),
    donors = c(paste(h5don, collapse = ";"), paste(p5don, collapse = ";"))
  ),
  file.path(OUT, "P5_RESTART4_FIX8_DONOR_ALIGNMENT.tsv"),
  sep = "\t"
)

if (!setequal(h5don, p5don)) {
  stop(
    "Exact H5AD donor_id still does not match P5 donor IDs.\nH5AD: ",
    paste(h5don, collapse = ", "),
    "\nP5: ", paste(p5don, collapse = ", ")
  )
}

# The original primary cNMF H5AD inputs predate the later analysis_condition
# field. After exact donor_id alignment, use the frozen P5/34A2 condition key.
condition_key <- unique(
  p5_score[, .(donor, condition_norm)]
)

if (nrow(condition_key) != 8L || anyDuplicated(condition_key$donor)) {
  stop("Frozen P5 condition key is not 8 unique donor-condition rows.")
}

audit <- merge(
  audit,
  condition_key,
  by = "donor",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(audit$condition_norm)) {
  stop("Some exact H5AD donor_id values could not be assigned frozen P5 conditions.")
}

fwrite(
  condition_key,
  file.path(OUT, "P5_RESTART4_FIX8_CONDITION_KEY_FROM_P5.tsv"),
  sep = "\t"
)

coverage <- audit[
  lineage %in% c(senders, "Oligodendrocytes"),
  .(
    donor,
    condition_norm,
    lineage,
    n_cells,
    total_counts,
    eligible = n_cells >= MIN_CELLS
  )
]

fwrite(
  coverage,
  file.path(OUT, "P5_RESTART4_FIX8_DONOR_LINEAGE_COVERAGE.tsv"),
  sep = "\t"
)

cw <- dcast(
  coverage,
  donor ~ lineage,
  value.var = "eligible",
  fill = FALSE
)

expected_lineages <- c(
  senders,
  "Oligodendrocytes"
)

lin_cols <- intersect(
  expected_lineages,
  names(cw)
)

if (!setequal(lin_cols, expected_lineages)) {

  fwrite(
    data.table(
      expected_lineage = expected_lineages,
      present_in_dcast = expected_lineages %in% names(cw)
    ),
    file.path(
      OUT,
      "P5_RESTART4_FIX8_ELIGIBILITY_COLUMN_AUDIT.tsv"
    ),
    sep = "\t"
  )

  stop(
    "Eligibility dcast is missing one or more required lineage columns.\n",
    "Expected: ", paste(expected_lineages, collapse = ", "), "\n",
    "Observed columns: ", paste(names(cw), collapse = ", "), "\n",
    "See P5_RESTART4_FIX8_ELIGIBILITY_COLUMN_AUDIT.tsv"
  )
}

eligibility_matrix <- as.matrix(
  cw[, ..lin_cols]
)

storage.mode(
  eligibility_matrix
) <- "logical"

joint_eligible <- rowSums(
  eligibility_matrix,
  na.rm = FALSE
) == length(lin_cols)

common_donors <- cw[
  joint_eligible,
  donor
]

if (length(common_donors) < 6L) {
  stop("Fewer than 6 donors jointly eligible across sender+receiver lineages.")
}

p5_use <- p5_score[donor %in% common_donors]

# -------------------------------------------------------------------------
# 4. Read exported pseudobulk counts -> TMM/logCPM
message("[RESTART4] 3/6 Building donor-level TMM/logCPM matrices...")
# -------------------------------------------------------------------------

read_pb <- function(lin) {
  f <- file.path(EXPORT, paste0(lin, "__DONOR_PSEUDOBULK_COUNTS.tsv"))
  z <- fread(f)
  genes <- toupper(trimws(as.character(z$gene)))
  if (anyDuplicated(genes)) stop("Duplicated genes in ", lin)

  if (!all(common_donors %in% names(z))) stop("Missing donor columns in ", lin)
  mat <- as.matrix(z[, ..common_donors])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes

  if (
    any(!is.finite(mat)) ||
    any(mat < 0) ||
    max(abs(mat - round(mat))) > 1e-6
  ) {
    stop("Non-count values in ", lin, " pseudobulk.")
  }

  y <- DGEList(counts = round(mat))
  y <- calcNormFactors(y, method = "TMM")

  list(
    cpm = cpm(y, log = FALSE),
    logcpm = cpm(y, log = TRUE, prior.count = 1)
  )
}

target_lineages <- c(senders, "Oligodendrocytes")
expr <- lapply(target_lineages, read_pb)
names(expr) <- target_lineages

# -------------------------------------------------------------------------
# 5. NicheNet resources
message("[RESTART4] 4/6 Loading NicheNet resources...")
# -------------------------------------------------------------------------

RESOURCE_DIR <- file.path(ROOT, "00_REFERENCE_CACHE", "NICHENET_V2")
dir.create(RESOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

find_resource <- function(pattern, min_size) {
  h <- list.files(ROOT, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  h <- h[file.exists(h) & is.finite(file.info(h)$size) & file.info(h)$size >= min_size]
  if (!length(h)) return(NA_character_)
  sc <- 50 * grepl("00_REFERENCE_CACHE", h, ignore.case = TRUE) - nchar(h)/10000
  h[which.max(sc)]
}

LR_FILE <- find_resource("^lr_network_human_21122021\\.rds$", 10000)
LTM_FILE <- find_resource("^ligand_target_matrix(_nsga2r_final)?\\.rds$", 1000000)

if (is.na(LR_FILE)) {
  LR_FILE <- file.path(RESOURCE_DIR, "lr_network_human_21122021.rds")
  url <- "https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"
  tmp <- paste0(LR_FILE, ".part")
  if (file.exists(tmp)) unlink(tmp)
  tryCatch(download.file(url, tmp, mode = "wb", quiet = FALSE), error = function(e) NULL)
  if (!file.exists(tmp) || file.info(tmp)$size < 10000) {
    stop("NicheNet LR network unavailable; retry with network access.")
  }
  file.rename(tmp, LR_FILE)
}

if (is.na(LTM_FILE)) {
  LTM_FILE <- file.path(RESOURCE_DIR, "ligand_target_matrix_nsga2r_final.rds")
  url <- "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
  tmp <- paste0(LTM_FILE, ".part")
  if (file.exists(tmp)) unlink(tmp)
  tryCatch(download.file(url, tmp, mode = "wb", quiet = FALSE), error = function(e) NULL)
  if (!file.exists(tmp) || file.info(tmp)$size < 1000000) {
    stop("NicheNet ligand-target matrix unavailable; retry with network access.")
  }
  file.rename(tmp, LTM_FILE)
}

lr_raw <- readRDS(LR_FILE)
if (all(c("from","to") %in% names(lr_raw))) {
  lr <- data.table(
    ligand = toupper(as.character(lr_raw$from)),
    receptor = toupper(as.character(lr_raw$to))
  )
} else if (all(c("ligand","receptor") %in% names(lr_raw))) {
  lr <- data.table(
    ligand = toupper(as.character(lr_raw$ligand)),
    receptor = toupper(as.character(lr_raw$receptor))
  )
} else stop("NicheNet LR schema not recognized.")

lr <- unique(lr[nzchar(ligand) & nzchar(receptor)])

ltm <- readRDS(LTM_FILE)
if (!is.matrix(ltm)) ltm <- as.matrix(ltm)
rownames(ltm) <- toupper(rownames(ltm))
colnames(ltm) <- toupper(colnames(ltm))

# -------------------------------------------------------------------------
# 6. Donor statistics helpers
# -------------------------------------------------------------------------

cond_assoc <- function(x, y, g) {
  ok <- is.finite(x) & is.finite(y) & !is.na(g)
  x <- x[ok]; y <- y[ok]; g <- g[ok]

  if (length(x) < 6L || uniqueN(g) < 2L) {
    return(list(r = NA_real_, lodo = NA_real_))
  }

  # Condition adjustment by group-centering; algebraically equivalent to
  # residualizing on condition for this two-group design, but much faster.
  rx <- x - ave(x, g, FUN = mean)
  ry <- y - ave(y, g, FUN = mean)

  r <- suppressWarnings(cor(rx, ry))
  s <- sign(r)

  same <- vapply(
    seq_along(x),
    function(i) {
      qx <- x[-i]
      qy <- y[-i]
      qg <- g[-i]

      if (length(qx) < 5L || uniqueN(qg) < 2L) {
        return(NA)
      }

      qrx <- qx - ave(qx, qg, FUN = mean)
      qry <- qy - ave(qy, qg, FUN = mean)

      rr <- suppressWarnings(cor(qrx, qry))

      if (!is.finite(rr) || !is.finite(s)) {
        return(NA)
      }

      if (s == 0) rr == 0 else sign(rr) == s
    },
    logical(1)
  )

  list(
    r = r,
    lodo = mean(same, na.rm = TRUE)
  )
}

group_effect <- function(v, g) {
  ok <- is.finite(v) & g %in% c("TLE", "Control")
  x <- v[ok]
  gg <- g[ok]

  if (
    sum(gg == "TLE") < 1L ||
    sum(gg == "Control") < 1L
  ) {
    return(NA_real_)
  }

  mean(x[gg == "TLE"]) -
    mean(x[gg == "Control"])
}

# -------------------------------------------------------------------------
# 7. Gene metrics
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 7. Targeted gene metrics only
# -------------------------------------------------------------------------

message("[RESTART4] 5/6 Computing targeted ligand/receptor donor metrics...")

donor_order <- sort(common_donors)
po <- p5_use[match(donor_order, donor)]

if (anyNA(po$P5_score)) {
  stop("P5 donor alignment failed.")
}

p5_vec <- po$P5_score
cond <- po$condition_norm

# Receiver background for NicheNet is calculated vectorially from all
# oligodendrocyte genes; expensive donor correlations are NOT computed genome-wide.
receiver_cpm_all <- expr[["Oligodendrocytes"]]$cpm[
  ,
  donor_order,
  drop = FALSE
]

receiver_available_n <- rowSums(
  receiver_cpm_all >= MIN_CPM_AVAILABLE,
  na.rm = TRUE
)

receiver_bg <- intersect(
  rownames(receiver_cpm_all)[
    receiver_available_n >= MIN_DONORS_AVAILABLE
  ],
  rownames(ltm)
)

p5_eval <- intersect(
  P5_TOP100,
  receiver_bg
)

fwrite(
  data.table(
    frozen_P5_top100 = 100L,
    evaluable_P5_targets = length(p5_eval),
    receiver_background = length(receiver_bg)
  ),
  file.path(
    OUT,
    "P5_RESTART4_FIX8_P5_TARGET_COVERAGE.tsv"
  ),
  sep = "\t"
)

if (length(p5_eval) < 40L) {
  stop(
    "Too few P5 Top100 genes evaluable: ",
    length(p5_eval),
    "/100"
  )
}

# Only genes that can participate in the frozen LR question receive the
# expensive donor-correlation/LODO calculation.
sender_candidate_genes <- lapply(
  senders,
  function(lin) {
    intersect(
      rownames(expr[[lin]]$logcpm),
      intersect(
        lr$ligand,
        colnames(ltm)
      )
    )
  }
)
names(sender_candidate_genes) <- senders

receiver_candidate_genes <- intersect(
  rownames(expr[["Oligodendrocytes"]]$logcpm),
  lr$receptor
)

calc_targeted_metrics <- function(lin, genes) {

  genes <- unique(genes)

  lc <- expr[[lin]]$logcpm[
    genes,
    donor_order,
    drop = FALSE
  ]

  rc <- expr[[lin]]$cpm[
    genes,
    donor_order,
    drop = FALSE
  ]

  rows <- vector("list", length(genes))

  for (i in seq_along(genes)) {

    vals <- as.numeric(lc[i, ])
    cpms <- as.numeric(rc[i, ])

    a <- cond_assoc(
      vals,
      p5_vec,
      cond
    )

    rows[[i]] <- data.table(
      lineage = lin,
      gene = genes[i],
      available_donors = sum(
        cpms >= MIN_CPM_AVAILABLE,
        na.rm = TRUE
      ),
      median_CPM = median(
        cpms,
        na.rm = TRUE
      ),
      P5_adjusted_r = a$r,
      P5_LODO = a$lodo,
      TLE_minus_Control_logCPM = group_effect(
        vals,
        cond
      )
    )
  }

  rbindlist(rows)
}

metrics <- list()

for (lin in senders) {

  message(
    "[RESTART4]   sender ",
    lin,
    ": ",
    length(sender_candidate_genes[[lin]]),
    " LR ligands"
  )

  metrics[[lin]] <- calc_targeted_metrics(
    lin,
    sender_candidate_genes[[lin]]
  )
}

message(
  "[RESTART4]   receiver Oligodendrocytes: ",
  length(receiver_candidate_genes),
  " LR receptors"
)

metrics[["Oligodendrocytes"]] <- calc_targeted_metrics(
  "Oligodendrocytes",
  receiver_candidate_genes
)

receiver <- metrics[["Oligodendrocytes"]]

# -------------------------------------------------------------------------
# 8. NicheNet compatibility with frozen P5 Top100
# -------------------------------------------------------------------------

sender_ligands <- unique(unlist(lapply(senders, function(lin) {
  metrics[[lin]][available_donors >= MIN_DONORS_AVAILABLE, gene]
})))

candidate_ligands <- intersect(sender_ligands, intersect(lr$ligand, colnames(ltm)))

lt_rows <- lapply(candidate_ligands, function(lig) {
  w <- ltm[receiver_bg, lig]
  w[!is.finite(w)] <- 0
  ord <- order(w, decreasing = TRUE)
  top <- receiver_bg[ord[seq_len(min(NICHENET_TOP_TARGETS, length(ord)))]]
  ov <- intersect(p5_eval, top)

  data.table(
    ligand = lig,
    receiver_background_N = length(receiver_bg),
    P5_targets_K = length(p5_eval),
    top_targets_m = length(top),
    P5_overlap_k = length(ov),
    P5_overlap_genes = paste(ov, collapse = ";"),
    hypergeom_p = phyper(
      length(ov) - 1,
      length(p5_eval),
      length(receiver_bg) - length(p5_eval),
      length(top),
      lower.tail = FALSE
    ),
    P5_target_weight_sum = sum(ltm[p5_eval, lig], na.rm = TRUE)
  )
})

ligand_target <- rbindlist(lt_rows)
ligand_target[, NicheNet_P5_FDR := p.adjust(hypergeom_p, method = "BH")]

fwrite(
  ligand_target,
  file.path(OUT, "P5_RESTART4_NICHENET_P5_TARGET_COMPATIBILITY.tsv"),
  sep = "\t"
)

# -------------------------------------------------------------------------
# 9. Build targeted LR edges
message("[RESTART4] 6/6 Integrating targeted LR evidence...")
# -------------------------------------------------------------------------

edge_list <- lapply(senders, function(lin) {
  sg <- metrics[[lin]][available_donors >= MIN_DONORS_AVAILABLE, gene]
  rg <- receiver[available_donors >= MIN_DONORS_AVAILABLE, gene]
  e <- lr[ligand %in% sg & receptor %in% rg]
  if (!nrow(e)) return(NULL)
  e[, sender_lineage := lin]
  e
})

edges <- unique(rbindlist(edge_list, fill = TRUE))
if (!nrow(edges)) stop("No targeted LR edges passed availability gates.")

sm <- rbindlist(lapply(senders, function(lin) {
  z <- copy(metrics[[lin]])
  z[, sender_lineage := lin]
  z
}))

setnames(
  sm,
  c(
    "gene",
    "available_donors",
    "median_CPM",
    "P5_adjusted_r",
    "P5_LODO",
    "TLE_minus_Control_logCPM"
  ),
  c(
    "ligand",
    "ligand_available_donors",
    "ligand_median_CPM",
    "ligand_P5_adjusted_r",
    "ligand_P5_LODO",
    "ligand_TLE_minus_Control_logCPM"
  )
)
sm[, lineage := NULL]

rmx <- copy(receiver)
setnames(
  rmx,
  c(
    "gene",
    "available_donors",
    "median_CPM",
    "P5_adjusted_r",
    "P5_LODO",
    "TLE_minus_Control_logCPM"
  ),
  c(
    "receptor",
    "receptor_available_donors",
    "receptor_median_CPM",
    "receptor_P5_adjusted_r",
    "receptor_P5_LODO",
    "receptor_TLE_minus_Control_logCPM"
  )
)
rmx[, lineage := NULL]

edges <- merge(edges, sm, by = c("sender_lineage","ligand"), all.x = TRUE, sort = FALSE)
edges <- merge(edges, rmx, by = "receptor", all.x = TRUE, sort = FALSE)
edges <- merge(edges, ligand_target, by = "ligand", all.x = TRUE, sort = FALSE)

edges[, positive_ligand_P5_coupling :=
  is.finite(ligand_P5_adjusted_r) &
  ligand_P5_adjusted_r >= MIN_CORR_ABS &
  is.finite(ligand_P5_LODO) &
  ligand_P5_LODO >= MIN_LODO
]

edges[, inverse_ligand_P5_coupling :=
  is.finite(ligand_P5_adjusted_r) &
  ligand_P5_adjusted_r <= -MIN_CORR_ABS &
  is.finite(ligand_P5_LODO) &
  ligand_P5_LODO >= MIN_LODO
]

edges[, receptor_P5_context_support :=
  is.finite(receptor_P5_adjusted_r) &
  receptor_P5_adjusted_r >= MIN_CORR_ABS &
  is.finite(receptor_P5_LODO) &
  receptor_P5_LODO >= MIN_LODO
]

edges[, NicheNet_P5_support :=
  is.finite(NicheNet_P5_FDR) &
  NicheNet_P5_FDR < NICHENET_FDR_SUPPORT &
  P5_overlap_k >= MIN_P5_TARGET_OVERLAP
]

edges[, ligand_TLE_direction_support :=
  is.finite(ligand_TLE_minus_Control_logCPM) &
  ligand_TLE_minus_Control_logCPM > 0
]

edges[, positive_support_count :=
  as.integer(positive_ligand_P5_coupling) +
  as.integer(receptor_P5_context_support) +
  as.integer(NicheNet_P5_support) +
  as.integer(ligand_TLE_direction_support)
]

edges[, evidence_tier := fcase(
  positive_ligand_P5_coupling & NicheNet_P5_support & ligand_TLE_direction_support,
  "TIER_A_POSITIVE_P5_COMPATIBLE",

  positive_ligand_P5_coupling & positive_support_count >= 2L,
  "TIER_B_POSITIVE_P5_COMPATIBLE",

  inverse_ligand_P5_coupling,
  "INVERSE_OR_PUTATIVE_SUPPRESSIVE_CONTEXT",

  NicheNet_P5_support | receptor_P5_context_support | ligand_TLE_direction_support,
  "TIER_C_EXPLORATORY",

  default = "LOW_PRIORITY"
)]

tier_ord <- c(
  TIER_A_POSITIVE_P5_COMPATIBLE = 1L,
  TIER_B_POSITIVE_P5_COMPATIBLE = 2L,
  TIER_C_EXPLORATORY = 3L,
  INVERSE_OR_PUTATIVE_SUPPRESSIVE_CONTEXT = 4L,
  LOW_PRIORITY = 5L
)

edges[, tier_rank := unname(tier_ord[evidence_tier])]
setorder(
  edges,
  tier_rank,
  -positive_support_count,
  NicheNet_P5_FDR,
  -ligand_P5_adjusted_r,
  -ligand_TLE_minus_Control_logCPM,
  sender_lineage,
  ligand,
  receptor
)
edges[, final_rank := seq_len(.N)]

fwrite(
  edges,
  file.path(OUT, "P5_RESTART4_ALL_TARGETED_LR_EDGES.tsv"),
  sep = "\t"
)

shortlist <- edges[
  evidence_tier %in% c(
    "TIER_A_POSITIVE_P5_COMPATIBLE",
    "TIER_B_POSITIVE_P5_COMPATIBLE"
  )
]

if (nrow(shortlist) < 10L) {
  fill <- edges[evidence_tier == "TIER_C_EXPLORATORY"][
    seq_len(min(.N, 10L - nrow(shortlist)))
  ]
  shortlist <- rbindlist(list(shortlist, fill), fill = TRUE)
}

shortlist <- unique(shortlist, by = c("sender_lineage","ligand","receptor"))
shortlist <- head(shortlist, 20L)

fwrite(
  shortlist,
  file.path(OUT, "P5_RESTART4_LR_SHORTLIST.tsv"),
  sep = "\t"
)

# -------------------------------------------------------------------------
# 10. Decision/report
# -------------------------------------------------------------------------

nA <- edges[evidence_tier == "TIER_A_POSITIVE_P5_COMPATIBLE", .N]
nB <- edges[evidence_tier == "TIER_B_POSITIVE_P5_COMPATIBLE", .N]
nInv <- edges[evidence_tier == "INVERSE_OR_PUTATIVE_SUPPRESSIVE_CONTEXT", .N]

FINAL_CLASS <- if (nA >= 1L) {
  "AT_LEAST_ONE_PRIORITY_POSITIVE_LR_AXIS_FOR_SPATIAL_QUALIFICATION"
} else if (nB >= 1L) {
  "SUPPORTIVE_POSITIVE_LR_AXES_ONLY_NO_PRIMARY_AXIS_YET"
} else {
  "NO_POSITIVE_LR_AXIS_PASSES_TARGETED_INTEGRATION_GATE"
}

NEXT_STAGE <- if (nA >= 1L || nB >= 1L) {
  "RESTART5_SPATIAL_OR_REFERENCE_QUALIFICATION_OF_FROZEN_LR_SHORTLIST"
} else {
  "RESTART5_P5_SPATIAL_REFERENCE_MAPPING_WITHOUT_MECHANISTIC_LR_CLAIM"
}

decision <- data.table(
  item = c(
    "input_strategy",
    "frozen_sender_lineages",
    "receiver_lineage",
    "n_common_donors",
    "n_LR_edges",
    "n_Tier_A",
    "n_Tier_B",
    "n_inverse_context",
    "FINAL_CLASS",
    "NEXT_STAGE"
  ),
  value = c(
    "EXACT_PRIMARY_CNMf_H5AD_INPUTS_DONOR_ID_ALIGNED",
    paste(senders, collapse = ";"),
    "Oligodendrocytes",
    length(common_donors),
    nrow(edges),
    nA,
    nB,
    nInv,
    FINAL_CLASS,
    NEXT_STAGE
  )
)

fwrite(decision, file.path(OUT, "P5_RESTART4_DECISION.tsv"), sep = "\t")

top_txt <- if (nrow(shortlist)) {
  paste(
    paste0(
      shortlist$sender_lineage, " -> ",
      shortlist$ligand, " -> ",
      shortlist$receptor, " | ",
      shortlist$evidence_tier,
      " | ligand-P5 r=", signif(shortlist$ligand_P5_adjusted_r, 4),
      " | LODO=", signif(shortlist$ligand_P5_LODO, 3),
      " | P5 overlap=", shortlist$P5_overlap_k,
      " | NicheNet FDR=", signif(shortlist$NicheNet_P5_FDR, 4),
      " | ligand TLE-Control=", signif(shortlist$ligand_TLE_minus_Control_logCPM, 4)
    ),
    collapse = "\n"
  )
} else "No LR candidate entered the shortlist."

writeLines(
  c(
    "P5 MATURE OLIGODENDROCYTE PURE-BIOINFORMATICS MANUSCRIPT",
    "RESTART-4 FIX8 FAST | TARGETED DONOR-AWARE LR",
    "=========================================",
    "",
    "INPUT CORRECTION",
    "----------------",
    "- The blinded CP1R B01-B08 object is not used for donor matching.",
    "- Exact primary cNMF H5AD inputs are used because their donor_id generated the P5/34A2 donor-score evidence. The original H5ADs need only donor_id; condition is attached from the frozen P5/34A2 donor-score table after exact donor alignment.",
    "- No donor identity is guessed from condition alone.",
    "",
    paste0("Frozen sender lineages: ", paste(senders, collapse = ", ")),
    "Receiver: Oligodendrocytes",
    paste0("Common eligible donors: ", length(common_donors)),
    paste0("Targeted LR edges: ", nrow(edges)),
    paste0("Tier A: ", nA),
    paste0("Tier B: ", nB),
    paste0("Inverse context: ", nInv),
    "",
    paste0("FINAL_CLASS=", FINAL_CLASS),
    paste0("NEXT_STAGE=", NEXT_STAGE),
    "",
    "TOP SHORTLIST",
    "-------------",
    top_txt,
    "",
    "INTERPRETATION BOUNDARY",
    "-----------------------",
    "- Donor is the biological replicate.",
    "- FIX8 computes expensive donor association only for LR-relevant ligand/receptor genes; whole-genome exact permutation tests were removed because they were not used by the LR decision rule.",
    "- Sender scope is frozen from RESTART3 FIX1.",
    "- NicheNet is ligand-target compatibility with the frozen P5 Top100, not a causal p-value.",
    "- P5 remains directional/high-priority exploratory rather than formally FDR-supported.",
    "- Tier A/B means computational mechanism candidate only.",
    "- No secretion, receptor activation, physical binding, or causal claim is authorized.",
    "- No TF is forced into the mechanism.",
    "- Exact primary cNMF lineage inputs are used; this is not whole-tissue communication inference."
  ),
  file.path(OUT, "P5_RESTART4_REPORT.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(OUT, "P5_RESTART4_SESSIONINFO.txt")
)

message("P5 RESTART-4 FIX8 FAST completed.")
message("Output: ", OUT)

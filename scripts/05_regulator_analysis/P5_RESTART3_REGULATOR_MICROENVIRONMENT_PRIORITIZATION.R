# P5_RESTART3_REGULATOR_MICROENVIRONMENT_PRIORITIZATION.R
# Pure-bioinformatics manuscript: Oligodendrocytes__P5
#
# Purpose:
#   Reuse already-completed donor-aware mechanism outputs to build ONE compact,
#   publication-oriented mechanism shortlist for P5.
#
# This script does NOT rerun TF inference, cNMF, QC, clustering, DE, external
# validation, LR inference, spatial analysis, or any new machine-learning model.
#
# It integrates:
#   A) existing 34A3 intrinsic TF/regulator evidence for Oligodendrocytes__P5
#   B) existing 34A2 donor-level cross-program correlations to identify
#      cross-lineage microenvironment programs that covary with P5
#
# Output:
#   - ranked P5 regulator candidates from existing 34A3 evidence
#   - ranked cross-lineage program partners from existing 34A2 donor evidence
#   - a small sender-lineage shortlist for the NEXT targeted LR analysis
#
# Important:
#   Cross-lineage program correlation is NOT ligand-receptor evidence.
#   TF association is NOT causal regulation.
#   P5 remains a high-priority directional/exploratory disease-linked program.

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- Sys.getenv("P5_WORK_ROOT")
RELEASE_ROOT <- file.path(ROOT, "99_frozen_releases")

PROJECT <- file.path(
  ROOT,
  "P5_MATURE_OLIGODENDROCYTE_PURE_BIOINFO_RESTART"
)

R2 <- file.path(
  PROJECT,
  "RESTART2_FUNCTIONAL_ARCHITECTURE"
)

OUT <- file.path(
  PROJECT,
  "RESTART3_REGULATOR_MICROENVIRONMENT_PRIORITIZATION"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

P5 <- "Oligodendrocytes__P5"

R2_DECISION <- file.path(R2, "P5_RESTART2_DECISION.tsv")
if (!file.exists(R2_DECISION)) {
  stop("Missing RESTART-2 decision file: ", R2_DECISION)
}

r2_decision <- fread(R2_DECISION)

get_r2 <- function(key) {
  x <- r2_decision[item == key, value]
  if (length(x) != 1L) return(NA_character_)
  as.character(x)
}

# Preserve frozen interpretation boundary.
functional_class <- get_r2("FUNCTIONAL_CLASS")
senescence_wording <- get_r2("SENESCENCE_WORDING")

# =============================================================================
# Helpers
# =============================================================================

latest_release <- function(prefix) {

  dirs <- list.dirs(
    RELEASE_ROOT,
    full.names = TRUE,
    recursive = FALSE
  )

  dirs <- dirs[
    grepl(paste0("^", prefix), basename(dirs))
  ]

  if (!length(dirs)) return(NA_character_)

  info <- file.info(dirs)
  dirs[order(info$mtime, decreasing = TRUE)][1]
}

find_one <- function(root, basename_target, required = TRUE) {

  if (is.na(root) || !dir.exists(root)) {
    if (required) stop("Missing release root for ", basename_target)
    return(NA_character_)
  }

  hits <- list.files(
    root,
    recursive = TRUE,
    full.names = TRUE
  )

  hits <- hits[basename(hits) == basename_target]

  if (length(hits) == 1L) return(hits)

  if (required) {
    stop(
      "Could not uniquely locate ", basename_target,
      "\nHits:\n", paste(hits, collapse = "\n")
    )
  }

  NA_character_
}

contains_p5 <- function(dt) {

  if (is.null(dt) || !nrow(dt)) return(rep(FALSE, 0))

  hit <- rep(FALSE, nrow(dt))

  for (cc in names(dt)) {
    if (is.character(dt[[cc]]) || is.factor(dt[[cc]])) {
      z <- as.character(dt[[cc]])
      hit <- hit |
        z == P5 |
        grepl("Oligodendrocytes__P5", z, fixed = TRUE)
    }
  }

  hit
}

pick_first_col <- function(nms, patterns) {

  low <- tolower(nms)

  for (p in patterns) {
    hit <- grep(p, low)
    if (length(hit)) return(nms[hit[1]])
  }

  NA_character_
}

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

to_bool <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z %in% c("true", "t", "1", "yes")
}

extract_lineage <- function(program_name) {
  z <- as.character(program_name)
  sub("__.*$", "", z)
}

# =============================================================================
# 1. Recover authoritative 34A2 / 34A3 releases
# =============================================================================

A2 <- latest_release(
  "34A2_DONOR_LEVEL_PROGRAM_MECHANISM_INPUT_AUDIT_"
)

A3 <- latest_release(
  "34A3_DONOR_LEVEL_INTRINSIC_TF_ACTIVITY_"
)

if (is.na(A2)) {
  stop("No 34A2 donor-level mechanism-input release found.")
}

if (is.na(A3)) {
  stop("No 34A3 intrinsic-TF release found.")
}

COR_FILE <- find_one(
  A2,
  "34A2_DONOR_PROGRAM_SCORE_CORRELATIONS.csv"
)

SCORE_FILE <- find_one(
  A2,
  "34A2_DONOR_PROGRAM_SCORE_TABLE.csv",
  required = FALSE
)

TF_CONC_FILE <- find_one(
  A3,
  "34A3_TF_PROGRAM_CROSS_NETWORK_CONCORDANCE.csv"
)

TF_ALL_FILE <- find_one(
  A3,
  "34A3_TF_PROGRAM_ASSOCIATIONS_ALL_NETWORKS.csv.gz",
  required = FALSE
)

provenance <- data.table(
  layer = c(
    "RESTART2_DECISION",
    "34A2_DONOR_PROGRAM_SCORE_CORRELATIONS",
    "34A2_DONOR_PROGRAM_SCORE_TABLE",
    "34A3_TF_CROSS_NETWORK_CONCORDANCE",
    "34A3_TF_ASSOCIATIONS_ALL_NETWORKS"
  ),
  path = c(
    R2_DECISION,
    COR_FILE,
    SCORE_FILE,
    TF_CONC_FILE,
    TF_ALL_FILE
  ),
  exists = file.exists(c(
    R2_DECISION,
    COR_FILE,
    ifelse(is.na(SCORE_FILE), "", SCORE_FILE),
    TF_CONC_FILE,
    ifelse(is.na(TF_ALL_FILE), "", TF_ALL_FILE)
  ))
)

fwrite(
  provenance,
  file.path(OUT, "P5_RESTART3_INPUT_PROVENANCE.tsv"),
  sep = "\t"
)

# =============================================================================
# 2. Existing 34A3 intrinsic regulator evidence
# =============================================================================

tf_conc <- fread(TF_CONC_FILE)
p5_tf <- tf_conc[contains_p5(tf_conc)]

if (!nrow(p5_tf)) {
  warning("No explicit P5 rows found in 34A3 TF concordance table.")
}

# Recover TF/regulator column.
tf_col <- pick_first_col(
  names(p5_tf),
  c("^tf$", "transcription.*factor", "^regulator$", "^source$", "source_gene")
)

if (nrow(p5_tf) && is.na(tf_col)) {
  stop(
    "P5 rows found in 34A3 TF concordance, but TF column could not be identified.\n",
    "Columns: ", paste(names(p5_tf), collapse = ", ")
  )
}

# Build generic ranking components from available columns.
if (nrow(p5_tf)) {

  p5_tf[, TF := as.character(get(tf_col))]

  # Evidence-support count: sum TRUE-like fields mentioning concordance/support.
  support_cols <- names(p5_tf)[
    grepl(
      "support|concord|same_direction|replic|network",
      tolower(names(p5_tf))
    )
  ]

  logical_support_cols <- support_cols[
    vapply(
      p5_tf[, ..support_cols],
      function(x) {
        is.logical(x) ||
          all(
            na.omit(tolower(as.character(x))) %in%
              c("true", "false", "t", "f", "1", "0", "yes", "no")
          )
      },
      logical(1)
    )
  ]

  if (length(logical_support_cols)) {
    p5_tf[, support_count_recovered :=
            rowSums(
              as.data.frame(
                lapply(.SD, to_bool)
              ),
              na.rm = TRUE
            ),
          .SDcols = logical_support_cols]
  } else {
    p5_tf[, support_count_recovered := 0L]
  }

  # Prefer an existing FDR/q-value if present.
  fdr_col <- pick_first_col(
    names(p5_tf),
    c("fdr", "qvalue", "q_value", "adj.*p", "padj")
  )

  if (!is.na(fdr_col)) {
    p5_tf[, recovered_FDR := to_num(get(fdr_col))]
  } else {
    p5_tf[, recovered_FDR := NA_real_]
  }

  # Prefer an existing correlation/effect statistic if present.
  assoc_col <- pick_first_col(
    names(p5_tf),
    c("rho", "correlation", "spearman", "effect", "estimate", "association")
  )

  if (!is.na(assoc_col)) {
    p5_tf[, recovered_association := to_num(get(assoc_col))]
  } else {
    p5_tf[, recovered_association := NA_real_]
  }

  p5_tf[, abs_association := abs(recovered_association)]

  # Rank without inventing a new significance gate.
  setorder(
    p5_tf,
    -support_count_recovered,
    recovered_FDR,
    -abs_association,
    TF
  )

  p5_tf[, RESTART3_regulator_rank := seq_len(.N)]

  p5_tf[, RESTART3_regulator_class := fifelse(
    support_count_recovered >= 2L &
      (is.na(recovered_FDR) | recovered_FDR < 0.10),
    "PRIORITY_REGULATOR_ASSOCIATED_CANDIDATE",
    fifelse(
      support_count_recovered >= 1L,
      "SUPPORTIVE_REGULATOR_ASSOCIATED_CANDIDATE",
      "EXPLORATORY_REGULATOR_ASSOCIATION"
    )
  )]

  fwrite(
    p5_tf,
    file.path(OUT, "P5_RESTART3_P5_TF_EVIDENCE_RANKED.tsv"),
    sep = "\t"
  )

} else {

  fwrite(
    data.table(
      status = "NO_P5_TF_ROWS_RECOVERED_FROM_EXISTING_34A3"
    ),
    file.path(OUT, "P5_RESTART3_P5_TF_EVIDENCE_RANKED.tsv"),
    sep = "\t"
  )
}

# Recover network-specific detail only as support, no new inference.
p5_tf_all <- data.table()

if (!is.na(TF_ALL_FILE) && file.exists(TF_ALL_FILE)) {

  tf_all <- fread(TF_ALL_FILE)
  p5_tf_all <- tf_all[contains_p5(tf_all)]

  fwrite(
    p5_tf_all,
    file.path(OUT, "P5_RESTART3_P5_TF_NETWORK_DETAILS.tsv"),
    sep = "\t"
  )
}

# =============================================================================
# 3. Existing 34A2 donor-level cross-program correlations
# =============================================================================

cor_dt <- fread(COR_FILE)
p5_cor <- cor_dt[contains_p5(cor_dt)]

if (!nrow(p5_cor)) {
  stop("No P5-related rows found in 34A2 donor program correlation table.")
}

# Identify columns that contain program identifiers.
program_cols <- names(p5_cor)[
  grepl("program", tolower(names(p5_cor)))
]

if (length(program_cols) < 2L) {

  # Fallback: character columns containing "__P" patterns.
  program_cols <- names(p5_cor)[
    vapply(
      p5_cor,
      function(x) {
        if (!(is.character(x) || is.factor(x))) return(FALSE)
        any(grepl("__P[0-9]+", as.character(x)))
      },
      logical(1)
    )
  ]
}

if (length(program_cols) < 2L) {
  stop(
    "Could not identify the two program columns in 34A2 correlations.\n",
    "Columns: ", paste(names(p5_cor), collapse = ", ")
  )
}

# Use first two program-like columns.
p1_col <- program_cols[1]
p2_col <- program_cols[2]

p5_cor[, program_1 := as.character(get(p1_col))]
p5_cor[, program_2 := as.character(get(p2_col))]

p5_cor[, partner_program := fifelse(
  program_1 == P5,
  program_2,
  fifelse(
    program_2 == P5,
    program_1,
    NA_character_
  )
)]

# If one of the program fields embeds P5 with extra text, recover partner safely.
p5_cor[
  is.na(partner_program),
  partner_program := fifelse(
    grepl(P5, program_1, fixed = TRUE),
    program_2,
    fifelse(
      grepl(P5, program_2, fixed = TRUE),
      program_1,
      NA_character_
    )
  )
]

p5_cor <- p5_cor[
  !is.na(partner_program) &
    partner_program != P5
]

# Correlation statistic.
corr_col <- pick_first_col(
  names(p5_cor),
  c("spearman.*rho", "^rho$", "correlation", "pearson", "spearman")
)

if (is.na(corr_col)) {
  stop(
    "Could not identify correlation statistic in 34A2 P5 rows.\n",
    "Columns: ", paste(names(p5_cor), collapse = ", ")
  )
}

p5_cor[, donor_program_correlation := to_num(get(corr_col))]
p5_cor[, abs_donor_program_correlation := abs(donor_program_correlation)]
p5_cor[, partner_lineage := extract_lineage(partner_program)]

# Optional P/FDR fields, if already present.
cor_p_col <- pick_first_col(
  names(p5_cor),
  c("^p$", "pvalue", "p_value", "p\\.value")
)

cor_fdr_col <- pick_first_col(
  names(p5_cor),
  c("fdr", "qvalue", "adj.*p", "padj")
)

p5_cor[, recovered_p := if (!is.na(cor_p_col)) to_num(get(cor_p_col)) else NA_real_]
p5_cor[, recovered_FDR := if (!is.na(cor_fdr_col)) to_num(get(cor_fdr_col)) else NA_real_]

# Cross-lineage context only for microenvironment prioritization.
crosslineage <- p5_cor[
  partner_lineage != "Oligodendrocytes"
]

setorder(
  crosslineage,
  -abs_donor_program_correlation,
  recovered_FDR,
  partner_program
)

crosslineage[, RESTART3_microenvironment_rank := seq_len(.N)]

fwrite(
  crosslineage,
  file.path(OUT, "P5_RESTART3_CROSS_LINEAGE_PROGRAM_CORRELATIONS.tsv"),
  sep = "\t"
)

# Keep only a compact top10 for manuscript-oriented review.
top10_context <- head(crosslineage, 10)

fwrite(
  top10_context,
  file.path(OUT, "P5_RESTART3_TOP10_MICROENVIRONMENT_PROGRAMS.tsv"),
  sep = "\t"
)

# Sender-lineage shortlist: maximum absolute P5 correlation within lineage,
# plus number of top10 programs contributed.
if (nrow(crosslineage)) {

  sender_lineages <- crosslineage[
    ,
    .(
      max_abs_P5_correlation = max(
        abs_donor_program_correlation,
        na.rm = TRUE
      ),
      strongest_partner_program = partner_program[
        which.max(abs_donor_program_correlation)
      ],
      strongest_signed_correlation = donor_program_correlation[
        which.max(abs_donor_program_correlation)
      ],
      n_programs_available = .N
    ),
    by = partner_lineage
  ]

  top10_counts <- top10_context[
    ,
    .(n_top10_programs = .N),
    by = partner_lineage
  ]

  sender_lineages <- merge(
    sender_lineages,
    top10_counts,
    by = "partner_lineage",
    all.x = TRUE
  )

  sender_lineages[is.na(n_top10_programs), n_top10_programs := 0L]

  setorder(
    sender_lineages,
    -n_top10_programs,
    -max_abs_P5_correlation,
    partner_lineage
  )

  sender_lineages[, sender_lineage_rank := seq_len(.N)]

  # Restrict the NEXT LR step to at most 3 sender lineages.
  sender_lineages[, next_LR_scope := sender_lineage_rank <= 3L]

  fwrite(
    sender_lineages,
    file.path(OUT, "P5_RESTART3_SENDER_LINEAGE_SHORTLIST.tsv"),
    sep = "\t"
  )
}

# =============================================================================
# 4. Integrated manuscript-oriented decision
# =============================================================================

n_priority_tf <- 0L
n_supportive_tf <- 0L

if (nrow(p5_tf)) {
  n_priority_tf <- p5_tf[
    RESTART3_regulator_class == "PRIORITY_REGULATOR_ASSOCIATED_CANDIDATE",
    .N
  ]
  n_supportive_tf <- p5_tf[
    RESTART3_regulator_class %in% c(
      "PRIORITY_REGULATOR_ASSOCIATED_CANDIDATE",
      "SUPPORTIVE_REGULATOR_ASSOCIATED_CANDIDATE"
    ),
    .N
  ]
}

n_crosslineage <- nrow(crosslineage)
n_sender_lineages <- if (exists("sender_lineages")) nrow(sender_lineages) else 0L
n_lr_lineages <- if (exists("sender_lineages")) sum(sender_lineages$next_LR_scope) else 0L

regulator_status <- if (n_priority_tf >= 1L) {
  "EXISTING_34A3_PRIORITY_REGULATOR_CANDIDATES_RECOVERED"
} else if (n_supportive_tf >= 1L) {
  "SUPPORTIVE_TF_ASSOCIATIONS_ONLY_DO_NOT_FORCE_MASTER_TF"
} else {
  "NO_CONVERGENT_INTRINSIC_REGULATOR_DO_NOT_ADD_NEW_TF_ALGORITHMS"
}

microenv_status <- if (n_crosslineage >= 1L) {
  "DONOR_LEVEL_CROSS_LINEAGE_CONTEXT_RECOVERED_FOR_TARGETED_LR"
} else {
  "NO_CROSS_LINEAGE_PROGRAM_CONTEXT_RECOVERED"
}

next_stage <- if (n_lr_lineages >= 1L) {
  "RESTART4_TARGETED_DONOR_AWARE_LR_TO_P5_USING_ONLY_SHORTLISTED_SENDER_LINEAGES"
} else {
  "RESTART4_SPATIAL_REFERENCE_MAPPING_NO_LR_EXPANSION"
}

decision <- data.table(
  item = c(
    "P5_FUNCTIONAL_CLASS",
    "P5_SENESCENCE_WORDING",
    "n_P5_TF_rows_recovered",
    "n_priority_TF_candidates",
    "n_supportive_or_priority_TF_candidates",
    "REGULATOR_STATUS",
    "n_cross_lineage_program_correlations",
    "n_sender_lineages_recovered",
    "n_sender_lineages_for_next_LR",
    "MICROENVIRONMENT_STATUS",
    "NEXT_STAGE"
  ),
  value = c(
    functional_class,
    senescence_wording,
    nrow(p5_tf),
    n_priority_tf,
    n_supportive_tf,
    regulator_status,
    n_crosslineage,
    n_sender_lineages,
    n_lr_lineages,
    microenv_status,
    next_stage
  )
)

fwrite(
  decision,
  file.path(OUT, "P5_RESTART3_DECISION.tsv"),
  sep = "\t"
)

# =============================================================================
# 5. Compact report
# =============================================================================

top_tf_text <- if (nrow(p5_tf)) {
  x <- head(p5_tf, 10)
  paste(
    paste0(
      x$TF,
      " [",
      x$RESTART3_regulator_class,
      "; support=",
      x$support_count_recovered,
      ifelse(
        is.finite(x$recovered_association),
        paste0("; assoc=", signif(x$recovered_association, 4)),
        ""
      ),
      ifelse(
        is.finite(x$recovered_FDR),
        paste0("; FDR=", signif(x$recovered_FDR, 4)),
        ""
      ),
      "]"
    ),
    collapse = "\n"
  )
} else {
  "No P5 TF rows recovered from existing 34A3."
}

top_context_text <- if (nrow(top10_context)) {
  paste(
    paste0(
      top10_context$partner_program,
      " [rho/correlation=",
      signif(top10_context$donor_program_correlation, 4),
      "]"
    ),
    collapse = "\n"
  )
} else {
  "No cross-lineage P5 program associations recovered."
}

writeLines(
  c(
    "P5 MATURE OLIGODENDROCYTE PURE-BIOINFORMATICS MANUSCRIPT",
    "RESTART-3 | REGULATOR + MICROENVIRONMENT PRIORITIZATION",
    "=======================================================",
    "",
    paste0("FUNCTIONAL_CLASS=", functional_class),
    paste0("SENESCENCE_WORDING=", senescence_wording),
    "",
    "INTRINSIC REGULATOR BRANCH",
    "--------------------------",
    paste0("P5 TF rows recovered from existing 34A3=", nrow(p5_tf)),
    paste0("priority TF candidates=", n_priority_tf),
    paste0("supportive+priority TF candidates=", n_supportive_tf),
    paste0("REGULATOR_STATUS=", regulator_status),
    "",
    "Top recovered P5 TF/regulator rows:",
    top_tf_text,
    "",
    "MICROENVIRONMENT PROGRAM CONTEXT",
    "--------------------------------",
    paste0("P5 cross-lineage donor-program correlations recovered=", n_crosslineage),
    paste0("sender lineages recovered=", n_sender_lineages),
    paste0("sender lineages retained for next LR=", n_lr_lineages),
    paste0("MICROENVIRONMENT_STATUS=", microenv_status),
    "",
    "Top cross-lineage program partners:",
    top_context_text,
    "",
    paste0("NEXT_STAGE=", next_stage),
    "",
    "INTERPRETATION BOUNDARY",
    "-----------------------",
    "- Existing 34A3 TF evidence is reused; no new TF algorithm is added.",
    "- TF association/concordance identifies regulator-associated candidates, not causal master regulators.",
    "- 34A2 cross-lineage donor correlations identify microenvironment context only; they are not ligand-receptor evidence.",
    "- The next LR analysis, if reached, is restricted to at most three sender lineages selected here.",
    "- P5 remains directional/high-priority exploratory rather than a formally FDR-supported TLE program.",
    "- Senescence remains partial transcriptional overlap only."
  ),
  file.path(OUT, "P5_RESTART3_REPORT.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(OUT, "P5_RESTART3_SESSIONINFO.txt")
)

message("P5 RESTART-3 completed.")
message("Output: ", OUT)

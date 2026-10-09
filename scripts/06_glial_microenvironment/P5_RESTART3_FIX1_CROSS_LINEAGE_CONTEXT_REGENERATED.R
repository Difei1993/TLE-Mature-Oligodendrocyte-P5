# P5_RESTART3_FIX1_CROSS_LINEAGE_CONTEXT.R
# Correction for RESTART-3 microenvironment branch.
#
# Why this fix is needed:
#   The original 34A2_DONOR_PROGRAM_SCORE_CORRELATIONS.csv is a WITHIN-LINEAGE
#   correlation table. For P5 (Oligodendrocytes__P5), all partners therefore
#   belong to Oligodendrocytes and the original RESTART-3 cross-lineage filter
#   necessarily returns zero rows.
#
# This FIX1 does NOT rerun cNMF, QC, clustering, DE, TF inference, LR inference,
# or external validation. It uses the already-existing donor-level program-score
# table to compute cross-lineage donor correlations directly.
#
# Statistical boundary:
#   - donor = biological replicate
#   - correlations are condition-adjusted/descriptive prioritization only
#   - they define the sender-lineage SEARCH SCOPE for a later LR analysis
#   - they are NOT ligand-receptor evidence and NOT disease-effect tests

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

OUT <- file.path(
  PROJECT,
  "RESTART3_REGULATOR_MICROENVIRONMENT_PRIORITIZATION"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

P5 <- "Oligodendrocytes__P5"

# =============================================================================
# 1. Locate authoritative 34A2 donor-level score table
# =============================================================================

dirs <- list.dirs(
  RELEASE_ROOT,
  recursive = FALSE,
  full.names = TRUE
)

dirs <- dirs[
  grepl(
    "^34A2_DONOR_LEVEL_PROGRAM_MECHANISM_INPUT_AUDIT_",
    basename(dirs)
  )
]

if (!length(dirs)) {
  stop("No 34A2 donor-level mechanism-input release found.")
}

info <- file.info(dirs)
A2 <- dirs[order(info$mtime, decreasing = TRUE)][1]

score_hits <- list.files(
  A2,
  recursive = TRUE,
  full.names = TRUE
)
score_hits <- score_hits[
  basename(score_hits) == "34A2_DONOR_PROGRAM_SCORE_TABLE.csv"
]

if (length(score_hits) != 1L) {
  stop(
    "Could not uniquely locate 34A2_DONOR_PROGRAM_SCORE_TABLE.csv.\nHits:\n",
    paste(score_hits, collapse = "\n")
  )
}

SCORE_FILE <- score_hits[1]
x <- fread(SCORE_FILE)

# Save schema immediately for audit.
schema <- data.table(
  column = names(x),
  class = vapply(x, function(z) paste(class(z), collapse = "/"), character(1)),
  n_unique = vapply(x, uniqueN, integer(1)),
  example = vapply(
    x,
    function(z) {
      u <- unique(as.character(z))
      u <- u[!is.na(u)]
      paste(head(u, 5), collapse = " | ")
    },
    character(1)
  )
)

fwrite(
  schema,
  file.path(OUT, "P5_RESTART3_FIX1_SCORE_TABLE_SCHEMA.tsv"),
  sep = "\t"
)

# =============================================================================
# 2. Robustly identify long-table columns
# =============================================================================

pick_exact <- function(nms, candidates) {
  low <- tolower(nms)
  for (cand in tolower(candidates)) {
    hit <- which(low == cand)
    if (length(hit)) return(nms[hit[1]])
  }
  NA_character_
}

pick_pattern <- function(nms, patterns) {
  low <- tolower(nms)
  for (pat in patterns) {
    hit <- grep(pat, low)
    if (length(hit)) return(nms[hit[1]])
  }
  NA_character_
}

DONOR_COL <- pick_exact(
  names(x),
  c("donor_id", "donor", "subject_id", "subject")
)
if (is.na(DONOR_COL)) {
  DONOR_COL <- pick_pattern(names(x), c("donor", "subject"))
}

PROGRAM_COL <- pick_exact(
  names(x),
  c("program_id", "program", "score_name", "program_name")
)
if (is.na(PROGRAM_COL)) {
  PROGRAM_COL <- pick_pattern(names(x), c("program"))
}

LINEAGE_COL <- pick_exact(
  names(x),
  c("lineage", "cell_type", "celltype", "broad_cell_type")
)
if (is.na(LINEAGE_COL)) {
  LINEAGE_COL <- pick_pattern(names(x), c("lineage", "cell.*type"))
}

CONDITION_COL <- pick_exact(
  names(x),
  c("condition", "disease", "group", "status")
)
if (is.na(CONDITION_COL)) {
  CONDITION_COL <- pick_pattern(names(x), c("condition", "disease", "group"))
}

# score column: prioritize known names, otherwise choose a numeric column
# containing score but not n_cells/count.
SCORE_COL <- pick_exact(
  names(x),
  c(
    "program_score",
    "mean_program_score",
    "mean_score",
    "score",
    "donor_score"
  )
)

if (is.na(SCORE_COL)) {
  candidates <- names(x)[
    grepl("score", tolower(names(x))) &
      !grepl("count|cell|n_", tolower(names(x)))
  ]
  numeric_candidates <- candidates[
    vapply(x[, ..candidates], is.numeric, logical(1))
  ]
  if (length(numeric_candidates)) SCORE_COL <- numeric_candidates[1]
}

required_cols <- c(DONOR_COL, PROGRAM_COL, SCORE_COL)
if (any(is.na(required_cols))) {
  stop(
    "Could not identify donor/program/score columns.\n",
    "Detected: donor=", DONOR_COL,
    "; program=", PROGRAM_COL,
    "; score=", SCORE_COL,
    "; lineage=", LINEAGE_COL,
    "; condition=", CONDITION_COL,
    "\nSee P5_RESTART3_FIX1_SCORE_TABLE_SCHEMA.tsv"
  )
}

dt <- data.table(
  donor = as.character(x[[DONOR_COL]]),
  program = as.character(x[[PROGRAM_COL]]),
  score = suppressWarnings(as.numeric(x[[SCORE_COL]]))
)

if (!is.na(LINEAGE_COL)) {
  dt[, lineage := as.character(x[[LINEAGE_COL]])]
} else {
  dt[, lineage := sub("__.*$", "", program)]
}

if (!is.na(CONDITION_COL)) {
  dt[, condition := as.character(x[[CONDITION_COL]])]
} else {
  dt[, condition := NA_character_]
}

dt <- dt[
  !is.na(donor) & nzchar(donor) &
    !is.na(program) & nzchar(program) &
    is.finite(score)
]

# Collapse duplicates only if the upstream table contains >1 row per donor-program.
dt <- dt[
  ,
  .(
    score = mean(score, na.rm = TRUE),
    lineage = {
      z <- unique(lineage[!is.na(lineage) & nzchar(lineage)])
      if (length(z)) z[1] else sub("__.*$", "", program[1])
    },
    condition = {
      z <- unique(condition[!is.na(condition) & nzchar(condition)])
      if (length(z)) z[1] else NA_character_
    }
  ),
  by = .(donor, program)
]

fwrite(
  dt,
  file.path(OUT, "P5_RESTART3_FIX1_DONOR_PROGRAM_SCORE_LONG.tsv"),
  sep = "\t"
)

# =============================================================================
# 3. Recover P5 donor vector
# =============================================================================

p5 <- dt[program == P5]

if (!nrow(p5)) {
  stop("Oligodendrocytes__P5 is absent from donor-level score table.")
}

if (anyDuplicated(p5$donor)) {
  stop("P5 has duplicated donor rows after collapse.")
}

if (nrow(p5) < 6L) {
  stop("Fewer than 6 donor scores available for P5: n=", nrow(p5))
}

# Condition map comes primarily from P5 rows.
cond_map <- unique(p5[, .(donor, condition)])

# If P5 condition is missing, recover donor condition from all rows.
if (all(is.na(cond_map$condition) | !nzchar(cond_map$condition))) {
  cond_map <- dt[
    !is.na(condition) & nzchar(condition),
    .(condition = unique(condition)[1]),
    by = donor
  ]
}

# normalize condition only for centering; preserve original strings.
normalize_condition <- function(z) {
  y <- tolower(trimws(as.character(z)))
  out <- y
  out[grepl("tle|epilep|case", y)] <- "TLE"
  out[grepl("control|normal|ctrl|non.?epilep", y)] <- "Control"
  out
}

if (nrow(cond_map)) {
  cond_map[, condition_norm := normalize_condition(condition)]
} else {
  cond_map <- data.table(
    donor = p5$donor,
    condition = NA_character_,
    condition_norm = NA_character_
  )
}

p5 <- merge(
  p5[, .(donor, P5_score = score)],
  cond_map[, .(donor, condition, condition_norm)],
  by = "donor",
  all.x = TRUE,
  sort = FALSE
)

# =============================================================================
# 4. Cross-lineage donor correlations
# =============================================================================

partner_programs <- unique(dt$program)
partner_programs <- setdiff(partner_programs, P5)

rows <- list()

for (pr in partner_programs) {

  z <- dt[
    program == pr,
    .(
      donor,
      partner_score = score,
      partner_lineage = lineage
    )
  ]

  z <- merge(
    p5,
    z,
    by = "donor",
    all = FALSE,
    sort = FALSE
  )

  z <- z[
    is.finite(P5_score) &
      is.finite(partner_score)
  ]

  if (nrow(z) < 6L) next

  partner_lineage <- {
    a <- unique(z$partner_lineage[!is.na(z$partner_lineage)])
    if (length(a)) a[1] else sub("__.*$", "", pr)
  }

  # We want cross-lineage context only.
  if (partner_lineage == "Oligodendrocytes") next

  global_spearman <- suppressWarnings(
    cor(
      z$P5_score,
      z$partner_score,
      method = "spearman",
      use = "complete.obs"
    )
  )

  # Condition-adjusted residual correlation.
  # If usable condition labels exist, regress both scores on condition and
  # correlate residuals. Otherwise fall back to donor-global Pearson and flag it.
  condition_usable <- (
    sum(!is.na(z$condition_norm) & nzchar(z$condition_norm)) == nrow(z) &&
      uniqueN(z$condition_norm) >= 2L
  )

  if (condition_usable) {

    f <- factor(z$condition_norm)

    r1 <- residuals(lm(P5_score ~ f, data = z))
    r2 <- residuals(lm(partner_score ~ f, data = z))

    adjusted_r <- suppressWarnings(
      cor(r1, r2, method = "pearson", use = "complete.obs")
    )

    adjusted_method <- "condition_residual_pearson"

  } else {

    adjusted_r <- suppressWarnings(
      cor(
        z$P5_score,
        z$partner_score,
        method = "pearson",
        use = "complete.obs"
      )
    )

    adjusted_method <- "global_pearson_condition_unavailable"
  }

  # Leave-one-donor-out sign stability for adjusted association.
  full_sign <- sign(adjusted_r)

  loo_same <- vapply(
    seq_len(nrow(z)),
    function(i) {

      q <- z[-i]

      if (nrow(q) < 5L) return(NA)

      cond_ok <- (
        sum(!is.na(q$condition_norm) & nzchar(q$condition_norm)) == nrow(q) &&
          uniqueN(q$condition_norm) >= 2L
      )

      if (cond_ok) {
        ff <- factor(q$condition_norm)
        rr1 <- residuals(lm(P5_score ~ ff, data = q))
        rr2 <- residuals(lm(partner_score ~ ff, data = q))
        rr <- suppressWarnings(cor(rr1, rr2, method = "pearson"))
      } else {
        rr <- suppressWarnings(cor(q$P5_score, q$partner_score, method = "pearson"))
      }

      if (!is.finite(rr) || !is.finite(full_sign)) return(NA)

      if (full_sign == 0) {
        rr == 0
      } else {
        sign(rr) == full_sign
      }
    },
    logical(1)
  )

  # TLE/control-specific Spearman are descriptive only because n is tiny.
  tle_rho <- NA_real_
  control_rho <- NA_real_

  if (condition_usable) {
    if (sum(z$condition_norm == "TLE") >= 3L) {
      tle_rho <- suppressWarnings(
        cor(
          z[condition_norm == "TLE", P5_score],
          z[condition_norm == "TLE", partner_score],
          method = "spearman"
        )
      )
    }

    if (sum(z$condition_norm == "Control") >= 3L) {
      control_rho <- suppressWarnings(
        cor(
          z[condition_norm == "Control", P5_score],
          z[condition_norm == "Control", partner_score],
          method = "spearman"
        )
      )
    }
  }

  rows[[length(rows) + 1L]] <- data.table(
    partner_program = pr,
    partner_lineage = partner_lineage,
    n_donors = nrow(z),
    global_spearman_rho = global_spearman,
    adjusted_r = adjusted_r,
    adjusted_method = adjusted_method,
    LODO_same_sign_fraction = mean(loo_same, na.rm = TRUE),
    TLE_spearman_rho_descriptive = tle_rho,
    Control_spearman_rho_descriptive = control_rho
  )
}

cross <- if (length(rows)) {
  rbindlist(rows, fill = TRUE)
} else {
  data.table(
    partner_program = character(),
    partner_lineage = character(),
    n_donors = integer(),
    global_spearman_rho = numeric(),
    adjusted_r = numeric(),
    adjusted_method = character(),
    LODO_same_sign_fraction = numeric(),
    TLE_spearman_rho_descriptive = numeric(),
    Control_spearman_rho_descriptive = numeric()
  )
}

if (nrow(cross)) {
  cross[, abs_adjusted_r := abs(adjusted_r)]
  setorder(
    cross,
    -LODO_same_sign_fraction,
    -abs_adjusted_r,
    partner_program
  )
  cross[, microenvironment_rank := seq_len(.N)]
}

fwrite(
  cross,
  file.path(OUT, "P5_RESTART3_FIX1_CROSS_LINEAGE_PROGRAM_CORRELATIONS.tsv"),
  sep = "\t"
)

top10 <- head(cross, 10)

fwrite(
  top10,
  file.path(OUT, "P5_RESTART3_FIX1_TOP10_MICROENVIRONMENT_PROGRAMS.tsv"),
  sep = "\t"
)

# =============================================================================
# 5. Sender-lineage shortlist
# =============================================================================

if (nrow(cross)) {

  sender <- cross[
    ,
    .(
      n_programs = .N,
      n_top10_programs = sum(microenvironment_rank <= 10L),
      strongest_partner_program = partner_program[
        which.max(abs(adjusted_r))
      ],
      strongest_adjusted_r = adjusted_r[
        which.max(abs(adjusted_r))
      ],
      strongest_LODO_same_sign_fraction = LODO_same_sign_fraction[
        which.max(abs(adjusted_r))
      ],
      max_abs_adjusted_r = max(abs(adjusted_r), na.rm = TRUE)
    ),
    by = partner_lineage
  ]

  setorder(
    sender,
    -n_top10_programs,
    -strongest_LODO_same_sign_fraction,
    -max_abs_adjusted_r,
    partner_lineage
  )

  sender[, sender_lineage_rank := seq_len(.N)]
  sender[, next_LR_scope := sender_lineage_rank <= 3L]

} else {

  # Always write the expected file, even if empty.
  sender <- data.table(
    partner_lineage = character(),
    n_programs = integer(),
    n_top10_programs = integer(),
    strongest_partner_program = character(),
    strongest_adjusted_r = numeric(),
    strongest_LODO_same_sign_fraction = numeric(),
    max_abs_adjusted_r = numeric(),
    sender_lineage_rank = integer(),
    next_LR_scope = logical()
  )
}

fwrite(
  sender,
  file.path(OUT, "P5_RESTART3_SENDER_LINEAGE_SHORTLIST.tsv"),
  sep = "\t"
)

# =============================================================================
# 6. Corrected microenvironment decision
# =============================================================================

n_cross <- nrow(cross)
n_sender <- nrow(sender)
n_scope <- if (n_sender) sum(sender$next_LR_scope) else 0L

status <- if (n_cross > 0L && n_scope > 0L) {
  "CROSS_LINEAGE_CONTEXT_RECOVERED_FOR_TARGETED_LR"
} else {
  "NO_USABLE_CROSS_LINEAGE_CONTEXT_FROM_DONOR_SCORE_TABLE"
}

decision <- data.table(
  item = c(
    "source_table",
    "n_P5_donors",
    "n_cross_lineage_program_correlations",
    "n_sender_lineages",
    "n_sender_lineages_for_next_LR",
    "MICROENVIRONMENT_STATUS",
    "NEXT_STAGE"
  ),
  value = c(
    SCORE_FILE,
    nrow(p5),
    n_cross,
    n_sender,
    n_scope,
    status,
    ifelse(
      n_scope > 0L,
      "RESTART4_TARGETED_DONOR_AWARE_LR_USING_SHORTLISTED_SENDER_LINEAGES",
      "RESTART4_SPATIAL_REFERENCE_MAPPING_NO_LR_EXPANSION"
    )
  )
)

fwrite(
  decision,
  file.path(OUT, "P5_RESTART3_FIX1_DECISION.tsv"),
  sep = "\t"
)

writeLines(
  c(
    "P5 RESTART-3 FIX1 | CROSS-LINEAGE MICROENVIRONMENT CONTEXT",
    "===========================================================",
    "",
    "Reason for correction:",
    "- Original 34A2 correlation table was within-lineage.",
    "- Therefore filtering P5 correlations for non-oligodendrocyte partners produced zero rows by construction.",
    "- FIX1 uses the existing donor-level program-score table and computes cross-lineage donor correlations directly.",
    "",
    paste0("P5 donors=", nrow(p5)),
    paste0("cross-lineage program correlations=", n_cross),
    paste0("sender lineages recovered=", n_sender),
    paste0("sender lineages retained for next LR=", n_scope),
    paste0("MICROENVIRONMENT_STATUS=", status),
    "",
    "Interpretation boundary:",
    "- These are donor-level context correlations, not ligand-receptor interactions.",
    "- Condition-adjusted residual correlations are used when condition is available.",
    "- TLE/control within-group correlations are descriptive only because donor n is small.",
    "- The shortlist defines search scope for the next LR step; it does not establish sender causality."
  ),
  file.path(OUT, "P5_RESTART3_FIX1_REPORT.txt")
)

message("P5 RESTART-3 FIX1 completed.")
message("Expected sender file: ", file.path(OUT, "P5_RESTART3_SENDER_LINEAGE_SHORTLIST.tsv"))

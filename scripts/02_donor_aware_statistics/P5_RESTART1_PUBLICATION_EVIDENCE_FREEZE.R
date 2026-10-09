# P5_RESTART1_PUBLICATION_EVIDENCE_FREEZE.R
# Restart of the pure-bioinformatics manuscript centered on:
#   Oligodendrocytes__P5
#
# Purpose:
#   Build ONE publication-oriented evidence freeze from already completed analyses.
#   This script does NOT rerun QC, clustering, cNMF, differential expression,
#   external validation, LR, TF, spatial analysis, or any new biological test.
#
# It answers:
#   1) What exactly is already established for P5?
#   2) What was only exploratory/directional?
#   3) What mature-OL specificity / cross-reference evidence already exists?
#   4) What are the true remaining gaps before mechanism/spatial work?
#
# Important:
#   Existing frozen decisions are preserved even if they are negative.
#   No post hoc redefinition of P5 is allowed here.

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- Sys.getenv("P5_WORK_ROOT")

PROJECT <- file.path(
  ROOT,
  "P5_MATURE_OLIGODENDROCYTE_PURE_BIOINFO_RESTART"
)

OUT <- file.path(
  PROJECT,
  "RESTART1_PUBLICATION_EVIDENCE_FREEZE"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

RELEASE_ROOT <- file.path(ROOT, "99_frozen_releases")
CROSSREF_ROOT <- file.path(ROOT, "35C_P5_CROSS_REFERENCE")

P5 <- "Oligodendrocytes__P5"

# =============================================================================
# Helpers
# =============================================================================

latest_complete_release <- function(prefix, required_names) {

  dirs <- list.dirs(
    RELEASE_ROOT,
    full.names = TRUE,
    recursive = FALSE
  )

  dirs <- dirs[
    grepl(
      paste0("^", prefix),
      basename(dirs),
      ignore.case = FALSE
    )
  ]

  if (!length(dirs)) {
    return(NA_character_)
  }

  info <- file.info(dirs)
  dirs <- dirs[order(info$mtime, decreasing = TRUE)]

  for (d in dirs) {

    files <- list.files(
      d,
      recursive = TRUE,
      full.names = TRUE
    )

    ok <- vapply(
      required_names,
      function(nm) sum(basename(files) == nm) == 1L,
      logical(1)
    )

    if (all(ok)) {
      return(normalizePath(d, winslash = "/", mustWork = TRUE))
    }
  }

  NA_character_
}

one_file <- function(root, name, required = TRUE) {

  if (is.na(root) || !dir.exists(root)) {
    if (required) stop("Root directory unavailable for: ", name)
    return(NA_character_)
  }

  hits <- list.files(
    root,
    recursive = TRUE,
    full.names = TRUE
  )

  hits <- hits[basename(hits) == name]

  if (length(hits) == 1L) {
    return(hits)
  }

  if (required) {
    stop(
      "Could not uniquely resolve file: ", name,
      "\nHits:\n", paste(hits, collapse = "\n")
    )
  }

  NA_character_
}

read_any_table <- function(path) {

  if (is.na(path) || !file.exists(path)) {
    return(NULL)
  }

  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("csv", "tsv", "txt")) {
    sep <- if (ext == "csv") "," else "\t"
    x <- tryCatch(
      fread(path, sep = sep),
      error = function(e) {
        tryCatch(fread(path), error = function(e2) NULL)
      }
    )
    return(x)
  }

  NULL
}

contains_p5 <- function(dt) {

  if (is.null(dt) || !nrow(dt)) return(logical(0))

  char_cols <- names(dt)[
    vapply(dt, function(x) is.character(x) || is.factor(x), logical(1))
  ]

  if (!length(char_cols)) {
    return(rep(FALSE, nrow(dt)))
  }

  hit <- rep(FALSE, nrow(dt))

  for (cc in char_cols) {
    v <- as.character(dt[[cc]])
    hit <- hit |
      v == P5 |
      grepl("Oligodendrocytes__P5", v, fixed = TRUE)
  }

  hit
}

extract_p5_rows <- function(path, source_name) {

  dt <- read_any_table(path)

  if (is.null(dt)) {
    return(data.table(
      source = source_name,
      file = path,
      extraction_status = "FILE_UNREADABLE_OR_NON_TABULAR"
    ))
  }

  idx <- contains_p5(dt)

  if (!length(idx) || !any(idx)) {
    return(data.table(
      source = source_name,
      file = path,
      extraction_status = "NO_P5_ROW_FOUND"
    ))
  }

  out <- copy(dt[idx])
  out[, source := source_name]
  out[, file := path]
  out[, extraction_status := "P5_ROW_RECOVERED"]

  setcolorder(
    out,
    c(
      "source",
      "file",
      "extraction_status",
      setdiff(names(out), c("source", "file", "extraction_status"))
    )
  )

  out
}

write_if <- function(x, path) {
  if (!is.null(x)) {
    fwrite(x, path, sep = "\t")
  }
}

safe_bool <- function(x) {

  if (length(x) == 0 || is.na(x[1])) return(NA)

  if (is.logical(x)) return(as.logical(x[1]))

  z <- tolower(trimws(as.character(x[1])))

  if (z %in% c("true", "t", "1", "yes")) return(TRUE)
  if (z %in% c("false", "f", "0", "no")) return(FALSE)

  NA
}

# =============================================================================
# 1. Resolve authoritative frozen releases
# =============================================================================

B3A <- latest_complete_release(
  "32B3A_UNBIASED_PROGRAM_FREEZE_",
  c(
    "32B3A_PROGRAM_TLE_EFFECT_SUMMARY.csv",
    "32B3A_FINAL_PROGRAM_LOCK.csv",
    "32B3A_TOP100_GENES.csv"
  )
)

B3B0 <- latest_complete_release(
  "32B3B0_PROGRAM_CONTENT_AUDIT_",
  c(
    "32B3B0_PROGRAM_CONTENT_AUDIT.csv"
  )
)

B3B1B <- latest_complete_release(
  "32B3B1B_PATHWAY_AND_26E_INTERPRETATION_",
  c(
    "32B3B1B_PROGRAM_INTERPRETATION.csv",
    "32B3B1B_26E_OVERLAP_WITH_DISEASE_DIRECTION.csv"
  )
)

B3C <- latest_complete_release(
  "32B3C_MANUSCRIPT_SAFE_PROGRAM_FREEZE_",
  c(
    "32B3C_MANUSCRIPT_SAFE_PROGRAM_LABELS.csv",
    "32B3C_CLAIM_MATRIX.csv",
    "32B3C_VALIDATION_PRIORITY.csv"
  )
)

A0 <- latest_complete_release(
  "34A0_GLIAL_STATE_CONSENSUS_INPUT_AUDIT_",
  c(
    "34A0_MASTER_PROGRAM_EVIDENCE_MATRIX.csv",
    "34A0_STATE_CONSENSUS_LOCK.json"
  )
)

resolved_releases <- data.table(
  evidence_layer = c(
    "32B3A_unbiased_program_freeze",
    "32B3B0_content_audit",
    "32B3B1B_pathway_interpretation",
    "32B3C_manuscript_safe_freeze",
    "34A0_integrated_consensus"
  ),
  path = c(B3A, B3B0, B3B1B, B3C, A0),
  available = !is.na(c(B3A, B3B0, B3B1B, B3C, A0))
)

fwrite(
  resolved_releases,
  file.path(OUT, "P5_RESTART1_RESOLVED_FROZEN_RELEASES.tsv"),
  sep = "\t"
)

# =============================================================================
# 2. Recover discovery / donor-aware P5 evidence
# =============================================================================

EFFECT_FILE <- one_file(
  B3A,
  "32B3A_PROGRAM_TLE_EFFECT_SUMMARY.csv"
)

LOCK_FILE <- one_file(
  B3A,
  "32B3A_FINAL_PROGRAM_LOCK.csv"
)

TOP100_FILE <- one_file(
  B3A,
  "32B3A_TOP100_GENES.csv"
)

DONOR_EFFECT_FILE <- one_file(
  B3A,
  "32B3A_ALL_DONOR_EFFECTS.csv",
  required = FALSE
)

DONOR_USAGE_FILE <- one_file(
  B3A,
  "32B3A_DONOR_USAGE_LONG.csv",
  required = FALSE
)

effects <- fread(EFFECT_FILE)
lock <- fread(LOCK_FILE)
top100 <- fread(TOP100_FILE)

p5_effect <- effects[contains_p5(effects)]
p5_lock <- lock[contains_p5(lock)]
p5_top100 <- top100[contains_p5(top100)]

if (nrow(p5_effect) != 1L) {
  stop(
    "Expected exactly one P5 row in 32B3A_PROGRAM_TLE_EFFECT_SUMMARY.csv; found ",
    nrow(p5_effect)
  )
}

if (nrow(p5_lock) != 1L) {
  stop(
    "Expected exactly one P5 row in 32B3A_FINAL_PROGRAM_LOCK.csv; found ",
    nrow(p5_lock)
  )
}

if (!nrow(p5_top100)) {
  stop("Could not recover P5 Top100 genes from 32B3A_TOP100_GENES.csv.")
}

# Preserve original ranking if available.
rank_candidates <- c(
  "rank", "gene_rank", "program_rank", "meta_rank", "loading_rank"
)
rank_col <- rank_candidates[rank_candidates %in% names(p5_top100)]

if (length(rank_col)) {
  setorderv(p5_top100, rank_col[1])
}

fwrite(
  p5_effect,
  file.path(OUT, "P5_RESTART1_DISCOVERY_EFFECT_ROW.tsv"),
  sep = "\t"
)

fwrite(
  p5_lock,
  file.path(OUT, "P5_RESTART1_PROGRAM_LOCK_ROW.tsv"),
  sep = "\t"
)

fwrite(
  p5_top100,
  file.path(OUT, "P5_RESTART1_FROZEN_TOP100_GENES.tsv"),
  sep = "\t"
)

p5_donor_effect <- NULL
if (!is.na(DONOR_EFFECT_FILE)) {
  donor_effect <- fread(DONOR_EFFECT_FILE)
  p5_donor_effect <- donor_effect[contains_p5(donor_effect)]
  write_if(
    p5_donor_effect,
    file.path(OUT, "P5_RESTART1_DONOR_EFFECTS.tsv")
  )
}

p5_donor_usage <- NULL
if (!is.na(DONOR_USAGE_FILE)) {
  donor_usage <- fread(DONOR_USAGE_FILE)
  p5_donor_usage <- donor_usage[contains_p5(donor_usage)]
  write_if(
    p5_donor_usage,
    file.path(OUT, "P5_RESTART1_DONOR_USAGE.tsv")
  )
}

# =============================================================================
# 3. Recover functional/content interpretation without re-enrichment
# =============================================================================

CONTENT_FILE <- one_file(
  B3B0,
  "32B3B0_PROGRAM_CONTENT_AUDIT.csv"
)

INTERPRET_FILE <- one_file(
  B3B1B,
  "32B3B1B_PROGRAM_INTERPRETATION.csv"
)

OVERLAP_FILE <- one_file(
  B3B1B,
  "32B3B1B_26E_OVERLAP_WITH_DISEASE_DIRECTION.csv"
)

LABEL_FILE <- one_file(
  B3C,
  "32B3C_MANUSCRIPT_SAFE_PROGRAM_LABELS.csv"
)

CLAIM_FILE <- one_file(
  B3C,
  "32B3C_CLAIM_MATRIX.csv"
)

PRIORITY_FILE <- one_file(
  B3C,
  "32B3C_VALIDATION_PRIORITY.csv"
)

p5_content <- extract_p5_rows(
  CONTENT_FILE,
  "32B3B0_CONTENT_AUDIT"
)

p5_interpretation <- extract_p5_rows(
  INTERPRET_FILE,
  "32B3B1B_PROGRAM_INTERPRETATION"
)

p5_overlap <- extract_p5_rows(
  OVERLAP_FILE,
  "32B3B1B_26E_OVERLAP"
)

p5_label <- extract_p5_rows(
  LABEL_FILE,
  "32B3C_MANUSCRIPT_SAFE_LABEL"
)

p5_claim <- extract_p5_rows(
  CLAIM_FILE,
  "32B3C_CLAIM_MATRIX"
)

p5_priority <- extract_p5_rows(
  PRIORITY_FILE,
  "32B3C_VALIDATION_PRIORITY"
)

fwrite(
  p5_content,
  file.path(OUT, "P5_RESTART1_CONTENT_AUDIT.tsv"),
  sep = "\t"
)

fwrite(
  p5_interpretation,
  file.path(OUT, "P5_RESTART1_FUNCTIONAL_INTERPRETATION.tsv"),
  sep = "\t"
)

fwrite(
  p5_overlap,
  file.path(OUT, "P5_RESTART1_26E_OVERLAP.tsv"),
  sep = "\t"
)

fwrite(
  p5_label,
  file.path(OUT, "P5_RESTART1_MANUSCRIPT_SAFE_LABEL.tsv"),
  sep = "\t"
)

fwrite(
  p5_claim,
  file.path(OUT, "P5_RESTART1_CLAIM_MATRIX.tsv"),
  sep = "\t"
)

fwrite(
  p5_priority,
  file.path(OUT, "P5_RESTART1_VALIDATION_PRIORITY.tsv"),
  sep = "\t"
)

# =============================================================================
# 4. Recover integrated 34A0 evidence boundary
# =============================================================================

p5_a0 <- data.table(
  extraction_status = "34A0_NOT_AVAILABLE"
)

if (!is.na(A0)) {

  A0_MASTER <- one_file(
    A0,
    "34A0_MASTER_PROGRAM_EVIDENCE_MATRIX.csv"
  )

  a0 <- fread(A0_MASTER)
  p5_a0 <- a0[contains_p5(a0)]

  if (nrow(p5_a0) != 1L) {
    stop(
      "Expected exactly one P5 row in 34A0 master matrix; found ",
      nrow(p5_a0)
    )
  }
}

fwrite(
  p5_a0,
  file.path(OUT, "P5_RESTART1_34A0_INTEGRATED_EVIDENCE.tsv"),
  sep = "\t"
)

# =============================================================================
# 5. Recover existing 35C mature-OL / platform-specificity evidence
# =============================================================================

SPEC_ROOT <- file.path(
  CROSSREF_ROOT,
  "09_PROGRAM_SPECIFICITY_AND_PLATFORM_SENSITIVITY"
)

P5_SPEC_FILE <- file.path(
  SPEC_ROOT,
  "02_results",
  "P5_PROGRAM_SPECIFICITY_SUMMARY.csv"
)

ALL_COMP_FILE <- file.path(
  SPEC_ROOT,
  "02_results",
  "ALL_PROGRAM_COMPARISON_RESULTS.csv"
)

MATURE_RANK_FILE <- file.path(
  SPEC_ROOT,
  "02_results",
  "MATURE_ADULT_3_VS_3_PROGRAM_SPECIFICITY_RANKING.csv"
)

DONOR_SCORE_35C <- file.path(
  SPEC_ROOT,
  "01_program_scores",
  "ALL_OLIGODENDROCYTE_PROGRAM_DONOR_SCORES.csv"
)

SUMMARY_35C <- file.path(
  SPEC_ROOT,
  "35C_06_PROGRAM_SPECIFICITY_SUMMARY.md"
)

p5_specificity <- extract_p5_rows(
  P5_SPEC_FILE,
  "35C_P5_PROGRAM_SPECIFICITY_SUMMARY"
)

p5_all_comparisons <- extract_p5_rows(
  ALL_COMP_FILE,
  "35C_ALL_PROGRAM_COMPARISONS"
)

p5_mature_ranking <- extract_p5_rows(
  MATURE_RANK_FILE,
  "35C_MATURE_ADULT_3_VS_3_RANKING"
)

p5_35c_donor_scores <- extract_p5_rows(
  DONOR_SCORE_35C,
  "35C_ALL_OLIGODENDROCYTE_PROGRAM_DONOR_SCORES"
)

fwrite(
  p5_specificity,
  file.path(OUT, "P5_RESTART1_35C_SPECIFICITY_SUMMARY.tsv"),
  sep = "\t"
)

fwrite(
  p5_all_comparisons,
  file.path(OUT, "P5_RESTART1_35C_ALL_COMPARISONS.tsv"),
  sep = "\t"
)

fwrite(
  p5_mature_ranking,
  file.path(OUT, "P5_RESTART1_35C_MATURE_RANKING.tsv"),
  sep = "\t"
)

fwrite(
  p5_35c_donor_scores,
  file.path(OUT, "P5_RESTART1_35C_DONOR_SCORES.tsv"),
  sep = "\t"
)

if (file.exists(SUMMARY_35C)) {
  file.copy(
    SUMMARY_35C,
    file.path(OUT, "P5_RESTART1_35C_ORIGINAL_SUMMARY.md"),
    overwrite = TRUE
  )
}

# =============================================================================
# 6. Inventory existing human external/reference resources
#    Presence only; NO reanalysis and NO claims are inferred from file existence.
# =============================================================================

resource_paths <- c(
  GSE140393_MATURE_OLIGODENDROCYTES = file.path(
    CROSSREF_ROOT,
    "06_GSE140393_FINAL_MANUAL_OLIGO_ANNOTATION",
    "objects",
    "GSE140393_MATURE_OLIGODENDROCYTES.rds"
  ),
  GSE140393_BONA_FIDE_OLIGODENDROGLIA = file.path(
    CROSSREF_ROOT,
    "06_GSE140393_FINAL_MANUAL_OLIGO_ANNOTATION",
    "objects",
    "GSE140393_BONA_FIDE_OLIGODENDROGLIA.rds"
  ),
  SEAAD_OLIGODENDROCYTE_FORMAL_RAWCOUNTS = file.path(
    CROSSREF_ROOT,
    "02_SEAAD_TARGET_DONOR_SUBSETS",
    "SEAAD_OLIGODENDROCYTE_FORMAL_RAWCOUNTS.h5ad"
  ),
  SEAAD_OLIGODENDROCYTE_TARGET5_RAWCOUNTS = file.path(
    CROSSREF_ROOT,
    "02_SEAAD_TARGET_DONOR_SUBSETS",
    "SEAAD_OLIGODENDROCYTE_TARGET5_RAWCOUNTS.h5ad"
  )
)

resource_inventory <- data.table(
  resource = names(resource_paths),
  path = unname(resource_paths),
  exists = file.exists(resource_paths),
  role = c(
    "independent_human_TLE_mature_OL_context",
    "independent_human_TLE_oligodendroglia_context",
    "adult_MTG_mature_OL_reference",
    "adult_MTG_target_donor_reference"
  )
)

fwrite(
  resource_inventory,
  file.path(OUT, "P5_RESTART1_EXISTING_RESOURCE_INVENTORY.tsv"),
  sep = "\t"
)

# =============================================================================
# 7. Publication-oriented master evidence registry
# =============================================================================

formal_flag <- NA
directional_flag <- NA

if ("formal_TLE_program_candidate" %in% names(p5_effect)) {
  formal_flag <- safe_bool(p5_effect$formal_TLE_program_candidate)
}

if ("directional_support_only" %in% names(p5_effect)) {
  directional_flag <- safe_bool(p5_effect$directional_support_only)
}

a0_status <- NA_character_
if ("program_use_status" %in% names(p5_a0)) {
  a0_status <- as.character(p5_a0$program_use_status[1])
}

safe_label <- NA_character_
if ("safe_label" %in% names(p5_a0)) {
  safe_label <- as.character(p5_a0$safe_label[1])
} else if ("safe_label" %in% names(p5_label)) {
  safe_label <- as.character(p5_label$safe_label[1])
}

master_evidence <- data.table(
  evidence_block = c(
    "P5_program_exists_and_is_frozen",
    "formal_TLE_program_status",
    "directional_TLE_support",
    "integrated_program_use_status",
    "manuscript_safe_label_available",
    "frozen_top100_available",
    "donor_effect_table_available",
    "donor_usage_table_available",
    "content_annotation_available",
    "35C_mature_specificity_results_available",
    "GSE140393_mature_OL_resource_available",
    "SEAAD_adult_MTG_reference_available"
  ),
  status = c(
    nrow(p5_effect) == 1L && nrow(p5_lock) == 1L,
    formal_flag,
    directional_flag,
    !is.na(a0_status),
    !is.na(safe_label) && nzchar(safe_label),
    nrow(p5_top100) > 0,
    !is.null(p5_donor_effect) && nrow(p5_donor_effect) > 0,
    !is.null(p5_donor_usage) && nrow(p5_donor_usage) > 0,
    any(p5_content$extraction_status == "P5_ROW_RECOVERED"),
    any(p5_specificity$extraction_status == "P5_ROW_RECOVERED") ||
      any(p5_all_comparisons$extraction_status == "P5_ROW_RECOVERED"),
    resource_inventory[
      resource == "GSE140393_MATURE_OLIGODENDROCYTES",
      exists
    ],
    resource_inventory[
      resource == "SEAAD_OLIGODENDROCYTE_FORMAL_RAWCOUNTS",
      exists
    ]
  ),
  detail = c(
    P5,
    ifelse(
      is.na(formal_flag),
      "not_recovered",
      ifelse(formal_flag, "formal_TLE_program_candidate=TRUE",
             "formal_TLE_program_candidate=FALSE")
    ),
    ifelse(
      is.na(directional_flag),
      "not_recovered",
      ifelse(directional_flag, "directional_support_only=TRUE",
             "directional_support_only=FALSE")
    ),
    ifelse(is.na(a0_status), "not_available", a0_status),
    ifelse(is.na(safe_label), "not_available", safe_label),
    paste0("rows_recovered=", nrow(p5_top100)),
    paste0("rows_recovered=", ifelse(is.null(p5_donor_effect), 0, nrow(p5_donor_effect))),
    paste0("rows_recovered=", ifelse(is.null(p5_donor_usage), 0, nrow(p5_donor_usage))),
    paste0("content_status=", p5_content$extraction_status[1]),
    paste0("35C_status=", p5_specificity$extraction_status[1]),
    "file_presence_only_not_result_claim",
    "file_presence_only_not_result_claim"
  )
)

fwrite(
  master_evidence,
  file.path(OUT, "P5_RESTART1_MASTER_EVIDENCE_REGISTRY.tsv"),
  sep = "\t"
)

# =============================================================================
# 8. Freeze the true remaining manuscript gaps
# =============================================================================

gaps <- data.table(
  module = c(
    "FOUNDATION",
    "FUNCTIONAL_ARCHITECTURE",
    "REGULATORY_MECHANISM",
    "MICROENVIRONMENT_LR",
    "SPATIAL_LOCALIZATION",
    "EXTERNAL_HUMAN_INTEGRATION",
    "MANUSCRIPT"
  ),
  current_status = c(
    ifelse(
      identical(formal_flag, TRUE),
      "FORMAL_DISEASE_ASSOCIATION_AVAILABLE",
      "EXPLORATORY_OR_DIRECTIONAL_PROGRAM; DO_NOT_OVERSTATE_DISEASE_ASSOCIATION"
    ),
    "PARTIAL_EXISTING_ANNOTATION; NEED_COMPACT_P5_FUNCTIONAL_DECOMPOSITION_FOR_PAPER",
    "NOT_YET_FROZEN_FOR_P5_MANUSCRIPT",
    "NOT_YET_FROZEN_FOR_P5_MANUSCRIPT",
    "NOT_YET_FROZEN_FOR_P5_MANUSCRIPT",
    "EXISTING_GSE140393_AND_SEAAD_RESOURCES; RESULTS_MUST_BE_INTEGRATED_WITH_CORRECT_ROLE",
    "NOT_STARTED_AS_STANDALONE_P5_MANUSCRIPT"
  ),
  next_action = c(
    "Preserve 32B3A/34A0 boundary; no redefinition of P5.",
    "Next analysis after review: one compact functional decomposition of frozen P5 Top100/loading.",
    "Only after functional axes are frozen, prioritize a small number of regulator-associated programs.",
    "Only after functional axes are frozen, perform donor-aware sender/receptor prioritization; keep 5-10 axes.",
    "Project frozen P5 program to appropriate adult cortical spatial/reference data; distinguish TLE from reference mapping.",
    "Use GSE140393 as TLE context and SEA-AD as adult MTG reference; do not treat SEA-AD as TLE replication.",
    "Build figures only after foundation + functional architecture are frozen."
  )
)

fwrite(
  gaps,
  file.path(OUT, "P5_RESTART1_PUBLICATION_GAP_MATRIX.tsv"),
  sep = "\t"
)

# =============================================================================
# 9. Compact restart report
# =============================================================================

status_line <- if (identical(formal_flag, TRUE)) {
  "FORMAL_TLE_ASSOCIATED_P5_PROGRAM"
} else {
  "P5_REMAINS_HIGH_PRIORITY_EXPLORATORY_DIRECTIONAL_PROGRAM"
}

report <- c(
  "P5 MATURE OLIGODENDROCYTE PURE-BIOINFORMATICS MANUSCRIPT",
  "RESTART-1 | PUBLICATION EVIDENCE FREEZE",
  "========================================================",
  "",
  paste0("Primary frozen program: ", P5),
  paste0("Foundation status: ", status_line),
  paste0(
    "32B3A formal_TLE_program_candidate=",
    ifelse(is.na(formal_flag), "NA", formal_flag)
  ),
  paste0(
    "32B3A directional_support_only=",
    ifelse(is.na(directional_flag), "NA", directional_flag)
  ),
  paste0(
    "34A0 program_use_status=",
    ifelse(is.na(a0_status), "NA", a0_status)
  ),
  paste0(
    "manuscript_safe_label=",
    ifelse(is.na(safe_label), "NA", safe_label)
  ),
  paste0("Frozen P5 Top100 rows recovered=", nrow(p5_top100)),
  paste0(
    "P5 donor-effect rows recovered=",
    ifelse(is.null(p5_donor_effect), 0, nrow(p5_donor_effect))
  ),
  paste0(
    "P5 donor-usage rows recovered=",
    ifelse(is.null(p5_donor_usage), 0, nrow(p5_donor_usage))
  ),
  paste0(
    "35C P5 specificity row recovered=",
    any(p5_specificity$extraction_status == "P5_ROW_RECOVERED")
  ),
  paste0(
    "GSE140393 mature OL object present=",
    resource_inventory[
      resource == "GSE140393_MATURE_OLIGODENDROCYTES",
      exists
    ]
  ),
  paste0(
    "SEA-AD adult MTG OL reference present=",
    resource_inventory[
      resource == "SEAAD_OLIGODENDROCYTE_FORMAL_RAWCOUNTS",
      exists
    ]
  ),
  "",
  "PUBLICATION INTERPRETATION BOUNDARY",
  "-----------------------------------",
  "- This stage performs no new biological inference.",
  "- Existing negative or exploratory frozen decisions are preserved.",
  "- P5 is not automatically a senescent-oligodendrocyte label.",
  "- GSE190452 tissue should remain described as neocortex unless a more specific frozen source supports otherwise.",
  "- GSE140393 is independent human TLE context; SEA-AD is adult MTG reference, not TLE replication.",
  "- Mouse hippocampal evidence is not required for this restarted pure-bioinformatics manuscript.",
  "",
  "NEXT DECISION",
  "-------------",
  "Review the recovered P5 foundation before any new mechanism analysis.",
  "If the frozen evidence supports manuscript continuation, the next substantive step is ONE compact functional decomposition of P5 (myelin/lipid, lysosome-proteostasis, metabolic/oxidative stress, senescence-associated remodeling), followed by integrated TF/LR prioritization rather than repeated validation."
)

writeLines(
  report,
  file.path(OUT, "P5_RESTART1_REPORT.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(OUT, "P5_RESTART1_SESSIONINFO.txt")
)

message("P5 RESTART-1 completed.")
message("Output: ", OUT)

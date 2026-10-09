# P5_RESTART2_FUNCTIONAL_ARCHITECTURE.R
# Pure-bioinformatics manuscript: Oligodendrocytes__P5
#
# This is the first NEW substantive analysis after RESTART-1.
# It performs one compact task:
#   1) repair/recover the existing 35C P5-specificity outputs;
#   2) decompose the frozen P5 program against the already-frozen 26E catalog;
#   3) freeze a manuscript-safe functional architecture.
#
# It does NOT rerun QC, clustering, cNMF, disease testing, external validation,
# TF analysis, LR analysis, spatial analysis, trajectory, or mouse validation.
#
# Statistical meaning:
# - Hypergeometric enrichment maps frozen gene-set content onto the frozen P5
#   top100 program. It is a PROGRAM-ANNOTATION test, not a disease-effect test.
# - Existing disease-association boundary from RESTART-1 is preserved:
#   P5 is directional/high-priority exploratory, not a formal TLE program.
# - Senescence remains "partial senescence-associated transcriptional features"
#   and must not be upgraded to a senescent-cell label.

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

R1 <- file.path(
  PROJECT,
  "RESTART1_PUBLICATION_EVIDENCE_FREEZE"
)

OUT <- file.path(
  PROJECT,
  "RESTART2_FUNCTIONAL_ARCHITECTURE"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

RELEASE_ROOT <- file.path(ROOT, "99_frozen_releases")
CROSSREF_ROOT <- file.path(ROOT, "35C_P5_CROSS_REFERENCE")

P5 <- "Oligodendrocytes__P5"

# =============================================================================
# 0. Essential frozen inputs
# =============================================================================

R1_EFFECT <- file.path(R1, "P5_RESTART1_DISCOVERY_EFFECT_ROW.tsv")
R1_A0 <- file.path(R1, "P5_RESTART1_34A0_INTEGRATED_EVIDENCE.tsv")
R1_TOP100 <- file.path(R1, "P5_RESTART1_FROZEN_TOP100_GENES.tsv")

required <- c(R1_EFFECT, R1_A0, R1_TOP100)
missing <- required[!file.exists(required)]

if (length(missing)) {
  stop(
    "Missing RESTART-1 input(s):\n",
    paste(missing, collapse = "\n")
  )
}

effect <- fread(R1_EFFECT)
a0 <- fread(R1_A0)
top100_raw <- fread(R1_TOP100)

if (nrow(effect) != 1L || nrow(a0) != 1L) {
  stop("RESTART-1 effect or integrated-evidence table is not one-row P5 evidence.")
}

# Preserve the disease-association boundary.
if ("formal_TLE_program_candidate" %in% names(effect)) {
  if (isTRUE(as.logical(effect$formal_TLE_program_candidate[1]))) {
    warning("RESTART-1 now reports formal_TLE_program_candidate=TRUE; review upstream freeze.")
  }
}

# =============================================================================
# Helpers
# =============================================================================

norm_gene <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z[is.na(z) | z == "" | z == "NA"] <- NA_character_
  z
}

pick_gene_col <- function(nms) {
  low <- tolower(nms)
  pref <- c("gene", "gene_symbol", "symbol", "feature", "human_gene")
  for (p in pref) {
    hit <- which(low == p)
    if (length(hit)) return(nms[hit[1]])
  }
  hit <- grep("gene|symbol|feature", low)
  if (length(hit)) return(nms[hit[1]])
  stop("Could not identify gene column in: ", paste(nms, collapse = ", "))
}

pick_score_col <- function(nms) {
  low <- tolower(nms)
  pref <- c(
    "meta_program_score",
    "loading",
    "score",
    "weight",
    "consensus_score",
    "program_score"
  )
  for (p in pref) {
    hit <- which(low == p)
    if (length(hit)) return(nms[hit[1]])
  }
  hit <- grep("score|loading|weight", low)
  if (length(hit)) return(nms[hit[1]])
  NA_character_
}

contains_p5_anywhere <- function(dt) {
  if (is.null(dt) || !nrow(dt)) return(rep(FALSE, 0))
  hit <- rep(FALSE, nrow(dt))
  for (cc in names(dt)) {
    if (is.character(dt[[cc]]) || is.factor(dt[[cc]])) {
      z <- as.character(dt[[cc]])
      hit <- hit |
        z == P5 |
        z == "P5" |
        grepl("Oligodendrocytes__P5", z, fixed = TRUE) |
        grepl("(^|[^A-Za-z0-9])P5([^A-Za-z0-9]|$)", z)
    }
  }
  hit
}

find_exact_files <- function(root, basename_target) {
  if (!dir.exists(root)) return(character())
  z <- list.files(root, recursive = TRUE, full.names = TRUE)
  z[basename(z) == basename_target]
}

select_latest_file <- function(paths) {
  paths <- unique(paths[file.exists(paths)])
  if (!length(paths)) return(NA_character_)
  info <- file.info(paths)
  paths[order(info$mtime, decreasing = TRUE)][1]
}

split_semicolon <- function(x) {
  if (length(x) == 0 || is.na(x[1]) || !nzchar(x[1])) return(character())
  z <- unlist(strsplit(as.character(x[1]), ";", fixed = TRUE))
  trimws(z[nzchar(trimws(z))])
}

# =============================================================================
# 1. Correctly recover 35C P5-specificity results
#
# RESTART-1 searched for the literal string "P5" inside a file whose filename
# already declares it to be P5-specific. Therefore "NO_P5_ROW_FOUND" is an
# extraction artifact, not a biological negative result.
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

DONOR_SCORE_FILE <- file.path(
  SPEC_ROOT,
  "01_program_scores",
  "ALL_OLIGODENDROCYTE_PROGRAM_DONOR_SCORES.csv"
)

if (!file.exists(P5_SPEC_FILE)) {
  stop("Existing 35C P5 summary file not found: ", P5_SPEC_FILE)
}

p5_spec <- fread(P5_SPEC_FILE)
p5_spec[, recovered_as := "WHOLE_FILE_IS_P5_SPECIFIC_BY_FROZEN_FILENAME"]

fwrite(
  p5_spec,
  file.path(OUT, "P5_RESTART2_35C_P5_SPECIFICITY_RECOVERED.tsv"),
  sep = "\t"
)

recover_p5_rows <- function(path, label) {

  if (!file.exists(path)) {
    return(data.table(
      source = label,
      recovery_status = "FILE_MISSING"
    ))
  }

  x <- fread(path)
  idx <- contains_p5_anywhere(x)

  if (length(idx) && any(idx)) {
    y <- copy(x[idx])
    y[, source := label]
    y[, recovery_status := "P5_ROWS_RECOVERED"]
    return(y)
  }

  # Do not pretend that no literal label means a negative result.
  # Preserve the full table for inspection if it is small; otherwise record schema.
  if (nrow(x) <= 200L) {
    y <- copy(x)
    y[, source := label]
    y[, recovery_status := "NO_LITERAL_P5_LABEL_FULL_TABLE_PRESERVED"]
    return(y)
  }

  data.table(
    source = label,
    recovery_status = "NO_LITERAL_P5_LABEL",
    n_rows = nrow(x),
    columns = paste(names(x), collapse = ";")
  )
}

p5_allcomp <- recover_p5_rows(
  ALL_COMP_FILE,
  "35C_ALL_PROGRAM_COMPARISON_RESULTS"
)

p5_mature_rank <- recover_p5_rows(
  MATURE_RANK_FILE,
  "35C_MATURE_ADULT_3_VS_3_PROGRAM_SPECIFICITY_RANKING"
)

p5_donor_scores <- recover_p5_rows(
  DONOR_SCORE_FILE,
  "35C_ALL_OLIGODENDROCYTE_PROGRAM_DONOR_SCORES"
)

fwrite(
  p5_allcomp,
  file.path(OUT, "P5_RESTART2_35C_ALL_COMPARISONS_RECOVERED.tsv"),
  sep = "\t"
)

fwrite(
  p5_mature_rank,
  file.path(OUT, "P5_RESTART2_35C_MATURE_RANKING_RECOVERED.tsv"),
  sep = "\t"
)

fwrite(
  p5_donor_scores,
  file.path(OUT, "P5_RESTART2_35C_DONOR_SCORES_RECOVERED.tsv"),
  sep = "\t"
)

# =============================================================================
# 2. Recover frozen P5 Top100 and full meta-spectrum
# =============================================================================

gene_col_top <- pick_gene_col(names(top100_raw))
top100_genes <- unique(na.omit(norm_gene(top100_raw[[gene_col_top]])))

if (length(top100_genes) != 100L) {
  stop(
    "Expected exactly 100 unique frozen P5 Top100 genes; found ",
    length(top100_genes)
  )
}

SPECTRUM_BASENAME <- "Oligodendrocytes__P5__meta_spectrum.csv"

spectrum_hits <- unique(c(
  find_exact_files(file.path(ROOT, "04_tables"), SPECTRUM_BASENAME),
  find_exact_files(RELEASE_ROOT, SPECTRUM_BASENAME),
  find_exact_files(ROOT, SPECTRUM_BASENAME)
))

SPECTRUM_FILE <- select_latest_file(spectrum_hits)

if (is.na(SPECTRUM_FILE)) {
  stop("Could not locate frozen P5 full meta-spectrum: ", SPECTRUM_BASENAME)
}

spectrum <- fread(SPECTRUM_FILE)

gene_col_spec <- pick_gene_col(names(spectrum))
score_col_spec <- pick_score_col(names(spectrum))

if (is.na(score_col_spec)) {
  stop("Could not identify P5 meta-spectrum loading/score column.")
}

spectrum[, gene := norm_gene(get(gene_col_spec))]
spectrum[, meta_program_score := suppressWarnings(as.numeric(get(score_col_spec)))]
spectrum <- spectrum[
  !is.na(gene) &
    nzchar(gene) &
    is.finite(meta_program_score)
]

setorder(spectrum, -meta_program_score)
spectrum <- spectrum[!duplicated(gene)]
spectrum[, spectrum_rank := seq_len(.N)]

if (nrow(spectrum) < 5000L) {
  stop("Unexpectedly small P5 meta-spectrum: ", nrow(spectrum), " genes.")
}

spectrum[, in_frozen_top100 := gene %in% top100_genes]

fwrite(
  spectrum[in_frozen_top100 == TRUE],
  file.path(OUT, "P5_RESTART2_TOP100_WITH_LOADINGS.tsv"),
  sep = "\t"
)

# =============================================================================
# 3. Recover exact frozen 26E catalog
# =============================================================================

CATALOG_BASENAME <- "32B3B1A_FROZEN_26E_GENESETS_LONG.csv"

catalog_hits <- find_exact_files(RELEASE_ROOT, CATALOG_BASENAME)

if (!length(catalog_hits)) {
  stop("Could not locate exact frozen 26E catalog: ", CATALOG_BASENAME)
}

# Prefer a release with the frozen-ready lock; otherwise use latest valid catalog.
catalog_info <- data.table(
  path = unique(catalog_hits)
)

catalog_info[, mtime := as.numeric(file.info(path)$mtime)]
catalog_info[, ready_lock := FALSE]

for (i in seq_len(nrow(catalog_info))) {
  parent <- dirname(catalog_info$path[i])
  lock_hits <- list.files(
    parent,
    pattern = "^32B3B1A_CATALOG_LOCK\\.json$",
    full.names = TRUE
  )
  if (length(lock_hits) == 1L) {
    txt <- paste(readLines(lock_hits, warn = FALSE), collapse = "\n")
    catalog_info$ready_lock[i] <- grepl(
      "READY_EXACT_26E_CATALOG_FROZEN",
      txt,
      fixed = TRUE
    )
  }
}

setorder(catalog_info, -ready_lock, -mtime)

catalog <- NULL
CATALOG_FILE <- NA_character_

for (p in catalog_info$path) {
  x <- tryCatch(fread(p), error = function(e) NULL)
  if (is.null(x)) next
  if (!all(c("program", "gene") %in% names(x))) next

  x[, program := as.character(program)]
  x[, gene := norm_gene(gene)]
  x <- x[!is.na(program) & nzchar(program) & !is.na(gene) & nzchar(gene)]

  if (uniqueN(x$program) >= 20L) {
    catalog <- x
    CATALOG_FILE <- p
    break
  }
}

if (is.null(catalog)) {
  stop("Frozen 26E catalog(s) found but none had a valid program/gene schema.")
}

fwrite(
  catalog_info,
  file.path(OUT, "P5_RESTART2_26E_CATALOG_SELECTION.tsv"),
  sep = "\t"
)

# =============================================================================
# 4. Program-level Top100 enrichment against all frozen 26E programs
# =============================================================================

universe <- spectrum$gene
N <- length(universe)

top100_use <- intersect(top100_genes, universe)
m <- length(top100_use)

if (m < 95L) {
  stop("Fewer than 95 frozen Top100 genes map to the full spectrum.")
}

programs <- sort(unique(catalog$program))

rows <- vector("list", length(programs))

for (i in seq_along(programs)) {

  pr <- programs[i]
  gs <- intersect(unique(catalog[program == pr, gene]), universe)

  K <- length(gs)
  overlap <- intersect(top100_use, gs)
  k <- length(overlap)

  p <- if (K > 0L) {
    phyper(k - 1, K, N - K, m, lower.tail = FALSE)
  } else {
    NA_real_
  }

  spec_sub <- spectrum[gene %in% gs]

  rows[[i]] <- data.table(
    program = pr,
    universe_N = N,
    frozen_set_K = K,
    top100_m = m,
    overlap_k = k,
    overlap_fraction_top100 = k / m,
    overlap_fraction_geneset = ifelse(K > 0, k / K, NA_real_),
    overlap_genes = paste(
      spectrum[
        gene %in% overlap
      ][order(spectrum_rank), gene],
      collapse = ";"
    ),
    hypergeom_p = p,
    median_P5_rank_of_overlap = ifelse(
      k > 0,
      median(spectrum[gene %in% overlap, spectrum_rank]),
      NA_real_
    ),
    P5_loading_sum_overlap = sum(
      spectrum[gene %in% overlap, pmax(meta_program_score, 0)],
      na.rm = TRUE
    )
  )
}

enrich <- rbindlist(rows)
enrich[, FDR_26E := p.adjust(hypergeom_p, method = "BH")]
setorder(enrich, FDR_26E, -overlap_k)

fwrite(
  enrich,
  file.path(OUT, "P5_RESTART2_26E_TOP100_ENRICHMENT.tsv"),
  sep = "\t"
)

# =============================================================================
# 5. Predefined publication macro-axes
#
# These are groupings of the already-frozen 26E programs, not new gene sets.
# MYELIN/OLIGODENDROCYTE identity is handled separately using the already
# frozen 34A0 content/pathway annotation because the 26E catalog was not
# designed as a complete myelination catalog.
# =============================================================================

AXIS_MAP <- rbindlist(list(

  data.table(
    macro_axis = "LYSOSOME_PROTEOSTASIS_LIPID",
    program = c(
      "Lysosome_autophagy_proteostasis_deconfounded",
      "Phagolysosome_lipid_processing"
    )
  ),

  data.table(
    macro_axis = "METABOLIC_OXIDATIVE_STRESS",
    program = c(
      "Mitochondrial_oxidative_stress_deconfounded"
    )
  ),

  data.table(
    macro_axis = "SENESCENCE_ASSOCIATED_REMODELING",
    program = c(
      "DDR_genome_instability_deconfounded",
      "Chromatin_remodelling_deconfounded",
      "Stable_cell_cycle_exit_deconfounded",
      "SASP_ECM_remodelling_deconfounded",
      "Anti_apoptosis_survival_deconfounded",
      "SAUL_SEN_MAYO",
      "FRIDMAN_SENESCENCE_UP",
      "REACTOME_CELLULAR_SENESCENCE",
      "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP"
    )
  ),

  data.table(
    macro_axis = "COMPETING_REACTIVE_PROGRAMS",
    program = c(
      "Apoptosis_cell_death",
      "IFN_response",
      "Immediate_early_response",
      "JAK_STAT_activation",
      "NFkB_inflammation_chemokine",
      "Proliferation"
    )
  )
))

axis_long <- merge(
  AXIS_MAP,
  enrich,
  by = "program",
  all.x = TRUE,
  sort = FALSE
)

axis_long[, evaluable := !is.na(FDR_26E)]

fwrite(
  axis_long,
  file.path(OUT, "P5_RESTART2_MACRO_AXIS_PROGRAMS.tsv"),
  sep = "\t"
)

axis_summary <- axis_long[
  evaluable == TRUE,
  .(
    n_programs_evaluable = .N,
    n_programs_FDR05 = sum(FDR_26E < 0.05, na.rm = TRUE),
    n_programs_FDR10 = sum(FDR_26E < 0.10, na.rm = TRUE),
    min_FDR = min(FDR_26E, na.rm = TRUE),
    total_top100_overlap_events = sum(overlap_k, na.rm = TRUE),
    supporting_programs_FDR05 = paste(
      program[FDR_26E < 0.05],
      collapse = ";"
    )
  ),
  by = macro_axis
]

fwrite(
  axis_summary,
  file.path(OUT, "P5_RESTART2_MACRO_AXIS_SUMMARY.tsv"),
  sep = "\t"
)

# =============================================================================
# 6. Existing frozen myelin/lineage annotation from 34A0
# =============================================================================

get_field <- function(dt, nm) {
  if (nm %in% names(dt)) as.character(dt[[nm]][1]) else NA_character_
}

myelin_context <- data.table(
  item = c(
    "dominant_general_module",
    "top_general_pathways",
    "reviewed_anchor_pathways",
    "top20_genes",
    "safe_label",
    "interpretation",
    "strong_competing_explanation"
  ),
  value = c(
    get_field(a0, "dominant_general_module"),
    get_field(a0, "top_general_pathways"),
    get_field(a0, "reviewed_anchor_pathways"),
    get_field(a0, "top20_genes"),
    get_field(a0, "safe_label"),
    get_field(a0, "interpretation"),
    get_field(a0, "strong_competing_explanation")
  )
)

fwrite(
  myelin_context,
  file.path(OUT, "P5_RESTART2_EXISTING_MYELIN_LINEAGE_CONTEXT.tsv"),
  sep = "\t"
)

# =============================================================================
# 7. Annotate the Top100 by frozen 26E memberships
# =============================================================================

membership <- catalog[
  gene %in% top100_use,
  .(
    frozen_26E_programs = paste(sort(unique(program)), collapse = ";"),
    n_frozen_26E_programs = uniqueN(program)
  ),
  by = gene
]

top100_annot <- merge(
  spectrum[
    in_frozen_top100 == TRUE,
    .(
      gene,
      meta_program_score,
      spectrum_rank
    )
  ],
  membership,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

setorder(top100_annot, spectrum_rank)
top100_annot[is.na(frozen_26E_programs), frozen_26E_programs := ""]
top100_annot[is.na(n_frozen_26E_programs), n_frozen_26E_programs := 0L]

fwrite(
  top100_annot,
  file.path(OUT, "P5_RESTART2_TOP100_FUNCTIONAL_ANNOTATION.tsv"),
  sep = "\t"
)

# =============================================================================
# 8. Publication-safe functional freeze
# =============================================================================

strong_competing <- FALSE
if ("strong_competing_explanation" %in% names(a0)) {
  z <- tolower(trimws(as.character(a0$strong_competing_explanation[1])))
  strong_competing <- z %in% c("true", "t", "1", "yes")
}

sen_row <- axis_summary[
  macro_axis == "SENESCENCE_ASSOCIATED_REMODELING"
]

lys_row <- axis_summary[
  macro_axis == "LYSOSOME_PROTEOSTASIS_LIPID"
]

ox_row <- axis_summary[
  macro_axis == "METABOLIC_OXIDATIVE_STRESS"
]

sen_sig_n <- if (nrow(sen_row)) sen_row$n_programs_FDR05[1] else 0L
lys_sig_n <- if (nrow(lys_row)) lys_row$n_programs_FDR05[1] else 0L
ox_sig_n <- if (nrow(ox_row)) ox_row$n_programs_FDR05[1] else 0L

functional_class <- if (
  lys_sig_n >= 1L &&
  ox_sig_n >= 1L &&
  sen_sig_n >= 2L
) {
  "MULTI_AXIS_MYELIN_STRESS_REMODELING_WITH_PARTIAL_SENESCENCE_ASSOCIATED_CONTENT"
} else if (
  lys_sig_n >= 1L &&
  ox_sig_n >= 1L
) {
  "MYELIN_LYSOSOME_METABOLIC_STRESS_REMODELING"
} else if (
  lys_sig_n >= 1L
) {
  "MYELIN_LYSOSOME_PROTEOSTASIS_REMODELING"
} else {
  "MIXED_OLIGODENDROGLIAL_PROGRAM_NO_SINGLE_FUNCTIONAL_AXIS_DOMINATES"
}

# Never upgrade to a senescent-cell label when the frozen interpretation already
# records a strong competing explanation.
senescence_wording <- if (strong_competing) {
  "PARTIAL_SENESCENCE_ASSOCIATED_TRANSCRIPTIONAL_FEATURES_ONLY"
} else if (sen_sig_n >= 3L) {
  "SENESCENCE_ASSOCIATED_PROGRAM_ENRICHMENT_WITHOUT_CELLULAR_SENESCENCE_CLAIM"
} else {
  "LIMITED_SENESCENCE_ASSOCIATED_CONTENT"
}

decision <- data.table(
  item = c(
    "P5_foundation_status",
    "formal_TLE_program_candidate",
    "directional_support_only",
    "primary_effect",
    "primary_hedges_g",
    "primary_exact_p",
    "primary_fdr",
    "cross_sampling_correlation",
    "35C_P5_specific_summary_recovered",
    "P5_full_spectrum_genes",
    "P5_top100_genes",
    "lysosome_proteostasis_lipid_FDR05_programs",
    "metabolic_oxidative_FDR05_programs",
    "senescence_associated_FDR05_programs",
    "strong_competing_explanation",
    "FUNCTIONAL_CLASS",
    "SENESCENCE_WORDING",
    "NEXT_STAGE"
  ),
  value = c(
    get_field(a0, "program_use_status"),
    if ("formal_TLE_program_candidate" %in% names(effect))
      as.character(effect$formal_TLE_program_candidate[1]) else NA_character_,
    if ("directional_support_only" %in% names(effect))
      as.character(effect$directional_support_only[1]) else NA_character_,
    if ("primary_effect" %in% names(effect))
      as.character(effect$primary_effect[1]) else NA_character_,
    if ("primary_hedges_g" %in% names(effect))
      as.character(effect$primary_hedges_g[1]) else NA_character_,
    if ("primary_exact_p" %in% names(effect))
      as.character(effect$primary_exact_p[1]) else NA_character_,
    if ("primary_fdr" %in% names(effect))
      as.character(effect$primary_fdr[1]) else NA_character_,
    get_field(a0, "cross_sampling_correlation"),
    TRUE,
    nrow(spectrum),
    length(top100_use),
    lys_sig_n,
    ox_sig_n,
    sen_sig_n,
    strong_competing,
    functional_class,
    senescence_wording,
    "RESTART3_INTEGRATED_REGULATOR_AND_MICROENVIRONMENT_PRIORITIZATION"
  )
)

fwrite(
  decision,
  file.path(OUT, "P5_RESTART2_DECISION.tsv"),
  sep = "\t"
)

# =============================================================================
# 9. Compact report
# =============================================================================

top_enrich <- enrich[
  is.finite(FDR_26E)
][
  order(FDR_26E, -overlap_k)
][
  seq_len(min(.N, 10))
]

report <- c(
  "P5 MATURE OLIGODENDROCYTE PURE-BIOINFORMATICS MANUSCRIPT",
  "RESTART-2 | FUNCTIONAL ARCHITECTURE",
  "========================================================",
  "",
  "FOUNDATION (preserved from RESTART-1)",
  "-------------------------------------",
  paste0(
    "formal_TLE_program_candidate=",
    if ("formal_TLE_program_candidate" %in% names(effect))
      effect$formal_TLE_program_candidate[1] else NA
  ),
  paste0(
    "directional_support_only=",
    if ("directional_support_only" %in% names(effect))
      effect$directional_support_only[1] else NA
  ),
  paste0(
    "primary_effect=",
    if ("primary_effect" %in% names(effect))
      signif(effect$primary_effect[1], 6) else NA
  ),
  paste0(
    "primary_Hedges_g=",
    if ("primary_hedges_g" %in% names(effect))
      signif(effect$primary_hedges_g[1], 6) else NA
  ),
  paste0(
    "primary_exact_p=",
    if ("primary_exact_p" %in% names(effect))
      signif(effect$primary_exact_p[1], 6) else NA
  ),
  paste0(
    "primary_FDR=",
    if ("primary_fdr" %in% names(effect))
      signif(effect$primary_fdr[1], 6) else NA
  ),
  paste0(
    "cross_sampling_correlation=",
    get_field(a0, "cross_sampling_correlation")
  ),
  "",
  "35C RECOVERY",
  "------------",
  paste0(
    "P5_PROGRAM_SPECIFICITY_SUMMARY rows recovered as P5-specific by frozen filename=",
    nrow(p5_spec)
  ),
  "The RESTART-1 NO_P5_ROW_FOUND flag was an extraction artifact, not a biological negative result.",
  "",
  "FUNCTIONAL ARCHITECTURE",
  "-----------------------",
  paste0("P5 full spectrum genes=", nrow(spectrum)),
  paste0("Frozen Top100 mapped=", length(top100_use), "/100"),
  paste0(
    "Lysosome/proteostasis/lipid FDR<0.05 programs=",
    lys_sig_n
  ),
  paste0(
    "Metabolic/oxidative FDR<0.05 programs=",
    ox_sig_n
  ),
  paste0(
    "Senescence-associated FDR<0.05 programs=",
    sen_sig_n
  ),
  paste0("strong_competing_explanation=", strong_competing),
  "",
  paste0("FUNCTIONAL_CLASS=", functional_class),
  paste0("SENESCENCE_WORDING=", senescence_wording),
  "",
  "TOP FROZEN-26E OVERLAPS",
  "-----------------------",
  paste(
    paste0(
      top_enrich$program,
      ": overlap=",
      top_enrich$overlap_k,
      ", FDR=",
      signif(top_enrich$FDR_26E, 4)
    ),
    collapse = "\n"
  ),
  "",
  "INTERPRETATION BOUNDARY",
  "-----------------------",
  "- This is functional annotation of the frozen P5 program, not a new TLE disease test.",
  "- P5 remains a high-priority directional/exploratory disease-linked program unless the frozen upstream disease gate is independently changed.",
  "- Mature-OL/myelin identity is supported by the existing frozen content/pathway annotation; the 26E catalog is used to quantify stress/senescence/reactive content.",
  "- Senescence-associated overlap must not be converted into a senescent-oligodendrocyte label.",
  "",
  "NEXT STAGE",
  "----------",
  "Proceed once to integrated regulator + microenvironment prioritization.",
  "Do not rerun validation, cNMF, QC, or multiple redundant functional-enrichment workflows."
)

writeLines(
  report,
  file.path(OUT, "P5_RESTART2_REPORT.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(OUT, "P5_RESTART2_SESSIONINFO.txt")
)

message("P5 RESTART-2 completed.")
message("Output: ", OUT)

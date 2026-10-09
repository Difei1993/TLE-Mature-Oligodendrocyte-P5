# P5_RESTART6_MANUSCRIPT_EVIDENCE_FREEZE.R
# Final evidence integration for the pure-bioinformatics manuscript centered on
# Oligodendrocytes__P5.
#
# This step DOES NOT run new biological analyses.
# It only integrates already-frozen RESTART1-5 and existing 35C evidence into:
#   1) a manuscript evidence matrix,
#   2) a claim matrix,
#   3) a figure plan,
#   4) a compact manuscript-readiness report.
#
# Scientific boundary:
# - P5 remains directional/high-priority exploratory, not formally FDR-supported.
# - Functional core: myelin–lysosome/proteostasis remodeling.
# - Senescence: partial transcriptional overlap only.
# - TF: no convergent priority regulator.
# - APOE–LRP1: exploratory mechanism candidate only because RESTART5 did not
#   provide consistent spatial/reference support.
# - No further dataset fishing is authorized in this freeze.

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}
suppressPackageStartupMessages(library(data.table))

ROOT <- Sys.getenv("P5_WORK_ROOT")

PROJECT <- file.path(
  ROOT,
  "P5_MATURE_OLIGODENDROCYTE_PURE_BIOINFO_RESTART"
)

R1 <- file.path(PROJECT, "RESTART1_PUBLICATION_EVIDENCE_FREEZE")
R2 <- file.path(PROJECT, "RESTART2_FUNCTIONAL_ARCHITECTURE")
R3 <- file.path(PROJECT, "RESTART3_REGULATOR_MICROENVIRONMENT_PRIORITIZATION")
R4 <- file.path(PROJECT, "RESTART4_TARGETED_DONOR_AWARE_LR")
R5 <- file.path(PROJECT, "RESTART5_GSE295237_SPATIAL_QUALIFICATION")

CROSSREF <- file.path(
  ROOT,
  "35C_P5_CROSS_REFERENCE",
  "09_PROGRAM_SPECIFICITY_AND_PLATFORM_SENSITIVITY"
)

OUT <- file.path(
  PROJECT,
  "RESTART6_MANUSCRIPT_EVIDENCE_FREEZE"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

must_read <- function(path, sep = "\t") {
  if (!file.exists(path)) stop("Missing required file:\n", path)
  fread(path, sep = sep)
}

optional_read <- function(path, sep = NULL) {
  if (!file.exists(path)) return(NULL)
  if (is.null(sep)) return(fread(path))
  fread(path, sep = sep)
}

kv_value <- function(dt, key) {
  if (is.null(dt) || !all(c("item","value") %in% names(dt))) return(NA_character_)
  z <- dt[item == key, value]
  if (length(z) != 1L) return(NA_character_)
  as.character(z)
}

num_value <- function(dt, key) {
  suppressWarnings(as.numeric(kv_value(dt, key)))
}

bool_from_any <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z %in% c("true","t","1","yes")
}

# -------------------------------------------------------------------------
# 1. Read frozen evidence
# -------------------------------------------------------------------------

R1_EFFECT <- must_read(
  file.path(R1, "P5_RESTART1_DISCOVERY_EFFECT_ROW.tsv")
)

R1_A0 <- must_read(
  file.path(R1, "P5_RESTART1_34A0_INTEGRATED_EVIDENCE.tsv")
)

R2_DEC <- must_read(
  file.path(R2, "P5_RESTART2_DECISION.tsv")
)

R2_SPEC <- optional_read(
  file.path(R2, "P5_RESTART2_35C_P5_SPECIFICITY_RECOVERED.tsv"),
  sep = "\t"
)

R3_DEC <- optional_read(
  file.path(R3, "P5_RESTART3_DECISION.tsv"),
  sep = "\t"
)

R3_FIX1_DEC <- optional_read(
  file.path(R3, "P5_RESTART3_FIX1_DECISION.tsv"),
  sep = "\t"
)

R4_DEC <- must_read(
  file.path(R4, "P5_RESTART4_DECISION.tsv")
)

R5_DEC <- must_read(
  file.path(R5, "P5_RESTART5_DECISION.tsv")
)

R5_STATS <- must_read(
  file.path(R5, "P5_RESTART5_SAMPLE_SPATIAL_STATISTICS.tsv")
)

R5_GROUP <- must_read(
  file.path(R5, "P5_RESTART5_GROUP_DIRECTION_SUMMARY.tsv")
)

# Existing 35C files are read directly; no reanalysis.
SPEC_FILE <- file.path(
  CROSSREF,
  "02_results",
  "P5_PROGRAM_SPECIFICITY_SUMMARY.csv"
)

RANK_FILE <- file.path(
  CROSSREF,
  "02_results",
  "MATURE_ADULT_3_VS_3_PROGRAM_SPECIFICITY_RANKING.csv"
)

DONOR_SCORE_FILE <- file.path(
  CROSSREF,
  "01_program_scores",
  "ALL_OLIGODENDROCYTE_PROGRAM_DONOR_SCORES.csv"
)

SPEC_EXIST <- file.exists(SPEC_FILE)
RANK_EXIST <- file.exists(RANK_FILE)
DONOR_SCORE_EXIST <- file.exists(DONOR_SCORE_FILE)

spec35 <- if (SPEC_EXIST) fread(SPEC_FILE) else NULL
rank35 <- if (RANK_EXIST) fread(RANK_FILE) else NULL

if (!is.null(spec35)) {
  fwrite(
    spec35,
    file.path(OUT, "P5_RESTART6_35C_SPECIFICITY_SNAPSHOT.tsv"),
    sep = "\t"
  )
}

if (!is.null(rank35)) {
  fwrite(
    rank35,
    file.path(OUT, "P5_RESTART6_35C_MATURE_RANKING_SNAPSHOT.tsv"),
    sep = "\t"
  )
}

# -------------------------------------------------------------------------
# 2. Recover key frozen values
# -------------------------------------------------------------------------

formal_flag <- if ("formal_TLE_program_candidate" %in% names(R1_EFFECT)) {
  bool_from_any(R1_EFFECT$formal_TLE_program_candidate[1])
} else {
  FALSE
}

directional_flag <- if ("directional_support_only" %in% names(R1_EFFECT)) {
  bool_from_any(R1_EFFECT$directional_support_only[1])
} else {
  NA
}

primary_effect <- if ("primary_effect" %in% names(R1_EFFECT)) {
  as.numeric(R1_EFFECT$primary_effect[1])
} else NA_real_

primary_g <- if ("primary_hedges_g" %in% names(R1_EFFECT)) {
  as.numeric(R1_EFFECT$primary_hedges_g[1])
} else NA_real_

primary_p <- if ("primary_exact_p" %in% names(R1_EFFECT)) {
  as.numeric(R1_EFFECT$primary_exact_p[1])
} else NA_real_

primary_fdr <- if ("primary_fdr" %in% names(R1_EFFECT)) {
  as.numeric(R1_EFFECT$primary_fdr[1])
} else NA_real_

cross_sampling <- if ("cross_sampling_correlation" %in% names(R1_A0)) {
  as.numeric(R1_A0$cross_sampling_correlation[1])
} else NA_real_

functional_class <- kv_value(R2_DEC, "FUNCTIONAL_CLASS")
senescence_wording <- kv_value(R2_DEC, "SENESCENCE_WORDING")

r4_final <- kv_value(R4_DEC, "FINAL_CLASS")
r5_final <- kv_value(R5_DEC, "FINAL_CLASS")

# -------------------------------------------------------------------------
# 3. Evidence matrix
# -------------------------------------------------------------------------

evidence <- data.table(
  evidence_block = c(
    "DISCOVERY_PROGRAM",
    "DISCOVERY_ROBUSTNESS",
    "MATURE_OL_SPECIFICITY_EXTERNAL_CONTEXT",
    "FUNCTIONAL_ARCHITECTURE",
    "SENESCENCE_INTERPRETATION",
    "INTRINSIC_REGULATOR",
    "MICROENVIRONMENT_CONTEXT",
    "TARGETED_LR",
    "SPATIAL_QUALIFICATION",
    "EXTERNAL_HUMAN_ROLE"
  ),

  current_status = c(
    if (formal_flag) "FORMAL" else "DIRECTIONAL_EXPLORATORY",
    "STRONG_INTERNAL_REPRODUCIBILITY",
    if (SPEC_EXIST) "EXISTING_35C_EVIDENCE_AVAILABLE" else "MISSING_35C_SNAPSHOT",
    functional_class,
    senescence_wording,
    "NO_CONVERGENT_PRIORITY_TF",
    if (!is.null(R3_FIX1_DEC)) kv_value(R3_FIX1_DEC, "MICROENVIRONMENT_STATUS") else "EXPLORATORY_CONTEXT",
    r4_final,
    r5_final,
    "GSE140393_TLE_CONTEXT_PLUS_SEAAD_ADULT_MTG_REFERENCE"
  ),

  manuscript_role = c(
    "MAIN",
    "MAIN",
    "MAIN",
    "MAIN",
    "SECONDARY",
    "SUPPLEMENTARY_NEGATIVE",
    "SECONDARY",
    "EXPLORATORY_SECONDARY",
    "SUPPLEMENTARY_NEGATIVE",
    "MAIN_SUPPORT"
  ),

  allowed_claim = c(
    "P5 is a high-priority mature-oligodendrocyte program with directional TLE association; do not call it a formally FDR-significant disease program.",
    "P5 shows large directional effect and strong cross-sampling reproducibility.",
    "Use existing 35C mature-OL specificity/context evidence exactly as originally reported; do not relabel it as definitive replication.",
    "P5 is primarily a myelin–lysosome/proteostasis remodeling program.",
    "Only partial senescence-associated transcriptional overlap is supported; do not call P5 senescent oligodendrocytes.",
    "No master TF/regulator claim.",
    "Astrocyte and myeloid programs covary with P5 at donor level and define microenvironment context only.",
    "APOE–LRP1 and related lipid-input axes are computational candidates; do not call them established mechanisms.",
    "GSE295237 does not provide consistent spatial/reference support for APOE–LRP1–P5; retain this as a negative boundary.",
    "GSE140393 supports independent human TLE context; SEA-AD is adult MTG reference, not TLE replication."
  )
)

fwrite(
  evidence,
  file.path(OUT, "P5_RESTART6_MASTER_EVIDENCE_MATRIX.tsv"),
  sep = "\t"
)

# -------------------------------------------------------------------------
# 4. Manuscript claim matrix
# -------------------------------------------------------------------------

claims <- data.table(
  claim_id = sprintf("C%02d", 1:9),
  proposed_claim = c(
    "A mature oligodendroglial P5 program is directionally increased in TLE.",
    "P5 is reproducible across the two frozen sampling schemes.",
    "P5 is enriched for myelin/lysosome/proteostasis remodeling.",
    "P5 contains partial senescence-associated transcriptional features.",
    "P5 represents senescent oligodendrocytes.",
    "A single master TF drives P5.",
    "Astrocytic and myeloid extracellular programs covary with P5.",
    "APOE–LRP1 is a candidate extracellular lipid-input axis linked to P5.",
    "APOE–LRP1–P5 is spatially validated in independent human epilepsy tissue."
  ),
  decision = c(
    "ALLOW_WITH_QUALIFIER",
    "ALLOW",
    "ALLOW",
    "ALLOW_WITH_QUALIFIER",
    "DO_NOT_CLAIM",
    "DO_NOT_CLAIM",
    "ALLOW_WITH_QUALIFIER",
    "ALLOW_AS_EXPLORATORY",
    "DO_NOT_CLAIM"
  ),
  rationale = c(
    paste0(
      "effect=", signif(primary_effect, 4),
      "; Hedges_g=", signif(primary_g, 4),
      "; exact_P=", signif(primary_p, 4),
      "; FDR=", signif(primary_fdr, 4),
      "; formal_TLE_program_candidate=", formal_flag
    ),
    paste0("cross_sampling_correlation=", signif(cross_sampling, 4)),
    paste0("FUNCTIONAL_CLASS=", functional_class),
    paste0("SENESCENCE_WORDING=", senescence_wording),
    "Frozen upstream explicitly prohibits a senescent-cell label.",
    "RESTART3 recovered no convergent priority TF.",
    "RESTART3 FIX1 recovered donor-level cross-lineage context; correlations are not LR/causal evidence.",
    paste0(
      "RESTART4=", r4_final,
      "; but spatial qualification must remain negative."
    ),
    paste0(
      "RESTART5=", r5_final,
      "; receiver-positive samples=",
      kv_value(R5_DEC, "receiver_positive_samples"),
      "/4; neighbor-APOE-positive samples=",
      kv_value(R5_DEC, "neighbor_APOE_positive_samples"),
      "/4."
    )
  )
)

fwrite(
  claims,
  file.path(OUT, "P5_RESTART6_CLAIM_MATRIX.tsv"),
  sep = "\t"
)

# -------------------------------------------------------------------------
# 5. Figure plan
# -------------------------------------------------------------------------

figplan <- data.table(
  figure = c(
    "Figure 1",
    "Figure 2",
    "Figure 3",
    "Figure 4",
    "Supplementary Figure S1",
    "Supplementary Figure S2"
  ),
  purpose = c(
    "Discovery and donor-aware P5 association",
    "Independent mature-OL specificity / human context",
    "Functional architecture of P5",
    "Microenvironment and exploratory LR candidates",
    "Regulator-negative / robustness evidence",
    "GSE295237 spatial qualification negative boundary"
  ),
  core_panels = c(
    "P5 program definition; donor-level primary/sensitivity effects; effect size; cross-sampling concordance",
    "35C mature-OL program specificity; GSE140393 TLE mature-OL context; SEA-AD adult MTG reference role",
    "P5 Top genes; myelin/lysosome/proteostasis; lipid-processing; limited senescence overlap",
    "Astrocyte/myeloid donor-program coupling; targeted LR shortlist; APOE–LRP1 shown as exploratory candidate",
    "TF screen summary; LODO/alternative sampling robustness",
    "Four GSE295237 spatial maps plus sample-level rho summary; clearly labeled FCDI temporal-neocortex reference, not TLE"
  ),
  main_or_supplement = c(
    "MAIN",
    "MAIN",
    "MAIN",
    "MAIN",
    "SUPPLEMENT",
    "SUPPLEMENT"
  )
)

fwrite(
  figplan,
  file.path(OUT, "P5_RESTART6_FIGURE_PLAN.tsv"),
  sep = "\t"
)

# -------------------------------------------------------------------------
# 6. Readiness decision
# -------------------------------------------------------------------------

spatial_negative <- identical(
  r5_final,
  "NO_CONSISTENT_SPATIAL_REFERENCE_SUPPORT"
)

external_ready <- SPEC_EXIST && RANK_EXIST && DONOR_SCORE_EXIST

FINAL_MANUSCRIPT_STATUS <- if (
  directional_flag &&
  is.finite(cross_sampling) &&
  cross_sampling > 0.8 &&
  identical(
    functional_class,
    "MYELIN_LYSOSOME_PROTEOSTASIS_REMODELING"
  ) &&
  external_ready
) {
  "READY_FOR_FIGURE_ASSEMBLY_AND_MANUSCRIPT_DRAFTING_WITH_EXPLORATORY_MECHANISM"
} else {
  "MANUSCRIPT_POSSIBLE_BUT_EXTERNAL_EVIDENCE_SNAPSHOT_REQUIRES_REVIEW"
}

decision <- data.table(
  item = c(
    "formal_TLE_program_candidate",
    "directional_support_only",
    "primary_effect",
    "primary_Hedges_g",
    "primary_exact_p",
    "primary_FDR",
    "cross_sampling_correlation",
    "functional_class",
    "senescence_wording",
    "35C_specificity_file_present",
    "35C_ranking_file_present",
    "35C_donor_score_file_present",
    "RESTART4_LR_status",
    "RESTART5_spatial_status",
    "FINAL_MANUSCRIPT_STATUS",
    "NEXT_ACTION"
  ),
  value = c(
    formal_flag,
    directional_flag,
    primary_effect,
    primary_g,
    primary_p,
    primary_fdr,
    cross_sampling,
    functional_class,
    senescence_wording,
    SPEC_EXIST,
    RANK_EXIST,
    DONOR_SCORE_EXIST,
    r4_final,
    r5_final,
    FINAL_MANUSCRIPT_STATUS,
    "FIGURE_ASSEMBLY_AND_RESULTS_METHODS_DRAFTING_NO_NEW_DISCOVERY_ANALYSIS"
  )
)

fwrite(
  decision,
  file.path(OUT, "P5_RESTART6_DECISION.tsv"),
  sep = "\t"
)

# -------------------------------------------------------------------------
# 7. Compact report
# -------------------------------------------------------------------------

writeLines(
  c(
    "P5 MATURE OLIGODENDROCYTE PURE-BIOINFORMATICS MANUSCRIPT",
    "RESTART-6 | FINAL MANUSCRIPT EVIDENCE FREEZE",
    "=================================================",
    "",
    paste0("FINAL_MANUSCRIPT_STATUS=", FINAL_MANUSCRIPT_STATUS),
    "",
    "CORE STORY",
    "----------",
    "- P5 is a high-priority mature-oligodendrocyte program with directional TLE association and strong internal reproducibility.",
    "- The dominant biological interpretation is myelin–lysosome/proteostasis remodeling.",
    "- Senescence is secondary/partial transcriptional overlap only.",
    "- No convergent master TF is supported.",
    "- Astrocyte/myeloid microenvironment coupling and APOE–LRP1 are exploratory candidate mechanisms.",
    "- GSE295237 did NOT consistently spatially qualify APOE–LRP1–P5; this negative result is preserved and prevents a mechanistic claim.",
    "",
    "KEY NUMBERS",
    "-----------",
    paste0("primary_effect=", signif(primary_effect, 6)),
    paste0("primary_Hedges_g=", signif(primary_g, 6)),
    paste0("primary_exact_p=", signif(primary_p, 6)),
    paste0("primary_FDR=", signif(primary_fdr, 6)),
    paste0("cross_sampling_correlation=", signif(cross_sampling, 6)),
    paste0("RESTART5_final=", r5_final),
    "",
    "EDITORIAL POSITIONING",
    "---------------------",
    "- Main paper should be a disease-associated/pathological mature-OL program paper, not a causal APOE–LRP1 mechanism paper.",
    "- APOE–LRP1 may appear as a hypothesis-generating microenvironmental axis in a main or secondary figure, with explicit caveats.",
    "- The negative spatial result belongs in Supplementary material or a limitation panel; it should not be hidden.",
    "- Do not add new external datasets, new TF algorithms, new LR algorithms, or alternative spatial thresholds to rescue the mechanism.",
    "",
    "NEXT ACTION",
    "-----------",
    "Assemble Figures 1-4 and Supplementary Figures S1-S2, then draft Results and Methods from the frozen claim matrix."
  ),
  file.path(OUT, "P5_RESTART6_REPORT.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(OUT, "P5_RESTART6_SESSIONINFO.txt")
)

message("P5 RESTART-6 completed.")
message("Output: ", OUT)

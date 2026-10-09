###############################################################################
# 32B-0
# Prepare donor-balanced raw-count inputs for unbiased cNMF discovery
#
# No senescence signatures are used.
# No TLE/Control label is used for feature or cell-state discovery.
###############################################################################

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

set.seed(20260728)

PROJECT_ROOT <- normalizePath(
  Sys.getenv("P5_WORK_ROOT"),
  winslash = "/",
  mustWork = TRUE
)

OBJECT_PATH <- normalizePath(
  file.path(
    PROJECT_ROOT,
    "02_objects",
    "human_harmony_r0p4_level1_annotated_with_clinical_metadata_repaired.rds"
  ),
  winslash = "/",
  mustWork = TRUE
)

STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

OUT_ROOT <- file.path(
  PROJECT_ROOT,
  "99_frozen_releases",
  paste0("32B0_BALANCED_CNMf_INPUTS_", STAMP)
)

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "cnmf_input"), showWarnings = FALSE)

###############################################################################
# 1. Packages
###############################################################################

packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "dplyr",
  "readr",
  "tibble",
  "jsonlite",
  "digest"
)

missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
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
  library(readr)
  library(tibble)
})

###############################################################################
# 2. Load and audit object
###############################################################################

human <- readRDS(OBJECT_PATH)

if (!inherits(human, "Seurat")) {
  stop("输入不是Seurat对象。")
}

required_metadata <- c(
  "donor_id",
  "analysis_condition",
  "annotation_level1_manual"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(human[[]])
)

if (length(missing_metadata) > 0) {
  stop(
    "缺少必要元数据：",
    paste(missing_metadata, collapse = ", ")
  )
}

meta <- as.data.frame(
  human[[]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(meta) <- colnames(human)

meta$donor_id <- as.character(meta$donor_id)
meta$condition <- as.character(meta$analysis_condition)
meta$cell_type <- as.character(meta$annotation_level1_manual)

if (!setequal(unique(meta$condition), c("Normal", "TLE"))) {
  stop("疾病分组不是预期的Normal/TLE。")
}

donor_map <- meta %>%
  dplyr::distinct(donor_id, condition)

if (
  nrow(donor_map) != 8L ||
  sum(donor_map$condition == "Normal") != 4L ||
  sum(donor_map$condition == "TLE") != 4L
) {
  stop("供者结构不是4 Normal + 4 TLE。")
}

###############################################################################
# 3. Lineage plan
###############################################################################

lineage_plan <- tibble::tribble(
  ~cell_type,               ~safe_name,               ~k_min, ~k_max,
  "Astrocytes",             "Astrocytes",                  4L,    12L,
  "Microglia/macrophage",   "Microglia_macrophage",        4L,    12L,
  "OPC",                    "OPC",                         3L,    10L,
  "Oligodendrocytes",       "Oligodendrocytes",             4L,    14L
)

readr::write_csv(
  lineage_plan,
  file.path(OUT_ROOT, "32B0_lineage_plan.csv")
)

###############################################################################
# 4. Extract one joined raw-count matrix
###############################################################################

extract_lineage_counts <- function(
    object,
    lineage_cells,
    lineage_name) {
  
  message("")
  message("提取谱系：", lineage_name)
  
  sub_object <- object[, lineage_cells]
  
  if (!"RNA" %in% Seurat::Assays(sub_object)) {
    stop(lineage_name, "：对象中不存在RNA assay。")
  }
  
  SeuratObject::DefaultAssay(sub_object) <- "RNA"
  
  if (inherits(sub_object[["RNA"]], "Assay5")) {
    
    layer_names <- SeuratObject::Layers(sub_object[["RNA"]])
    
    count_layers <- grep(
      "^counts($|\\.)",
      layer_names,
      value = TRUE
    )
    
    if (length(count_layers) == 0L) {
      stop(lineage_name, "：没有counts layer。")
    }
    
    if (length(count_layers) > 1L) {
      message(
        lineage_name,
        "：合并",
        length(count_layers),
        "个counts layers"
      )
      
      sub_object <- SeuratObject::JoinLayers(
        object = sub_object,
        assay = "RNA",
        layers = count_layers,
        new = "counts"
      )
    }
    
    counts <- SeuratObject::LayerData(
      object = sub_object,
      assay = "RNA",
      layer = "counts"
    )
    
  } else {
    
    counts <- SeuratObject::GetAssayData(
      object = sub_object,
      assay = "RNA",
      layer = "counts"
    )
  }
  
  if (
    nrow(counts) == 0L ||
    ncol(counts) == 0L
  ) {
    stop(lineage_name, "：原始计数矩阵为空。")
  }
  
  if (any(counts@x < 0)) {
    stop(lineage_name, "：原始计数中出现负值。")
  }
  
  if (any(abs(counts@x - round(counts@x)) > 1e-8)) {
    stop(lineage_name, "：counts不是整数型原始计数。")
  }
  
  counts
}

###############################################################################
# 5. Write 10x-style MTX input
###############################################################################

write_10x_input <- function(
    counts,
    cell_metadata,
    output_directory) {
  
  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if (!identical(
    colnames(counts),
    rownames(cell_metadata)
  )) {
    stop("写出前counts与metadata细胞顺序不一致。")
  }
  
  Matrix::writeMM(
    counts,
    file.path(output_directory, "counts.mtx")
  )
  
  write.table(
    colnames(counts),
    file.path(output_directory, "barcodes.tsv"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )
  
  gene_table <- data.frame(
    gene_id = rownames(counts),
    gene_name = rownames(counts),
    stringsAsFactors = FALSE
  )
  
  write.table(
    gene_table,
    file.path(output_directory, "genes.tsv"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )
  
  readr::write_csv(
    tibble::rownames_to_column(
      cell_metadata,
      var = "cell_id"
    ),
    file.path(output_directory, "cell_metadata.csv")
  )
}

###############################################################################
# 6. Balanced sampling and input export
###############################################################################

technical_pattern <- paste0(
  "^MT-|",
  "^RPL[0-9]|",
  "^RPS[0-9]|",
  "^HBA[0-9]|",
  "^HBB$|",
  "^MALAT1$"
)

sex_linked_genes <- c(
  "XIST",
  "RPS4Y1",
  "KDM5D",
  "DDX3Y",
  "EIF1AY",
  "UTY",
  "ZFY",
  "USP9Y",
  "TMSB4Y",
  "NLGN4Y"
)

seed_plan <- c(
  primary = 20260728L,
  sensitivity = 20260729L
)

summary_list <- list()
command_lines <- character()

for (i in seq_len(nrow(lineage_plan))) {
  
  lineage <- lineage_plan$cell_type[i]
  lineage_safe <- lineage_plan$safe_name[i]
  
  lineage_cells <- rownames(meta)[
    meta$cell_type == lineage
  ]
  
  lineage_meta <- meta[
    lineage_cells,
    c("donor_id", "condition", "cell_type"),
    drop = FALSE
  ]
  
  donor_counts <- lineage_meta %>%
    dplyr::count(
      donor_id,
      condition,
      name = "available_cells"
    ) %>%
    arrange(condition, donor_id)
  
  if (nrow(donor_counts) != 8L) {
    stop(
      lineage,
      "：不是所有8名供者均有该谱系细胞。"
    )
  }
  
  minimum_donor_cells <- min(
    donor_counts$available_cells
  )
  
  # Use at most 85% of the smallest donor to preserve an unused reserve.
  donor_cap <- floor(
    0.85 * minimum_donor_cells / 50
  ) * 50
  
  donor_cap <- min(
    donor_cap,
    1000L
  )
  
  if (donor_cap < 150L) {
    stop(
      lineage,
      "：最小供者细胞数不足以进行正式平衡cNMF。"
    )
  }
  
  counts_full <- extract_lineage_counts(
    object = human,
    lineage_cells = lineage_cells,
    lineage_name = lineage
  )
  
  common_cells <- intersect(
    colnames(counts_full),
    rownames(lineage_meta)
  )
  
  counts_full <- counts_full[
    ,
    common_cells,
    drop = FALSE
  ]
  
  lineage_meta <- lineage_meta[
    common_cells,
    ,
    drop = FALSE
  ]
  
  if (!identical(
    colnames(counts_full),
    rownames(lineage_meta)
  )) {
    stop(lineage, "：counts与元数据未对齐。")
  }
  
  minimum_detected_cells <- max(
    10L,
    ceiling(0.005 * ncol(counts_full))
  )
  
  detected_cells <- Matrix::rowSums(
    counts_full > 0
  )
  
  eligible_genes <- rownames(counts_full)[
    detected_cells >= minimum_detected_cells
  ]
  
  eligible_genes <- eligible_genes[
    !grepl(
      technical_pattern,
      eligible_genes,
      ignore.case = TRUE
    )
  ]
  
  eligible_genes <- setdiff(
    eligible_genes,
    sex_linked_genes
  )
  
  if (length(eligible_genes) < 2000L) {
    stop(
      lineage,
      "：过滤后可用基因少于2000。"
    )
  }
  
  counts_full <- counts_full[
    eligible_genes,
    ,
    drop = FALSE
  ]
  
  lineage_root <- file.path(
    OUT_ROOT,
    "cnmf_input",
    lineage_safe
  )
  
  dir.create(
    lineage_root,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  readr::write_csv(
    donor_counts,
    file.path(
      lineage_root,
      "donor_cell_availability.csv"
    )
  )
  
  writeLines(
    eligible_genes,
    file.path(
      lineage_root,
      "eligible_genes_after_technical_filter.txt"
    )
  )
  
  for (seed_name in names(seed_plan)) {
    
    sampling_seed <- seed_plan[[seed_name]]
    set.seed(sampling_seed)
    
    sampled_cells <- unlist(
      lapply(
        split(
          rownames(lineage_meta),
          lineage_meta$donor_id
        ),
        function(cells_one_donor) {
          sample(
            cells_one_donor,
            size = donor_cap,
            replace = FALSE
          )
        }
      ),
      use.names = FALSE
    )
    
    sampled_cells <- sample(
      sampled_cells,
      length(sampled_cells),
      replace = FALSE
    )
    
    sampled_counts <- counts_full[
      ,
      sampled_cells,
      drop = FALSE
    ]
    
    sampled_meta <- lineage_meta[
      sampled_cells,
      ,
      drop = FALSE
    ]
    
    input_directory <- file.path(
      lineage_root,
      paste0("balanced_", seed_name)
    )
    
    write_10x_input(
      counts = sampled_counts,
      cell_metadata = sampled_meta,
      output_directory = input_directory
    )
    
    readr::write_csv(
      sampled_meta %>%
        dplyr::count(
          donor_id,
          condition,
          name = "sampled_cells"
        ),
      file.path(
        input_directory,
        "sampling_audit.csv"
      )
    )
    
    rm(sampled_counts, sampled_meta)
    invisible(gc())
  }
  
  k_values <- seq(
    lineage_plan$k_min[i],
    lineage_plan$k_max[i]
  )
  
  summary_list[[lineage_safe]] <- tibble(
    cell_type = lineage,
    safe_name = lineage_safe,
    total_available_cells = ncol(counts_full),
    minimum_donor_cells = minimum_donor_cells,
    balanced_cells_per_donor = donor_cap,
    balanced_total_cells = donor_cap * 8L,
    eligible_genes = length(eligible_genes),
    k_values = paste(k_values, collapse = " "),
    pilot_iterations = 30L,
    final_iterations = 100L
  )
  
  input_mtx <- normalizePath(
    file.path(
      lineage_root,
      "balanced_primary",
      "counts.mtx"
    ),
    winslash = "/",
    mustWork = TRUE
  )
  
  run_root <- file.path(
    PROJECT_ROOT,
    "04_tables",
    paste0("32B_CNMf_WORK_", STAMP),
    lineage_safe
  )
  
  command_lines <- c(
    command_lines,
    paste0("# ", lineage),
    paste(
      "cnmf prepare",
      "--output-dir",
      shQuote(gsub("\\\\", "/", run_root)),
      "--name",
      paste0(lineage_safe, "_primary"),
      "-c",
      shQuote(input_mtx),
      "-k",
      paste(k_values, collapse = " "),
      "--n-iter 30",
      "--seed 20260728"
    ),
    paste(
      "cnmf factorize",
      "--output-dir",
      shQuote(gsub("\\\\", "/", run_root)),
      "--name",
      paste0(lineage_safe, "_primary"),
      "--worker-index 0",
      "--total-workers 1"
    ),
    paste(
      "cnmf combine",
      "--output-dir",
      shQuote(gsub("\\\\", "/", run_root)),
      "--name",
      paste0(lineage_safe, "_primary")
    ),
    paste(
      "cnmf k_selection_plot",
      "--output-dir",
      shQuote(gsub("\\\\", "/", run_root)),
      "--name",
      paste0(lineage_safe, "_primary")
    ),
    ""
  )
  
  rm(counts_full)
  invisible(gc())
}

cnmf_plan <- dplyr::bind_rows(summary_list)

readr::write_csv(
  cnmf_plan,
  file.path(OUT_ROOT, "32B0_CNMf_PLAN.csv")
)

writeLines(
  command_lines,
  file.path(OUT_ROOT, "32B1_PILOT_COMMANDS.txt")
)

###############################################################################
# 7. Report and manifest
###############################################################################

report <- c(
  "# 32B-0 Balanced cNMF Input Preparation",
  "",
  "## Status: READY_FOR_CNMf_PILOT",
  "",
  "- No senescence gene set was used.",
  "- No disease label was used for gene or program discovery.",
  "- Cells were sampled equally from all eight donors.",
  "- Two independent balanced cell samples were generated.",
  "- Mitochondrial, ribosomal, hemoglobin, MALAT1 and major sex-linked genes were excluded from discovery.",
  "- Raw integer RNA counts were retained.",
  "",
  "## Next stage",
  "",
  "Run the primary balanced pilot for K selection. Do not run consensus until K has been reviewed."
)

writeLines(
  report,
  file.path(OUT_ROOT, "32B0_EXECUTIVE_SUMMARY.md")
)

output_files <- list.files(
  OUT_ROOT,
  recursive = TRUE,
  full.names = TRUE
)

output_files <- output_files[
  !file.info(output_files)$isdir
]

manifest <- tibble(
  relative_path = substring(
    output_files,
    nchar(OUT_ROOT) + 2
  ),
  size_bytes = file.info(output_files)$size,
  sha256 = vapply(
    output_files,
    function(path) {
      digest::digest(
        file = path,
        algo = "sha256",
        serialize = FALSE
      )
    },
    character(1)
  )
)

readr::write_csv(
  manifest,
  file.path(OUT_ROOT, "32B0_FILE_MANIFEST.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(OUT_ROOT, "32B0_sessionInfo.txt")
)

message("")
message("==========================================================")
message("32B-0运行完成")
message("==========================================================")
message("状态：READY_FOR_CNMf_PILOT")
message("")
print(cnmf_plan)
message("")
message("结果目录：")
message(OUT_ROOT)
message("")
message("下一步命令：")
message(file.path(OUT_ROOT, "32B1_PILOT_COMMANDS.txt"))
message("==========================================================")
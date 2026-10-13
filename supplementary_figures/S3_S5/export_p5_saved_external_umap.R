# Read-only coordinate/label conversion. No graphics device or analysis functions.
options(stringsAsFactors=FALSE)
# SPDX-License-Identifier: MIT
# CLI/configuration only; stored UMAP and metadata extraction below is unchanged.
cli <- commandArgs(trailingOnly=TRUE)
usage <- paste(
    "Usage: Rscript export_p5_saved_external_umap.R --objects-dir PATH --output-dir PATH [--check-inputs-only]",
    "Reads two existing frozen RDS objects without calculating UMAP, clustering, or expression statistics.",
    "Coordinate/cell-key exports are private local inputs and are not licensed for redistribution by this code.",
    sep="\n"
)
if ("--help" %in% cli || "-h" %in% cli) {cat(usage,"\n");quit(status=0)}
read_flag <- function(flag) {
    where <- which(cli==flag)
    if (length(where)!=1L || where>=length(cli) || startsWith(cli[where+1L],"--")) stop(usage)
    cli[where+1L]
}
root <- normalizePath(path.expand(read_flag("--objects-dir")),winslash="/",mustWork=FALSE)
output <- normalizePath(path.expand(read_flag("--output-dir")),winslash="/",mustWork=FALSE)
check_only <- "--check-inputs-only" %in% cli
required <- c('GSE140393_ALL_CANDIDATE_CELLS_MANUALLY_ANNOTATED.rds','GSE140393_BONA_FIDE_OLIGODENDROGLIA.rds')
missing <- required[!file.exists(file.path(root,required))]
if (check_only) {
    cat("AVAILABILITY_CHECK_ONLY: ",length(required)-length(missing),"/2 objects; data not read, exports not written.\n",sep="")
    if (length(missing)) cat("Missing: ",paste(missing,collapse=", "),"\n",sep="")
    quit(status=if(length(missing)) 2 else 0)
}
if (length(missing)) stop("Required archived RDS inputs are unavailable: ",paste(missing,collapse=", "))
if (output==root || startsWith(paste0(root,"/"),paste0(output,"/"))) stop("Output must not equal or contain the input directory")
if (file.exists(output) && (!dir.exists(output) || length(list.files(output,all.files=TRUE,no..=TRUE))>0)) stop("Output must be a new or empty local directory")
dir.create(output,recursive=TRUE,showWarnings=FALSE)

jobs <- c(all_candidates='GSE140393_ALL_CANDIDATE_CELLS_MANUALLY_ANNOTATED.rds',bona_fide='GSE140393_BONA_FIDE_OLIGODENDROGLIA.rds')
for (name in names(jobs)) {
    obj <- readRDS(file.path(root,jobs[[name]]))
    metadata <- methods::slot(obj,'meta.data')
    reductions <- methods::slot(obj,'reductions')
    coordinates <- methods::slot(reductions[['umap']],'cell.embeddings')
    stopifnot(identical(rownames(metadata),rownames(coordinates)))
    label <- as.character(metadata[['manual_label_35c']])
    stopifnot(!anyNA(label),ncol(coordinates)==2L,all(is.finite(coordinates)))
    result <- data.frame(cell_key=rownames(metadata),umap_1=coordinates[,1],umap_2=coordinates[,2],manual_label=label,check.names=FALSE)
    write.csv(result,file.path(output,paste0('GSE140393_',name,'_saved_umap.csv')),row.names=FALSE)
    writeBin(as.numeric(coordinates),file.path(output,paste0('GSE140393_',name,'_coordinates_double.bin')),size=8,endian='little')
    levels <- sort(unique(label))
    colors <- grDevices::hcl(h=seq(15,375,length.out=length(levels)+1L)[seq_along(levels)],c=100,l=65)
    write.csv(data.frame(manual_label=levels,color=colors),file.path(output,paste0('GSE140393_',name,'_original_discrete_colors.csv')),row.names=FALSE)
    cat(name,':',nrow(result),'stored coordinates; labels:',length(levels),'\n')
    rm(obj,metadata,reductions,coordinates,result);invisible(gc())
}
writeLines(paste('R conversion only:',getRversion()),file.path(output,'coordinate_conversion_runtime.txt'))

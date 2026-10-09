suppressPackageStartupMessages(library(xml2))
root <- paste0(Sys.getenv("P5_WORK_ROOT"), "/P5_MATURE_OLIGODENDROCYTE_PURE_BIOINFO_RESTART/MOLECULAR_NEUROBIOLOGY_SUBMISSION/FINAL_SUBMISSION_FILES/FIGURES")
targets <- list(
 Figure1='Discovery-family FDR = 0.60',
 FigureS1=c('615 evaluated regulator associations;','no regulator met the priority criterion'),
 FigureS2='GSE295237: FCD type I and controls; FCD I is not TLE.'
)
for(name in names(targets)){
 path <- file.path(root,paste0(name,'.svg'))
 doc <- read_xml(path)
 nodes <- xml_find_all(doc,'.//*[local-name()="text"]')
 for(target in targets[[name]]){
  ix <- which(xml_text(nodes)==target);stopifnot(length(ix)==1)
  node <- nodes[[ix]]
  stopifnot(xml_attr(node,'textLength')=='NA',xml_attr(node,'lengthAdjust')=='NA')
  xml_set_attr(node,'textLength',NULL);xml_set_attr(node,'lengthAdjust',NULL)
 }
 write_xml(doc,path)
 cat(name,': annotation glyph-length attributes removed.\n',sep='')
}

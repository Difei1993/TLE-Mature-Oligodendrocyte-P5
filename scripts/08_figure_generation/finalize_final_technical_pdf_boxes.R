suppressPackageStartupMessages(library(jsonlite))
p <- paste0(Sys.getenv("P5_WORK_ROOT"), "/P5_MATURE_OLIGODENDROCYTE_PURE_BIOINFO_RESTART/MOLECULAR_NEUROBIOLOGY_SUBMISSION/FINAL_SUBMISSION_FILES")
metadata<-fromJSON(file.path(p,'QC/WORK/pdf_page_box_metadata.json'),simplifyVector=FALSE)
dimensions<-fromJSON(file.path(p,'QC/WORK/figure_export_log.json'),simplifyVector=FALSE)
for(m in metadata){
 name<-sub('\\.pdf$','',basename(m$path));h<-dimensions[[name]]$height_mm/25.4*72;w<-174/25.4*72
 dict<-sub('/MediaBox \\[.*?\\]',sprintf('/MediaBox [0 0 %.12f %.12f]',w,h),m$dictionary,perl=TRUE)
 stopifnot(dict!=m$dictionary)
 bytes<-readBin(m$path,'raw',n=file.info(m$path)$size);offset<-length(bytes)
 object<-sprintf('%d 0 obj\n%s\nendobj\n',m$object_id,dict)
 xrefpos<-offset+nchar(object,type='bytes')
 append<-paste0(object,sprintf('xref\n%d 1\n%010d 00000 n \ntrailer\n<< /Size %d /Root %d 0 R /Info %d 0 R /Prev %d >>\nstartxref\n%d\n%%%%EOF\n',m$object_id,offset,m$size,m$root_id,m$info_id,m$prev_xref,xrefpos))
 con<-file(m$path,'wb');writeBin(bytes,con);writeBin(charToRaw(append),con);close(con)
 cat(name,': 174 mm page box; content streams unchanged.\n',sep='')
}

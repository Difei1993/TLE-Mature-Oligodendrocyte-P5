options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(xml2);library(grid);library(ragg);library(png);library(base64enc);library(jsonlite);library(svglite)})
p <- paste0(Sys.getenv("P5_WORK_ROOT"), "/P5_MATURE_OLIGODENDROCYTE_PURE_BIOINFO_RESTART")
base <- file.path(p,'MOLECULAR_NEUROBIOLOGY_SUBMISSION')
out <- file.path(base,'FINAL_SUBMISSION_FILES');figdir<-file.path(out,'FIGURES');work<-file.path(out,'QC/WORK')
width_mm<-174; dpi<-600; family<-'Arial'
# Contract: preserve all frozen quantitative grids and the spatial composite.
# Only three specified annotation replacements require new rendering. Figure
# 2/3/4 are copied byte-for-byte because their local PNGs already meet 600 dpi.
# Reuse only the R primitive renderer functions, never the previous source
# selection or export loop. Every input SVG comes from the current submission.
helper<-parse(file.path(base,'WORK/format_submission_figures.R'))
for(expr in helper){
 if(is.call(expr) && as.character(expr[[1]])=='<-' && as.character(expr[[2]])%in%c('attrnum','style_list','svg_grob'))eval(expr)
}
replacements<-list(
 Figure1=c('Primary FDR = 0.60'='Discovery-family FDR = 0.60'),
 FigureS1=c('615 evaluated regulator-associated rows;'='615 evaluated regulator associations;','no priority regulator'='no regulator met the priority criterion'),
 FigureS2=c('GSE295237: FCD type I and controls; FCDI is not TLE.'='GSE295237: FCD type I and controls; FCD I is not TLE.')
)
log<-list()
for(name in names(replacements)){
 doc<-read_xml(file.path(base,'FIGURES',paste0(name,'.svg')))
 texts<-xml_find_all(doc,'.//*[local-name()="text"]');txt<-xml_text(texts)
 for(old in names(replacements[[name]])){
  pos<-which(txt==old);stopifnot(length(pos)==1)
  xml_text(texts[[pos]])<-unname(replacements[[name]][old])
  xml_set_attr(texts[[pos]],'textLength',NULL);xml_set_attr(texts[[pos]],'lengthAdjust',NULL)
 }
 write_xml(doc,file.path(figdir,paste0(name,'.svg')))
 vb<-scan(text=xml_attr(xml_root(doc),'viewBox'),quiet=TRUE);W<-vb[3];H<-vb[4]
 info<-list(doc=doc,name=name,W=W,H=H,scale=174/25.4*72/W,height_mm=174*H/W)
 g<-svg_grob(info)
 svglite::svglite(file.path(work,paste0(name,'_render_check.svg')),width=width_mm/25.4,height=info$height_mm/25.4,bg='white')
 grid.newpage();grid.draw(g);dev.off()
 grDevices::cairo_pdf(file.path(figdir,paste0(name,'.pdf')),width=width_mm/25.4+1/72,height=info$height_mm/25.4+1/72,family=family,bg='white')
 grid.newpage();grid.draw(g);dev.off()
 pw<-round(width_mm/25.4*dpi);ph<-round(pw*H/W)
 ragg::agg_png(file.path(figdir,paste0(name,'.png')),width=pw,height=ph,units='px',res=dpi,background='white')
 grid.newpage();grid.draw(g);dev.off()
 ragg::agg_png(file.path(work,paste0(name,'_preview.png')),width=1200,height=round(1200*H/W),units='px',res=600*1200/pw,background='white')
 grid.newpage();grid.draw(g);dev.off()
 log[[name]]<-info[c('name','W','H','scale','height_mm')]
 cat(name,': only specified annotation text changed; all geometry preserved.\n',sep='')
}
write_json(log,file.path(work,'figure_export_log.json'),pretty=TRUE,auto_unbox=TRUE,digits=12)

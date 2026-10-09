options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(xml2);library(grid);library(ragg);library(png);library(base64enc);library(jsonlite);library(svglite)})
p<-paste0(Sys.getenv("P5_WORK_ROOT"), "/P5_MATURE_OLIGODENDROCYTE_PURE_BIOINFO_RESTART")
f<-file.path(p,'MOLECULAR_NEUROBIOLOGY_SUBMISSION/FINAL_SUBMISSION_FILES')
o<-file.path(p,'MOLECULAR_NEUROBIOLOGY_FINAL_REVIEW_PACKAGE')
w<-file.path(p,'MOLECULAR_NEUROBIOLOGY_REVIEW_BUILD_WORK')
dir.create(file.path(o,'02_MAIN_FIGURES'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(o,'04_SOURCE_PREVIEWS'),recursive=TRUE,showWarnings=FALSE)
width_mm<-174
dpi<-600
family<-'Arial'
# Exact reuse of frozen vector primitives: no input data or statistics are recomputed.
helper<-parse(file.path(p,'MOLECULAR_NEUROBIOLOGY_SUBMISSION/WORK/format_submission_figures.R'))
for(expr in helper){if(is.call(expr)&&as.character(expr[[1]])=='<-'&&as.character(expr[[2]])%in%c('attrnum','style_list','svg_grob'))eval(expr)}
log<-list()
for(i in 1:4){
 name<-paste0('Figure',i);prefix<-file.path(o,'02_MAIN_FIGURES',paste0('Fig',i))
 doc<-read_xml(file.path(f,'FIGURES',paste0(name,'.svg')))
 vb<-scan(text=xml_attr(xml_root(doc),'viewBox'),quiet=TRUE);W<-vb[3];H<-vb[4]
 info<-list(doc=doc,name=name,W=W,H=H,scale=174/25.4*72/W,height_mm=174*H/W)
 stopifnot(info$height_mm<=234)
 g<-svg_grob(info);pw<-round(width_mm/25.4*dpi);ph<-round(pw*H/W)
 # Internal vector QA exports complement the user-requested TIFF/EPS files.
 svglite::svglite(file.path(w,paste0('Fig',i,'_export_check.svg')),width=width_mm/25.4,height=info$height_mm/25.4,bg='white')
 grid.newpage();grid.draw(g);dev.off()
 grDevices::cairo_pdf(file.path(w,paste0('Fig',i,'_export_check.pdf')),width=width_mm/25.4+1/72,height=info$height_mm/25.4+1/72,family=family,bg='white')
 grid.newpage();grid.draw(g);dev.off()
 ragg::agg_tiff(paste0(prefix,'.tif'),width=pw,height=ph,units='px',res=dpi,bitsize=8,compression='lzw',background='white')
 grid.newpage();grid.draw(g);dev.off()
 grDevices::cairo_ps(paste0(prefix,'.eps'),width=width_mm/25.4+1/72,height=info$height_mm/25.4+1/72,onefile=FALSE,family=family,fallback_resolution=600,bg='white')
 grid.newpage();grid.draw(g);dev.off()
 # Precise bounding box declares the original fixed-size drawing viewport.
 eps<-readLines(paste0(prefix,'.eps'),warn=FALSE)
 bw<-width_mm/25.4*72;bh<-info$height_mm/25.4*72
 eps[grepl('^%%BoundingBox:',eps)]<-sprintf('%%%%BoundingBox: 0 0 %d %d',ceiling(bw),ceiling(bh))
 eps<-append(eps,sprintf('%%%%HiResBoundingBox: 0 0 %.8f %.8f',bw,bh),after=which(grepl('^%%BoundingBox:',eps))[1])
 writeLines(eps,paste0(prefix,'.eps'),useBytes=TRUE)
 if(!file.exists(paste0(prefix,'_preview.png')))file.copy(file.path(f,'FIGURES',paste0(name,'.png')),paste0(prefix,'_preview.png'))
 ragg::agg_png(file.path(w,paste0('Fig',i,'_review.png')),width=1200,height=round(1200*H/W),units='px',res=dpi*1200/pw,background='white')
 grid.newpage();grid.draw(g);dev.off()
 log[[name]]<-list(width_px=pw,height_px=ph,width_mm=width_mm,height_mm=info$height_mm,dpi=dpi)
 cat(name,': frozen SVG -> TIFF/EPS; geometry retained.\n',sep='')
}
for(name in c('FigureS1','FigureS2'))if(!file.exists(file.path(o,'04_SOURCE_PREVIEWS',paste0(name,'_preview.png'))))file.copy(file.path(f,'FIGURES',paste0(name,'.png')),file.path(o,'04_SOURCE_PREVIEWS',paste0(name,'_preview.png')))
write_json(log,file.path(w,'figure_export_log.json'),pretty=TRUE,auto_unbox=TRUE,digits=12)

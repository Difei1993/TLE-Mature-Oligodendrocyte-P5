options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(xml2);library(grid);library(ragg);library(png);library(base64enc);library(jsonlite)})
p <- paste0(Sys.getenv("P5_WORK_ROOT"), "/P5_MATURE_OLIGODENDROCYTE_PURE_BIOINFO_RESTART")
source_dir <- file.path(p,'MANUSCRIPT_FINAL_POLISH/FIGURES_FINAL_POLISHED')
out <- file.path(p,'MOLECULAR_NEUROBIOLOGY_SUBMISSION')
figdir <- file.path(out,'FIGURES');work <- file.path(out,'WORK')
width_mm <- 174; dpi <- 600; family <- 'Arial'
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
dir.create(work,recursive=TRUE,showWarnings=FALSE)
# Contract: reuse the frozen quantitative grids and spatial composite exactly.
# Existing claims, panel order, vector coordinates, axis ranges and numeric text
# are immutable. Only final size, panel tags and Figure 4 marker shape change.
# No statistical functions or analysis pipelines are called by this script.
attrnum <- function(n,k,default=0){v<-xml_attr(n,k);if(is.na(v))return(default);as.numeric(sub('[a-z%]+$','',v))}
style_list <- function(n){
 s<-xml_attr(n,'style');if(is.na(s))return(list())
 entries<-trimws(strsplit(s,';',fixed=TRUE)[[1]]);entries<-entries[nzchar(entries)]
 setNames(lapply(entries,function(x)trimws(sub('^[^:]+:','',x))),vapply(entries,function(x)trimws(sub(':.*$','',x)),character(1)))
}
to_triangle <- function(n){
 cx<-attrnum(n,'cx');cy<-attrnum(n,'cy');r<-attrnum(n,'r')
 # Equal-area equilateral triangle, with centroid at the original data point.
 side<-sqrt(4*pi*r*r/sqrt(3));alt<-side*sqrt(3)/2
 pts<-sprintf('%.9f,%.9f %.9f,%.9f %.9f,%.9f',cx,cy-2*alt/3,cx-side/2,cy+alt/3,cx+side/2,cy+alt/3)
 xml_name(n)<-'polygon';xml_set_attr(n,'points',pts)
 xml_set_attr(n,'cx',NA_character_);xml_set_attr(n,'cy',NA_character_);xml_set_attr(n,'r',NA_character_)
}
modify_svg <- function(name){
 doc<-read_xml(file.path(source_dir,paste0(name,'.svg')));root<-xml_root(doc)
 vb<-scan(text=xml_attr(root,'viewBox'),quiet=TRUE);W<-vb[3];H<-vb[4]
 scale<-width_mm/25.4*72/W
 xml_set_attr(root,'width','174mm');xml_set_attr(root,'height',sprintf('%.10fmm',width_mm*H/W))
 texts<-xml_find_all(doc,'.//*[local-name()="text"]')
 tags<-texts[xml_text(texts)%in%LETTERS[1:6]]
 expected<-if(name=='FigureS2')letters[1:6] else letters[1:4]
 stopifnot(identical(tolower(xml_text(tags)),expected))
 for(t in tags){
  xml_text(t)<-tolower(xml_text(t));s<-xml_attr(t,'style')
  s<-sub('font-size: [0-9.]+px',sprintf('font-size: %.8fpx',11/scale),s)
  xml_set_attr(t,'style',s);xml_set_attr(t,'textLength',NA_character_);xml_set_attr(t,'lengthAdjust',NA_character_)
 }
 changed<-list()
 if(name=='Figure4'){
  circles<-xml_find_all(doc,'.//*[local-name()="circle"]')
  for(c in circles){
   s<-xml_attr(c,'style');x<-attrnum(c,'cx');y<-attrnum(c,'cy')
   is_b<-x>300 && y<275 && grepl('#4C918B',s,fixed=TRUE)
   is_c<-y>300 && y<480 && grepl('#BD694B',s,fixed=TRUE)
   if(is_b || is_c){changed[[length(changed)+1]]<-list(panel=if(is_b)'b' else 'c',cx=x,cy=y,style=s);to_triangle(c)}
  }
  # Existing legends stay in place: their matching Myeloid and TLE symbols
  # receive the same shape as the plotted points. Astrocyte/Control stay circles.
  stopifnot(sum(vapply(changed,function(z)z$panel=='b',logical(1)))==8,
            sum(vapply(changed,function(z)z$panel=='c',logical(1)))==13)
 }
 svgfile<-file.path(figdir,paste0(name,'.svg'));write_xml(doc,svgfile)
 list(doc=doc,name=name,W=W,H=H,scale=scale,height_mm=width_mm*H/W,shape_changes=changed)
}
# Render the frozen SVG primitives using R grid graphics. This preserves the
# existing geometry rather than reconstructing plots or recomputing summaries.
# PDF text uses native Cairo text and embedded fonts, not rasterized glyphs.
svg_grob <- function(info){
 doc<-info$doc;W<-info$W;H<-info$H;sc<-info$scale
 clips<-xml_find_all(doc,'.//*[local-name()="clipPath"]')
 clipmap<-setNames(lapply(clips,function(z)xml_find_first(z,'./*[local-name()="rect"]')),xml_attr(clips,'id'))
 make_gp<-function(n,text=FALSE){
  s<-style_list(n);get<-function(k,d)if(is.null(s[[k]]))d else s[[k]]
  col<-get(if(text)'fill' else 'stroke',if(text)'#000000' else '#000000');fill<-get('fill','none')
  alpha<-as.numeric(get('opacity','1'))
  acol<-function(c,a){if(c=='none')return(NA);grDevices::adjustcolor(c,alpha.f=a)}
  dash<-get('stroke-dasharray',NULL);lty<-1
  if(!is.null(dash)){v<-as.numeric(strsplit(dash,',')[[1]]);lw<-as.numeric(get('stroke-width','1'));lty<-paste(format(as.hexmode(pmax(1,pmin(15,round(v/lw))))),collapse='')}
  gpar(col=acol(col,alpha*as.numeric(get(if(text)'fill-opacity' else 'stroke-opacity','1'))),
       fill=acol(fill,alpha*as.numeric(get('fill-opacity','1'))),
       lwd=as.numeric(get('stroke-width','1'))*sc/.75,lty=lty,lineend='round',linejoin='round',
       fontsize=as.numeric(sub('px','',get('font-size','8')))*sc,fontfamily=family,
       fontface=if(get('font-weight','normal')=='bold')'bold' else 'plain')
 }
 parse_xy<-function(n){v<-scan(text=gsub(',',' ',xml_attr(n,'points'),fixed=TRUE),quiet=TRUE);matrix(v,ncol=2,byrow=TRUE)}
 make_node<-function(n){
  tag<-xml_name(n)
  if(tag%in%c('defs','style','clipPath'))return(NULL)
  if(tag%in%c('svg','g')){
   children<-xml_children(n);gs<-list();i<-1
   while(i<=length(children)){
    ch<-children[[i]]
    # Batch consecutive identical circle marks; all coordinates are retained.
    if(xml_name(ch)=='circle'){
     s<-xml_attr(ch,'style');j<-i
     while(j<length(children) && xml_name(children[[j+1]])=='circle' && identical(xml_attr(children[[j+1]],'style'),s))j<-j+1
     z<-children[i:j]
     gs[[length(gs)+1]]<-circleGrob(x=unit(vapply(z,attrnum,numeric(1),k='cx'),'native'),y=unit(H-vapply(z,attrnum,numeric(1),k='cy'),'native'),r=unit(vapply(z,attrnum,numeric(1),k='r'),'native'),gp=make_gp(ch))
     i<-j+1;next
    }
    child<-make_node(ch);if(!is.null(child))gs[[length(gs)+1]]<-child;i<-i+1
   }
   vp<-NULL;cp<-xml_attr(n,'clip-path')
   if(!is.na(cp)){
    id<-sub('^url\\(#(.*)\\)$','\\1',cp);r<-clipmap[[id]]
    stopifnot(!is.null(r));x<-attrnum(r,'x');y<-attrnum(r,'y');w<-attrnum(r,'width');h<-attrnum(r,'height')
    vp<-viewport(x=unit(x,'native'),y=unit(H-y-h,'native'),width=unit(w,'native'),height=unit(h,'native'),just=c('left','bottom'),xscale=c(x,x+w),yscale=c(H-y-h,H-y),clip='on')
   }
   return(gTree(children=do.call(gList,gs),vp=vp))
  }
  gp<-make_gp(n,tag=='text')
  if(tag=='rect'){
   w<-if(xml_attr(n,'width')=='100%')W else attrnum(n,'width');h<-if(xml_attr(n,'height')=='100%')H else attrnum(n,'height')
   return(rectGrob(x=unit(attrnum(n,'x'),'native'),y=unit(H-attrnum(n,'y')-h,'native'),width=unit(w,'native'),height=unit(h,'native'),just=c('left','bottom'),gp=gp))
  }
  if(tag=='line')return(segmentsGrob(x0=unit(attrnum(n,'x1'),'native'),y0=unit(H-attrnum(n,'y1'),'native'),x1=unit(attrnum(n,'x2'),'native'),y1=unit(H-attrnum(n,'y2'),'native'),gp=gp))
  if(tag%in%c('polygon','polyline')){
   pts<-parse_xy(n)
   if(tag=='polygon')return(polygonGrob(x=unit(pts[,1],'native'),y=unit(H-pts[,2],'native'),gp=gp))
   return(linesGrob(x=unit(pts[,1],'native'),y=unit(H-pts[,2],'native'),gp=gp))
  }
  if(tag=='text'){
   x<-attrnum(n,'x');y<-attrnum(n,'y');angle<-0;tr<-xml_attr(n,'transform')
   if(!is.na(tr)){
    vals<-as.numeric(regmatches(tr,gregexpr('-?[0-9]+(?:\\.[0-9]+)?',tr,perl=TRUE))[[1]])
    stopifnot(length(vals)==3);x<-vals[1];y<-vals[2];angle<--vals[3]
   }
   anchor<-xml_attr(n,'text-anchor');just<-if(is.na(anchor))'left' else switch(anchor,start='left',middle='centre',end='right')
   g<-textGrob(xml_text(n),x=unit(0,'pt'),y=unit(0,'pt'),just=c(just,'bottom'),gp=gp)
   g<-editGrob(g,y=unit(0,'pt')-grobDescent(g))
   return(gTree(children=gList(g),vp=viewport(x=unit(x,'native'),y=unit(H-y,'native'),width=unit(0,'pt'),height=unit(0,'pt'),angle=angle)))
  }
  if(tag=='image'){
   href<-xml_attr(n,'href');if(is.na(href))href<-xml_attr(n,'href',ns=xml_ns(doc))
   if(is.na(href)){a<-xml_attrs(n);href<-unname(a[grepl('href$',names(a))][1])}
   raster<-readPNG(base64decode(sub('^data:image/png;base64,','',href)))
   return(rasterGrob(raster,x=unit(attrnum(n,'x'),'native'),y=unit(H-attrnum(n,'y')-attrnum(n,'height'),'native'),width=unit(attrnum(n,'width'),'native'),height=unit(attrnum(n,'height'),'native'),just=c('left','bottom'),interpolate=FALSE))
  }
  stop(paste('Unsupported frozen SVG element:',tag))
 }
 g<-make_node(xml_root(doc));g$vp<-viewport(x=unit(0,'mm'),y=unit(0,'mm'),width=unit(width_mm,'mm'),height=unit(info$height_mm,'mm'),just=c('left','bottom'),xscale=c(0,W),yscale=c(0,H),clip='on');g
}
allinfo<-list()
for(name in c('Figure1','Figure2','Figure3','Figure4','FigureS1','FigureS2')){
 info<-modify_svg(name);g<-svg_grob(info)
 # Independent editable-text rendering check. The delivery SVG retains the
 # source nodes exactly; this working render checks the R graphics pathway.
 svglite::svglite(file.path(work,paste0(name,'_render_check.svg')),width=width_mm/25.4,height=info$height_mm/25.4,bg='white')
 grid.newpage();grid.draw(g);dev.off()
 # Cairo rounds its page box to integer points. Add one blank point so the
 # exact-size viewport is never clipped; trim only the page box afterwards.
 grDevices::cairo_pdf(file.path(figdir,paste0(name,'.pdf')),width=width_mm/25.4+1/72,height=info$height_mm/25.4+1/72,family=family,bg='white',onefile=TRUE)
 grid.newpage();grid.draw(g);dev.off()
 pw<-round(width_mm/25.4*dpi);ph<-round(pw*info$H/info$W)
 ragg::agg_png(file.path(figdir,paste0(name,'.png')),width=pw,height=ph,units='px',res=dpi,background='white')
 grid.newpage();grid.draw(g);dev.off()
 ragg::agg_png(file.path(work,paste0(name,'_preview.png')),width=1200,height=round(1200*info$H/info$W),units='px',res=600*1200/pw,background='white')
 grid.newpage();grid.draw(g);dev.off()
 allinfo[[name]]<-info[c('name','W','H','scale','height_mm','shape_changes')]
 cat(name,' exported at 174 mm, PNG ',pw,' x ',ph,' pixels.\n',sep='')
}
write_json(allinfo,file.path(work,'figure_export_log.json'),pretty=TRUE,auto_unbox=TRUE,digits=12)

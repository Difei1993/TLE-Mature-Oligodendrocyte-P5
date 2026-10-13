# SPDX-License-Identifier: MIT
from pathlib import Path
import argparse
import json
import sys

def load_display_paths():
    parser = argparse.ArgumentParser(description="Display S3-S5 from matching archived inputs; no biological analysis.")
    parser.add_argument("--config", type=Path, help="JSON configuration; relative paths resolve against its directory")
    parser.add_argument("--analysis-root", type=Path, help="Local MANUAL_R_ANALYSIS archive root")
    parser.add_argument("--frozen-table-dir", type=Path, help="Directory with the four archived display TSV inputs")
    parser.add_argument("--umap-export-dir", type=Path, help="Private, existing saved-coordinate/label exports")
    parser.add_argument("--output-dir", type=Path, help="New empty local output directory; no input directory may be overwritten")
    parser.add_argument("--check-inputs-only", action="store_true", help="Check file availability only; do not open data, create outputs, or draw")
    args = parser.parse_args()
    config = {}
    if args.config:
        try:
            config = json.loads(args.config.read_text(encoding="utf-8-sig"))
        except (OSError, ValueError) as error:
            parser.error("Cannot read configuration: " + str(error))
        if not isinstance(config, dict):
            parser.error("Configuration must be a JSON object")
    values = {}
    for key in ("analysis_root", "frozen_table_dir", "umap_export_dir", "output_dir"):
        direct = getattr(args, key)
        value = direct if direct is not None else config.get(key)
        if not isinstance(value, (str, Path)) or not str(value).strip():
            parser.error("Missing path: " + key)
        path = Path(value).expanduser()
        if direct is None and args.config and not path.is_absolute():
            path = args.config.resolve().parent / path
        values[key] = path.resolve()
    return parser, args, values

_parser, _args, _paths = load_display_paths()
B = _paths["analysis_root"]
O = _paths["output_dir"]
TABLES = _paths["frozen_table_dir"]
UMAP_EXPORTS = _paths["umap_export_dir"]
_analysis_files = [
    "03_figures/human_harmony_level1_annotation_umap.png",
    "03_figures/human_harmony_umap_resolution_grid.png",
    "03_figures/human_level1_canonical_marker_dotplot.png",
    "03_figures/human_harmony_r0p4_dotplot_glial.png",
    "35C_P5_CROSS_REFERENCE/06_GSE140393_FINAL_MANUAL_OLIGO_ANNOTATION/figures/GSE140393_FINAL_OLIGO_MARKER_DOTPLOT.png",
] + [
    "04_tables/32B2C_FILTERED_CONSENSUS_AUDIT_20260728_215801/_FILTERED_CONSENSUS_REVIEW/"
    + "Oligodendrocytes/" + scheme + "/dt_0_01/Oligodendrocytes_" + scheme
    + "_final_k6.clustering.k_6.dt_0_01.png"
    for scheme in ("primary", "sensitivity")
]
_table_files = ["S3C_display_label_to_public_sample.tsv", "historical_K_selection_statistics.tsv",
                "S4C_frozen_top8_OL_genes.tsv", "S4F_V3_TableS3_21_primary_core_effects.tsv"]
_export_files = ["GSE140393_" + view + suffix for view in ("all_candidates", "bona_fide")
                 for suffix in ("_saved_umap.csv", "_coordinates_double.bin", "_original_discrete_colors.csv")]
_required = [("analysis/" + name, B / name) for name in _analysis_files]
_required += [("frozen_tables/" + name, TABLES / name) for name in _table_files]
_required += [("private_umap_exports/" + name, UMAP_EXPORTS / name) for name in _export_files]
_missing = [name for name, path in _required if not path.is_file()]
if _args.check_inputs_only:
    print(json.dumps({"mode": "AVAILABILITY_CHECK_ONLY", "required": len(_required),
                      "available": len(_required) - len(_missing), "missing": _missing,
                      "data_opened": False, "output_created": False, "render_executed": False}, indent=2))
    raise SystemExit(2 if _missing else 0)
if _missing:
    _parser.error("Archived inputs are unavailable: " + ", ".join(_missing))
for _input_root in (B, TABLES, UMAP_EXPORTS):
    if O == _input_root or O in _input_root.parents:
        _parser.error("Output directory cannot equal or contain an input directory")
if O.exists() and (not O.is_dir() or any(O.iterdir())):
    _parser.error("Output directory must be new or empty; existing archived/submission results will not be overwritten")
FIG = O / "02_SUPPLEMENTARY_FIGURES"
QA = O / "_INTERNAL_QA"
DATA = O / "05_SOURCE_AND_CHANGE_RECORDS/source_data"
for _directory in (FIG, QA, DATA):
    _directory.mkdir(parents=True, exist_ok=True)

"""Display only: unchanged archived pixel layers, saved coordinates, and frozen tables.
No clustering, normalization, fitting, enrichment or statistical testing is performed.
"""
from pathlib import Path
from collections import Counter
import csv,json,hashlib,sys
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from scipy import ndimage
from p5_display_image_layers import dotplot_layers,plot_geometry,groups
sys.stdout.reconfigure(encoding='utf8')
WIDTH=174/25.4
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'lines.linewidth':1,'legend.frameon':False})
def rows(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig'),delimiter='\t' if Path(p).suffix=='.tsv' else ','))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def tsv(p,rs):
    with Path(p).open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
records=[];fig_pages={}
def new(height):return plt.figure(figsize=(WIDTH,height/25.4),facecolor='white')
def label(fig,letter,x,y):fig.text(x,y,letter,fontsize=10,weight='bold',ha='left',va='top')
def output(fig,snum,part):
    root=FIG/f'Figure_S{snum}_{part}';fig.savefig(root.with_suffix('.pdf'),metadata={'CreationDate':None});fig.savefig(root.with_suffix('.svg'))
    fig.savefig(root.with_suffix('.png'),dpi=600)
    fig.savefig(root.with_suffix('.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    # Submission export: opaque eight-bit RGB rather than an unused alpha channel.
    for ext in ['.png','.tiff']:
        with Image.open(root.with_suffix(ext)) as im:
            if im.mode=='RGBA':assert im.getchannel('A').getextrema()==(255,255)
            rgb=im.convert('RGB')
        opts={'dpi':(600,600)}
        if ext=='.tiff':opts['compression']='tiff_lzw'
        rgb.save(root.with_suffix(ext),**opts)
    fig_pages.setdefault(snum,[]).append(root.with_suffix('.pdf'));plt.close(fig)
def record(s,panel,source,operation,details):
    records.append({'figure':f'S{s}','panel':panel,'source':str(source),'source_sha256':sha(source),'operation':operation,'details':details,'scientific_reanalysis':'NONE'})

def archived_umap(fig,image,box,xvals,yvals,letter,source):
    geo=plot_geometry(image);x0,y0,x1,y1=geo['box'];roi=image[y0:y1,x0:x1]
    left,bottom,width,height=box
    actual=width*fig.get_figwidth()*roi.shape[0]/roi.shape[1]/fig.get_figheight()
    bottom=bottom+height-actual;height=actual
    assert bottom>0 and bottom+height<1
    ax=fig.add_axes([left,bottom,width,height]);ax.imshow(roi,extent=(0,roi.shape[1],roi.shape[0],0),interpolation='nearest',aspect='auto')
    assert len(geo['xticks'])==len(xvals) and len(geo['yticks'])==len(yvals),(source,len(geo['xticks']),len(geo['yticks']))
    ax.set_xticks(geo['xticks'],xvals);ax.set_yticks(geo['yticks'],yvals);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.tick_params(length=2,pad=2)
    label(fig,letter,max(.012,left-.075),min(.995,bottom+height+.032))
    record(3,letter,source,'Preserved archived UMAP pixel layer; axis labels re-typeset',json.dumps({'pixel_crop':geo['box'],'no_point_or_cluster_change':True}))

def original_dotplot(fig,source,features,annotations,percents,box,letter,snum):
    d=dotplot_layers(source,len(features),len(annotations),len(percents));im=d['image'];x0,y0,x1,y1=d['plot_box'];roi=im[y0:y1,x0:x1]
    left,bottom,width,height=box
    actual=width*fig.get_figwidth()*roi.shape[0]/roi.shape[1]/fig.get_figheight()
    bottom=bottom+height-actual;height=actual
    ax=fig.add_axes([left,bottom,width,height])
    ax.imshow(roi,extent=(0,roi.shape[1],roi.shape[0],0),interpolation='nearest',aspect='auto')
    ax.set_xticks(d['xticks'],features,rotation=70,ha='right',rotation_mode='anchor');ax.set_yticks(d['yticks'],annotations)
    ax.tick_params(length=2,pad=2);ax.spines[['top','right']].set_visible(False)
    legx=left+width+.028
    fig.text(legx,bottom+height,'Percent\nexpressed',fontsize=8,va='top')
    pt_per_pixel=(width*fig.get_figwidth()*72)/roi.shape[1]
    for i,(part,percent) in enumerate(zip(d['percent'],percents)):
        x,y,xx,yy=part['box'];patch=im[max(0,y-2):yy+2,max(0,x-2):xx+2]
        wi=patch.shape[1]*pt_per_pixel/72/fig.get_figwidth();hi=patch.shape[0]*pt_per_pixel/72/fig.get_figheight()
        center=bottom+height-.052-i*.025
        a=fig.add_axes([legx+.015-wi/2,center-hi/2,wi,hi]);a.imshow(patch,interpolation='nearest');a.axis('off')
        fig.text(legx+.055,center,str(percent),fontsize=8,va='center')
    bx,by,bxx,byy=d['bar_box'];gradient=im[by:byy,bx:bxx]
    gh=.07 if fig.get_figheight()>8 else .18
    cy=bottom+height-.052-(len(percents)-1)*.025-.034
    fig.text(legx,cy,'Scaled mean\nexpression',fontsize=8,va='top')
    ticks=['2','1','0'] if snum==3 else ['1.0','0.5','0.0','-0.5','-1.0']
    assert len(d['bar_tick_y'])==len(ticks),(source,d['bar_tick_y'])
    bar=fig.add_axes([legx,cy-.045-gh,.025,gh]);bar.imshow(gradient,extent=(0,gradient.shape[1],gradient.shape[0],0),interpolation='nearest',aspect='auto')
    bar.set_yticks(d['bar_tick_y'],ticks);bar.yaxis.tick_right();bar.tick_params(length=0,pad=3);bar.set_xticks([])
    for sp in bar.spines.values():sp.set_visible(False)
    label(fig,letter,.012,bottom+height+.026)
    record(snum,letter,source,'Unchanged original expression-dot and scale pixel layers; typography/legend repositioned',json.dumps({'plot_crop':d['plot_box'],'color_scale_crop':d['bar_box'],'percent_marker_scale_matches_data':True,'no_expression_summary_recalculation':True}))

global_source=B/'03_figures/human_harmony_level1_annotation_umap.png'
grid_source=B/'03_figures/human_harmony_umap_resolution_grid.png'
global_im=np.asarray(Image.open(global_source).convert('RGB'))
# Exclude the existing legend solely to replace it with readable typography; retain every plotted point.
global_plot=global_im[:,:2915]
grid_im=np.asarray(Image.open(grid_source).convert('RGB'));cluster_plot=grid_im[:1800,2250:4500]
fig=new(234)
archived_umap(fig,global_plot,[.085,.57,.69,.385],[-15,-10,-5,0,5,10],[10,0,-10],'a',global_source)
names=['Excitatory neurons','Inhibitory neurons','Astrocytes','Oligodendrocytes','OPC','Microglia/macrophage','Vascular cells','T lymphocytes']
colors=['#F8766D','#CD9600','#7CAE00','#00BE67','#00BFC4','#00A9FF','#C77CFF','#F564E3']
handles=[Line2D([],[],ls='',marker='o',markersize=4,color=c,label=n) for c,n in zip(colors,names)]
fig.legend(handles=handles,loc='center left',bbox_to_anchor=(.775,.76),fontsize=8,borderaxespad=0,labelspacing=.8,handletextpad=.4)
archived_umap(fig,cluster_plot,[.105,.048,.855,.47],[-15,-10,-5,0,5,10],[10,0,-10],'b',grid_source)
output(fig,3,1)

fig=new(234)
canonical=['SLC17A7','SATB2','GAD1','GAD2','AQP4','ALDH1L1','PLP1','MOBP','PDGFRA','VCAN','P2RY12','INPP5D','CLDN5','PECAM1','RGS5','PDGFRB','CD3D','CD3E']
original_dotplot(fig,B/'03_figures/human_level1_canonical_marker_dotplot.png',canonical,list(reversed(names)),[0,25,50,75],[.225,.722,.58,.235],'c',3)
glial=['AQP4','ALDH1L1','SLC1A2','SLC1A3','GFAP','GJA1','PLP1','MBP','MOG','MAG','MOBP','SOX10','PDGFRA','CSPG4','VCAN','OLIG1','OLIG2']
original_dotplot(fig,B/'03_figures/human_harmony_r0p4_dotplot_glial.png',glial,[str(i) for i in range(18,-1,-1)],[0,25,50,75],[.225,.35,.58,.285],'d',3)
source=TABLES/'S3C_display_label_to_public_sample.tsv';ds=rows(source)
ax=fig.add_axes([.15,.045,.81,.155]);vals=[int(r['OL_nuclei']) for r in ds]
ax.bar(range(8),vals,color=['#4477AA']*4+['#CC6677']*4,width=.62);ax.axhline(450,color='#333333',ls='--',lw=.8)
for i,v in enumerate(vals):ax.text(i,v+100,f'{v:,}',ha='center',fontsize=8)
ax.set_xticks(range(8),[r['display_label'].replace('Control ','C').replace('TLE ','T') for r in ds]);ax.set_ylim(0,3850);ax.set_ylabel('OL nuclei');ax.tick_params(length=2)
ax.text(7.4,510,'450',ha='left',fontsize=8);label(fig,'e',.012,.225)
record(3,'e',source,'Frozen Table S1 counts; existing 450/donor sampling annotation','No sampling or abundance test');tsv(DATA/'S3e_donor_availability.tsv',ds);output(fig,3,2)

stats_source=TABLES/'historical_K_selection_statistics.tsv';ks=rows(stats_source)
primary=[r for r in ks if r['sampling']=='primary'];sens=[r for r in ks if r['sampling']=='sensitivity']
assert all(float(r['local_density_threshold'])==.5 for r in ks)
fig=new(234);a=fig.add_axes([.12,.675,.75,.265]);a2=a.twinx();a2.spines['top'].set_visible(False)
a.plot([float(r['k']) for r in primary],[float(r['silhouette']) for r in primary],'o-',color='#4477AA',ms=3)
a2.plot([float(r['k']) for r in primary],[float(r['prediction_error'])/1e6 for r in primary],'s-',color='#BB5566',ms=3)
a.set_xticks([4,6,8,10,12,14]);a.set_xlabel('Number of components K');a.set_ylabel('Stability (silhouette)',color='#4477AA');a2.set_ylabel('Reconstruction error ($10^6$)',color='#BB5566');a.tick_params(length=2);a2.tick_params(length=2)
label(fig,'a',.012,.982)
for box,key,color,letter in [([.13,.413,.335,.155],'silhouette','#4477AA','b'),([.62,.413,.335,.155],'prediction_error','#BB5566','c')]:
    ax=fig.add_axes(box);values=[float(r[key]) for r in sens];disp=values if key=='silhouette' else [x/1e6 for x in values]
    ax.plot([6,7],disp,'o-',color=color,ms=3);ax.set_xlim(5.85,7.15);ax.set_xticks([6,7]);ax.set_xlabel('K')
    ax.set_ylabel('Stability' if key=='silhouette' else 'Error ($10^6$)')
    pad=(max(disp)-min(disp))*.30;ax.set_ylim(min(disp)-pad,max(disp)+pad*2)
    for x,v in zip([6,7],disp):ax.annotate(f'{v:.6f}',(x,v),xytext=(0,5),textcoords='offset points',ha='center',fontsize=8)
    label(fig,letter,box[0]-.10,.611)
top_source=TABLES/'S4C_frozen_top8_OL_genes.tsv';top=rows(top_source)
assert len(top)==48
ax=fig.add_axes([.04,.04,.925,.245]);ax.axis('off')
data=[[next(r['gene'] for r in top if r['program_id']==f'Oligodendrocytes__P{i}' and int(r['rank'])==j) for i in range(1,7)] for j in range(1,9)]
tab=ax.table(cellText=data,colLabels=[f'P{i}' for i in range(1,7)],rowLabels=[str(i) for i in range(1,9)],loc='center',cellLoc='center',bbox=[.025,0,.975,1]);tab.auto_set_font_size(False);tab.set_fontsize(8)
for (r,c),cell in tab.get_celld().items():
    cell.set_linewidth(.4);cell.set_edgecolor('#CCCCCC')
    if r==0:cell.set_facecolor('#EBEFF3');cell.set_text_props(weight='bold')
    if c==4:cell.set_facecolor('#FFF1DF')
label(fig,'d',.012,.322)
record(4,'a-c',stats_source,'Display saved K metrics; split sensitivity metrics into independent axes','Pilot threshold0.5; raw saved numbers unchanged')
record(4,'d',top_source,'Preserve all six Top8 lists, including neuronal-associated transcripts','48/48 genes retained; no reranking or enrichment');tsv(DATA/'S4a_c_archived_K_metrics.tsv',ks);tsv(DATA/'S4d_frozen_top8_genes.tsv',top);output(fig,4,1)

for index,scheme,panel,ticks in [(2,'primary','e',['0.000','0.628','1.257']),(3,'sensitivity','f',['0.000','0.642','1.285'])]:
    source=B/f'04_tables/32B2C_FILTERED_CONSENSUS_AUDIT_20260728_215801/_FILTERED_CONSENSUS_REVIEW/Oligodendrocytes/{scheme}/dt_0_01/Oligodendrocytes_{scheme}_final_k6.clustering.k_6.dt_0_01.png'
    im=np.asarray(Image.open(source).convert('RGB'));assert im.shape[:2]==(2375,3750)
    fig=new(142);label(fig,panel,.012,.978)
    matrix=im[43:2360,35:2345];ax=fig.add_axes([.025,.095,.615,.82]);ax.imshow(matrix,interpolation='nearest');ax.axis('off')
    # Exact historical histogram/colorbar pixels retained; labels re-typeset at 8 pt.
    histbox=(2462,169,3436,896);x0,y0,x1,y1=histbox;hist=im[y0:y1,x0:x1].copy()
    # Replace only the archived black text in a verified blank label region, not histogram bars/threshold line.
    oldlabel=hist[30:108,70:248]
    assert np.max(oldlabel.max(axis=2).astype(int)-oldlabel.min(axis=2).astype(int))==0
    hist[30:108,70:248]=255
    ax=fig.add_axes([.665,.51,.265,.365]);ax.imshow(hist,extent=(0,hist.shape[1],hist.shape[0],0),aspect='auto',interpolation='nearest')
    dark=im.max(axis=2)<100
    xt=groups(dark[y1+2:y1+16,x0+1:x1-1].sum(axis=0)>=6)
    yt=groups(dark[y0:y1,x1+3:x1+15].sum(axis=1)>=5)
    xt=[g.mean()+1 for g in xt];yt=[g.mean() for g in yt]
    assert len(xt)==6 and len(yt)==7,(scheme,len(xt),len(yt))
    ax.set_xticks(xt,['0.0','0.2','0.4','0.6','0.8','1.0']);ax.set_yticks(yt,[600,500,400,300,200,100,0]);ax.yaxis.tick_right();ax.tick_params(length=2,pad=2)
    ax.set_xlabel('Mean distance to\nnearest spectra',fontsize=8)
    ax.text(.14,.94,'Filtering\nthreshold',transform=ax.transAxes,va='top',fontsize=8)
    fig.text(.69,.93,'Local density',fontsize=8)
    fig.text(.69,.365,('0' if scheme=='primary' else '7')+'/600 spectra\nremoved',fontsize=8)
    gradient=im[1262:1336,2463:3436];cb=fig.add_axes([.665,.17,.265,.04]);cb.imshow(gradient,aspect='auto',interpolation='nearest');cb.axis('off')
    fig.text(.68,.238,'Euclidean distance',fontsize=8)
    for x,t,ha in zip([.665,.7975,.93],ticks,['left','center','right']):fig.text(x,.13,t,ha=ha,fontsize=8)
    record(4,panel,source,'Unchanged final consensus matrix/histogram/color strip; readable labels','K6 dt0.01; spectra counts and metric ticks copied from archive; no consensus recalculated')
    output(fig,4,index)

source=TABLES/'S4F_V3_TableS3_21_primary_core_effects.tsv';effects=rows(source);assert len(effects)==21
fig=new(200);ax=fig.add_axes([.22,.115,.63,.815]);abbr={'Astrocytes':'Astro','Microglia_macrophage':'Myeloid','OPC':'OPC','Oligodendrocytes':'OL'}
for i,r in enumerate(effects):
    v=float(r['effect_TLE_minus_Normal']);color='#D55E00' if r['program_id']=='Oligodendrocytes__P5' else '#4477AA'
    ax.plot([0,v],[i,i],color=color,lw=1);ax.scatter(v,i,color=color,s=16);ax.text(1.045,i,f"{float(r['formal_primary_fdr']):.2f}",transform=ax.get_yaxis_transform(),va='center',fontsize=8)
ax.axvline(0,color='#777777',lw=.7);ax.set_yticks(range(21),[abbr[r['lineage']]+' '+r['program_id'].split('__')[-1] for r in effects]);ax.set_ylim(20.7,-.7);ax.set_xlabel('Primary CLR-usage difference\n(TLE - control)');ax.text(1.02,1.025,'BH FDR',transform=ax.transAxes,fontsize=8);label(fig,'g',.012,.982)
assert next(r for r in effects if r['program_id']=='Oligodendrocytes__P5')['formal_primary_fdr']=='0.6'
record(4,'g',source,'All21 frozen primary effects/FDR; zero-to-effect segments','Segments are not intervals; no recalculated tests/BH; CLR effects are not absolute gene expression');tsv(DATA/'S4g_21_program_frozen_effects_FDR.tsv',effects);output(fig,4,4)

aliases={'B_cell':'B cell','Cycling_lineage_unresolved':'Cycling, unresolved','Endothelial_cell':'Endothelial cell','ENPP6_maturing_oligodendrocyte':'ENPP6-maturing OL','Fibroblast':'Fibroblast','GPR17_transitional_oligodendrocyte':'GPR17-transitional OL','Mature_oligodendrocyte':'Mature OL','Microglia_macrophage':'Microglia/macrophage','Microglia_vascular_mixed':'Microglia/vascular mixed','Mural_fibroblast_cell':'Mural/fibroblast cell','OPC':'OPC','Plasma_cell':'Plasma cell','Plasmacytoid_dendritic_cell':'Plasmacytoid dendritic cell','T_NK_cell':'T/NK cell'}
fig=new(234);stage_data={}
for name,box,panel in [('all_candidates',[.105,.557,.52,.375],'a'),('bona_fide',[.105,.095,.52,.375],'b')]:
    inp=UMAP_EXPORTS/f'GSE140393_{name}_saved_umap.csv';rs=rows(inp)
    xy=np.fromfile(UMAP_EXPORTS/f'GSE140393_{name}_coordinates_double.bin',dtype='<f8').reshape((len(rs),2),order='F')
    assert np.all(np.isfinite(xy))
    assert np.all((xy[:,0]>-16.2)&(xy[:,0]<19.8)&(xy[:,1]>-21)&(xy[:,1]<21)), 'No saved point may be clipped'
    labels=[r['manual_label'] for r in rs];counts=Counter(labels);stage_data[name]=dict(counts)
    rawcols=rows(UMAP_EXPORTS/f'GSE140393_{name}_original_discrete_colors.csv')
    # Restore the archived ggplot legend's case-insensitive order, independent of current R locale.
    cols=[{'manual_label':lab,'color':r['color']} for lab,r in zip(sorted(counts,key=str.casefold),rawcols)]
    tsv(DATA/f'S5{panel}_retained_discrete_palette.tsv',cols)
    ax=fig.add_axes(box)
    for r in cols:
        idx=np.array([v==r['manual_label'] for v in labels]);assert idx.any()
        ax.scatter(xy[idx,0],xy[idx,1],s=1.5,color=r['color'],linewidths=0,rasterized=True)
    ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.set_xlim(-16.2,19.8);ax.set_ylim(-21,21);ax.set_xticks([-10,0,10]);ax.set_yticks([-20,-10,0,10,20]);ax.tick_params(length=2,pad=2)
    handles=[Line2D([],[],ls='',marker='o',markersize=3.5,color=r['color'],label=aliases[r['manual_label']]+f" ({counts[r['manual_label']]:,})") for r in cols]
    fig.legend(handles=handles,loc='center left',bbox_to_anchor=(.65,box[1]+box[3]/2),fontsize=8,borderaxespad=0,handletextpad=.35,labelspacing=.57)
    label(fig,panel,.012,box[1]+box[3]+.035)
    record(5,panel,inp,'Replot only stored double-precision UMAP coordinates and original manual labels','No point removed, no UMAP recalculated; original per-view discrete palette retained; labels moved to legend')
assert stage_data['all_candidates']['Mature_oligodendrocyte']==1381 and sum(stage_data['bona_fide'].values())==2138
tsv(DATA/'S5_saved_label_counts.tsv',[{'view':name,'original_label':k,'display_label':aliases[k],'cells':v} for name,data in stage_data.items() for k,v in data.items()]);output(fig,5,1)
fig=new(130)
features=['SOX10','OLIG1','OLIG2','PDGFRA','CSPG4','VCAN','PTPRZ1','GPR17','BCAS1','TCF7L2','ENPP6','PLP1','CNP','MOBP','MOG','MBP','MAG','CLDN11','ERMN','CRYAB','APOD','PTGDS','SCD','SELENOP','FTH1','FTL']
original_dotplot(fig,B/'35C_P5_CROSS_REFERENCE/06_GSE140393_FINAL_MANUAL_OLIGO_ANNOTATION/figures/GSE140393_FINAL_OLIGO_MARKER_DOTPLOT.png',features,['OPC','Mature OL','GPR17-transitional OL','ENPP6-maturing OL'],[0,25,50,75,100],[.225,.245,.58,.58],'c',5);output(fig,5,2)

from pypdf import PdfReader,PdfWriter
for snum,paths in fig_pages.items():
    w=PdfWriter()
    for p in paths:w.add_page(PdfReader(p).pages[0])
    w.write(FIG/f'Figure_S{snum}.pdf')
tsv(O/'05_SOURCE_AND_CHANGE_RECORDS/P5_FIGURE_PANEL_SOURCE_AND_CHANGES.tsv',records)
(QA/'figure_leaf_manifest.json').write_text(json.dumps({str(k):[str(p) for p in v] for k,v in fig_pages.items()},indent=2),encoding='utf8')
print(json.dumps({'figures':list(fig_pages),'leaves':{k:len(v) for k,v in fig_pages.items()},'scientific_analysis':'NOT EXECUTED','all_frozen_genes_and_effects_preserved':True},indent=2))

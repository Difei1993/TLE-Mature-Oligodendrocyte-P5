# SPDX-License-Identifier: MIT
"""Presentation-only reuse of unchanged image data layers and labels."""
from pathlib import Path
from collections import Counter
import numpy as np
from PIL import Image
from scipy import ndimage

def groups(v):
    idx=np.flatnonzero(v)
    if not len(idx):return []
    return [g for g in np.split(idx,np.flatnonzero(np.diff(idx)>1)+1) if len(g)]

def plot_geometry(image):
    dark=image.max(axis=2)<100
    yline=int(np.argmax(dark.sum(axis=1)))
    segment=max(groups(dark[yline]),key=len);x0=int(segment[0]);xr=int(segment[-1])
    vertical=max(groups(dark[:,x0]),key=len);y0=int(vertical[0])
    xt=groups(dark[yline+3:yline+16,x0+1:xr].sum(axis=0)>=7)
    yt=groups(dark[y0:yline,x0-12:x0-2].sum(axis=1)>=5)
    return {'box':(x0,y0,xr+1,yline+1),'xticks':np.array([g.mean()+x0+1 for g in xt])-x0,'yticks':np.array([g.mean()+y0 for g in yt])-y0}

def dotplot_layers(path,nfeatures,ngroups,npercent):
    im=np.asarray(Image.open(path).convert('RGB'))
    dark=im.max(axis=2)<100
    yline=int(np.argmax(dark.sum(axis=1)))
    segment=max(groups(dark[yline]),key=len);x0=int(segment[0]);xr=int(segment[-1])
    vertical=max(groups(dark[:,x0]),key=len);y0=int(vertical[0]);yb=int(vertical[-1])
    if abs(yb-yline)>5: raise ValueError(('axis endpoints',path,yb,yline))
    xt=groups(dark[yline+3:yline+16,x0+1:xr].sum(axis=0)>=7)
    yt=groups(dark[y0:yline,x0-12:x0-2].sum(axis=1)>=5)
    xticks=np.array([g.mean()+x0+1 for g in xt]);yticks=np.array([g.mean()+y0 for g in yt])
    if len(xticks)!=nfeatures or len(yticks)!=ngroups:raise ValueError(('ticks',Path(path).name,len(xticks),len(yticks),nfeatures,ngroups))
    labels,n=ndimage.label(dark[:,xr+20:]);boxes=ndimage.find_objects(labels);candidates=[]
    for k,sl in enumerate(boxes,1):
        if sl is None:continue
        sy,sx=sl;w=sx.stop-sx.start;h=sy.stop-sy.start
        if .75<=w/max(h,1)<=1.3 and w>=2:
            fill=np.count_nonzero(labels[sl]==k)/(w*h)
            if fill>.67:candidates.append({'box':(sx.start+xr+20,sy.start,sx.stop+xr+20,sy.stop),'x':(sx.start+sx.stop)/2+xr+20,'y':(sy.start+sy.stop)/2,'w':w,'h':h,'area':w*h})
    largest=max(candidates,key=lambda x:x['area']);xc=largest['x']
    aligned=sorted([q for q in candidates if abs(q['x']-xc)<8 and q['area']<=largest['area']],key=lambda x:x['y'])
    aligned=[q for q in aligned if y0<q['y']<largest['y']+5]
    aligned=aligned[-npercent:]
    if len(aligned)!=npercent:raise ValueError(('percent legend',path,len(aligned)))
    colorful=(im.max(axis=2).astype(float)-im.min(axis=2)>25)&(np.indices(im.shape[:2])[1]>xr)
    lab,n=ndimage.label(colorful);bs=ndimage.find_objects(lab)
    bars=[]
    for k,sl in enumerate(bs,1):
        if sl is not None:
            sy,sx=sl;w=sx.stop-sx.start;h=sy.stop-sy.start
            if h>100 and w>15 and h/w>2:bars.append((w*h,sl))
    _,(sy,sx)=max(bars,key=lambda x:x[0])
    barbox=(sx.start,sy.start,sx.stop,sy.stop)
    scan0=max(0,sy.start-25);scan1=min(im.shape[0],sy.stop+25)
    tickgroups=groups(dark[scan0:scan1,sx.stop+5:min(im.shape[1],sx.stop+125)].sum(axis=1)>2)
    tickcenters=[float(g.mean()+scan0-sy.start) for g in tickgroups if len(g)>5 and g.mean()+scan0>=sy.start-5]
    return {'image':im,'plot_box':(x0,y0,xr+1,yline+1),'xticks':xticks-x0,'yticks':yticks-y0,'percent':aligned,'bar_box':barbox,'bar_tick_y':tickcenters}

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Inspect geometry of three archived DotPlots; no expression statistics are recalculated.')
    parser.add_argument('--analysis-root', type=Path, required=True, help='Local analysis archive root')
    b = parser.parse_args().analysis_root.expanduser().resolve()
    for p,nx,ny,npct in [('03_figures/human_level1_canonical_marker_dotplot.png',18,8,4),('03_figures/human_harmony_r0p4_dotplot_glial.png',17,19,4),('35C_P5_CROSS_REFERENCE/06_GSE140393_FINAL_MANUAL_OLIGO_ANNOTATION/figures/GSE140393_FINAL_OLIGO_MARKER_DOTPLOT.png',26,4,5)]:
        d=dotplot_layers(b/p,nx,ny,npct)
        print(p,d['plot_box'],'ticks',len(d['xticks']),len(d['yticks']),'percent',[(q['w'],q['h']) for q in d['percent']], 'bar',d['bar_box'])

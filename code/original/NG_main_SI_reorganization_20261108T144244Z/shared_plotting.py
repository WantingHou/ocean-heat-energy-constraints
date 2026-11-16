"""Offline, create-only figures from frozen inputs and display-only postprocessing.
No downloads, fits, optimizers, source-history changes or circulation simulations.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,LogNorm,ListedColormap,BoundaryNorm,FuncNorm,LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch,Rectangle,Polygon
from matplotlib.projections.geo import GeoAxes
from PIL import Image,ImageOps,ImageDraw
ap=argparse.ArgumentParser();ap.add_argument('--source-data',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--figure',choices=['all','2','3','4','5'],default='all');ap.add_argument('--dpi',type=int,choices=[450,600],default=450);args=ap.parse_args()
S=args.source_data;R=S/'retained';O=args.output_dir
if O.exists():raise FileExistsError(O)
O.mkdir(parents=True)
BLUE='#2a78d6';ORANGE='#eb6834';GREEN='#1baf7a';INK='#263440';GRAY='#65727b';LAND='#e4e1da';UNKNOWN='#c3c9cf'
COLORS=[BLUE,ORANGE,GREEN];OPS=['direct','loss_100yr','loss_30yr'];OPN=['Complete retention','100-year loss','30-year loss'];ALLOC=['geological_power','geological_area'];ALN=['Geological power','Geological area']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.0,'axes.labelsize':7.5,'xtick.labelsize':7.0,'ytick.labelsize':7.0,'legend.fontsize':6.7,'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':'#a7b0b7','axes.linewidth':.5,'pdf.fonttype':42,'svg.fonttype':'none','svg.hashsalt':'NG-five-20261106','figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','hatch.linewidth':.35})
GeoAxes.RESOLUTION=2
manifest=[];checks=[];files=[]
def read(name):return pd.read_csv(R/name,float_precision='round_trip')
def record(panel,source,fields,window,depth,units,condition,operator,denominator,screen,identity,notes=''):
 manifest.append(dict(panel=panel,source_files=source,fields=fields,time_window=window,depth_m=depth,units=units,source_conditions=condition,response_operator=operator,denominator=denominator,joint_budget_screening=screen,candidate_model_observation_ids=identity,interpretation=notes))
def clean(ax):
 ax.spines[['top','right']].set_visible(False);ax.grid(color='#e6e9ec',lw=.4);ax.set_axisbelow(True);ax.tick_params(length=2.5,width=.5)
def letter(ax,label):ax.text(-.065,1.04,label,transform=ax.transAxes,weight='bold',fontsize=10,va='bottom')
def save(fig,name):
 import io,html,xml.etree.ElementTree as ET
 for ax in fig.axes:
  for patch in ax.patches:
   if patch.get_gid() in ['ppt-Label — contact','ppt-Label — fluid rise','ppt-Label — flow']:patch.set_gid(patch.get_gid().replace('ppt-Label — ','ppt-mask-'))
 fig.canvas.draw();rend=fig.canvas.get_renderer();outside=[]
 for t in fig.findobj(matplotlib.text.Text):
  if t.get_visible() and t.get_text():
   bb=t.get_window_extent(rend)
   if bb.x0 < -2 or bb.y0 < -2 or bb.x1>fig.bbox.width+2 or bb.y1>fig.bbox.height+2:outside.append(t.get_text())
 checks.append(dict(figure=name,text_outside_canvas=outside,minimum_font_pt=min(t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text())))
 for ext in ['pdf','svg','png']:
  p=O/(name+'.'+ext)
  if ext=='svg':
   buffer=io.StringIO();fig.savefig(buffer,format='svg',dpi=args.dpi);svg=buffer.getvalue()
   for gid,title in CELL_TITLES.items():svg=svg.replace('<g id="'+gid+'">','<g id="'+gid+'"><title>'+html.escape(title)+'</title>')
   ids=[e.attrib['id'] for e in ET.fromstring(svg).iter() if 'id' in e.attrib];assert len(ids)==len(set(ids)),name
   with p.open('x') as f:f.write(svg)
  else:
   with p.open('xb') as f:fig.savefig(f,format=ext,dpi=args.dpi)
 with (O/(name+'_preview.png')).open('xb') as f:fig.savefig(f,format='png',dpi=160)
 plt.close(fig)
geo=np.load(R/'F1a_geographic_grids.npz',allow_pickle=False)
# Pacific-centred longitude, registered identically for data and coastlines.
x0=geo['terrain_lon']%360-180;tx,ix=np.unique(x0,return_index=True);tz=geo['terrain_m'][:,ix];tx=np.r_[tx,tx[0]+360];tz=np.column_stack([tz,tz[:,0]])
def coast(ax):
 ax.contourf(np.deg2rad(tx),np.deg2rad(geo['terrain_lat']),tz,levels=[0,10000],colors=[LAND],zorder=4)
 ax.contour(np.deg2rad(tx),np.deg2rad(geo['terrain_lat']),tz,levels=[0],colors=['#7e8b94'],linewidths=.22,zorder=5)
 ax.set_xticks([]);ax.set_yticks([]);ax.spines['geo'].set_linewidth(.4)
def mapplot(fig,rect,values,lon,lat,cmap,norm,label):
 ax=fig.add_axes(rect,projection='mollweide');order=np.argsort(lon%360-180);dx=float(np.diff(lon)[0]);dy=float(np.diff(lat)[0]);xx=np.r_[lon[order]%360-180-dx/2,180];yy=np.r_[lat-dy/2,90]
 cm=(plt.get_cmap(cmap) if isinstance(cmap,str) else cmap).copy();cm.set_bad(UNKNOWN)
 im=ax.pcolormesh(np.deg2rad(xx),np.deg2rad(yy),np.ma.masked_invalid(values[:,order]),cmap=cm,norm=norm,shading='flat',rasterized=True,antialiased=False)
 coast(ax);letter(ax,label)
 return ax,im


W=180/25.4
GC=[BLUE,GREEN,ORANGE]
GL=['Thin, 1–<25 Ma','Thin, 25–<75 Ma','Thick, 1–<25 Ma']
GKEY=['thin_1_25Ma','thin_25_75Ma','thick_1_25Ma']
OPS=['direct','loss_100yr','loss_30yr'];OPN=['Complete retention','100-year loss','30-year loss']
ALLOC=['geological_power','geological_area'];ALN=['Geological power','Geological area']
st=pd.read_csv(S/'SO240_stations_work_areas.csv');comp=read('external_validation/regional_comparison.csv')
phys=pd.read_csv(S/'physical_schematic_conditions.csv').set_index('parameter').value
ab=read('evidence_closure/constraint_ablation.csv');sup=read('physical_limits/support_attribution.csv')
cur=xr.load_dataset(R/'geological_support/conditional_budget/geology_conditioned_source_fields.nc');area=cur.area_m2.values
unit=lambda e:np.divide(e*1e12,area,out=np.full(area.shape,np.nan),where=area>0)
div=LinearSegmentedColormap.from_list('signed_heat',[BLUE,'#f0efec',ORANGE])
obs=read('F2_observation_summary.csv');force=read('F2_forcing_summary.csv');force=force[force.component=='A']
den=float((ab.outer_ZJ/(ab.percent_fixed_mean/100)).median())
assert np.allclose(ab.percent_fixed_mean,100*ab.outer_ZJ/den)

MODELS={'GISS-E2-1-G':BLUE,'IPSL-CM6A-LR':ORANGE}
exec(compile((Path(__file__).parent/'colour_panels_publication.py').read_text(),'colour_panels_publication.py','exec'))


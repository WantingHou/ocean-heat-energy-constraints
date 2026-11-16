"""Display-only Figure 5 revision from frozen full-precision inputs.
Computes only algebraic display curves, state masks and cap intersections.
"""
from pathlib import Path
import argparse, json, hashlib, html, io, xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
ap=argparse.ArgumentParser();ap.add_argument('--source-data',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args()
S=args.source_data;O=args.output_dir
if O.exists():raise FileExistsError(O)
O.mkdir(parents=True)
read=lambda n:pd.read_csv(S/n,float_precision='round_trip')
front=read('F5a_saved_frontiers.csv');crossings=read('F5a_saved_intersections.csv');req=read('F5b_necessary_additional_power.csv');caps=read('F5b_additional_power_caps.csv');gain=read('response_summaries.csv').set_index('loss');ab=read('constraint_ablation.csv')
den=float((ab.outer_ZJ/(ab.percent_fixed_mean/100)).median())
assert np.allclose(ab.percent_fixed_mean,100*ab.outer_ZJ/den,rtol=1e-12)
assert np.allclose(req.target_ZJ,req.target_percent_fixed_mean_O*den/100,rtol=1e-12)
assert set(req.allocation)=={'geological_power'}
cap=caps[caps.allocation=='geological_power'].set_index('release_condition').additional_anomaly_cap_TW
OPS=['direct','loss_100yr','loss_30yr'];NAMES=['Complete retention','100-year loss','30-year loss'];COLORS=['#2a78d6','#1baf7a','#eb6834'];STYLES=['-','--',':'];ALLOC=['geological_power','geological_area']
INK='#263440';GRAY='#65727b';BOUND='#424a50';PALEBLUE='#dceaf5';HATCH='#939ca3'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.2,'axes.labelsize':7.2,'xtick.labelsize':6.8,'ytick.labelsize':6.8,'legend.fontsize':6.9,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'axes.edgecolor':'#a7b0b7','axes.linewidth':.5,'pdf.fonttype':42,'svg.fonttype':'none','svg.hashsalt':'Figure5-states-caps-20261109','figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','hatch.linewidth':.35})
W=180/25.4;H=8.6;fig=plt.figure(figsize=(W,H));titles={};numeric_labels=[];records=[];bounds=[];curves=[];intersections=[];saved=[];a_axes=[]
def tag(a,gid,text):a.set_gid(gid);titles[gid]=text;return a
def clean(ax):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(length=2.5,width=.5);ax.set_axisbelow(True)
def label(ax,text,xy,offset=(0,6),ha='center',leader=False,fontsize=7):
 kw={'arrowprops':dict(arrowstyle='-',color=GRAY,lw=.45,shrinkA=2,shrinkB=3)} if leader else {}
 t=ax.annotate(text,xy,xytext=offset,textcoords='offset points',ha=ha,va='bottom' if offset[1]>=0 else 'top',fontsize=fontsize,color=INK,bbox=dict(fc='white',ec='none',pad=.35,alpha=.96),**kw)
 numeric_labels.append(t);return t
red=LinearSegmentedColormap.from_list('positive_gap_red',['#fee1dd','#de7068','#a51f32']);norm=Normalize(0,26)
rgb=red(np.linspace(0,1,256))[:,:3];linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4);lum=linear@np.array([.2126,.7152,.0722]);assert np.all(np.diff(lum)<0)
fig.text(.04,8.28/H,'a',fontsize=10,weight='bold');fig.text(.09,8.28/H,'Exposure–target bounds · 1955–2023, 0–6000 m',fontsize=8.3)
fig.text(.09,8.05/H,'3 MW actual / unit · monthly exposure bounds',fontsize=7.0)
for i,(op,name) in enumerate(zip(OPS,NAMES)):
 ax=fig.add_axes([.11+i*.295,5.85/H,.255,1.85/H]);a_axes.append(ax)
 qs=[front[(front.loss==op)&(front.allocation==al)].sort_values('equivalent_full_power_years') for al in ALLOC]
 x=qs[0].equivalent_full_power_years.to_numpy();raw=np.stack([q.outer_response_ZJ.to_numpy() for q in qs]);assert len(x)==277 and x[0]==0 and x[-1]==69 and np.array_equal(x,qs[1].equivalent_full_power_years)
 lo=raw.min(axis=0);hi=raw.max(axis=0);targets=np.linspace(0,5,101);gap=targets[:,None]*den/100-hi[None,:]
 ax.pcolormesh(x,targets,np.maximum(gap,0),cmap=red,norm=norm,shading='gouraud',rasterized=True,zorder=0)
 tag(ax.fill_between(x,0,100*lo/den,fc=PALEBLUE,ec='none',zorder=1),'below-both-'+op,'Uniform pale blue: below both upper bounds; not excluded by either, not demonstrated achievable and not cooling.')
 tag(ax.fill_between(x,100*lo/den,100*hi/den,fc='white',ec=HATCH,hatch='//',lw=0,zorder=2),'between-'+op,'Between upper bounds; white with sparse gray hatching; not a confidence interval.')
 for j,(al,q) in enumerate(zip(ALLOC,qs)):
  line=ax.plot(x,100*raw[j]/den,c=BOUND,ls='-' if j==0 else '--',lw=1.15,marker='o' if j==0 else 's',markevery=[92,184,276],ms=3.4,mfc='white',mec=BOUND,mew=.7,zorder=4)[0]
  assert np.array_equal(line.get_xdata(),q.equivalent_full_power_years.to_numpy()) and np.array_equal(line.get_ydata(),100*q.outer_response_ZJ.to_numpy()/den)
  tag(line,f'frontier-{op}-{al}',f'Original saved frontier: {op}; {al}; all 277 coordinates unchanged; zero-exposure bound {raw[j,0]!r} ZJ')
  c=crossings[(crossings.loss==op)&(crossings.allocation==al)&(crossings.target_percent_fixed_mean_O==1)].iloc[0]
  if np.isfinite(c.minimum_equivalent_full_power_years):
   xc=float(c.minimum_equivalent_full_power_years);ax.plot(xc,1,marker='o' if j==0 else 's',ms=4.2,mfc='white',mec=BOUND,mew=.8,ls='none',zorder=6)
   label(ax,f'{xc:.1f} y',(xc,1),(0,19) if j==0 else (-3,-18),ha='center',leader=True,fontsize=6.4)
  bounds.extend({'loss':op,'allocation':al,'equivalent_years':float(xx),'upper_B_ZJ':float(yy),'target_percent':float(100*yy/den)} for xx,yy in zip(x,raw[j]))
 for t in [1,5]:ax.axhline(t,c=GRAY,ls=':',lw=.55,zorder=3)
 ax.axvline(69,c=GRAY,lw=.7,zorder=5)
 if op=='loss_30yr':ax.text(.52,.255,'1% not reached',transform=ax.transAxes,ha='center',fontsize=6.6,bbox=dict(fc='white',ec='none',pad=.5,alpha=.9))
 ax.set(xlim=(0,69),ylim=(0,5),xticks=[0,30,60,69],yticks=[0,1,3,5],xlabel='Equivalent full-power years');clean(ax);ax.set_title(name,fontsize=7.6,pad=7);ax.set_xlabel('Equivalent full-power years',fontsize=6.8)
 if i==0:ax.set_ylabel('Target / fixed mean O (%)',fontsize=7.1)
 else:ax.set_yticklabels([])
 for k,xx in enumerate(x):
  for m,t in enumerate(targets):records.append({'loss':op,'equivalent_years':float(xx),'target_percent':float(t),'smaller_upper_ZJ':float(lo[k]),'larger_upper_ZJ':float(hi[k]),'positive_gap_ZJ':float(max(gap[m,k],0)),'state':'below_both' if t*den/100<lo[k] else 'above_both' if t*den/100>hi[k] else 'between_or_on_bounds'})
fig.legend(handles=[Line2D([],[],c=BOUND,ls='-',marker='o',ms=4,mfc='white',label='Power allocation'),Line2D([],[],c=BOUND,ls='--',marker='s',ms=4,mfc='white',label='Area allocation')],loc='center',bbox_to_anchor=(.53,5.25/H),ncol=2,frameon=False,handlelength=2.5,columnspacing=2.5)
fig.legend(handles=[Patch(fc=PALEBLUE,ec='#9badba',lw=.4,label='Below both upper bounds\nNot excluded by either'),Patch(fc='white',ec=HATCH,hatch='//',lw=.4,label='Between upper bounds\nAllocation-dependent exclusion'),Patch(fc='#d56865',ec='none',label='Above both upper bounds\nExcluded by both upper bounds')],loc='center',bbox_to_anchor=(.535,4.91/H),ncol=3,frameon=False,fontsize=6.5,handlelength=1.8,handleheight=1.4,columnspacing=1.6)
cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=red),cax=fig.add_axes([.26,4.49/H,.53,.10/H]),orientation='horizontal');cb.set_ticks([0,5,10,15,20,26]);cb.set_label('Red region only: target energy − larger upper bound (ZJ)',fontsize=7);cb.ax.tick_params(labelsize=6.5,length=2);cb.outline.set_visible(False)
fig.text(.04,3.96/H,'b',fontsize=10,weight='bold');fig.text(.09,3.96/H,'Necessary additional peak power · power allocation',fontsize=8.3)
fig.legend(handles=[Line2D([],[],c=co,ls=ls,lw=1.3,label=name) for name,co,ls in zip(NAMES,COLORS,STYLES)],loc='center',bbox_to_anchor=(.53,3.72/H),ncol=3,frameon=False,handlelength=2.5,columnspacing=2.3)
fig.legend(handles=[Line2D([],[],c=INK,marker='o',ls='none',ms=4,label='Saved target value'),Line2D([],[],c=INK,marker='D',mfc='white',ls='none',ms=4,label='Calculated cap intersection')],loc='center',bbox_to_anchor=(.53,3.47/H),ncol=2,frameon=False,columnspacing=2.5)
ax=fig.add_axes([.11,.66/H,.85,2.66/H]);inset=fig.add_axes([.165,2.02/H,.285,1.04/H],zorder=20)
ax.set(xlim=(0,20),ylim=(0,125),xticks=[0,5,10,15,20],yticks=[0,25,50,75,100,125],xlabel='Target / fixed mean O (%)',ylabel='Necessary extra peak power (TW)');clean(ax)
ax.grid(axis='y',c='#e6e9ec',lw=.4)
for t in [5,10]:ax.axvline(t,c='#abb4bb',lw=.55,zorder=0)
for condition in ['3 MW','14.4 MW']:
 v=float(cap.loc[condition]);ax.axhline(v,c=GRAY,ls='-.' if condition=='3 MW' else '--',lw=.9,zorder=1)
 ax.annotate(f'{v:.2f} TW · {condition} actual / unit',(19.8,v),xytext=(0,4 if condition=='3 MW' else -5),textcoords='offset points',ha='right',va='bottom' if condition=='3 MW' else 'top',fontsize=6.8,color=INK,bbox=dict(fc='white',ec='none',pad=.3))
errors=[]
for op,name,co,ls in zip(OPS,NAMES,COLORS,STYLES):
 g=float(gain.loc[op,'epoch_constant_gain_ZJ_per_TW']);q=req[req.loss==op]
 errors.extend(abs(q.target_ZJ/g-q.epoch_necessary_anomaly_peak_TW).tolist())
 assert np.allclose(q.target_ZJ/g,q.epoch_necessary_anomaly_peak_TW,rtol=1e-12)
 for a,lim in [(ax,20),(inset,2.6)]:
  xx=np.linspace(0,lim,201);yy=xx*den/100/g;a.plot(xx,yy,c=co,ls=ls,lw=1.3,zorder=3)
  if a is ax:curves.extend({'loss':op,'target_percent':float(x),'necessary_extra_peak_TW':float(y)} for x,y in zip(xx,yy))
 for row in q.itertuples():
  t=float(row.target_percent_fixed_mean_O);v=float(row.epoch_necessary_anomaly_peak_TW)
  tag(ax.plot(t,v,marker='o',ms=4.3,mfc=co,mec='white',mew=.5,ls='none',zorder=5,clip_on=False)[0],f'saved-{op}-{t:g}',f'Saved target requirement, not an observation: target {t!r}%; necessary extra peak {v!r} TW')
  saved.append({'loss':op,'target_percent':t,'necessary_extra_peak_TW':v})
  if t in [5,10]:
   offsets={(5,'direct'):(6,-2),(5,'loss_100yr'):(-5,5),(5,'loss_30yr'):(-5,11),(10,'direct'):(5,-9),(10,'loss_100yr'):(-7,8),(10,'loss_30yr'):(6,5)}
   label(ax,f'{v:.2f}',(t,v),offsets[(t,op)],ha='right' if offsets[(t,op)][0]<0 else 'left',leader=True)
 for condition in ['3 MW','14.4 MW']:
  v=float(cap.loc[condition]);t=100*v*g/den;assert np.isclose(t*den/100/g,v,rtol=1e-14)
  destination=inset if condition=='3 MW' else ax
  tag(destination.plot(t,v,marker='D',ms=5.5,mfc='white',mec=co,mew=1.05,ls='none',zorder=6)[0],f'cap-intersection-{op}-{condition.replace(" ","-")}',f'Calculated intersection, not an observation: {op}; {condition} actual unit cap; target {t!r}%; extra power {v!r} TW')
  offset=(0,(8.65-v)/8*1.04*72) if condition=='3 MW' else (11,13)
  label(destination,f'{t:.2f}%',(t,v),offset,ha='center' if condition=='3 MW' else 'left',leader=True,fontsize=6.8)
  intersections.append({'loss':op,'actual_unit_condition':condition,'gain_ZJ_per_TW':g,'fixed_mean_O_ZJ':den,'additional_power_cap_TW':v,'target_percent':t,'display_panel':'inset' if condition=='3 MW' else 'main','kind':'algebraic cap intersection; not an observation'})
inset.axhline(float(cap.loc['3 MW']),c=GRAY,ls='-.',lw=.8,zorder=1)
inset.set(xlim=(0,2.6),ylim=(0,8),xticks=[0,1,2],yticks=[0,4,8],xlabel='Target (%)',ylabel='Extra peak (TW)');clean(inset);inset.tick_params(labelsize=6);inset.set_xlabel('Target (%)',fontsize=6.2,labelpad=2);inset.set_ylabel('Extra peak (TW)',fontsize=6.2,labelpad=2);inset.set_title('3 MW cap intersections · 5.07 TW',fontsize=6.6,pad=5)
fig.text(.11,.17/H,'Necessary, not sufficient: being below a cap does not establish a feasible source history.',fontsize=7)
fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
for t in fig.findobj(matplotlib.text.Text):
 if t.get_visible() and t.get_text():
  bb=t.get_window_extent(renderer)
  if bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
collisions=[]
for i,t in enumerate(numeric_labels):
 for u in numeric_labels[i+1:]:
  if t.axes is u.axes and t.get_bbox_patch().get_window_extent(renderer).overlaps(u.get_bbox_patch().get_window_extent(renderer)):collisions.append([t.get_text(),u.get_text()])
assert not collisions,collisions
assert len(bounds)==1662 and len(saved)==12 and len(intersections)==6
for rows,n in [(bounds,'F5a_plotted_frontiers.csv'),(records,'F5a_state_and_gap_display.csv'),(saved,'F5b_saved_target_values.csv'),(curves,'F5b_fixed_gain_curves.csv'),(intersections,'F5b_cap_intersections.csv')]:pd.DataFrame(rows).to_csv(O/n,index=False,mode='x')
name='Figure_5_power_persistence'
for ext in ['pdf','png']:
 with (O/(name+'.'+ext)).open('xb') as f:fig.savefig(f,format=ext,dpi=600)
stream=io.StringIO();fig.savefig(stream,format='svg',dpi=600);svg=stream.getvalue()
for gid,title in titles.items():svg=svg.replace('<g id="'+gid+'">','<g id="'+gid+'"><title>'+html.escape(title)+'</title>')
root=ET.fromstring(svg);ids=[e.attrib['id'] for e in root.iter() if 'id' in e.attrib];assert len(ids)==len(set(ids));assert any(e.tag.endswith('}text') for e in root.iter())
with (O/(name+'.svg')).open('x') as f:f.write(svg)
with (O/(name+'_preview.png')).open('xb') as f:fig.savefig(f,format='png',dpi=180)
fig.set_size_inches(W*170/180,H*170/180)
for t in fig.findobj(matplotlib.text.Text):t.set_fontsize(t.get_fontsize()*170/180)
with (O/'Figure_5_print_size_170mm.pdf').open('xb') as f:fig.savefig(f,format='pdf',dpi=600)
checks={'fixed_mean_O_ZJ':den,'window':'1955–2023','depth_m':'0–6000','original_frontier_values':len(bounds),'frontier_coordinates_unchanged':True,'zero_year_offsets_preserved':front[front.equivalent_full_power_years==0][['allocation','loss','outer_response_ZJ']].to_dict('records'),'original_one_percent_intersections':crossings[(crossings.target_percent_fixed_mean_O==1)&crossings.minimum_equivalent_full_power_years.notna()][['allocation','loss','minimum_equivalent_full_power_years']].to_dict('records'),'saved_target_values':len(saved),'cap_intersections':intersections,'saved_values_vs_gain_max_abs_difference_TW':max(errors),'red_ramp_luminance_strictly_decreasing':True,'uniform_blue_not_a_magnitude_scale':True,'state_legend_separate_from_positive_red_scale':True,'text_outside_canvas':outside,'numeric_label_collisions':collisions,'svg_text_editable':True,'svg_ids_unique':True,'inset_count':1,'png_dpi':600,'width_mm':180,'print_check_width_mm':170,'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in S.iterdir() if p.is_file()},'scientific_changes':False}
with (O/'figure_checks.json').open('x') as f:json.dump(checks,f,indent=2,ensure_ascii=False)
print(json.dumps(checks,indent=2,ensure_ascii=False))

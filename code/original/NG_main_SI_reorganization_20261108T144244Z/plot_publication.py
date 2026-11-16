"""Offline display reorganization using frozen source tables only.
No candidate screening, source optimization or model fitting is performed.
"""
from pathlib import Path
exec(compile((Path(__file__).parent/'shared_plotting.py').read_text(),'shared_plotting.py','exec'))
from matplotlib.transforms import Bbox, blended_transform_factory
from matplotlib.patches import ConnectionPatch
RC=dict(zip(OPS,[BLUE,GREEN,ORANGE]));LS=dict(zip(OPS,['-','--',':']))
def table(name):return pd.read_csv(S/'panels'/name,float_precision='round_trip')
def csv(frame,name):frame.to_csv(O/name,index=False,mode='x')
def tagged(artist,gid,title):artist.set_gid(gid);CELL_TITLES[gid]=title;return artist
def heading(fig,y,letter,title):fig.text(.065,y,letter,fontsize=10,weight='bold');fig.text(.11,y,title,fontsize=8)
class PanelCanvas:
 def __init__(self,canvas,height,bottom=0,target=0):self.canvas=canvas;self.height=height;self.bottom=bottom;self.target=target
 def y(self,y):return ((y-self.bottom)*self.height+self.target)/self.canvas.get_figheight()
 def add_axes(self,rect,**kw):
  x,y,w,h=rect;return self.canvas.add_axes([x,self.y(y),w,h*self.height/self.canvas.get_figheight()],**kw)
 def text(self,x,y,*args,**kw):return self.canvas.text(x,self.y(y),*args,**kw)
 def legend(self,*args,**kw):
  if 'bbox_to_anchor' in kw:
   x,y=kw['bbox_to_anchor'];kw['bbox_to_anchor']=(x,self.y(y))
  return self.canvas.legend(*args,**kw)
 def colorbar(self,*args,**kw):return self.canvas.colorbar(*args,**kw)
fragment=json.loads((Path(__file__).parent/'layout_fragments.json').read_text())
fragment['S7_body']=fragment['S7_body'].replace("f'{v:.1f}%'","f'{v:.1f}% outside'")

# Main Figure 2: retained regional observations, then original-allocation area test.
exec(compile(fragment['F2_upper'],'retained_F2ab','exec'),globals());fig=canvas
heading(fig,.489,'c','Same-basin collection-area constraints')
audit=read('geological_support/reference_geometry_audit.csv');audit['R_required_over_eligible']=audit.required_collection_million_km2/audit.known_eligible_flank_million_km2;audit['allocation']='retained'
assert np.allclose(audit.R_required_over_eligible*audit.known_area_coverage_ratio,1,rtol=1e-12)
assert np.isfinite(audit.R_required_over_eligible).all()
basins=read('basin_geometry_summary.csv').basin.tolist();assert basins==['Pacific','Atlantic','Indian','Southern','Arctic','Mediterranean']
central=audit[audit.reference=='central'].set_index('basin');geo_audit=read('basin_mask_assignment_audit.csv');assert set(geo_audit.basin_index)==set(range(6))
# Aggregate only saved cell geometry, never the response fields of an older branch.
lon=((geo_audit.longitude.to_numpy()+180)%360)-180;lat=geo_audit.latitude.to_numpy();ix=np.clip(((lon+180)//2).astype(int),0,179);iy=np.clip(((lat+90)//2).astype(int),0,89);bins=np.zeros((6,90,180));np.add.at(bins,(geo_audit.basin_index.to_numpy(int),iy,ix),geo_audit.area_m2.to_numpy());wet=bins.sum(axis=0)>0;dominant=bins.argmax(axis=0)
ratios=np.array([central.loc[b,'R_required_over_eligible'] for b in basins]);mapped=np.where(wet,ratios[dominant],np.nan)
cn=BoundaryNorm([0,1,2,4,8],4);cm=ListedColormap([GREEN,'#f9d2b7','#ef9a60','#c75b24'])
ax,im=mapplot(fig,[.055,.267,.49,.18],mapped,np.arange(-179,180,2),np.arange(-89,90,2),cm,cn,'')
fig.text(.07,.458,'Original allocation · central ordinary reference',fontsize=6.7)
for b,lo,la in [('Pacific',-150,-3),('Atlantic',-30,10),('Indian',75,-25),('Southern',20,-65),('Arctic',0,76)]:
 ax.text(np.deg2rad(lo%360-180),np.deg2rad(la),b,ha='center',fontsize=5.7,bbox=dict(fc='white',ec='none',alpha=.8,pad=.25),zorder=10)
cb=fig.colorbar(im,cax=fig.add_axes([.10,.225,.37,.010]),orientation='horizontal',spacing='proportional');cb.set_ticks([0,1,2,4,8]);cb.ax.tick_params(labelsize=6,length=2);cb.set_label('R = required / eligible same-basin area',fontsize=6.2);cb.outline.set_visible(False)
fig.text(.07,.179,'Basin-level support only; not modern warming.',fontsize=6.0)
ax=fig.add_axes([.735,.265,.215,.177]);range_rows=[]
for i,b in enumerate(basins):
 q=audit[audit.basin==b].set_index('reference');vals=q.loc[['low','central','high'],'R_required_over_eligible'].to_numpy();ax.plot([vals[0],vals[-1]],[i,i],c='#b5bec4',lw=1)
 for value,marker,ref in zip(vals,['o','D','^'],['low','central','high']):
  tagged(ax.scatter(value,i,s=21,c=[cm(cn(value))],marker=marker,edgecolors=INK,lw=.35,zorder=3),'R-'+b+'-'+ref,f'Original allocation; {b}; {ref}; R={value}')
  range_rows.append({'allocation':'retained','basin':b,'reference':ref,'R':value})
ax.axvline(1,c=INK,ls='--',lw=.75);ax.set(xlim=(0,8),xticks=[0,1,4,8],ylim=(5.6,-.6),yticks=range(6),yticklabels=basins,xlabel='Required / eligible area');clean(ax);ax.tick_params(labelsize=6.0);ax.set_xlabel('Required / eligible area',fontsize=6.1)
fig.text(.61,.458,'Original allocation · reference range',fontsize=6.4)
fig.legend(handles=[Line2D([],[],ls='',marker=m,mfc='white',mec=INK,ms=3.5,label=r) for m,r in [('o','Low'),('D','Central'),('^','High')]],loc='center',bbox_to_anchor=(.775,.201),ncol=3,frameon=False,fontsize=6.0,columnspacing=.65,handlelength=.65)
feas=pd.read_csv(S/'allocation_collection_feasibility.csv',float_precision='round_trip');f=feas.set_index('allocation').loc[['retained']+ALLOC]
ax=fig.add_axes([.31,.067,.64,.091]);low=float(f.ordinary_F_low_TW.iloc[0]);high=float(f.ordinary_F_high_TW.iloc[0]);ax.axvspan(low,high,fc='#e2e8ed',zorder=0)
for i,(al,r) in enumerate(f.iterrows()):
 ax.plot([0,r.Fcap_TW],[i,i],c=GRAY,lw=1.3);ax.plot(r.Fcap_TW,i,marker='|',ms=7,c=INK)
 if r.ordinary_reference_nonempty:ax.plot([r.admitted_reference_low_TW,r.admitted_reference_high_TW],[i,i],c=GREEN,lw=3,solid_capstyle='butt')
 else:assert np.isnan(r.admitted_reference_low_TW);ax.text(2.4,i,'No intersection',va='center',ha='center',fontsize=6.1)
ax.set(xlim=(0,8.5),xticks=[0,2,4,6,8],ylim=(2.6,-.6),yticks=[0,1,2],yticklabels=['Original','Geological power','Geological area'],xlabel='Ordinary flank power, F (TW)');clean(ax);ax.tick_params(labelsize=6.1);ax.set_xlabel('Ordinary flank power, F (TW)',fontsize=6.4)
fig.text(.065,.164,'Allocation-specific power intersection',fontsize=6.6)
fig.legend(handles=[Patch(fc='#e2e8ed',ec='none',label='Ordinary reference band'),Line2D([],[],c=GREEN,lw=3,label='Admitted intersection')],loc='center',bbox_to_anchor=(.70,.176),ncol=2,frameon=False,fontsize=5.9,columnspacing=.8,handlelength=1.2)
save(fig,'Figure_2_paired_budget');csv(audit,'F2c_basin_area_ratios.csv');csv(pd.DataFrame(range_rows),'F2c_reference_ranges.csv');csv(f.reset_index(),'F2c_allocation_power_intersections.csv')
yy,xx=np.where(wet);csv(pd.DataFrame({'latitude':-89+2*yy,'longitude':-179+2*xx,'basin_index':dominant[wet],'R':mapped[wet]}),'F2c_basin_display_grid.csv')

# Main Figure 3: original points, adjusted label, frozen model-deletion support.
exec(compile(fragment['F3_abc'],'retained_F3abc','exec'),globals())
fig.text(.61,.443,'d',fontsize=10,weight='bold');fig.text(.65,.443,'Model dependence',fontsize=7.7);fig.text(.61,.409,'Power allocation · 3 MW · complete retention',fontsize=6.2)
paired=pd.read_csv(S/'paired_current_joint_support.csv',float_precision='round_trip');primary=paired[paired.allocation==ALLOC[0]];assert len(primary)==4 and set(primary.loss)=={'direct'} and set(primary.family)=={'S2_P3'}
deldata=read('physical_limits/deletion_diagnostics.csv');groups=[('All retained\nhistories',primary,None),('Without GISS',primary[primary.model!='GISS-E2-1-G'],'GISS-E2-1-G'),('Without IPSL',primary[primary.model!='IPSL-CM6A-LR'],'IPSL-CM6A-LR')]
ax=fig.add_axes([.755,.19,.195,.18]);support=[]
for i,(label,q,deleted) in enumerate(groups):
 lo=q.B_over_O_percent.min();hi=q.B_over_O_percent.max();assert len(q)==[4,3,1][i]
 if deleted:
  check=deldata[(deldata.allocation==ALLOC[0])&(deldata.family=='S2_P3')&(deldata.loss=='direct')&(deldata.deleted_field=='model')&(deldata.deleted_value==deleted)].iloc[0]
  assert int(check.count_remaining)==len(q) and np.allclose([lo,hi],[check.lower_percent,check.upper_percent])
 if len(q)>1:ax.plot([lo,hi],[i,i],c='#ccd4da',lw=1.6,zorder=1)
 for r in q.itertuples():
  tagged(ax.scatter(r.B_over_O_percent,i,s=23,c=MODELS[r.model],edgecolors='white',lw=.5,zorder=3),'deletion-'+str(i)+'-'+str(r.candidate_index),f'{label}; {r.model}; candidate {r.candidate_index}; B/O={r.B_over_O_percent}%')
  support.append(dict(group=i,label=label.replace('\n',' '),candidate_index=r.candidate_index,model=r.model,paired_O_ZJ=r.O_1955_2023_ZJ,B_over_O_percent=r.B_over_O_percent))
 if i<2:ax.annotate(f'{hi:.2f}%',(hi,i),xytext=(-1,7),textcoords='offset points',ha='right',fontsize=6.2)
 ax.text(1.10,i,str(len(q)),transform=ax.get_yaxis_transform(),ha='center',va='center',fontsize=6.6)
ax.axvline(0,c=GRAY,lw=.6);ax.set(xlim=(-1.48,2.0),xticks=[-1,0,1,2],ylim=(2.55,-.55),yticks=range(3),yticklabels=[g[0] for g in groups],xlabel='Paired B/O (%)');clean(ax);ax.tick_params(labelsize=6.2);ax.set_xlabel('Paired B/O (%)',fontsize=6.5);ax.text(1.10,1.07,'n',transform=ax.transAxes,ha='center',fontsize=6.4)
fig.text(.61,.120,'Upper endpoint: 1.79% → 0.42% without GISS',fontsize=6.1);fig.text(.61,.075,'Finite-support deletion sensitivity,\nnot independent prediction validation.',fontsize=6.0,linespacing=1.25)
save(fig,'Figure_3_energy_compatibility');csv(pd.DataFrame(labels),'F3c_paired_label_audit.csv');csv(pd.DataFrame(support),'F3d_model_dependence_points.csv')

# Main Figure 4: restored full native schematic and focused source upper bounds.
exec(compile((Path(__file__).parent/'ppt_full.py').read_text(),'ppt_full.py','exec'),globals())
fig=plt.figure(figsize=(W,10.2));a=fig.add_axes([0,5.20/10.2,1,(W/1.5*1.04)/10.2]);ppt_audit=draw_supplied_ppt(a,Path(__file__).parents[1]/'assets/Figure4a_supplied_editable.pptx');a.set_ylim(9753600*1.04,0)
heading(fig,4.84/10.2,'b','Release-condition upper bounds')
fig.legend(handles=[Patch(fc='#2565a8',ec='none',label='3 MW actual / unit'),Patch(fc='#e3edf6',ec='#2565a8',lw=.8,label='14.4 MW expanded')],loc='center',bbox_to_anchor=(.72,4.88/10.2),ncol=2,frameon=False,fontsize=6.0,handlelength=1.2)
pairs=table('F4b_controlled_pairs.csv');shown=pairs[pairs.allocation.isin(ALLOC)];assert len(shown)==6
for i,(op,title) in enumerate(zip(OPS,OPN)):
 ax=fig.add_axes([.14+i*.29,3.65/10.2,.23,.80/10.2]);z=fig.add_axes([.14+i*.29,2.65/10.2,.23,.55/10.2]);q=shown[shown.loss==op].set_index('allocation')
 for j,al in enumerate(ALLOC):
  r=q.loc[al]
  ax.barh(j-.16,r.restricted_outer_ZJ,height=.25,color='#2565a8',ec='white',lw=.4)
  ax.barh(j+.16,r.broad_outer_ZJ,height=.25,fc='#e3edf6',ec='#2565a8',lw=.7)
  if i==0:ax.text(r.broad_outer_ZJ+1,j+.16,f'{r.broad_outer_ZJ:.2f}',fontsize=5.9,va='center')
  z.barh(j,r.restricted_outer_ZJ,height=.35,color='#2565a8');z.text(r.restricted_outer_ZJ+.3,j,f'{r.restricted_outer_ZJ:.2f}',fontsize=6.0,va='center')
 ax.set(xlim=(0,100),xticks=[0,50,100],ylim=(1.6,-.6),yticks=[0,1],yticklabels=['Power','Area'] if i==0 else ['','']);clean(ax);ax.set_title(title,fontsize=7,pad=6);ax.tick_params(labelsize=6.0)
 z.set(xlim=(0,15),xticks=[0,5,10,15],ylim=(1.6,-.6),yticks=[0,1],yticklabels=['Power','Area'] if i==0 else ['',''],xlabel='Upper B (ZJ)');clean(z);z.tick_params(labelsize=5.9);z.set_title('Low-value zoom · 3 MW only',fontsize=6,pad=5);z.set_xlabel('Upper B (ZJ)',fontsize=6.3)
heading(fig,2.20/10.2,'c','Three-stage tightening · power allocation / 3 MW')
stock=read('geological_support/conditional_budget/source_only_bounds.csv');stock=stock[(stock.allocation==ALLOC[0])&(stock.family=='S2_P3')&(stock.loss=='direct')&(stock.first_year==1955)&(stock.last_year==2023)&(stock.endpoint=='high')];assert len(stock)==1;stock=stock.iloc[0]
fig.text(.14,2.02/10.2,f'Direct source-upper witness: {stock.old_remaining_2023_ZJ:.2f} ZJ old stock remains; not a claim about all histories.',fontsize=6.0)
stages=table('F4c_nested_bounds.csv');stage_rows=[]
for i,(op,title) in enumerate(zip(OPS,OPN)):
 ax=fig.add_axes([.14+i*.29,.58/10.2,.23,1.10/10.2]);r=stages[stages.loss==op].iloc[0];vals=r[['monthly_outer_ZJ','epoch_source_upper_ZJ','joint_upper_ZJ']].to_numpy(float);assert np.all(np.diff(vals)<1e-10)
 if op=='direct':assert abs(vals[0]-vals[1])<1e-10 and np.isclose(stock.B_ZJ,vals[1])
 ax.stairs(vals,[-.5,.5,1.5,2.5],fill=True,color=RC[op],alpha=.10);ax.stairs(vals,[-.5,.5,1.5,2.5],color=RC[op],lw=1.3);ax.plot(range(3),vals,'o',c=RC[op],ms=3,mec='white',mew=.4)
 for j in [0,2]:ax.text(j,vals[j]+.6,f'{vals[j]:.2f}',ha='center',fontsize=6.1)
 ax.set(xlim=(-.5,2.5),ylim=(0,13),yticks=[0,4,8,12],xticks=[0,1,2],xticklabels=['Timing\nrelaxed','Complete\nsource','Joint\nbudget']);clean(ax);ax.tick_params(labelsize=6);ax.set_title(title,fontsize=7,pad=6)
 if i==0:ax.set_ylabel('Upper B (ZJ)',fontsize=6.4)
 else:ax.set_yticklabels([])
 for stage,value in zip(['timing_relaxed','complete_source','joint_budget'],vals):stage_rows.append({'loss':op,'stage':stage,'upper_B_ZJ':value})
fig.text(.065,.13/10.2,'Stage differences tighten constraints; they do not measure physical heat loss.',fontsize=6.1)
for text in fig.findobj(matplotlib.text.Text):
 if text.get_visible() and text.get_text() and text.get_fontsize()<6.8:text.set_fontsize(6.8)
save(fig,'Figure_4_release_controls');csv(shown,'F4b_primary_allocation_bounds.csv');csv(pd.DataFrame(stage_rows),'F4c_stage_display.csv');csv(stock.to_frame().T,'F4c_direct_source_upper_witness.csv')
with (O/'F4a_native_object_audit.json').open('x') as f:json.dump(ppt_audit,f,indent=2)

# Main Figure 5: two panels, positive gaps start orange, power zoom repeats saved points.
exec(compile(fragment['F5_a'],'retained_F5a','exec'),globals());fig=canvas
heading(fig,.483,'b','Necessary additional peak power · power allocation')
req=table('F5b_necessary_additional_power.csv');caps=table('F5b_additional_power_caps.csv');gains=read('physical_limits/response_summaries.csv').set_index('loss');assert set(req.allocation)=={ALLOC[0]}
ax=fig.add_axes([.12,.12,.49,.31]);zoom=fig.add_axes([.76,.21,.18,.16]);curve=[]
for op,title in zip(OPS,OPN):
 q=req[req.loss==op];gain=float(gains.loc[op,'epoch_constant_gain_ZJ_per_TW']);assert np.allclose(q.target_ZJ/gain,q.epoch_necessary_anomaly_peak_TW,rtol=1e-12)
 for a,limit in [(ax,20),(zoom,2.5)]:
  t=np.linspace(0,limit,201);power=t/100*den/gain;a.plot(t,power,c=RC[op],ls=LS[op],lw=1.2,label=title)
  visible=q[q.target_percent_fixed_mean_O<=limit];a.scatter(visible.target_percent_fixed_mean_O,visible.epoch_necessary_anomaly_peak_TW,s=15,c=RC[op],edgecolors='white',lw=.3,zorder=4)
  if a is ax:curve.extend({'loss':op,'target_percent':x,'necessary_extra_peak_TW':y} for x,y in zip(t,power))
for condition,style in [('3 MW','-.'),('14.4 MW','--')]:
 cap=float(caps[(caps.allocation==ALLOC[0])&(caps.release_condition==condition)].additional_anomaly_cap_TW.iloc[0]);ax.axhline(cap,c=GRAY,ls=style,lw=.8)
 ax.annotate(f'{cap:.2f} TW\n{condition} units',(20,cap),xytext=(5,0),textcoords='offset points',fontsize=6,va='bottom' if condition=='3 MW' else 'center',annotation_clip=False)
 if condition=='3 MW':zoom.axhline(cap,c=GRAY,ls=style,lw=.8);zoom.annotate(f'{cap:.2f}',(2.5,cap),xytext=(2,0),textcoords='offset points',va='center',fontsize=6.0,annotation_clip=False)
ax.add_patch(Rectangle((0,0),2.5,16,fill=False,ec=GRAY,ls=':',lw=.8));ax.text(2.8,17,'Zoom',fontsize=5.6)
ax.set(xlim=(0,20),ylim=(0,125),xticks=[0,5,10,20],yticks=[0,40,80,120],xlabel='Target / fixed mean O (%)',ylabel='Necessary extra peak (TW)');clean(ax);ax.tick_params(labelsize=6.3);ax.set_xlabel('Target / fixed mean O (%)',fontsize=6.4);ax.set_ylabel('Necessary extra peak (TW)',fontsize=6.4);ax.legend(loc='upper left',frameon=False,fontsize=6.1)
zoom.set(xlim=(0,2.5),ylim=(0,16),xticks=[0,1,2],yticks=[0,5,10,15],xlabel='Target (%)');clean(zoom);zoom.tick_params(labelsize=5.8);zoom.set_xlabel('Target (%)',fontsize=6.1);zoom.set_title('Low-target zoom\n0–2.5%; 0–16 TW',fontsize=6.1,pad=7)
fig.text(.065,.035,'Below a power cap is necessary, not sufficient: stock, timing and joint budgets still apply.',fontsize=6.2)
save(fig,'Figure_5_power_persistence');csv(pd.DataFrame(gaprows),'F5a_positive_gap_display.csv');csv(pd.DataFrame(curve),'F5b_algebraic_display_curves.csv');csv(gains.reset_index(),'F5b_frozen_response_gains.csv')

# S6 replaces the promoted area diagnostic with the regional paired budgets.
fig=plt.figure(figsize=(W,4.4));fig.text(.065,.945,'Regional paired heat budgets: SO240',fontsize=9,weight='bold');fig.text(.065,.88,'The same observations under three age conditions; not independent samples.',fontsize=7)
budget=table('F2c_paired_heat_flow.csv')
for i,age in enumerate([18,20,22]):
 ax=fig.add_axes([.18+i*.26,.34,.23,.40]);q=budget[budget.age_Ma==age].set_index('summary');ref=float(q.cooling_reference_mW_m2.iloc[0])
 for j,key in enumerate(['primary_station','primary_profile_equal','model']):
  r=q.loc[key];assert np.isclose(r.conduction_mW_m2+r.deficit_mW_m2,ref);observed=j<2
  ax.barh(j,r.conduction_mW_m2,height=.60,color=BLUE if observed else ORANGE,ec='none')
  ax.barh(j,r.deficit_mW_m2,left=r.conduction_mW_m2,height=.60,fc='#e6eff9' if observed else '#fff0e4',ec='#8ba6c4' if observed else '#c58f67',hatch='///' if observed else '..',lw=.3)
 ax.axvline(ref,c=INK,lw=.8);ax.text(ref,2.65,f'{ref:.1f}',ha='center',fontsize=6.2);ax.set(xlim=(0,130),xticks=[0,40,80,120],ylim=(2.9,-.7),yticks=[0,1,2],yticklabels=['Station equal','Profile equal','Model'] if i==0 else ['','',''],xlabel='Heat flux (mW m⁻²)');clean(ax);ax.set_title(f'{age} Ma',fontsize=8,pad=8);ax.tick_params(labelsize=6.3);ax.set_xlabel('Heat flux (mW m⁻²)',fontsize=6.3)
fig.legend(handles=[Patch(fc=BLUE,label='Observed conduction'),Patch(fc=ORANGE,label='Model conduction'),Patch(fc='#e6eff9',ec='#8ba6c4',hatch='///',label='Inferred deficit'),Patch(fc='#fff0e4',ec='#c58f67',hatch='..',label='Model deficit'),Line2D([],[],c=INK,lw=.8,label='Cooling reference')],loc='center',bbox_to_anchor=(.55,.16),ncol=3,frameon=False,fontsize=6.1,columnspacing=.9)
fig.text(.065,.043,'Conduction + inferred deficit = cooling reference; inferred deficit is not measured discharge.',fontsize=6.3)
save(fig,'Supplementary_Figure_S6_regional_paired_budgets')

# S7 receives the six primary marginal rows, explicitly under the expanded envelope.
fig=plt.figure(figsize=(W,4.2));fig.text(.065,.945,'Marginal requirements under the expanded 14.4 MW release envelope',fontsize=8,weight='bold');fig.text(.065,.877,'1.5× → 7.2× ordinary supply · power allocation · independent marginal minima',fontsize=6.7)
exec(compile(fragment['S7_body'],'relocated_marginal_requirements','exec'),globals())
fig.text(.065,.028,'Enhanced allowance is ordinary-power weighted, not field count or seafloor-area fraction.',fontsize=6.1)
save(fig,'Supplementary_Figure_S7_expanded_marginal_requirements');csv(pd.DataFrame(rows),'S7_expanded_primary_marginal_rows.csv')

mapping=[
 ('2a-b','SO240_stations_work_areas.csv','2015','regional observations; original station values; descriptive violins'),
 ('2c-map/ranges','retained/geological_support/reference_geometry_audit.csv;retained/basin_mask_assignment_audit.csv','ordinary reference','ORIGINAL retained allocation only; R=required/eligible; area-dominant 2-degree display bins, no smoothing'),
 ('2c-power','allocation_collection_feasibility.csv','ordinary reference','three distinct allocations; empty admitted interval remains missing, not measured zero'),
 ('3a','retained/F2_observation_summary.csv;retained/F2_forcing_summary.csv','1955–2023;0–6000m','unchanged means and descriptive ranges; no band outlines'),
 ('3b','retained/geological_support/conditional_budget/geology_conditioned_source_fields.nc;native_A_latitude_integral.csv','1970–2014;0–6000m','A spatial diagnostic only; no B field used'),
 ('3c','panels/F3c_primary_budgets.csv','1955–2023;0–6000m','five original tuples; original paired O; label position only changed'),
 ('3d','paired_current_joint_support.csv;retained/physical_limits/deletion_diagnostics.csv','1955–2023;0–6000m','power/S2_P3/direct; frozen support counts4/3/1; finite support, not CI'),
 ('4a','assets/Figure4a_supplied_editable.pptx;physical_schematic_conditions.csv','conceptual','restored optimized native schematic; not measured connectivity'),
 ('4b','panels/F4b_controlled_pairs.csv','1955–2023;0–6000m','12 values for power/area; actual unit cap3/14.4MW, not extra power; low zoom only3MW'),
 ('4c','panels/F4c_nested_bounds.csv;retained/geological_support/conditional_budget/source_only_bounds.csv','1955–2023;0–6000m','nine frozen heights; direct source-max witness stock only'),
 ('5a','panels/F5a_saved_frontiers.csv;panels/F5a_saved_intersections.csv','1955–2023;0–6000m','fixed mean O; pure white below larger bound, between-bound hatching, any positive gap starts orange'),
 ('5b','panels/F5b_necessary_additional_power.csv;panels/F5b_additional_power_caps.csv;retained/physical_limits/response_summaries.csv','1955–2023;0–6000m','frozen constant gain; geological-power extra caps; 12 saved points, inset repeats existing1% points'),
 ('S6','panels/F2c_paired_heat_flow.csv','2015;18/20/22Ma conditions','same data; paired conduction/deficit/cooling bars'),
 ('S7','panels/F5c_marginal_requirements.csv','1955–2023;0–6000m','expanded1.5×to7.2×/14.4MW envelope; two independent marginal minima; missing duration not zero')]
csv(pd.DataFrame(mapping,columns=['panel','sources','window','scope']),'panel_source_mapping.csv')
with (O/'display_contract.json').open('x') as f:json.dump({'scientific_results_changed':False,'no_model_or_optimizer_run':True,'main_window':'1955–2023;0–6000m','independent_A_window':'1970–2014;0–6000m','fixed_mean_O_ZJ_for_F5':den,'Figure3_uses_paired_O':True,'S6':'regional paired heat budgets','S7':'expanded-release marginal requirements','S8_S9':'unchanged complete conditions','F2c_ratio':'required collection / eligible same-basin area; original allocation','F5a_positive_gap_colours':['#ffe1b8','#d77a35','#8b351c'],'F5a_unexcluded_region':'pure white, not proven feasible','F5b_points':'saved values only; zoom repeats existing1% values','F4_direct_source_stock_ZJ':float(stock.old_remaining_2023_ZJ)},f,indent=2)
with (O/'figure_checks.json').open('x') as f:json.dump(checks,f,indent=2)
print(json.dumps(checks,indent=2))

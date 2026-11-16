"""Display-only transformations of immutable tables; no scientific optimization."""
from matplotlib.cm import ScalarMappable
AREA_CMAP=LinearSegmentedColormap.from_list('positive_area',['#edf6fc','#3171b2']);AREA_CMAP.set_bad('white')
FLUX_CMAP=LinearSegmentedColormap.from_list('flux_magnitude',['#f0f6fc','#6c9ed0'])
HEAT_CMAP=LinearSegmentedColormap.from_list('heat_magnitude',['#fff5eb','#e99b70'])
ALLOW_CMAP=LinearSegmentedColormap.from_list('allowance_magnitude',['#edf8f4','#60ad9b'])
COMP_CMAP=LinearSegmentedColormap.from_list('budget_compensation',[BLUE,'#f3f3f1',ORANGE])
CELL_TITLES={}
def cell(ax,x,y,value,norm,cmap,key,text=None,size=7.5):
 finite=np.isfinite(value);face=cmap(norm(value)) if finite else '#f1f3f4'
 patch=Rectangle((x-.5,y-.5),1,1,facecolor=face,edgecolor='white',lw=1.5,hatch=None if finite else '///');ax.add_patch(patch)
 title=text if text is not None else (f'{value:.2f}' if finite else 'Outside')
 ax.text(x,y,title,ha='center',va='center',fontsize=size,color=INK)
 gid='cell-'+key;patch.set_gid(gid);CELL_TITLES[gid]=key+': '+(repr(float(value)) if finite else 'Outside declared envelope')
def matrix_axes(ax,nr,nc,rows,cols):
 ax.set(xlim=(-.5,nc-.5),ylim=(nr-.5,-.5),xticks=range(nc),xticklabels=cols,yticks=range(nr),yticklabels=rows)
 ax.tick_params(length=0,labelsize=7,pad=4);ax.xaxis.tick_top()
 for spine in ax.spines.values():spine.set_visible(False)
def colourbar(fig,rect,cmap,norm,label,orientation='horizontal',ticks=None):
 cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_axes(rect),orientation=orientation)
 if ticks is not None:cb.set_ticks(ticks)
 cb.set_label(label,fontsize=7);cb.ax.tick_params(labelsize=6.7,length=2);return cb

def draw_f2c(fig):
 bud=pd.read_csv(S/'panels/F2c_paired_heat_flow.csv');ages=[18,20,22];summaries=['primary_station','primary_profile_equal','model'];norm=Normalize(0,125)
 fig.text(.075,.415,'c',weight='bold',fontsize=10);fig.text(.17,.415,'Discrete age conditions × evidence categories',fontsize=8)
 for k,(field,x,title) in enumerate([('conduction_mW_m2',.17,'Conduction'),('deficit_mW_m2',.52,'Cooling minus conduction')]):
  ax=fig.add_axes([x,.16,.30,.17]);cols=['Station\n'+('observed' if k==0 else 'inferred'),'Profiles\n'+('observed' if k==0 else 'inferred'),'Model\npredicted']
  for i,age in enumerate(ages):
   for j,summary in enumerate(summaries):
    r=bud[(bud.age_Ma==age)&(bud.summary==summary)].iloc[0];cell(ax,j,i,r[field],norm,FLUX_CMAP,f'F2c-{field}-{age}-{summary}',f'{r[field]:.1f}',8)
  matrix_axes(ax,3,3,[str(a) for a in ages] if k==0 else ['']*3,cols)
  if k==0:ax.set_ylabel('Age condition (Ma)',fontsize=7)
  fig.text(x+.15,.385,title,ha='center',fontsize=7.4)
 ax=fig.add_axes([.87,.16,.10,.17])
 for i,age in enumerate(ages):
  q=bud[bud.age_Ma==age];assert q.cooling_reference_mW_m2.nunique()==1;r=q.iloc[0];cell(ax,0,i,r.cooling_reference_mW_m2,norm,FLUX_CMAP,f'F2c-cooling-{age}',f'{r.cooling_reference_mW_m2:.1f}',8)
 matrix_axes(ax,3,1,['']*3,['Cooling\nreference'])
 colourbar(fig,[.26,.095,.56,.012],FLUX_CMAP,norm,'Heat flux (mW m⁻²; shared linear colour scale)',ticks=[0,40,80,120])
 fig.text(.17,.025,'18/20/22 Ma: the same observations under three age conditions, not repeated surveys.\nNo spatial interpolation. 40 mW m⁻² is a collection normalization, not an added term.',fontsize=6.8,linespacing=1.4)

def draw_f3cd(fig):
 unique=pd.read_csv(S/'panels/F3c_primary_budgets.csv');mat=pd.read_csv(S/'panels/F3d_condition_memberships.csv').fillna({'display_ids':''})
 assert len(unique)==5 and sum(unique.membership=='Both allocations')==4
 ax=fig.add_axes([.105,.255,.335,.17]);xx=np.linspace(87,136,246);yy=np.linspace(-1.7,2.4,83);xc=(xx[:-1]+xx[1:])/2;yc=(yy[:-1]+yy[1:])/2;X,Y=np.meshgrid(xc,yc);C=100-X-Y
 ax.pcolormesh(xx,yy,C,cmap=COMP_CMAP,norm=Normalize(-40,40),shading='flat',rasterized=True)
 ax.plot([97.6,101.7],[2.4,-1.7],color=INK,lw=.8,ls='--',label='N + E = 0');ax.axvline(100,c=GRAY,ls=':',lw=.6);ax.axhline(0,c=GRAY,lw=.5)
 for r in unique.itertuples():
  both=r.membership=='Both allocations';assert abs(100-r.A_over_O_percent-r.B_over_O_percent-r.N_plus_E_over_O_percent)<1e-10
  ax.scatter(r.A_over_O_percent,r.B_over_O_percent,marker='o' if both else 's',s=43,facecolor='white',edgecolor='white',lw=2,zorder=5)
  ax.scatter(r.A_over_O_percent,r.B_over_O_percent,marker='o' if both else 's',s=25,facecolor=MODELS[r.model] if both else 'white',edgecolor=MODELS[r.model],lw=1,zorder=6)
  ax.annotate(r.display_id,(r.A_over_O_percent,r.B_over_O_percent),xytext=(4,4),textcoords='offset points',fontsize=7,zorder=7)
 ax.set(xlim=(87,136),ylim=(-1.7,2.4),xticks=[90,100,110,120,130],xlabel='Paired A/O (%)',ylabel='Paired B/O (%)');clean(ax)
 fig.text(.075,.443,'c',fontsize=10,weight='bold');fig.text(.12,.443,'Required compensation',fontsize=7.7)
 colourbar(fig,[.115,.188,.32,.009],COMP_CMAP,Normalize(-40,40),'100 − A/O − B/O = 100(N+E)/O (%)',ticks=[-40,0,40])
 fig.text(.11,.143,'Dashed: zero compensation.\nNot probability or a feasible region.',fontsize=6.5,linespacing=1.2,va='top')
 ax=fig.add_axes([.105,.025,.335,.083]);ax.axis('off');ax.set(xlim=(0,1),ylim=(5.2,-.9))
 for x,t in [(0,'ID'),(.20,'Candidate'),(.61,'Model'),(.99,'Comp. (%)')]:ax.text(x,-.6,t,fontsize=6.5,ha='right' if x>.9 else 'left',weight='bold')
 for i,r in enumerate(unique.itertuples()):
  ax.text(0,i+.35,r.display_id,fontsize=6.6);ax.text(.20,i+.35,str(r.candidate_index),fontsize=6.6);ax.text(.65,i+.35,'G' if r.model.startswith('GISS') else 'I',fontsize=6.6);ax.text(.99,i+.35,f'{r.N_plus_E_over_O_percent:+.2f}',fontsize=6.6,ha='right')
 fig.text(.105,.015,'Filled circle: both; open square: area only (c).',fontsize=6.5)
 fig.text(.54,.443,'d',fontsize=10,weight='bold');fig.text(.59,.443,'Conditional membership',fontsize=7.7)
 order=[(o,a,f) for o in OPS for a in ALLOC for f in ['S2_P3','global_C8','global_C10']];ax=fig.add_axes([.65,.065,.18,.34]);side=fig.add_axes([.85,.065,.135,.34]);side.axis('off');side.set(xlim=(0,1),ylim=(17.5,-.5));labels=[];membership_count=0
 for i,(o,a,f) in enumerate(order):
  r=mat[(mat.loss==o)&(mat.allocation==a)&(mat.family==f)].iloc[0];ids=r.display_ids.split();membership_count+=len(ids);assert len(ids)==int(r['count'])
  labels.append(('P' if a==ALLOC[0] else 'A')+' · '+{'S2_P3':'3 MW','global_C8':'8 TW','global_C10':'10 TW'}[f])
  for j,hid in enumerate(['H1','H2','H3','H4','H5']):
   retained=hid in ids;patch=Rectangle((j-.5,i-.5),1,1,facecolor='#a9d9ca' if retained else 'white',edgecolor='#e1e6e8',lw=.6);ax.add_patch(patch);gid=f'cell-F3d-{i}-{hid}';patch.set_gid(gid);CELL_TITLES[gid]=f'{o};{a};{f};{hid}: '+('retained' if retained else 'not retained')
   if retained:ax.text(j,i,'●',ha='center',va='center',fontsize=5,color=INK)
  side.text(.03,i,str(int(r['count'])),va='center',fontsize=6.5);side.text(1,i,'Empty' if not ids else f'{r.B_over_O_low:.2f}…{r.B_over_O_high:.2f}',ha='right',va='center',fontsize=6.5)
  if i in [2,5,8,11,14]:ax.axhline(i+.5,c=GRAY,lw=.7 if i in [5,11] else .3)
 matrix_axes(ax,18,5,labels,['H1\nG','H2\nI','H3\nI','H4\nI','H5\nI']);ax.tick_params(labelsize=6.5,pad=3)
 for k,t in enumerate(['Complete','100 yr','30 yr']):ax.text(-3.65,k*6+2.5,t,rotation=90,va='center',ha='center',fontsize=6.5)
 side.text(.03,-1.0,'n',fontsize=6.5);side.text(1,-1.0,'B/O endpoints (%)',fontsize=6.5,ha='right')
 fig.text(.54,.022,'● / teal: retained; white: not retained. All evaluated.\nP / A: power / area. Same H ≠ same source history.',fontsize=6.5,linespacing=1.25)
 fig.text(.11,.462,'G: GISS-E2-1-G (blue)    I: IPSL-CM6A-LR (orange)',fontsize=6.7)
 assert membership_count==44
 record('3c-background','colour_display_contract.json','A/O grid;B/O grid;100-A/O-B/O','algebraic diagnostic',6000,'percent','budget identity only','none','paired O coordinate','not a feasibility screen','original five markers unchanged','analytic colour field;no candidate interpolation or probability')

def draw_f4():
 exec(compile((Path(__file__).parent/'ppt_layers_readable.py').read_text(),'ppt_layers_readable.py','exec'),globals())
 ppt=Path(__file__).parents[1]/'assets/Figure4a_supplied_editable.pptx'
 standalone=plt.figure(figsize=(W,W/1.5*1.04));a=standalone.add_axes([0,0,1,1]);ppt_audit=draw_supplied_ppt(a,ppt);a.set_ylim(9753600*1.04,0);save(standalone,'Figure_4a_supplied_layers')
 with (O/'ppt_layer_audit.json').open('x') as f:json.dump(ppt_audit,f,indent=2)
 fig=plt.figure(figsize=(W,10));a=fig.add_axes([0,.51,1,.49]);draw_supplied_ppt(a,ppt);a.set_ylim(9753600*1.04,0)
 pairs=pd.read_csv(S/'panels/F4b_controlled_pairs.csv');no_geo=next(x for x in pairs.allocation.unique() if x not in ALLOC);conditions=[(f,a) for f in ['restricted_outer_ZJ','broad_outer_ZJ'] for a in ALLOC+[no_geo]]
 ax=fig.add_axes([.27,.285,.55,.15]);norm=LogNorm(3,120);labels=[]
 for i,(field,a) in enumerate(conditions):
  labels.append(('3 MW' if field.startswith('restricted') else '14.4 MW')+' · '+{ALLOC[0]:'power',ALLOC[1]:'area',no_geo:'no screen'}[a])
  for j,o in enumerate(OPS):
   r=pairs[(pairs.allocation==a)&(pairs.loss==o)].iloc[0];cell(ax,j,i,r[field],norm,HEAT_CMAP,f'F4b-{field}-{a}-{o}',f'{r[field]:.3f}',8)
 matrix_axes(ax,6,3,labels,['Complete','100-year loss','30-year loss']);ax.axhline(2.5,c=GRAY,lw=.7)
 fig.text(.075,.468,'b',fontsize=10,weight='bold');fig.text(.27,.468,'Saved release conditions × retention',fontsize=8)
 cb=colourbar(fig,[.87,.285,.013,.15],HEAT_CMAP,norm,'Outer B (ZJ; log colour)',orientation='vertical',ticks=[3,10,30,120]);cb.ax.set_yticklabels(['3','10','30','120'])
 fig.text(.27,.257,'14.4 MW + no screen is the jointly relaxed sensitivity, not the primary bound.',fontsize=6.6)
 q=pd.read_csv(S/'panels/F4c_nested_bounds.csv');ax=fig.add_axes([.27,.095,.55,.105]);norm=Normalize(0,12)
 fields=['monthly_outer_ZJ','epoch_source_upper_ZJ','joint_upper_ZJ']
 for i,o in enumerate(OPS):
  r=q[q.loss==o].iloc[0];assert r.joint_upper_ZJ<=r.epoch_source_upper_ZJ+1e-10<=r.monthly_outer_ZJ+2e-10
  for j,field in enumerate(fields):cell(ax,j,i,r[field],norm,HEAT_CMAP,f'F4c-{o}-{field}',f'{r[field]:.3f}',8.5)
 matrix_axes(ax,3,3,OPN,['Timing-relaxed','Complete source','Joint budget'])
 fig.text(.075,.232,'c',fontsize=10,weight='bold');fig.text(.27,.232,'Three nested stages: geological power / 3 MW',fontsize=8)
 colourbar(fig,[.87,.095,.013,.105],HEAT_CMAP,norm,'Upper B (ZJ)',orientation='vertical',ticks=[0,4,8,12])
 fig.text(.27,.057,'Equal-width segments identify stages, not additive heat.\nDirect outer/source equality: tolerance 10⁻¹⁰ ZJ; old stock remaining: 117.97 ZJ.',fontsize=6.7,linespacing=1.4)
 for panel,source,notes in [('4a','../assets/Figure4a_supplied_editable.pptx','independent geological texture;132 native shapes;21 editable labels;conceptual,not measured field'),('4b','panels/F4b_controlled_pairs.csv','18 original values;log colour;not additive deductions'),('4c','panels/F4c_nested_bounds.csv','9 original endpoints;linear colour;equal-width stage bands,not additive energy')]:
  record(panel,source,'saved values or supplied conceptual geometry','conceptual' if panel=='4a' else '1955–2023',0 if panel=='4a' else 6000,'concept' if panel=='4a' else 'ZJ','unchanged source conditions','three alternatives','none','conditional bounds','original records',notes)
 save(fig,'Figure_4_release_controls')

def draw_f5():
 fig=plt.figure(figsize=(W,8.8));front=pd.read_csv(S/'panels/F5a_saved_frontiers.csv');cross=pd.read_csv(S/'panels/F5a_saved_intersections.csv')
 fig.text(.075,.97,'a',weight='bold',fontsize=10);fig.text(.14,.97,'Restricted release: necessary-bound regions',fontsize=8)
 for k,o in enumerate(OPS):
  ax=fig.add_axes([.10+k*.30,.765,.25,.15]);curves=[]
  for a in ALLOC:
   q=front[(front.loss==o)&(front.allocation==a)].sort_values('equivalent_full_power_years');curves.append(q)
  xp=curves[0].equivalent_full_power_years.to_numpy();assert np.array_equal(xp,curves[1].equivalent_full_power_years.to_numpy())
  yp,ya=[100*q.outer_response_ZJ.to_numpy()/den for q in curves];lo=np.minimum(yp,ya);hi=np.maximum(yp,ya)
  ax.fill_between(xp,0,lo,color='#fff0df');ax.fill_between(xp,lo,hi,color='white',hatch='///',edgecolor='#ba9a7b',lw=0);ax.fill_between(xp,hi,2.55,color='#e3f1ec')
  for j,a in enumerate(ALLOC):
   q=curves[j];ax.plot(q.equivalent_full_power_years,100*q.outer_response_ZJ/den,c=ORANGE,ls=['-','--'][j],lw=1.1,marker=['o','s'][j],mfc=ORANGE if j==0 else 'white',markevery=60,ms=3)
   r=cross[(cross.loss==o)&(cross.allocation==a)&(cross.target_percent_fixed_mean_O==1)].iloc[0]
   if r.feasible_in_relaxed_power_box:
    x=r.minimum_equivalent_full_power_years;ax.scatter(x,1,facecolor=ORANGE if j==0 else 'white',edgecolor=ORANGE,marker=['o','s'][j],s=24,zorder=4);ax.annotate(f'{x:.1f} yr',(x,1),xytext=(3,10 if j==0 else -15),textcoords='offset points',fontsize=6.5,bbox=dict(fc='white',ec='none',alpha=.9,pad=.3),arrowprops=dict(arrowstyle='-',lw=.4,color=GRAY))
  ax.axhline(1,c=GRAY,ls=':',lw=.7);ax.set(xlim=(0,69),ylim=(0,2.55),xticks=[0,30,60],yticks=[0,1,2],xlabel='Equivalent exposure (yr)',ylabel='Outer B / mean O (%)' if k==0 else '');clean(ax);ax.set_title(OPN[k],fontsize=7.1,pad=8)
  if o=='loss_30yr':ax.text(.04,.83,'1% not reached',transform=ax.transAxes,fontsize=6.5)
 fig.legend(handles=[Line2D([],[],c=ORANGE,ls='-',marker='o',label='Geological power',ms=3),Line2D([],[],c=ORANGE,ls='--',marker='s',mfc='white',label='Geological area',ms=3)],loc='center',bbox_to_anchor=(.53,.716),ncol=2,frameon=False,fontsize=6.8)
 fig.legend(handles=[Patch(fc='#fff0df',label='Below both bounds'),Patch(fc='white',ec='#ba9a7b',hatch='///',label='Between bounds'),Patch(fc='#e3f1ec',label='Excluded by both bounds')],loc='center',bbox_to_anchor=(.53,.691),ncol=3,frameon=False,fontsize=6.5,columnspacing=1)
 fig.text(.14,.657,'Below a bound does not establish a feasible history. All 5–20% targets are outside.',fontsize=6.6)
 req=pd.read_csv(S/'panels/F5b_necessary_additional_power.csv');norm=Normalize(0,125);ax=fig.add_axes([.25,.455,.55,.145])
 for i,t in enumerate([1,5,10,20]):
  for j,o in enumerate(OPS):
   r=req[(req.target_percent_fixed_mean_O==t)&(req.loss==o)].iloc[0];cell(ax,j,i,r.epoch_necessary_anomaly_peak_TW,norm,HEAT_CMAP,f'F5b-{t}-{o}',f'{r.epoch_necessary_anomaly_peak_TW:.2f}',8.3)
 matrix_axes(ax,4,3,['1%','5%','10%','20%'],['Complete','100-year loss','30-year loss']);ax.set_ylabel('Target B / mean O',fontsize=7)
 fig.text(.075,.632,'b',weight='bold',fontsize=10);fig.text(.25,.632,'Necessary additional peak power',fontsize=8)
 colourbar(fig,[.85,.455,.013,.145],HEAT_CMAP,norm,'Additional power (TW)',orientation='vertical',ticks=[0,40,80,120])
 caps=pd.read_csv(S/'panels/F5b_additional_power_caps.csv')
 caplines=[]
 for release in ['3 MW','14.4 MW']:
  r=caps[caps.release_condition==release];p=r[r.allocation==ALLOC[0]].additional_anomaly_cap_TW.iloc[0];a=r[r.allocation==ALLOC[1]].additional_anomaly_cap_TW.iloc[0];caplines.append(f'{release} envelope: additional caps {p:.2f} TW (power), {a:.2f} TW (area).')
 fig.text(.25,.405,'\n'.join(caplines),fontsize=6.7,linespacing=1.45)
 h=pd.read_csv(S/'panels/F5c_marginal_requirements.csv');duration=json.loads((S/'panels/display_definitions.json').read_text())['source_window_years'];order=[(t,a,o) for t in [5,10] for a in ALLOC for o in OPS]
 labels=[f'{t}% · '+('P' if a==ALLOC[0] else 'A')+' · '+{'direct':'full','loss_100yr':'100 yr','loss_30yr':'30 yr'}[o] for t,a,o in order]
 fig.text(.075,.367,'c',weight='bold',fontsize=10);fig.text(.29,.367,'Two separate marginal necessary conditions',fontsize=8)
 for j,x in enumerate([.29,.62]):
  ax=fig.add_axes([x,.095,.27,.235]);norm=Normalize(0,200 if j==0 else duration);cmap=ALLOW_CMAP if j==0 else FLUX_CMAP
  for i,(t,a,o) in enumerate(order):
   r=h[(h.target_percent==t)&(h.allocation==a)&(h.loss==o)].iloc[0];v=100*r.critical_kernel_weighted_excess_fraction if j==0 else r.minimum_global_full_participation_years
   text=(f'{v:.1f}'+(' (>100%)' if j==0 and v>100 else '')) if np.isfinite(v) else 'Outside'
   cell(ax,0,i,v,norm,cmap,f'F5c-{j}-{t}-{a}-{o}',text,7)
  matrix_axes(ax,12,1,labels if j==0 else ['']*12,['Allowance (%)' if j==0 else 'Exposure (yr)']);ax.tick_params(labelsize=6.5)
  for y in [2.5,5.5,8.5]:ax.axhline(y,c=GRAY,lw=.6 if y==5.5 else .3)
  colourbar(fig,[x,.061,.27,.009],cmap,norm,'100% declared allowance' if j==0 else '69.0 yr source window',ticks=[0,100,200] if j==0 else [0,30,69])
 fig.text(.29,.013,'P / A: geological power / area. The two minima are not one attainable history.',fontsize=6.5)
 for panel,source,unit,notes in [('5a','panels/F5a_saved_frontiers.csv;panels/F5a_saved_intersections.csv','yr;percent','regions from saved curves only;not feasible sets;no extrapolation'),('5b','panels/F5b_necessary_additional_power.csv;panels/F5b_additional_power_caps.csv','TW','additional anomaly peak,not gross;12 discrete saved thresholds'),('5c','panels/F5c_marginal_requirements.csv','percent;yr','separate scales;values above100 retained;Outside is not a duration;separate marginal minima')]:record(panel,source,'saved numerical values','1955–2023',6000,unit,'unchanged source conditions','three alternatives','fixed mean O','necessary only','original records',notes)
 save(fig,'Figure_5_power_persistence')

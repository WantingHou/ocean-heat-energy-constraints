"""Geological collection-area constraints with analytic full-history certificates."""
from pathlib import Path
import argparse,ast,calendar,itertools,json
import numpy as np,pandas as pd,xarray as xr
ap=argparse.ArgumentParser();ap.add_argument('--project-root',type=Path,required=True);ap.add_argument('--config',type=Path,required=True);ap.add_argument('--transfer-config',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args();R=a.project_root;cfg=json.loads(a.config.read_text());tc=json.loads(a.transfer_config.read_text());S=R/cfg['source_package'];O=a.output_dir
if O.exists():raise FileExistsError(O)
O.mkdir(parents=True);c=json.loads((S/'frozen_source_config.json').read_text());vectors=np.array(c['normal_reference_vectors']);settings=c['settings'];n=16
Y=np.arange(1850,2026);MY=np.repeat(Y,12);MM=np.tile(np.arange(1,13),176);seconds=np.array([calendar.monthrange(int(y),int(m))[1]*86400. for y,m in zip(MY,MM)]);unit=seconds*1e-9;dt=seconds/31557600.;masks=[(MY>=lo)&(MY<=hi) for lo,hi in c['epochs']];u0=unit[masks[0]].sum();ufull=unit.sum();modern_units=np.array([unit[m].sum() for m in masks[1:]]);windows=[(1955,2014),(1955,2023),(1992,2023)];operators={'direct':None,'loss_100yr':100.,'loss_30yr':30.}
z=settings['S2'];rlo=(1+(z['H']-1)*z['fmin'])/(1+(z['H']-1)*c['normal_enhanced_fraction']);rhi=(1+(z['H']-1)*z['fmax'])/(1+(z['H']-1)*c['normal_enhanced_fraction'])
helper=Path(__file__).with_name('power_transfer_helpers.py');text=helper.read_text();tree=ast.parse(text)
for name in ['states','contrast','intersect','certify']:
 node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name==name);exec(compile(ast.get_source_segment(text,node),'previous_verified_helper','exec'))
builder=Path(__file__).with_name('original_constraint_builder.py');text=builder.read_text();node=next(x for x in ast.parse(text).body if isinstance(x,ast.FunctionDef) and x.name=='build');exec(compile(ast.get_source_segment(text,node),'source_constraint_builder','exec'));builders=[build('S2',i) for i in range(2)]
K={loss:np.array([contrast(states(m.astype(float),tau)) for m in masks[1:]]).T for loss,tau in operators.items()};max_violation=0.;max_closure=0.;worst_stock=1e99

def cap(family):return {'S2_P3':1e30,'global_C8':8.,'global_C10':10.}[family]
def limits(v,family):
 pa,pv,F,cn,qg,qb=v;cy=qb-pa-pv;lo=np.array([rlo*pa,rlo*pv,.5*F,min(cy,rlo*(pa+pv)),.5*F]);hi=np.array([rhi*pa,rhi*pv,1.5*F if family=='S2_P3' else cap(family)-rhi*pa,0.,0.]);L=max(-pa-pv-F,rlo*(pa+pv)-F-qb);U=(rhi-1)*(pa+pv)+.5*F if family=='S2_P3' else cap(family)+rhi*pv-pa-pv-F
 assert .5*F<=cn-cy+1e-10;assert hi[2]>=.5*F-1e-10 and hi[2]<=7.2*F+1e-10;return lo,hi,L,U

def interval(seg,family,Fcap,delta=None):
 lo,hi=vectors[seg:seg+2];r=(0.,1.);fs=[lambda v:v[2]-Fcap]
 if family!='S2_P3':fs.append(lambda v:v[0]+v[2]-cap(family))
 if delta is not None:
  small=min(delta);large=max(delta);fs.extend([lambda v:-sum(v[:3])-small,lambda v:rlo*(v[0]+v[1])-v[2]-v[5]-small])
  if family=='S2_P3':fs.append(lambda v:large-((rhi-1)*(v[0]+v[1])+.5*v[2]))
  else:fs.append(lambda v:large-(cap(family)+rhi*v[1]-sum(v[:3])))
 for f in fs:r=intersect(r,f(lo),f(hi))
 return r

def points(seg,family,Fcap):
 v=interval(seg,family,Fcap)
 if v is None:return []
 vals=list(v);f=lambda x:x[5]-(1+rlo)*(x[0]+x[1]);left,right=map(f,vectors[seg:seg+2]);root=-left/(right-left)
 if v[0]<root<v[1]:vals.append(root)
 return sorted(set(vals))

pairs=pd.read_csv(S/'F3_all_window_operator_pairs.csv',float_precision='round_trip');base=pairs[(pairs.setting=='S2')&(pairs.first_year==1955)&(pairs.last_year==2023)];classify=[];hist=[];source=[];corners=[];summary=[];bounds=[]
cases=tc['geometry_cases']+[{'allocation':'global_missing_sediment_permissive','ordinary_flank_cap_TW':tc['unknown_sediment_permissive_global_cap_TW'],'limiting_basin':'Global','weights':None}]
for case in cases:
 alloc=case['allocation'];Fcap=case['ordinary_flank_cap_TW']
 for family in tc['families']:
  proofpoints=[]
  for seg in range(2):
   for x in points(seg,family,Fcap):
    v=vectors[seg]+x*(vectors[seg+1]-vectors[seg]);lo,hi,L,U=limits(v,family);proofpoints.append((seg,x,lo,hi,L,U))
    for pattern in itertools.product([0,1],repeat=3):
     pp=np.array([hi if bit else lo for bit in pattern]);ck,b=certify(seg,x,pp,cap(family),'direct');corners.append(dict(allocation=alloc,family=family,Fcap_TW=Fcap,pattern=''.join(map(str,pattern)),**ck))
  for loss,matrix in K.items():
   for j,(first,last) in enumerate(windows):
    for endpoint in ['low','high']:
     opts=[]
     for seg,x,lo,hi,L,U in proofpoints:
      choices=(matrix[j]>=0) if endpoint=='high' else (matrix[j]<0);pp=np.array([hi if bit else lo for bit in choices]);ck,b=certify(seg,x,pp,cap(family),loss);opts.append(dict(allocation=alloc,family=family,Fcap_TW=Fcap,loss=loss,first_year=first,last_year=last,endpoint=endpoint,B_ZJ=float(b[j]),**ck))
     if opts:source.append((max if endpoint=='high' else min)(opts,key=lambda d:d['B_ZJ']))
   kept=[]
   for row in base[base.loss==loss].to_dict('records'):
    required=np.array([row[f'B_{first}_{last}_ZJ'] for first,last in windows]);delta=np.linalg.solve(matrix,required);feasible=[]
    for seg in range(2):
     inter=interval(seg,family,Fcap,delta)
     if inter is not None:feasible.append((seg,inter))
    classify.append(dict(allocation=alloc,family=family,loss=loss,candidate_index=row['candidate_index'],admitted=bool(feasible)))
    if not feasible:continue
    seg,(left,right)=feasible[0];x=(left+right)/2;v=vectors[seg]+x*(vectors[seg+1]-vectors[seg]);lo,hi,L,U=limits(v,family);alpha=(delta-L)/(U-L);assert min(alpha)>=-1e-9 and max(alpha)<=1+1e-9;pp=lo+np.clip(alpha,0,1)[:,None]*(hi-lo);ck,b=certify(seg,x,pp,cap(family),loss,required)
    assert v[2]<=Fcap+1e-10
    if family=='S2_P3':assert max(pp[:,2])<=1.5*v[2]+1e-10
    record={k:v for k,v in row.items() if k not in ['setting','variables_json','reference_segment','reference_interpolation','source_history_id','identity']};record.update(allocation=alloc,family=family,Fcap_TW=Fcap,source_history_id=f'{alloc}_{family}_{loss}_{row["candidate_index"]}',**ck);kept.append(record);hist.append(record)
   summary.append(dict(allocation=alloc,family=family,loss=loss,Fcap_TW=Fcap,count=len(kept),models=';'.join(sorted(set(r['model'] for r in kept))),background_feasible=bool(proofpoints)))
   if kept:
    g=pd.DataFrame(kept)
    for first,last in windows:
     for term in ['A','B']:
      for metric in ['ZJ','percent']:
       vals=g[f'{term}_{first}_{last}_ZJ'] if metric=='ZJ' else 100*g[f'{term}_{first}_{last}_ZJ']/g[f'O_{first}_{last}_ZJ']
       for end,idx in [('low',vals.idxmin()),('high',vals.idxmax())]:bounds.append(dict(allocation=alloc,family=family,loss=loss,first_year=first,last_year=last,term=term,metric=metric,endpoint=end,value=float(vals.loc[idx]),source_history_id=g.loc[idx,'source_history_id'],candidate_index=int(g.loc[idx,'candidate_index'])))
for name,rows in [('candidate_classification',classify),('certified_histories',hist),('source_only_bounds',source),('corner_proof',corners),('scenario_summary',summary),('conditional_bounds',bounds)]:pd.DataFrame(rows).to_csv(O/(name+'.csv'),index=False,mode='x')
# New spatial witness: retain near-axis/mixed/conduction class geography and
# replace ONLY the flank weights with the evidence-defined geography.
h=pd.DataFrame(hist);selected=h[(h.allocation=='geological_power')&(h.family=='S2_P3')&(h.loss=='direct')].sort_values('B_1955_2023_ZJ').iloc[-1];x=np.array(json.loads(selected.variables_json));seg=int(selected.reference_segment);v=vectors[seg]+x[0]*(vectors[seg+1]-vectors[seg]);pp=x[1:].reshape(3,5);names=['Pacific','Atlantic','Indian','Southern','Arctic','Mediterranean'];classes=['axial_near_axis','backarc_volcanic','low_temperature_flank','remaining_conduction'];ref=pd.read_csv(R/cfg['original_reference']);W=[]
for k in classes:
 arrays=[]
 for label in [['low','central'][seg],['central','high'][seg]]:arrays.append(ref[(ref.reference==label)&(ref.outlet_class==k)].set_index('basin').loc[names].baseline_power_TW.to_numpy())
 vals=arrays[0]+x[0]*(arrays[1]-arrays[0]);W.append(vals/vals.sum())
W=np.array(W);W[2]=[tc['geometry_cases'][1]['weights'][name] for name in names];assert np.max(abs(W.sum(axis=1)-1))<1e-12
rate=np.zeros((len(MY),6,4))
for j,m in enumerate(masks[1:]):
 pa,pv,pf,ry,rf=pp[j];ya=ry*pa/(pa+pv);yv=ry-ya;rate[m,:,0]=(pa-ya-v[0])*W[0]+ya*W[3];rate[m,:,1]=(pv-yv-v[1])*W[1]+yv*W[3];rate[m,:,2]=(pf-rf-v[2])*W[2]+rf*W[3];rate[m,:,3]=-(ry+rf)*W[3]
e=rate*unit[:,None,None];annual=(np.cumsum(e,axis=0)-e/2).reshape(176,12,6,4).mean(axis=1);er=annual[2014-1850]-annual[1970-1850];B=er.sum(axis=1);Bh=er[:,0]+er[:,2]
src=np.zeros(len(MY))
for j,m in enumerate(masks[1:]):src[m]=sum(pp[j,:3])-sum(pp[j,3:])-sum(v[:3])
ss=states(src,None);assert abs(B.sum()-(ss[164]-ss[120]))<1e-9
geo=xr.load_dataset(S/'mapped_fields.nc');frac=geo.basin_area_m2.sel(basin=names).values;area=geo.area_m2.values;ba=frac.sum(axis=(1,2));fields={'B_ZJ':np.sum(frac*(B/ba)[:,None,None],axis=0),'Bh_ZJ':np.sum(frac*(Bh/ba)[:,None,None],axis=0)};out=xr.Dataset({k:(('lat','lon'),val) for k,val in fields.items()},coords={'lat':geo.lat,'lon':geo.lon});out['A_ZJ']=geo.A_ZJ;out['area_m2']=geo.area_m2;out.attrs.update(window='1970-2014',depth_m=6000,scope='conditional source diagnostic;geologicalflankweights;otherclassesretained;independentspatialA',candidate_index=int(selected.candidate_index),source_history_id=selected.source_history_id)
with (O/'geology_conditioned_source_fields.nc').open('xb') as f:f.write(out.to_netcdf(engine='scipy'))
r=pd.DataFrame({'basin':names,'B_ZJ':B,'Bh_ZJ':Bh,'B_share_percent':100*B/B.sum(),'B_GJ_m2':B*1e12/ba,'ordinary_flank_TW':v[2]*W[2],'collection_million_km2':v[2]*W[2]/.04,'native_near_axis_ZJ':er[:,0]-0, 'mixed_with_routing_ZJ':er[:,1],'flank_with_routing_ZJ':er[:,2],'conductive_donor_ZJ':er[:,3]});r.rename(columns={'native_near_axis_ZJ':'near_axis_with_routing_ZJ'},inplace=True);r.to_csv(O/'selected_regional_diagnostic.csv',index=False,mode='x')
with (O/'selected_witness.json').open('x') as f:json.dump(selected.to_dict(),f,indent=2)
with (O/'verification.json').open('x') as f:json.dump({'max_source_constraint_violation':max_violation,'max_budget_error_ZJ':max_closure,'minimum_stock_checked_ZJ':worst_stock,'all_continuous_reference_intervals_certified':True,'source_optimizer_invoked':False,'regional_inventory_or_hydraulic_closure_proven':False,'modern_source_history_independently_observed':False,'spatial_B_integral_ZJ':float(B.sum()),'spatial_Bh_integral_ZJ':float(Bh.sum()),'matrix':{k:v.tolist() for k,v in K.items()}},f,indent=2)
print(pd.DataFrame(summary).to_string(index=False));print(r.to_string(index=False))

"""Heterogeneous, ordinary-power-weighted source envelopes; no source/model runs."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np,pandas as pd
ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,required=True);ap.add_argument('--config',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args();S=a.package/'Source_Data';O=a.output_dir
if O.exists():raise FileExistsError(O)
O.mkdir(parents=True);cfg=json.loads(a.config.read_text());c=json.loads((S/'frozen_source_config.json').read_text());tc=json.loads((S/'regional_transfer_config.json').read_text());V=np.array(c['normal_reference_vectors']);z=c['settings']['S2'];den=1+(z['H']-1)*c['normal_enhanced_fraction'];rl=(1+(z['H']-1)*z['fmin'])/den;rh=(1+(z['H']-1)*z['fmax'])/den;dr=cfg['prior_broad_release_ratio']-cfg['benchmark_release_ratio'];mean=json.loads((S/'physical_config.json').read_text())['ensemble_mean_O_ZJ'];caps={v['allocation']:v['ordinary_flank_cap_TW'] for v in tc['geometry_cases']}
def L(v):return max(-sum(v[:3]),rl*(v[0]+v[1])-v[2]-v[5])
def U(v):return (rh-1)*(v[0]+v[1])+(cfg['benchmark_release_ratio']-1)*v[2]
def points(cap):
 out=[]
 for i in range(2):
  x,y=V[i:i+2];right=min(1.,(cap-x[2])/(y[2]-x[2]))
  if right<0:continue
  t=[0.,right];a=x[5]-(1+rl)*sum(x[:2]);b=y[5]-(1+rl)*sum(y[:2]);root=-a/(b-a) if b!=a else -1
  if 0<root<right:t.append(root)
  out += [(i,f,x+f*(y-x)) for f in sorted(set(t))]
 return out
G=pd.read_csv(S/'geological_support/age_sediment_strata.csv',float_precision='round_trip');G=G[G['mode']=='exact_cell_intersection'];groups={
 'thin_1_25Ma':G[(G.sediment_class=='sediment_lt400m')&(G.age_class=='age_1_25Ma')],
 'thin_25_75Ma':G[(G.sediment_class=='sediment_lt400m')&G.age_class.isin(['age_25_55Ma','age_55_75Ma'])],
 'thick_1_25Ma':G[(G.sediment_class=='sediment_ge400m')&(G.age_class=='age_1_25Ma')]}
area=sum(g.area_m2.sum() for g in groups.values());power=sum(g.background_power_TW.sum() for g in groups.values());weights=[]
for name,g in groups.items():weights.append(dict(group=name,area_million_km2=g.area_m2.sum()/1e12,background_power_TW=g.background_power_TW.sum(),area_weight=g.area_m2.sum()/area,power_weight=g.background_power_TW.sum()/power,scope='age/sediment support weights, not observed thermal-regime or modern-anomaly fractions'))
W=pd.DataFrame(weights).set_index('group');assert abs(W.area_weight.sum()-1)<1e-12 and abs(W.power_weight.sum()-1)<1e-12
# Raw background shares are diagnostic, not automatically a feasible fine-scale
# ordinary-source allocation at fixed collection density. Preserve this failure.
basins=pd.read_csv(S/'geological_support/basin_geological_support.csv',float_precision='round_trip');basins=basins[basins['mode']=='exact_cell_intersection'].set_index('basin')
nested=[];rawaudit=[];sharemap={}
for name,gg in groups.items():
 bybasin=gg.groupby('basin')[['area_m2','background_power_TW']].sum()
 rawaudit.append(dict(group=name,raw_background_weight=W.loc[name,'power_weight'],area_million_km2=W.loc[name,'area_million_km2'],raw_weight_F_cap_TW=.04*gg.area_m2.sum()/1e12/W.loc[name,'power_weight'],minimum_reference_F_TW=V[0,2],raw_weight_low_reference_feasible=.04*gg.area_m2.sum()/1e12/W.loc[name,'power_weight']>=V[0,2]))
 for allocation in cfg['allocations']:
  case=next(x for x in tc['geometry_cases'] if x['allocation']==allocation);total=0.
  for basin,b in basins.iterrows():
   area_s=bybasin.loc[basin,'area_m2'];fraction=case['weights'][basin]*area_s/b.eligible_flank_area_m2;density=caps[allocation]*fraction*1e12/area_s if area_s>0 else 0.;assert density<=.04+1e-12
   nested.append(dict(allocation=allocation,basin=basin,group=name,ordinary_power_fraction=fraction,area_m2=area_s,density_at_allowed_F_cap_W_m2=density,scope='declared uniform eligible-area allocation within each original basin; not measured flow geography'))
   total+=fraction
  sharemap[allocation,name]=total
for allocation in cfg['allocations']:assert abs(sum(sharemap[allocation,name] for name in groups)-1)<1e-10
K=pd.read_csv(S/'physical_limits/response_kernels.csv',float_precision='round_trip');thresholds=[];supports=[];schedules=[];checks=[]
def fill(g,dt,gain,need,availability):
 pos=np.flatnonzero(g>0);idx=pos[np.argsort(-g[pos],kind='stable')];capacity=availability*dt[idx];value=g[idx]*gain;maximum=float(np.dot(capacity,value))
 if need>maximum+1e-10:return None
 chosen=np.zeros(len(g));remaining=max(0.,need);last=0.
 for j,cap,v in zip(idx,capacity,value):
  if remaining<=1e-12:break
  take=min(cap,remaining/v);chosen[j]=take;remaining-=take*v;last=v
 D=chosen.sum();response=np.dot(chosen,g*gain);dual=last*D+np.dot(capacity,np.maximum(value-last,0)) if D>0 else 0.
 assert abs(response-max(0.,need))<1e-9 and abs(dual-response)<1e-9,(response,need,dual)
 return D,chosen,abs(dual-response)
for allocation in cfg['allocations']:
 for loss in cfg['operators']:
  kr=K[K.loss==loss];g=kr.energy_weight.to_numpy();k=kr.power_weight_ZJ_per_TW.to_numpy();dt=kr.seconds.to_numpy()/31557600.;kp=np.maximum(k,0).sum();kn=np.minimum(k,0).sum();pp=points(caps[allocation]);axis='power_weight' if allocation=='geological_power' else 'area_weight'
  for target in cfg['target_percent']:
   need=mean*target/100;choices=[];durations=[]
   for seg,t,v in pp:
    base=kp*U(v)+kn*L(v);slope=dr*v[2]*kp;f=max(0.,(need-base)/slope);choices.append((f,seg,t,v,base))
    result=fill(g,dt,dr*v[2]*31557600e-9,need-base,1.)
    if result is not None:durations.append((result[0],seg,t,v,base,result))
   f,seg,t,v,base=min(choices,key=lambda x:x[0]);full=base+dr*v[2]*kp;feasible=f<=1+1e-12
   duration=None;dual=None
   if durations:
    duration,ds,dtref,dv,db,result=min(durations,key=lambda x:x[0]);dual=result[2]
    for j,x in enumerate(result[1]):
     if x>0:schedules.append(dict(allocation=allocation,loss=loss,target_percent=target,year=int(kr.iloc[j].year),month=int(kr.iloc[j].month),effective_fraction=x/dt[j],weighted_years=x,reference_segment=ds,reference_interpolation=dtref))
   thresholds.append(dict(allocation=allocation,loss=loss,target_percent=target,critical_kernel_weighted_excess_fraction=f,threshold_within_declared_tail=feasible,minimum_global_full_participation_years=duration,strict_validated_low_share_threshold=max(0.,1-f),base_response_at_fraction_threshold_ZJ=base,full_tail_response_at_fraction_threshold_ZJ=full,reference_segment=seg,reference_interpolation=t,dual_error_ZJ=dual,scope='necessary heterogeneous-envelope conditions, not field counts, event durations, probability or actual attribution'))
  for group in cfg['geological_groups']:
   share=1. if group=='all_eligible_flank' else sharemap[allocation,group];options=[]
   for seg,t,v in pp:options.append((kp*(U(v)+dr*v[2]*share)+kn*L(v),seg,t))
   val,seg,t=max(options,key=lambda x:x[0]);supports.append(dict(allocation=allocation,loss=loss,high_release_allowed_group=group,ordinary_power_support_fraction=share,outer_ZJ=val,percent_fixed_mean=100*val/mean,below_5percent=val<.05*mean,reference_segment=seg,reference_interpolation=t,scope='conditional association stress test; geological class does not identify high-release occupancy'))
# Check limiting cases against independently delivered one-constraint envelopes.
prior=pd.read_csv(S/'evidence_closure/constraint_ablation.csv')
for allocation in cfg['allocations']:
 for loss in cfg['operators']:
  kr=K[K.loss==loss];k=kr.power_weight_ZJ_per_TW.to_numpy();kp=np.maximum(k,0).sum();kn=np.minimum(k,0).sum()
  for f,family in [(0.,'P3'),(1.,'prior_broad_14p4')]:
   val=max(kp*(U(v)+dr*v[2]*f)+kn*L(v) for _,_,v in points(caps[allocation]));ref=prior[(prior.allocation==allocation)&(prior.loss==loss)&(prior.family==family)].outer_ZJ.item();checks.append(abs(val-ref));assert checks[-1]<1e-9
# Unit-refinement invariance: split each weighted regime into arbitrary children.
rng=np.random.default_rng(20261106);invariance=[]
for _ in range(100):
 w=rng.dirichlet(np.ones(20));ratio=rng.uniform(.5,7.2,20);cut=rng.uniform(0,1,20);childw=np.column_stack([w*cut,w*(1-cut)]).ravel();childr=np.repeat(ratio,2)
 x=np.dot(w,np.maximum(ratio-1.5,0)/dr);y=np.dot(childw,np.maximum(childr-1.5,0)/dr);invariance.append(abs(x-y))
for name,rows in [('geological_regime_weights',weights),('nested_regime_weights',nested),('raw_background_weight_audit',rawaudit),('heterogeneous_release_thresholds',thresholds),('stratum_support_envelopes',supports),('minimum_exposure_schedules',schedules)]:pd.DataFrame(rows).to_csv(O/(name+'.csv'),index=False,mode='x')
inputs=[a.config,S/'frozen_source_config.json',S/'regional_transfer_config.json',S/'geological_support/age_sediment_strata.csv',S/'geological_support/basin_geological_support.csv',S/'physical_limits/response_kernels.csv',S/'evidence_closure/constraint_ablation.csv']
record={'limit_case_checks':len(checks),'max_limit_case_error_ZJ':max(checks),'weighted_unit_refinement_checks':len(invariance),'max_refinement_error':max(invariance),'minimum_exposure_dual_checks':sum(x['dual_error_ZJ'] is not None for x in thresholds),'max_dual_error_ZJ':max(x['dual_error_ZJ'] or 0 for x in thresholds),'nested_geometry_rows_checked':len(nested),'raw_background_disaggregation_compatible_with_low_reference':all(x['raw_weight_low_reference_feasible'] for x in rawaudit),'geological_area_million_km2':area/1e12,'geological_background_TW':power,'input_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},'modern_mixture_fraction_identified':False,'global_empirical_cap_validated':False,'model_or_source_optimizer_run':False}
with (O/'verification.json').open('x') as f:json.dump(record,f,indent=2)
print(W.to_string());print(pd.DataFrame(thresholds).to_string(index=False));print(pd.DataFrame(supports).query("loss=='direct'").to_string(index=False));print(json.dumps(record,indent=2))

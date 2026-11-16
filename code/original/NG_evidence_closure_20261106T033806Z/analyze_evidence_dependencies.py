"""Create-only analytical dependency audit; no climate/source optimizer execution."""
from pathlib import Path
import argparse,json,hashlib,calendar
import numpy as np
import pandas as pd
ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,required=True);ap.add_argument('--config',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args()
P=a.package;S=P/'Source_Data';O=a.output_dir
if O.exists():raise FileExistsError(O)
O.mkdir(parents=True)
cfg=json.loads(a.config.read_text());c=json.loads((S/'frozen_source_config.json').read_text());tc=json.loads((S/'regional_transfer_config.json').read_text());pc=json.loads((S/'physical_config.json').read_text())
V=np.array(c['normal_reference_vectors']);z=c['settings']['S2'];den=1+(z['H']-1)*c['normal_enhanced_fraction'];lo=(1+(z['H']-1)*z['fmin'])/den;hi=(1+(z['H']-1)*z['fmax'])/den
caps={x['allocation']:x['ordinary_flank_cap_TW'] for x in tc['geometry_cases']};caps['no_geometry_screen']=float('inf')
K=pd.read_csv(S/'physical_limits/response_kernels.csv',float_precision='round_trip');mean=pc['ensemble_mean_O_ZJ']
def low(v):return max(-sum(v[:3]),lo*(v[0]+v[1])-v[2]-v[5])
def offset(v):return v[0]+v[2]-(hi-1)*v[1]
def points(cap,C=None):
 out=[]
 for seg in range(2):
  v,w=V[seg:seg+2];left,right=0.,1.
  constraints=[(v[2]-cap,w[2]-cap)] if np.isfinite(cap) else []
  if C is not None:constraints.append((v[0]+v[2]-C,w[0]+w[2]-C))
  for f0,f1 in constraints:
   slope=f1-f0
   if abs(slope)<1e-14:
    if f0>0:left,right=1.,0.
   elif slope>0:right=min(right,-f0/slope)
   else:left=max(left,-f0/slope)
  if left>right:continue
  ts=[left,right];f0=v[5]-(1+lo)*(v[0]+v[1]);f1=w[5]-(1+lo)*(w[0]+w[1])
  if abs(f1-f0)>1e-14:
   root=-f0/(f1-f0)
   if left<root<right:ts.append(root)
  out.extend((seg,t,v+t*(w-v)) for t in sorted(set(ts)))
 return out
def upper(v,P=None,C=None):return C-offset(v) if C is not None else (hi-1)*(v[0]+v[1])+(P/2-1)*v[2]
rows=[];thresholds=[];forward=[]
def forward_contrast(q,label):
 tau=cfg['operators'][label];H=0.;annual={}
 for year in range(1955,2024):
  mm=[]
  for m in range(1,13):
   j=(year-1955)*12+m-1;dt=calendar.monthrange(year,m)[1]*86400/cfg['year_seconds'];energy=q[j]*dt*cfg['year_seconds']*1e-9
   if tau is None:average=H+energy/2;H+=energy
   else:
    x=dt/tau;aa=np.exp(-x);bb=-np.expm1(-x)/x;dd=(1-bb)/x;average=bb*H+dd*energy;H=aa*H+bb*energy
   mm.append(average)
  annual[year]=np.mean(mm)
 return annual[2023]-annual[1955]
for alloc in cfg['allocations']:
 for label in cfg['operators']:
  kk=K[K.loss==label].power_weight_ZJ_per_TW.to_numpy();kp=np.maximum(kk,0).sum();kn=np.minimum(kk,0).sum()
  for family in cfg['families']:
   C={'C8':8.,'C10':10.}.get(family);Pmax={'P3':3.,'prior_broad_14p4':14.4}.get(family);pp=points(caps[alloc],C)
   vals=[(kp*upper(v,Pmax,C)+kn*low(v),seg,t,v) for seg,t,v in pp];B,seg,t,v=max(vals,key=lambda x:x[0]);u=upper(v,Pmax,C);l=low(v)
   q=np.where(kk>=0,u,l);actual=forward_contrast(q,label);assert abs(actual-B)<1e-9
   rows.append(dict(allocation=alloc,family=family,loss=label,outer_ZJ=B,percent_fixed_mean=100*B/mean,below_5percent=B<.05*mean,reference_segment=seg,reference_interpolation=t,U_TW=u,L_TW=l,forward_error_ZJ=abs(actual-B),scope='power-only envelope; stock and chronology relaxed; not attained history or actual attribution'))
  for pct in cfg['targets_percent']:
   target=pct/100*mean
   # On each reference/lower-envelope branch this is linear-fractional with F>0.
   choices=[(2*(1+(target-kn*low(v)-kp*(hi-1)*(v[0]+v[1]))/(kp*v[2])),seg,t,v) for seg,t,v in points(caps[alloc])]
   Pcrit,seg,t,v=min(choices,key=lambda x:x[0]);achieved=kp*upper(v,Pcrit)+kn*low(v);assert abs(achieved-target)<1e-10
   # Independently solve the monotone envelope equation by scalar bisection.
   left,right=0.,1000.
   for _ in range(80):
    mid=(left+right)/2;B=max(kp*upper(v,mid)+kn*low(v) for _,_,v in points(caps[alloc]))
    if B<target:left=mid
    else:right=mid
   assert abs((left+right)/2-Pcrit)<1e-10
   thresholds.append(dict(allocation=alloc,loss=label,target_percent=pct,target_ZJ=target,necessary_equivalent_flank_peak_MW=Pcrit,multiplier_relative_ordinary=Pcrit/2,reference_segment=seg,reference_interpolation=t,bisection_error_MW=abs((left+right)/2-Pcrit),scope='necessary only; fixed 2 MW normalization and nonflank envelope; no 8/10 TW cap imposed; stocks relaxed'))
df=pd.DataFrame(rows);baseline=pd.read_csv(S/'physical_limits/monthly_outer_bounds.csv');errs=[]
for _,old in baseline.iterrows():
 family={'S2_P3':'P3','global_C8':'C8','global_C10':'C10'}[old.family];row=df[(df.allocation==old.allocation)&(df.loss==old.loss)&(df.family==family)].iloc[0];errs.append(abs(row.outer_ZJ-old.monthly_timing_outer_ZJ))
assert max(errs)<1e-10
# The direct source ceiling contains neither a climate-member path nor a fitted O coefficient.
# A missing coupled response need not be geological energy; quantify a scale, not a cause.
channels=[]
for label in cfg['operators']:
 subset=df[(df.loss==label)&(df.allocation!='no_geometry_screen')&(df.family!='prior_broad_14p4')];bound=subset.outer_ZJ.max()
 for pct in cfg['targets_percent']:
  gap=max(0,pct/100*mean-bound)
  equivalent_TW=gap/(68*cfg['year_seconds']*1e-9)
  channels.append(dict(loss=label,target_percent=pct,represented_source_outer_ZJ=bound,additional_response_gap_ZJ=gap,constant_net_exchange_equivalent_TW=equivalent_TW,earth_area_equivalent_W_m2=equivalent_TW*1e12/cfg['global_area_m2'],scope='68-year contrast-equivalent scale only; not an inferred surface flux, forcing or loss-operator input; coupled response unmeasured'))
density=[]
for alloc in ['retained','geological_power','geological_area']:
 for ref,v in zip(['low','central','high'],V):density.append(dict(allocation=alloc,reference=ref,ordinary_flank_TW=v[2],minimum_uniform_collection_density_mW_m2=40*v[2]/caps[alloc],declared_density_mW_m2=40,scope='threshold for fixed weights and mask, not measured collection density; changing density changes unit footprint'))
for name,records in [('constraint_ablation',rows),('critical_flank_power',thresholds),('coupled_response_requirements',channels),('geometry_density_thresholds',density)]:pd.DataFrame(records).to_csv(O/(name+'.csv'),index=False,mode='x')
inputs=[a.config,S/'frozen_source_config.json',S/'regional_transfer_config.json',S/'physical_config.json',S/'physical_limits/response_kernels.csv',S/'physical_limits/monthly_outer_bounds.csv']
verification={'analytic_envelopes':len(rows),'forward_monthly_checks':len(rows),'max_forward_error_ZJ':max(r['forward_error_ZJ'] for r in rows),'critical_thresholds':len(thresholds),'max_threshold_bisection_error_MW':max(r['bisection_error_MW'] for r in thresholds),'prior_envelopes_reproduced':len(errs),'max_prior_envelope_error_ZJ':max(errs),'source_model_or_optimizer_run':False,'new_attainable_histories_certified':False,'empirical_global_upper_bound_established':False,'input_hashes':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs}}
with (O/'verification.json').open('x') as f:json.dump(verification,f,indent=2)
print(df[df.loss=='direct'].to_string(index=False));print(pd.DataFrame(thresholds).to_string(index=False));print(json.dumps(verification,indent=2))

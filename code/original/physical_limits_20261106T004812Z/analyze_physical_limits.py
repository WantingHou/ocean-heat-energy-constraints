"""Analytical physical limits; no source optimizer or climate model execution."""
from pathlib import Path
import argparse, calendar, json, hashlib
import numpy as np
import pandas as pd

ap=argparse.ArgumentParser()
ap.add_argument('--source-package',type=Path,required=True)
ap.add_argument('--config',type=Path,required=True)
ap.add_argument('--output-dir',type=Path,required=True)
a=ap.parse_args(); S=a.source_package; O=a.output_dir
if O.exists(): raise FileExistsError(O)
O.mkdir(parents=True)
cfg=json.loads(a.config.read_text()); c=json.loads((S/'Source_Data/frozen_source_config.json').read_text()); tc=json.loads((S/'regional_transfer_config.json').read_text())
vectors=np.array(c['normal_reference_vectors']); z=c['settings']['S2']
rlo=(1+(z['H']-1)*z['fmin'])/(1+(z['H']-1)*c['normal_enhanced_fraction'])
rhi=(1+(z['H']-1)*z['fmax'])/(1+(z['H']-1)*c['normal_enhanced_fraction'])
years=np.repeat(np.arange(1850,2026),12); months=np.tile(np.arange(1,13),176)
seconds=np.array([calendar.monthrange(int(y),int(m))[1]*86400. for y,m in zip(years,months)])
unit=seconds*1e-9; dt=seconds/31557600.; modern=(years>=1955)&(years<=2023)
epochs=[(years>=lo)&(years<=hi) for lo,hi in c['epochs'][1:]]
windows=[(1955,2014),(1955,2023),(1992,2023)]
h=pd.read_csv(S/'Source_Data/certified_histories.csv',float_precision='round_trip')
sb=pd.read_csv(S/'Source_Data/source_only_bounds.csv',float_precision='round_trip')
pairs=pd.read_csv(S/'Source_Data/F3_all_window_operator_pairs.csv',float_precision='round_trip')
base=pairs[(pairs.setting=='S2')&(pairs.first_year==1955)&(pairs.last_year==2023)]
Omin=float(base.O_1955_2023_ZJ.min()); Omax=float(base.O_1955_2023_ZJ.max())
caps={r['allocation']:r['ordinary_flank_cap_TW'] for r in tc['geometry_cases']}

def kernel(tau, first=1955, last=2023):
    w=np.zeros(len(years)); w[years==first]=-1/12; w[years==last]=1/12
    if tau is None: aa=np.ones(len(dt)); bb=np.ones(len(dt)); dd=np.full(len(dt),.5)
    else:
        x=dt/tau; aa=np.exp(-x); bb=-np.expm1(-x)/x; dd=(1-bb)/x
    g=np.zeros(len(dt)); future=0.
    for j in range(len(dt)-1,-1,-1):
        g[j]=w[j]*dd[j]+future*bb[j]
        future=w[j]*bb[j]+future*aa[j]
    return g, g*unit

def low(v):
    pa,pv,F,cn,qg,qb=v
    return max(-pa-pv-F, rlo*(pa+pv)-F-qb)

def offset(v):
    pa,pv,F,*_=v
    return pa+F-(rhi-1)*pv

def upper(v,family):
    if family=='S2_P3': return (rhi-1)*(v[0]+v[1])+.5*v[2]
    return {'global_C8':8.,'global_C10':10.}[family]-offset(v)

def segment_interval(seg,Fcap,C=None,delta=None,family=None):
    v0,v1=vectors[seg:seg+2]; left,right=0.,1.
    ff=[lambda v:v[2]-Fcap]
    if C is not None:ff.append(lambda v:v[0]+v[2]-C)
    if delta is not None:
        ff += [lambda v:-sum(v[:3])-min(delta),lambda v:rlo*(v[0]+v[1])-v[2]-v[5]-min(delta)]
        if family=='S2_P3':ff.append(lambda v:max(delta)-upper(v,'S2_P3'))
        elif C is not None:ff.append(lambda v:max(delta)+offset(v)-C)
    for f in ff:
        f0,f1=f(v0),f(v1); slope=f1-f0
        if abs(slope)<1e-14:
            if f0>1e-10:return None
        elif slope>0:right=min(right,-f0/slope)
        else:left=max(left,-f0/slope)
    if left>right+1e-10:return None
    return left,right

def points(Fcap,C=None,delta=None,family=None):
    result=[]
    for seg in range(2):
        inter=segment_interval(seg,Fcap,C,delta,family)
        if inter is None:continue
        vals=list(inter);v0,v1=vectors[seg:seg+2]
        f0=v0[5]-(1+rlo)*(v0[0]+v0[1]); f1=v1[5]-(1+rlo)*(v1[0]+v1[1])
        if abs(f1-f0)>1e-14:
            root=-f0/(f1-f0)
            if inter[0]<root<inter[1]:vals.append(root)
        for t in sorted(set(vals)):result.append((seg,t,v0+t*(v1-v0)))
    return result

kr=[]; response=[]; kernels={}; matrices={}
for label,tau in cfg['operators'].items():
    g,k=kernel(tau); kernels[label]=(g,k)
    K=np.array([[kernel(tau,*w)[1][m].sum() for m in epochs] for w in windows]); matrices[label]=K
    for j in np.flatnonzero(modern):kr.append(dict(loss=label,year=int(years[j]),month=int(months[j]),seconds=seconds[j],energy_weight=g[j],power_weight_ZJ_per_TW=k[j]))
    response.append(dict(loss=label,epoch_constant_gain_ZJ_per_TW=float(k[modern].sum()),monthly_positive_gain_ZJ_per_TW=float(np.maximum(k[modern],0).sum()),monthly_negative_gain_ZJ_per_TW=float(np.minimum(k[modern],0).sum()),maximum_positive_energy_weight=float(g[modern].max()),minimum_energy_weight=float(g[modern].min()),negative_months=int((k[modern]<-1e-14).sum()),K_condition_number=float(np.linalg.cond(K))))
res=pd.DataFrame(response).set_index('loss')
outer=[]; requirements=[]; references=[]
for alloc in cfg['allocations']:
    Fcap=caps[alloc]
    for seg,t,v in points(Fcap):references.append(dict(allocation=alloc,segment=seg,interpolation=t,Pa0_TW=v[0],Pv0_TW=v[1],F_TW=v[2],L_TW=low(v),S2_P3_U_TW=upper(v,'S2_P3'),gross_power_offset_TW=offset(v)))
    for label in cfg['operators']:
        rr=res.loc[label]; kp=rr.monthly_positive_gain_ZJ_per_TW; kn=rr.monthly_negative_gain_ZJ_per_TW
        for family in cfg['families']:
            C=None if family=='S2_P3' else {'global_C8':8.,'global_C10':10.}[family]
            pp=points(Fcap,C)
            choices=[(kp*upper(v,family)+kn*low(v),seg,t,v) for seg,t,v in pp]
            val,seg,t,v=max(choices,key=lambda x:x[0]); old=sb[(sb.allocation==alloc)&(sb.family==family)&(sb.loss==label)&(sb.first_year==1955)&(sb.last_year==2023)&(sb.endpoint=='high')].iloc[0]
            outer.append(dict(allocation=alloc,family=family,loss=label,monthly_timing_outer_ZJ=val,epoch_source_upper_ZJ=old.B_ZJ,timing_relaxation_ZJ=val-old.B_ZJ,reference_segment=seg,reference_interpolation=t,anomaly_upper_TW=upper(v,family),anomaly_lower_TW=low(v),percent_fixed_mean_O=100*val/cfg['ensemble_mean_O_ZJ'],percent_min_retained_O=100*val/Omin,scope='power-only outer bound; stock relaxed; no pre1955 anomaly; not actual attribution'))
        # Relax the gross cap to an arbitrary peak C; retain background/mixed power
        # and geology. Dropping the ordinary C check can only enlarge the envelope.
        pp=points(Fcap); intercept,seg,t,v=max([(-kp*offset(v)+kn*low(v),s,x,v) for s,x,v in pp],key=lambda x:x[0])
        minoff=min(offset(v) for s,x,v in pp); Lmin=min(low(v) for s,x,v in pp)
        recovery=kn*Lmin
        for frac in cfg['target_percent']:
            target=cfg['ensemble_mean_O_ZJ']*frac/100
            E=max(0.,target-recovery)/rr.maximum_positive_energy_weight
            area=E*1e21/(3330*1171*60*300)/1e6
            requirements.append(dict(allocation=alloc,loss=label,target_percent_fixed_mean_O=frac,target_ZJ=target,epoch_necessary_anomaly_peak_TW=target/rr.epoch_constant_gain_ZJ_per_TW,epoch_necessary_axis_flank_peak_TW=target/rr.epoch_constant_gain_ZJ_per_TW+minoff,monthly_necessary_axis_flank_peak_TW=(target-intercept)/kp,nonnegative_anomaly_peak_lower_TW=target/kp,maximum_baseline_recovery_allowance_ZJ=recovery,necessary_positive_extra_input_ZJ=E,storage_only_equivalent_area_km2=area,storage_only_fraction_eligible_percent=100*area/151825323.990614,scope='necessary only; fixed-mean target; storage-only equivalent excludes ordinary supply and replenishment'))

support=[]; deletion=[]; slack=[]; critical=[]
for alloc in cfg['allocations']:
    for family in cfg['families']:
        for label in cfg['operators']:
            g=h[(h.allocation==alloc)&(h.family==family)&(h.loss==label)].copy()
            outer_row=next(r for r in outer if r['allocation']==alloc and r['family']==family and r['loss']==label)
            if g.empty:
                support.append(dict(allocation=alloc,family=family,loss=label,count=0,epoch_source_upper_ZJ=outer_row['epoch_source_upper_ZJ'],monthly_outer_ZJ=outer_row['monthly_timing_outer_ZJ']))
                continue
            g['ratio']=100*g.B_1955_2023_ZJ/g.O_1955_2023_ZJ
            top=g.loc[g.ratio.idxmax()]
            support.append(dict(allocation=alloc,family=family,loss=label,count=len(g),model_count=g.model.nunique(),upper_candidate=int(top.candidate_index),upper_model=top.model,joint_upper_percent=top.ratio,joint_upper_ZJ=top.B_1955_2023_ZJ,paired_O_ZJ=top.O_1955_2023_ZJ,epoch_source_upper_ZJ=outer_row['epoch_source_upper_ZJ'],monthly_outer_ZJ=outer_row['monthly_timing_outer_ZJ'],source_headroom_ZJ=outer_row['epoch_source_upper_ZJ']-top.B_1955_2023_ZJ))
            for field in ['candidate_index','model','O_member']:
                for value in g[field].unique():
                    remaining=g[g[field]!=value]
                    deletion.append(dict(allocation=alloc,family=family,loss=label,deleted_field=field,deleted_value=value,count_remaining=len(remaining),upper_percent=None if remaining.empty else remaining.ratio.max(),lower_percent=None if remaining.empty else remaining.ratio.min(),scope='finite support dependence; not confidence interval'))
            for kind,row in [('joint_ratio_upper',top),('source_energy_upper',sb[(sb.allocation==alloc)&(sb.family==family)&(sb.loss==label)&(sb.first_year==1955)&(sb.last_year==2023)&(sb.endpoint=='high')].iloc[0])]:
                x=np.array(json.loads(row.variables_json));seg=int(row.reference_segment);v=vectors[seg]+x[0]*(vectors[seg+1]-vectors[seg]);pp=x[1:].reshape(3,5);q=pp[:,:3].sum(axis=1)-pp[:,3:].sum(axis=1)-sum(v[:3]);delta=matrices[label]@q
                slack.append(dict(allocation=alloc,family=family,loss=label,kind=kind,F_TW=v[2],geometry_slack_TW=caps[alloc]-v[2],anomaly_peak_TW=q.max(),scalar_upper_slack_TW=upper(v,family)-q.max(),source_B_ZJ=delta[1],minimum_young_stock_ZJ=row.minimum_young_stock_ZJ,minimum_old_stock_ZJ=row.minimum_old_stock_ZJ,minimum_ordinary_old_stock_ZJ=row.minimum_ordinary_old_stock_ZJ,required_connected_thickness_m=row.required_connected_thickness_m,thickness_headroom_m=300-row.required_connected_thickness_m,gross_axis_flank_peak_TW=row.maximum_gross_axis_flank_TW))
# Determine the reference-support thresholds of the direct ratio-maximizing tuple.
top=h[(h.allocation=='geological_power')&(h.family=='S2_P3')&(h.loss=='direct')].sort_values('B_over_O_percent').iloc[-1]
b=np.array([top[f'B_{f}_{l}_ZJ'] for f,l in windows])
for label,K in matrices.items():
    d=np.linalg.solve(K,b)
    for alloc in cfg['allocations']:
        pp=points(caps[alloc],delta=d)
        vals=[]
        for seg in range(2):
            inter=segment_interval(seg,caps[alloc],delta=d)
            if inter is None:continue
            v0,v1=vectors[seg:seg+2]; vals_t=list(inter)
            f=lambda v:v[0]+v[2]-max(d)-offset(v)
            f0,f1=f(v0),f(v1)
            if abs(f1-f0)>1e-14:
                root=-f0/(f1-f0)
                if inter[0]<root<inter[1]:vals_t.append(root)
            for t in vals_t:
                v=v0+t*(v1-v0);vals.append((max(v[0]+v[2],max(d)+offset(v)),seg,t,v))
        if vals:
            C,seg,t,v=min(vals,key=lambda z:z[0])
            critical.append(dict(allocation=alloc,loss=label,candidate_index=int(top.candidate_index),delta_1955_1985_TW=d[0],delta_1986_2010_TW=d[1],delta_2011_2025_TW=d[2],necessary_peak_axis_flank_cap_TW=C,reference_segment=seg,reference_interpolation=t,scope='exact scalar cap threshold before inventory and broad flank-ceiling checks; necessary only'))

for name,rows in [('response_kernels',kr),('response_summaries',response),('reference_breakpoints',references),('monthly_outer_bounds',outer),('target_requirements',requirements),('support_attribution',support),('deletion_diagnostics',deletion),('witness_slacks',slack),('selected_tuple_thresholds',critical)]:
    pd.DataFrame(rows).to_csv(O/(name+'.csv'),index=False,mode='x')
inputs=[a.config,S/'Source_Data/frozen_source_config.json',S/'regional_transfer_config.json',S/'Source_Data/certified_histories.csv',S/'Source_Data/source_only_bounds.csv',S/'Source_Data/F3_all_window_operator_pairs.csv']
with (O/'analysis_record.json').open('x') as f:json.dump(dict(config=cfg,rlo=rlo,rhi=rhi,retained_pool_O_min_ZJ=Omin,retained_pool_O_max_ZJ=Omax,retained_pool_observation_members=int(base.O_member.nunique()),retained_pool_unique_candidates=int(base.candidate_index.nunique()),input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},matrix={k:v.tolist() for k,v in matrices.items()},optimizer_invoked=False,model_invoked=False),f,indent=2)
print(pd.DataFrame(response).to_string(index=False))
print(pd.DataFrame(outer).to_string(index=False))
print(pd.DataFrame(requirements).query("allocation=='geological_power'").to_string(index=False))
print(pd.DataFrame(critical).to_string(index=False))

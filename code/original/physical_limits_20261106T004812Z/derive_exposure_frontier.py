"""Exact rearrangement exposure frontiers for the retained response operators."""
from pathlib import Path
import argparse,json,calendar
import pandas as pd,numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--config',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args()
if a.output_dir.exists():raise FileExistsError(a.output_dir)
a.output_dir.mkdir(parents=True);cfg=json.loads(a.config.read_text());r=a.results
kernels=pd.read_csv(r/'response_kernels.csv');refs=pd.read_csv(r/'reference_breakpoints.csv');outer=pd.read_csv(r/'monthly_outer_bounds.csv');rows=[];inversions=[];gaps=[];verify=[]
for loss,tau in cfg['operators'].items():
    k=kernels[kernels.loss==loss];g=k.energy_weight.to_numpy();u=k.seconds.to_numpy()*1e-9;positive=g>0
    idx=np.flatnonzero(positive)[np.argsort(-g[positive],kind='stable')];gg=g[idx];uu=u[idx];cum_u=np.r_[0.,np.cumsum(uu)];cum_w=np.r_[0.,np.cumsum(uu*gg)];negative_weight=float(np.minimum(g*u,0).sum())
    def weighted_duration(D):return float(np.interp(D*31557600e-9,cum_u,cum_w))
    for alloc in cfg['allocations']:
        vr=refs[refs.allocation==alloc]
        for family in cfg['families']:
            candidates=[]
            for v in vr.itertuples():
                C=None if family=='S2_P3' else {'global_C8':8.,'global_C10':10.}[family]
                if C is not None and v.Pa0_TW+v.F_TW>C+1e-10:continue
                U=v.S2_P3_U_TW if C is None else C-v.gross_power_offset_TW
                candidates.append((U,v.L_TW,v.segment,v.interpolation))
            # All reference intervals have ordinary axis+flank below8TW here;
            # no extra cap-intersection breakpoints are required.
            assert all(v.Pa0_TW+v.F_TW<8 for v in vr.itertuples())
            for D in np.linspace(0,69,277):
                W=weighted_duration(D);B,U,L,seg,t=max([(U*W+L*negative_weight,U,L,s,t) for U,L,s,t in candidates],key=lambda x:x[0])
                rows.append(dict(allocation=alloc,family=family,loss=loss,equivalent_full_power_years=D,outer_response_ZJ=B,reference_segment=seg,reference_interpolation=t,anomaly_cap_TW=U,baseline_recovery_ZJ=L*negative_weight))
            full=max(U*cum_w[-1]+L*negative_weight for U,L,s,t in candidates)
            previous=outer[(outer.allocation==alloc)&(outer.family==family)&(outer.loss==loss)].iloc[0].monthly_timing_outer_ZJ
            verify.append(dict(allocation=alloc,family=family,loss=loss,full_duration_error_ZJ=abs(full-previous)))
            for f in cfg['target_percent']:
                T=f/100*cfg['ensemble_mean_O_ZJ'];options=[]
                for U,L,s,t in candidates:
                    targetW=(T-L*negative_weight)/U
                    if targetW<=cum_w[-1]+1e-12:
                        D=float(np.interp(max(0,targetW),cum_w,cum_u))/(31557600e-9)
                        options.append((D,U,L,s,t))
                if options:D,U,L,s,t=min(options,key=lambda x:x[0])
                inversions.append(dict(allocation=alloc,family=family,loss=loss,target_percent_fixed_mean_O=f,target_ZJ=T,feasible_in_relaxed_power_box=bool(options),minimum_equivalent_full_power_years=D if options else None,scope='not uninterrupted duration; power-only relaxation; not sufficient stock or transport'))
    # Initial anomaly is defined at the start of January1955. Its decay is not
    # represented by the source calculations; quantify the excluded legacy.
    secs=k.seconds.to_numpy();dy=secs/31557600.;ends=np.cumsum(dy);starts=ends-dy
    if tau is None:legacy=0.
    else:
        mean=np.exp(-starts/tau)*(-np.expm1(-dy/tau))/(dy/tau)
        legacy=float(mean[k.year.to_numpy()==2023].mean()-mean[k.year.to_numpy()==1955].mean())
    ceiling=float(outer[outer.loss==loss].monthly_timing_outer_ZJ.max())
    for f in cfg['target_percent']:
        T=f/100*cfg['ensemble_mean_O_ZJ'];gap=max(0,T-ceiling)
        gaps.append(dict(loss=loss,target_percent_fixed_mean_O=f,target_ZJ=T,most_permissive_represented_outer_ZJ=ceiling,additional_retained_response_needed_ZJ=gap,initial1955_legacy_coefficient=legacy,required_negative_initial1955_anomaly_ZJ=gap/(-legacy) if legacy<0 else None,scope='separate missing-source/negative-legacy alternatives; not observed or jointly inferred'))
for name,data in [('exposure_frontier',rows),('exposure_thresholds',inversions),('unrepresented_response_requirements',gaps)]:pd.DataFrame(data).to_csv(a.output_dir/(name+'.csv'),index=False,mode='x')
assert max(x['full_duration_error_ZJ'] for x in verify)<1e-9
with (a.output_dir/'verification.json').open('x') as f:json.dump(dict(full_duration_checks=verify,proof='Exchange argument: higher energy-weight months receive capped positive input before lower ones. Piecewise-linear cumulative gain is inverted without an optimizer.',no_retention_or_temperature_scan=True),f,indent=2)
print(pd.DataFrame(inversions).query('target_percent_fixed_mean_O==1').to_string(index=False));print(pd.DataFrame(gaps).to_string(index=False))

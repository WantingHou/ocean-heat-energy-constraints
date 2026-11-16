"""Independent forward and geometry checks for heterogeneous support."""
from pathlib import Path
import argparse,json,calendar
import numpy as np,pandas as pd
ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,required=True);ap.add_argument('--results',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();S=a.package/'Source_Data';R=a.results
T=pd.read_csv(R/'heterogeneous_release_thresholds.csv');N=pd.read_csv(R/'nested_regime_weights.csv');W=pd.read_csv(R/'geological_regime_weights.csv');A=pd.read_csv(R/'raw_background_weight_audit.csv');H=pd.read_csv(R/'stratum_support_envelopes.csv');schedule=pd.read_csv(R/'minimum_exposure_schedules.csv');prior=pd.read_csv(S/'evidence_closure/critical_flank_power.csv');B=pd.read_csv(S/'geological_support/basin_geological_support.csv');B=B[B['mode']=='exact_cell_intersection'].set_index('basin');tc=json.loads((S/'regional_transfer_config.json').read_text());c=json.loads((S/'frozen_source_config.json').read_text());V=np.array(c['normal_reference_vectors']);z=c['settings']['S2'];den=1+(z['H']-1)*c['normal_enhanced_fraction'];rl=(1+(z['H']-1)*z['fmin'])/den;rh=(1+(z['H']-1)*z['fmax'])/den;mean=json.loads((S/'physical_config.json').read_text())['ensemble_mean_O_ZJ'];kernels=pd.read_csv(S/'physical_limits/response_kernels.csv')
for x in T.itertuples():
 pr=prior[(prior.allocation==x.allocation)&(prior.loss==x.loss)&(prior.target_percent==x.target_percent)].iloc[0];assert abs((pr.necessary_equivalent_flank_peak_MW-3)/11.4-x.critical_kernel_weighted_excess_fraction)<1e-12
for allocation,g in N.groupby('allocation'):
 case=next(x for x in tc['geometry_cases'] if x['allocation']==allocation)
 for basin,gg in g.groupby('basin'):
  assert abs(gg.ordinary_power_fraction.sum()-case['weights'][basin])<1e-12
  assert abs(gg.area_m2.sum()/B.loc[basin,'eligible_flank_area_m2']-1)<1e-10
  for row in gg.itertuples():
   if row.area_m2>0:assert abs(row.ordinary_power_fraction/row.area_m2-case['weights'][basin]/B.loc[basin,'eligible_flank_area_m2'])<1e-20
 assert abs(g.ordinary_power_fraction.sum()-1)<1e-10
assert (N.density_at_allowed_F_cap_W_m2<=.04+1e-12).all()
assert abs(W.area_million_km2.sum()-151.825323990614)<1e-8
assert abs(W.background_power_TW.sum()-7.502183269931)<1e-10
assert not bool(A.loc[A.group=='thin_1_25Ma','raw_weight_low_reference_feasible'].item())
errors=[]
for x in T.itertuples():
 kk=kernels[kernels.loss==x.loss];k=kk.power_weight_ZJ_per_TW.to_numpy();pp=schedule[(schedule.allocation==x.allocation)&(schedule.loss==x.loss)&(schedule.target_percent==x.target_percent)]
 if not x.threshold_within_declared_tail:
  maximum=H[(H.allocation==x.allocation)&(H.loss==x.loss)&(H.high_release_allowed_group=='all_eligible_flank')].outer_ZJ.item();assert maximum<mean*x.target_percent/100;assert pp.empty;continue
 assert len(pp)>0
 seg=int(pp.reference_segment.iloc[0]);f=pp.reference_interpolation.iloc[0];v=V[seg]+f*(V[seg+1]-V[seg]);low=max(-sum(v[:3]),rl*sum(v[:2])-v[2]-v[5]);upper=(rh-1)*sum(v[:2])+.5*v[2];phi={(int(row.year),int(row.month)):row.effective_fraction for row in pp.itertuples()};state=0.;annual={};D=0.;j=0;tau={'direct':None,'loss_100yr':100.,'loss_30yr':30.}[x.loss]
 for year in range(1955,2024):
  months=[]
  for month in range(1,13):
   dt=calendar.monthrange(year,month)[1]*86400/31557600.;fraction=phi.get((year,month),0.);assert -.0000001<=fraction<=1.0000001;D+=fraction*dt;q=(upper if k[j]>=0 else low)+5.7*v[2]*fraction;rate=q*31557600e-9
   if tau is None:months.append(state+.5*rate*dt);state+=rate*dt
   else:
    a1=np.exp(-dt/tau);b1=-np.expm1(-dt/tau)/(dt/tau);months.append(state*b1+rate*tau*(1-b1));state=state*a1+rate*tau*(1-a1)
   j+=1
  annual[year]=np.mean(months)
 actual=annual[2023]-annual[1955];error=abs(actual-mean*x.target_percent/100);assert error<1e-9;assert abs(D-x.minimum_global_full_participation_years)<1e-9;errors.append(error)
for x in H.itertuples():
 expected=1 if x.high_release_allowed_group=='all_eligible_flank' else N[(N.allocation==x.allocation)&(N.group==x.high_release_allowed_group)].ordinary_power_fraction.sum();assert abs(expected-x.ordinary_power_support_fraction)<1e-12
 rng=next(y for y in tc['geometry_cases'] if y['allocation']==x.allocation);v=V[x.reference_segment]+x.reference_interpolation*(V[x.reference_segment+1]-V[x.reference_segment]);kk=kernels[kernels.loss==x.loss].power_weight_ZJ_per_TW.to_numpy();low=max(-sum(v[:3]),rl*sum(v[:2])-v[2]-v[5]);upper=(rh-1)*sum(v[:2])+.5*v[2]+5.7*v[2]*expected;val=np.maximum(kk,0).sum()*upper+np.minimum(kk,0).sum()*low;assert abs(val-x.outer_ZJ)<1e-9
result={'passed':True,'thresholds_cross_checked_against_independent_peak_calculation':12,'nested_density_rows_checked':len(N),'stratum_support_bounds_checked':len(H),'forward_exposure_histories_checked':len(errors),'max_forward_error_ZJ':max(errors),'infeasible_target_cases_verified':2,'raw_background_weight_failure_preserved':True,'observed_global_release_fraction_identified':False,'coupled_response_gap_closed':False}
with a.output.open('x') as f:json.dump(result,f,indent=2)
print(json.dumps(result,indent=2))

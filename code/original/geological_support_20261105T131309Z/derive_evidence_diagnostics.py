"""Derive capacity/density requirements and the explicit retained-evidence time audit."""
from pathlib import Path
import argparse,json
import pandas as pd
ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--source-data',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args();D=a.results;Q=D/'conditional_budget';S=a.source_data;O=a.output_dir
if O.exists():raise FileExistsError(O)
O.mkdir(parents=True);ledger=pd.read_csv(S/'retained_local_evidence.csv')
roles={'Baby_Bare':('local_gross_power','no','snapshot','power_calibration'),'Dorado':('patch_gross_power_and_temperature','temperature_only','2013–2014','retrospective_external_local'),'Grotto':('one_plume_power_series','yes_one_plume','2011–2013','retrospective_external_local'),'Axial_Baker':('relative_plume_enhancement','relative_not_absolute_total','1985–2017','activity_calibration_not_independent'),'Lucky_Strike':('imagery_and_field_power_range','not_quantified_power_difference','1996–2009','contextual'),'ASHES_Jackson':('heat_flux_density','yes_density_not_full_field_power','2019–2021','nested_Axial'),'EPR_PNAS':('temperature_history','temperature_not_power','1991–2025','external_temperature_only'),'Wu_segment':('focused_plus_diffuse_power_estimates','not_comparable_full_power_history','2018/2019/2021','contextual'),'Tagoro':('multiple_field_power_estimators','not_independent_time_endpoints','2014–2023','contextual')}
rows=[]
for key,(quantity,repeat,period,role) in roles.items():
 r=ledger[ledger.id==key].iloc[0];rows.append({'evidence_id':key,'quantity':quantity,'repeat_status':repeat,'retained_period':period,'role':role,'matched_1955_baseline_identified':False,'global_energy_weight_identified':False,'eligible_for_global_1955_2023_extra_power_reconstruction':False,'unknown_global_increment_assigned_zero':False,'reason':r.modern_extra_status,'source':r.source,'scope':'audit of retained evidence, not proof that no other measurements exist'})
pd.DataFrame(rows).to_csv(O/'temporal_evidence_bridge_audit.csv',index=False,mode='x')
b=pd.read_csv(D/'basin_geological_support.csv');b=b[b['mode']=='exact_cell_intersection'].set_index('basin');g=pd.read_csv(D/'reference_geometry_audit.csv');rows=[]
for r in g[g.basin!='Global'].itertuples():
 ar=b.loc[r.basin,'eligible_flank_area_m2'];p=b.loc[r.basin,'power_flank_W'];rows.append({'reference':r.reference,'basin':r.basin,'required_ordinary_density_mW_m2':r.ordinary_flank_TW*1e15/ar,'geological_background_flank_density_mW_m2':p*1000/ar,'fixed_equivalent_density_mW_m2':40,'interpretation':'failed fixed-density geometry can indicate wrong footprint, weight or screening mask; not absence of hydrothermal activity'})
pd.DataFrame(rows).to_csv(O/'collection_density_requirements.csv',index=False,mode='x')
h=pd.read_csv(Q/'source_only_bounds.csv');h=h[(h.allocation=='geological_power')&(h.loss=='direct')&(h.endpoint=='high')&(h.first_year==1955)&(h.last_year==2023)];v=json.loads((S/'frozen_source_config.json').read_text())['normal_reference_vectors'];rows=[]
for r in h.itertuples():
 F=v[r.reference_segment][2]+r.reference_interpolation*(v[r.reference_segment+1][2]-v[r.reference_segment][2]);N=F/2e-6;conn=N*1e6;need=conn*r.required_connected_thickness_m/300;A=b.eligible_flank_area_m2.sum();rows.append({'family':r.family,'source_only_upper_B_ZJ':r.B_ZJ,'ordinary_flank_TW':F,'required_collection_million_km2':F/.04,'assumed_connected_million_km2':conn/1e12,'minimum_connected_million_km2_at300m60K':need/1e12,'minimum_connected_fraction_of_eligible_flank_percent':100*need/A,'required_rock_excess_at300m_K':r.required_connected_thickness_m*60/300,'old_remaining_2023_ZJ':r.old_remaining_2023_ZJ,'scope':'global inventory threshold, not independently measured rock temperature or per-basin hydraulic reachability'})
pd.DataFrame(rows).to_csv(O/'geology_capacity_requirements.csv',index=False,mode='x')

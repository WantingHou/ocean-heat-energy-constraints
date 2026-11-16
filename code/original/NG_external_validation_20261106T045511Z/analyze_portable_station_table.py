"""Create-only, unfitted SO240 regional transfer check; no area-weighted inference."""
from pathlib import Path
import argparse, io, json, hashlib, math
import pandas as pd
ap=argparse.ArgumentParser();inputs=ap.add_mutually_exclusive_group(required=True);inputs.add_argument('--raw',type=Path);inputs.add_argument('--stations',type=Path);ap.add_argument('--support',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
if a.output.exists():raise FileExistsError(a.output)
a.output.mkdir(parents=True)
source=a.raw if a.raw else a.stations
if a.raw:
 raw=a.raw.read_text();d=pd.read_csv(io.StringIO(raw.split('*/',1)[1].strip()),sep='\t')
else:
 d=pd.read_csv(a.stations,float_precision='round_trip')
support=json.loads(a.support.read_text())
d=d.rename(columns={'Station':'station','Event':'event','Latitude':'latitude','Longitude':'longitude','HF [mW/m**2]':'conductive_mW_m2','k [W/m/K]':'conductivity_W_m_K','T grad [mK/m]':'gradient_mK_m','Sensor c [#]':'sensor_count','Date/Time':'date_time','Elevation [m a.s.l.]':'archive_elevation_m'})
d['profile']=d.station.str.extract(r'(HF\d+)');d['primary']=d.profile.isin(support['primary_profiles']);d['assumed_conductivity']=d.conductivity_W_m_K.eq(.83)
assert len(d)==109 and d.station.is_unique and d.primary.sum()==105
assert not d[['conductive_mW_m2','conductivity_W_m_K','gradient_mK_m']].isna().any().any()
p=d.groupby('profile').agg(count=('station','size'),mean_conductive_mW_m2=('conductive_mW_m2','mean'),min_conductive_mW_m2=('conductive_mW_m2','min'),max_conductive_mW_m2=('conductive_mW_m2','max'),assumed_conductivity_count=('assumed_conductivity','sum'),primary=('primary','first')).reset_index()
means={'primary_station':float(d.loc[d.primary,'conductive_mW_m2'].mean()),'primary_profile_equal':float(p.loc[p.primary,'mean_conductive_mW_m2'].mean()),'all_cruise_station':float(d.conductive_mW_m2.mean()),'published_2017':48.1}
rows=[]
for agg in ['primary_station','primary_profile_equal']:
 for age in [18.,20.,22.]:
  qref=506.7/math.sqrt(age);ratio=.66-.0111*age;pred=qref*ratio;obs=qref-means[agg]
  rows.append(dict(aggregation=agg,age_Ma=age,conductive_mean_mW_m2=means[agg],cooling_reference_mW_m2=qref,inferred_deficit_mW_m2=obs,frozen_age_deficit_mW_m2=pred,ordinary_density_mW_m2=40.,deficit_minus_density_mW_m2=obs-40.,deficit_to_density=obs/40.,deficit_minus_age_prediction_mW_m2=obs-pred,predicted_conductive_mW_m2=qref-pred,deficit_fraction=obs/qref))
comparison=pd.DataFrame(rows)
d.to_csv(a.output/'stations.csv',index=False,mode='x');p.to_csv(a.output/'profiles.csv',index=False,mode='x');comparison.to_csv(a.output/'regional_comparison.csv',index=False,mode='x')
summary=dict(input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),input_kind='original_PANGAEA_table' if a.raw else 'preserved_station_CSV',n_all=109,n_primary=105,n_primary_profiles=13,n_all_profiles=14,means=means,primary_assumed_conductivity_count=int(d.loc[d.primary,'assumed_conductivity'].sum()),all_assumed_conductivity_count=int(d.assumed_conductivity.sum()),published_minus_archive_mean=48.1-means['all_cruise_station'],erratum_text_verified=False,statistical_CI_computed=False,area_weights_available=False,refit=False,modern_change_identified=False,conclusion='Fixed-density and frozen age-deficit regional mean predictions understate the inferred deficit under the declared cooling/closure interpretation; no universal-cap validation, no area-integrated measured extraction, no global parameter update.')
with (a.output/'summary.json').open('x') as f:json.dump(summary,f,indent=2)
print(json.dumps(summary,indent=2));print(comparison.to_string(index=False))

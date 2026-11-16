"""Classify existing evidence and sampling support; no parameter fitting."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np,pandas as pd,xarray as xr
ap=argparse.ArgumentParser();ap.add_argument('--project-root',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args();R=a.project_root;O=a.output_dir
if O.exists():raise FileExistsError(O)
O.mkdir(parents=True)
P=R/'docs/manuscript/NG_panel_color_revision_20261108T045452Z/Source_Data'
paths={'stations':P/'SO240_stations_work_areas.csv','catalogue':R/'data/raw/ocean_heat_1955_2025_package_20261029T134654Z/Ocean_Warming_1955_2025/seafloor/interridge_v34_20200325.csv','catalogue_ids':P/'hydrothermal_catalogue_locations.csv','vent_power':R/'data/processed/global_recorded_heat_budget_20261030T051153Z/reused_vent_power_measurements.csv','crosswalk':R/'data/processed/global_recorded_heat_budget_20261030T051153Z/global_source_crosswalk_basins_checked.csv','local_evidence':R/'data/processed/global_seafloor_heat_input_20261103T101537Z/local_power_time_evidence.csv','age':R/'data/raw/geological_support_20261105T131309Z/Seton_age_GTS2012/age.nc','sediment':R/'data/interim/geological_support_20261105T131309Z/GlobSed_extracted/GlobSed-v3.nc','weights':R/'data/processed/global_release_support_20261106T040811Z/nested_mixture_constraints/nested_regime_weights.csv'}
st=pd.read_csv(paths['stations']);cat=pd.read_csv(paths['catalogue']);ids=pd.read_csv(paths['catalogue_ids']);vent=pd.read_csv(paths['vent_power']);cross=pd.read_csv(paths['crosswalk'])
assert cat['Name.ID'].tolist()==ids['Name.ID'].tolist();assert len(st)==109 and len(cat)==721 and len(vent)==120
identity=st.event.astype(str)+'|'+st.station.astype(str)+'|'+st.date_time.astype(str)
assert identity.is_unique and vent.measurement_id.is_unique and ids.catalogue_id.is_unique
basin={'N. Pacific':'Pacific','S. Pacific':'Pacific','N. Atlantic':'Atlantic','S. Atlantic':'Atlantic','Indian':'Indian','Arctic':'Arctic','Mediterranean':'Mediterranean','Southern':'Southern'}
tectonic={'mid-ocean ridge':'mid_ocean_ridge_context','arc volcano':'arc_volcano_context','back-arc spreading center':'back_arc_context','intra-plate volcano':'intraplate_volcano_context','other':'other_context','NotProvided':'unknown'}
c=pd.DataFrame({'record_id':ids.catalogue_id,'source_type':'hydrothermal_occurrence','source':'InterRidge_v3.4','basin':cat.Ocean.map(basin),'latitude':cat.Latitude,'longitude':cat.Longitude,'physical_context':cat['Tectonic.setting'].map(tectonic),'sampling_unit':'catalogued_field_not_station','cruise_id':pd.NA,'station_id':pd.NA,'quality':cat.Activity,'q_mW_m2':np.nan,'power_kW':np.nan})
vbase=cross[cross.source=='VentPower'][['measurement_id','basin']];assert vbase.measurement_id.is_unique
v=vent.merge(vbase,on='measurement_id',how='left',validate='one_to_one');assert v.basin.notna().all()
vclass=np.where(v.site.str.contains('Lau Basin'),'back_arc_context','mid_ocean_ridge_context')
v=pd.DataFrame({'record_id':v.measurement_id,'source_type':'local_vent_power','source':'Germanovich2015_compilation','basin':v.basin,'latitude':v.latitude,'longitude':v.longitude,'physical_context':vclass,'sampling_unit':'measurement_not_independent_vent','cruise_id':v.cruise,'station_id':pd.NA,'site_label':v.site,'quality':v.calibration_status,'q_mW_m2':np.nan,'power_kW':v.reported_heat_power,'primary_direct_device_record':v.primary_direct_device_record})
s=pd.DataFrame({'record_id':identity,'source_type':'conductive_station_heat_flux','source':'SO240_2024_archive','basin':'Pacific','latitude':st.latitude,'longitude':st.longitude,'physical_context':'ridge_flank_survey','sampling_unit':'station_record','cruise_id':'SO240','station_id':st.station,'profile':st.profile,'work_area':st.work_area,'primary':st.primary,'quality':np.where(st.assumed_conductivity,'assumed_k','measured_k'),'q_mW_m2':st.conductive_mW_m2,'power_kW':np.nan})
records=pd.concat([s,c,v],ignore_index=True);assert len(records)==950
# Context lookup only: nearest native grid cell, not measured geology or interpolation.
for name,out in [('age','grid_age_Ma'),('sediment','grid_sediment_m')]:
 with xr.open_dataset(paths[name]) as ds:
  vals=ds.z.sel(lon=xr.DataArray(records.longitude.to_numpy(),dims='record'),lat=xr.DataArray(records.latitude.to_numpy(),dims='record'),method='nearest').values.astype(float)
 records[out]=np.where(np.isfinite(vals)&(vals>=0),vals,np.nan)
records['age_class']=pd.cut(records.grid_age_Ma,[-np.inf,1,25,75,np.inf],right=False,labels=['<1 Ma','1–<25 Ma','25–<75 Ma','>=75 Ma']).astype('object').fillna('unknown')
records['sediment_class']=np.where(records.grid_sediment_m.isna(),'unknown',np.where(records.grid_sediment_m<400,'thin <400 m','thick >=400 m'))
records['geology_role']='nearest-native-grid context; not measured station geology or sampled-area weight'
records.to_csv(O/'record_physical_classification.csv',index=False,mode='x')
cols=['source','source_type','basin','physical_context','age_class','sediment_class'];groups=[]
for key,g in records.groupby(cols,dropna=False,sort=True):
 row=dict(zip(cols,key));row.update(records=len(g),station_count=int(g.station_id.nunique()) if g.source_type.iloc[0]=='conductive_station_heat_flux' else pd.NA,identified_cruises=int(g.cruise_id.nunique()) if g.cruise_id.notna().any() else pd.NA,records_without_cruise=int(g.cruise_id.isna().sum()),independent_statistical_samples='not established',global_weight='not derived from record counts');groups.append(row)
pd.DataFrame(groups).to_csv(O/'basin_physical_sampling_support.csv',index=False,mode='x')
quality=records.groupby(['source','basin','quality'],dropna=False).size().reset_index(name='records');quality.to_csv(O/'quality_support_counts.csv',index=False,mode='x')
source_summary=[]
for source,g in records.groupby('source'):
 for b,bb in g.groupby('basin'):
  source_summary.append({'source':source,'basin':b,'records':len(bb),'identified_cruises':bb.cruise_id.nunique() if bb.cruise_id.notna().any() else pd.NA,'cruise_unknown_records':int(bb.cruise_id.isna().sum()),'station_count':bb.station_id.nunique() if source=='SO240_2024_archive' else pd.NA,'site_labels_not_stations':bb.site_label.nunique() if source=='Germanovich2015_compilation' else pd.NA,'support_note':'Distinct cruise IDs are acquisition support, not independent spatial samples; counts across strata must not be summed as independent cruises.'})
pd.DataFrame(source_summary).to_csv(O/'source_basin_support.csv',index=False,mode='x')
# Preserve existing local evidence, including unknown weights and non-admitted integrals.
local=pd.read_csv(paths['local_evidence']);local['n_independent_cruises']=pd.NA;local['n_standardized_heatflow_stations']=pd.NA;local['sampling_count_status']='not established by this summary ledger; not set to zero';local.to_csv(O/'local_power_evidence_inventory.csv',index=False,mode='x')
with (O/'retained_nested_regime_weights.csv').open('xb') as f:f.write(paths['weights'].read_bytes())
rows=[
('SO240_2024_archive','conductive observations','Pacific',109,1,109,'18–22 Ma field support; native-grid context separately recorded','thin sediment field context','56/105 primary k assumed; uncorrected sedimentation/topography; unresolved erratum','regional external test, not global calibration'),
('InterRidge_v3.4','occurrence catalogue','six adopted basins',721,None,None,'native-grid lookup, missing retained','native-grid lookup, missing retained','confirmed/inferred/inactive; uneven discovery; cruises not resolved','location and discovery support, not heat measurements'),
('Germanovich2015_compilation','local vent power','Pacific',120,int(vent.cruise.nunique()),None,'native-grid context only','native-grid context only','65 within calibration;39 below range;16 legacy; repeated outlet measurements','local source-power evidence, not conductive q or independent field sample'),
('Davies2013_grid','global background product','global',None,None,None,'product-supported background','not station sediment metadata','2-degree gridded product; underlying station/cruise independence not reconciled','background reference, not modern excess heat'),
('Hasterok2013','published global heat-flow deficit diagnostic','global',None,None,None,'age-stratified published relation','thin <400 m / thick >=400 m','underlying station data not harmonized here; overlaps older cooling/heat-flow evidence','frozen physical alternative, not independent global validation'),
('Winslow_Fisher2015','local literature/model transfer','Pacific',None,None,None,'Baby Bare warm flank end-member','field-specific outcrop geometry','used to define condition; cannot independently validate it','3 MW declared equivalent-unit actual-power ceiling, not universal measured cap'),
('equivalent_unit_normalization','accounting convention','not a sampled population',None,None,None,'not inferred from station ages','50 km2 declared collection area','2 MW divided by50 km2 =40 mW m-2','normalization only; not global measured ceiling'),
('IHFC_2024_attempt','unavailable station compilation','global',None,None,None,'unavailable','unavailable','retained xlsx file is zero bytes; prior retrieval failed','not admitted; no station or cruise counts invented'),
('Coggon2025','figure-design reference: crustal carbon','not heat-flow sampling',None,None,None,'geological ages, not modern warming years','talus breccia / upper crust','not a thermal observation dataset','classification and sampling-support presentation only')]
columns=['source','source_type','basin','records','identified_cruises','station_count','age_support','sediment_support','quality_limit','admissible_role']
pd.DataFrame(rows,columns=columns).to_csv(O/'evidence_source_inventory.csv',index=False,mode='x')
summary={'records_classified':len(records),'SO240_records':109,'SO240_primary':105,'SO240_profiles':14,'SO240_primary_profiles':13,'SO240_identified_cruises':1,'SO240_identity_duplicates':int(identity.duplicated().sum()),'InterRidge_records':721,'vent_power_records':120,'vent_identified_cruises':int(vent.cruise.nunique()),'vent_records_unknown_cruise':int(vent.cruise.isna().sum()),'vent_site_labels_not_standardized_stations':int(vent.site.nunique()),'vent_quality':vent.calibration_status.value_counts().to_dict(),'grid_context_missing':{c:int(records[c].isna().sum()) for c in ['grid_age_Ma','grid_sediment_m']},'multi_region_heatflow_validation_ready':False,'reason':'Occurrence, background grids, vent power and published area diagnostics do not constitute a harmonized multi-region station heat-flux sample with independent sampling and uncertainty. Local power comparisons remain heterogeneous evidence, not substituted q.','figure2_selection':'regional SO240 external test retained; global evidence and coverage inventory in SI','station_counts_used_for_global_weights':False,'source_parameters_changed':False,'main_window':'1955–2023','2025':'separate auxiliary comparisons; no splice','input_sha256':{k:{'path':str(p.relative_to(R)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in paths.items()}}
with (O/'evidence_audit.json').open('x') as f:json.dump(summary,f,indent=2)
print(json.dumps({k:v for k,v in summary.items() if k!='input_sha256'},indent=2))

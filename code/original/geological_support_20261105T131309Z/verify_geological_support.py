"""Independent sparse-overlap spatial integration and affine-root budget checks."""
from pathlib import Path
import argparse,calendar,json
import numpy as np,pandas as pd,xarray as xr
from scipy.sparse import csr_matrix
ap=argparse.ArgumentParser();ap.add_argument('--project-root',type=Path,required=True);ap.add_argument('--config',type=Path,required=True);ap.add_argument('--results',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();R=a.project_root;D=a.results;C=json.loads(a.config.read_text())
if a.output.exists():raise FileExistsError(a.output)
ag=xr.load_dataset(R/C['age_grid']);sd=xr.load_dataset(R/C['sediment_grid']);age=ag.z.values.astype(float);sed=sd.z.values;rad=C['earth_radius_m']
def edges(x,ends):return np.concatenate(([ends[0]],(x[1:]+x[:-1])/2,[ends[1]]))
def overlap(x,y,latitude=False):
 rows=[];cols=[];vals=[];i=j=0
 while i<len(x)-1 and j<len(y)-1:
  left=max(x[i],y[j]);right=min(x[i+1],y[j+1])
  if right>left:
   val=np.sin(np.deg2rad(right))-np.sin(np.deg2rad(left)) if latitude else np.deg2rad(right-left);rows.append(i);cols.append(j);vals.append(val)
  if x[i+1]<y[j+1]:i+=1
  elif x[i+1]>y[j+1]:j+=1
  else:i+=1;j+=1
 return csr_matrix((vals,(rows,cols)),shape=(len(x)-1,len(y)-1))
ae=edges(ag.lat.values,[-90,90]);ax=edges(ag.lon.values,[-180,180]);se=edges(sd.lat.values,[-90,90]);sx=edges(sd.lon.values,[-180,180]);lat=overlap(ae,se,True);lon=overlap(ax,sx)
def area(mask):return rad**2*(lon@(lat@mask.astype(float)).T).T
thin_area=area(np.isfinite(sed)&(sed<400));thick_area=area(np.isfinite(sed)&(sed>=400));known=thin_area+thick_area;full=rad**2*np.diff(np.sin(np.deg2rad(ae)))[:,None]*np.deg2rad(np.diff(ax))[None,:];valid=np.isfinite(age);t=np.where(valid,age,1.)
thin_eligible=valid&(t<75);thick_eligible=valid&(t<25);eligible=thin_area*thin_eligible+thick_area*thick_eligible;flank=eligible*(t>=1);unknown=(full-known)*valid
q=np.where(t<=48.1,506.7*t**-.5,53+106*np.exp(-.034607*t));rt=np.where(t<55,.66-.0111*t,np.where(t<75,.186-.0025*t,0.));rk=np.where(t<25,.40-.016*t,0.);power=(thin_area*np.maximum(rt,0)+thick_area*np.maximum(rk,0))*q*.001*valid
numbers={'age_area_m2':float(full[valid].sum()),'known_area_m2':float(known[valid].sum()),'unknown_sediment_area_m2':float(unknown.sum()),'eligible_area_m2':float(eligible.sum()),'eligible_flank_area_m2':float(flank.sum()),'power_all_W':float(power.sum()),'power_flank_W':float((power*(t>=1)).sum())};reported=pd.read_csv(D/'global_geological_support.csv').set_index('mode').loc['exact_cell_intersection'];relative={k:abs(v-reported[k])/max(abs(reported[k]),1) for k,v in numbers.items()};assert max(relative.values())<1e-9,relative
print('Sparse-overlap global verification passed',flush=True)
# Independent response coefficients by backward adjoint.
S=R/C['source_package'];conf=json.loads((S/'frozen_source_config.json').read_text());vv=np.array(conf['normal_reference_vectors']);z=conf['settings']['S2'];rlo=1/(1+(z['H']-1)*conf['normal_enhanced_fraction']);rhi=z['H']*rlo;ys=np.repeat(np.arange(1850,2026),12);ms=np.tile(np.arange(1,13),176);sec=np.array([calendar.monthrange(int(y),int(m))[1]*86400. for y,m in zip(ys,ms)]);u=sec*1e-9;dt=sec/31557600.;masks=[(ys>=lo)&(ys<=hi) for lo,hi in conf['epochs'][1:]];windows=[(1955,2014),(1955,2023),(1992,2023)];K={}
for loss,tau in [('direct',None),('loss_100yr',100),('loss_30yr',30)]:
 if tau is None:aa=np.ones(len(u));bb=aa.copy();dd=np.full(len(u),.5)
 else:
  v=dt/tau;aa=np.exp(-v);bb=-np.expm1(-v)/v;dd=(1-bb)/v
 rows=[]
 for first,last in windows:
  w=((ys==last).astype(float)-(ys==first).astype(float))/12;future=0.;g=np.zeros(len(u))
  for j in range(len(u)-1,-1,-1):g[j]=w[j]*dd[j]+future*bb[j];future=w[j]*bb[j]+future*aa[j]
  rows.append([sum(g[m]*u[m]) for m in masks])
 K[loss]=np.array(rows)
BD=D/'conditional_budget';ss=pd.read_csv(BD/'scenario_summary.csv');Fcap=ss.groupby('allocation').Fcap_TW.first().to_dict();pairs=pd.read_csv(S/'F3_all_window_operator_pairs.csv',float_precision='round_trip');pairs=pairs[(pairs.setting=='S2')&(pairs.first_year==1955)&(pairs.last_year==2023)].set_index(['loss','candidate_index']);d=pd.read_csv(BD/'candidate_classification.csv');count=0
for row in d.itertuples():
 p=pairs.loc[(row.loss,row.candidate_index)];delta=np.linalg.solve(K[row.loss],[p[f'B_{lo}_{hi}_ZJ'] for lo,hi in windows]);found=False
 for seg in range(2):
  values=[]
  for v in vv[seg:seg+2]:
   pa,pv,F,cn,qg,qb=v;g=[F-Fcap[row.allocation],-pa-pv-F-min(delta),rlo*(pa+pv)-F-qb-min(delta)]
   if row.family=='S2_P3':g.append(max(delta)-((rhi-1)*(pa+pv)+.5*F))
   else:
    cap=8 if row.family=='global_C8' else 10;g+=[pa+F-cap,max(delta)-(cap+rhi*pv-pa-pv-F)]
   values.append(g)
  left=np.array(values[0]);slope=np.array(values[1])-left;points=[0.,1.]
  for f,change in zip(left,slope):
   if abs(change)>1e-13:
    root=-f/change
    if 0<=root<=1:points.append(root)
  points.sort();points+=list((np.array(points[:-1])+np.array(points[1:]))/2)
  if any(np.max(left+x*slope)<1e-9 for x in points):found=True
 assert found==bool(row.admitted),(row.allocation,row.family,row.loss,row.candidate_index);count+=1
worst=0.;stockdiff=0.;checks=0
for name in ['certified_histories','source_only_bounds','corner_proof']:
 df=pd.read_csv(BD/(name+'.csv'),float_precision='round_trip')
 for row in df.itertuples():
  xx=np.array(json.loads(row.variables_json));v=vv[row.reference_segment]+xx[0]*(vv[row.reference_segment+1]-vv[row.reference_segment]);pa,pv,F,cn,qg,qb=v;cy=qb-pa-pv;p=xx[1:].reshape(3,5);assert F<=row.Fcap_TW+1e-10
  assert np.min(p)>=-1e-10 and np.all(p[:,0]>=rlo*pa-1e-10) and np.all(p[:,0]<=rhi*pa+1e-10) and np.all(p[:,1]>=rlo*pv-1e-10) and np.all(p[:,1]<=rhi*pv+1e-10)
  assert np.all(p[:,2]>=.5*F-1e-10) and np.all(p[:,2]<=7.2*F+1e-10) and np.all(p[:,3]<=cy+1e-10) and np.all(p[:,3]<=p[:,0]+p[:,1]+1e-10)
  assert np.all(p[:,4]<=cn-cy+1e-10) and np.all(p[:,4]<=3.5*F+1e-10) and np.all(p[:,4]<=p[:,2]+1e-10)
  if row.family=='S2_P3':assert max(p[:,2])<=1.5*F+1e-10
  else:
   cap=8 if row.family=='global_C8' else 10;assert pa+F<=cap+1e-10 and max(p[:,0]+p[:,2])<=cap+1e-10
  y0=qg*50*31557600e-9;old0=F/2e-6*1e6*300*3330*1171*60/1e21-qg*u.sum()-y0;young=y0;old=old0;ylist=[young];olist=[old]
  for k,(lo,hi) in enumerate(conf['epochs']):
   duration=sum(calendar.monthrange(y,m)[1]*86400e-9 for y in range(lo,hi+1) for m in range(1,13));yr=qb if k==0 else p[k-1,0]+p[k-1,1]+cy-p[k-1,3];fr=F if k==0 else p[k-1,2]-p[k-1,4];young+=(qg-yr)*duration;old-=fr*duration;ylist.append(young);olist.append(old)
  assert min(ylist+olist+[old0-F*u.sum()])>=-1e-8;stockdiff=max(stockdiff,abs(min(olist)-row.minimum_old_stock_ZJ))
  delta=p[:,:3].sum(axis=1)-p[:,3:].sum(axis=1)-sum(v[:3])
  if name=='certified_histories':
   required=np.array([getattr(row,f'B_{lo}_{hi}_ZJ') for lo,hi in windows]);err=float(abs(K[row.loss]@delta-required).max());assert err<1e-8;worst=max(worst,err)
  elif name=='source_only_bounds':assert abs(K[row.loss][windows.index((row.first_year,row.last_year))]@delta-row.B_ZJ)<1e-8
  checks+=1
assert stockdiff<1e-8
# Area constraints verified individually for the selected geography, not just min cap.
reg=pd.read_csv(BD/'selected_regional_diagnostic.csv').set_index('basin');area=pd.read_csv(D/'basin_geological_support.csv');area=area[area['mode']=='exact_cell_intersection'].set_index('basin');assert np.all(reg.collection_million_km2<=area.eligible_flank_area_m2/1e12+1e-10)
geo=xr.load_dataset(BD/'geology_conditioned_source_fields.nc');assert abs(float(geo.B_ZJ.sum())-reg.B_ZJ.sum())<1e-10
result={'independent_sparse_overlap_values':numbers,'relative_spatial_differences':relative,'classified_combinations':count,'all_classifications_match':True,'full_source_history_records_checked':checks,'max_independent_budget_error_ZJ':worst,'max_independent_stock_difference_ZJ':stockdiff,'selected_regional_collection_areas_pass':True,'selected_spatial_integral_matches':True,'geological_mask_does_not_prove_hydraulic_connection':True,'modern_source_anomaly_not_identified_from_static_grids':True}
with a.output.open('x') as f:json.dump(result,f,indent=2)
print(json.dumps(result,indent=2))

"""Area-conservative geological support from age/sediment cell intersections."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np,pandas as pd,xarray as xr
from scipy.spatial import cKDTree
ap=argparse.ArgumentParser();ap.add_argument('--project-root',type=Path,required=True);ap.add_argument('--config',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args();R=a.project_root;C=json.loads(a.config.read_text());O=a.output_dir
if O.exists():raise FileExistsError(O)
O.mkdir(parents=True)
age=xr.load_dataset(R/C['age_grid']);sed=xr.load_dataset(R/C['sediment_grid']);native=xr.load_dataset(R/C['native_geometry']);av=age.z.values.astype(float);sv=sed.z.values.astype(float)
assert av.shape==(1801,3601) and sv.shape==(2161,4321);assert np.nanmin(av)>0 and np.nanmin(sv)>=0
ax=age.lon.values;ay=age.lat.values;sx=sed.lon.values;sy=sed.lat.values;rad=C['earth_radius_m'];names=['Pacific','Atlantic','Indian','Southern','Arctic','Mediterranean']
def edges(x,lo,hi):return np.r_[lo,(x[:-1]+x[1:])/2,hi]
aex=edges(ax,-180,180);aey=edges(ay,-90,90);sex=edges(sx,-180,180);sey=edges(sy,-90,90)
assert np.all(np.diff(aex)>0) and np.all(np.diff(sex)>0)
def xyz(lat,lon):
 la=np.deg2rad(lat);ll=np.deg2rad(lon);return np.column_stack([np.cos(la)*np.cos(ll),np.cos(la)*np.sin(ll),np.sin(la)])
wet=(native.volume_0_6000_m3.values>0)&(native.basin_index.values>=0)&(native.basin_index.values<6);tree=cKDTree(xyz(native.latitude.values[wet],native.longitude.values[wet]));labels=native.basin_index.values[wet];basin=np.full(av.shape,-1,dtype=np.int8);distance=np.full(av.shape,np.nan,dtype=np.float32)
for start in range(0,len(ay),90):
 block=av[start:start+90];jj,ii=np.nonzero(np.isfinite(block));dist,idx=tree.query(xyz(ay[start+jj],ax[ii]),workers=2);basin[start+jj,ii]=labels[idx];distance[start+jj,ii]=(2*rad*np.arcsin(np.clip(dist/2,0,1))/1000).astype(np.float32)
print('Native basin assignment complete',flush=True)
# Each nodal value represents its clipped Voronoi control area. Shared dateline
# nodes receive half widths, rather than each receiving an entire cell.
xedge=np.unique(np.round(np.r_[aex,sex,np.arange(-180,181,2)],10));yedge=np.unique(np.round(np.r_[aey,sey,np.arange(-90,91,2)],10));xm=(xedge[:-1]+xedge[1:])/2;ym=(yedge[:-1]+yedge[1:])/2
ix=np.searchsorted(aex,xm,side='right')-1;jx=np.searchsorted(sex,xm,side='right')-1;iy=np.searchsorted(aey,ym,side='right')-1;jy=np.searchsorted(sey,ym,side='right')-1
metrics=['age_area_m2','known_area_m2','unknown_sediment_area_m2','eligible_area_m2','eligible_flank_area_m2','thin_flank_area_m2','thick_flank_area_m2','unknown_possible_flank_area_m2','power_all_W','power_flank_W','negative_fit_clipping_W','basin_distance_over_200km_area_m2']
mapkeys=['known_area_m2','eligible_area_m2','eligible_flank_area_m2','power_all_W','power_flank_W','unknown_sediment_area_m2'];maps={k:np.zeros((90,180)) for k in mapkeys};sums={};stratum=[]

def process(mode,latmid,lonmid,latidx,lonidx,sedlatidx,sedlonidx,lat_edges,lon_edges,make_maps):
 total=np.zeros((6,len(metrics)));hist=np.zeros((6,3,5,2));maxdist=0.;sphere=0.;lonwidth=np.deg2rad(np.diff(lon_edges));lonbin=np.clip(((lonmid+180)//2).astype(int),0,179)
 for start in range(0,len(latmid),48):
  end=min(len(latmid),start+48);al=latidx[start:end];sl=sedlatidx[start:end];t=av[np.ix_(al,lonidx)];th=sv[np.ix_(sl,sedlonidx)];b=basin[np.ix_(al,lonidx)];dist=distance[np.ix_(al,lonidx)];ar=rad**2*np.diff(np.sin(np.deg2rad(lat_edges[start:end+1])))[:,None]*lonwidth[None,:];sphere+=float(ar.sum())
  ageok=np.isfinite(t);known=ageok&np.isfinite(th);unknown=ageok&~np.isfinite(th);thin=known&(th<400);thick=known&(th>=400);eligible=(thin&(t<75))|(thick&(t<25));flank=eligible&(t>=1);missingpossible=unknown&(t>=1)&(t<75)
  tref=np.where(ageok,t,1.);qref=np.where(tref<=48.1,506.7/np.sqrt(tref),53+106*np.exp(-.034607*tref));ratio=np.zeros(t.shape);ratio[thin&(t<55)]=.66-.0111*t[thin&(t<55)];ratio[thin&(t>=55)&(t<75)]=.186-.0025*t[thin&(t>=55)&(t<75)];ratio[thick&(t<25)]=.40-.016*t[thick&(t<25)]
  clipping=np.maximum(-ratio,0)*qref*1e-3*ar;power=np.maximum(ratio,0)*qref*1e-3*ar
  vals=[ar*ageok,ar*known,ar*unknown,ar*eligible,ar*flank,ar*(thin&flank),ar*(thick&flank),ar*missingpossible,power,power*(t>=1),clipping,ar*(ageok&(dist>200))]
  for k,val in enumerate(vals):total[:,k]+=np.bincount(b[ageok],weights=val[ageok],minlength=6)
  maxdist=max(maxdist,float(np.nanmax(dist)) if np.isfinite(dist).any() else 0.)
  sb=np.where(thin,0,np.where(thick,1,2));ab=np.searchsorted([1,25,55,75],t,side='right');code=((b*3+sb)*5+ab).astype(int)
  hist[:,:,:,0]+=np.bincount(code[ageok],weights=ar[ageok],minlength=90).reshape(6,3,5);hist[:,:,:,1]+=np.bincount(code[ageok],weights=power[ageok],minlength=90).reshape(6,3,5)
  if make_maps:
   latbin=np.clip(((latmid[start:end]+90)//2).astype(int),0,89);bins=(latbin[:,None]*180+lonbin[None,:]).ravel()
   for key in mapkeys:maps[key]+=np.bincount(bins,weights=vals[metrics.index(key)].ravel(),minlength=90*180).reshape(90,180)
 assert abs(sphere-4*np.pi*rad**2)/(4*np.pi*rad**2)<1e-12
 records=[]
 for k,name in enumerate(names):records.append(dict(mode=mode,basin=name,**dict(zip(metrics,total[k])),maximum_assignment_distance_km=maxdist))
 for b,name in enumerate(names):
  for si,label in enumerate(['sediment_lt400m','sediment_ge400m','sediment_unknown']):
   for ai,agelabel in enumerate(['age_lt1Ma','age_1_25Ma','age_25_55Ma','age_55_75Ma','age_ge75Ma']):stratum.append(dict(mode=mode,basin=name,sediment_class=label,age_class=agelabel,area_m2=hist[b,si,ai,0],background_power_TW=hist[b,si,ai,1]/1e12))
 return records,sphere

records,sphere=process('exact_cell_intersection',ym,xm,iy,ix,jy,jx,yedge,xedge,True);sums['exact_cell_intersection']=records
# Independent registration diagnostic, not another independent evidence source.
sni=np.clip(np.searchsorted(sex,ax,side='right')-1,0,len(sx)-1);snj=np.clip(np.searchsorted(sey,ay,side='right')-1,0,len(sy)-1)
near,_=process('nearest_sediment_on_age_nodes',ay,ax,np.arange(len(ay)),np.arange(len(ax)),snj,sni,aey,aex,False);sums['nearest_sediment_on_age_nodes']=near
basins=pd.DataFrame(records+near);basins['eligible_flank_million_km2']=basins.eligible_flank_area_m2/1e12;basins['background_flank_TW']=basins.power_flank_W/1e12;basins['background_all_TW']=basins.power_all_W/1e12;basins.to_csv(O/'basin_geological_support.csv',index=False,mode='x');pd.DataFrame(stratum).to_csv(O/'age_sediment_strata.csv',index=False,mode='x')
globalrows=[]
for mode in sums:
 g=basins[basins['mode']==mode];tot=g[metrics].sum();globalrows.append(dict(mode=mode,**tot.to_dict(),eligible_million_km2=tot.eligible_area_m2/1e12,eligible_flank_million_km2=tot.eligible_flank_area_m2/1e12,background_all_TW=tot.power_all_W/1e12,background_flank_TW=tot.power_flank_W/1e12,global_reference_flank_cap_TW=.04*tot.eligible_flank_area_m2/1e12,missing_sediment_permissive_flank_cap_TW=.04*(tot.eligible_flank_area_m2+tot.unknown_possible_flank_area_m2)/1e12,area_difference_from_2013_million_km2=tot.eligible_area_m2/1e12-150.9,power_difference_from_2013_TW=tot.power_all_W/1e12-8))
globaldf=pd.DataFrame(globalrows);globaldf.to_csv(O/'global_geological_support.csv',index=False,mode='x')
refs=pd.read_csv(R/C['original_reference']);primary=basins[basins['mode']=='exact_cell_intersection'].set_index('basin');geoaudit=[]
for reference in ['low','central','high']:
 f=refs[(refs.reference==reference)&(refs.outlet_class=='low_temperature_flank')].set_index('basin').loc[names].baseline_power_TW
 for region in names+['Global']:
  F=float(f.sum()) if region=='Global' else float(f.loc[region]);area=float(primary.eligible_flank_area_m2.sum()) if region=='Global' else float(primary.loc[region,'eligible_flank_area_m2']);unknown=float(primary.unknown_possible_flank_area_m2.sum()) if region=='Global' else float(primary.loc[region,'unknown_possible_flank_area_m2']);req=F*1e12/.04
  geoaudit.append(dict(reference=reference,basin=region,ordinary_flank_TW=F,required_collection_million_km2=req/1e12,known_eligible_flank_million_km2=area/1e12,known_area_coverage_ratio=area/req,missing_sediment_permissive_coverage_ratio=(area+unknown)/req,known_area_sufficient=bool(area>=req),area_deficit_million_km2=max(req-area,0)/1e12,scope='conditional nonoverlap collection-area test; age/sediment criterion not an empirical hard exclusion of all circulation'))
pd.DataFrame(geoaudit).to_csv(O/'reference_geometry_audit.csv',index=False,mode='x')
# Derived geographic power shares are diagnostics, never fitted to modern residuals.
weights=[]
for name in names:
 p=primary.loc[name];weights.append(dict(basin=name,eligible_flank_area_share_percent=100*p.eligible_flank_area_m2/primary.eligible_flank_area_m2.sum(),background_flank_power_share_percent=100*p.power_flank_W/primary.power_flank_W.sum(),retained_flank_weight_percent=100*refs[(refs.reference=='high')&(refs.outlet_class=='low_temperature_flank')].set_index('basin').loc[name,'baseline_power_TW']/refs[(refs.reference=='high')&(refs.outlet_class=='low_temperature_flank')].baseline_power_TW.sum()))
pd.DataFrame(weights).to_csv(O/'geological_vs_retained_weights.csv',index=False,mode='x')
for key in mapkeys:assert abs(maps[key].sum()-primary[key].sum())<=1e-10*max(abs(primary[key].sum()),1)
out=xr.Dataset({k:(('lat','lon'),v) for k,v in maps.items()},coords={'lat':np.arange(-89,90,2),'lon':np.arange(-179,180,2)})
for name in out.data_vars:out[name].attrs['units']='W' if name.endswith('_W') else 'm2'
out.attrs.update(scope='present-day geological background and eligibility, not modern extra heat B',age_source='Seton2020GTS2012 6arcmin',sediment_source='GlobSedv3PANGAEA982339 5arcmin',area_operator='exact spherical intersections of piecewise-constant nodal control areas; conservative2degreeaggregation',no_scalar_integral_rescaling='true',age_units='Ma',sediment_units='m')
with (O/'geological_background_maps.nc').open('xb') as f:f.write(out.to_netcdf(engine='scipy'))
check={'sphere_area_relative_error':abs(sphere-4*np.pi*rad**2)/(4*np.pi*rad**2),'common_refinement_shape':[len(ym),len(xm)],'mapped_integrals_match_basin_integrals':True,'prior_integrals_forced':False,'maximum_nearest_native_assignment_distance_km':max(r['maximum_assignment_distance_km'] for r in records),'age_grid_minimum_Ma':float(np.nanmin(av)),'published_Hasterok_area_million_km2':150.9,'published_Hasterok_background_TW':8,'interpretation':'eligibility is a geological screening mask; no independent hydraulic connectivity, inventory temperature or modern anomaly inferred','input_sha256':{key:hashlib.sha256((R/C[key]).read_bytes()).hexdigest() for key in ['age_grid','sediment_grid','native_geometry','original_reference']}}
with (O/'spatial_verification.json').open('x') as f:json.dump(check,f,indent=2)
print(globaldf.to_string(index=False));print(primary[['eligible_flank_million_km2','background_flank_TW']].to_string());print(pd.DataFrame(geoaudit).to_string(index=False));print(json.dumps(check,indent=2))

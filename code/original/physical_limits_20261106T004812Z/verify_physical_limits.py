"""Independent explicit impulse propagation and envelope checks."""
from pathlib import Path
import argparse,calendar,json
import numpy as np,pandas as pd
ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--source-package',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
R=a.results;S=a.source_package;c=json.loads((S/'Source_Data/frozen_source_config.json').read_text());tc=json.loads((S/'regional_transfer_config.json').read_text())
vectors=np.array(c['normal_reference_vectors']);z=c['settings']['S2'];den=1+(z['H']-1)*c['normal_enhanced_fraction'];rhi=(1+(z['H']-1)*z['fmax'])/den;rlo=(1+(z['H']-1)*z['fmin'])/den
Y=np.repeat(np.arange(1850,2026),12);M=np.tile(np.arange(1,13),176);days=np.array([calendar.monthrange(int(y),int(m))[1] for y,m in zip(Y,M)]);duration=days/365.25;end=np.cumsum(duration);start=end-duration;units=days*86400e-9;modern=(Y>=1955)&(Y<=2023)
weights=pd.read_csv(R/'response_kernels.csv');outer=pd.read_csv(R/'monthly_outer_bounds.csv');req=pd.read_csv(R/'target_requirements.csv');records=pd.read_csv(R/'support_attribution.csv');delrows=pd.read_csv(R/'deletion_diagnostics.csv');h=pd.read_csv(S/'Source_Data/certified_histories.csv');mat=json.loads((R/'analysis_record.json').read_text())['matrix']
max_kernel=0.;max_outer=0.;max_forward=0.;max_grid_excess=-1e99;max_requirement=0.;kernel_data={}
for label,tau in [('direct',None),('loss_100yr',100.),('loss_30yr',30.)]:
    g=np.zeros(len(Y))
    # Explicit impulse survival into each endpoint-month mean, no backward adjoint.
    for i in np.flatnonzero((Y==1955)|(Y==2023)):
        sign=(-1 if Y[i]==1955 else 1)/12
        if tau is None:g[:i]+=sign;g[i]+=.5*sign
        else:
            bi=tau/duration[i]*(-np.expm1(-duration[i]/tau))
            bj=tau/duration[:i]*(-np.expm1(-duration[:i]/tau))
            g[:i]+=sign*bj*bi*np.exp(-(start[i]-end[:i])/tau)
            g[i]+=sign*tau/duration[i]*(1-bi)
    got=weights[weights.loss==label].energy_weight.to_numpy();max_kernel=max(max_kernel,float(np.max(abs(g[modern]-got))))
    kernel_data[label]=(g,g*units)
    # Independent forward recurrence for arbitrary signed, not just epoch sources.
    rng=np.random.default_rng(73519)
    for n in range(10):
        q=np.zeros(len(Y));q[modern]=rng.uniform(-8,6,modern.sum());e=q*units;state=0.;hh=[]
        for j in range(len(Y)):
            if tau is None:hh.append(state+e[j]/2);state+=e[j]
            else:
                A=np.exp(-duration[j]/tau); b=tau/duration[j]*(1-A)
                hh.append(b*state+tau/duration[j]*(1-b)*e[j]);state=A*state+b*e[j]
        hh=np.array(hh);b=hh[Y==2023].mean()-hh[Y==1955].mean();max_forward=max(max_forward,abs(b-g@e))

def L(v):return max(-sum(v[:3]),rlo*sum(v[:2])-v[2]-v[5])
def U(v,f):return (rhi-1)*sum(v[:2])+.5*v[2] if f=='S2_P3' else {'global_C8':8,'global_C10':10}[f]+rhi*v[1]-sum(v[:3])
caps={r['allocation']:r['ordinary_flank_cap_TW'] for r in tc['geometry_cases']}
for row in outer.itertuples():
    g,k=kernel_data[row.loss];v=vectors[row.reference_segment]+row.reference_interpolation*np.diff(vectors,axis=0)[row.reference_segment]
    q=np.zeros(len(Y));q[modern]=np.where(k[modern]>=0,U(v,row.family),L(v));value=k@q
    max_outer=max(max_outer,abs(value-row.monthly_timing_outer_ZJ))
    # Check every affine segment boundary and donor branch, independently found,
    # plus a dense interior set. The grid is a numerical check, not the proof.
    check=[]
    for seg in range(2):
        v0,v1=vectors[seg:seg+2];dv=v1-v0;ts=[0.,1.]
        for f in [lambda v:v[2]-caps[row.allocation],lambda v:v[5]-(1+rlo)*sum(v[:2])]+([] if row.family=='S2_P3' else [lambda v:v[0]+v[2]-{'global_C8':8,'global_C10':10}[row.family]]):
            x0,x1=f(v0),f(v1)
            if abs(x1-x0)>1e-14:
                root=-x0/(x1-x0)
                if 0<=root<=1:ts.append(root)
        ts+=list(np.linspace(0,1,1001))
        for t in ts:
            vv=v0+t*dv
            if vv[2]>caps[row.allocation]+1e-10:continue
            if row.family!='S2_P3' and vv[0]+vv[2]>{'global_C8':8,'global_C10':10}[row.family]+1e-10:continue
            check.append(np.maximum(k[modern],0).sum()*U(vv,row.family)+np.minimum(k[modern],0).sum()*L(vv))
    max_grid_excess=max(max_grid_excess,max(check)-value)
    assert abs(max(check)-value)<1e-9
for row in req.itertuples():
    g,k=kernel_data[row.loss];kp=np.maximum(k[modern],0).sum();kn=np.minimum(k[modern],0).sum()
    # The relaxed peak thresholds attain their envelope at the low reference.
    v=vectors[0];C=row.monthly_necessary_axis_flank_peak_TW
    B=kp*(C+rhi*v[1]-sum(v[:3]))+kn*L(v)
    max_requirement=max(max_requirement,abs(B-row.target_ZJ))
    assert abs(row.necessary_positive_extra_input_ZJ*g[modern].max()+row.maximum_baseline_recovery_allowance_ZJ-row.target_ZJ)<1e-9
for row in records.itertuples():
    hh=h[(h.allocation==row.allocation)&(h.family==row.family)&(h.loss==row.loss)]
    assert len(hh)==row.count
    if len(hh):assert abs((100*hh.B_1955_2023_ZJ/hh.O_1955_2023_ZJ).max()-row.joint_upper_percent)<1e-10
for row in delrows.itertuples():
    hh=h[(h.allocation==row.allocation)&(h.family==row.family)&(h.loss==row.loss)]
    value=int(row.deleted_value) if row.deleted_field in ['candidate_index','O_member'] else row.deleted_value
    remain=hh[hh[row.deleted_field]!=value];assert len(remain)==row.count_remaining
    if len(remain):assert abs((100*remain.B_1955_2023_ZJ/remain.O_1955_2023_ZJ).max()-row.upper_percent)<1e-10
assert max_kernel<1e-10 and max_forward<1e-9 and max_outer<1e-9 and max_requirement<1e-9
result=dict(passed=True,explicit_impulse_vs_adjoint_max_error=max_kernel,random_signed_forward_checks=30,random_signed_forward_max_error_ZJ=max_forward,outer_witness_checks=len(outer),outer_max_error_ZJ=max_outer,reference_boundary_and_grid_max_excess_ZJ=max_grid_excess,power_requirement_checks=len(req),power_requirement_max_error_ZJ=max_requirement,support_rows=len(records),deletion_rows=len(delrows),notes='Numerical verification of analytical envelopes, not empirical validation. Dense interior checks supplement the piecewise-affine breakpoint proof.')
with a.output.open('x') as f:json.dump(result,f,indent=2)
print(json.dumps(result,indent=2))

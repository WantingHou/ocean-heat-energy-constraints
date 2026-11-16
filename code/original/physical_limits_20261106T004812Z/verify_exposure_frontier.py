"""Verify the duration frontier with its independent scalar dual bound."""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();R=a.results
k=pd.read_csv(R/'response_kernels.csv');refs=pd.read_csv(R/'reference_breakpoints.csv');f=pd.read_csv(R/'exposure/exposure_frontier.csv');inv=pd.read_csv(R/'exposure/exposure_thresholds.csv');maxerr=0.;inv_err=0.;checked=0
for loss in k.loss.unique():
    x=k[k.loss==loss];g=x.energy_weight.to_numpy();u=x.seconds.to_numpy()*1e-9;neg=np.minimum(g*u,0).sum();gp=np.maximum(g,0);lam=np.unique(np.r_[0,gp]);penalty=(np.maximum(gp[None,:]-lam[:,None],0)*u).sum(axis=1)
    def W(D):return float(np.min(lam*D*31557600e-9+penalty))
    for alloc in refs.allocation.unique():
        vr=refs[refs.allocation==alloc]
        for family in f.family.unique():
            Us=vr.S2_P3_U_TW.to_numpy() if family=='S2_P3' else {'global_C8':8.,'global_C10':10.}[family]-vr.gross_power_offset_TW.to_numpy();Ls=vr.L_TW.to_numpy()
            def bound(D):return float(np.max(Us*W(D)+Ls*neg))
            rows=f[(f.loss==loss)&(f.allocation==alloc)&(f.family==family)]
            for row in rows.itertuples():maxerr=max(maxerr,abs(bound(row.equivalent_full_power_years)-row.outer_response_ZJ));checked+=1
            for row in inv[(inv.loss==loss)&(inv.allocation==alloc)&(inv.family==family)].itertuples():
                if row.feasible_in_relaxed_power_box:
                    D=row.minimum_equivalent_full_power_years;inv_err=max(inv_err,abs(bound(D)-row.target_ZJ));assert bound(D-1e-5)<row.target_ZJ
                else:assert bound(100)<row.target_ZJ
assert maxerr<1e-9 and inv_err<1e-9
result=dict(passed=True,independent_method='dual minimization over observed kernel breakpoints; no source optimizer',frontier_points_checked=checked,max_frontier_error_ZJ=maxerr,inversion_rows_checked=len(inv),max_inversion_error_ZJ=inv_err)
with a.output.open('x') as z:json.dump(result,z,indent=2)
print(json.dumps(result,indent=2))

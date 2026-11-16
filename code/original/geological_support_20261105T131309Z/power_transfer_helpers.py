def states(power,tau):
 e=np.asarray(power)*unit
 if tau is None:return (np.cumsum(e)-e/2).reshape(176,12).mean(axis=1)
 zz=dt/tau;aa=np.exp(-zz);bb=-np.expm1(-zz)/zz;dd=(1-bb)/zz;h=0.;out=[]
 for j,ee in enumerate(e):out.append(bb[j]*h+dd[j]*ee);h=aa[j]*h+bb[j]*ee
 return np.array(out).reshape(176,12).mean(axis=1)

def contrast(h):return np.array([h[hi-1850]-h[lo-1850] for lo,hi in windows])

def intersect(interval,f0,f1):
 if interval is None:return None
 lo,hi=interval;slope=f1-f0
 if abs(slope)<1e-13:return interval if f0<=1e-10 else None
 boundary=-f0/slope
 if slope>0:hi=min(hi,boundary)
 else:lo=max(lo,boundary)
 if lo>hi+1e-10:return None
 if lo>hi:lo=hi=(lo+hi)/2
 return (max(0.,lo),min(1.,hi))

def certify(seg,x,p,C,loss,required=None):
 global max_violation,max_closure,worst_stock
 v=vectors[seg]+x*(vectors[seg+1]-vectors[seg]);g=builders[seg];xx=np.r_[x,np.asarray(p).ravel()];res=g['A']@xx-g['b'];viol=max(float(res.max()),float(-xx.min()),float(x-1),float((p[:,0]+p[:,2]-C).max()),float(v[0]+v[2]-C));assert viol<1e-8,('constraint',viol);max_violation=max(max_violation,viol)
 cy=v[5]-v[0]-v[1];fy=np.full(len(MY),v[5]);ff=np.full(len(MY),v[2]);power=np.zeros(len(MY))
 for j,m in enumerate(masks[1:]):
  pa,pv,pf,ry,rf=p[j];fy[m]=pa+pv+cy-ry;ff[m]=pf-rf;power[m]=pa+pv+pf-ry-rf-sum(v[:3])
 initial_y=v[4]*50*31557600e-9;reserved=v[4]*ufull+initial_y;N=v[2]/2e-6;per_m=N*1e6*3330*1171*60/1e21;initial_old=300*per_m-reserved
 ys=initial_y+np.cumsum((v[4]-fy)*unit);os=initial_old-np.cumsum(ff*unit);ns=initial_old-v[2]*np.cumsum(unit);stock=min(float(ys.min()),float(os.min()),float(ns.min()));assert stock>=-1e-8;worst_stock=min(worst_stock,stock)
 energy=g['extra']@xx+g['extra_c'];assert np.max(abs(energy-power*unit))<1e-10
 h=states(power,operators[loss]);b=contrast(h)
 if required is not None:
  err=float(np.max(abs(b-required)));assert err<1e-8;max_closure=max(max_closure,err)
 h_required=(reserved+max(float(np.cumsum(ff*unit).max()),float(v[2]*ufull)))/per_m
 return dict(variables_json=json.dumps(xx.tolist()),reference_segment=seg,reference_interpolation=float(x),minimum_young_stock_ZJ=float(ys.min()),minimum_old_stock_ZJ=float(os.min()),minimum_ordinary_old_stock_ZJ=float(ns.min()),old_remaining_2023_ZJ=float(os[(MY==2023)&(MM==12)][0]),required_connected_thickness_m=h_required,max_native_flank_TW=float(max(ff)),maximum_gross_axis_flank_TW=float(max(p[:,0]+p[:,2])),maximum_equivalent_flank_MW=float(max(p[:,2])/N*1e6)),b

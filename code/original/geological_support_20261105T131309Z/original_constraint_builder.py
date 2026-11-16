def build(st,segment,h_override=None):
 z=settings[st];v0=vectors[segment];dv=vectors[segment+1]-v0;h=z['h_m'] if h_override is None else h_override
 # Affine background vector: row * x + constant.
 bg=np.zeros((6,n));bg[:,0]=dv;bc=v0.copy();cy=bg[5]-bg[0]-bg[1];cyc=bc[5]-bc[0]-bc[1]
 H=z['H'];ratio_lo=(1+(H-1)*z['fmin'])/(1+(H-1)*10/33);ratio_hi=(1+(H-1)*z['fmax'])/(1+(H-1)*10/33)
 A=[];b=[];names=[]
 def le(a,c,name):A.append(a);b.append(c);names.append(name)
 def pos(k):r=np.zeros(n);r[k]=1;return r
 for j in range(3):
  k=1+5*j
  for c,lo,hi in [(0,ratio_lo,ratio_hi),(1,ratio_lo,ratio_hi),(2,.5,z['flank_max_ratio'])]:
   le(pos(k+c)-hi*bg[c],hi*bc[c],f'power_high_{j}_{c}');le(-pos(k+c)+lo*bg[c],-lo*bc[c],f'power_low_{j}_{c}')
  ry=pos(k+3);rf=pos(k+4)
  le(ry-cy,cyc,f'young_donor_{j}');le(rf-bg[3]+cy,bc[3]-cyc,f'old_donor_{j}');le(rf-3.5*bg[2],3.5*bc[2],f'catchment_donor_{j}');le(rf-pos(k+2),0,f'flank_route_{j}');le(ry-pos(k)-pos(k+1),0,f'young_route_{j}')
 # All possible young cohorts reserved before accessing old connected rock; no full-plate maximization.
 ecoef=1e6*h*3330*1171*60/1e21/(2e-6)
 yreserve_coef=z['young_years']*31557600*1e-9
 young_initial=bg[4]*yreserve_coef;young_initial_c=bc[4]*yreserve_coef
 old_initial=ecoef*bg[2]-(ufull+yreserve_coef)*bg[4];old_initial_c=ecoef*bc[2]-(ufull+yreserve_coef)*bc[4]
 # Native flank release is actual flank minus already-counted conductive routing.
 ycum=u0*bg[5];ycumc=u0*bc[5];fcum=u0*bg[2];fcumc=u0*bc[2];elapsed=u0
 for j,uj in enumerate(modern_units):
  k=1+5*j;ycum=ycum+uj*(pos(k)+pos(k+1)-pos(k+3)+cy);ycumc+=uj*cyc;fcum=fcum+uj*(pos(k+2)-pos(k+4));elapsed+=uj
  le(ycum-young_initial-elapsed*bg[4],young_initial_c+elapsed*bc[4]-ycumc,f'young_stock_{j}')
  le(fcum-old_initial,old_initial_c-fcumc,f'old_stock_{j}')
 le(ufull*bg[2]-old_initial,old_initial_c-ufull*bc[2],'background_old_stock')
 le(bg[5]-bg[4],bc[4]-bc[5],'background_young_supply')
 # Calendar source matrices for whole model support; background and actual share the same reference.
 normal=np.tile(bg[:4].sum(axis=0),(len(MY),1))*unit[:,None];normal_c=np.full(len(MY),bc[:4].sum())*unit
 actual=normal.copy();actual_c=normal_c.copy()
 for j,mask in enumerate(masks[1:]):
  k=1+5*j;power=bg[3]+pos(k)+pos(k+1)+pos(k+2)-pos(k+3)-pos(k+4)
  actual[mask]=unit[mask,None]*power;actual_c[mask]=unit[mask]*bc[3]
 return dict(A=np.array(A),b=np.array(b),names=names,bounds=[(0,1)]+[(0,None)]*15,bg=bg,bc=bc,actual=actual,actual_c=actual_c,extra=actual-normal,extra_c=actual_c-normal_c,old_initial=old_initial,old_initial_c=old_initial_c,young_initial=young_initial,young_initial_c=young_initial_c,cy=cy,cyc=cyc,segment=segment,st=st,h=h)

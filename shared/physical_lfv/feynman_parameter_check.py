"""Independent Euclidean/below-cut Feynman-parameter quadrature check of hard constants and finite-q expansion."""
from pathlib import Path
import numpy as np,json
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss
import physical_threebody as p
zz,ww=leggauss(96);zs=(zz+1)/2;ws=ww/2
rows=[];M=5000.
for mk in [p.ME,p.MMU,p.MTAU]:
 for r in [.05,1.,3.9]:
  q2=r*mk*mk;tmax=np.log(M*M/mk**2)+25
  ff=hh=0.
  for z,weight in zip(zs,ws):
   w=z*(1-z)
   def kernels(t,which):
    a=np.exp(-t);u=-np.expm1(-t);xy=u*u*w;d0=a*M*M+u*mk*mk;dq=d0-xy*q2
    if which=='f':return M*M*u*a*(xy/dq-np.log1p(-xy*q2/d0)/q2+mk*mk*xy/(dq*d0))
    return M*M*u*u*a/dq
   ff+=weight*quad(lambda t:kernels(t,'f'),0,tmax,epsabs=1e-10,epsrel=1e-10)[0]
   hh+=weight*quad(lambda t:kernels(t,'h'),0,tmax,epsabs=1e-10,epsrel=1e-10)[0]
  I,J=p.logmoments(np.array([r]));fl=(np.log(M*M/mk**2)-4/3)/3-2*J[0].real;hl=np.log(M*M/mk**2)-1.5-I[0].real
  rows.append({'mk':mk,'r':r,'f_FP':ff,'f_leading':fl,'f_relative_difference':ff/fl-1,'H_FP':hh,'H_leading':hl,'H_relative_difference':hh/hl-1})
assert max(abs(t['f_relative_difference']) for t in rows)<1e-5
assert max(abs(t['H_relative_difference']) for t in rows)<1e-6
out={'cases':rows,'note':'Exact Feynman denominators omit external masses consistently; below-cut quadrature validates finite constants and q dependence independently of the phase-space code.'}
Path(__file__).with_name('feynman_check_v10.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS: independent finite-q Feynman-parameter check; max f and H relative errors',max(abs(t['f_relative_difference']) for t in rows),max(abs(t['H_relative_difference']) for t in rows))

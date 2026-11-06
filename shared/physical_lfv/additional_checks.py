import json,sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import physical_threebody as p
R=Path(__file__).parent;H=np.load(R/'gram_320.npy');Q=p.fields();rd=p.old.ALPHA/(3*np.pi)*(np.log(p.MMU**2/p.ME**2)-11/4)
def coeff(d,M=5000):return p.coefficients(p.fields(M),M*np.tile([1+d,1-d],4))
roots={}
for name,target in [('R_gt_10RD',10*rd),('MEG_90CL_at5TeV',1.5e-13)]:
 f=(lambda d:p.rate(coeff(d),H)/p.radiative(coeff(d))-target) if name.startswith('R') else (lambda d:p.radiative(coeff(d))-target)
 roots[name]=[brentq(f,-.03,0),brentq(f,0,.03)]
errors=[]
for r in [1e-7,.05,1.,3.9,4.,4.1,10.,100.,10000.]:
 rootsz=([.5] if r==4 else []) if r<=4 else [(1-np.sqrt(1-4/r))/2,(1+np.sqrt(1-4/r))/2]
 i,j=p.logmoments(np.array([r]));ivals=[]
 for weight in [lambda z:1,lambda z:z*(1-z)]:
  re=quad(lambda z:weight(z)*np.log(abs(1-r*z*(1-z))),0,1,points=rootsz,epsabs=1e-10)[0]
  im=0 if r<=4 else -np.pi*quad(weight,rootsz[0],rootsz[1],epsabs=1e-13)[0]
  ivals.append(re+1j*im)
 errors.append(max(abs(i[0]-ivals[0]),abs(j[0]-ivals[1])))
# Retained exact x-dependent on-shell functions check leading heavy mass expansion.
errs=[]
for M in [3500.,5000.,6500.]:
 for mk in [p.ME,p.MMU,p.MTAU]:
  x=mk*mk/M**2;errs.extend([abs(p.old.Ffun(x)/(1/6)-1),abs(p.old.Hfun(x)/(np.log(1/x)-1.5)-1)])
out={'common_splitting_roots':roots,'finite_q_log_moment_max_absolute_error_against_independent_quadrature':max(errors),'on_shell_scalar_functions_max_relative_heavy_expansion_error':max(errs),'rescaled_split_counterexample':{}}
for M in [10000,11000]:
 c=coeff(.05,M);out['rescaled_split_counterexample'][str(M)]={'Beg':p.radiative(c),'B3e':p.rate(c,H),'R':p.rate(c,H)/p.radiative(c)}
(R/'additional_checks_v10.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

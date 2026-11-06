"""Declared scans; no probability over UV theories or electroweak-viable fits."""
import numpy as np,json,time
from pathlib import Path
from scipy.optimize import minimize,brentq
import physical_threebody as p
ROOT=Path(__file__).parent
H=np.load(ROOT/'gram_320.npy');Q=p.fields();m0=np.full(8,5000.);c0=p.coefficients(Q,m0);B0=p.rate(c0,H);G0=p.radiative(c0);RD=p.old.ALPHA/(3*np.pi)*(np.log(p.MMU**2/p.ME**2)-11/4)
def report(c):
 b=p.rate(c,H);g=p.radiative(c);return {'B3e':b,'Beg':g,'R':b/g,'R_over_RD':b/g/RD}
def rotation(n,i,j,a):
 O=np.eye(n);O[i,i]=O[j,j]=np.cos(a);O[i,j]=np.sin(a);O[j,i]=-np.sin(a);return O
def rotate(Q,O):return np.einsum('bi,ijk->bjk',O,Q)
def minmax(rows):
 a=np.array(rows);return {'n':len(a),'B3e_minmax':a[:,0].min(axis=0).tolist() if False else [float(a[:,0].min()),float(a[:,0].max())],'Beg_minmax':[float(a[:,1].min()),float(a[:,1].max())],'R_minmax':[float(a[:,2].min()),float(a[:,2].max())],'R_percentiles_5_50_95':np.percentile(a[:,2],[5,50,95]).tolist(),'fraction_R_gt_10RD':float(np.mean(a[:,2]>10*RD))}
def trip(c):return [p.rate(c,H),p.radiative(c),p.rate(c,H)/p.radiative(c)]
def matching(eta):
    # Independent LEFT basis: scalar RR/LL, four vectors, two dipoles.
    # Mixed scalar contacts obey S_RL=-V_RL/2 and S_LR=-V_LR/2 after antisymmetrization.
    ids=np.array([0,1,4,5,6,7,8,9]);nv=len(ids)
    scales=eta*np.r_[np.full(6,max(abs(c0[:4]))),np.full(2,max(abs(c0[8:10])))]
    def make(x):
        c=c0.copy();c[ids]+=scales*(x[:nv]+1j*x[nv:]);return c
    def objective(x,what,sign):
        c=make(x);b=p.rate(c,H)/B0;g=p.radiative(c)/G0
        return sign*np.log(b if what=='B' else b/g)
    def jac(x,what,sign):
        c=make(x);hc=H@c;bv=p.rate(c,H)
        gb=np.r_[2*scales*hc[ids].real,2*scales*hc[ids].imag]/bv
        if what=='R':
            gd=np.zeros(2*nv);norm=np.sum(abs(c[8:10])**2)
            gd[nv-2:nv]=2*scales[-2:]*c[8:10].real/norm
            gd[-2:]=2*scales[-2:]*c[8:10].imag/norm
            gb-=gd
        return sign*gb
    cons={'type':'ineq','fun':lambda x:1-x[:nv]**2-x[nv:]**2,
          'jac':lambda x:np.c_[-np.diag(2*x[:nv]),-np.diag(2*x[nv:])]}

    out={}
    for what in ['B','R']:
        for sign in [1,-1]:
            sols=[minimize(objective,np.zeros(2*nv),args=(what,sign),method='SLSQP',jac=jac,constraints=cons,bounds=[(-1,1)]*(2*nv),options={'ftol':1e-9,'maxiter':1000})]
            for seed in range(2):
                x=np.random.default_rng(900+seed).normal(size=2*nv);x=x.reshape(2,nv);x*=.5/np.maximum(np.sqrt(np.sum(x*x,axis=0)),1)[None,:]
                sols.append(minimize(objective,x.ravel(),args=(what,sign),method='SLSQP',jac=jac,constraints=cons,bounds=[(-1,1)]*(2*nv),options={'ftol':1e-9,'maxiter':1000}))
            good=[r for r in sols if r.success and np.min(cons['fun'](r.x))>-1e-7]
            if not good:raise RuntimeError('No feasible converged matching optimum')
            best=min(good,key=lambda r:r.fun)
            out[f'{what}_{"min" if sign==1 else "max"}']={'normalized':float(np.exp(sign*best.fun)),'physical':float(np.exp(sign*best.fun)*(B0 if what=='B' else B0/G0)),'success':bool(best.success),'disk_violation':float(max(0,-np.min(cons['fun'](best.x)))),'shifts_normalized':best.x.tolist(),'basis_indices':ids.tolist()}
    return out

def main():
 out={'scope':'LO scalar tree plus finite one-loop photon amplitude at leading inverse pole mass squared. Independent UV local Wilson terms varied explicitly. Not inclusive QED NLO.','baseline':report(c0),'RD':RD,'validation':json.loads((ROOT/'validation_320.json').read_text())}
 cd=c0.copy();cd[4:8]=0;cd[10:]=0;out['baseline_contact_dipole_exact']=report(cd)
 cq=c0.copy();cq[10:]=0;out['q0_photon_diagnostic']=report(cq)
 rng=np.random.default_rng(271828);OO=np.linalg.qr(rng.normal(size=(8,8)))[0]
 ci=p.coefficients(rotate(Q,OO),m0);out['exact_degeneracy_rotation_max_relative_coefficient_change']=float(np.max(abs(ci-c0))/np.max(abs(c0)))
 out['common_pair_splitting']={}
 for d in [-.3,-.2,-.1,-.05,0,.05,.1,.2,.3]:out['common_pair_splitting'][str(d)]=report(p.coefficients(Q,5000*np.tile([1+d,1-d],4)))
 out['label2_pair_splitting']={}
 for d in [-.2,-.1,-.05,0,.05,.1,.2]:
  m=m0.copy();m[2:4]=5000*np.array([1+d,1-d]);out['label2_pair_splitting'][str(d)]=report(p.coefficients(Q,m))
 out['restricted_rotation']={}
 for spread in [.1,.2]:
  m=m0.copy();m[[0,2]]=5000*np.array([1-spread,1+spread]);vals=[];best=None
  for a in np.linspace(-np.pi/2,np.pi/2,721):
   cc=p.coefficients(rotate(Q,rotation(8,0,2,a)),m);v=trip(cc);vals.append(v)
   if best is None or v[2]<best[0]:best=(v[2],a,cc)
  out['restricted_rotation'][str(spread)]={**minmax(vals),'minimum_angle_radians':float(best[1]),'minimum_witness':report(best[2])}
 out['independent_poles']={};out['mixed_poles']={}
 for spread in [.1,.2]:
  nomix=[];mixed=[];best=None
  for _ in range(3000):
   m=5000*(1+rng.uniform(-spread,spread,8));nomix.append(trip(p.coefficients(Q,m)))
   O=np.linalg.qr(rng.normal(size=(8,8)))[0];cc=p.coefficients(rotate(Q,O),m);v=trip(cc);mixed.append(v)
   if best is None or v[2]<best[0]:best=(v[2],O,m,cc)
  out['independent_poles'][str(spread)]=minmax(nomix);out['mixed_poles'][str(spread)]={**minmax(mixed),'minimum_witness':{'O':best[1].tolist(),'masses_GeV':best[2].tolist(),**report(best[3])}}
  np.save(ROOT/f'pole_scan_{int(spread*100)}.npy',np.stack([nomix,mixed]))
 out['matching']={str(eta):matching(eta) for eta in [1/(16*np.pi**2),.1,.3,.5]}
 # Local tree-sized coefficient freedom destroys any uniform lower ratio: cancel all C,V with explicit matching and retain dipole.
 dc=c0.copy();dc[:8]=0;out['unrestricted_local_contact_cancellation_diagnostic']=report(dc)
 # Common-scale reach, finite momentum amplitude included.
 def bs(M):
  cc=p.coefficients(p.fields(M*1000),np.full(8,M*1000));return p.rate(cc,H)
 out['common_scale_reach_TeV']={label:float(brentq(lambda M:bs(M)-target,2,20)) for label,target in [('Mu3e_I_SES',2e-15),('Mu3e_I_zero_background_90CL',-np.log(.1)*2e-15),('Mu3e_declared_1e-16_target',1e-16)]}
 out['scalar_box_norm_budget_GeV_minus2']={'sum_spectral_norms_squared':float(sum(np.linalg.norm(q,2)**2 for q in Q)),'conservative_dimensional_envelope':float(sum(np.linalg.norm(q,2)**2 for q in Q)**2/(16*np.pi**2*min(m0)**2)),'note':'Power-counting envelope, not a calculated box Wilson coefficient or a rigorous diagram-count bound.'}
 (ROOT/'physical_results_v10.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k not in ['matching','mixed_poles']},indent=2))
 print('matching extrema', {k:{j:x['physical'] for j,x in v.items()} for k,v in out['matching'].items()})
 print('mixed extrema',{k:{j:x for j,x in v.items() if j!='minimum_witness'} for k,v in out['mixed_poles'].items()})
if __name__=='__main__':main()

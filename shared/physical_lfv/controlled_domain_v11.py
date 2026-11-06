"""Conditional amplitude enclosure; neutral-Yukawa box majorant, not UV matching.
Every budget is deterministic. No fit likelihood or inclusive NLO is implied.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import physical_threebody as p
import higgs_mixing_diagnostic as hm

R=Path(__file__).parent
H=np.load(R/'gram_320.npy');Q=p.fields();masses=np.full(8,5000.)
c=p.coefficients(Q,masses);B=p.rate(c,H);Bg=p.radiative(c)
N=np.sqrt(B);dn=np.linalg.norm(c[8:10]);C=max(abs(c[:4]));D=max(abs(c[8:10]))
ix=[0,1,4,5,6,7];eta=1/(16*np.pi**2)
# Entrywise sum includes BOTH chiral vertices, without relying on cancellations.
A=abs(Q)+abs(Q.swapaxes(1,2).conj())
F=sum((a@b)[0,1]*(b@a)[0,0] for a in A for b in A)
W=16*F/(16*np.pi**2*min(masses)**2)
box=W*sum(np.sqrt(H.diagonal().real[ix]))/N
local=eta*C*sum(np.sqrt(H.diagonal().real[ix]))/N
ld3=eta*D*sum(np.sqrt(H.diagonal().real[8:10]))/N
ldg=eta*D*np.sqrt(2)/dn

def mixed(i,t):
    qq=np.concatenate([Q,hm.qh[None]],axis=0);O=np.eye(9)
    O[i,i]=O[8,8]=np.cos(t);O[8,i]=np.sin(t);O[i,8]=-np.sin(t)
    q=np.einsum('bi,ijk->bjk',O,qq)
    return p.coefficients(q,np.r_[masses,125.])

cases={};limits=[];neutral3=neutralg=0.;higgs3=0.
for i,label in [(0,'s1'),(2,'s2'),(4,'s4'),(6,'s5')]:
    # Sum absolute top and W contributions: no accidental t/W cancellation used.
    ft,_=hm.fg(hm.MT**2/125**2);fth,_=hm.fg(hm.MT**2/5000**2)
    top=-6*(2/3)**2*p.old.V/hm.MT*hm.YTT*(ft-fth)
    total=hm.baz(125)-hm.baz(5000);weak=total-top
    slope=p.E2*hm.kappa*(abs(top)+abs(weak))/(4*np.pi**2)*np.linalg.norm([Q[i,0,1],Q[i,1,0].conjugate()])
    t=.0075*dn/slope;limits.append(t)
    delta=mixed(i,t)-c
    # Triangle bound in both norms, including every finite-photon coefficient.
    n3=sum(abs(delta[j])*np.sqrt(max(H[j,j].real,0)) for j in range(22))/N
    ng=np.linalg.norm(delta[8:10])/dn
    neutral3+=n3;neutralg+=ng
    dvec=p.E2*hm.kappa*t*total/(4*np.pi**2)*np.array([Q[i,0,1],Q[i,1,0].conjugate()])
    h3=sum(abs(dvec[j])*np.sqrt(H[8+j,8+j].real) for j in range(2))/N;higgs3+=h3
    vev=np.sqrt(2)*p.old.EPS[[1,2,4,5][i//2]]*5000
    cases[label]={'angle_abs_cap_rad':t,'portal_abs_cap_if_only_quartic_Hchi':t*(5000**2-125**2)/(p.old.V*vev),'neutral_amplitude_3e_triangle_bound':n3,'neutral_dipole_relative_bound':ng,'topW_photon_3e_triangle_bound':h3,'topW_dipole_budget_per_channel':.0075}

# Residual norms contain QED corrections, omitted weak diagrams, parent thresholds,
# higher orders, and interference remainders. These are MATCHING REQUIREMENTS.
res3=.075;resg=.05
eps3=box+local+ld3+neutral3+higgs3+res3
epsg=ldg+neutralg+.03+resg
assert eps3<.25 and epsg<.1
R0=B/Bg
lo=R0*(.75/1.1)**2;hi=R0*(1.25/.9)**2
# Validate underlying Euclidean inequalities numerically, independent of flavour sums.
rng=np.random.default_rng(711)
for _ in range(2000):
    a,b=10**rng.uniform(-5,1,2);x=10**rng.uniform(-12,6)
    assert x/((x+a*a)*(x+b*b)) <= (1+1e-13)/x
    assert a*b/((x+a*a)*(x+b*b)) <= (1+1e-13)/(4*x)
loop_checks=[]
for u,v in [(1.,1.),(.8,1.2),(.95,1.05)]:
    integral=quad(lambda x:1/((x+u*u)*(x+v*v)),0,np.inf,epsabs=1e-12)[0]
    analytic=1/u**2 if u==v else np.log(u*u/(v*v))/(u*u-v*v)
    assert abs(integral/analytic-1)<1e-10
    loop_checks.append({'mass_ratios':[u,v],'integral':integral,'analytic':analytic})
rows=[]
for M in [3000,4000,5000,6000,7500,10000]:
    cc=p.coefficients(p.fields(M),np.full(8,M));rows.append({'M_TeV':M/1000,'Beg':p.radiative(cc),'B3e':p.rate(cc,H),'R':p.rate(cc,H)/p.radiative(cc)})
reach={}
for name,t in [('SINDRUM',1e-12),('Mu3e_SES',2e-15),('Mu3e_90CL',2.302585093*2e-15),('Mu3e_target',1e-16)]:
    reach[name]=brentq(lambda M:p.rate(p.coefficients(p.fields(M),np.full(8,M)),H)-t,2000,14000)/1000
out={'scope':'Neutral-scalar leading-dimension-six box majorant; conditional matching-domain enclosure. Residual budgets are imposed, not computed UV/NLO predictions.','SM_top_Ytt':hm.YTT,'central':{'Beg':Bg,'B3e':B,'R':R0},'box':{'entrywise_flavour_sum':F,'coefficient_majorant_GeV_minus2':W,'coefficient_majorant_over_Cstar':W/C,'amplitude_norm_majorant':box,'deliberate_overcount_factor':16},'amplitude_budget_3e':{'box':box,'local_six_coefficient_disks':local,'local_dipole_disks':ld3,'neutral_Higgs_admixture':neutral3,'topW_photon':higgs3,'residual_matching_requirement':res3,'sum':eps3,'declared_cap':.25},'amplitude_budget_dipole':{'local_disks':ldg,'neutral_Higgs_admixture':neutralg,'topW_photon':.03,'residual_matching_requirement':resg,'sum':epsg,'declared_cap':.1},'Higgs_angle_caps':cases,'enclosure':{'R_lower':lo,'R_upper':hi,'B3e_lower':B*.75**2,'B3e_upper':B*1.25**2,'Beg_lower':Bg*.9**2,'Beg_upper':Bg*1.1**2},'loop_integral_checks':loop_checks,'common_scale_rows':rows,'central_reach_TeV':reach}
(R/'controlled_domain_v11.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['loop_integral_checks','common_scale_rows']},indent=2))

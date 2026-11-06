import json,numpy as np
from pathlib import Path
import physical_threebody as p
R=Path(__file__).parent;H=np.load(R/'gram_320.npy');H160=np.load(R/'gram_160.npy');c=p.coefficients(p.fields(),np.full(8,5000.));v=p.validate(H)
assert abs(v['contact']/v['contact_analytic']-1)<.001
assert abs(v['dipole']/v['dipole_analytic']-1)<.0002
assert abs(v['radiative']/v['radiative_v9']-1)<2e-6
assert abs(p.rate(c,H)/p.rate(c,H160)-1)<1e-7
assert min(v['gram_contact_vector_diag_relative'])>-.001
assert np.linalg.eigvalsh(H)[0]/np.linalg.norm(H,2)>-1e-14
s,t,w,p0,p1,p2,p3=p.dalitz(8)
for pp,m,anti in [(p1,p.ME,False),(p3,p.ME,True),(p0,p.MMU,False)]:
 u=p.spinor(pp,m,anti);rhs=np.einsum('nv,vab->nab',pp*p.metric,p.gamma)+(-1 if anti else 1)*m*np.eye(4)
 assert np.max(abs(np.einsum('nai,nib->nab',u,p.bar(u))-rhs))<1e-14
A=p.kernels(p0,p1,p2,p3)
assert np.max(abs(A[:,2]+A[:,6]/2))/np.max(abs(A[:,2]))<1e-12
assert np.max(abs(A[:,3]+A[:,7]/2))/np.max(abs(A[:,3]))<1e-12
q=p2+p3;b=p.bar(p.spinor(p2,p.ME));vv=p.spinor(p3,p.ME,True)
ward=np.einsum('nm,nmab->nab',q*p.metric,np.stack([p.bil(b,g,vv) for g in p.gamma],axis=1));assert np.max(abs(ward))<1e-14
res=json.loads((R/'physical_results_v10.json').read_text());checks=json.loads((R/'additional_checks_v10.json').read_text())
assert checks['finite_q_log_moment_max_absolute_error_against_independent_quadrature']<2e-8
assert res['exact_degeneracy_rotation_max_relative_coefficient_change']<1e-14
for x in res['matching'].values():
 for y in x.values():assert y['success'] and y['disk_violation']<1e-7
print('PASS: spin sums, antiparticle completeness, pair-current Ward identity, Fierz identities, contact/dipole normalizations, finite-q quadrature, convergence, degenerate mixing and matching constraints.')

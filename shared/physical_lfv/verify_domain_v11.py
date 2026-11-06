"""Check box-reduction Fierz identities and the shared conditional enclosure."""
import json
from pathlib import Path
import numpy as np
import physical_threebody as p
R=Path(__file__).parent
s,t,w,p0,p1,p2,p3=p.dalitz(5)
u0=p.spinor(p0,p.MMU);u1=p.spinor(p1,p.ME);u2=p.spinor(p2,p.ME)
v3=p.spinor(p3,p.ME,True);b1=p.bar(u1);b2=p.bar(u2)
K=p.kernels(p0,p1,p2,p3)
errors=[]
for idx,P in enumerate(p.projectors):
    T=np.zeros((len(s),16),complex)
    for mu in range(4):
        for nu in range(4):
            G=-1j*p.isig[mu,nu]@P
            T+=p.metric[mu]*p.metric[nu]*(p.product(p.bil(b1,G,u0),p.bil(b2,G,v3))-p.swapped(p.bil(b2,G,u0),p.bil(b1,G,v3)))
    error=np.max(abs(T+12*K[:,idx]))/np.max(abs(T))
    assert error<1e-12
    errors.append(float(error))
z=json.loads((R/'controlled_domain_v11.json').read_text())
assert z['amplitude_budget_3e']['sum']<.25
assert z['amplitude_budget_dipole']['sum']<.10
assert z['enclosure']['R_lower']>10*p.old.ALPHA/(3*np.pi)*(np.log(p.MMU**2/p.ME**2)-11/4)
out={'status':'PASS','same_chirality_antisymmetrized_tensor_equals_minus12_scalar_relative_errors':errors,'neutral_box_topology_scope':'two fermion-line contractions, two heavy-line orderings; coefficient majorant is not a finite UV box calculation','conditional_ratio_enclosure':z['enclosure'],'unknown_matching_status':'Residual amplitude budgets are requirements, not computed predictions.'}
(R/'validation_domain_v11.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))

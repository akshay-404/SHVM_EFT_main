"""Optional Higgs admixture diagnostic: neutral loop and finite top/W-photon Barr--Zee.
Not a complete two-loop EW matching: Z diagrams, set C and parent-specific thresholds omitted.
"""
import numpy as np
import json
from scipy.integrate import quad
from scipy.optimize import brentq
from pathlib import Path
import physical_threebody as p
R = Path(__file__).parent
H = np.load(R/'gram_320.npy')
Q = p.fields()
qh = np.diag([p.ME, p.MMU, p.MTAU])/p.old.V
MW = 80.379
MT = 172.5
YTT = MT/p.old.V
kappa = p.old.ALPHA*np.sqrt(2)*p.old.GF*p.old.V/(4*np.pi*p.MMU)


def fg(z):
    def fun(x, kind):
        a = x*(1-x)
        t = a/z
        base = np.log(t)/(a-z) if abs(t-1) > 1e-7 else (1-(t-1)/2+(t-1)**2/3)/z
        return z/2*base*(1-2*a if kind == 'f' else 1)
    return [quad(lambda x: fun(x, kind), 0, 1, epsabs=1e-11)[0] for kind in ['f', 'g']]


def baz(M):
    ft, _ = fg(MT**2/M**2)
    f, g = fg(MW**2/M**2)
    z = MW**2/M**2
    # h=z dg/dz-g; symmetric numerical derivative with convergence checked separately.
    eps = 1e-4
    gp = fg(z*np.exp(eps))[1]
    gm = fg(z*np.exp(-eps))[1]
    h = (gp-gm)/(2*eps)-g
    return -6*(2/3)**2*p.old.V/MT*YTT*ft+3*f+5*g+.75*g+.75*h+(f-g)/(2*z)


def evaluate(i, theta):
    qq = np.concatenate([Q, qh[None]], axis=0)
    O = np.eye(9)
    O[i, i] = O[8, 8] = np.cos(theta)
    O[8, i] = np.sin(theta)
    O[i, 8] = -np.sin(theta)
    q = p.rotate(qq, O) if hasattr(
        p, 'rotate') else np.einsum('bi,ijk->bjk', O, qq)
    c = p.coefficients(q, np.r_[np.full(8, 5000), 125])
    c1 = c.copy()
    delta = kappa*np.sin(theta)*np.cos(theta)*(baz(125)-baz(5000)) * \
        np.array([Q[i, 0, 1], Q[i, 1, 0].conjugate()])
    c[8:10] -= p.E2*delta/(4*np.pi**2)
    y = q[8]
    bhlfv = 125/(8*np.pi*.00407)*(abs(y[0, 1])**2+abs(y[1, 0])**2)
    return {'B3e_neutral_one_loop': p.rate(c1, H), 'Beg_neutral_one_loop': p.radiative(c1), 'R_neutral_one_loop': p.rate(c1, H)/p.radiative(c1), 'B3e_plus_topW_photon': p.rate(c, H), 'Beg_plus_topW_photon': p.radiative(c), 'R_plus_topW_photon': p.rate(c, H)/p.radiative(c), 'B_h_e_mu_fixed_SM_width': float(bhlfv), 'kappa_h_V': float(np.cos(theta)), 'kappa_h_mu_complex': [float((y[1, 1]/qh[1, 1]).real), float((y[1, 1]/qh[1, 1]).imag)], 'note': 'Diagnostic only; no global Higgs fit or complete EW two-loop result.'}


out = {'topW_photon_b_functions': {'125': baz(125), '5000': baz(
    5000)}, 'kappa': kappa, 'neutral_scalar_external_mass_expansion_m_mu_squared_over_m_h_squared': p.MMU**2/125**2, 'cases': {}}
for i, label in [(0, 's1'), (2, 's2'), (4, 's4'), (6, 's5')]:
    out['cases'][label] = {str(t): evaluate(i, t)
                           for t in [-.01, -.001, 0, .001, .01]}
(R/'higgs_mixing_v11.json').write_text(json.dumps(out, indent=2)+'\n')
if __name__ == '__main__':
    print(json.dumps(out, indent=2))

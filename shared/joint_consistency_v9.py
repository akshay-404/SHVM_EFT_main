"""Reproducible v9 cross-manuscript checks. Units: GeV except top matrix (TeV).
Run from any working directory: python3 shared/joint_consistency_v9.py
Requires numpy, scipy. This verifies the supplied witnesses and interface, not a global fit.
"""
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
ROOT = Path(__file__).resolve().parents[1]
# The shared archive includes the same LFV implementation in shared/lfv_letter.py.
spec = importlib.util.spec_from_file_location(
    'lfv', Path(__file__).with_name('lfv_letter.py'))
L = importlib.util.module_from_spec(spec)
spec.loader.exec_module(L)
P = dict(MT=6.398573980858272, MD=6.398573980858272, l=6.398573980858272,
         r=3.4857798273934257, yT=2.3736052677477124, yD=2.3736052677477124, FT=1., FD=1.)
f, v, chi3 = 2., .246, 4.35
th = math.asin(v/f)


def matrix(p, angle=th):
    m = np.zeros((6, 6))
    m[:2, :2] = p['l']*np.eye(2)
    m[:2, 2:4] = p['r']*np.eye(2)
    m[2:4, :2] = p['MT']*np.eye(2)
    m[4:6, 2:4] = p['MD']*np.eye(2)
    m[2:4, 4] = p['yT']*p['FT']*np.array([math.sin(angle), math.cos(angle)])
    m[4:6, 5] = p['yD']*p['FD']*np.array([-math.sin(angle), math.cos(angle)])
    return m


U, m, Vh = np.linalg.svd(matrix(P))
UR = Vh.conj().T
j = np.argmin(m)
Dm = np.zeros((6, 6))
Dm[:2, 2:4] = P['r']/(math.sqrt(2)*chi3)*np.eye(2)
gs = (U.conj().T@Dm@UR)[j, j]
ga = (U.conj().T@(1j*Dm)@UR)[j, j].imag
h = 1e-5
def mt(p): return np.linalg.svd(matrix(p), compute_uv=False).min()


eta = (math.log(mt(dict(P, r=P['r']*math.exp(h)))) -
       math.log(mt(dict(P, r=P['r']*math.exp(-h)))))/(2*h)
assert abs(gs-ga) < 1e-13 and abs(gs-eta*m[j]/(math.sqrt(2)*chi3)) < 1e-10
ye = L.yukawa(L.YE)
me = L.mass_matrix(ye)
ue, se, vhe = np.linalg.svd(me)
ue = ue[:, np.argsort(se)]
YN = {(1, 1): (.821553, 5.464722), (1, 2): (.994617, 5.971534), (1, 3): (1.028904, 3.854641), (2, 2): (
    1.143514, 5.627429), (2, 3): (1.187142, 1.396783), (3, 2): (1.031627, 1.386256), (3, 3): (1.000453, 4.705352)}
NT = {(1, 1): 1, (1, 2): 4, (1, 3): 4, (2, 2)      : 4, (2, 3): 4, (3, 2): 5, (3, 3): 5}
mn = np.zeros((3, 3), complex)
for (i, k), (a, ph) in YN.items():
    mn[i-1, k-1] = L.V/math.sqrt(2)*7.18e-10*a * \
        np.exp(1j*ph)*L.EPS[NT[(i, k)]]*1e9
un, sn, vhn = np.linalg.svd(mn)
order = np.argsort(sn)
sn = sn[order]
un = un[:, order]
pmns = ue.conj().T@un
s13 = abs(pmns[0, 2])**2
s12 = abs(pmns[0, 1])**2/(1-s13)
s23 = abs(pmns[1, 2])**2/(1-s13)
J = (pmns[0, 0]*pmns[1, 1]*np.conj(pmns[0, 1]*pmns[1, 0])).imag
sd = J/math.sqrt(s12*(1-s12)*s23*(1-s23)*s13)/(1-s13)
cd = (abs(pmns[1, 0])**2-s12*(1-s23)-(1-s12)*s23*s13) / \
    (2*math.sqrt(s12*(1-s12)*s23*(1-s23)*s13))
delta = math.degrees(math.atan2(sd, cd)) % 360
neutrino = dict(masses_eV=sn.tolist(), s12sq=s12, s23sq=s23, s13sq=s13,
                delta_deg=delta, J=J, dm21sq=sn[1]**2-sn[0]**2, dm31sq=sn[2]**2-sn[0]**2)
# Integral representations test the finite dipole functions over the physical hierarchy.
loop_errors = []
for x in (1e-14, 1e-10, 1e-7, .01, .3):
    fi = quad(lambda z: .5*z*(1-z)**2/(z+x*(1-z)), 0, 1, epsabs=1e-11)[0]
    hi = quad(lambda z: (1-z)**2/(z+x*(1-z)), 0, 1, epsabs=1e-10)[0]
    loop_errors.append([x, abs(fi-L.Ffun(x)), abs(hi-L.Hfun(x))])
assert max(a[1] for a in loop_errors) < 1e-7 and max(a[2]
                                                     for a in loop_errors) < 1e-7
_, cp, cl, cr, b, C, bc = L.rates(ye, 5000.)
sub = L.contact_dipole(cl, cr, C, bc)
assert sub['B_subset'] >= sub['independent_phase_envelope'][0]-1e-25
assert sub['B_subset'] <= sub['independent_phase_envelope'][1]+1e-25
# Mass-basis rephasing must preserve all rates, including fixed interference.
D = np.diag(np.exp(1j*np.array([.2, 1.1, -.7])))
cp2 = {r: (D.conj().T@S@D, D.conj().T@A@D) for r, (S, A) in cp.items()}
cl2, cr2, b2 = L.dipole(cp2, 5000.)
C2, bc2 = L.contact(cp2, 5000.)
sub2 = L.contact_dipole(cl2, cr2, C2, bc2)
rephase = max(abs(b2/b-1), abs(bc2/bc-1),
              abs(sub2['B_subset']/sub['B_subset']-1))
assert rephase < 1e-12
# Orientation identity independently checked before making cancellation claims.
orientation_errors = []
for S, A in cp.values():
    gp = (S+A)/2
    gm = (S-A)/2
    orientation_errors.append(np.max(abs(S@S-A@A-2*(gp@gm+gm@gp))))
assert max(orientation_errors) < 1e-15
# Explicit formal one-loop field-content extrapolations. No UV exclusion beyond thresholds.
aY = (1/127.951)/(1-.23122)


def landau(N, ND, Lam):
    strong = 2*N/3+2*N/3+20*ND/9+2*N*ND/27
    thresholds = [(1749.5, 16/9), (6398.6, 2.), (6398.6, 2.), (Lam, strong)]
    inv = 1/aY
    mu = 91.1876
    bY = 41/6
    for t, db in sorted(thresholds):
        inv -= bY/(2*math.pi)*math.log(t/mu)
        mu = t
        bY += db
    return mu*math.exp(inv*2*math.pi/bY), strong


links = []
for N, ND in [(2, 3), (3, 4)]:
    for ld in (2000., 3000., 5000.):
        ly, bs = landau(N, ND, ld)
        X = ND/(4*math.pi**2)/7.18e-10
        links.append(dict(N=N, ND=ND, LambdaD_GeV=ld, bY_strong=bs, LambdaY_formal_GeV=ly,
                          link_dim3_GeV=ld*math.sqrt(X), link_dim2_GeV=ld*X,
                          gamma_for_2e16=2-math.log(X)/math.log(2e16/ld),
                          gamma_max_formal=2-math.log(X)/math.log(ly/ld)))
assert abs(links[2]['bY_strong']-88/9) < 1e-12
out = dict(status='PASS', top=dict(params_TeV=P, mt_GeV=m[j]*1000, eta3=eta, g_s3=gs, g_a3=ga, kappa_chi=P['r']/chi3),
           charged_masses_MeV=(np.sort(se)*1e3).tolist(), neutrino=neutrino,
           loop_function_integral_errors=loop_errors, rephasing_relative_error=rephase,
           orientation_identity_max_error=float(max(orientation_errors)), contact_dipole=sub,
           formal_links=links,
           limitations=['No global-minimum certificate or full flavour likelihood', 'No physical q-dependent charge-radius, Z/box or LFV QED running', 'No finite two-loop dynamical-singlet threshold or complete UV selector assignment', 'Formal hypercharge extrapolation assumes no additional thresholds'])
Path(__file__).with_name('joint_consistency_v9.json').write_text(
    json.dumps(out, indent=2)+'\n')
print(json.dumps({k: v for k, v in out.items() if k not in [
      'top', 'limitations', 'formal_links']}, indent=2))

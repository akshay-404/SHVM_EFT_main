"""Independent algebra checks for the sufficient isolated dynamical bridge.
This is a trace, contraction, charge and tree-boundary audit, not computation
of unlisted finite loops or automatic proof of the full UV flavour model.
"""
import json
from pathlib import Path
import sympy as s
ell, k, MT, MD, yt, yd, cr, ci = s.symbols(
    'ell k MT MD yt yd cr ci', real=True)
tr = s.symbols('tr0:3', real=True)
ti = s.symbols('ti0:3', real=True)
dr = s.symbols('dr0:3', real=True)
di = s.symbols('di0:3', real=True)
pt = s.Matrix([a+s.I*b for a, b in zip(tr, ti)])
pd = s.Matrix([a+s.I*b for a, b in zip(dr, di)])
chi = cr+s.I*ci
r2 = k*k*(cr*cr+ci*ci)
M = s.zeros(9, 8)
M[:3, :3] = ell*s.eye(3)
M[:3, 3:6] = k*chi*s.eye(3)
M[3:6, :3] = MT*s.eye(3)
M[6:9, 3:6] = MD*s.eye(3)
M[3:6, 6] = yt*pt
M[6:9, 7] = yd*pd
rt = (pt.conjugate().T*pt)[0]
rd = (pd.conjugate().T*pd)[0]
expected = 3*((MT**2+ell**2)**2+(MD**2+r2)**2+2*ell**2*r2)+2 * \
    MT**2*yt**2*rt+2*MD**2*yd**2*rd+yt**4*rt**2+yd**4*rd**2
H = M.conjugate().T*M
assert s.expand(s.trace(H*H)-expected) == 0
# Complete two-real-component singlet contraction.
Ys = k/s.sqrt(2)*s.eye(3)
Ya = s.I*Ys
assert s.simplify(Ys*Ys.conjugate().T+Ya*Ya.conjugate().T -
                  k*k*s.eye(3)) == s.zeros(3)
X = dict(chi=1, Q=-1, FT=-1, t1=-1, FD=-2, t2=-2, PhiT=0, PhiD=0)
invariants = {'left_bridge': -X['Q']+X['FT'], 'right_bridge': -X['Q']+X['chi']+X['FD'],
              'yT': -X['FT']+X['PhiT']+X['t1'], 'yD': -X['FD']+X['PhiD']+X['t2']}
assert all(q == 0 for q in invariants.values())
forbidden = {'conjugate_bridge': -X['Q']-X['chi']+X['FD'],
             'bare_cross_messenger_mass': -X['FT']+X['FD']}
assert all(q != 0 for q in forbidden.values())
# Left-handed Weyl sums; vectorlike F triplets cancel pairwise.
fields = [(-1, s.Rational(1, 6), 6, s.Rational(1), s.Rational(3, 2)),
          (-1, s.Rational(2, 3), 3, s.Rational(1, 2), 0),
          (1, s.Rational(-2, 3), 3, s.Rational(1, 2), 0),
          (2, s.Rational(-2, 3), 3, s.Rational(1, 2), 0)]
a = dict(c2X=0, W2X=0, Y2X=0, grav2X=0, X3=0, X2Y=0)
for q, Y, n, Tc, Tw in fields:
    a['c2X'] += q*Tc
    a['W2X'] += q*Tw
    a['Y2X'] += n*q*Y**2
    a['grav2X'] += n*q
    a['X3'] += n*q**3
    a['X2Y'] += n*q**2*Y
assert a == dict(c2X=0, W2X=s.Rational(-3, 2),
                 Y2X=s.Rational(5, 2), grav2X=0, X3=18, X2Y=-7)
# Explicit tree singlet boundary, with chi=(sqrt(2)*w+s+i*a)/sqrt(2).
w = s.Rational(435, 100)
lam = s.Rational(25)/(4*w*w)
B = s.Rational(-25, 4)
ss, aa = s.symbols('s a', real=True)
chi = (s.sqrt(2)*w+ss+s.I*aa)/s.sqrt(2)
V = s.expand(lam*(chi*s.conjugate(chi))**2+B*(chi**2+s.conjugate(chi)**2))
grad = [s.simplify(s.diff(V, x).subs({ss: 0, aa: 0})) for x in [ss, aa]]
hess = s.hessian(V, [ss, aa]).subs({ss: 0, aa: 0})
assert grad == [0, 0] and s.simplify(hess) == 25*s.eye(2)
result = dict(status='PASS', fermion_trace_identity=str(expected), bridge_contraction='|kappa|^2 identity',
              bridge_phase_charges=X, bridge_phase_invariant_vertices=invariants, forbidden_without_soft_phase_matching=forbidden,
              bridge_phase_anomalies={k: str(v) for k, v in a.items()}, singlet_tree_boundary=dict(w_TeV=float(w), lambda_chi=float(lam), B_chi_TeV2=float(B), Hessian_TeV2=[[25, 0], [0, 25]]),
              hypotheses=['isolated ungauged parent plus QCD', 'soft universal left bridge', 'universal hard right bridge with single orientation',
                          'complete radial fields and norm potential', 'hard bridge-phase symmetry preserved; B_chi soft', 'no extra hard T connector'],
              not_established=['finite dynamical-singlet loops', 'gauged or complete flavour closure', 'UV selector and constituent bridge-phase assignment', 'full RG or global scalar stability'])
Path(__file__).with_name('dynamical_bridge_audit_v9.json').write_text(
    json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))

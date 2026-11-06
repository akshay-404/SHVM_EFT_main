"""Paper I (v7 radiatively closed top sector) -> Paper II interface. Self-contained; units TeV.
Run: python3 interface/paper2_interface_v7.py  -> interface/interface_v7.json
Benchmark parameters are those of Paper I's common v7 boundary (f = 2 TeV, on-shell M_t = 172.5 GeV)."""
import json, math
from pathlib import Path
import numpy as np
P = dict(MT=6.398573980858272, MD=6.398573980858272, l=6.398573980858272, r=3.4857798273934257,
         yT=2.3736052677477124, yD=2.3736052677477124, FT=1.0, FD=1.0)
F, V, MT_OS, CHI3 = 2.0, 0.246, 0.1725, 4.35            # <chi_3> = eps_3 Lambda_D = 0.87 x 5 TeV
def up_matrix(th, p):
    m = np.zeros((6, 6)); s, c = math.sin(th), math.cos(th)
    m[0:2, 0:2] = p['l']*np.eye(2); m[0:2, 2:4] = p['r']*np.eye(2)
    m[2:4, 0:2] = p['MT']*np.eye(2); m[4:6, 2:4] = p['MD']*np.eye(2)
    m[2:4, 4] = p['yT']*p['FT']*np.array([s, c]); m[4:6, 5] = p['yD']*p['FD']*np.array([-s, c])
    return m
mt = lambda p, th: np.linalg.svd(up_matrix(th, p), compute_uv=False).min()
xi = (V/F)**2; th = math.asin(math.sqrt(xi)); h = 1e-5
eta3 = (math.log(mt(dict(P, r=P['r']*math.exp(h)), th))-math.log(mt(dict(P, r=P['r']*math.exp(-h)), th)))/(2*h)
J = 4*xi*(1-xi)/F**2; Vpp = 0.1250637**2/J; alpha_unit = F*F*CHI3*CHI3/2
out = dict(light_eigenvalue=mt(P, th), eta3=eta3, g_s3=eta3*MT_OS/(math.sqrt(2)*CHI3), g_a3=eta3*MT_OS/(math.sqrt(2)*CHI3),
           lambda_H3_order_one_mh_shift=0.125**2/(J*2*alpha_unit/(1-2*xi)), lambda_H3_vacuum_bound=xi*Vpp/alpha_unit)
# one-loop hypercharge running, Q = T3 + Y; d alpha^-1/d ln mu = -b/(2 pi)
aY = (1/127.951)/(1-0.23122); MZ = 0.0911876
def landau(thr, b=41/6):
    inv, mu = 1/aY, MZ
    for s, db in sorted(thr):
        inv -= b/(2*math.pi)*math.log(s/mu); mu = s; b += db
    return mu*math.exp(inv*2*math.pi/b)*1e3
out['LambdaY_parent_GeV'] = landau([(1.7495, 16/9), (6.3986, 2.0), (6.3986, 2.0), (5.0, 88/9)])
out['LambdaY_parent_messengers_only_GeV'] = landau([(1.7495, 16/9), (6.3986, 2.0), (6.3986, 2.0)])
out['LambdaY_sextet_chain_validation_GeV'] = landau([(2.756, 44/3), (7.722, 44/3), (11.159, 44/3), (5.0, 88/9)])
X = (5.14e7/5e3)**2      # (Lambda_GUT^(3)/Lambda_D)^2 at Lambda_D = 5 TeV
out['Delta_nu_min'] = 1+math.log(X)/math.log(out['LambdaY_parent_GeV']/5e3)
out['colour_slope_above_messengers'] = -7+2+2+2/3
Path(__file__).with_name('interface_v7.json').write_text(json.dumps(out, indent=1)+'\n')
print(json.dumps(out, indent=1))

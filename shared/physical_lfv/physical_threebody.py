"""Leading 1/M^2 neutral-real-scalar QED amplitude, finite q^2, exact 3e phase space.
No UV electroweak matching or inclusive NLO QED rate is implied.
Run this file to build/check integrated operator Gram matrices.
"""
import numpy as np
import math
import json
from pathlib import Path
from numpy.polynomial.legendre import leggauss
import lfv_letter as old
ME, MMU, MTAU = old.ME, old.MMU, old.MTAU
E2 = 4*np.pi*old.ALPHA
I2 = np.eye(2, dtype=complex)
Z2 = np.zeros((2, 2), complex)
pauli = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                 [[1, 0], [0, -1]]], complex)
gamma = np.array([np.block([[I2, Z2], [Z2, -I2]])] +
                 [np.block([[Z2, s], [-s, Z2]]) for s in pauli])
g5 = 1j*gamma[0]@gamma[1]@gamma[2]@gamma[3]
PL = (np.eye(4)-g5)/2
PR = (np.eye(4)+g5)/2
projectors = [PR, PL]
metric = np.array([1, -1, -1, -1])
isig = -.5*(np.einsum('mij,njk->mnik', gamma, gamma) -
            np.einsum('nij,mjk->mnik', gamma, gamma))
CH = [(0, 0), (1, 1), (0, 1), (1, 0)]  # RR LL RL LR
NAMES = ['C_RR', 'C_LL', 'C_RL', 'C_LR']+['V_RR', 'V_LL', 'V_RL', 'V_LR']+['D_R', 'D_L'] + \
    [f'J_{x}_{k}' for x in ('R', 'L') for k in (
        'e', 'mu', 'tau')]+[f'I_{x}_{k}' for x in ('R', 'L') for k in ('e', 'mu', 'tau')]


def logmoments(r):
    """I=integral log(1-r z(1-z)-i0); J=integral z(1-z) log(...)."""
    r = np.asarray(r, float)
    i = np.empty(r.shape, complex)
    j = np.empty(r.shape, complex)
    sm = r < .01
    i[sm] = 0
    j[sm] = 0
    # convergent exact expansion, beta(n+1,n+1), beta(n+2,n+2)
    for n in range(1, 9):
        i[sm] -= r[sm]**n/math.comb(2*n, n)/(n*(2*n+1))
        j[sm] -= r[sm]**n*math.factorial(n+1)**2/math.factorial(2*n+3)/n
    low = (r >= .01) & (r <= 4)
    a = np.sqrt(np.maximum(4/r[low]-1, 0))
    i[low] = -2+2*a*np.arctan2(1, a)
    high = r > 4
    b = np.sqrt(1-4/r[high])
    i[high] = -2+b*(np.log1p(b)-np.log1p(-b))-1j*np.pi*b
    j[~sm] = 1/18+(1/6+1/(3*r[~sm]))*i[~sm]
    return i, j


def spinor(p, m, anti=False):
    # Columns are the two spin states. Standard Dirac normalization.
    ep = np.sqrt(p[:, 0]+m)
    sp = np.einsum('ni,ijk->njk', p[:, 1:], pauli)/ep[:, None, None]
    hi = ep[:, None, None]*np.broadcast_to(I2, (len(p), 2, 2))
    return np.concatenate([sp, hi] if anti else [hi, sp], axis=1)


def bar(u): return np.einsum('nai,ab->nib', u.conj(), gamma[0])
def bil(b, G, u): return np.einsum('nai,ij,njb->nab', b, G, u, optimize=True)
def bilv(b, G, u): return np.einsum('nai,nij,njb->nab', b, G, u, optimize=True)
def product(a, b): return np.einsum('nad,nbc->nabcd', a, b).reshape(len(a), 16)
def swapped(a, b): return np.einsum('nbd,nac->nabcd', a, b).reshape(len(a), 16)


def dalitz(n, me=ME):
    z, w = leggauss(n)
    v = (z+1)/2
    w = w/2
    m = MMU
    smin = 4*me**2
    smax = (m-me)**2
    s = smin*np.exp(v*np.log(smax/smin))
    ws = w*s*np.log(smax/smin)
    lam = np.sqrt(np.maximum((m*m-me*me-s)**2-4*me*me*s, 0))
    mid = 2*me*me+(m*m-me*me-s)/2
    half = lam*np.sqrt(1-4*me*me/s)/2
    tmin = mid-half
    tmax = mid+half
    t = tmin[:, None]*np.exp(v[None, :]*np.log(tmax/tmin)[:, None])
    wt = w[None, :]*t*np.log(tmax/tmin)[:, None]
    s = np.repeat(s, n)
    wts = (ws[:, None]*wt).ravel()
    t = t.ravel()
    u = m*m+3*me*me-s-t
    e1 = (m*m+me*me-s)/(2*m)
    e2 = (m*m+me*me-t)/(2*m)
    pz = np.sqrt(np.maximum(e1*e1-me*me, 0))
    p2 = np.sqrt(np.maximum(e2*e2-me*me, 0))
    c = (e1*e2-(u-2*me*me)/2)/(pz*p2)
    c = np.clip(c, -1, 1)
    p1 = np.stack([e1, np.zeros(len(s)), np.zeros(len(s)), pz], axis=1)
    p2 = np.stack([e2, p2*np.sqrt(1-c*c), np.zeros(len(s)), p2*c], axis=1)
    p3 = np.column_stack([m-e1-e2, -p1[:, 1:]-p2[:, 1:]])
    p0 = np.tile([m, 0, 0, 0], (len(s), 1))
    return s, t, wts, p0, p1, p2, p3


def kernels(p0, p1, p2, p3):
    n = len(p0)
    u0 = spinor(p0, MMU)
    u1 = spinor(p1, ME)
    u2 = spinor(p2, ME)
    v3 = spinor(p3, ME, anti=True)
    b1 = bar(u1)
    b2 = bar(u2)
    q1 = p2+p3
    q2 = p1+p3
    s = np.sum(q1*q1*metric, axis=1)
    t = np.sum(q2*q2*metric, axis=1)
    C = []
    V = []
    D = []
    vfirst = []
    vsecond = []
    dfirst = []
    dsecond = []
    for x, y in CH:
        C.append(product(bil(b1, projectors[x], u0), bil(
            b2, projectors[y], v3))-swapped(bil(b2, projectors[x], u0), bil(b1, projectors[y], v3)))
        z = np.zeros((n, 16), complex)
        for mu in range(4):
            z += metric[mu]*(product(bil(b1, gamma[mu]@projectors[x], u0), bil(b2, gamma[mu]@projectors[y], v3)) -
                             swapped(bil(b2, gamma[mu]@projectors[x], u0), bil(b1, gamma[mu]@projectors[y], v3)))
        V.append(z)
    for P in projectors:
        a = np.zeros((n, 16), complex)
        b = a.copy()
        d = a.copy()
        f = a.copy()
        for mu in range(4):
            a += metric[mu]*product(bil(b1, gamma[mu]@P, u0),
                                    bil(b2, gamma[mu], v3))
            b -= metric[mu]*swapped(bil(b2, gamma[mu]@P, u0),
                                    bil(b1, gamma[mu], v3))
            g1 = MMU*np.einsum('vij,nv->nij',
                               isig[mu]@P, q1*metric)/s[:, None, None]
            g2 = MMU*np.einsum('vij,nv->nij',
                               isig[mu]@P, q2*metric)/t[:, None, None]
            d += metric[mu]*product(bilv(b1, g1, u0), bil(b2, gamma[mu], v3))
            f -= metric[mu]*swapped(bilv(b2, g2, u0), bil(b1, gamma[mu], v3))
        vfirst.append(a)
        vsecond.append(b)
        dfirst.append(d)
        dsecond.append(f)
        D.append(d+f)
    Js = []
    Is = []
    for x in range(2):
        for mk in [ME, MMU, MTAU]:
            i, j = logmoments(s/mk**2)
            ii, jj = logmoments(t/mk**2)
            Js.append(vfirst[x]*j[:, None]+vsecond[x]*jj[:, None])
            Is.append(dfirst[x]*i[:, None]+dsecond[x]*ii[:, None])
    return np.stack(C+V+D+Js+Is, axis=1)


def gram(n=100):
    s, t, w, p0, p1, p2, p3 = dalitz(n)
    H = np.zeros((22, 22), complex)
    for start in range(0, len(s), 1000):
        ix = slice(start, start+1000)
        a = kernels(p0[ix], p1[ix], p2[ix], p3[ix])
        weights = w[ix]/(1024*np.pi**3*MMU**3*old.GAMMA_MU)
        H += np.einsum('nki,nli,n->kl', a.conj(), a, weights, optimize=True)
    return (H+H.conj().T)/2


def fields(Lam=5000):
    _, cp = old.couplings(old.yukawa(old.YE), Lam)
    return np.array([Q for S, A in cp.values() for Q in (S, 1j*A)])


def coefficients(Q, masses):
    """Tree and finite photon Wilson terms in basis NAMES; Q_b general real-field coupling."""
    masses = np.asarray(masses)
    invm = 1/masses**2
    c = np.zeros(22, complex)
    a2 = np.zeros(2, complex)
    a1 = np.zeros(2, complex)
    mk = np.array([ME, MMU, MTAU])
    ll = np.empty((2, 3), complex)
    cross = ll.copy()
    for b, q in enumerate(Q):
        wr = invm[b]
        c[:4] += wr*np.array([q[0, 1]*q[0, 0], q[1, 0].conjugate()*q[0, 0].conjugate(),
                             q[0, 1]*q[0, 0].conjugate(), q[1, 0].conjugate()*q[0, 0]])
        v = np.array([q[0, :]*q[1, :].conj(), q[:, 0].conj()*q[:, 1]])
        h = np.array([q[0, :]*q[:, 1], q[:, 0].conj()*q[1, :].conj()])
        # vector current PR uses right-handed external leptons, i.e. v[1]; dipole PR uses h[0].
        for x in range(2):
            a1[x] += wr * \
                np.sum(v[1-x]*(np.log(mk**2/masses[b]**2)+4/3))/(48*np.pi**2)
            a2[x] += wr*np.sum((v[x]+ME/MMU*v[1-x])/6+mk/MMU*h[x]
                               * (np.log(masses[b]**2/mk**2)-1.5))/(16*np.pi**2)
            c[10+x*3:13+x*3] += E2*wr*v[1-x]/(8*np.pi**2)
            c[16+x*3:19+x*3] += -E2*wr*mk/MMU*h[x]/(16*np.pi**2)
    c[4:8] = E2*np.array([a1[0], a1[1], a1[0], a1[1]])
    c[8:10] = E2*a2
    # HKZ Eq. A13--A16: physical vector = -old A1; dipole = -m i sigma q d.
    # Scalar C amplitude has the positive tree sign.
    c[4:] *= -1
    return c


def rate(c, H): return float(np.real(c.conj()@H@c))


def radiative(c):
    a2 = c[8:10]/E2
    return old.ALPHA*MMU**5/4/old.GAMMA_MU*np.sum(abs(a2)**2)


def validate(H):
    _, cp = old.couplings(old.yukawa(old.YE), 5000)
    CL, CR, Bg = old.dipole(cp, 5000)
    C, Bc = old.contact(cp, 5000)
    c = coefficients(fields(), np.full(8, 5000.))
    cc = c.copy()
    cc[4:] = 0
    dd = c.copy()
    dd[:8] = 0
    dd[10:] = 0
    cd = cc+dd
    pref = old.MU_NORM/old.GF**2
    diagtargets = np.array([1/64, 1/64, 1/32, 1/32, 1/4, 1/4, 1/8, 1/8])*pref
    result = {'contact': rate(cc, H), 'contact_analytic': Bc, 'dipole': rate(dd, H), 'dipole_analytic': (old.ALPHA/(3*np.pi)*(np.log(MMU**2/ME**2)-11/4))*Bg, 'cd': rate(cd, H), 'radiative': radiative(
        c), 'radiative_v9': Bg, 'physical_photon_tree': rate(c, H), 'gram_contact_vector_diag_relative': (np.diag(H).real[:8]/diagtargets-1).tolist(), 'min_gram_eigenvalue': float(np.linalg.eigvalsh(H)[0])}
    return result


if __name__ == '__main__':
    import sys
    import time
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    t = time.time()
    H = gram(n)
    np.save(Path(__file__).with_name(f'gram_{n}.npy'), H)
    out = validate(H)
    out.update(n=n, seconds=time.time()-t)
    print(json.dumps(out, indent=2))
    Path(__file__).with_name(f'validation_{n}.json').write_text(
        json.dumps(out, indent=2))

"""Legacy v9 contact and on-shell diagnostics; v10 uses physical_threebody.py.

Self-contained; units GeV. Inputs are those of Paper II (Table "Lepton fit", Eq. epsilon-benchmark,
texture Eq. charged-matrices, orientation of chi_2, couplings Eq. fiveTeV-couplings,
complete dipole Eq. fiveTeV-full-dipole, contact Eqs. mu3e-contact/mu3e-rate).

Run: python3 lfv_letter.py  -> lfv_letter.json
"""
import json, math, cmath
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, brentq

# ---------------------------------------------------------------- constants
ALPHA = 1/137.035999       # on-shell photon
GF = 1.1663788e-5
HBAR = 6.582119569e-25     # GeV s
TAU_MU = 2.1969811e-6      # s
ME, MMU, MTAU = 0.51099895e-3, 0.1056583755, 1.77686
GAMMA_MU = HBAR/TAU_MU
MU_NORM = GF**2*MMU**5/(192*math.pi**3*GAMMA_MU)  # common measured lifetime for both channels
V = 246.0
EPS = {1: 3.159e-6, 2: 3.10000e-3, 4: 6.1001e-5, 5: 2.70002e-4}     # from the benchmark lepton VEVs
# Lepton fit (printed-input witness; no certificate of global minimality), magnitudes and phases
YE = {(1, 1): (0.907574, 0.915208), (1, 2): (1.008157, 1.857883), (1, 3): (1.001456, 6.112102),
      (2, 2): (2.781819, 5.489036), (2, 3): (1.988401, 1.814003), (3, 3): (2.626645, 1.421511)}
# texture: which epsilon multiplies each entry, and orientation (+1 chi, -1 chi^dagger)
TEXTURE = {(1, 1): (1, +1), (1, 2): (4, +1), (1, 3): (5, +1), (2, 2): (5, +1), (2, 3): (2, -1), (3, 3): (2, +1)}
def yukawa(ye):
    return {k: m*cmath.exp(1j*ph) for k, (m, ph) in ye.items()}


def mass_matrix(y):
    M = np.zeros((3, 3), complex)
    for (i, j), val in y.items():
        r, _ = TEXTURE[(i, j)]
        M[i-1, j-1] = V/math.sqrt(2)*val*EPS[r]
    return M


def couplings(y, LamD):
    """S_r, A_r in the mass basis (Eq. fiveTeV-couplings); returns masses and dict r->(S,A)."""
    M = mass_matrix(y)
    U, s, Vh = np.linalg.svd(M)            # M = U diag(s) Vh
    order = np.argsort(s); U = U[:, order]; s = s[order]; UR = Vh.conj().T[:, order]
    out = {}
    for r in EPS:
        Y = np.zeros((3, 3), complex); Ya = np.zeros((3, 3), complex)
        for (i, j), val in y.items():
            rr, o = TEXTURE[(i, j)]
            if rr == r:
                Y[i-1, j-1] = val; Ya[i-1, j-1] = o*val
        pref = V/(2*LamD)
        out[r] = (pref*U.conj().T@Y@UR, pref*U.conj().T@Ya@UR)
    return s, out


def Ffun(x):
    if x < 1e-6:
        pass  # evaluate the stable exact expression below
    return (x*x-5*x-2)/(12*(x-1)**3)+x*math.log(x)/(2*(x-1)**4)


def Hfun(x):
    if x < 1e-6:
        pass  # evaluate the stable exact expression below
    return (x-3)/(2*(x-1)**2)+math.log(x)/(x-1)**3


def dipole(cpl, Mscal, channels=None):
    """Complete neutral-scalar one-loop dipole (Eq. fiveTeV-full-dipole). Returns c_L, c_R, B."""
    mk = [ME, MMU, MTAU]; e, mu = 0, 1
    cL = cR = 0j
    for r, (S, A) in cpl.items():
        if channels and r not in channels:
            continue
        for Q in (S, 1j*A):
            for k in range(3):
                x = mk[k]**2/Mscal**2
                pre = 1/(4*Mscal**2)
                cR += pre*(Ffun(x)*(Q[e, k]*np.conj(Q[mu, k])+ME/MMU*np.conj(Q[k, e])*Q[k, mu])
                           + mk[k]/MMU*Hfun(x)*Q[e, k]*Q[k, mu])
                cL += pre*(Ffun(x)*(np.conj(Q[k, e])*Q[k, mu]+ME/MMU*Q[e, k]*np.conj(Q[mu, k]))
                           + mk[k]/MMU*Hfun(x)*np.conj(Q[k, e])*np.conj(Q[mu, k]))
    G = ALPHA*MMU**5/(64*math.pi**4)*(abs(cL)**2+abs(cR)**2)
    return cL, cR, G/GAMMA_MU


def contact(cpl, Mscal, channels=None):
    """Tree contact coefficients (Eq. mu3e-contact) and the contact rate (Eq. mu3e-rate)."""
    e, mu = 0, 1
    CRR = CLL = CRL = CLR = 0j
    for r, (S, A) in cpl.items():
        if channels and r not in channels:
            continue
        m2 = Mscal**2
        CRR += (S[e, mu]*S[e, e]-A[e, mu]*A[e, e])/m2
        CLL += (np.conj(S[mu, e])*np.conj(S[e, e])-np.conj(A[mu, e])*np.conj(A[e, e]))/m2
        CRL += (S[e, mu]*np.conj(S[e, e])+A[e, mu]*np.conj(A[e, e]))/m2
        CLR += (np.conj(S[mu, e])*S[e, e]+np.conj(A[mu, e])*A[e, e])/m2
    B = (abs(CRR)**2+abs(CLL)**2)/(64*GF**2)+(abs(CRL)**2+abs(CLR)**2)/(32*GF**2)
    return dict(CRR=CRR, CLL=CLL, CRL=CRL, CLR=CLR), MU_NORM*B


# Kitano-Koike-Okada dipole overlap integrals (units m_mu^{5/2}) and capture rates (1/s)
KKO = {'Al': (0.0362, 0.7054e6), 'Ti': (0.0864, 2.59e6), 'Au': (0.189, 13.07e6)}


def conversion_over_meg(nucleus):
    """Dipole-dominated B(mu N -> e N)/B(mu -> e gamma) = 2 G_F^2 D^2/(Gamma_capt 384 pi^2)."""
    D, w = KKO[nucleus]
    Dg = D*MMU**2.5
    return 2*GF**2*Dg**2/(HBAR*w*384*math.pi**2*MU_NORM)


def rates(y, M):
    s, cpl = couplings(y, M)
    cL, cR, Beg = dipole(cpl, M)
    C, B3 = contact(cpl, M)
    return s, cpl, cL, cR, Beg, C, B3


# ---------------------------------------------------------------- robustness sampling
def sample_mass_preserving(rng, n, M, targets):
    """Random order-one charged-lepton coefficients that reproduce the three charged-lepton masses:
    off-diagonal magnitudes log-uniform in [1/3,3], all phases uniform; diagonal magnitudes solved;
    accepted if the solved diagonals lie in [1/3,3]. Neutrino constraints are NOT imposed, so this
    class is broader than the actual fit (a declared-prior diagnostic, not a probability over viable theories)."""
    Rs, Begs, B3s = [], [], []
    tries = 0
    while len(Rs) < n and tries < 40*n:
        tries += 1
        ph = rng.uniform(0, 2*math.pi, 6)
        off = {k: 10**rng.uniform(-math.log10(3), math.log10(3)) for k in [(1, 2), (1, 3), (2, 3)]}
        def build(d):
            y = {(1, 1): (d[0], ph[0]), (2, 2): (d[1], ph[1]), (3, 3): (d[2], ph[2]),
                 (1, 2): (off[(1, 2)], ph[3]), (1, 3): (off[(1, 3)], ph[4]), (2, 3): (off[(2, 3)], ph[5])}
            return yukawa(y)
        def resid(ld):
            sv = np.sort(np.linalg.svd(mass_matrix(build(np.exp(ld))), compute_uv=False))
            return np.log(sv/targets)
        sol = least_squares(resid, np.log([0.9, 2.5, 2.5]), xtol=1e-13, ftol=1e-13)
        d = np.exp(sol.x)
        if max(abs(resid(sol.x))) > 1e-8 or d.min() < 1/3 or d.max() > 3:
            continue
        _, cp, cL_, cR_, Beg, _, B3 = rates(build(d), M)
        Rs.append(B3/Beg); Begs.append(Beg); B3s.append(B3)
    return np.array(Rs), np.array(Begs), np.array(B3s), tries



def contact_dipole(cL, cR, C, Bcontact):
    """Kuno-Okada Eq.117, with the HKZ dipole L = +e m_mu c/(8 pi^2).
    All rates divided by the same measured muon width. AR*, AL* map to cR,cL;
    g6*,g5* map to CLR,CRL after the -1/2 Fierz identity.
    This is a contact+dipole subset, not the full physical three-body rate.
    """
    e = math.sqrt(4*math.pi*ALPHA)
    AR = -e*np.conj(cR)/(16*math.sqrt(2)*math.pi**2*GF)
    AL = -e*np.conj(cL)/(16*math.sqrt(2)*math.pi**2*GF)
    g6 = np.conj(C['CLR'])/(4*math.sqrt(2)*GF)
    g5 = np.conj(C['CRL'])/(4*math.sqrt(2)*GF)
    term = AR*np.conj(g6)+AL*np.conj(g5)
    I = MU_NORM*8*e*term.real
    Imax = MU_NORM*8*e*(abs(AR*g6)+abs(AL*g5))
    photon = MU_NORM*32*e**2*(math.log(MMU**2/ME**2)-11/4)*(abs(AR)**2+abs(AL)**2)
    return dict(B_contact=Bcontact, B_photon=photon, interference=I,
                phase_envelope_interference=Imax, B_subset=Bcontact+photon+I,
                independent_phase_envelope=[Bcontact+photon-Imax, Bcontact+photon+Imax])


def main():
    y = yukawa(YE)
    s, cp, cL, cR, Beg, C, B3 = rates(y, 5000.)
    out = dict(scope='Tree contact and neutral-scalar on-shell one-loop dipole. Contact+dipole subset includes exact fixed-fit interference. No q=0 monopole extrapolation.',
               normalization=dict(Gamma_mu_GeV=GAMMA_MU, measured_over_tree_normalization=MU_NORM),
               masses_MeV=(s*1e3).tolist(),
               benchmark_5TeV=dict(B_meg=Beg,B_3e_contact=B3,R_contact=B3/Beg,
                  chi2_only_dipole=dipole(cp,5000.,[2])[2],chi1_only_contact=contact(cp,5000.,[1])[1],
                  contact_dipole=contact_dipole(cL,cR,C,B3)),
               dipole_dominated_R3e=ALPHA/(3*math.pi)*(math.log(MMU**2/ME**2)-11/4),
               dipole_only_conversion_ratio_Al=conversion_over_meg('Al'))
    out['scale_scan']={}
    for M in (3.,4.,5.,6.,7.5,10.):
        _,_,cl,cr,b,C_,bc=rates(y,M*1e3)
        out['scale_scan'][str(M)]=dict(B_meg=b,B_3e_contact=bc,R_contact=bc/b,
                                       B_contact_dipole=contact_dipole(cl,cr,C_,bc)['B_subset'])
    fB=lambda m,idx: rates(y,m*1e3)[idx]
    out['common_scale_TeV']=dict(
       MEGII_90CL=brentq(lambda m:fB(m,4)-1.5e-13,2,20),
       SINDRUM_contact_diagnostic_90CL=brentq(lambda m:fB(m,6)-1e-12,1,20),
       MEGII_design_90CL=brentq(lambda m:fB(m,4)-6e-14,2,20),
       Mu3eI_contact_SES=brentq(lambda m:fB(m,6)-2e-15,2,30),
       Mu3eI_contact_zero_background_90CL=brentq(lambda m:fB(m,6)-(-math.log(.1))*2e-15,2,30),
       Mu3eII_contact_declared_BR_target=brentq(lambda m:fB(m,6)-1e-16,2,30))
    targets=s
    Rs,Begs,B3s,tries=sample_mass_preserving(np.random.default_rng(2026),1500,5000.,targets)
    out['robustness_contact_only']=dict(n=len(Rs),tries=tries,seed=2026,
        prior='Off-diagonal magnitudes log uniform [1/3,3]; all phases uniform; diagonal magnitudes fitted to printed-input singular values and accepted in [1/3,3]; no neutrino or experimental cuts.',
        R_percentiles_5_16_50_84_95=np.percentile(Rs,[5,16,50,84,95]).tolist(),
        frac_above_10x_dipole=float((Rs>10*out['dipole_dominated_R3e']).mean()),
        frac_below_dipole=float((Rs<out['dipole_dominated_R3e']).mean()))
    np.save(Path(__file__).with_name('robustness_v9.npy'),np.vstack([Rs,Begs,B3s]))
    Path(__file__).with_name('lfv_letter.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='scale_scan'},indent=2))

if __name__=='__main__':main()

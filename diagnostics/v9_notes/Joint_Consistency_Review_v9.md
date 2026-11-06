# Joint consistency review and submission scope — v9

Gauhar Abbas — 11 October 2026

## Files and version boundary

This revision continues the four supplied v8 packages. The two newly attached v5 archives are earlier versions: their main sources were compared with the working drafts, and they were not used to roll back later construction, matching or audit material. The deliverables include four revised manuscripts, a one-page PRL cover letter, a 100-word PRL justification, vector figures, numerical inputs, Python sources and JSON checks. The Letter uses REVTeX 4.2; the three longer manuscripts remain companion full articles.

The review verifies the declared EFT calculations and printed witnesses. It does not certify a complete microscopic TC–DTC–DQCD completion, every historical numerical diagnostic, a global phenomenological fit, or publication acceptance.

## Common physical boundary

| Quantity | Shared value or convention |
|---|---|
| Electroweak VEV | 246 GeV |
| Composite decay constant | 2 TeV; xi = 0.015129 |
| Higgs pole | 125 GeV; curvature mass 125.063701 GeV |
| Top boundary | 172.5 GeV on-shell input in the Higgs calculation |
| Higgs matching scale | 1 TeV, canonical Zh = 1 |
| Top parent | MT = MD = ell = 6.398573981 TeV; r = 3.485779827 TeV; yT = yD = 2.373605268; FT = FD = 1 TeV |
| Top subtraction | The positive top curvature contribution and twice-subtracted pole response are separate; they are not counted twice |
| Fermion pole response | −15.929393 GeV² |
| Spectral matching | L = 3.315732847 TeV; moment product 2.201908452 |
| LFV witness | LambdaD = Mr = 5 TeV, canonical unmixed singlets and degenerate scalar/pseudoscalar poles |
| Flavour top target | yt(1 TeV) = 0.8616 is a distinct running input, not a second pole subtraction |

Conversion of the running flavour target to the top on-shell input requires QCD, scheme and new-sector threshold matching. This conversion has not been computed; the manuscripts explicitly identify the interface rather than assert a completed joint top fit. The LFV normalization scale and scalar poles are distinct from the 1-TeV Higgs matching scale.

The common precision comparison uses Table 10.8 of the PDG 2025 electroweak review, carrying the 15 April 2026 update: means (0.021, 0.040, 0.008), standard deviations (0.096, 0.120, 0.092), and correlations (0.91, −0.62, −0.83) for ST, SU and TU. The full three-parameter covariance is used, not a U = 0 ellipse.

The calculated point is

**(S,T,U) = (0.071522659, 0.063359564, 0.001349597), chi² = 0.915201.**

It lies inside the stated three-dimensional 95% ellipsoid, whose threshold is 7.814728. This is compatibility of a conditional benchmark. Nonuniversal vertices, finite dynamical-singlet thresholds, direct-search likelihoods and a complete flavour/SMEFT fit are not thereby tested. Strong-sector S = 0.066467357; the total includes the declared infrared/finite terms.

## Collective protection and operator closure

The original Gaussian chain is not promoted to a radiatively closed complete theory by its determinant alone. The revised collective paper owns a soft-mixing parent with complete triplets, radial fields, a general norm potential, universal bridges and a relative phase symmetry. The full 9-by-8 mass matrix, including the down spectator, independently reproduces the one-loop trace with no angular scalar pole. Kinetic, Yukawa and norm counterterms still run.

The proof is an isolated-parent statement in a mass-independent scheme. Local ultraviolet residues cannot contain inverse masses: soft bridge insertions therefore do not produce the proposed hard shortcuts or an independent renormalizable angular quartic. Finite mediator operators remain allowed and are listed separately. Direct full-matrix tests give orientation-independent determinants on the four restoring cuts, while an explicit nonuniversal bridge breaks this property. Equal low-energy coefficients, equal norms, scalar degeneracy and matching choices are not upgraded to symmetry consequences.

A sufficient dynamical-singlet extension is now specified. It retains the soft left bridge and introduces a universal hard single-orientation right bridge, the complete singlet/norm potential and an accidental global U(1)chi. The genuine cuts become ell = 0, kappa_chi = 0, yT = 0 and yD = 0. **Zero singlet VEV with nonzero hard coupling is not a restoring cut**, because virtual singlets remain.

The soft phase-mass coefficient Bchi must carry selector charge (0,0,12,6,0) and zero Dirac-protector charge. Fixing it explicitly breaks the full selector. Its charge is already in D37, so it does not enlarge the condensate dressing subgroup or change its quotient. Finite flavour thresholds and the microscopic origin of this soft spurion remain open.

A tree singlet boundary with VEV 4.35 TeV, lambda_chi = 0.330294623, m_chi² = 0 and Bchi = −6.25 TeV² has a positive Hessian diag(25,25) TeV². This supplies a consistent tree witness for 5-TeV scalar and phase masses, not a prediction of pole degeneracy or full RG stability. The singlet norm portals must be retained under running. The physical Higgs–singlet portal belongs to the explicit breaking budget.

Gauge and flavour interactions modify the isolated argument. Electroweak gauging generates angular operators; hard flavour tensors can generate shortcuts; finite thresholds need matching. The bridge U(1)chi has vanishing colour and gravitational sums but nonzero electroweak and cubic sums for the displayed top fields. It is distinct from the gauged Dirac protector U(1)X = 3(B−L), and is not claimed to be an anomaly-free gauge symmetry. Relative-phase constituent cancellation is a sufficient displayed assignment, not a complete UV selector construction.

The common Higgs numbers remain the frozen-bridge calculation. Additional dynamical-singlet loop contributions have not been silently folded into them.

## Spectral contribution and novelty

The spectral manuscript credits established sum-rule inequalities, generalized Stieltjes positivity and the known three-resonance counterexample. Its claimed addition is the observable certificate that combines different current directions, normalization repair, weighted spectral defects and physical pole allowances under explicitly declared hypotheses.

The new worked comparison demonstrates an actual change of inference: imposing a same-current, zero-defect bound gives a false exclusion; keeping the measured overlap and budgets gives a nontrivial surviving constraint. An independent script verifies the mapping to the known counterexample and the two-current example. The near-equality with the earlier Marzocca–Serone–Shu resonance coefficient is explained by the two-pole shape factor; it is not presented as an unrelated new universal inequality.

Asymptotic freedom and a finite-energy ansatz do not prove a global positivity class. Tail matching, its cumulative defect and the required infrared hypothesis remain separate. The application to TC–DTC–DQCD is conditional. The collective paper owns the construction and radiative analysis; the spectral paper owns the counterexample, hierarchy and corrected conditional bounds. Repeated numerical applications are identified as shared interfaces.

## Flavour, LFV and numerical corrections

For a single-orientation top bridge, the exact mass-basis derivatives give

**eta3 = 0.648395241 and gs3tt = ga3tγ5t = 0.018181289.**

The older unsuppressed pseudoscalar vertex was incorrect. A phase redefinition generates kinetic currents and cannot remove the heavy–light mixing from the physical vertex. The diagonal top result is not a universal rescaling of all off-diagonal flavour vertices. Trading the singlet VEV against its coupling preserves a fitted mass matrix; it does not prove a flat scalar potential.

The charge audit independently enumerates the leading charged and degree-two neutrino operators. It verifies the sequential counts 67 → 14 → 2, scalar bilinear classes 10 → 2 → 0, selector order 5376, D37 order 56 and quotient order 96. All seven VEVs generate the full selector group, so an independent Dirac protector is necessary. The printed colour, weak and gravitational measure sums are reproduced; these congruences are not a complete modern discrete-anomaly classification or constituent completion.

Direct SVD of the printed six-decimal coefficients and precise insertions gives charged masses (0.495758796, 104.651014790, 1778.840675297) MeV and neutrino masses (0.000266933, 0.008669828, 0.050496214) eV. The displayed oscillation angles and splittings have been aligned with these values. The fit is a supplied witness, not a reproduced global optimizer. The advertised original Mathematica N1–N10 suite, unrounded fit vectors and multi-start optimizer were not supplied; unsupported execution and uncertainty-band claims have been removed.

Both LFV rates now use the measured muon width. At 5 TeV:

| Quantity | Result |
|---|---:|
| On-shell neutral-scalar B(mu → e gamma) | 2.338860475 × 10⁻¹⁴ |
| Tree contact-only B(mu → 3e) | 4.031753833 × 10⁻¹⁵ |
| Contact-to-radiative ratio RC | 0.172381118 |
| Photon-dipole three-body term | 1.433014297 × 10⁻¹⁶ |
| Fixed-witness contact–dipole interference | −2.435257062 × 10⁻¹⁶ |
| Contact-plus-dipole subset | 3.931529557 × 10⁻¹⁵ |

Opposite initial-muon chiralities do not eliminate interference. The normalization and conjugation are explicitly matched to the Kuno–Okada operator convention. Channel-only rates are coherent diagnostics, not additive branching fractions. Fixed orientation guarantees cancellation of the internal-mass-enhanced dipole for degenerate poles; mixed orientation permits a surviving term but need not prevent alignment or cancellation. External flips survive.

Neither RC nor the contact-plus-dipole subset is the full physical three-body prediction. Momentum-dependent photon charge-radius terms, allowed Z/box contributions and low-energy running remain to be matched. A q² = 0 approximation with light internal leptons is insufficient for physical three-body decay or nuclear conversion. Unsupported inclusive predictions and unique-signature claims were removed. The aluminium number retained in the Letter is explicitly dipole-only.

The deterministic 1500-point scan reproduces median RC = 0.146411486 and 16th/84th percentiles 0.049761726/0.441040714. It uses a declared prior, charged-mass preservation and no neutrino or experimental cuts. It is not a probability distribution over viable theories.

At matched 90% conventions the conditional MEG-II design reach is 4.42 TeV and the background-free contact-only Mu3e-I reach is 4.92 TeV. The SES crossing is instead 5.46 TeV. Two predicted events near 5 TeV establish neither discovery nor 90% exclusion. All three-body reaches are operator-subset sensitivity estimates pending full matching and experimental acceptance/background treatment.

The canonical neutrino dimension is a benchmark choice. A two-loop coefficient ratio does not prove absence of walking. The 5.29 × 10¹¹-GeV walking example lies below the stated formal hypercharge pole, as does the quoted MFV ceiling; the earlier contradictory exclusion language has been corrected. Formal extrapolations remain field-content and threshold dependent.

## Checks and practical submission assessment

The shared audit scripts pass for the corrected top derivatives, finite dipole integral representations, mass-basis rephasing, orientation identity, lepton/neutrino witness, exact charge enumeration, singlet trace/contraction and tree Hessian. The primary physical Higgs, closure and radial checks pass in their declared scopes. The spectral source/analytic audit passes. There are no unresolved labels/citations or duplicate labels in the four revised sources. The final manuscripts were compiled and their page layouts inspected; vector plots were regenerated and the crowded Letter legend corrected.

The Letter is three pages in PRL format; the cover is one page. A separate exactly 100-word significance statement is included. Data availability and substantive AI-assistance disclosures have been added. The cover distinguishes the Letter from the proposed longer PRD companion and requests coordinated handling if both are submitted. Submission history and author approval remain matters for the actual submission; no manuscript has been submitted by this workflow.

**Editorial assessment:** formatting and coherent presentation do not establish the significance threshold for PRL. The Letter's defensible result is a narrow orientation/operator comparison, and part of its central physical LFV matching remains unfinished. PRL editorial risk is therefore substantial. Completing the physical momentum-dependent amplitude would strengthen its experimental claim more than stronger novelty language. The longer papers have independent conditional theoretical content, but claims of a fully viable UV completion or a complete global fit would require additional calculations. The revision deliberately states these limits rather than guarantee PRL or PRD acceptance.

## Primary literature/data links used for the submission conventions

- APS PRL author guidance: https://journals.aps.org/prl/authors
- PDG electroweak review: https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf
- MEG-II limit: https://arxiv.org/abs/2504.15711
- Mu3e technical design: https://arxiv.org/abs/2009.11690
- Kuno–Okada operator conventions: https://arxiv.org/abs/hep-ph/9909265
- Lavoura finite dipole functions: https://arxiv.org/abs/hep-ph/0302221
- Pascoli–Zhou flavon LFV predecessor: https://arxiv.org/abs/1607.05599
- Koivunen–Raidal flavon analysis: https://arxiv.org/abs/2305.00014

Reproduction commands, scope notes and dependencies are supplied in each source archive's README. Earlier notes and examples are labelled under diagnostics; they are not the authority for the v9 benchmark.

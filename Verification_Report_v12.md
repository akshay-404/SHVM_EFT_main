# Four-manuscript consistency and reproduction review — v12

Gauhar Abbas · 12 October 2026

The spectral and collective manuscripts have been synchronized with the physical finite-photon LFV calculation in the v11 Letter and flavour paper. All four manuscripts are now v12. Fresh execution reproduces the supplied primary Higgs, LFV, spectral and construction calculations. This review establishes consistency within the declared EFT truncations; it does not establish the missing microscopic matching assumptions or a complete experimental likelihood.

## Changes made

**Spectral manuscript.** Appendix A now quotes the physical three-body rate, identifies the earlier contact and contact-plus-dipole results as separately labelled diagnostics, and references the current shared data. A new subsection next to the precision application lists the independently required conditions for a numerical exclusion. The physical plot now displays the known finite fermion pole allowance explicitly and marks the stronger retained-top bound only at the common point. Its f-curves hold matching data fixed and are labelled diagnostic slices.

**Collective manuscript.** The dispersive definition of `S_str` now consistently denotes the pure current moment, `4*pi*xi*I`. Finite matching terms are added to total S separately. The shortened flavour discussion quotes the current finite-photon rates. The order-three proposition specifies the direction of its infrared sign crossing. Its order-two obstruction distinguishes a tail with a positive exponent extending indefinitely from a threshold-dependent tail, which requires its own integral-sign test. The precision–mass corollary now states finite moments, equal embeddings, f=f_h, xi<1/2, positive kinetic/gauge matching, and the complete physical pole remainder explicitly.

**Long SHVM flavour manuscript.** The statement equating the running top texture target to the pole input is corrected. Stationarity eliminates a fitted first derivative; it does not cancel the physical top contribution to Higgs curvature. The remaining determinant convexity contribution is stated at the declared accuracy. The three-parameter 95% cutoff is consistently rounded to 7.815. The active common-scale LFV plot uses the same physical amplitude implementation as the Letter.

**Letter and supplement.** Their independent low-energy calculation and substantive v11 results are retained. They now identify the common v12 boundary and explicitly distinguish its Higgs-parent assumptions from the LFV pole-alignment and residual-matching requirements. The supplement records the different threshold and electromagnetic-coupling conventions. The four-page Letter, nine-page supplement and PLB cover draft rebuild successfully.

## One physical boundary, with explicit conventions

Identical copies of `shared/joint_boundary_v12.json` are included in all four source trees. Its LFV section is also available as `shared/lfv_boundary_v12.json`.

| Input or convention | Common declaration |
|---|---|
| Higgs and top pole inputs | 125 GeV and 172.5 GeV |
| Higgs angular scale and electroweak VEV | f=2 TeV; v=246 GeV |
| Higgs parent matching | mu_H=1 TeV; one SM-top subtraction; Z_h=1 |
| Retained Higgs pole response | Twice-subtracted finite fermion self-energy, counted once |
| Charged-singlet poles and flavour normalization | All eight canonical real poles M_b=Lambda_D=5 TeV |
| Lepton texture | Printed charged-lepton SVD witness; physical lepton masses in loops and phase space |
| Electromagnetic inputs | alpha_EW^-1=128 for the Higgs/current calculation; alpha_QED^-1=137.035999 for low-energy LFV |
| Running top texture target | y_t(1 TeV)=0.8616; its conversion to the pole input remains uncomputed |
| Precision comparison | PDG 2025 EW review, April 2026 release, Table 10.8, freely floating S,T,U |
| Additional matching | Explicit scenarios or residual requirements; not UV-derived error distributions |

The Higgs matching scale and flavour normalization/pole scale have different roles. No uncomputed running between them is asserted. Canonical pole residues are absorbed into the LFV Yukawa matrices. The aligned LFV poles are a physical matching boundary, not an equality of bare masses enforced by the collective top symmetry.

## Fresh numerical results

| Quantity | Reproduced result | Status |
|---|---:|---|
| Top input realized by the retained parent | 172.5000000000003 GeV | Leading on-shell matching |
| Higgs pole | 124.9999999999996 GeV | Declared finite-fermion truncation |
| Higgs curvature | 125.0637013424 GeV | Before the retained pole correction |
| Finite fermion pole remainder | -15.9293934629 GeV² | Counted once |
| Reduced gauge mass square | 10305.4319452 GeV² | Matched parent calculation |
| Positive determinant convexity term | 5335.49744824 GeV² | Fixed parent and truncation |
| Moment product I*Omega/f⁴ | 2.20190845242 | Stated aggregate spectrum |
| Pure current S | 0.06646735749 | Excludes other precision terms |
| Total retained (S,T,U) | (0.07152265860, 0.06335956435, 0.00134959673) | Oblique truncation |
| Correlated chi²_STU | 0.9152012150 | Three-parameter Gaussian comparison |
| Maximum chi² at the eight declared threshold corners | 2.6339548933 | Chosen ±0.01, ±0.01, ±0.005 scenarios |
| B(mu→e gamma) | 2.33885809159e-14 | Scalar one-loop dipole |
| B(mu→3e) | 4.34749992785e-15 | Scalar tree + finite-momentum photon |
| R=B_3e/B_eg | 0.185881304363 | Same amplitude and normalization |
| Conditional LFV ratio enclosure | [0.0864117634, 0.3585673309] | Imposed Higgs/local/residual budgets |

The covariance was checked against the [official PDG review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf), Table 10.8: means (0.021,0.040,0.008), standard deviations (0.096,0.120,0.092), and correlations (0.91,-0.62,-0.83). It is positive definite. The exact three-parameter 95% chi-square cutoff is 7.81472790325; the manuscripts use 7.815. The reference fit is an oblique-only test, and these numbers do not include uncomputed nonuniversal vertices or a full joint likelihood.

## Conditional spectral exclusions

The fresh numerator is N=466.623506817 GeV². With the matched global order-three hypothesis, coincident currents, the calculated positive gauge factor and the known finite fermion allowance E_f=15.929393463 GeV²:

- Discarding additional positive top information gives `S_str >= N/(125²+E_f) = 0.0298334897549`.
- Retaining P=5335.49744824 GeV² at the fixed parent gives `S_str >= N/(125²+E_f-P) = 0.0452793739551`.
- An independent moment route gives `S_spec*(3/2)/(I*Omega/f⁴) = 0.0452793739549`.
- The separately imposed Z_h≤1.05 and total remainder Delta≥-0.06*m_h² scenario gives 0.0268318997631. It is not derived from the limited matching scans.

The two-current repair diagnostic reproduces the surviving nonzero bound of approximately 0.00761696. Incorrectly identifying its two current embeddings gives the false 0.2 exclusion. The known three-state counterexample mapping, the two-pole coefficient comparison, noncommuting moment tests, normalization repairs and basis covariance reproduce.

An exclusion still requires independent proof of the global positivity class or weighted repair, physical current overlap, gauge allowance and entire one-sided pole remainder. The infrared partial-area condition and third-moment/OPE matching are independent of asymptotic freedom. A finite fractional-order scan does not prove global positivity; its near-1.968 coefficient is not used in the rigorous conditional bound. Scalar mixing without a controlled gap can invalidate a diagonal-curvature argument. The known positive term is a floor only within the stated truncation, with every extra negative correction included in the remainder.

## Checks actually executed

Fresh logs and JSON are included under `verification/`; current parent outputs are in `checks/v12/`.

1. Solved the shared parent for the top input, stationary vacuum and Higgs pole; recomputed fermion currents, leading scalar projection, Higgs couplings, matching corners and vacuum diagnostics.
2. Executed the symbolic full-matrix UV trace, restoring-cut determinants, gauge quartic polynomial and stated Weyl anomaly sums. The nonuniversal-bridge counterexample violates its restoring cut as expected. These checks verify the isolated parent and assigned fields, not every gauge/flavour completion.
3. Executed the norm-sector RG stability diagnostic. It does not establish the uncomputed angular, portal or global scalar stability conditions.
4. Executed the spectral analytic comparison, canonical-Hessian/basis-covariance checks, two-regime threshold calculation and finite fractional-order diagnostic scan.
5. Executed the printed charged/neutrino witness, discrete-charge/operator counts, dynamical-bridge audit and top-interface calculation. No unavailable Mathematica fit suite or absent multistart global fit is claimed.
6. Reintegrated the physical three-body Gram matrices at 160×160 and 320×320 quadrature resolution. The central rate changes by 4.8318e-8 relatively. Spin sums, antiparticle completeness, photon Ward identity, antisymmetrized Fierz identities, contact/dipole normalizations, finite timelike loop quadrature and exact-degeneracy mixing invariance pass.
7. Reexecuted pole/mixing/matching scans, the leading neutral-box majorant, Higgs-angle budget and deterministic enclosure checks. The missing inclusive-QED and full-parent thresholds remain imposed residual requirements.
8. Regenerated active numerical figures in all four manuscripts. The compact-spectrum figure is regenerated from explicitly historical data; static construction diagrams are supplied inputs. The spectral first figure is regenerated by the embedded PGFPlots code during the build.
9. Rebuilt all six PDFs through three LaTeX passes. Source references and citations resolve; no duplicate labels or overfull boxes remain. The collective build emits only the benign amsmath warning for retained TeX fraction syntax.
10. Rendered and reviewed all 164 pages using contact sheets, with detailed inspection of changed pages and figures. The geometry audit reports no off-page text.

## Reproduction

From the combined source root:

```bash
python3 -m pip install -r requirements.txt
python3 reproduce_all_v12.py
```

Optional all-page rendering and geometry audit:

```bash
python3 -m pip install PyMuPDF Pillow
python3 render_audit_v12.py
```

Individual source packages also include manuscript-specific commands in their README files. Version suffixes on unchanged implementations identify historical filenames, not older authoritative benchmark data. The current shared boundary and verification outputs take precedence over historical diagnostic notes.

## Remaining scientific scope

The revised package supports a conditional spectral article, a construction/radiative-analysis article, a flavour-EFT construction article and a focused LFV Letter with supplement. It does not establish a microscopic TC–DTC–DQCD spectral class, complete gauge/flavour finite matching, pole-to-running top conversion, inclusive LFV NLO, full scalar stability, cosmology or a global viability fit. Those limitations are retained in the manuscripts and must accompany any submission claim. No manuscript has been submitted by this workflow.

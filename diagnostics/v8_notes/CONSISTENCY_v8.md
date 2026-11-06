# Three-paper consistency revision (v8, 11 October 2026)

Series: Paper I (collective Higgs protection, top sector), companion spectral analysis
(current-correlator bounds), Paper II (flavour, LFV, Dirac neutrinos). All three now share one
benchmark boundary: f = 2 TeV, on-shell M_t = 172.5 GeV (equivalent to the running
y_t(1 TeV) = 0.8616 of Paper II's flavour fits), M_h = 125 GeV pole, PDG 2025 three-parameter
covariance; (S,T,U) = (0.0715, 0.0634, 0.0013), chi2_STU = 0.92.

## Inconsistencies found and fixed
Paper II still described Paper I's superseded top sector. The following were all re-derived for
Paper I v7's radiatively closed soft-mixing parent (script: PaperII interface/paper2_interface_v7.py).

| quantity | Paper II before | now (all three papers) |
|---|---|---|
| top sector | four-spurion chain of three SO(6) sextets | radiatively closed soft-mixing parent; chi_3 in the right bridge r = kappa_chi <chi_3> |
| eta_3 = dln m_t/dln<chi_3> | 0.764 | 0.648 |
| g_{s3 tt}, g_{a3 tt} | 0.0186, 0.0244 | 0.0182, 0.0280 |
| chi_3 portal bound | 4.2e-4 | 4.2e-4 (unchanged; depends only on m_h, xi, <chi_3>) |
| hypercharge strong-coupling scale | 1.1e8 GeV (Delta b_Y = 44/3 per sextet) | 1.4e15 GeV (parent Delta b_Y = 52/9 in total) |
| link consistency | Delta_nu >= 2.85, gamma_m <= 0.15 | Delta_nu >= 1.70, gamma_m <= 1.30 |
| link table, max gamma_m column | 0.12-0.30 | recomputed with the new Lambda_Y |
| epsilon_7 floor from Lambda_Y | 1.7e-10 (dominant) | ~1e-24; the tree-level floor is now operative |
| electroweak section | GHLS, Gfitter U=0, chi2_ST = 0.87 | correlator-based, PDG 2025 three-parameter, chi2_STU = 0.92 |
| series | "second of two companion papers" | three-paper series with mutual citations |

- Method validation: the same running code reproduces Paper II's old Lambda_Y = 1.1e8 GeV
  (1.0e8) and 7.3e8 GeV for messengers alone (6.7e8).
- Physical consequence: hypercharge perturbativity no longer excludes walking-branch link
  scales. The non-walking branch now rests on the supercritical reduced DQCD spectrum alone,
  as Paper II states.
- The parent's selector-charge assignment is flagged in Paper II as an open input; Paper I
  audits only the relative-phase anomalies.
- Paper I gains a subsection "Interface with the flavour paper" that derives eta_3, the
  couplings, the portal bound and the running, so Paper II cites derived numbers.
- Bibliography: companion titles corrected in all three; uncited GHLS, Gfitter and other
  stale entries removed.

## Presentation
Each paper opens its introduction with a highlighted "Main results" box (numbered key results
with numbers and section links, a scope line and a common series line). All dates are
harmonized.

## Checks
All three compile with no undefined references, no uncited or missing entries and no
duplicate labels. Overfull boxes: Paper I 0, spectral 0, Paper II 1 (9.7 pt, pre-existing).

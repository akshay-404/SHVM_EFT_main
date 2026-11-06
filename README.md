# SHVM flavour manuscript for PRD, v12

Primary source: Paper2_SHVM_flavour.tex. This is the revised longer flavour paper,
with its operator, anomaly, messenger and neutrino construction retained independently
of the PLB Letter. The LFV application now uses the same physical finite-photon
headline, pole/residue conventions, Higgs/top inputs and sensitivity definitions.

Current LFV boundary: shared/lfv_boundary_v12.json.
The old shared/joint_boundary_v9.json remains the unchanged collective-Higgs/flavour
construction audit; its LFV operator subsets are historical diagnostics. In case of
LFV differences the v11 boundary is primary. No broad viability or full global fit is asserted.

Run:
```
python3 shared/joint_consistency_v9.py
python3 shared/flavour_charge_audit_v9.py
python3 shared/dynamical_bridge_audit_v9.py
python3 interface/paper2_interface_v7.py
python3 shared/physical_lfv/verify_physical_v10.py
python3 shared/physical_lfv/controlled_domain_v11.py
python3 shared/physical_lfv/verify_domain_v11.py
python3 code/flavour_figures_v9.py
```
The v7/v9 filenames on unchanged construction/interface checks are historical.
Dependencies are in requirements.txt. The v11 LFV code includes the neutral-box
majorant, conservative amplitude enclosure and Higgs-angle caps. Unknown residual
matching budgets remain imposed requirements, not computed UV or inclusive-QED errors.
The figure generator retains its v9 filename but panel (a) now uses the physical
finite-photon rate; panel (b) and the separate-scale curve are labelled diagnostics.

Compile twice from this directory:
```
pdflatex -interaction=nonstopmode -halt-on-error Paper2_SHVM_flavour.tex
```
Original absent Mathematica N1-N10 fits or global likelihoods are not claimed to have run.

All four papers share `shared/joint_boundary_v12.json`. The v11 suffixes on unchanged LFV-domain scripts identify implementations, not a stale scientific boundary. Current execution logs are in the combined package.

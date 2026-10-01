# STAR cuts applied
#
# mikael.mieskolainen@cern.ch (2026)

1. Basic fiducial cuts are applied in `gencard.json`.
2. More involved fiducial cuts are applied in `src/Graniitti/MUserCuts.cc` at the code level, referenced by `"USERCUTS" : 1792394010(1)(2)` under the corresponding `gencard.json`.
3. Categorical fiducial cuts are applied over the central system mass `M` and forward proton `dPhi` at the Python code level.

All fiducial acceptances are automatically treated.

Alternatively, one could apply all (complicated) cuts at the Python level, but this would be computationally more heavy (CPU wasteful) for parameter tuning purposes.

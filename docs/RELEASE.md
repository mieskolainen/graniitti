## GRANIITTI release notes

mikael.mieskolainen@cern.ch

> *Version 1.5 beta (September 2026)*

This release includes new and updated physics models, more processes, large architecture changes, technology improvements, automated analysis & testing and bugfixes.

> This is a **beta** version, which means some (soft QCD) processes may be in R&D stage and features unstable, at least in terms of physics or technical parameters, or steering validation scripts. Run `icepacks` to see comparisons with the data. See also `docs/underwork.md`.


### General

- **[New]** General machinery for $k_T$-EPA $\gamma \gamma$ + MG5 amplitudes (matrix elements) and photoproduction processes across `pp`, `ee`, `ep`, `pA`, `eA` and `AA` (when applicable).

- **[New]** **UPC initial states** (advanced prototype): coherent & incoherent, EMD, screening convolutions, forward emissions at particle level, with `icepack/UPC` LHC data comparisons for validation and development. [R&D-level]


### Soft processes

- **[New]** End-to-end helicity coupled soft Pomeron models, with process syntax `(GP|XP|MP|TP)[RES|CON|RES+CON]`.

`GP`: Exchange with a continuous angular momentum $J = \alpha(t)$ with analytically continued LS basis production coupling structures (or direct helicity amplitudes) per resonance.

`XP`: Fixed spin exchange (adjustable), with a covariant basis LS operators and reduced local Jacob-Wick tensors.

`MP`: Fixed spin exchange (adjustable), with a 'minimal coupling' principle modulated with adjustable spin polarization $J_z$-vector (coherent sum) or density matrix (incoherent sum) per resonance in the chosen Lorentz frame.

`TP`: Tensor Pomeron model, with covariant rank-2 current coupling structures per resonance.

- **[New]** Multi-channel (1,2,3,...) and multi-exchange (P,O,...) matrix unitarization based screening (absorption) loop models & elastic scattering, including CNI (Coulomb Nuclear Interference) for low-$|t|$ and non-commutative Good-Walker structure.

- **[Update]** Tensor Pomeron amplitudes are now 10-100x faster than in the previous version, due to optimized tensor contractions and cache code.

- **[New]** Automated and fully modelcard adjustable (Pomeron, Reggeon, Odderon, Gamma)-pair exchange combinatorics for central production, with arbitrary production channels in superposition (e.g. Pomeron-Pomeron + Reggeon-Pomeron).

- **[New]** Sub-amplitude form factors `FF_transfer`, `FF_prod`, `FF_decay` can be adjusted in a modular way, explicitly per model, resonance and decay.

- **[New]** Direct helicity amplitude and LS-basis coupling steering options (production and decays).

- **[New]** Improved gamma-Pomeron photoproduction model (MP/XP/GP) with HERA vector meson data parameter extractions.

- **[New]** Bose and Fermi statistics coherent symmetrization option for Jacob-Wick cascade decays.

- **[New]** Forward system N* excitation q(diquark) string skeleton output for Pythia fragmentation.

- **[Extended]** Tensor Pomeron model amplitudes (axial vectors) and cascaded decays, Drell-Söding + rho/omega amplitude `TP[PHOTO]` added, all parameters made adjustable.

- **[Extended]** New soft form factors and adjustable control per final state meson/baryon pair family.

- **[New]** Entropy and quantum entanglement diagnostic metrics (optional) computed during the cross-section integration stage [MP/XP/GP-models] [R&D-level].


### Hard processes

- **[New]** Rewritten, fully automated C++ HELAS amplitude export from MadGraph 5 to GRANIITTI [R&D-level].

- **[New]** Kinematically more accurate on-shell MG5 amplitude matching procedure for kT-EPA, and Durham QCD partonic processes.

- **[New]** In-library loop amplitude processes $\gamma\gamma \mapsto \gamma\gamma$ (kT-EPA) and $gg \mapsto \gamma\gamma$ (Durham QCD).

- **[Finished]** Generic $SU(3)$ color algebra projector for interfacing MG5 amplitudes with Durham QCD color singlet requirement, with color flow sampling in ($gg$, $ggg$, $q\bar{q}$ ...) for Pythia.

- **[Improved]** Durham process screening gluon loop numerical integration, Sudakov veto and Shuvaev transform numerics.

- **[Finished]** Durham $\chi_c(0,1,2)$ and perturbative meson pair amplitudes implemented end-to-end.

- **[New]** Z-boson and heavy vector mesons via hard-photoproduction (gamma-UGD).

- **[New]** Hard-DPDF single/double diffraction with MG5 amplitudes and Pythia interfacing.


### Model parameters

The baseline soft model parameter inference strategy has been improved:

0. Extract coupling and trajectory `starting` values based on Donnachie-Landshoff fits `develop/tools/DL_couplings.py`.

1. Extract vector meson photoproduction parameters based on HERA data `develop/tools/HERA_couplings.py`.

2. Fit the eikonal Pomeron (screening model) parameters against elastic proton-(anti)proton scattering HEPData.

3. Fit central production parameters against STAR and CMS exclusive meson and baryon pair production HEPData [R&D-level].

Where steps 2. and 3. are implemented using the new HPC-distributed `icetune` machinery.

Other model parameters are from the literature or in their ansatz, or in placeholder values. The current parameter optimization effort is concentrated on `GP` model, then later `XP`, `MP` and `TP`.

### Technology

- **[New]** Neural spline coupling normalizing flow MC integrator, enabled with `-g NEUROJAC`, with vectorized standalone code and optional `libtorch` backend. The unweighted event generation efficiency, with `<C>` phase space in particular, can improve significantly (process dependent), because the network goes beyond dim-by-dim factorized importance sampling marginals of VEGAS [R&D-level].

- **[Updated]** The phase space `<C>` class uses now the system invariant mass as the user facing generation range like `<F>` option, which is easier, instead of intermediate transverse momenta.

- **[Updated]** Improved floating point and algorithmic treatments of numerical loop integrals and MC integrators.

- **[Updated]** Improved sampling efficiency for multi-body cascaded amplitudes.

- **[New/Extended]** Python tools (`icetune`, `iceplot` ...) now CERN Condor compatible for HPC-distributed simulation-based global optimization inference against HEPData.

- **[New]** Central production amplitude fitter `ampfit` based on interpolating the complex scattering amplitudes with autograd based gradient descent. Available as an optimizer option in `icetune`, to be used once the global optimization has been run [R&D-level].

- **[New]** The icetune based parameter inference handles covariance structure of the observables of experimental correlated systematics [R&D-level].

- **[Extended]** Spherical Harmonics moment expansion and inversion tool `fitharmonic` now supports both central only measurements (ALICE,LHCb) and forward proton tagged (ATLAS,CMS,STAR), including a new steering card format. Future interfaces for `neural networks` (e.g. with `libtorch`) beyond N-dim hypercell histograms, and improved uncertainty propagation [R&D-level].

- **[New]** Any .json file parameter can be changed via cli `--set` command, very practical for changing e.g. production or decay helicity / LS coupling values.

- **[Update]** The particular soft production channels, such as Gamma-Pomeron, are now defined under `modeldata/TUNE0/RES/*.json` per resonance (e.g. old syntax `yP` deprecated).

- **[Extended]** All technical parameters are now available under `NUMERICS.json`, model parameters under `GENERAL.json`, and other json-files under the modeldata.

- **[Extended]** More generic fiducial cuts machine via steering cards.

- **[Extended]** Automated, re-organized unit, integration and physics process validation test suite with `icepacks`. [R&D-level]

- **[Update]** Multiple bugfixes (e.g. amplitude signs, kinematic and coupling normalizations, logistics issues).

- **[Update]** CMake compilation, docs, steering and syntax updates.

- **[Update]** Case study analysis scripts (e.g. eikonal/elastic) changed to Python.

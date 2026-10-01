# *DOCS*
## MadGraph 5 to GRANIITTI conversion [R&D]

#### 30/09/2026
mikael.mieskolainen@cern.ch

`MG2GRA` generates and installs amplitudes from the [process registry](../develop/MG2GRA/processes.json). Do not edit generated C++ manually.

## Contents

[Build](#build) · [Add processes](#add-processes) · [Physics settings](#physics-settings) · [Remove processes](#remove-processes)

## Build

```bash
conda activate graniitti
source install/setenv.sh
./develop/MG2GRA/run.sh --mg5 /path/to/MG5_aMC/bin/mg5_aMC
cmake -S . -B build
cmake --build build -j4
```

Generated code goes to `include/Graniitti/Amplitude/MG5/` and `src/Amplitude/MG5/`, cards to `MG5cards/`. Changed outputs are preserved as `._old`. Rebuild after regeneration.

```bash
./develop/MG2GRA/run.sh --validate-only
```

This checks installed consistency, not numerical or physics agreement.

## Add processes

Standalone Durham process:

```bash
./develop/MG2GRA/run.sh --mg5 /path/to/MG5_aMC/bin/mg5_aMC \
  --add gg_ttbar --process "g g > t t~ QED=0" --model sm --projection durham
```

Use `--projection photon` for a concrete `a a` subprocess. Standalone final-state PDG signatures must be unique within each projection. Reordered particles are accepted when decay topology and PDGs match unambiguously.

Subprocess family:

```bash
./develop/MG2GRA/run.sh --mg5 /path/to/MG5_aMC/bin/mg5_aMC \
  --add-family MG5_PP_ZJJ --model sm --projection parton --channel Zjj \
  --definition "define p = g u d s c b u~ d~ s~ c~ b~" \
  --definition "define j = g u d s c b u~ d~ s~ c~ b~" \
  --family-process "p p > z j j, z > mu+ mu-"
```

Repeat `--family-process` for channels. Parton families serve `IPp` and `IPIP`, photon families serve `yy`, `yy_DZ` and `yy_LUX`. For new UFO models, add `--model-import /path/to/model`. Particle names, PDG codes and pole parameters are read from that UFO. Add new external particles to the GRANIITTI particle card with their spin, charge and color representation.

## Physics settings

| Setting | Requirement |
|---|---|
| Model | Choose the UFO restriction before generation. Cards cannot restore removed couplings. Supplied SM models retain charged lepton masses |
| Masses | Change independent SLHA parameters. The evaluated model controls amplitudes and phase space. Signed Majorana masses remain signed in HELAS and use positive physical poles in phase space. Conflicting `@PDG` overrides are rejected |
| Partons | Match `p` and `j` to the PDF flavour scheme. Incoming partons must be massless |
| Orders | `QED=2` is an upper bound and `QED==2` is exact. Coupling powers come from generated diagrams. Tree amplitudes retain interference between different orders. Include `QCD=0` for pure electroweak processes |
| Decays | Specify them in MG5. Decay chains select resonant diagrams, unrestricted final states can include nonresonant interference |
| Finite W width | The supplied $\gamma\gamma\to W^+W^-\to4\text{ leptons}$ amplitudes use `sm_cms`, our Standard Model configuration with MadGraph's complex-mass scheme enabled |
| Couplings | `QED_alpha: ZERO` requests the optional $\alpha(0)$ override configured with `--alpha-charge`. Charge aliases and mixed couplings are updated together. Without this setting, the UFO electromagnetic scheme is retained. `MG` and `LL` retain the generated scheme. PDFs supply running $\alpha_s$ |

The [complex-mass scheme](https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/ComplexMassScheme) uses $\mu_W^2=M_W^2-iM_W\Gamma_W$ consistently in propagators and electroweak couplings to preserve gauge identities. It is already enabled for the supplied WW family. When registering a new family and model, `--complex-mass-scheme` enables it during MadGraph generation.

Exact finite $N_C=3$ color contractions determine matrix elements. Durham projects the incoming gluons onto a singlet. Shower flows only assign color tags, including the signed LHE convention for sextets.

Loop generation requires `--type loop --projection durham`, an explicit `noborn` selector and a Fortran compiler for MadLoop. This adapter supports colorless final states with one incoming gluon trace and one QCD/QED amplitude order. Coupling rescaling requires positive reference couplings. The current registry contains tree processes only.


## Remove processes

```bash
./develop/MG2GRA/run.sh --mg5 /path/to/MG5_aMC/bin/mg5_aMC --remove gg_ttbar
./develop/MG2GRA/run.sh --mg5 /path/to/MG5_aMC/bin/mg5_aMC --remove-family MG5_PP_ZJJ
```

Old outputs become `._old`. Shared routines are regenerated for remaining processes. Rebuild afterward.

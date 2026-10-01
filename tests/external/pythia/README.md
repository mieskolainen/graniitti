# Pythia showering and hadronization [R&D]

## Build and run

After the [common setup](../../README.md), run from the repository root:

```bash
bash tests/external/pythia/build_pythia.sh
NEVENTS=10000 LOOPSCREEN=0 WEIGHTED=1 bash tests/external/pythia/cases/hard_pomeron_z/run.sh
bash tests/external/pythia/cases/hard_pomeron_z/analyze.sh
```

`hard_pomeron_z`, `yy_mumu_dissociative` and `durham_qcd` provide generation (`run.sh`) and analysis (`analyze.sh`). Scripts resolve repository paths, build drivers as needed and save LHE and HepMC3 files under `output/`. Keep generation and analysis tags consistent (`TAG`, `GG_TAG`, `UUBAR_TAG`, `IPP_TAG`, `IPIP_TAG` and native Pythia tags).

`PYTHIA_DIR` selects the installation. `JOBS` sets build concurrency (1 to 4, default 4). Drivers store absolute library and XML paths, with `PYTHIA8DATA` overriding the XML directory at runtime.

## LHE conversion

```bash
bash tests/external/pythia/drivers/lhe_converter/run.sh CARD_OR_LHE TAG [GRANIITTI_OPTIONS]
```

Cards default to 20 events, uncompressed LHE files to full conversion. `NEVENTS` overrides the count (first N events for LHE input). `SEED` controls both generators. `PYTHIA_CMND` selects settings, or Pythia defaults when empty. Streaming input accepts only weight strategies $\pm3$ and $\pm4$ because others can resample source events.

| `CONVERTER_MODE` | Behavior |
| --- | --- |
| `auto` (default) | Full showers for diffractive records with hard metadata and reduced hard records, isolated fragmentation for complete records without hard metadata |
| `fragment` | Preserve the supplied color state and disable ISR, the default for EPA and Durham cases |
| `shower` | Rebuild remnants and follow the card's ISR and primordial transverse momentum settings |

Both modes disable MPI and preserve realized nuclear breakup products.

Low mass remnant strings use cluster collapse with tagged proton recoil, preserving the hard system and event weight. Proton attributes `graniitti_recoil_dpx`, `graniitti_recoil_dpy`, `graniitti_recoil_dpz` and `graniitti_recoil_de` record the recoil. Diffractive variables refer to the input proton.

Hard fractions outside Pythia's remnant support trigger fragmentation of the complete source record without ISR, flagged by `graniitti_pythia_remnant_fallback` and `graniitti_low_remnant_isr_disabled`. Native Pythia SD excludes this region, whose cross section contribution must be distinguished.

## Cross section comparisons

iceplot normalizes converted LHE events using HepMC3 cross sections and relative weights. Native GRANIITTI weights retain their normalization. Absolute SD comparisons require common muon and tagged proton cuts, matched PDFs and flux normalization, and consistent MPI and screening.

Native Pythia covers wider proton momentum loss than the GRANIITTI cards and includes first emission matrix element matching for Z production. `shower.cmnd` limits LHE showers to the supplied hard scale (`SpaceShower:pTmaxMatch = 1`). Matching integrated rates does not ensure agreement at high $p_T$.

## Z merging

```bash
bash tests/external/pythia/cases/z_merging/run.sh
```

`icepack/HARDPOM/z_merging` combines SD Z and Z+jet matrix elements using Pythia CKKW-L kT merging, a 5 GeV parton cut and merging scales of 10, 15 and 20 GeV. iceplot compares absolute dimuon spectra and hadronic kT jet resolution against an unmerged shower and native Pythia SD with common proton and muon cuts. Outputs: `figs/iceplot/z_merging` and `tmp/z_merging/iceplot.json`.

| Option | Effect |
| --- | --- |
| `--tag NAME` | Separate events, logs and plots |
| `--events N` | Set all sample sizes |
| `--native-events N` | Override the native Pythia sample size |
| `--reuse` | Regenerate only changed or incomplete samples, checking arguments, executables, steering, model cards and output contents |

### Weights and normalization

Events vetoed by CKKW-L keep zero weight and are not reshowered to obtain acceptance. `graniitti_merging_weight` multiplies the preserved `graniitti_lhe_weight`. The final HepMC3 cross section includes weighted merging acceptance. Events outside native remnant support receive zero merging weight and `graniitti_merging_remnant_veto=1`, with excluded source weight reported separately.

Multiplicities use independent shower seeds and sum absolute cross sections, including samples with all merging weights zero. Normalization follows the [Pythia CKKW-L prescription](https://pythia.org/latest-manual/CKKWLMerging.html). Other prescriptions lack implemented weights and are rejected. This pilot uses diffractive full showers without MPI, screening or NLO accuracy.

Cross section uncertainty combines independent source integration error and weighted acceptance variance. For source weights $w_i$, merged weights $m_i$, $W=\sum_i w_i$ and $R=\sum_i m_i/W$, $\mathrm{Var}(R)=N\sum_i(m_i-Rw_i)^2/[(N-1)W^2]$ for $N>1$. Merging scales share source events, so their errors are correlated.

## Checks

After building both drivers:

```bash
pytest -q -s -rP tests/technical/drivers/test_pythia_lhe.py tests/technical/drivers/test_pythia_zmumu.py tests/technical/drivers/test_merging.py
```

The matching checker accepts continuum dimuons and merged weights. `--require-z-peak` adds a Z resonance check.

# Development tools

Run these commands from the repository root after setting up the environment:

```bash
conda activate graniitti
source install/setenv.sh
```

| Tool | Purpose |
| --- | --- |
| [DL_couplings.py](DL_couplings.py) | Derive MP, XP and GP continuum couplings from DL hadron ratios and the active SOFT model proton couplings |
| [HERA_couplings.py](HERA_couplings.py) | Fit vector meson photoproduction data, optionally including LHCb, and derive elastic and proton dissociation parameters together |
| [resonance_card_builder.py](resonance_card_builder.py) | Construct GP or XP resonance vertices in the LS or helicity basis |
| [analytic_Jz.py](analytic_Jz.py) | Derive GP couplings for pure $\lvert J_z\rvert$ sectors at an integer exchange pole |
| [ls_helicity_conversion.py](ls_helicity_conversion.py) | Compute symbolic Jacob-Wick LS to helicity transformations |
| [jpsi_gamma_helicity.py](jpsi_gamma_helicity.py) | Convert $\chi_{cJ}\to J/\psi\gamma$ multipoles into helicity rows |
| [tuneplot.py](tuneplot.py) | Combine icetune linear and logarithmic plots into a titled PDF, using `pdflatex` |

Every tool accepts `--help`. The scripts also support `python -m develop.tools.<name>` without the `.py` suffix.

## Update TUNE0 couplings

Derive continuum couplings from the current SOFT model fit:

```bash
python develop/tools/DL_couplings.py --tune-general modeldata/TUNE0/GENERAL.json --push
```

DL preserves the relative amplitudes and phases by default. Use `--shape lowest` to initialize the lowest allowed LS tensor. The active SOFT model fixes the absolute proton couplings and the DL coefficients fix the hadron to proton ratios. TP continuum couplings use a separate convention.

Fit HERA and LHCb measurements and update photoproduction parameters:

```bash
python develop/tools/HERA_couplings.py \
  --lhcb-fit \
  --tune-general modeldata/TUNE0/GENERAL.json \
  --output-dir "tmp/HERA_$(date -u +%Y-%m-%d_%H-%M-%S_UTC)" \
  --push
```

Both commands preview changes and ask for confirmation before writing, with backups of changed cards. Run both again whenever the physical SOFT model proton coupling or form factor changes. These downstream calculations are currently separate from `icetune --push`.

HERA writes `fit.json`, `cards.json` and comparison plots in `fit_comparisons.pdf`. It updates MP, XP and GP photoproduction couplings and the associated `photoprod_diss` rows together. Omit `--lhcb-fit` for HERA inputs alone. Replace `--push` with `--write-tune <new-directory>` to create a separate candidate tune, or with `--format json` to write reports only.

The direct HERA and LHCb fit uses measured bin integrals and fixed single channel SOFT model absorption, without amplitude banks or event generation. It assumes small photon momentum fractions and an on shell vector meson. The phi energy dependence, excited vector ratios and Upsilon shrinkage retain the assumptions recorded in `fit.json`. TP parameters are not refitted by this command. Native generator validation is a separate icepack run.

HERA physics inputs, quadrature, fit controls and dataset selections are in the `hera` section of [settings.json](settings.json). `--lhcb-fit` uses this file by default. Use `--lhcb-fit <settings-file.json>` to supply another file with the same structure. Measurements are read through the icepack readers from unmodified files under `HEPData/`. Install missing original papers using `python HEPData/download.py --papers 0903.4205v1 1111.2133v2 hep-ex/0205107v1 1505.08139v2`.

## Other commands

```bash
python develop/tools/resonance_card_builder.py --model GP --fuse 22 990 --spinX2 2 --P -1 --C -1 --basis helicity --helicity-row -1 0 1.0 0.0 --format card
python develop/tools/analytic_Jz.py --J 0 --mmax 1 --pole-spin 2 --format cards
python develop/tools/ls_helicity_conversion.py --J 1 --s1 1 --s2 0
python develop/tools/jpsi_gamma_helicity.py --spin 2 --preset pure-e1 --json
python develop/tools/tuneplot.py figs/icetune/<run> -o tmp/tuneplots.pdf
```

The radiative decay tool also accepts explicit normalized `--m2` and `--e3` amplitudes. `tuneplot.py` pairs `hplot__*.pdf` files under `linear/` and `log/` directories and excludes `default/` predictions.

## Package and checks

Shared functions live in the [lib](lib/) Python package, with `common.py`, `push.py` and `soft_exchange.py` for formatting, card updates and SOFT model proton couplings and form factors. The [lib/hera](lib/hera/) subpackage contains the measurement readers, amplitudes, fits, reports and candidate tune writer. Import these modules directly, for example `from develop.tools.lib.hera import model`.

```bash
ruff check develop/tools
pytest -q -s -rP tests/physics/unit/test_repository_tools.py tests/physics/unit/test_hera_pp.py tests/physics/unit/test_jspi_gamma_helicity.py tests/technical/data/test_hera_upsilon.py
```

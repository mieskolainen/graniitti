# Multichannel eikonal analysis

Compare single and double channel eikonal models for $pp$ and $p\bar p$. Run from the repository root:

```bash
bash tests/physics/studies/multichannel/run.sh
```

Plot existing outputs without running the generator:

```bash
STUDY_ACTION=analyze bash tests/physics/studies/multichannel/run.sh
```

The script activates the environment and scans $\sqrt{s}=200,7000,13000,60000$ GeV for $pp$ and $546,1960$ GeV for $p\bar p$. `SQRTS="7000 13000"` selects the $pp$ energies. Matrix tables go to `eikonal/`, figures to `figs/multichannel/pp/` and `ppbar/`, comparing models, channels and energies.

The analysis accepts `--sqrts`, `--models`, `--channels` and `--ordered-channels`. See `python tests/physics/studies/multichannel/analysis/analyze.py --help`.

Elastic scans also produce the elastic comparison plots. Matrix files are selected from completed scan logs under `tmp/multichannel/<model>/<beam>/`.

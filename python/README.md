# core

Python tools for simulation tuning, inference and HEP data analysis, with GRANIITTI and Pandora drivers. Programs live directly under `core`, with reusable modules under `tune`, `numerics`, `stats`, `kinematics`, `analysis`, `inference`, `io` and `plot`. Simulator drivers live under `tune/drivers`.

Run directly from the repository root after environment setup:

```bash
conda activate graniitti
source install/setenv.sh
python -m core.icetune --help
python -m core.iceplot --help
```

Use the same `python -m core.<tool>` form for `iceproxy`, `icescape`, `iceviz` and `iceweb`. Run the optimizer benchmark with `python -m core.tune.optimizers.icebo.cli`. `install/setenv.sh` adds `python/src` to `PYTHONPATH`. Python dependencies are listed in the root `requirements.txt`. Simulator binaries and HepMC3 bindings are supplied separately. Repository campaigns and CERN submission setup live under `submit/`. Pass an explicit JSON tuning definition, optionally `file.json#/pointer`, to `python -m core.icetune --tunesetup`.

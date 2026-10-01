# Gaussian process uncertainties

Study Gaussian process regression with input and target uncertainties. Run from the repository root:

```bash
conda activate graniitti
source install/setenv.sh
bash tests/physics/studies/inference_gproc/run.sh
```

`run.sh` runs `gaussian_process.py` and saves figures under the repository root `figs/`. It is included in `run_all.sh`. There is no separate analysis stage, so `STUDY_ACTION=analyze` skips this study.

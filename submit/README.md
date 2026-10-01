# icetune submission

Fit definitions and computing resources are set in [campaigns.yml](campaigns.yml). Physics parameters, fit ranges and datasets are selected by the [GRANIITTI](tunecards/graniitti/) and [Pandora](tunecards/pandora/) tune cards. Run all commands from the repository root.

## Step 1. Submit the fit

### Option 1: Single machine

```bash
conda activate graniitti
source install/setenv.sh
python -m submit --campaign tune-gpom-res-con --scheduler local --run-name GPOM_RES_CON
```

### Option 2: CERN lxbatch

The submission setup defaults to `graniitti`. The batch environment is selected by `runtime.conda_env` in `campaigns.yml`.

```bash
module load lxbatch/eossubmit # if using EOS
export CONDA_EXE="/eos/user/<PATH>/miniconda3/bin/conda"
source install/setconda_lxplus.sh
python -m submit --campaign tune-gpom-res-con --scheduler lxplus --run-name GPOM_RES_CON
```

### Option 3: Custom computing cluster

<details>
<summary>Start Ray and submit the fit</summary>

Use the same Conda environment and shared repository path on every node. Replace `<head-ip>` with the head node address. Start the head without trial CPUs:

```bash
ray start --head --node-ip-address="<head-ip>" --port=6379 --num-cpus=0 --num-gpus=0
```

For two workers with 16 CPUs each, run this on each worker:

```bash
ray start --address="<head-ip>:6379" --num-cpus=16 --num-gpus=0
```

Submit once from the head:

```bash
python -m submit --campaign tune-gpom-res-con --run-name GPOM_RES_CON \
  --scheduler external --address "<head-ip>:6379" \
  --workers 2 --worker-cpu 16 --worker-gpu 0 --set ray.cpu_per_trial=4
```

Match resources to the cluster. This example provides CPU capacity for up to eight concurrent trials.

</details>

Replace the campaign name to select another fit, for example `tune-gpom-star-cms-ampfit`. Add `--no-submit` to check the fit and prepare submission files without starting jobs. Other schedulers are `condor`, `slurm` and `pbs`. See `python -m submit --help` for options.

Outputs default to `runs/icetune/` under the repository. Shared paths must be accessible from all computing nodes. Use the resulting `<RUN>` and `<FINGERPRINT>` in Steps 2 and 3.

## Step 2. Fit the likelihood surface

Use icescape to fit a surrogate of the likelihood surface for further minimization and parameter uncertainty estimates:

```bash
python -m core.icescape --input "runs/icetune/<RUN>/campaigns/<FINGERPRINT>/history.json"
```

## Step 3. Fit observable distributions

Use iceproxy to predict histogram bins and MC uncertainties as functions of model parameters, enabling further fits without generating new events. Set `fit.full_output: 1` in the campaign before Step 1 to retain the required trial histograms.

```bash
python -m core.iceproxy --input "runs/icetune/<RUN>/campaigns/<FINGERPRINT>"
```

## Step 4. Update model parameters

Apply the selected fit result to the model cards:

```bash
python -m core.icetune --push figs/icetune/<PATH>/summary.json modeldata/TUNE0
```

## Event samples for simulation based inference

Save complete MC samples for surrogate training with `fit.save_events: 1` in the campaign, or:

```bash
python -m submit --campaign tune-gpom-res-con --scheduler lxplus \
  --run-name GENTRAIN --set fit.save_events=1
```

Successful trials retain HepMC3 events, model cards, generator settings and random seeds. Samples are transferred directly from workers to shared storage, verified by checksums and retried on transfer failure. They are stored separately from fit checkpoints.

For amortized Bayesian inference requiring a known parameter sampling density, use uniform random trials in Step 1. Adaptive proposals generally have no simple analytic density and complicate the required correction. See [AIMS25](https://github.com/mieskolainen/AIMS25) for generative model training.

## Settings

Set these fields in `campaigns.yml`, or override them with `--set SECTION.FIELD=VALUE`.

| Field | Purpose |
| --- | --- |
| `runtime.shared_output_dir` | Shared output directory, defaulting to `--repo-dir` or the current directory |
| `runtime.conda_env` | Conda environment name or absolute path |
| `runtime.pfa_dir` | Pandora installation, default `../PandoraPFA` |
| `fit.full_output` | Set to `1` to save trial histograms for iceproxy |
| `fit.save_events` | Set to `1` to retain complete event samples from successful GRANIITTI trials |
| `ray.init_jobs` | Maximum concurrent sample and amplitude bank preparation jobs |
| `ray.init_cpu` | CPUs per preparation job |
| `ray.workers` | Worker allocations for the fit |

Relative storage and Pandora paths are resolved from `--repo-dir`. On lxplus, storage and the Conda installation must be accessible from worker nodes. Conda is found through `CONDA_EXE` or `PATH`, and jobs receive the resolved absolute environment path. Set `CONDA_EXE=/path/to/bin/conda` if needed.

## Amplitude fits on lxplus

For each active sample, SAMPLE generates events and BANK constructs amplitude banks. With `ray.init_jobs > 1`, banks are built in independent jobs and different samples can be prepared concurrently. INIT combines the banks before STEER starts the fit. The head proposes fit parameters, while workers evaluate amplitudes and gradients. The head does not run trial jobs.

Prepared code and amplitude banks are stored under the run's `ray/runtime/` directory and copied to local worker storage.

See [icetune](../docs/icetune.md) for fit settings and [the main README](../README.md) for output locations and supported optimizers.

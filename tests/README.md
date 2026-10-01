# Tests

## Local tests

Run from the repository root:

```bash
conda activate graniitti
source install/setenv.sh
cmake -S . -B build -DWITH_TEST=ON
cmake --build build -j4
ctest --test-dir build --output-on-failure
pytest -q -s -rP tests/technical tests/physics/unit --run-integration
```

## Condor

Submit the full suite or select one:

```bash
bash tests/condor/submit.sh all      # Full build and test suite
bash tests/condor/submit.sh cpp      # C++ tests
bash tests/condor/submit.sh pytest   # Python tests
bash tests/condor/submit.sh studies  # Physics studies
```

The program path and conda environment must be visible to workers. Builds use eight compiler jobs by default. Resource settings are in [condor/SETTINGS.json](condor/SETTINGS.json). Results are saved under `runs/tests/`.

## icepacks

Submit all active measurement comparisons:

```bash
NEVENTS=50000 LOOPSCREEN=1 WEIGHTED=1 bash tests/condor/submit.sh icepacks
```

Submit pure MC studies, for example the SPIN family:

```bash
NEVENTS=50000 LOOPSCREEN=1 WEIGHTED=1 bash tests/condor/submit.sh icepacks --name SPIN
```

Each job generates events and runs iceplot. Results are saved under `figs/icepack/`. Open the run's `report.md` for comparison results. See [icepack](../icepack/README.md) for details.

## Test directories

| Directory | Contents |
| --- | --- |
| [cpp/](cpp/) | C++ library and amplitude tests |
| [technical/](technical/) | Python, integration and basic operation tests |
| [physics/unit/](physics/unit/) | Physics functions, observables and cuts |
| [physics/validation/](physics/validation/) | Published measurement comparisons |
| [physics/symbolic/](physics/symbolic/) | Symbolic calculations |
| [physics/studies/](physics/studies/README.md) | Generation and analysis studies |
| [external/](external/README.md) | Pythia and Rivet comparisons |
| [references/](references/) | Published reference tables and generator cards |
| [data/](data/) | Test input files |
| [condor/](condor/) | Submission scripts and resource settings |

# External tests

After the [common setup](../README.md), run the Pythia and Rivet workflows from the repository root with their external environments installed:

```bash
NEVENTS=10000 LOOPSCREEN=0 WEIGHTED=1 bash tests/external/pythia/cases/hard_pomeron_z/run.sh
NEVENTS=10000 bash tests/external/rivet/generic/run.sh
NEVENTS=10000 bash tests/external/rivet/cms/run.sh
NEVENTS=10000 bash tests/external/rivet/star/run.sh
```

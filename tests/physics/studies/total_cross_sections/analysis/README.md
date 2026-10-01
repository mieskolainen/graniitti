# Minimum bias cross sections

Generate, fit and plot total, elastic and diffractive cross sections from the repository root:

```bash
conda activate graniitti
source install/setenv.sh
bash tests/physics/studies/total_cross_sections/run.sh
```

CSV tables with integration errors for all four processes are saved in this analysis directory, with plots under its `figs/`. Regenerate tables lacking error columns. `ENERGIES=7000` selects a single collision energy in GeV. Fits use only measurements with matching beams and energies within the scan, with logarithmic interpolation. `scan.json` records the coupling and beam used for generation. If absent, supply `--g3p-orig` and `--beam`.

Comparison results are saved in `figs/measurements.json` before reporting failure for disagreement or insufficient MC precision.

Repeat the analysis with `python tests/physics/studies/total_cross_sections/analysis/analyze.py --quiet-fit`. Use `--scan` for another input table and `--help` for options.

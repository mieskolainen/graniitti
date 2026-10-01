# Sudakov and Shuvaev plots

Run a Durham process to generate interpolation tables under `sudakov/`, then plot the selected JSON files:

```bash
python tests/physics/studies/sudakov/analyze.py path/to/SUDA_file.json path/to/SHUV_file.json
```

The JSON files store compressed double precision arrays, as for eikonal and nuclear tables. Coordinates, logarithmic scales, collision energy, PDF, grid sizes, ranges and curve positions are read from these files. Each input has a separate figure directory under this study's `figs/`. `--fig-dir` changes the location and `--slices` sets the number of curves. ASCII tables must be regenerated.

# Elastic energy scans

Generate and plot a scan from the repository root:

```bash
conda activate graniitti
source install/setenv.sh
EIKONAL_MODEL=double bash tests/physics/studies/elastic/run.sh
```

`EIKONAL_MODEL` selects `single`, `double` or `triple`. `STUDY_ACTION=analyze` plots completed scans from `tmp/elastic/<model>/`. Figures are saved below this analysis directory in `figs/<beam>/N1`, `N2` or `N3`. Momentum space plots include measurements selected by `icepack/ELASTIC` cards.

See `python -m tests.physics.studies.elastic.analysis.analyze_momentum --help`, or `analyze_impact`, for beam, energy and channel selections. Channel `0,0` is physical elastic scattering.

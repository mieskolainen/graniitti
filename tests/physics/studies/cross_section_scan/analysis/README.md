# Central production cross sections

Compare central production with intact protons, single dissociation and double dissociation with and without Pomeron loop screening. After the [common setup](../../../../README.md), run from the repository root:

```bash
bash tests/physics/studies/cross_section_scan/run.sh
```

`run.sh` generates and plots both scans. `STUDY_ACTION=analyze` reuses saved tables. `ENERGIES=7000,13000` selects collision energies in GeV. Tables are saved in this analysis directory, with plots in `figs/xsec.pdf` and `figs/xsecratios.pdf` beneath it. The analysis accepts other tables through `--false-scan` and `--true-scan`. See `--help` for plotting options.

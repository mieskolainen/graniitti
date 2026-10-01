# Physics icepacks

Run icepacks and check results with the [test commands](../tests/README.md). Dataset and generator cards define process parameters and random seeds. Use `bash icepack/run.sh --help` for generation and plotting options. For Pythia comparisons, see [external tests](../tests/external/README.md).

## Example:

```bash
bash icepack/run.sh PHOTOPROD/LHCb* --NEVENTS 10000 --LOOPSCREEN 1 --WEIGHTED 1
```

Outputs go to `figs/icepack/<entry>/`. For example, `DURHAM/partons` writes to `figs/icepack/DURHAM__partons/`:

| Directory | Contents |
| --- | --- |
| `generation/<run>/<sample>/` | Generator logs, cards and commands to repeat generation |
| `logs/` | Complete run logs and commands to repeat iceplot |
| `reports/` | Analysis reports |
| `plots/` | Plots, including density and stacked variants |

Use `--output-dir PATH` to change the output root or `--report-dir PATH` for reports only. Events stay in `output/` and integration grids in `vgrid/`.

**Automated comparisons.** Comparisons check distribution shapes, integrated cross sections and MC statistical precision using the [test criteria](../tests/README.md) and [_common/SETTINGS.json](_common/SETTINGS.json) thresholds. Readers process original versioned HEPData JSON tables without modifying them. Histogram `covariance` fields select HEPData covariance tables. When HEPData omits bin correlations, document the assumed correlations in the reader. Measurement comparisons without original HEPData input are inactive. Cited numerical results and theory curves are stored in [tests/references](../tests/references).

**MC normalization.** `mc_scale` multiplies the MC prediction to match the measured quantity, including any decay branching fraction or photon flux conversion. Document its formula, PDG branching fraction and citations in card comments. The ZEUS 582237 and ALICE 2658375 cards include published photon fluxes in `mc_scale`. Readers apply fluxes tabulated in HEPData.

**UPC neutron classes.** Analysis cuts select final HepMC3 beam remnant neutrons. Generate all neutron multiplicities and use internal or external nuclear fragmentation to write the final remnant particles. Combined `0nXn` selections accept zero neutrons from one beam and at least one from the other, in either order. If a remnant is stored as a placeholder or an intermediate nucleus without decay products, its neutron multiplicity is unknown and neutron classification raises an error.

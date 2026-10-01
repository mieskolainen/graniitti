# Rivet tests

Install Rivet with the three analyses used here:

```bash
RIVET_ANALYSES="CMS_2011_I954992 STAR_2020_I1792394 MC_FSPARTICLES" bash install/installRIVET.sh
```

Omit `RIVET_ANALYSES` to build all bundled analyses. The installer uses HepMC3 and at most four compiler jobs, with build files and logs under `tmp/rivet-install/`. After the [common setup](../../README.md), select `all` (default), `generate` or `analyze`:

```bash
NEVENTS=10000 bash tests/external/rivet/cms/run.sh all
bash tests/external/rivet/cms/run.sh analyze
NEVENTS=10000 bash tests/external/rivet/star/run.sh all
NEVENTS=10000 bash tests/external/rivet/generic/run.sh all
```

For weighted generation with a cross section estimate, first save a VEGAS proposal, then generate with an independent seed:

```bash
NEVENTS=-1 CORES=2 WEIGHTED=1 bash tests/external/rivet/cms/run.sh generate
NEVENTS=10000 CORES=2 WEIGHTED=1 SEED=67890 VGRID=vgrid/rivet_cms.vgrid bash tests/external/rivet/cms/run.sh all
```

Replace `cms` with `star` for STAR. Rivet uses event weights and the HepMC3 cross section. GRANIITTI histograms are disabled. Generation alone requires only GRANIITTI. For analysis, Rivet must be on `PATH` or configured through `RIVET_ENV=/path/to/rivetenv.sh`, defaulting to `.local/rivet/rivetenv.sh` when present. The `all` action checks Rivet before generation. Rivet and generator environments are loaded separately.

Events go to `output/rivet_CASE.hepmc3`, histograms and plots to `runs/tests/rivet/rivet_CASE/`, relative to the automatically located repository root. `TAG` selects the sample for generation and analysis. `NEVENTS`, `SEED` and `CORES` control generation. `VGRID` selects a compatible saved VEGAS proposal. Screening and event weighting follow the card unless overridden by `LOOPSCREEN=0|1` or `WEIGHTED=0|1`.

## Fiducial cuts

The CMS card describes exclusive dimuons at 7 TeV with $p_T>4$ GeV, $|\eta|<2.1$ and $m_{\mu\mu}>11.5$ GeV, as in [CMS, arXiv:1111.5536](https://arxiv.org/abs/1111.5536). The CMS Rivet analysis further requires exactly two charged particles within $|\eta|<2.4$, a muon opening angle below $0.95\pi$, $\Delta\phi>0.9\pi$ and $|\Delta p_T|<1$ GeV. Its accepted cross section can therefore differ from the generator integral.

The STAR card describes pion pairs at 200 GeV with $p_T>0.2$ GeV, $|\eta|<0.7$ and the measured proton momentum acceptance through `USERCUTS=1792394000`, as in [STAR, arXiv:2004.11078](https://arxiv.org/abs/2004.11078). Compare the pion distributions within these cuts. The sample contains no kaons or protons, although the Rivet analysis includes their histograms. `MC_FSPARTICLES` provides general event distributions.

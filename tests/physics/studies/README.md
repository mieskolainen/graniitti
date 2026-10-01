# Studies

| Study | Description | Run or analyze |
| --- | --- | --- |
| [alice_harmonic](alice_harmonic/) | Spherical harmonic reconstruction with ALICE cuts | [run.sh](alice_harmonic/run.sh) |
| [alice_multi](alice_multi/) | Pion, kaon and proton pair spectra with ALICE cuts | [run.sh](alice_multi/run.sh) |
| [alice_single](alice_single/) | Pion pair comparison with external ALICE data | [run.sh](alice_single/run.sh) |
| [atlas_multi](atlas_multi/) | Pion pairs with ATLAS proton tagging and transverse momentum selections | [run.sh](atlas_multi/run.sh) |
| [cdf_single](cdf_single/) | Pion pair production with CDF cuts | [run.sh](cdf_single/run.sh) |
| [cms_harmonic](cms_harmonic/) | Spherical harmonic reconstruction with CMS central and proton tagging cuts | [run.sh](cms_harmonic/run.sh) |
| [cms_multi](cms_multi/) | Pion, kaon and proton pair spectra with CMS cuts | [run.sh](cms_multi/run.sh) |
| [cms_mumu](cms_mumu/) | EPA and QED dimuon predictions at low and high masses | [run.sh](cms_mumu/run.sh) |
| [cross_section_scan](cross_section_scan/analysis/README.md) | Central production energy scans with screening and proton dissociation | [run.sh](cross_section_scan/run.sh), [analyze.py](cross_section_scan/analysis/analyze.py) |
| [decay_chains](decay_chains/) | Resonance decay chains with stable and decaying daughters | [run.sh](decay_chains/run.sh) |
| [durham_chic0](durham_chic0/) | Durham $\chi_{c0}\to\pi^+\pi^-$ production with different PDFs | [run.sh](durham_chic0/run.sh) |
| [durham_mmbar](durham_mmbar/) | Durham continuum production of pion, kaon, eta and eta prime pairs | [run.sh](durham_mmbar/run.sh) |
| [elastic](elastic/analysis/README.md) | Elastic momentum and impact parameter space amplitudes across beam energies | [run.sh](elastic/run.sh) |
| [excitation](excitation/) | Pion pair production with zero, one or two dissociating protons | [run.sh](excitation/run.sh) |
| [flux](flux/) | Photon EPA and Durham gluon flux comparisons | [run.sh](flux/run.sh) |
| [inference_ampfit](inference_ampfit/) | GPOM amplitude refinement, interference analysis and independent STAR/CMS simulation checks | [refine.py](inference_ampfit/refine.py), [inspect.py](inference_ampfit/inspect.py), [validate.py](inference_ampfit/validate.py) |
| [inference_gproc](inference_gproc/README.md) | Gaussian process regression with input and target uncertainties | [run.sh](inference_gproc/run.sh) |
| [jw_frames](jw_frames/) | Tensor resonance decay angles in different polarization frames | [run.sh](jw_frames/run.sh) |
| [jw_polarization](jw_polarization/) | Jacob-Wick decay distributions for different spin projections | [run.sh](jw_polarization/run.sh) |
| [lux](lux/) | Inclusive photon fusion lepton pair production with LUX photon PDFs | [run.sh](lux/run.sh) |
| [minbias](minbias/) | Minimum bias event generation and Rivet comparisons | [run.sh](minbias/run.sh) |
| [multichannel](multichannel/analysis/README.md) | Eikonal model, channel and beam energy comparisons | [run.sh](multichannel/run.sh) |
| [odderon_photoprod](odderon_photoprod/) | Photon and Odderon production of $\phi(1020)$ with proton tagging selections | [run.sh](odderon_photoprod/run.sh) |
| [optimal_transport](optimal_transport/) | Optimal transport comparison of continuum and resonance pion pair samples | [run.sh](optimal_transport/run.sh) |
| [phase_space](phase_space/) | Two, four and six pion phase space constructions | [run.sh](phase_space/run.sh) |
| [pileup_kron](pileup_kron/README.md) | Delphes jet reconstruction with pileup, PUPPI and graph Kron comparisons | [run.sh](pileup_kron/run.sh) |
| [rhorho_phiphi](rhorho_phiphi/) | $\rho\rho$ and $\phi\phi$ production and decays in MP and TP models | [run.sh](rhorho_phiphi/run.sh) |
| [screening](screening/) | Kaon pair continuum with and without Pomeron loop screening | [run.sh](screening/run.sh) |
| [sommerfeld](sommerfeld/) | Sommerfeld diffraction fields and Feynman path integral distributions | [plotdiffraction.py](sommerfeld/plotdiffraction.py), [plotpath.py](sommerfeld/plotpath.py) |
| [sudakov](sudakov/README.md) | Sudakov and Shuvaev interpolation plots from generated JSON caches | [analyze.py](sudakov/analyze.py) |
| [tensor0](tensor0/) | Tensor Pomeron scalar and pseudoscalar coupling comparisons | [run.sh](tensor0/run.sh) |
| [tensor2](tensor2/) | Seven tensor Pomeron production couplings for $f_2(1270)$ | [run.sh](tensor2/run.sh) |
| [tensor_spectrum](tensor_spectrum/) | Tensor Pomeron pion, kaon and proton pair spectra | [run.sh](tensor_spectrum/run.sh) |
| [total_cross_sections](total_cross_sections/analysis/README.md) | Total, elastic and diffractive cross section scans and measurement comparisons | [run.sh](total_cross_sections/run.sh) |

Run all studies with a `run.sh` script:

```bash
bash tests/physics/studies/run_all.sh
```

Run the selected plotting studies:

```bash
bash tests/physics/studies/run_plots.sh
```

Run one study:

```bash
NEVENTS=10000 WEIGHTED=true LOOPSCREEN=false bash tests/physics/studies/screening/run.sh
```

`STUDY_ACTION=both` generates and analyzes by default. `STUDY_ACTION=analyze` reuses saved results where supported. Screening comparisons and elastic studies set screening explicitly. Each study script or README lists its inputs and outputs.

Enable optional studies with their required inputs and `RUN_HARMONIC=1` for `alice_harmonic` and `cms_harmonic`, `RUN_RIVET=1` for `minbias`, `RUN_DELPHES=1` for `pileup_kron`, and `RUN_EXTERNAL_DATA=1` for `alice_single`. Studies with missing inputs are skipped and reported.

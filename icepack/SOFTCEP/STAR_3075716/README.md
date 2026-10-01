# STAR 510 GeV central exclusive production

Original tables: [HEPData 167048, version 1](https://www.hepdata.net/record/ins3075716?version=1). Measurement: [JHEP 07 (2026) 143](https://doi.org/10.1007/JHEP07(2026)143), [arXiv:2510.27482](https://arxiv.org/abs/2510.27482). The 36 tables in Figures 7 to 13 cover $\pi^+\pi^-$, $K^+K^-$ and $p\bar p$. Run `bash icepack/run.sh SOFTCEP/STAR_3075716` after `conda activate graniitti` and `source install/setenv.sh`.

The generator cards apply $|\eta|<0.9$ and the species dependent lower $p_T$ cuts. `src/MUserCuts.cc` hardcodes the four Roman Pot regions and the pair upper $p_T$ cuts with identifiers `3075716000`, `3075716010` and `3075716020`. West is positive $z$. Equation (4.1) uses the second circle for ED and WU, and the upper $|p_y|$ bound for EU and WD. Python selections apply the published mass and proton azimuth regions.

The reader preserves asymmetric errors and masks unpublished bins. It assumes uncorrelated experimental errors because no bin covariance is supplied. Only uncertainties supplied in the original HEPData JSON are included. These tables do not supply a separate luminosity component. Cross sections retain their absolute normalization. The original HEPData files remain unchanged.

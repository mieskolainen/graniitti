# *DOCS*
## Under work

#### 30/09/2026
mikael.mieskolainen@cern.ch

-----

#### Parameters

- `(GP|MP|XP|TP)[RES+CON]` processes require global multi-energy dataset icetune runs, with current focus on the `GP` model, currently MP/XP with a placeholder values, TP with paper defaults.


#### Processes

- UPC all around, especially forward emissions.

- Dissociation in semi-exclusive CEP (inelastic vertices, fragmentation).

- (Minimum bias processes, in the future ...).

- Understand GP model (resonance) / continued trajectory amplitude behavior with higher $|t_{1(2)}|$ CEP + SD/DD dissociative processes (see below related).


#### Numerics and amplitudes

- GP model LS coupling mode angular singularities: the singularity of the continued scalar fusion coefficient at $|\alpha_P(t_1)-\alpha_P(t_2)|=1$ is outside the current STAR/CMS forward proton fiducial cuts for nominal `TUNE0` parameters. Wider cuts, changed trajectory parameters and screening integration ranges require separate checks. The angular floor `GP_alpha_min=-0.4900` enforces $2\alpha_{\mathrm{ang}}+1\geq0.0200$, avoiding zero denominators of this form while propagators retain the full trajectory.

- Inspect complete vertex products and sums to distinguish removable singularities from those requiring changes to the amplitude or its LS analytic continuation. Simplify factors before numerical evaluation, then examine the remaining continued Wigner functions for singularities.

- Study helicity structure for `GP[CON]`, currently `m=0` used only in `CON_GP.json`. E.g. "tilt" in $|t_1+t_2|$-distribution of data vs GRANIITTI can be `|m| != 0`.

- Screening loop integral discretization, under `NUMERICS.json`, find optimal tradeoff between accuracy and speed (process dependent).

- Durham QCD loop integral discretization, under `NUMERICS.json`, find optimal tradeoff between accuracy and speed (subprocess dependent).


#### Tuning and validation

- Improve icepack driven autotesting.
- Improve condor submission logic for icepacks and tests.
- Improve unit tests.
- Improve automated tuning.

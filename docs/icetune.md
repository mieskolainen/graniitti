# *DOCS*
## icetune [R&D]

#### 30/09/2026
mikael.mieskolainen@cern.ch

Execution and submission: [submit/README.md](../submit/README.md).


## Programs and configuration

| Program | Input | Calculation |
| --- | --- | --- |
| icetune | [Tune card](../submit/tunecards/), [campaign](../submit/campaigns.yml), [optimizer settings](../python/src/core/tune/settings/) | Parameter fit using simulation or amplitude reweighting |
| icescape | `history.json`, including completed trials during optimization | Surrogate of the recorded likelihood, best fit, profiles and parameter uncertainties |
| iceproxy | Saved trial histograms | Histogram and MC uncertainty surrogates with quadratic likelihood. Requires `fit.full_output: 1` before tuning. Gaussian input is unsupported |

| Setting or interface | Meaning |
| --- | --- |
| Direct `core.icetune` | Options use underscores. Campaign settings are not loaded automatically |
| `fit.plot: 1` | Initial and improving best fit comparisons |
| `fit.full_output: 1` | Trial histograms for iceproxy |
| `fit.save_events: 1` | Trial events and cards for generative model training |
| Surrogate `--posterior` | Posterior samples with uniform priors within optimizer bounds |

## Parametrizations

For a complex coupling, $g=R e^{i\phi}$. For a real coupling vector, $\boldsymbol g=R\boldsymbol n$ with $\boldsymbol n^T\boldsymbol n=1$.

| Coupling `mode` | Coordinates | Use |
| --- | --- | --- |
| `none` | None | Fixed coupling |
| `magnitude` | $R$ | Vary magnitude at fixed phase. With vector geometry, also vary the real direction |
| `phase_raw` | $\phi\in[-\pi,\pi]$ | One periodic phase coordinate |
| `phase_cayley` | $u,v\in[-1,1]$, $\phi=2(\arctan u+\arctan v)$ modulo $2\pi$ | Two redundant coordinates for one phase |
| `mag_phase_raw` | $R,\phi$ | Magnitude and periodic phase |
| `mag_phase_cayley` | $R,u,v$ | Magnitude and Cayley phase |
| `mag_phase_cartesian` | $\operatorname{Re}g,\operatorname{Im}g$ | Complex coefficient without an angular coordinate at zero. Component bounds define a rectangle |

| Geometry | Physical identification | When to use |
| --- | --- | --- |
| `direct` | Individual coupling components | Independent rows, including row phases where supported |
| `projective` | $\boldsymbol n\sim-\boldsymbol n$, direction in $\mathbb{RP}^{N-1}$ | Real relative couplings with an irrelevant common sign, or with the sign carried by a separate overall coefficient. Uses $N-1$ angles and removes the redundant norm |
| `spherical` | Oriented direction in $S^{N-1}$, retaining $\boldsymbol n$ and $-\boldsymbol n$ | Real amplitudes whose common sign matters for interference. Uses $N-1$ angles |
| Simplex | $p_i\geq0$, $\sum_i p_i=1$ | Incoherent spin populations. Squared spherical coordinates enforce positivity and normalization |

### Selection in tune cards

Paths below are relative to an entry in `models`. The builders generate the angle coordinates and topology metadata.

| Tune card fields | Effect and restrictions |
| --- | --- |
| `production: {"mode": "magnitude", "geometry": "projective"}` | Real XP/GP production direction and overall norm. TP also supports real tensor directions. XP/GP rows must be real |
| `production.geometry: "spherical"` | Retain the orientation of real production vectors. XP/GP requires `production.mode: "magnitude"` |
| `production.mode: "mag_phase_cartesian"`, `production.geometry: "projective"` | XP/GP real direction multiplied by one complex coefficient. Select the resonances with `phi.mode: "phase_raw"` and optional `phi.resonances`. The builder replaces their separate overall phase with Cartesian coefficient coordinates |
| `phi.mode: "phase_raw"` | Periodic overall resonance phase. `phi` supports only `none` and `phase_raw` |
| `decay.mode: "magnitude"` | Normalized real projective decay vector. No independent complex row phases. `decay.zeta.mode` controls the separate decay phase |
| `production.coherence: "coherent"` | MP coherent spin amplitudes. `production.geometry` selects projective or spherical directions. `direct` selects the projective default for this spin basis |
| `production.coherence: "incoherent"` | MP parity symmetric spin populations on a simplex. Production mode is `none` or `magnitude`, decay amplitudes and phases must be fixed |
| `continuum_options.production_geometry: "projective"` | Real continuum direction with `production_mode` set to `magnitude` or `mag_phase_cartesian`. Direct geometry is also supported, spherical is not |
| `production.active_only`, `continuum_options.production_active_only` | Select initially nonzero rows. Including zero continuum rows requires projective geometry and a nonzero total norm |
| Top level `parameters` entries | Explicit `type`, `lower`, `upper`. Finite bounds, continuous uniform or integer domains. Gradient fits require continuous parameters |

Examples: [GP amplitude fit](../submit/tunecards/graniitti/tune-gpom-star-cms-ampfit.json), [MP coherent amplitude fit](../submit/tunecards/graniitti/tune-mpom-star-cms-ampfit.json).

### Topology in optimization and surrogates

| Representation | Topology features |
| --- | --- |
| Raw phase | Periodic features $(\cos\phi,\sin\phi)$ identify phases differing by $2\pi$. Use raw modes for automatic phase topology |
| Cayley phase | Decodes to the same physical phase, but the two coordinates do not generate an automatic phase topology group |
| Projective direction | Features from $\boldsymbol n\boldsymbol n^T$ identify opposite directions. Spherical features retain $\boldsymbol n$ itself |
| Overall coefficient and direction | Joint features represent the physical complex vector, preserving interference and removing direction dependence at zero magnitude |
| Simplex | Features use the physical populations, not unconstrained angles |

| Optimization method or surrogate | Topology treatment |
| --- | --- |
| HEBO | Topology kernel controlled by `topology.enabled` in [hebo.json](../python/src/core/tune/settings/hebo.json), enabled by default |
| ICEBO | Uses topology metadata for covariance, distances and direction sampling |
| ampfit | Differentiates the decoded physical amplitudes. Requires continuous parameters and precomputation, uses no surrogate kernel |
| icescape | Uses the saved topology metadata automatically |
| iceproxy | Uses the saved topology metadata automatically |

Definitions: [coordinates](../python/src/core/tune/parameters/tools.py), [topology features](../python/src/core/tune/parameters/topology.py), [model builders](../python/src/core/tune/drivers/graniitti/tunesetup/).

## Costs and likelihoods

GRANIITTI options. Let $r=m-d$ and $C=C_{\mathrm{data}}+C_{\mathrm{MC}}$.

| Campaign setting | Direct option | Direct default |
| --- | --- | --- |
| `fit.cost` | `--cost` | `chi2` |
| `fit.cost_avg` | `--cost_avg` | `global-mean` |
| `fit.data_covariance_mode` | `--data_covariance_mode` | `diagonal` |
| Direct option only | `--cost_rho` | `quadratic` |

### Cost and residuals

| `cost` | Definition and use | Restrictions |
| --- | --- | --- |
| `chi2` | $r^T C^+r$, with observable fit weights. Additive deviations relative to data and MC uncertainties | Always quadratic. Only `dataset-mean` changes the unscaled weighted sum |
| `ratio2` | Sum of $\rho(z_i)$, where $z_i$ is the difference $\operatorname{asinh}(m_i/s_i)-\operatorname{asinh}(d_i/s_i)$ divided by its propagated uncertainty, with $s_i^2=\sigma_{m,i}^2+\sigma_{d,i}^2$ | Approximately logarithmic at large positive values, defined also at zero and for negative bins. Correlations propagate through the asinh derivatives |
| `wasserstein` | Integral of the absolute cumulative cross section difference, averaged over uncertainty replicas. Histograms are not normalized to unit area | Nonnegative bins. No `cost_rho` dependence. Not supported by ampfit |
| `gaussian` | $Q=r^T C^+r+\log\operatorname{pdet} C$, Gaussian $-2\log L$ including covariance normalization | Requires ampfit, `sum`, `quadratic` and unit positive fit weights |

| `cost_rho` | $\rho(z)$ | Behaviour |
| --- | --- | --- |
| `quadratic` | $z^2$ | Standard least squares |
| `absolute` | $\lvert z\rvert$ | Linear penalty |
| `cauchy` | $\log(1+z^2)$ | Logarithmic penalty at large residuals |
| `huber` | $z^2/2$ for $\lvert z\rvert\leq1$, otherwise $\lvert z\rvert-1/2$ | Quadratic near zero, linear beyond one standard deviation |
| Full covariance | Penalty on uncorrelated covariance modes | ampfit requires `quadratic`, nonquadratic modes fail during evaluation |

### Averaging

$H$ and $D$ count active histograms and datasets. $c_h$ is the histogram cost, $n_h$ its informative bin count or retained covariance contribution, and $\bar w_h=Hw_h/\sum_h w_h$.

| `cost_avg` | Combined cost | Emphasis |
| --- | --- | --- |
| `sum` | $\sum_h w_h c_h$ | Unscaled weighted sum |
| `global-mean` | $\sum_h\bar w_h c_h/\sum_h n_h$ | Per bin |
| `local-mean` | $H^{-1}\sum_h\bar w_h c_h/n_h$ | Per histogram |
| `dataset-mean` | $D^{-1}\sum_a (\sum_{h\in a}w_h c_h)/(\sum_{h\in a}w_h n_h)$ | Equal dataset weight, weighted mean per bin within each dataset |
| `global-mean-unweighted` | $\sum_h c_h/\sum_h n_h$ | Per bin, ignores positive weight magnitudes |
| `local-mean-unweighted` | $H^{-1}\sum_h c_h/n_h$ | Per histogram, ignores positive weight magnitudes |

| Weight or covariance convention | Implementation |
| --- | --- |
| `fitw` | $w_h>0$, valid histogram and $n_h>0$ required. Fitted parameters are not subtracted from $n_h$ |
| Full covariance contributions | Histogram costs partition the joint cost, including correlations with other histograms |
| Full covariance `ratio2`, `sum` | Weights normalized to unit histogram mean. Differs from diagonal `sum` for nonunit mean weight |
| Full covariance `ratio2`, robust residuals | Weights enter residuals before the penalty. Diagonal mode weights the resulting penalty. Nonunit weights can therefore give different costs even without correlations |

### Covariance and constructed likelihood

| Setting | Choice | Meaning |
| --- | --- | --- |
| `fit.data_covariance_mode` | `diagonal` | Measurement and MC variances per bin |
| `fit.data_covariance_mode` | `full` | Measurement and MC correlations. Requires precomputation. Wasserstein uses data covariance within each histogram but only diagonal MC fluctuations |
| `optimizer.ampfit.bank.mc_stat` | `source` | Fixed MC uncertainties from the source sample |
| `optimizer.ampfit.bank.mc_stat` | `reweighted` | MC uncertainties recomputed under amplitude reweighting |

$C^+$ is the pseudoinverse on the retained covariance subspace, including normalization constraints. The pseudodeterminant $\operatorname{pdet} C$ is the product of the retained covariance eigenvalues. It is constant for fixed covariance and contributes to the Gaussian fit, gradients and Hessian otherwise. Parameter independent constants are omitted.

| Selected cost or quantity | Likelihood convention |
| --- | --- |
| `chi2`, `ratio2`, `wasserstein` | Record the unscaled quadratic $\chi^2$ separately from the ranking cost, with $Q\simeq\chi^2$ omitting the determinant |
| `gaussian` | Record $Q=\chi^2+\log\operatorname{pdet} C$. Per histogram diagnostics remain chi square. MC uncertainty on $Q$ is not estimated |
| `two_nll`, `nll`, `logL` | $Q$, $Q/2$, $-Q/2$, respectively. Representations, not additional cost choices |
| Nonunit `fitw`, averaged costs, robust penalties | Change the statistical interpretation. Ranking costs cannot be interpreted directly as $-2\log L$ |

Implementation: [objectives](../python/src/core/stats/objective.py), [averaging and likelihood](../python/src/core/tune/likelihood.py).

## Optimization and validation

| Method or check | Treatment |
| --- | --- |
| HEBO, ICEBO | Bayesian optimization |
| ampfit | Complex amplitude reweighting with bounded L-BFGS and automatic differentiation. Retains interference, generates no events per gradient step. Selects the best completed point across independent starts |
| Trial execution | Workers evaluate trials. The head optimizes and collects results, never runs trials |
| Random seed | Reproduces starting points. Asynchronous completion can change the optimization trajectory |
| ampfit closure | Initialization checks supported parameters and closure, including reused banks. Check interpolation with `core.tune.drivers.graniitti.ampfit.ampcheck` and independent simulation |
| Parameter intervals | Use the recorded unscaled objective. Finite MC uncertainty remains a Gaussian approximation. Check surrogate coverage and errors on withheld trials |
| Cross section validation | After `xsmode: sample` scans, repeat VEGAS integration with `xsmode: reset`, independent events and identical fiducial cuts |

## Results

| Output or setting | Contents |
| --- | --- |
| `runs/icetune/RUN/`, `figs/icetune/RUN/` | Run outputs. lxplus adds `campaigns/FINGERPRINT/` beneath each |
| `history.json`, `status.json` | Trials, timing and status |
| `summary.json`, comparison PDFs | Selected fit and data comparisons |
| `cost_evolution*.png`, `cost_sorted*.png`, `parameter_evolution_*.pdf` | Cost and parameter histories, updated during optimization |
| `results/`, `failures/` | Optional histograms, covariance and failure reports |
| Best fit plots | Use saved predictions |
| Model card export | Preview and confirmation through [Step 4](../submit/README.md#step-4-update-model-parameters), followed by [independent validation](../tests/physics/studies/inference_ampfit/validate.py) |

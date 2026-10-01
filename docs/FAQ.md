# *DOCS*
## GRANIITTI MC Event Generator
**Frequently Asked Questions (FAQ)**

#### 30/09/2026
mikael.mieskolainen@cern.ch

## CONTENTS

- [Quick start](#quick-start)
- [Compile](#compile)
- [First run](#first-run)
- [Source](#source)
- [Programs](#programs)
- [Automated tests](#automated-tests)
- [Installation](#installation)
- [Basic operation](#basic-operation)
- [Monte Carlo techniques](#monte-carlo-techniques)
- [Analysis](#analysis)
- [Advanced](#advanced)
- [Other](#other)

## QUICK START

See [README.md](../README.md) for installation and first run instructions on Ubuntu Linux.

## COMPILE

See [README.md](../README.md)


## FIRST RUN

Run `./bin/gr` and read the output.


## SOURCE

```
├── bin/                - Compiled binary program files
├── develop/            - Developer scripts
├── docs/               - Documentation
├── eikonal/            - Eikonal screening amplitude and proton density cache arrays
├── figs/               - Figure output directory
├── gencard/            - Steering cards for event generation
├── icepack/            - Measurement data and MC comparison encapsulations
├── HEPData/            - Measurement data
|   ├── CD/             - Central Diffraction
|   └── EL/             - Elastic scattering
├── include/            - Header files (C++)
|   └── Graniitti/
|      ├── Amplitude/MG5/ - MadGraph amplitude headers
|      ├── Analysis/    - ROOT dependent analysis header files
|      └── Math/, Kinematics/, Spin/, Particle/, Eikonal/, PDF/, Regge/, Photon/, QCD/, Process/, Sampling/, Tech/
├── install/            - Installation scripts and packages
├── libs/               - External libraries (C++)
├── MG5cards/           - MadGraph amplitude parameter cards
├── modeldata/          - Global model parameters
├── build/              - CMake build tree
├── output/             - Simulation output folder
├── python/             - Python driven tools such as iceplot, HPC distributed tuning via icetune
├── src/                - Source files (C++)
|   ├── Amplitude/MG5/  - MadGraph amplitudes
|   ├── Analysis/       - ROOT dependent analysis source files
|   ├── Math/, Kinematics/, Spin/, Particle/, Eikonal/, PDF/, Regge/, Photon/, QCD/, Process/, Sampling/, Tech/
|   └── Program/        - Main program source files
|      └── Analysis/    - ROOT dependent main program source files
├── sudakov/            - Special transformed gluon pdf cache arrays
├── nuclear/            - Nuclear (UPC) process cache arrays
├── tests/              - Test suite: automated tests, simulations & analysis benchmark examples
├── runs/               - Parameter inference runs with icetune, test outputs
├── vgrid/              - MC integration grid outputs, useful for e.g. distributed HPC computation
```

## PROGRAMS

### Simulation

| Executable | Description |
|---|---|
| `./bin/gr` | Main simulation program |
| `./bin/xscan` | Process cross section energy dependence scanner |
| `./bin/minbias` | Minimum bias processes (SD,DD,ND) simulation (for minimalistic rapidity gap fluctuation studies) |
| `./bin/hepmc3tolhe` | Convert `.hepmc3` to `.lhe` event output format |

### Analysis and Simulation Based Inference

#### New series (Python)

| Program | Description |
|---|---|
| `./python/*` | Analysis and parameter inference (tuning) with Python 3 libraries |
| `python -m core.iceplot` | Fast histogramming and plotting of fiducial observables against HEPData |
| `python -m core.icetune` | Distributed MC model parameter inference via Bayesian optimization, Active Learning, AI etc. |
| `python -m core.icescape` | Surrogate likelihood function post-fitting based on the icetune trials |
| `python -m core.iceproxy` | Surrogate histogram predictor fitting |

#### Older series (C++)

| Program | Description |
|---|---|
| `./bin/fitharmonic` | Spherical harmonic inverse expansion and analysis |
| `./bin/analyze` | Fast histogramming and plotting of fiducial observables |

#### Other

```
./bin/data2hepmc3    Push events from external sources to .hepmc3 format for a quick analysis (code template)
```

## Automated tests


From the repository root, compile with `cmake -S . -B build -DWITH_TEST=ON && cmake --build build -j4`

```bash
ctest --test-dir build --output-on-failure
```

Run the dedicated physics studies with:

```bash
bash tests/physics/studies/run_all.sh
```

To loop over predefined fiducial measurements and produce a comparison table:

```bash
pytest -q -s -rP tests/physics/validation/test_icepacks.py::test_integrated_xs_table --run-physics --LOOPSCREEN 1 --WEIGHTED 1 --NEVENTS 10000
```



## INSTALLATION

**Q**:  It does not compile :(

*A*:  If you use macOS, try to obtain the latest GCC toolchain. The software is always tested and developed on the latest Ubuntu Linux (+ ROOT 6).


**Q**:  I use lxplus.cern.ch. What is the right environment setup?

*A*:  Try this before compiling anything:

```bash
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_110 x86_64-el9-gcc13-opt
```

Try the Conda setup described in README.md, because that allows you to install the necessary Python libraries.

**Q**:  I use a standalone system, for example Ubuntu Linux. How do I set up the ROOT libraries?

*A*:  If you have compiled ROOT or extracted precompiled binaries from <https://root.cern.ch> to a directory such as `$HOME/local/root`, add these three lines to your environment setup:

```bash
export ROOTSYS=$HOME/local/root
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
```

And make sure that these paths make sense, so that you find binaries/objects/headers under `$ROOTSYS/bin`, `$ROOTSYS/lib`, `$ROOTSYS/include`.

> [!TIP]
> Also try the default ROOT environment setup macro
>
> ```bash
> source $ROOTSYS/bin/thisroot.sh
> ```

> [!TIP]
> Use a ROOT build compatible with the compiler environment and installed Python when enabling ROOT-dependent analysis tools.

**Q**: Compilation with ROOT fails with missing C++20 features such as `std::span` or `requires`.

*A*: GRANIITTI requires C++20 for every target, including ROOT analysis targets. Reconfigure the build with `cd build && cmake .. && make -j4`. A stale cache created with `-DCXX_STANDARD=17` must be reconfigured because current CMake rejects standards below C++20. ROOT may itself have been built as C++17 because its compile flags are not imported into GRANIITTI targets.

**Q**: How do I turn off compilation of tools that depend on ROOT completely?

*A*: Compile with `cd build && cmake -DWITH_ROOT=OFF .. && make -j4`. The simulation part of GRANIITTI is not dependent on ROOT, but some analysis tools depend on it.


## BASIC OPERATION

**Q**: Some processes take a long time to generate.

*A*: The code is maximally multithreaded and most code is reasonably optimized. However, many processes contain very heavy spin algebra. Try decreasing, for example, the number of integrand calls in the VEGAS adaptation stage `ncalls`, the number of adaptation `rounds`, the minimum number of integrand samples in the integration stage `min_samples`, or relaxing the MC integration uncertainty target `precision` and the max weight search quality criteria, which precede event generation. However, the price to pay can be e.g. maximum weight overflows, which technically can invalidate unweighted event generation. Typically, it is better to generate weighted events (and do weighted histogramming) when first exploring the processes, then do final unweighted runs.

> Also try `-DMARCH=native` when compiling, especially with new CPUs and new g++ compilers, to enable the latest instruction sets.

**Q**:  Your generator does not match (my) data.

*A*:  GRANIITTI includes a spectrum of different amplitudes and model variations. Please specify which ones you are using. Every effort is made to match the available data, but, for example, low mass central production processes require parameter optimization on HPC systems, with hundreds of CPU cores running for days, given all the model variations available in GRANIITTI (and not all models are tuned yet).

Try turning the screening (absorption) loop, `LOOPSCREEN`, on and off. For example, some UPC processes or hard interaction processes in $pp$ basically need it always on (very big impact).

First check the process generated and the generation and fiducial cuts used. For strong interaction (Pomeron driven) amplitudes, you may start with the simplest one, `GP[CON]`, and try changing the intermediate form factor types in `CON_GP.json` -- large effects in e.g. $\Delta\phi_{pp}$ and other forward observables. Once this is understood, try adding resonances if those are relevant, `GP[RES+CON]`, and perhaps adjust the resonance spin polarization parameters and so on.

> [!TIP]
> Solid fiducial measurements should be executed in a way that minimizes arbitrary generator model dependence, e.g. by using techniques such as *DeepEfficiency* [https://arxiv.org/abs/1809.06101].

**Q**: What do the ALICE coherent rho neutron classes measure, and how are nuclear fluctuations selected?

*A*: Coherent refers to rho production on the whole nucleus. Additional photon exchanges can excite either ion and emit neutrons. The 0n0n class has no forward neutrons, 0nXn includes either side alone, and XnXn requires neutrons on both sides. The inclusive result sums neutron classes within coherent production. The ALICE 5.02 TeV comparison uses the fitted rho Breit-Wigner contribution over $2m_\pi < m_{\pi\pi} < m_\rho + 5\Gamma_\rho$ and $|y|<0.8$, with the published $d\sigma/dy$ folded to $d\sigma/d|y|$. Detector selection cuts should not be reapplied to this corrected cross section. The inclusive result overlaps the exclusive classes and is excluded from their combined fit. See [the ALICE measurement](https://arxiv.org/abs/2002.10897).

`SCATTERING.LOOPSCREEN` enables ion survival. The `optical_ggcf` and `mc_ggcf` survival modes select the multichannel model named by `PARAM_NUCLEAR.GGCF.eikonal` in `GENERAL.json`, independently of the active production model. Plain `optical` uses the active soft model. Photonuclear `target_model: "glauber"` uses ordinary Glauber attenuation. Agreement in one neutron class does not validate the inclusive rate or all class fractions.

**Q**:  I get "error while loading shared libraries: libHepMC3.so.1: cannot open shared object file"

*A*:  Always remember the environment setup: `conda activate graniitti`, then `source ./install/setenv.sh`.

**Q**:  Where are the examples of .json steering cards?

*A*: Under `/gencard` and the physics family folders below `/icepack`. Make a copy and modify it.

**Q**: How do I choose the final state?

*A*: After the `->` arrow, write down the final state in terms of PDG numbers such as `22 22` or their names such as `gamma gamma`. Events are distributed according to the phase space if no amplitude is explicitly constructed for that particular decay. Any (cascade) final state is possible within the phase space.

> **Example**: `-> A > {A1 > {A11 A12} A2 A3} B > {B1 B2 > {B21 B22} B3}`.

**Q**: How do I select generic outgoing quark and gluon states?

*A*: Use `j` or `~j` in the final state and add `@j={u,d,s,c,b,g}`. The positive numeric forms `1,2,3,4,5,21` are also accepted, and one trailing comma is allowed. A lone alias includes both quark signs and a selected gluon. Immediate siblings in `j ~j` order select a correlated $q\bar q$ species or $gg$ pair. Generated pair backends reject the reverse order. Selection occurs before phase space generation, so heavy partons use the external mass scheme of the active amplitude. `PARAM_HARDPOMERON.parton_flavours` controls incoming hard Pomeron species only.

**Q**: How do I know the PDG numbers/names?

*A*: Write a dummy entry after the `->` arrow and try to generate events. You will get a list of all particles.

**Q**: What are the intermediate state (propagators, central systems, excited forward systems) PDG ID codes and the status codes in general?

*A*: Check `/include/Graniitti/Particle/MPDG.h`, `modeldata/mass_width_2026.mcd`, and `modeldata/TUNE0/PDG_EXTRA.json`.

**Q**: How do I get the minimum bias process cross sections (SD,DD,ND) and generate events?

*A*: Use `./bin/minbias`. The direct selectors are `X[SD]<Q>`, `X[DD]<Q>`, and `X[ND]<Q>`. Only `minbias` supplies the physical normalization for `X[ND]`, whose event generator amplitude is unit normalized. These events should normally be fragmented with Pythia. `SCATTERING.BEAMFRAG` in the generator card selects the forward fragmentation mode.

**Q**: How do I steer forward dissociation?

*A*: See the parameter block `PARAM_NSTAR` in `GENERAL.json` and use the generator card `SCATTERING.NSTARS` for activating the dissociative forward processes in central production.

**Q**: What do generation cuts versus fiducial cuts mean?

*A*: Generation cuts are applied at the level of phase space sampling and basic kinematic construction and are necessary. Fiducial cuts are optional and are applied after the event kinematics has been constructed.

**Q**:  For analysis or reference purposes, I want to generate a pure $J=0$ amplitude ($S$-wave) with a smooth invariant mass distribution, flat over rapidity and peripheral (exponential) $t$-distributions.

*A*: Add `@FLATAMP:X` at the end of your process string with `X = 1,2,3,4`, where the amplitudes are

- mode 1: $|A|^2 \propto e^{b(t_1+t_2)}$
- mode 2: $|A|^2 \propto e^{b(t_1+t_2)} / \sqrt{\hat{s}}$
- mode 3: $|A|^2 \propto e^{b(t_1+t_2)} / \hat{s}$
- mode 4: $|A|^2 = 1$

Modes 1 to 3 include a common overall $s^2$ factor. Set the cross section slope $b$ with `PARAM_FLAT.B` in `modeldata/<TUNE>/GENERAL.json`. This diagnostic amplitude replaces the selected physical matrix element and is unrelated to the phase space integrator.

**Q**:  I want to generate a tensor $J=2$ particle decaying to vector daughters with sequential spin correlations.

*A*: This is supported e.g. via Jacob-Wick helicity amplitudes (MP/XP/GP) or the TP model. First set up the mother resonance with a spin coupling or polarization parameters under `/modeldata/TUNE/RES/xxx.json` (create a new file). The command line `@MP_FRAME` choice affects only MP. XP and GP use their fixed center-of-mass basis. Second, adjust the decay $(l,s)$ couplings for the intermediate daughters under `/modeldata/TUNE/DECAYS.json`. Then create a process steering card by following the examples.


**Q**: Why does the integrated cross section seem very high?

*A*: Turn on the screening integral (event-by-event loop). The absorptive effect can be very high, depending on the process -- less for gamma-gamma processes in $pp$. Also, take care of the branching ratios and phase space volumes (see below). Some processes have couplings which may not be well tuned or have only placeholder values.

**Q**: Do forward protons have unlimited kinematic generation cuts by default?

*A*: No. Their transverse momentum is limited to a few GeV to improve CPU efficiency. You see all the cuts in the terminal ASCII output. You can change this by adding/modifying a field `"Pt"` under `"GENCUTS"`.

**Q**: Can I simulate processes as low in transverse momentum and mass as possible?

*A*: In principle, yes. The processes are numerically constructed to be floating point and kinematically safe. However, note that kt-EPA processes coupled with (MadGraph based) on-shell matrix elements may have some corners of phase space where the approximation becomes more prominent with respect to the full QED computation, which is available for lepton pair production.

For Durham QCD processes, the skewed unintegrated gluon densities $f_g(x,x',Q_\perp^2,\mu^2)$ probe $x'\ll x$, often at small $x$, while the loop may enter $Q_\perp^2<Q_0^2$ below the PDF matching scale. GRANIITTI therefore applies configurable loop-virtuality and hard-process cutoffs together with an infrared prescription. Predictions at small $x$, low $Q_\perp^2$, or low central mass are correspondingly model dependent. If you find unstable or unexpected behavior there, send an email.

**Q**: I want to know more details about the algorithms and calculations.

*A*: Very good! You can find my heavily commented code in the `/src` and `/include` folders. If something is unclear, get in contact. I am interested in what is right in terms of physics, what is wrong and in improving the code quality. All original calculations have inline arXiv references for more information or, if derived here, should be described in the code.

**Q**: Which features or processes are "under construction"?

*A*: See e.g. `/docs/underwork.md`. Many things are always in progress.


## MONTE CARLO TECHNIQUES

**Q**: Which integrator (importance sampler) should I use, NEUROJAC or VEGAS?

*A*: Start with `VEGAS` and try the `NEUROJAC` neural network (spline flows) when `VEGAS`, which factorizes dimension by dimension, does not converge easily (even when e.g. increasing `ncall`). The neural network can learn the correlation structure for more efficient importance sampling, but is slower to initialize and slower to sample from.

**Q**: I am running it on a GRID / distributed computing system. What do I need to do?

*A*: **Remember** to set the random number seed differently on separate instances/machines when generating events. Example: `./bin/gr -i mycard.json -r 12345`.

Also remember that different threads of GRANIITTI run with different seeds "grown" from the given master seed. Their values are shown on the screen. Thus, use significantly spaced master seeds on different HPC nodes. Also, use mode `"increment"` under `/modeldata/TUNE/NUMERICS.json`.

> [!TIP]
> It is useful to reuse the same VEGAS or NEUROJAC trained importance sampling proposal function data on every node with `-d vgrid/<name>.vgrid`. This is especially useful for NEUROJAC, because training a high quality neural network can take time. Keep a different, well-separated `-r` master seed on each node.

```bash

# First initialize on master node
./bin/gr -i card.json -g NEUROJAC -r 1234 -n 0 -l true -o myflow

# Reuse on independent HPC nodes
./bin/gr -i card.json -g NEUROJAC -r 1001 \
  -d vgrid/myflow.vgrid -n 10000 -l true -o node_001

./bin/gr -i card.json -g NEUROJAC -r 2002 \
  -d vgrid/myflow.vgrid -n 10000 -l true -o node_002
```

> Note. `-n -1` mode does not run the cross section integration stage or the max weight estimation for unweighting. It only computes the importance sampling adaptation grid. This can be practical to save time if only **weighted** events are to be generated later (e.g. in tuning applications).

**Q**: The MC integral does not converge (stat. error estimate is very high)!

*A*: This can happen if the generation cuts are too loose, the process is localized in a specific region of phase space, and VEGAS cannot find suitable binning. Manually tighten the generation cuts or try NEUROJAC. For example, if you generate an isolated resonance with mass $M$, use `<F>` phase space and put the mass boundaries around the resonance peak $M \pm 3 \Gamma$ full widths.

**Q**: What are `<C>` and `<F>` phase space classes?

*A*: They are different coordinate representations for the same Lorentz invariant phase space.

- `<C>` samples the true rapidity of every central particle and the invariant mass of the central system. The treatment of cascaded decays depends on whether the process has full end-to-end amplitudes or uses Jacob-Wick decay amplitudes. All integration variables are treated by the adaptive importance sampling routine, either by `NEUROJAC` or `VEGAS`.

- `<F>` samples the true rapidity and invariant mass of the complete central system before constructing its decay phase space, which is further treated by the `Rambo` algorithm. Thus, fewer variables undergo adaptive importance sampling than in the `<C>` case. In principle, this is less optimal, but it may yield greater stability if the learned proposal of `<C>`, which has more dimensions, is hard to optimize.

Both `Rap` and `M` ranges are required generation cuts in the active class. With equivalent physical coverage and identical fiducial cuts, the two classes should give the same integral, although their numerical efficiency can differ.

**Q**: I get different cross-section integrals with `<C>` and `<F>` algorithms.

*A*: First check that the generation domains describe the same physical region. `<C>.M` and `<F>.M` both bound the invariant mass of the complete central system, but their rapidity cuts are not the same observable. `<C>.Rap` bounds the rapidity of each particle, while `<F>.Rap` bounds the total central system rapidity. The system rapidity satisfies $\min_i y_i \leq Y_X \leq \max_i y_i$. Therefore a common `<C>.Rap` interval also contains the system rapidity, but the same `<F>.Rap` interval does not constrain each particle's rapidity. Rapidity cuts on later decay daughters must typically be imposed through fiducial cuts. Use generation ranges that fully contain the requested fiducial region in each coordinate system and then compare with identical fiducial cuts.

**Q**:  What is the difference between the `->` and `&>` operators?

*A*:  The `->` operator uses the full $2 \rightarrow N$ phase space volume in the cross section integration. The `&>` operator (with `<F>`, `<P>` phase space) decouples the central system phase space, and the process will be $[2 \rightarrow 3]$ + [unweighted $1 \rightarrow (N-2)$ decays].

In essence, the `&>` operator is useful to use with isolated scalar resonance processes, for example `yy[Higgs]<F>`, which is not a full SM end-to-end amplitude of Higgs production and its decays (such amplitudes can be generated with `MG2GRA`, though). It is intended for processes which do not contain any decay amplitude / matrix element decay part [check the code if unsure]. Typically, with these processes, the user then applies decay branching ratios manually. In general:

------------------

| Description               | `<C> ->` | `<F> ->` | `<F> &>` |
|---------------------------|----------|----------|----------|
| Amp with no decay algebra | -        | -        | yes      |
| Amp with decay algebra*   | yes      | yes      | -        |
------------------

*Decay algebra means relevant decay coupling constants, spin algebra etc.

**Q**: Can I generate an arbitrary multibody decay which is not described in `DECAYS.json`?

*A*: Yes, as an isolated phase-space decay, using `<F>` with the `&>` arrow. For example:

```bash
./bin/gr -i gencard/test.json -p "MP[RES]<F> &> pi+ pi- pi0 @RES{eta:1}" --set 'GENCUTS.<F>.M=[0.54784,0.547884]'
```

This produces the $\eta$ through the normal resonance amplitude and subsequently distributes its daughters uniformly in Lorentz invariant three-body phase space. The narrow mass interval is set according to the $\eta$ width. No three-body decay matrix element or branching fraction is inferred. Apply the branching fraction manually. A physical Dalitz structure would require an explicit two-body cascade or a dedicated $1 \rightarrow 3$ amplitude. The isolated `&>` arrow is not supported with `<C>`. This pure phase-space sampler does not validate decay quantum numbers, other than energy-momentum conservation.

**Q**: How are two-body decay couplings and LS coefficients normalized?

*A*: The decay coefficients `alpha_ls` are normalized to a unit squared Frobenius norm, and the overall absolute coupling is computed from, for example, the branching fraction and total decay width. All the production couplings are absolute.

**Q**: What do you mean by weighted vs unweighted events?

*A*: **Weighted** events have a positive definite weight associated with them. This comes naturally from the Monte Carlo sampling. Any histogramming etc. needs to be done by using these weights. You can find the weights in each HepMC3 event header. In the last event header, you can also find the total number of MC trials. The integrated fiducial cross section and its sampling uncertainty are obtained as

$$\widehat{\sigma}_{\mathrm{fid}} = \frac{1}{N_{\mathrm{trials}}}\sum_{i=1}^{N_{\mathrm{events}}} w_i$$

$$\delta\widehat{\sigma}_{\mathrm{fid}} = \frac{1}{N_{\mathrm{trials}}}\sqrt{\frac{N_{\mathrm{trials}}\sum_i w_i^2 - \left(\sum_i w_i\right)^2}{N_{\mathrm{trials}}-1}}$$

where $N_{\mathrm{trials}}$ need not equal $N_{\mathrm{events}}$.

**Unweighted** events have $w=1$ and generating them is based on acceptance-rejection sampling. The total fiducial cross section is provided in the HepMC3 event header.

**Q**: Should I use weighted events?

*A*: While generating weighted events is fast and can be useful e.g. in tuning or training statistical inference surrogates, the sample may have a (very) low effective statistical sample size (ESS). The effective sample size for $N$ weighted events is given by

ESS $= (\sum_i w_i)^2 / \sum_i w_i^2 \leq N$,

which is printed to the terminal.

Thus, typically **unweighted** events are used, because they carry unit statistical power per event and e.g. GEANT4 CPU hours & disk space are allocated equally to each event, rather than wasted on events with very small statistical power.

**Q**: After unweighted event generation, I get maximum weight overflow warnings. What does this mean?

*A*: The phase before event generation samples the phase space to obtain an **estimate** of the maximum integrand weight needed for acceptance-rejection sampling in unweighted event generation. If you get maximum weight overflows, it simply means that the event generation phase has found larger integrand values. This is not good, even if in practice it may not always be extremely serious. To find a better maximum weight, you should adjust VEGAS parameters, for example, by using ten times as many samples per iteration and a relative error target ten times smaller.

**Q**: How many CPU threads can GRANIITTI handle?

*A*: As many as your system can. You can try out which number gives the best performance.


## ANALYSIS

**Q**: How do I quickly plot the fiducial observables of a given MC sample and stack different MC samples in the same plot?

*A*: PYTHON: Use `python -m core.iceplot --hepmc3 file1 file2`

ROOT:   With ROOT installed, just try some examples under `/tests`. The program `./bin/analyze` calculates a large collection of automated publication quality plots for arbitrary central production processes.

RIVET:  We also recommend trying RIVET.

**Q**: How do I use RIVET?

*A*: Install it with `/install/installRIVET.sh` and see examples under `./tests/external/rivet/`.

**Q**:  I want to do a full spherical harmonic expansion.

*A*: Use `./bin/fitharmonic --card <analysis.json>`. The response is factorized into generated space with a flat angular distribution, particle level fiducial space, and reconstructed detector space.

**Q**: I get different spherical harmonic *acceptance inversion* results in different Lorentz rest frames.

*A*: This is expected and is due to a finite spherical harmonic expansion order, say `L_max=4`. However, in the limit of infinite MC statistics and infinite harmonic expansion order, different frames will expand the acceptance identically and the data inversion results should then also be identical in that respect. However, the harmonic coefficients themselves can (will) naturally differ in different frames.

> The angular expansion is a harmonic fit in $(\cos\theta,\phi)$ conditional on the selected kinematic coordinates (hyperbin) of the event. In some eventwise rotated Lorentz rest frames, the geometric acceptance is flatter and requires fewer harmonic components. Comparing two detectors at the same $\sqrt{s}$, the detector with larger $(\eta,p_T)$ acceptance generally requires fewer acceptance harmonics. The production dynamics still determine the angular order needed for the physics distribution.


## ADVANCED

**Q**: How do I see quantum entanglement metrics for resonance decays?

*A*: Try this example:

```bash
./bin/gr -i gencard/test.json -p "GP[RES]<F> -> rho(770)0 rho(770)0 @RES{f0_1710:1} @RES{f2_2150:1} @QMETRICS"
```

**Q**: Which one of the Pomeron models is the most flexible? Which one is the simplest?

*A*: `GP` and `XP/MP` offer broad flexibility through independent helicity or $(l,s)$ couplings that control spin correlations. `GP` uses analytic Regge angular momenta $\alpha_i(t_i)$, introducing trajectory dependence into the spin structure, while `XP/MP` use fixed exchange spins. These reduced couplings do not uniquely specify the off-shell interaction. `TP` describes the Pomeron as an effective rank-2 tensor exchange with explicit covariant vertices and propagators. The chosen Lorentz structures relate momentum dependence and helicity amplitudes. `MP` with minimal production couplings and prescribed resonance polarization is usually the simplest phenomenological starting point.

**Q**: What is the difference between `GP[RES]`, `XP[RES]`, `MP[RES]`, and `TP[RES]`?

*A*:

- `GP[RES]` is the analytic Regge helicity process and reads the `GP` model block. The exchange spins are $\alpha_i(t_i)$, while the `MMAX` parameter truncates the retained integer helicity labels. The coupling option `helicity` describes independent amplitude couplings and completes their (parity) symmetry orbits. `g_ls` uses analytically continued angular coefficient functions.

- `MP[RES]` is the Minimal Pomeron resonance process and reads the `MP` model block. In addition to various "Minimal Coupling" operator structure modes, the resonance spin state can be dressed with coherent `a_Jz` amplitudes or an incoherent spin-density matrix `rho`.

- `XP[RES]` is a fixed-spin Pomeron process with a reduced (simplified) covariant structure and reads the `XP` model block, with options `helicity` and `g_ls` for the coupling structure.

- `TP[RES]` is the covariant Tensor Pomeron process and reads the `TP` model block (see the original papers for reference).

**Q**: What is the polarization steering option in `MP` and how does it relate to the helicity or $(l,s)$ couplings?

*A*: The production helicity or $(l,s)$ couplings already determine the resonance polarization dynamically, event by event, through the amplitude and kinematics. Optional `MP` steering applies an additional phenomenological spin filter, using coherent `a_Jz` amplitudes or a general density matrix `rho` in the chosen `MP_FRAME`. In `rho` mode, each resonance is added incoherently to the continuum and other resonances, while its internal coupling terms remain coherent. The resulting event polarization depends on both the production amplitude and the filter.

**Q**: What are the `auto` coupling basis modes in `MP`?

*A*: The overall complex coupling is set by `g, phi` and the "minimal" coupling structures chosen by `auto_` modes are:

- `auto_min_L`: Select the allowed $(l,s)$ pair with the smallest $l$, then the smallest $s$.

- `auto_min_S`: Select the allowed $(l,s)$ pair with the smallest $s$, then the smallest $l$.

- `auto_equal_ls`: Assign equal magnitudes and phases to all allowed $(l,s)$ couplings.

- `auto_equal_helicity`: Construct equal helicity magnitudes with the required parity symmetry phases, projected onto the allowed coupling space.

With production spin correlations enabled, the terms are combined coherently, $C=\sum_{l,s}g_{ls}C_{ls}$, and optional polarization steering gives $C_{\mathrm{pol}}=C S^T$, with $S=\sqrt{(2J+1)\rho}$. Here $C_{ls}$ are production matrices with resonance spin columns, $J$ is the resonance spin, and $\rho$ is defined in `MP_FRAME`. For normalized coherent `a_Jz` amplitudes, $\rho=aa^\dagger$. The square root is the positive matrix square root, and $S=I$ without steering.

**Q**: Are GP/XP/MP/TP continuum and resonance amplitudes coherent?

*A*: Contributions in the same coherent sector are added before squaring and screening, e.g. `GP[CON]` and `GP[RES]` in `GP[RES+CON]`. In `MP`, each resonance using `rho` (density matrix) steering is added incoherently to the continuum and other resonances.

**Q**: Why does the production rate of $f_0(980)$ or another resonance not match data (or a previous version of GRANIITTI), e.g. in $K^+K^-$?

*A*: Branching fractions for some resonances, such as $f_0(980)$, are uncertain and are read from `modeldata/TUNE0/DECAYS.json`. Production normalizations are explicit tune parameters. GP, MP, and XP store their absolute complex couplings directly in `g`, `g_ls` or `helicity` under each pair key. TP uses its own covariant couplings. Equal numerical inputs provide a model comparison but do not constitute a common fit because the central operators differ.

**Q**: Why do distributions match, for example, STAR data at $\sqrt{s} = 0.2$ TeV but not CMS data at $\sqrt{s} = 13$ TeV?

*A*: This is most likely an issue with the underlying model, for example, missing Reggeon contributions, the eikonal screening loop in central production, or the global TUNE settings (couplings, continuum form factors ...).
 
The differential efficiency (or background) corrections in the data should also be carefully scrutinized across experiments, e.g. for forward protons and central final states. For these, we recommend a fully differential (multidimensional) treatment such as DeepEfficiency, which is a strategy designed to be as independent of the event generator as possible.

In principle, this allows efficiency corrections in the full phase space, taking into account both the central final states (e.g. $2 \times 3$-momentum) and the kinematic variables of forward protons. This is not possible with traditional histogramming techniques.

**Q**: How are spin correlations configured in sequential decays?

*A*: Set intermediate decay couplings in `modeldata/<TUNE>/DECAYS.json`. Resonance production is in `RES/*.json`. Continuum spin couplings are separated into `CON_MP.json`, `CON_XP.json`, and `CON_GP.json`.

**Q**: How do full cascade amplitudes differ from generic Jacob-Wick decays?

*A*:

| Amplitude | Decay treatment |
| --- | --- |
| Full MG5 and Tensor Pomeron cascades | Complete amplitudes end-to-end (except some special cases for TP). |
| Generic Jacob-Wick processes, including Higgs and monopolium | Helicity decay matrices through binary chains of arbitrary depth. |

**Q**: How do I adjust the numerical parameters of the screening loop integral?

*A*: This depends on the tradeoff between accuracy and computation time and can have a major impact on certain observables. Briefly, this belongs to the category of `sign problems` in quantum mechanics / QFT because of the very high sensitivity to the numerical matching of complex phases.

In general, one should try increasingly fine discretizations for the loop momentum magnitude and angle under `/modeldata/TUNE/NUMERICS.json`, and check when the differential results for various observables stop changing.

**Q**: Can I add some arbitrary new process / amplitude?

*A*: Yes, you can. The kinematics is constructed to be fully generic and obey standard QFT normalization rules. What you need is some knowledge of C++ and then you need to study a little bit how the existing amplitudes are implemented in the machinery. Basically, you can add them under the `MSubProc` 'umbrella' class.

**Q**: How do I add new amplitudes / matrix elements based on MadGraph?

*A*: See the automatic import script under `./develop/MG2GRA` and `docs/madgraph.md`.

**Q**: Are the soft central production parameters independent of the screening loop?

*A*: No, the screening loop changes e.g. the phase structure of amplitudes. Strictly speaking, one should have two sets of parameters (or tunes), one with screening turned on and one without. Alternatively, use just one set with screening always turned on (computationally very demanding for some processes, though).

**Q**: How well are the soft parameters tuned/fitted?

*A*: The soft central production amplitude parameters are such that an approximate match with LHC and RHIC data can be obtained. More high precision data are needed for both fully exclusive and semi-exclusive proton dissociative processes. The same applies to many PDG resonance branching ratios. See `modeldata/<TUNE>/DECAYS.json` and the model specific `CON_*.json` cards.

Elastic (eikonal Pomeron) parameters are tuned with better accuracy to high energy elastic $d\sigma/dt$ and provide a reasonable starting point for their use in screening corrections. However, the available multichannel eikonals probably require simultaneous fits to inelastic diffractive data (either inclusive or exclusive) to be most useful in practice. Similarly, photoproduction model parameters are more accurate thanks to HERA and LHC data.

**Q**: What is the difference between Jacob-Wick LS couplings and covariant current models?

*A*: LS amplitudes provide a general angular basis, while covariant currents add model-dependent full dynamics, in short.

**Q**: Can I make my own tune?

*A*:  MANUAL: Sure, just copy `/modeldata/TUNE0` to e.g. `/modeldata/TUNE1`, update the reference in your .json steering card and start tuning. Then go beyond and put in your own scattering amplitudes.

AUTOMATED: For distributed HPC tuning via `icetune`, see `/docs/icetune.md`. Basically, the tuning algorithms modify the input parameter JSON files, generate new MC samples, construct the MC histograms from HepMC3 files and compare these with HEPData histograms, and then repeat this parameter optimization process by minimizing global chi2 or some other objective in an efficient way. Uncertainties can be extracted later with `icescape`.

> How to do this optimally is a research topic in its own right: active learning, Bayesian optimization, black-box optimization, surrogate models, likelihood free and simulation-based inference.

## OTHER

**Q**:  Your .json files do not follow the JSON standard.

*A*:  I use a subset of the extended `JSON5` format which allows commenting and other user friendly properties. You can read them with JSON5 reader libraries or with standard JSON libraries after some regexp cleaning. In addition, GRANIITTI uses custom `$ref` field references (links) within and across cards.

**Q**: I want to install it into our collaboration software platform.

*A*: Great! However, do NOT embed or create custom interfaces to the generator at the level of GRANIITTI classes. Instead, use it only through .json files and the command line interface, and take the HepMC3/2 output. First of all, this allows very streamlined and fast updates. Second, the generator engine is fully multithreaded and a high level of expertise is required for any safe modifications or interfacing at the level of code. Just keep the package factorized!

**Q**: My icetune Condor job exits immediately with `GLIBCXX_3.4.31 not found`.

*A*: The worker is loading the system `libstdc++`, while `bin/gr` was built e.g. against the newer Conda runtime. After activating the `graniitti` environment, use

```bash
conda activate graniitti
CONDA_CXX_ABI=1 source install/setenv.sh
```

This gives `$CONDA_PREFIX/lib` priority over `/lib64`. If the error persists, rebuild `bin/gr` with the target worker toolchain.


**Q**: Something seems weird.

*A*: Send an email. It might be a bug or just a confusing user interface.


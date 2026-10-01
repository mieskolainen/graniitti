# H1 elastic pion pair photoproduction

This icepack compares `TP[PHOTO]` with the elastic H1 $\gamma p \to \pi^+\pi^-p$. The generator uses 27.6 GeV positrons and 920 GeV protons ($e^+p$), as specified in [Section 3.2 of the H1 paper](https://arxiv.org/pdf/2005.14471). The reader obtains the effective photon flux $\Phi_{\gamma/e}$ from the original HEPData qualifier and supplies the conversion of the generated $ep$ rate to the $\gamma p$ cross section. The inputs are original HEPData JSON tables for [Tables 17, 18-19 and 22-23 with the available statistical correlations](https://www.hepdata.net/record/ins1798511?version=1), including the published signed systematic variations. The mass spectra include the measured $W$ and $|t|$ intervals, each compared with DS, DS plus $\rho/\omega$, and DS plus $\rho/\omega/\rho_3$.

The current coherent `TP[PHOTO]` calculation is restricted to $Q^2<0.5$ GeV $^2$.

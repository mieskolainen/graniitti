# Durham dijets with GRANIITTI and Pythia

GRANIITTI generates exclusive Durham $gg \to gg$ and $gg \to u\bar{u}$ events in LHE format. Pythia 8 showers and hadronizes the events and writes HepMC3.

Use the common [icepack instructions](../../../tests/README.md) and [Pythia setup and comparisons](../../../tests/external/pythia/README.md).

## Fiducial phase space

The generator cuts define a loose hard process phase space, not the particle level selection. Hard partons require $p_T>10$ GeV with no generator level $\eta$ cut. The central system range is $40<M_X<400$ GeV and $|y_X|<3.5$. After Pythia, stable visible particles are clustered with anti-$k_T$, $R=0.6$. The fiducial selection requires two jets with $p_T>20$ GeV and $|\eta|<2.5$, followed by $40<m_{jj}<200$ GeV. The loose parton threshold and rapidity range retain shower and hadronization migration. The generator upper mass is twice the fiducial $m_{jj}$ limit because radiation outside the leading jets gives $m_{jj}<M_X$.

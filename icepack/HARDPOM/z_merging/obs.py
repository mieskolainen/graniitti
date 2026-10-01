# Dimuon spectra and hadronic kT jet resolution for the SD merging study
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import fastjet
import numpy as np

from ..z_with_pythia import obs as dimuon


# Compute the final hadronic kT transition from one jet to zero within |eta| < 5
def kt01(event):
    particles = []
    for p in event.evt.particles():
        if p.status() != 1 or abs(p.pid()) in (12, 13, 14, 16):
            continue
        q = p.momentum()
        if abs(q.eta()) < 5:
            particles.append(fastjet.PseudoJet(q.px(), q.py(), q.pz(), q.e()))
    if not particles:
        return 0.0
    sequence = fastjet.ClusterSequence(particles, fastjet.JetDefinition(fastjet.kt_algorithm, 1.0))
    return np.sqrt(max(0.0, sequence.exclusive_dmerge(0)))


obs_m_mumu = dict(dimuon.obs_m_mumu, xlim=(80, 100), bins=np.linspace(80, 100, 21))
obs_pt_mumu = dict(dimuon.obs_pt_mumu, xlim=(0, 150), bins=np.linspace(0, 150, 31))
obs_y_mumu = dict(dimuon.obs_y_mumu)
obs_kt01 = dict(
    obs_pt_mumu,
    tag="kt01",
    func=kt01,
    xlabel=r"$k_{T,01}$",
    ylabel=r"$d\sigma/dk_{T,01}$",
    label="Hadronic kT jet resolution",
)

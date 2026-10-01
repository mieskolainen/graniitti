# Jacob-Wick polarization observables for central resonances
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numpy as np
from core.analysis import obs
from core.analysis.observables.default import obs_costheta_CS as default_costheta_CS
from core.analysis.observables.default import obs_dPhi_pp as default_dPhi_pp
from core.analysis.observables.default import obs_phi_CS as default_phi_CS


# Compute the analyzer glueball filter from the forward proton momentum difference
def project_glueball_filter(event):
    _, final_protons = obs.proj_init_final_protons(event)
    return (final_protons[0] - final_protons[1]).pt


obs_costheta_CS = copy.deepcopy(default_costheta_CS)
obs_costheta_CS["bins"] = np.linspace(-1.0, 1.0, 41)

obs_phi_CS = copy.deepcopy(default_phi_CS)
obs_phi_CS["xlim"] = (-180.0, 180.0)
obs_phi_CS["bins"] = np.linspace(-180.0, 180.0, 41)

obs_dPhi_pp = copy.deepcopy(default_dPhi_pp)
obs_dPhi_pp["bins"] = np.linspace(0.0, 180.0, 41)

obs_glueball_filter = {
    "tag": "glueball_filter",
    "xlim": (0.0, 2.0),
    "ylim": None,
    "xlabel": r"$|\Delta \vec{p}_{T,pp}|$",
    "ylabel": r"$1/\sigma\,d\sigma/d|\Delta \vec{p}_{T,pp}|$",
    "units": {"x": r"GeV", "y": r"GeV$^{-1}$"},
    "label": r"Glueball filter",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 2.0, 41),
    "density": False,
    "func": project_glueball_filter,
}

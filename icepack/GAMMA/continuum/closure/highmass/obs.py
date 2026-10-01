# High mass dimuon continuum closure observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numpy as np
from core.analysis import obs
from core.analysis.observables.default import obs_dPhi_pp as default_dPhi_pp
from core.analysis.observables.default import obs_M as default_M
from core.analysis.observables.default import obs_Rap as default_Rap

from ....ATLAS_1377585.mumu.obs import proj_acoplanarity


# Compute the dimuon system transverse momentum
def project_pair_pt(event):
    return sum(obs.proj_central_particles(event)).pt


obs_M = copy.deepcopy(default_M)
obs_M["xlim"] = (20.0, 200.0)
obs_M["bins"] = np.linspace(20.0, 200.0, 37)

obs_Rap = copy.deepcopy(default_Rap)
obs_Rap["xlim"] = (-2.4, 2.4)
obs_Rap["bins"] = np.linspace(-2.4, 2.4, 33)

obs_Pt = {
    "tag": "Pt",
    "xlim": (0.0, 5.0),
    "ylim": None,
    "xlabel": r"$p_{T,\mu\mu}$",
    "ylabel": r"$d\sigma/dp_{T,\mu\mu}$",
    "units": {"x": r"GeV", "y": r"pb"},
    "label": r"Dimuon transverse momentum",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 5.0, 41),
    "density": False,
    "func": project_pair_pt,
}

obs_acoplanarity = {
    "tag": "acoplanarity",
    "xlim": (0.0, 0.052),
    "ylim": None,
    "xlabel": r"$1-|\Delta\phi_{\mu\mu}|/\pi$",
    "ylabel": r"$d\sigma/dA_\phi$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Dimuon acoplanarity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 0.052, 27),
    "density": False,
    "func": proj_acoplanarity,
}

obs_dPhi_pp = copy.deepcopy(default_dPhi_pp)

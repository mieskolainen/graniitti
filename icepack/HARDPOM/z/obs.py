# Compact dimuon observables for standalone Z production
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numpy as np
from core.analysis import obs
from core.analysis.observables.default import obs_M as default_M
from core.analysis.observables.default import obs_Rap as default_Rap


# Compute the transverse momentum of the selected dimuon system
def project_dimuon_pt(event):
    return sum(obs.proj_central_particles(event)).pt

obs_m_mumu = copy.deepcopy(default_M)
obs_m_mumu["tag"] = "m_mumu"
obs_m_mumu["xlim"] = (80.0, 100.0)
obs_m_mumu["bins"] = np.linspace(80.0, 100.0, 41)
obs_m_mumu["xlabel"] = r"$m_{\mu\mu}$"

obs_y_mumu = copy.deepcopy(default_Rap)
obs_y_mumu["tag"] = "y_mumu"
obs_y_mumu["xlim"] = (-6.0, 6.0)
obs_y_mumu["bins"] = np.linspace(-6.0, 6.0, 49)
obs_y_mumu["xlabel"] = r"$y_{\mu\mu}$"

obs_pt_mumu = {
    "tag": "pt_mumu",
    "xlim": (0.0, 100.0),
    "ylim": None,
    "xlabel": r"$p_{T}^{\mu\mu}$",
    "ylabel": r"$d\sigma/dp_{T}^{\mu\mu}$",
    "units": {"x": r"GeV", "y": r"pb"},
    "label": r"Dimuon transverse momentum",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 100.0, 51),
    "density": False,
    "func": project_dimuon_pt,
}

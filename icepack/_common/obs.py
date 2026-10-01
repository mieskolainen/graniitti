# Common event observables for icepack physics validations
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numpy as np
from core.analysis.observables.default import obs_costheta_CS as default_costheta_CS
from core.analysis.observables.default import obs_dPhi_pp as default_dPhi_pp
from core.analysis.observables.default import obs_M as default_M
from core.analysis.observables.default import obs_Rap as default_Rap


# Map every accepted event into one unit width cross section bin
def proj_cross_section(_event):
    return 0.5

obs_cross_section = {
    "tag": "cross_section",
    "xlim": (0.0, 1.0),
    "ylim": None,
    "xlabel": r"Fiducial selection",
    "ylabel": r"$\sigma$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Fiducial cross section",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.array([0.0, 1.0]),
    "density": False,
    "func": proj_cross_section,
}

obs_M = copy.deepcopy(default_M)
obs_Rap = copy.deepcopy(default_Rap)
obs_dPhi_pp = copy.deepcopy(default_dPhi_pp)
obs_costheta_CS = copy.deepcopy(default_costheta_CS)

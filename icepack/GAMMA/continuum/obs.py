# Photon continuum observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numpy as np
from core.analysis.observables.default import obs_costheta_CS as default_costheta_CS
from core.analysis.observables.default import obs_dPhi_pp as default_dPhi_pp
from core.analysis.observables.default import obs_M as default_M
from core.analysis.observables.default import obs_Rap as default_Rap

obs_M = copy.deepcopy(default_M)
obs_M["xlim"] = (20.0, 200.0)
obs_M["bins"] = np.linspace(20.0, 200.0, 46)

obs_Rap = copy.deepcopy(default_Rap)
obs_costheta_CS = copy.deepcopy(default_costheta_CS)
obs_dPhi_pp = copy.deepcopy(default_dPhi_pp)

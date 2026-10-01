# Photon fusion Higgs observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numpy as np
from core.analysis.observables.default import obs_M as default_M
from core.analysis.observables.default import obs_Rap as default_Rap

obs_M = copy.deepcopy(default_M)
obs_M["xlim"] = (124.0, 126.0)
obs_M["bins"] = np.linspace(124.0, 126.0, 41)

obs_Rap = copy.deepcopy(default_Rap)

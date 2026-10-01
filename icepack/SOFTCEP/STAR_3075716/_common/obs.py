# STAR 510 GeV central pair and tagged proton observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from copy import deepcopy

from core.analysis.observables import default

obs_M = deepcopy(default.obs_M)
obs_M.update(xlim=None, bins=None)
obs_Rap = deepcopy(default.obs_Rap)
obs_Rap.update(xlim=None, bins=None)
obs_Abs_t1t2 = deepcopy(default.obs_Abs_t1t2)
obs_Abs_t1t2.update(xlim=None, bins=None)
obs_dPhi_pp = deepcopy(default.obs_dPhi_pp)
obs_dPhi_pp.update(xlim=None, bins=None)

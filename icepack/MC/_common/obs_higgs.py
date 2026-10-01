# Parton level Higgs observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np

from .observables import observable, production

globals().update(production(115.0, 135.0))
obs_costheta_CS = observable("costheta_CS", -1.0, 1.0)

# Resolve the Higgs pole with 1 MeV bins and retain both wide mass tails
obs_M = observable("M", 125.10, 125.16, 60)
obs_M["bins"] = np.concatenate(([115.0], obs_M["bins"], [135.0]))

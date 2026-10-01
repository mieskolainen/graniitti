# Dijet mass and rapidity over the configured generator phase space
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
from pathlib import Path

import numpy as np
import pyjson5
from core.analysis.observables import default
from core.io.serialize import load_json_file

phase_space = load_json_file(Path(__file__).with_name("gencard_IPIP_jj.json"), loader=pyjson5.load)["GENCUTS"]["<C>"]
obs_M = copy.deepcopy(default.obs_M)
obs_Rap = copy.deepcopy(default.obs_Rap)
for observable in (obs_M, obs_Rap):
    observable["xlim"] = tuple(phase_space[observable["tag"]])
    observable["bins"] = np.linspace(*observable["xlim"], len(observable["bins"]))

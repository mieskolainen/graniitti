# ALICE coherent rho0 production rapidity selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

cut_param = {"abs_rap_max": 0.8}


# Apply the published ALICE rapidity interval
def cut_func(event):
    return abs(obs.proj_1D_Rap(event)) < event.cut_param["abs_rap_max"]

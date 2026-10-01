# ALICE incoherent J/psi central rapidity and momentum transfer selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs as observables

from . import obs

cut_param = {"abs_rap_max": 0.8, "abs_t": [0.04, 1.0]}


# Apply the published rapidity and momentum transfer intervals
def cut_func(event):
    rapidity = abs(observables.proj_1D_Rap(event))
    abs_t = obs.proj_abs_t(event)
    lower, upper = event.cut_param["abs_t"]
    return rapidity < event.cut_param["abs_rap_max"] and lower < abs_t < upper

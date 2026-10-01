# ALICE coherent J/psi production rapidity selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

cut_param = {"rap": [[-0.2, 0.8], [2.5, 4.0]]}


# Apply the published ALICE central and forward rapidity intervals
def cut_func(event):
    rapidity = obs.proj_1D_Rap(event)
    return any(lower < rapidity < upper for lower, upper in event.cut_param["rap"])

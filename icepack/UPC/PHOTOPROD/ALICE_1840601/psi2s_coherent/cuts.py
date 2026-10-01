# ALICE coherent psi(2S) production rapidity selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

cut_param = {"rap": [0.0, 0.8]}


# Apply one side of the published symmetric ALICE rapidity interval
def cut_func(event):
    rapidity = obs.proj_1D_Rap(event)
    lower, upper = event.cut_param["rap"]
    return lower < rapidity < upper

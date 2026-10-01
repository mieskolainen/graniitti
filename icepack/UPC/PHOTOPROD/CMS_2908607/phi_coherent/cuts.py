# CMS coherent phi production rapidity selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

cut_param = {"abs_rap": [0.3, 1.0]}


# Apply the published CMS absolute rapidity interval
def cut_func(event):
    rapidity = abs(obs.proj_1D_Rap(event))
    lower, upper = event.cut_param["abs_rap"]
    return lower < rapidity < upper

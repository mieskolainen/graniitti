# CMS exclusive rho production selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

cut_param = {"rap": [-2.0, 2.0], "pt2": [0.025, 1.0]}


# Select meson laboratory rapidity and reconstructed momentum transfer
def cut_func(event):
    low, high = event.cut_param["rap"]
    tmin, tmax = event.cut_param["pt2"]
    return low < obs.proj_1D_Rap(event) < high and tmin < obs.proj_1D_Pt(event)**2 < tmax

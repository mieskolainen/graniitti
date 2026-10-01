# LHCb coherent rho production selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

cut_param = {"rap": [2.05, 4.90]}


# Select the published meson rapidity interval after acceptance correction
def cut_func(event):
    low, high = event.cut_param["rap"]
    return low < obs.proj_1D_Rap(event) < high

# CMS coherent J/psi rapidity observable
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache


# Compute the absolute J/psi rapidity
@cache
def proj_abs_rap(event):
    return abs(obs.proj_1D_Rap(event))


obs_abs_rap = {
    "tag": "abs_rap",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$|y_{J/\psi}|$",
    "ylabel": r"$d\sigma/dy$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"CMS coherent $J/\psi$",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_abs_rap,
}

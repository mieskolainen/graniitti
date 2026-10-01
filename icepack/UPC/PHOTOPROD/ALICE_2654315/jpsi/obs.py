# ALICE pPb J/psi rapidity observable
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache


# Compute the J/psi rapidity along the proton beam direction
@cache
def proj_rap(event):
    return obs.proj_1D_Rap(event)


obs_rap = {
    "tag": "rap",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$y_{J/\psi}$",
    "ylabel": r"$d\sigma/dy$",
    "units": {"x": r"unit", "y": r"ub"},
    "label": r"ALICE pPb $J/\psi$",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_rap,
}

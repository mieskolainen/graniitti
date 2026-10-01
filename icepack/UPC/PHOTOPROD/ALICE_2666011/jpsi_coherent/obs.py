# ALICE coherent J/psi rapidity observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache


# Compute the J/psi rapidity
@cache
def proj_rap(event):
    return obs.proj_1D_Rap(event)


obs_rap = {
    "tag": "rap",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$y_{J/\psi}$",
    "ylabel": r"$d\sigma_{J/\psi}/dy$",
    "units": {"x": r"unit", "y": r"mb"},
    "label": r"ALICE coherent $J/\psi$",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_rap,
}

obs_rap_prod = {
    "tag": "rap_prod",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$y_{J/\psi}$",
    "ylabel": r"$d\sigma_{J/\psi}/dy_{J/\psi}$",
    "units": {"x": r"unit", "y": r"mb"},
    "label": r"ALICE coherent $J/\psi$",
    "ylim_ratio": None,
    "bins": None,
    "density": False,
    "func": proj_rap,
}

# CMS coherent Upsilon(1S) rapidity observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache


# Compute the absolute Upsilon(1S) rapidity
@cache
def proj_abs_rap(event):
    return abs(obs.proj_1D_Rap(event))


obs_abs_rap = {
    "tag": "abs_rap",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$|y_{\Upsilon(1S)}|$",
    "ylabel": r"$d\sigma_{\Upsilon(1S)}/dy$",
    "units": {"x": r"unit", "y": r"ub"},
    "label": r"CMS coherent $\Upsilon(1S)$",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_abs_rap,
}

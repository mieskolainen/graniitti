# CMS coherent phi rapidity observable
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache


# Compute the absolute phi rapidity
@cache
def proj_abs_rap(event):
    return abs(obs.proj_1D_Rap(event))


obs_abs_rap = {
    "tag": "abs_rap",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$|y_{\phi}|$",
    "ylabel": r"$d\sigma_{\phi}/dy$",
    "units": {"x": r"unit", "y": r"mb"},
    "label": r"CMS coherent $\phi$",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_abs_rap,
}

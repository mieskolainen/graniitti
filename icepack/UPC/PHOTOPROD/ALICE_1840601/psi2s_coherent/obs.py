# ALICE coherent psi(2S) rapidity observable
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache


# Compute the psi(2S) rapidity
@cache
def proj_rap(event):
    return obs.proj_1D_Rap(event)


obs_rap = {
    "tag": "rap",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$y_{\psi(2S)}$",
    "ylabel": r"$d\sigma_{\psi(2S)}/dy$",
    "units": {"x": r"unit", "y": r"mb"},
    "label": r"ALICE coherent $\psi(2S)$",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_rap,
}

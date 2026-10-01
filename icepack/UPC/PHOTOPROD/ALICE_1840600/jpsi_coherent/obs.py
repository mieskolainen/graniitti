# ALICE coherent J/psi UPC transverse momentum observable
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache


# Compute the squared transverse momentum of the J/psi decay system
@cache
def proj_pt2(event):
    value = obs.proj_1D_Pt(event)
    return value * value


obs_pt2 = {
    "tag": "pt2",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$p_{T,J/\psi}^{2}$",
    "ylabel": r"$d^2\sigma_{J/\psi}/(dy\,dp_T^2)$",
    "units": {"x": r"GeV$^2$", "y": r"pb"},
    "label": r"ALICE coherent $J/\psi$",
    "ylim_ratio": (0.0, 2.5),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_pt2,
}

# ALICE pPb exclusive dimuon mass observable
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache


# Compute the dimuon invariant mass
@cache
def proj_mass(event):
    return obs.proj_1D_M(event)


obs_mass = {
    "tag": "mass",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$m_{\mu\mu}$",
    "ylabel": r"$d\sigma/dm_{\mu\mu}$",
    "units": {"x": r"GeV", "y": r"ub"},
    "label": r"ALICE pPb exclusive dimuons",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_mass,
}

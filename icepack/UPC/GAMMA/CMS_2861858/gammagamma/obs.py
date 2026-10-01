# CMS PbPb light by light Figure 7 observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.io.cache import cache

from .cuts import diphoton


# Compute the signed diphoton rapidity
@cache
def proj_rapidity(event):
    return diphoton(event).rapidity


# Compute the diphoton invariant mass
@cache
def proj_mass(event):
    return diphoton(event).m


obs_rapidity = {
    "tag": "rapidity",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$y^{\gamma\gamma}$",
    "ylabel": r"$d\sigma^{\gamma\gamma}/dy^{\gamma\gamma}$",
    "units": {"x": r"unit", "y": r"nb"},
    "label": r"CMS UPC light by light rapidity",
    "ylim_ratio": (0.0, 3.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_rapidity,
}

obs_mass = {
    "tag": "mass",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$m^{\gamma\gamma}$",
    "ylabel": r"$d\sigma^{\gamma\gamma}/dm^{\gamma\gamma}$",
    "units": {"x": r"GeV", "y": r"nb/GeV"},
    "label": r"CMS UPC light by light invariant mass",
    "ylim_ratio": (0.0, 3.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_mass,
}

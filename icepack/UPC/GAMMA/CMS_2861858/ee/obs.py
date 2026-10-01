# CMS PbPb Breit-Wheeler Figure 4 observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

from core.io.cache import cache

from .cuts import dielectron, selected_electrons


# Compute the dielectron transverse momentum
@cache
def proj_pair_pt(event):
    return dielectron(event).pt


# Compute the signed dielectron rapidity
@cache
def proj_rapidity(event):
    return dielectron(event).rapidity


# Compute the dielectron invariant mass
@cache
def proj_mass(event):
    return dielectron(event).m


# Compute the absolute Collins Soper scattering angle from the leading electron
@cache
def proj_abs_costheta(event):
    system = dielectron(event)
    electron = selected_electrons(event)[0]["p4"]
    transverse_mass = math.sqrt(max(0.0, system.e**2 - system.pz**2))
    denominator = system.m * transverse_mass
    if denominator <= 0.0:
        return 0.0
    numerator = 2.0 * abs(system.e * electron.pz - system.pz * electron.e)
    return min(1.0, numerator / denominator)


obs_pair_pt = {
    "tag": "pair_pt",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$p_{T}^{ee}$",
    "ylabel": r"$d\sigma^{ee}/dp_{T}^{ee}$",
    "units": {"x": r"GeV", "y": r"ub/GeV"},
    "label": r"CMS UPC dielectron transverse momentum",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_pair_pt,
}

obs_rapidity = {
    "tag": "rapidity",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$y^{ee}$",
    "ylabel": r"$d\sigma^{ee}/dy^{ee}$",
    "units": {"x": r"unit", "y": r"ub"},
    "label": r"CMS UPC dielectron rapidity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_rapidity,
}

obs_mass = {
    "tag": "mass",
    "xscale": "log",
    "yscale": "log",
    "xlim": (5.0, 100.0),
    "ylim": None,
    "xlabel": r"$m^{ee}$",
    "ylabel": r"$d\sigma^{ee}/dm^{ee}$",
    "units": {"x": r"GeV", "y": r"ub/GeV"},
    "label": r"CMS UPC dielectron invariant mass",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_mass,
}

obs_abs_costheta = {
    "tag": "abs_costheta",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$|\cos\theta^*|^{ee}$",
    "ylabel": r"$d\sigma^{ee}/d|\cos\theta^*|^{ee}$",
    "units": {"x": r"unit", "y": r"ub"},
    "label": r"CMS UPC dielectron Collins Soper angle",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_abs_costheta,
}

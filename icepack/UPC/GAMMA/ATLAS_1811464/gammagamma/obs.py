# ATLAS PbPb light-by-light observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

from core.io.cache import cache

from ....._common.obs import obs_cross_section as obs_cross_section
from .cuts import diphoton, selected_photons


# Compute the absolute diphoton rapidity
@cache
def proj_abs_rap(event):
    return abs(diphoton(event).rapidity)


# Compute the mean transverse momentum of the two photons
@cache
def proj_mean_pt(event):
    photons = selected_photons(event)
    return 0.5 * (photons[0]["pt"] + photons[1]["pt"])


# Compute the absolute two photon scattering angle proxy
@cache
def proj_abs_costheta(event):
    photons = selected_photons(event)
    return abs(math.tanh(0.5 * (photons[0]["eta"] - photons[1]["eta"])))


obs_abs_costheta = {
    "tag": "abs_costheta",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$|\cos\theta^*|$",
    "ylabel": r"$d\sigma/d|\cos\theta^*|$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"ATLAS UPC light-by-light scattering angle",
    "ylim_ratio": (0.0, 3.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_abs_costheta,
}

obs_mean_pt = {
    "tag": "mean_pt",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$(p_{T,1}^{\gamma}+p_{T,2}^{\gamma})/2$",
    "ylabel": r"$d\sigma/d\langle p_T^{\gamma}\rangle$",
    "units": {"x": r"GeV", "y": r"pb/GeV"},
    "label": r"ATLAS UPC light-by-light mean photon transverse momentum",
    "ylim_ratio": (0.0, 3.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_mean_pt,
}

obs_abs_rap = {
    "tag": "abs_rap",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$|y_{\gamma\gamma}|$",
    "ylabel": r"$d\sigma/d|y_{\gamma\gamma}|$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"ATLAS UPC light-by-light rapidity",
    "ylim_ratio": (0.0, 3.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_abs_rap,
}

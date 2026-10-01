# ATLAS PbPb exclusive dimuon observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

from core.io.cache import cache

from .cuts import dimuon, selected_muons


# Compute the absolute dimuon rapidity
@cache
def proj_abs_rap(event):
    return abs(dimuon(event).rapidity)

# Compute the dimuon invariant mass
@cache
def proj_mass(event):
    return dimuon(event).m

# Compute the ATLAS dimuon acoplanarity
@cache
def proj_acoplanarity(event):
    muons = selected_muons(event)
    value = 1.0 - muons[0]["p4"].abs_delta_phi(muons[1]["p4"]) / math.pi
    return min(1.0, max(0.0, value))

obs_mass = {
    "tag": "mass",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$m_{\mu\mu}$",
    "ylabel": r"$d\sigma/dm_{\mu\mu}$",
    "units": {"x": r"GeV", "y": r"pb"},
    "label": r"ATLAS UPC dimuon mass",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_mass,
}

obs_abs_rap = {
    "tag": "abs_rap",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$|y_{\mu\mu}|$",
    "ylabel": r"$d\sigma/d|y_{\mu\mu}|$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"ATLAS UPC dimuon rapidity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_abs_rap,
}

obs_acoplanarity = {
    "tag": "acoplanarity",
    "xscale": "log",
    "xlim": (1.0e-4, 0.2),
    "ylim": None,
    "xlabel": r"$1-|\Delta\phi_{\mu\mu}|/\pi$",
    "ylabel": r"$d\sigma/dA_{\phi}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"ATLAS UPC dimuon acoplanarity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_acoplanarity,
}

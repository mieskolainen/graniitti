# HERA vector meson decay angles in the helicity frame
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
from core.analysis import obs
from core.io.cache import cache

from .exclusive_jpsi import unique_momentum


# Project the first selected daughter with z opposite the target recoil and y along p cross q
# [REFERENCE: ZEUS, arXiv:hep-ex/0201043, Sections 3 and 8]
@cache
def daughter(event):
    proton = unique_momentum(event, 2212, True)
    photon = unique_momentum(event, 11, True) - unique_momentum(event, 11, False)
    records = obs.proj_central_particle_records(event)
    first = obs.flatten_pid_selection(event.pid)[0]
    selected = [record["p4"] for record in records if record["pid"] == first]
    if len(selected) != 1:
        raise ValueError("HERA helicity angles require one selected decay daughter")
    beams = obs.LorentFramePrepare(proton, photon, selected, obs.proj_central_system(event))
    return obs.LorentzFrame(*beams, frametype="HX")[0]


# Compute the helicity polar cosine of the selected decay daughter
def costheta(event):
    return daughter(event).costheta


# Compute the helicity azimuth on the measured interval from zero to two pi
def phi(event):
    return daughter(event).phi % (2.0 * np.pi)


obs_costheta_HX = {
    "tag": "costheta_HX", "xlim": None, "ylim": None,
    "xlabel": r"$\cos\theta_h$", "ylabel": r"$(1/N)\,dN/d\cos\theta_h$",
    "units": {"x": "", "y": ""}, "label": "Helicity polar angle",
    "ylim_ratio": (0.0, 2.0), "ytick_ratio_step": 0.5,
    "bins": None, "func": costheta,
}

obs_phi_HX = {
    "tag": "phi_HX", "xlim": None, "ylim": None,
    "xlabel": r"$\phi_h$", "ylabel": r"$(1/N)\,dN/d\phi_h$",
    "units": {"x": "rad", "y": r"rad$^{-1}$"}, "label": "Helicity azimuth",
    "ylim_ratio": (0.0, 2.0), "ytick_ratio_step": 0.5,
    "bins": None, "func": phi,
}

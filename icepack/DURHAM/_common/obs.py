# Event observables for Durham QCD SuperChic comparisons
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numpy as np
from core.analysis import obs
from core.analysis.observables.default import obs_M as default_M
from core.io.cache import cache


# Map every accepted event into one unit width cross section bin
def proj_cross_section(_event):
    return 0.5


# Compute the stable muon pair in charge order
@cache
def chic_muon_pair(event):
    records = [
        record
        for record in obs.proj_event_particles(event)
        if record["is_final"] and abs(record["pid"]) == 13
    ]
    positive = [record["p4"] for record in records if record["pid"] == -13]
    negative = [record["p4"] for record in records if record["pid"] == 13]
    if len(positive) != 1 or len(negative) != 1:
        return None
    return positive[0], negative[0]


# Compute the mu+ helicity angle used in SuperChic Figure 8a
@cache
def proj_costheta_jpsi_boost(event):
    pair = chic_muon_pair(event)
    if pair is None:
        return np.nan
    mu_plus, mu_minus = pair
    jpsi_lab = mu_plus + mu_minus
    mu_plus_rest = mu_plus.copy()
    mu_plus_rest.boost(b=jpsi_lab, sign=-1)
    denominator = mu_plus_rest.p3mod * jpsi_lab.p3mod
    if denominator <= 0.0:
        return np.nan
    value = np.dot(mu_plus_rest.p3, jpsi_lab.p3) / denominator
    return float(np.clip(value, -1.0, 1.0))


obs_cross_section = {
    "tag": "cross_section",
    "xlim": (0.0, 1.0),
    "ylim": None,
    "xlabel": r"Fiducial selection",
    "ylabel": r"$\sigma$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Fiducial cross section",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.array([0.0, 1.0]),
    "density": False,
    "func": proj_cross_section,
}

obs_M = copy.deepcopy(default_M)
obs_M["xlim"] = (3.5, 200.0)
obs_M["bins"] = np.linspace(3.5, 200.0, 50)
obs_M["xlabel"] = r"$M_X$"

obs_M_figure3 = copy.deepcopy(default_M)
obs_M_figure3["tag"] = "M_figure3"
obs_M_figure3["xlim"] = (40.0, 200.0)
obs_M_figure3["bins"] = np.linspace(40.0, 200.0, 31)
obs_M_figure3["xlabel"] = r"$M_X$"

obs_costheta_jpsi_boost = {
    "tag": "costheta_jpsi_boost",
    "xlim": (-1.0, 1.0),
    "ylim": None,
    "xlabel": r"$\cos\theta$",
    "ylabel": r"$1/\sigma\,d\sigma/d\cos\theta$",
    "units": {"x": r"unit", "y": r"unit"},
    "label": r"$\mu^+$ angle in the $J/\psi$ rest frame",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(-1.0, 1.0, 21),
    "density": False,
    "func": proj_costheta_jpsi_boost,
}

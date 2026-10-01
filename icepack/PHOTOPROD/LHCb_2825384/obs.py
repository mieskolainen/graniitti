# LHCb 13 TeV exclusive charmonium histogram definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
from core.analysis import obs


# Project every selected event into the single integrated cross-section bin
def proj_XS(event):
    return 0.5

obs_Rap = {
    "tag": "Rap",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$Y$",
    "ylabel": r"$d\sigma/dY$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Charmonium rapidity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,  # Supplied by the HEPData reader
    "density": False,
    "func": obs.proj_1D_Rap,
}

obs_XS = {
    "tag": "XS",
    "xlim": (0.0, 1.0),
    "ylim": None,
    "xlabel": r"Fiducial region",
    "ylabel": r"$\sigma$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Integrated fiducial cross section",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.asarray([0.0, 1.0]),
    "density": False,
    "func": proj_XS,
}

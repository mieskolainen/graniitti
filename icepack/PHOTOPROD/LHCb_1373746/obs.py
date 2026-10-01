# LHCb shared exclusive Upsilon histogram definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

obs_Rap = {
    "tag": "Rap",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$Y$",
    "ylabel": r"$d\sigma/dY$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Upsilon rapidity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,  # Supplied by the HEPData reader
    "density": False,
    "func": obs.proj_1D_Rap,
}

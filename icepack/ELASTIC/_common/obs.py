# Shared elastic histogram definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

obs_mandelstam_t = {
    "tag": "mandelstam_t",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$|t|$",
    "ylabel": r"$d\sigma/d|t|$",
    "units": {"x": r"GeV$^2$", "y": r"pb"},
    "label": r"Mandelstam $|t|$",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": obs.proj_1D_mandelstam_t,
}

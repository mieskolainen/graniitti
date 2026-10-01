# Shared exclusive HERA photoproduction observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from . import exclusive_jpsi

obs_abs_t = {
    "tag": "abs_t",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$|t|$",
    "ylabel": r"$d\sigma/d|t|$",
    "units": {"x": r"GeV$^2$", "y": r"pb GeV$^{-2}$"},
    "label": r"Proton momentum transfer",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,  # Supplied by the HEPData reader
    "density": False,
    "func": exclusive_jpsi.abs_t,
}

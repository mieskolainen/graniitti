# HERA proton dissociation momentum transfer and energy observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from . import dissociative

obs_abs_t = {
    "tag": "abs_t",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$|t|$",
    "ylabel": r"$d\sigma(\gamma p\to VY)/d|t|$",
    "units": {"x": r"GeV$^2$", "y": r"pb GeV$^{-2}$"},
    "label": "Proton momentum transfer",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,  # Supplied by the HEPData reader
    "density": False,
    "func": dissociative.abs_t,
}

obs_w = {
    "tag": "w",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$W_{\gamma p}$",
    "ylabel": r"$\sigma(\gamma p\to VY)$",
    "units": {"x": "GeV", "y": "pb"},
    "label": "Photon-proton energy",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,  # Supplied by the HEPData reader
    "density": False,
    "func": dissociative.w,
}

# ALICE incoherent J/psi momentum transfer observable
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache


# Compute the measured |t| estimator from the dimuon transverse momentum
@cache
def proj_abs_t(event):
    # [REFERENCE: arXiv:2305.06169, event selection uses |t| = pT^2 without photon pT unfolding]
    pt = obs.proj_1D_Pt(event)
    return pt * pt


obs_abs_t = {
    "tag": "abs_t",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$|t|$",
    "ylabel": r"$d\sigma_{\gamma\mathrm{Pb}}/d|t|$",
    "units": {"x": r"GeV$^2$", "y": r"pb"},
    "label": r"ALICE incoherent $J/\psi$",
    "ylim_ratio": (0.0, 2.5),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_abs_t,
}

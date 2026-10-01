# CMS reconstructed rho momentum transfer observable
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache


# Compute the published approximation |t| = pT(rho)^2
@cache
def proj_pt2(event):
    return obs.proj_1D_Pt(event)**2


obs_pt2 = {
    "tag": "pt2", "xlim": None, "ylim": None, "bins": None, "density": False,
    "xlabel": r"$|t|\simeq p_{T,\rho}^{2}$", "ylabel": r"$d\sigma_{\gamma p}/d|t|$",
    "units": {"x": r"GeV$^2$", "y": "ub"}, "label": r"CMS exclusive $\rho^0$",
    "ylim_ratio": (0.0, 2.0), "ytick_ratio_step": 0.5, "func": proj_pt2,
}

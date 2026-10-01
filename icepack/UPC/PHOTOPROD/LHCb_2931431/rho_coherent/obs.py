# LHCb signed rho rapidity observable
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

obs_rap = {
    "tag": "rap", "xlim": None, "ylim": None, "bins": None, "density": False,
    "xlabel": r"$y_{\rho^0}$", "ylabel": r"$d\sigma/dy$",
    "units": {"x": "unit", "y": "mb"}, "label": r"LHCb coherent $\rho^0$",
    "ylim_ratio": (0.0, 2.0), "ytick_ratio_step": 0.5, "func": obs.proj_1D_Rap,
}

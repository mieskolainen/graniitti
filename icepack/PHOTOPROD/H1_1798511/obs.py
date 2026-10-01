# H1 elastic pion pair photoproduction histogram definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

obs_M = {
    "tag": "M",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$M_{\pi^+\pi^-}$",
    "ylabel": r"$d\sigma/dM_{\pi^+\pi^-}$",
    "units": {"x": r"GeV", "y": r"pb"},
    "label": r"Elastic pion pair invariant mass",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,  # Supplied by the HEPData reader
    "density": False,
    "func": obs.proj_1D_M,
}

# Compute the dipion mass density per unit proton momentum transfer
obs_M_t = {**obs_M, 'tag':'M_t', 'ylabel':r'$d^2\sigma/(dM_{\pi\pi}\,d|t|)$',
           'units':{'x':'GeV', 'y':'pb', 'yden':r'GeV$^3$'}}

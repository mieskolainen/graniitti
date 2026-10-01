# CMS 13 TeV pion-pair histogram definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

obs_fp1pt_fp2pt_M = {
    "tag": "fp1pt_fp2pt_M",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$M_{\pi^+\pi^-}$",
    "ylabel": r"$d^3\sigma/(dM\,dp_{1T}\,dp_{2T})$",
    "units": {"x": r"GeV", "y": r"pb"},
    "category_labels": (r"$p_{1T}$", r"$p_{2T}$"),
    "category_units": (r"GeV", r"GeV"),
    "label": r"Pion-pair invariant mass",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "symmetrized_fill": True,
    "func": obs.proj_3D_fp1pt_fp2pt_M,
}

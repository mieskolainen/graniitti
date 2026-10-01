# Default histogram definitions for the observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np

from core.analysis import obs

obs_M = {
    "tag": "M",
    "xlim": (0, 3.0),
    "ylim": None,
    "xlabel": r"$M$",
    "ylabel": r"$d\sigma/dM$",
    "units": {"x": r"GeV", "y": r"pb"},
    "label": r"Invariant mass",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0, 3, 50),
    "density": False,
    "func": obs.proj_1D_M,
}


obs_Abs_t1t2 = {
    "tag": "Abs_t1t2",
    "xlim": (0, 2),
    "ylim": None,
    "xlabel": r"$|t_1+t_2|$",
    "ylabel": r"$d\sigma/d|t_1+t_2|$",
    "units": {"x": r"GeV$^{2}$", "y": r"pb"},
    "label": r"Sum of momentum transfers squared",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0, 2, 50),
    "density": False,
    "func": obs.proj_1D_Abs_t1t2,
}


obs_Rap = {
    "tag": "Rap",
    "xlim": (-4, 4),
    "ylim": None,
    "xlabel": r"$Y$",
    "ylabel": r"$d\sigma/dY$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Rapidity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(-3, 3, 50),
    "density": False,
    "func": obs.proj_1D_Rap,
}


obs_dPhi_pp = {
    "tag": "dPhi_pp",
    "xlim": (0, 180),
    "ylim": None,
    "xlabel": r"$\Delta\phi_{pp}$",
    "ylabel": r"$d\sigma/d\Delta\phi_{pp}$",
    "units": {"x": r"deg", "y": r"pb"},
    "label": r"Transverse angle separation",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0, 180, 50),
    "density": False,
    "func": obs.proj_1D_dPhi_pp,
}


obs_costheta_CS = {
    "tag": "costheta_CS",
    "xlim": (-1, 1),
    "ylim": None,
    "xlabel": r"cos$(\theta)_{CS}$",
    "ylabel": r"$d\sigma/dcos(\theta)_{CS}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Angle in Collins-Soper frame",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(-1, 1, 50),
    "density": False,
    "func": obs.proj_1D_costheta_CS,
}


obs_phi_CS = {
    "tag": "phi_CS",
    "xlim": (-180, 180),
    "ylim": None,
    "xlabel": r"$\phi_{CS}$",
    "ylabel": r"$d\sigma/d\phi_{CS}$",
    "units": {"x": r"deg", "y": r"pb"},
    "label": r"Angle in Collins-Soper frame",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(-180, 180, 50),
    "density": False,
    "func": obs.proj_1D_phi_CS,
}


obs_mandelstam_t = {
    "tag": "mandelstam_t",
    "xlim": (0, 2.5),
    "ylim": None,
    "xlabel": r"$|t|$",
    "ylabel": r"$d \sigma/d|t|$",
    "units": {"x": r"GeV$^2$", "y": r"pb"},
    "label": r"Mandelstam $|t|$",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0, 2.5, 50),
    "density": False,
    "func": obs.proj_1D_mandelstam_t,
}


obs_fp1pt_fp2pt_dPhi_pp = {
    "tag": "fp1pt_fp2pt_dPhi_pp",
    "xlim": {"cat0": (0.0, 3.14159)},  # (dict of categories)
    "ylim": None,
    "xlabel": r"$\Delta\phi_{pp}$",
    "ylabel": r"$d^3 \sigma/d\Delta\phi_{pp} dp_{p_{1T}} dp_{p_{2T}}$",
    "units": {"x": r"rad", "y": r"pb"},
    "category_labels": (r"$p_{1T}$", r"$p_{2T}$"),
    "category_units": (r"GeV", r"GeV"),
    "label": r"Transverse angle separation",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": {"cat0": {0: np.array([0.0, 0.8]), 1: np.array([0.0, 0.8]), 2: np.linspace(0.0, 3.14159, 50)}},
    "density": False,
    "symmetrized_fill": True,  # Fill both MC events pt1 < pt2 & pt1 > pt2 into the same plot
    "func": obs.proj_3D_fp1pt_fp2pt_dPhi_pp,
}

obs_fp1pt_fp2pt_M = {
    "tag": "fp1pt_fp2pt_M",
    "xlim": {"cat0": (0.0, 2.5)},  # (dict of categories)
    "ylim": None,
    "xlabel": r"$M$",
    "ylabel": r"$d^3 \sigma/dM dp_{p_{1T}} dp_{p_{2T}}$",
    "units": {"x": r"GeV", "y": r"pb"},
    "category_labels": (r"$p_{1T}$", r"$p_{2T}$"),
    "category_units": (r"GeV", r"GeV"),
    "label": r"Invariant mass",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": {"cat0": {0: np.array([0.0, 0.8]), 1: np.array([0.0, 0.8]), 2: np.linspace(0.0, 2.5, 50)}},
    "density": False,
    "symmetrized_fill": True,  # Fill both MC events pt1 < pt2 & pt1 > pt2 into the same plot
    "func": obs.proj_3D_fp1pt_fp2pt_M,
}

obs_fp1pt_fp2pt_max_hat_tu = {
    "tag": "fp1pt_fp2pt_max_hat_tu",
    "xlim": {"cat0": (-2.5, 0.0)},  # (dict of categories)
    "ylim": None,
    "xlabel": r"$\max(\hat{t}, \hat{u})$",
    "ylabel": r"$d^3 \sigma/d\max(\hat{t}, \hat{u}) dp_{p_{1T}} dp_{p_{2T}} $",
    "units": {"x": r"GeV$^2$", "y": r"pb"},
    "category_labels": (r"$p_{1T}$", r"$p_{2T}$"),
    "category_units": (r"GeV", r"GeV"),
    "label": r"Sub-Mandelstam",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": {"cat0": {0: np.array([0.0, 0.8]), 1: np.array([0.0, 0.8]), 2: np.linspace(-2.5, 0.0, 50)}},
    "density": False,
    "symmetrized_fill": True,  # Fill both MC events pt1 < pt2 & pt1 > pt2 into the same plot
    "func": obs.proj_3D_fp1pt_fp2pt_max_hat_tu,
}

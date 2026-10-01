# Shared STAR 200 GeV central-production histogram factories
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

from . import projectors


# Compute the STAR central-system invariant-mass definition
def invariant_mass():
    return {
        "tag": "M",
        "xlim": None,
        "ylim": None,
        "xlabel": r"$M$",
        "ylabel": r"$d\sigma/dM$",
        "units": {"x": r"GeV", "y": r"pb"},
        "label": r"Central-system invariant mass",
        "ylim_ratio": (0.0, 2.0),
        "ytick_ratio_step": 0.5,
        "bins": None,
        "density": False,
        "func": projectors.invariant_mass,
    }


# Compute the STAR summed momentum-transfer definition
def absolute_momentum_transfer_sum():
    return {
        "tag": "Abs_t1t2",
        "xlim": None,
        "ylim": None,
        "xlabel": r"$|t_1+t_2|$",
        "ylabel": r"$d\sigma/d|t_1+t_2|$",
        "units": {"x": r"GeV$^2$", "y": r"pb"},
        "label": r"Sum of momentum transfers",
        "ylim_ratio": (0.0, 2.0),
        "ytick_ratio_step": 0.5,
        "bins": None,
        "density": False,
        "func": obs.proj_1D_Abs_t1t2,
    }


# Compute the STAR central-system rapidity definition
def rapidity():
    return {
        "tag": "Rap",
        "xlim": None,
        "ylim": None,
        "xlabel": r"$Y$",
        "ylabel": r"$d\sigma/dY$",
        "units": {"x": r"unit", "y": r"pb"},
        "label": r"Central-system rapidity",
        "ylim_ratio": (0.0, 2.0),
        "ytick_ratio_step": 0.5,
        "bins": None,
        "density": False,
        "func": projectors.rapidity,
    }


# Compute the STAR outgoing-proton azimuthal-separation definition
def proton_azimuthal_separation():
    return {
        "tag": "dPhi_pp",
        "xlim": None,
        "ylim": None,
        "xlabel": r"$\Delta\phi_{pp}$",
        "ylabel": r"$d\sigma/d\Delta\phi_{pp}$",
        "units": {"x": r"deg", "y": r"pb"},
        "label": r"Outgoing-proton azimuthal separation",
        "ylim_ratio": (0.0, 2.0),
        "ytick_ratio_step": 0.5,
        "bins": None,
        "density": False,
        "func": obs.proj_1D_dPhi_pp,
    }


# Compute the STAR Collins-Soper polar-angle definition
def collins_soper_polar_angle():
    return {
        "tag": "costheta_CS",
        "xlim": None,
        "ylim": None,
        "xlabel": r"$\cos\theta_{\mathrm{CS}}$",
        "ylabel": r"$d\sigma/d\cos\theta_{\mathrm{CS}}$",
        "units": {"x": r"unit", "y": r"pb"},
        "label": r"Collins-Soper polar angle",
        "ylim_ratio": (0.0, 2.0),
        "ytick_ratio_step": 0.5,
        "bins": None,
        "density": False,
        "func": projectors.collins_soper_costheta,
    }


# Compute the STAR Collins-Soper azimuthal-angle definition
def collins_soper_azimuthal_angle():
    return {
        "tag": "phi_CS",
        "xlim": None,
        "ylim": None,
        "xlabel": r"$\phi_{\mathrm{CS}}$",
        "ylabel": r"$d\sigma/d\phi_{\mathrm{CS}}$",
        "units": {"x": r"deg", "y": r"pb"},
        "label": r"Collins-Soper azimuthal angle",
        "ylim_ratio": (0.0, 2.0),
        "ytick_ratio_step": 0.5,
        "bins": None,
        "density": False,
        "func": projectors.collins_soper_phi,
    }

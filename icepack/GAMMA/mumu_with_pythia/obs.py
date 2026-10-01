# Forward excitation dimuon and charged particle observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
#
# [REFERENCE: ATLAS, Phys. Lett. B 749 (2015) 242, arXiv:1506.07098]

import math

import numba
import numpy as np
from core.analysis import obs
from core.io.cache import cache

from . import cuts as GAMMA_CUTS

CHARGED_ABS_PDGS = np.asarray(
    [11, 13, 15, 211, 321, 2212, 3112, 3222, 3312, 3334],
    dtype=np.int64,
)


# Test whether one stable particle PDG code carries electric charge
@numba.njit
def is_charged_pid(pid):
    return np.any(abs(pid) == CHARGED_ABS_PDGS)


# Count stable charged particles in one absolute pseudorapidity interval
@numba.njit
def count_charged(pid, eta, pt, is_final, eta_min, eta_max, pt_min):
    count = 0
    for index in range(pid.size):
        abs_eta = abs(eta[index])
        if (
            is_final[index]
            and is_charged_pid(pid[index])
            and eta_min < abs_eta < eta_max
            and pt[index] > pt_min
        ):
            count += 1
    return count


# Compute compact arrays for stable charged particle observables
@cache
def particle_arrays(event):
    records = obs.proj_event_particles(event)
    return (
        np.asarray([record["pid"] for record in records], dtype=np.int64),
        np.asarray([record["eta"] for record in records], dtype=np.float64),
        np.asarray([record["pt"] for record in records], dtype=np.float64),
        np.asarray([record["is_final"] for record in records], dtype=np.bool_),
    )


# Compute the selected dimuon four-vector pair
def selected_pair(event):
    first, second = GAMMA_CUTS.selected_muons(event)
    return first["p4"], second["p4"]


# Compute the selected dimuon invariant mass
def proj_m_mumu(event):
    first, second = selected_pair(event)
    return (first + second).m


# Compute the selected dimuon transverse momentum
def proj_pt_mumu(event):
    first, second = selected_pair(event)
    return (first + second).pt


# Compute the selected dimuon absolute rapidity
def proj_abs_y_mumu(event):
    first, second = selected_pair(event)
    return abs((first + second).rapidity)


# Compute the selected dimuon acoplanarity
def proj_acoplanarity(event):
    first, second = selected_pair(event)
    value = 1.0 - first.abs_delta_phi(second) / math.pi
    return min(1.0, max(0.0, value))


# Compute central stable charged multiplicity after removing the selected muons
def proj_n_charged_central(event):
    pid, eta, pt, is_final = particle_arrays(event)
    value = count_charged(pid, eta, pt, is_final, 0.0, 2.5, 0.1)
    return max(value - 2, 0)


# Compute forward stable charged multiplicity from the dissociative system
def proj_n_charged_forward(event):
    return count_charged(*particle_arrays(event), 2.5, 5.0, 0.1)


# Construct one standard cross section histogram definition
def histogram(tag, xlim, xlabel, ylabel, units, label, bins, function):
    return {
        "tag": tag,
        "xlim": xlim,
        "ylim": None,
        "xlabel": xlabel,
        "ylabel": ylabel,
        "units": {"x": units, "y": r"pb"},
        "label": label,
        "ylim_ratio": (0.0, 2.0),
        "ytick_ratio_step": 0.5,
        "bins": bins,
        "density": False,
        "func": function,
    }


obs_m_mumu = histogram(
    "m_mumu",
    (20.0, 200.0),
    r"$m_{\mu\mu}$",
    r"$d\sigma/dm_{\mu\mu}$",
    r"GeV",
    r"Dimuon invariant mass",
    np.linspace(20.0, 200.0, 37),
    proj_m_mumu,
)
obs_pt_mumu = histogram(
    "pt_mumu",
    (0.0, 1.5),
    r"$p_T^{\mu\mu}$",
    r"$d\sigma/dp_T^{\mu\mu}$",
    r"GeV",
    r"Dimuon transverse momentum",
    np.linspace(0.0, 1.5, 31),
    proj_pt_mumu,
)
obs_abs_y_mumu = histogram(
    "abs_y_mumu",
    (0.0, 2.5),
    r"$|y_{\mu\mu}|$",
    r"$d\sigma/d|y_{\mu\mu}|$",
    r"unit",
    r"Dimuon absolute rapidity",
    np.linspace(0.0, 2.5, 26),
    proj_abs_y_mumu,
)
obs_acoplanarity = histogram(
    "acoplanarity",
    (0.0, 0.05),
    r"$1-|\Delta\phi_{\mu\mu}|/\pi$",
    r"$d\sigma/dA_\phi$",
    r"unit",
    r"Dimuon acoplanarity",
    np.linspace(0.0, 0.05, 26),
    proj_acoplanarity,
)
obs_n_charged_central = histogram(
    "n_charged_central",
    (0.0, 30.0),
    r"$N_{\mathrm{ch}}(|\eta|<2.5)$ excluding muons",
    r"$d\sigma/dN_{\mathrm{ch}}$",
    r"unit",
    r"Central stable charged multiplicity",
    np.arange(-0.5, 30.5, 1.0),
    proj_n_charged_central,
)
obs_n_charged_forward = histogram(
    "n_charged_forward",
    (0.0, 50.0),
    r"$N_{\mathrm{ch}}(2.5<|\eta|<5.0)$",
    r"$d\sigma/dN_{\mathrm{ch}}$",
    r"unit",
    r"Forward stable charged multiplicity",
    np.arange(-0.5, 50.5, 1.0),
    proj_n_charged_forward,
)

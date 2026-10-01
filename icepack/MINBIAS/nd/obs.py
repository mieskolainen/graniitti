# Inclusive charged particle observables for nondiffractive events
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numba
import numpy as np
from core.analysis import obs

CHARGED_ABS_PDGS = np.asarray(
    [11, 13, 15, 211, 321, 2212, 3112, 3222, 3312, 3334],
    dtype=np.int64,
)


# Test whether one stable-particle PDG code carries electric charge
@numba.njit
def is_charged_pid(pid):
    return np.any(abs(pid) == CHARGED_ABS_PDGS)


# Count stable charged particles inside one pseudorapidity interval
@numba.njit
def count_charged_in_eta(pid, eta, is_final, eta_max):
    count = 0
    for index in range(pid.size):
        if is_final[index] and abs(eta[index]) < eta_max:
            count += int(is_charged_pid(pid[index]))
    return count


# Convert cached HepMC particle records into compact Numba arrays
def particle_arrays(event):
    records = obs.proj_event_particles(event)
    return (
        np.asarray([record["pid"] for record in records], dtype=np.int64),
        np.asarray([record["eta"] for record in records], dtype=np.float64),
        np.asarray([record["is_final"] for record in records], dtype=np.bool_),
    )


# Compute the stable charged multiplicity in |eta| below 2.5
def proj_n_charged(event):
    pid, eta, is_final = particle_arrays(event)
    return count_charged_in_eta(pid, eta, is_final, 2.5)


# Compute the stable charged multiplicity in |eta| below 0.5
def proj_n_charged_midrapidity(event):
    pid, eta, is_final = particle_arrays(event)
    return count_charged_in_eta(pid, eta, is_final, 0.5)


obs_n_charged = {
    "tag": "n_charged",
    "xlim": (0.0, 100.0),
    "ylim": None,
    "xlabel": r"$N_{\mathrm{ch}}(|\eta|<2.5)$",
    "ylabel": r"$1/N_{\mathrm{evt}}\,dN_{\mathrm{evt}}/dN_{\mathrm{ch}}$",
    "units": {"x": r"unit", "y": r"unit"},
    "label": r"Stable charged multiplicity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.arange(-0.5, 100.5, 1.0),
    "density": True,
    "func": proj_n_charged,
}

obs_n_charged_midrapidity = {
    "tag": "n_charged_midrapidity",
    "xlim": (0.0, 40.0),
    "ylim": None,
    "xlabel": r"$N_{\mathrm{ch}}(|\eta|<0.5)$",
    "ylabel": r"$1/N_{\mathrm{evt}}\,dN_{\mathrm{evt}}/dN_{\mathrm{ch}}$",
    "units": {"x": r"unit", "y": r"unit"},
    "label": r"Midrapidity stable charged multiplicity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.arange(-0.5, 40.5, 1.0),
    "density": True,
    "func": proj_n_charged_midrapidity,
}

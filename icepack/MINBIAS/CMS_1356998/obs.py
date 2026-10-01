# CMS 7 TeV diffraction observable definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import numpy as np
from core.analysis import obs
from core.io.cache import cache

PROTON_MASS_GEV = 0.9382720813
UNDERFLOW = -100.0


# Reconstruct X and Y from stable particles on either side of the largest rapidity gap
@cache
def diffractive_systems(event):
    # [REFERENCE: CMS Phys. Rev. D 92 (2015) 012003, Sec. 3, arXiv:1503.08689]
    particles = sorted([
        record["p4"]
        for record in obs.proj_event_particles(event)
        if record["is_final"]
    ], key=lambda momentum: momentum.rapidity)
    if len(particles) < 2:
        return ()
    split = int(np.argmax(np.diff([momentum.rapidity for momentum in particles]))) + 1
    return sum(particles[split:]), sum(particles[:split])


# Compute the exact squared collision energy from the two beam records
@cache
def collision_s(event):
    beams = [record["p4"] for record in obs.proj_event_particles(event) if record["status"] == 4]
    if len(beams) != 2:
        raise ValueError("cms_diffraction: expected exactly two beam particles")
    return (beams[0] + beams[1]).m2


# Convert one dissociation-system mass to log10 xi
def system_log10_xi(momentum, s):
    if momentum.m2 <= 0.0 or s <= 0.0:
        raise ValueError("cms_diffraction: nonpositive diffractive mass or s")
    return math.log10(momentum.m2 / s)


# Compute log10 xi_X in the SD dominated region, including low mass DD
def proj_log10_xi_sd(event):
    systems = diffractive_systems(event)
    # [REFERENCE: CMS Phys. Rev. D 92 (2015) 012003, Table 1 FG2 selection]
    if len(systems) != 2 or systems[1].m <= 0.0 or math.log10(systems[1].m) >= 0.5:
        return UNDERFLOW
    value = system_log10_xi(systems[0], collision_s(event))
    return value if -5.5 < value < -2.5 else UNDERFLOW


# Compute log10 xi_X in the CMS CASTOR double-dissociation region
def proj_log10_xi_dd_forward(event):
    systems = diffractive_systems(event)
    if len(systems) != 2:
        return UNDERFLOW
    # [REFERENCE: CMS Phys. Rev. D 92 (2015) 012003, Table 1 FG2 selection]
    system_x, system_y = systems
    s = collision_s(event)
    log_xi_x = system_log10_xi(system_x, s)
    log_m_y = math.log10(system_y.m)
    if -5.5 < log_xi_x < -2.5 and 0.5 < log_m_y < 1.1:
        return log_xi_x
    return UNDERFLOW


# Compute the CMS rapidity-gap estimator for one double-dissociation event
def dd_delta_eta(event):
    systems = diffractive_systems(event)
    if len(systems) != 2:
        return UNDERFLOW
    masses = [system.m for system in systems]
    xi_product = masses[0] ** 2 * masses[1] ** 2 / (collision_s(event) * PROTON_MASS_GEV**2)
    if xi_product <= 0.0:
        raise ValueError("cms_diffraction: nonpositive DD xi product")
    return -math.log(xi_product)


# Compute the CMS central rapidity-gap estimator for double dissociation
def proj_delta_eta_dd_central(event):
    systems = diffractive_systems(event)
    if len(systems) != 2:
        return UNDERFLOW
    if any(math.log10(system.m) <= 1.1 for system in systems):
        return UNDERFLOW
    delta_eta = dd_delta_eta(event)
    return delta_eta if delta_eta > 3.0 else UNDERFLOW


# Compute every generated DD event in the extrapolated CMS gap region
def proj_delta_eta_dd_all(event):
    delta_eta = dd_delta_eta(event)
    return delta_eta if delta_eta > 3.0 else UNDERFLOW


obs_log10_xi_sd = {
    "tag": "log10_xi_sd",
    "xlim": (-5.5, -2.5),
    "ylim": None,
    "xlabel": r"$\log_{10}\xi_X$",
    "ylabel": r"$d\sigma/d\log_{10}\xi_X$",
    "units": {"x": "unit", "y": "mb"},
    "label": "SD dominated mass fraction",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.arange(-5.5, -2.0, 0.5),
    "density": False,
    "func": proj_log10_xi_sd,
}

obs_log10_xi_dd_forward = {
    "tag": "log10_xi_dd_forward",
    "xlim": (-5.5, -2.5),
    "ylim": None,
    "xlabel": r"$\log_{10}\xi_X$",
    "ylabel": r"$d\sigma/d\log_{10}\xi_X$",
    "units": {"x": "unit", "y": "mb"},
    "label": "DD dominated mass fraction",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.arange(-5.5, -2.0, 0.5),
    "density": False,
    "func": proj_log10_xi_dd_forward,
}

obs_delta_eta_dd_central = {
    "tag": "delta_eta_dd_central",
    "xlim": (3.0, 7.5),
    "ylim": None,
    "xlabel": r"$\Delta\eta$",
    "ylabel": r"$d\sigma/d\Delta\eta$",
    "units": {"x": "unit", "y": "mb"},
    "label": "Double-dissociation central rapidity gap",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.asarray([3.0, 4.5, 6.0, 7.5]),
    "density": False,
    "func": proj_delta_eta_dd_central,
}

obs_delta_eta_dd_all = {
    "tag": "delta_eta_dd_all",
    "xlim": (3.0, 40.0),
    "ylim": None,
    "xlabel": r"$\Delta\eta$",
    "ylabel": r"$d\sigma/d\Delta\eta$",
    "units": {"x": "unit", "y": "mb"},
    "label": "Double-dissociation rapidity gap",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.asarray([3.0, 40.0]),
    "density": False,
    "func": proj_delta_eta_dd_all,
}

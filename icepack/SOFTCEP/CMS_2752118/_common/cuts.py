# Shared CMS pion-pair fiducial and proton angular selections
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numba
import numpy as np
from core.analysis import obs

cut_param = {}


# Apply one strict open interval
@numba.njit
def rcut(x, r):
    return r[0] < x and x < r[1]


# Apply one strict histogram interval expressed in degrees
@numba.njit
def degree_rcut(x, lower, upper):
    degree = np.pi / 180.0
    return rcut(x, (lower * degree, upper * degree))


# Apply the shared published CMS pion and proton fiducial selection
def fiducial_cuts(event, mass_range):
    records = obs.proj_event_particles(event)
    pions = [record for record in records if record["is_final"] and abs(record["pid"]) == 211]

    if sum(record["pid"] == 211 for record in pions) != 1:
        return False
    if sum(record["pid"] == -211 for record in pions) != 1:
        return False
    if not all(rcut(record["p4"].rapidity, (-2.0, 2.0)) for record in pions):
        return False

    _, final_protons = obs.proj_init_final_protons(event)
    if not all(rcut(proton.pt, (0.2, 0.8)) for proton in final_protons):
        return False

    mass = sum(record["p4"] for record in pions).m
    return bool(rcut(mass, mass_range))


# Reject the missing bins of the published Delta phi distributions
@numba.njit
def dPhi_cuts(pt1: float, pt2: float, dphi: float):
    """Return false in angular cells missing from CMS HEPData Table 1"""

    cut = np.zeros(41, dtype=np.bool_)

    cut[0] = not (rcut(pt1, (0.2, 0.25)) and rcut(pt2, (0.2, 0.25)) and degree_rcut(dphi, 70, 130))
    cut[1] = not (rcut(pt1, (0.25, 0.3)) and rcut(pt2, (0.2, 0.25)) and degree_rcut(dphi, 80, 110))
    cut[2] = not (rcut(pt1, (0.2, 0.25)) and rcut(pt2, (0.25, 0.3)) and degree_rcut(dphi, 80, 110))
    cut[3] = not (rcut(pt1, (0.25, 0.3)) and rcut(pt2, (0.25, 0.3)) and degree_rcut(dphi, 90, 100))
    cut[4] = not (rcut(pt1, (0.3, 0.35)) and rcut(pt2, (0.2, 0.25)) and degree_rcut(dphi, 90, 100))
    cut[5] = not (rcut(pt1, (0.2, 0.25)) and rcut(pt2, (0.3, 0.35)) and degree_rcut(dphi, 90, 100))
    cut[6] = not (rcut(pt1, (0.35, 0.4)) and rcut(pt2, (0.2, 0.25)) and degree_rcut(dphi, 90, 100))
    cut[7] = not (rcut(pt1, (0.2, 0.25)) and rcut(pt2, (0.35, 0.4)) and degree_rcut(dphi, 90, 100))
    cut[8] = not (rcut(pt1, (0.6, 0.65)) and rcut(pt2, (0.6, 0.65)) and degree_rcut(dphi, 60, 70))
    cut[9] = not (rcut(pt1, (0.65, 0.7)) and rcut(pt2, (0.65, 0.7)) and degree_rcut(dphi, 60, 80))
    cut[10] = not (rcut(pt1, (0.7, 0.75)) and rcut(pt2, (0.45, 0.5)) and degree_rcut(dphi, 60, 70))
    cut[11] = not (rcut(pt1, (0.45, 0.5)) and rcut(pt2, (0.7, 0.75)) and degree_rcut(dphi, 60, 70))
    cut[12] = not (rcut(pt1, (0.7, 0.75)) and rcut(pt2, (0.5, 0.55)) and degree_rcut(dphi, 60, 70))
    cut[13] = not (rcut(pt1, (0.5, 0.55)) and rcut(pt2, (0.7, 0.75)) and degree_rcut(dphi, 60, 70))
    cut[14] = not (rcut(pt1, (0.7, 0.75)) and rcut(pt2, (0.55, 0.6)) and degree_rcut(dphi, 60, 70))
    cut[15] = not (rcut(pt1, (0.55, 0.6)) and rcut(pt2, (0.7, 0.75)) and degree_rcut(dphi, 60, 70))
    cut[16] = not (rcut(pt1, (0.7, 0.75)) and rcut(pt2, (0.6, 0.65)) and degree_rcut(dphi, 60, 70))
    cut[17] = not (rcut(pt1, (0.6, 0.65)) and rcut(pt2, (0.7, 0.75)) and degree_rcut(dphi, 60, 70))
    cut[18] = not (rcut(pt1, (0.7, 0.75)) and rcut(pt2, (0.65, 0.7)) and degree_rcut(dphi, 70, 80))
    cut[19] = not (rcut(pt1, (0.65, 0.7)) and rcut(pt2, (0.7, 0.75)) and degree_rcut(dphi, 70, 80))
    cut[20] = not (rcut(pt1, (0.7, 0.75)) and rcut(pt2, (0.7, 0.75)) and degree_rcut(dphi, 50, 80))
    cut[21] = not (
        rcut(pt1, (0.75, 0.8)) and rcut(pt2, (0.2, 0.25)) and degree_rcut(dphi, 170, 180)
    )
    cut[22] = not (
        rcut(pt1, (0.2, 0.25)) and rcut(pt2, (0.75, 0.8)) and degree_rcut(dphi, 170, 180)
    )
    cut[23] = not (rcut(pt1, (0.75, 0.8)) and rcut(pt2, (0.35, 0.4)) and degree_rcut(dphi, 60, 70))
    cut[24] = not (rcut(pt1, (0.35, 0.4)) and rcut(pt2, (0.75, 0.8)) and degree_rcut(dphi, 60, 70))
    cut[25] = not (rcut(pt1, (0.75, 0.8)) and rcut(pt2, (0.4, 0.45)) and degree_rcut(dphi, 50, 60))
    cut[26] = not (rcut(pt1, (0.4, 0.45)) and rcut(pt2, (0.75, 0.8)) and degree_rcut(dphi, 50, 60))
    cut[27] = not (rcut(pt1, (0.75, 0.8)) and rcut(pt2, (0.4, 0.45)) and degree_rcut(dphi, 70, 80))
    cut[28] = not (rcut(pt1, (0.4, 0.45)) and rcut(pt2, (0.75, 0.8)) and degree_rcut(dphi, 70, 80))
    cut[29] = not (rcut(pt1, (0.75, 0.8)) and rcut(pt2, (0.45, 0.5)) and degree_rcut(dphi, 60, 80))
    cut[30] = not (rcut(pt1, (0.45, 0.5)) and rcut(pt2, (0.75, 0.8)) and degree_rcut(dphi, 60, 80))
    cut[31] = not (rcut(pt1, (0.75, 0.8)) and rcut(pt2, (0.55, 0.6)) and degree_rcut(dphi, 60, 80))
    cut[32] = not (rcut(pt1, (0.55, 0.6)) and rcut(pt2, (0.75, 0.8)) and degree_rcut(dphi, 60, 80))
    cut[33] = not (rcut(pt1, (0.75, 0.8)) and rcut(pt2, (0.6, 0.65)) and degree_rcut(dphi, 70, 80))
    cut[34] = not (rcut(pt1, (0.6, 0.65)) and rcut(pt2, (0.75, 0.8)) and degree_rcut(dphi, 70, 80))
    cut[35] = not (rcut(pt1, (0.75, 0.8)) and rcut(pt2, (0.65, 0.7)) and degree_rcut(dphi, 50, 80))
    cut[36] = not (rcut(pt1, (0.65, 0.7)) and rcut(pt2, (0.75, 0.8)) and degree_rcut(dphi, 50, 80))
    cut[37] = not (rcut(pt1, (0.75, 0.8)) and rcut(pt2, (0.7, 0.75)) and degree_rcut(dphi, 50, 70))
    cut[38] = not (rcut(pt1, (0.7, 0.75)) and rcut(pt2, (0.75, 0.8)) and degree_rcut(dphi, 50, 70))
    cut[39] = not (rcut(pt1, (0.75, 0.8)) and rcut(pt2, (0.75, 0.8)) and degree_rcut(dphi, 30, 90))
    cut[40] = not (
        rcut(pt1, (0.75, 0.8)) and rcut(pt2, (0.75, 0.8)) and degree_rcut(dphi, 140, 170)
    )

    return np.sum(cut) == len(cut)


# Apply the missing Delta phi bin veto to one event
def cut_func(event):
    _, final_proton = obs.proj_init_final_protons(event)

    pt1 = final_proton[0].pt
    pt2 = final_proton[1].pt
    dphi = final_proton[0].abs_delta_phi(final_proton[1])

    return bool(dPhi_cuts(pt1, pt2, dphi))

# Fiducial Z to dimuon selection for Pythia showered events
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

from core.analysis import obs
from core.io.cache import cache

CHARGED_FINAL_STATE_ABS_PDGS = {11, 13, 15, 211, 321, 2212, 3112, 3222, 3312, 3334}

cut_param = {
    "Z_MASS": 91.1876,
    # ----------------------------------------------
    # The event must pass these cuts exclusively
    # Muons
    "MUON_PT_MIN": 25.0,
    "MUON_ETA_MAX": 2.5,
    # Di-muon system
    "M_MUMU_MIN": 80.0,
    "M_MUMU_MAX": 100.0,
    "DETA_MUMU_MAX": 1000.0,
    "PT_MUMU_MIN": 0.0,
    "PT_MUMU_MAX": 1000.0,
    "ACOPLANARITY_MAX": 1.0,
    # Total event multiplicity
    "N_CH_MAX": 1000,
    # ----------------------------------------------
    # Once event is selected, these parameters define
    # which (charged) final states are used for N-body
    # (energy-gap-flow -- EGF) observables
    "EGF_CH_ETA_MAX": 2.5,
    "EGF_CH_ETA_MIN": -2.5,
    "EGF_CH_PT_MIN": 0.5,  # GeV
    # Algorithm parameters
    "EGF_BETA": 0.1,
    "EGF_Q0": 2.5,  # GeV
    "G2_EF_Q0": 0.15,  # GeV
    "MGF_P0": 2.5,  # GeV
    # Forward leading-proton fiducial definition
    "FORWARD_PROTON_ETA_MIN": 6.0,
    "FORWARD_PROTON_PT_MAX": 2.0,  # GeV
    "FORWARD_PROTON_XI_FLOOR": 1.0e-12,
    # ----------------------------------------------
}


# Compute true for stable charged-particle PDG IDs used in central multiplicity
def is_charged_final_state_pid(pid):
    return abs(pid) in CHARGED_FINAL_STATE_ABS_PDGS


# Compute total stable charged final-state multiplicity used by the event cut
@cache
def total_charged_multiplicity(event):
    records = obs.proj_event_particles(event)

    return sum(1 for r in records if r["is_final"] and is_charged_final_state_pid(r["pid"]))


# Compute stable fiducial muon records from one HepMC3 event
def final_state_muons(event):
    records = obs.proj_event_particles(event)
    pt_min = event.cut_param["MUON_PT_MIN"]
    eta_max = event.cut_param["MUON_ETA_MAX"]

    return [
        r
        for r in records
        if abs(r["pid"]) == 13 and r["is_final"] and r["pt"] > pt_min and abs(r["eta"]) < eta_max
    ]


# Compute dimuon acoplanarity from two four-vectors
def acoplanarity_from_p4(first, second):
    acop = 1.0 - first.abs_delta_phi(second) / math.pi
    return max(0.0, min(1.0, acop))


# Compute absolute dimuon pseudorapidity separation from two four-vectors
def deta_mumu_from_p4(first, second):
    return abs(first.eta - second.eta)


# Compute dimuon acoplanarity 1 - |Delta phi| / pi
def dimuon_acoplanarity(first, second):
    return acoplanarity_from_p4(first["p4"], second["p4"])


# Compute absolute dimuon pseudorapidity separation from two muon records
def dimuon_delta_eta(first, second):
    return deta_mumu_from_p4(first["p4"], second["p4"])


# Compute the opposite-sign dimuon pair closest to the Z pole
@cache
def select_mumu_pair(event):
    muons = final_state_muons(event)
    mass_min = event.cut_param["M_MUMU_MIN"]
    mass_max = event.cut_param["M_MUMU_MAX"]
    deta_max = event.cut_param["DETA_MUMU_MAX"]
    pt_min = event.cut_param["PT_MUMU_MIN"]
    pt_max = event.cut_param["PT_MUMU_MAX"]
    acop_max = event.cut_param["ACOPLANARITY_MAX"]
    z_mass = event.cut_param["Z_MASS"]

    best = None
    best_score = None

    for i, first in enumerate(muons):
        for second in muons[i + 1 :]:
            if first["pid"] * second["pid"] > 0:
                continue

            pair = first["p4"] + second["p4"]
            if (pair.m < mass_min) or (pair.m > mass_max):
                continue
            if dimuon_delta_eta(first, second) >= deta_max:
                continue
            if (pair.pt < pt_min) or (pair.pt > pt_max):
                continue
            if dimuon_acoplanarity(first, second) >= acop_max:
                continue

            score = abs(pair.m - z_mass)
            if best_score is None or score < best_score:
                ordered = sorted([first["p4"], second["p4"]], key=lambda p: p.pt, reverse=True)
                best = {
                    "pair": pair,
                    "leading": ordered[0],
                    "subleading": ordered[1],
                    "records": (first, second),
                }
                best_score = score

    return best


# Compute true when an event contains one selected fiducial dimuon pair
def cut_func(event):
    return (
        select_mumu_pair(event) is not None
        and total_charged_multiplicity(event) <= event.cut_param["N_CH_MAX"]
    )

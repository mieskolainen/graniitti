# ATLAS PbPb light-by-light fiducial selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

from core.analysis import obs
from core.io.cache import cache

cut_param = {
    "photon_abs_eta_max": 2.4,
    "photon_pt_min": 2.5,
    "mass_min": 5.0,
    "system_pt_max": 1.0,
    "acoplanarity_max": 0.01,
}


# Compute the two final state central photons
@cache
def selected_photons(event):
    records = [
        record for record in obs.proj_central_particle_records(event) if record["pid"] == 22
    ]
    return records if len(records) == 2 else []


# Compute the selected diphoton four momentum
@cache
def diphoton(event):
    photons = selected_photons(event)
    if len(photons) != 2:
        return None
    return photons[0]["p4"] + photons[1]["p4"]


# Compute the diphoton acoplanarity
@cache
def acoplanarity(event):
    photons = selected_photons(event)
    if len(photons) != 2:
        return math.inf
    return 1.0 - photons[0]["p4"].abs_delta_phi(photons[1]["p4"]) / math.pi


# Apply the published ATLAS UPC light-by-light fiducial selection
def fiducial_cut(event):
    photons = selected_photons(event)
    if len(photons) != 2:
        return False
    for record in photons:
        abs_eta = abs(record["eta"])
        if abs_eta >= event.cut_param["photon_abs_eta_max"]:
            return False
        if record["pt"] <= event.cut_param["photon_pt_min"]:
            return False
    system = diphoton(event)
    return (
        system.m > event.cut_param["mass_min"]
        and system.pt < event.cut_param["system_pt_max"]
        and acoplanarity(event) < event.cut_param["acoplanarity_max"]
    )


cut_func = fiducial_cut

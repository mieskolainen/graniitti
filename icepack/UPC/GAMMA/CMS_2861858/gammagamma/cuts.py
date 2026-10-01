# CMS PbPb light by light fiducial selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

from core.analysis import obs
from core.io.cache import cache

from .. import cuts as common

cut_param = {
    **common.cut_param,
    "photon_abs_eta_max": 2.2,
    "photon_pt_min": 2.0,
    "mass_min": 5.0,
    "system_pt_max": 1.0,
    "acoplanarity_max": 0.01,
}


# Compute the two final state central photon candidates
@cache
def selected_photons(event):
    records = [
        record for record in obs.proj_central_particle_records(event)
        if record["pid"] == 22
        and record["pt"] > event.cut_param["photon_pt_min"]
        and abs(record["eta"]) < event.cut_param["photon_abs_eta_max"]
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


# Apply the published CMS UPC light by light fiducial selection
def fiducial_cut(event):
    photons = selected_photons(event)
    if len(photons) != 2:
        return False
    system = diphoton(event)
    return (
        system.m > event.cut_param["mass_min"]
        and system.pt < event.cut_param["system_pt_max"]
        and acoplanarity(event) < event.cut_param["acoplanarity_max"]
        and common.exclusivity_cut(event, photons)
    )


cut_func = fiducial_cut

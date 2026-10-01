# CMS PbPb Breit-Wheeler fiducial selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

from core.analysis import obs
from core.io.cache import cache

from .. import cuts as common

cut_param = {
    **common.cut_param,
    "electron_abs_eta_max": 2.2,
    "electron_et_min": 2.0,
    "mass_min": 5.0,
    "system_pt_max": 1.0,
    "acoplanarity_max": 0.01,
}


# Compute the opposite sign central electron pair ordered by transverse momentum
@cache
def selected_electrons(event):
    records = [
        record
        for record in obs.proj_central_particle_records(event)
        if abs(record["pid"]) == 11
    ]
    if len(records) != 2 or {record["pid"] for record in records} != {-11, 11}:
        return []
    return sorted(records, key=lambda record: record["pt"], reverse=True)


# Compute the selected dielectron four momentum
@cache
def dielectron(event):
    electrons = selected_electrons(event)
    if len(electrons) != 2:
        return None
    return electrons[0]["p4"] + electrons[1]["p4"]


# Compute the dielectron acoplanarity
@cache
def acoplanarity(event):
    electrons = selected_electrons(event)
    if len(electrons) != 2:
        return math.inf
    return 1.0 - electrons[0]["p4"].abs_delta_phi(electrons[1]["p4"]) / math.pi


# Apply the published CMS UPC dielectron fiducial selection
def fiducial_cut(event):
    electrons = selected_electrons(event)
    if len(electrons) != 2:
        return False
    if not all(
        abs(record["eta"]) < event.cut_param["electron_abs_eta_max"]
        and record["p4"].mt > event.cut_param["electron_et_min"]
        for record in electrons
    ):
        return False
    system = dielectron(event)
    return (
        system.m > event.cut_param["mass_min"]
        and system.pt < event.cut_param["system_pt_max"]
        and acoplanarity(event) < event.cut_param["acoplanarity_max"]
        and common.exclusivity_cut(event, electrons)
    )


cut_func = fiducial_cut

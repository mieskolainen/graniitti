# ATLAS PbPb exclusive dimuon fiducial selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache

cut_param = {
    "muon_abs_eta_max": 2.4,
    "muon_pt_min": 4.0,
    "mass_min": 10.0,
    "mass_max": 200.0,
    "system_pt_max": 2.0,
}


# Compute the two final state central muons
@cache
def selected_muons(event):
    records = [
        record for record in obs.proj_central_particle_records(event) if abs(record["pid"]) == 13
    ]
    if len(records) != 2 or {record["pid"] for record in records} != {-13, 13}:
        return []
    return records


# Compute the selected dimuon four momentum
@cache
def dimuon(event):
    muons = selected_muons(event)
    if len(muons) != 2:
        return None
    return muons[0]["p4"] + muons[1]["p4"]


# Apply the ATLAS UPC dimuon fiducial selection [REFERENCE: arXiv:2011.12211, Section 1]
def fiducial_cut(event):
    muons = selected_muons(event)
    if len(muons) != 2:
        return False
    if not all(
        abs(record["eta"]) < event.cut_param["muon_abs_eta_max"]
        and record["pt"] > event.cut_param["muon_pt_min"]
        for record in muons
    ):
        return False
    system = dimuon(event)
    return (
        event.cut_param["mass_min"] < system.m < event.cut_param["mass_max"]
        and system.pt < event.cut_param["system_pt_max"]
    )


cut_func = fiducial_cut

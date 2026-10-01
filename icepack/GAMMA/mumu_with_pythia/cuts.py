# Fiducial forward excitation dimuon selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache

cut_param = {
    "MUON_ABS_ETA_MAX": 2.4,
    "MUON_PT_MIN": 10.0,
    "M_MUMU_MIN": 20.0,
    "M_MUMU_MAX": 200.0,
    "M_Z_VETO_MIN": 70.0,
    "M_Z_VETO_MAX": 105.0,
    "PT_MUMU_MAX": 1.5,
}


# Compute the highest transverse momentum opposite sign fiducial muon pair
@cache
def selected_muons(event):
    muons = [
        record
        for record in obs.proj_event_particles(event)
        if record["is_final"]
        and abs(record["pid"]) == 13
        and record["pt"] > event.cut_param["MUON_PT_MIN"]
        and abs(record["eta"]) < event.cut_param["MUON_ABS_ETA_MAX"]
    ]
    pairs = [
        (first, second)
        for index, first in enumerate(muons)
        for second in muons[index + 1 :]
        if first["pid"] * second["pid"] < 0
    ]
    if len(pairs) == 0:
        return None
    pair = max(pairs, key=lambda records: records[0]["pt"] + records[1]["pt"])
    return sorted(pair, key=lambda record: record["pt"], reverse=True)


# Apply the ATLAS inspired exclusive dimuon fiducial selection
def cut_func(event):
    muons = selected_muons(event)
    if muons is None:
        return False

    pair = muons[0]["p4"] + muons[1]["p4"]
    mass = pair.m
    return (
        event.cut_param["M_MUMU_MIN"] < mass < event.cut_param["M_MUMU_MAX"]
        and not event.cut_param["M_Z_VETO_MIN"] < mass < event.cut_param["M_Z_VETO_MAX"]
        and pair.pt < event.cut_param["PT_MUMU_MAX"]
    )

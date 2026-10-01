# ATLAS 7 TeV exclusive dimuon event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

cut_param = {
    "MUON_ABS_ETA_MAX": 2.4,
    "MUON_PT_MIN": 10.0,
    "M_MUMU_MIN": 20.0,
    "M_SYSTEM_VETO": [70.0, 105.0],
    "SYSTEM_PT_MAX": 1.5,
}


# Apply the ATLAS fiducial muon cuts, mass veto and pair pT selection
def cut_func(event):
    muons = [
        record
        for record in obs.proj_event_particles(event)
        if record["is_final"] and abs(record["pid"]) == 13
    ]
    if len(muons) != 2 or muons[0]["pid"] != -muons[1]["pid"]:
        return False

    eta_max = event.cut_param["MUON_ABS_ETA_MAX"]
    pt_min = event.cut_param["MUON_PT_MIN"]
    if not all(abs(record["eta"]) < eta_max and record["pt"] > pt_min for record in muons):
        return False

    system = muons[0]["p4"] + muons[1]["p4"]
    system_mass = system.m
    veto_min, veto_max = event.cut_param["M_SYSTEM_VETO"]
    return (
        system_mass > event.cut_param["M_MUMU_MIN"]
        and not (veto_min < system_mass < veto_max)
        and system.pt < event.cut_param["SYSTEM_PT_MAX"]
    )

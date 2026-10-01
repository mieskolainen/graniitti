# ATLAS 7 TeV exclusive dielectron event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

cut_param = {
    "ELECTRON_ABS_ETA_MAX": 2.4,
    "ELECTRON_PT_MIN": 12.0,
    "M_EE_MIN": 24.0,
    "M_SYSTEM_VETO": [70.0, 105.0],
    "SYSTEM_PT_MAX": 1.5,
}


# Apply the ATLAS fiducial electron cuts, mass veto and pair pT selection
def cut_func(event):
    electrons = [
        record
        for record in obs.proj_event_particles(event)
        if record["is_final"] and abs(record["pid"]) == 11
    ]
    if len(electrons) != 2 or electrons[0]["pid"] != -electrons[1]["pid"]:
        return False

    eta_max = event.cut_param["ELECTRON_ABS_ETA_MAX"]
    pt_min = event.cut_param["ELECTRON_PT_MIN"]
    if not all(abs(record["eta"]) < eta_max and record["pt"] > pt_min for record in electrons):
        return False

    system = electrons[0]["p4"] + electrons[1]["p4"]
    system_mass = system.m
    veto_min, veto_max = event.cut_param["M_SYSTEM_VETO"]
    return bool(
        system_mass > event.cut_param["M_EE_MIN"]
        and not (veto_min < system_mass < veto_max)
        and system.pt < event.cut_param["SYSTEM_PT_MAX"]
    )

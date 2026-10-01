# CDF chi_c to J/psi gamma fiducial event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

cut_param = {"MUON_ABS_ETA_MAX": 0.6, "MUON_PT_MIN": 1.4}


# Apply the published CDF cuts to both final state muons
def cut_func(event):
    muons = [
        record
        for record in obs.proj_event_particles(event)
        if record["is_final"] and abs(record["pid"]) == 13
    ]
    if len(muons) != 2:
        return False
    eta_max = event.cut_param["MUON_ABS_ETA_MAX"]
    pt_min = event.cut_param["MUON_PT_MIN"]
    return all(abs(record["eta"]) < eta_max and record["pt"] > pt_min for record in muons)

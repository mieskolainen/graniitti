# LHCb shared exclusive Upsilon event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

cut_param = {
    "CENTRAL_RAPIDITY": [2.0, 4.5],
    "MUON_ETA": [2.0, 4.5],
}


# Compute the stable muons in one HepMC3 event
def final_state_muons(event):
    return [
        record
        for record in obs.proj_event_particles(event)
        if abs(record["pid"]) == 13 and record["is_final"]
    ]


# Apply the published open central-system rapidity interval
def central_system_cut(event):
    rapidity_min, rapidity_max = event.cut_param["CENTRAL_RAPIDITY"]
    rapidity = obs.proj_1D_Rap(event)
    return rapidity_min < rapidity < rapidity_max


# Apply the published rapidity and opposite-sign dimuon fiducial selection
def fiducial_cut(event):
    if not central_system_cut(event):
        return False
    muons = final_state_muons(event)
    if len(muons) != 2:
        return False
    if sum(record["pid"] == 13 for record in muons) != 1:
        return False
    if sum(record["pid"] == -13 for record in muons) != 1:
        return False

    eta_min, eta_max = event.cut_param["MUON_ETA"]
    return all(eta_min < record["eta"] < eta_max for record in muons)


cut_func = fiducial_cut

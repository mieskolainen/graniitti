# H1 elastic pion pair photoproduction fiducial selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs, pdg

cut_param = {
    "W": [20.0, 80.0],
    "ABS_T_MAX": 1.5,
    "M": [0.5, 2.2],
    "Q2_MAX": 2.5,
}


# Compute one unique particle momentum selected by PDG and event status
def unique_momentum(event, abs_pid, initial):
    records = obs.proj_event_particles(event)
    selected = [
        record["p4"]
        for record in records
        if abs(record["pid"]) == abs_pid
        and ((record["status"] == pdg.INITIAL_STATE) if initial else record["is_final"])
    ]
    if len(selected) != 1:
        raise ValueError(
            f"H1 elastic selection requires one {'initial' if initial else 'final'} PDG {abs_pid}"
        )
    return selected[0]


# Compute the unique stable charged pion pair
def pion_pair(event):
    records = obs.proj_event_particles(event)
    plus = [record["p4"] for record in records if record["pid"] == 211 and record["is_final"]]
    minus = [record["p4"] for record in records if record["pid"] == -211 and record["is_final"]]
    if len(plus) != 1 or len(minus) != 1:
        raise ValueError("H1 elastic selection requires one stable pi+ pi- pair")
    return plus[0] + minus[0]


# Apply the H1 measured electron-proton phase space
def cut_func(event):
    electron_in = unique_momentum(event, 11, True)
    electron_out = unique_momentum(event, 11, False)
    proton_in = unique_momentum(event, 2212, True)
    proton_out = unique_momentum(event, 2212, False)
    pair = pion_pair(event)
    photon = electron_in - electron_out
    q2 = -photon.m2
    w2 = (photon + proton_in).m2
    transfer = (proton_in - proton_out).m2
    if w2 <= 0.0:
        return False
    w_min, w_max = event.cut_param["W"]
    m_min, m_max = event.cut_param["M"]
    return (
        w_min < w2**0.5 < w_max
        and abs(transfer) < event.cut_param["ABS_T_MAX"]
        and m_min < pair.m < m_max
        and 0.0 <= q2 < event.cut_param["Q2_MAX"]
    )

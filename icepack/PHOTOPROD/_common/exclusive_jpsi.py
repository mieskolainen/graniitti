# Shared exclusive HERA vector meson event reconstruction
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs, pdg
from core.io.cache import cache


# Compute one unique initial or stable final momentum selected by absolute PDG
def unique_momentum(event, abs_pid, initial, records=None):
    records = obs.proj_event_particles(event) if records is None else records
    central = set() if initial else {record["id"] for record in obs.proj_central_particle_records(event)}
    selected = [
        record["p4"]
        for record in records
        if abs(record["pid"]) == abs_pid and record["id"] not in central
        and ((record["status"] == pdg.INITIAL_STATE) if initial else record["is_final"])
    ]
    if len(selected) != 1:
        state = "initial" if initial else "final"
        raise ValueError(
            f"Exclusive HERA selection requires one {state} PDG {abs_pid}, found {len(selected)}"
        )
    return selected[0]


# Reconstruct exact photon virtuality, gamma-proton mass and proton transfer
@cache
def kinematics(event):
    records = obs.proj_event_particles(event)
    lepton_in = unique_momentum(event, 11, True, records)
    lepton_out = unique_momentum(event, 11, False, records)
    proton_in = unique_momentum(event, 2212, True, records)
    proton_out = unique_momentum(event, 2212, False, records)
    photon = lepton_in - lepton_out
    w2 = (photon + proton_in).m2
    return {
        "q2": -photon.m2,
        "w": w2**0.5 if w2 > 0.0 else 0.0,
        "abs_t": abs((proton_in - proton_out).m2),
    }


# Apply one published exclusive HERA photoproduction phase-space region
def accepted(event, cuts=None):
    cuts = event.cut_param if cuts is None else cuts
    selected = sorted(record["pid"] for record in obs.proj_central_particle_records(event))
    if selected != sorted(obs.flatten_pid_selection(event.pid)):
        return False
    values = kinematics(event)
    w_min, w_max = cuts["W"]
    return (
        w_min < values["w"] < w_max
        and 0.0 <= values["q2"] < cuts["Q2_MAX"]
        and cuts.get("ABS_T_MIN", 0.0) <= values["abs_t"]
        < cuts.get("ABS_T_MAX", float("inf"))
        and ("M" not in cuts
             or cuts["M"][0] < obs.proj_central_system(event).m < cuts["M"][1])
    )


# Compute the absolute proton-vertex momentum transfer
def abs_t(event):
    return kinematics(event)["abs_t"]

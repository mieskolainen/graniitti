# STAR 200 GeV fiducial central pair observable projectors
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs, pdg
from core.io.cache import cache

# See generator card cuts and MUserCuts.cc also
CENTRAL_PT_REGION = {
    211: (0.2, None),
    321: (0.3, 0.7),
    2212: (0.4, 1.1),
}
CENTRAL_ABS_ETA_MAX = 0.7


# Compute the STAR transverse momentum region for one central pair species
def _central_pt_region(event):
    species = {abs(int(pid_value)) for pid_value in event.pid}
    if len(species) != 1:
        raise ValueError(f"STAR central pair requires one particle species, found {event.pid}")
    abs_pid = species.pop()
    if abs_pid not in CENTRAL_PT_REGION:
        raise ValueError(f"STAR central pair has unsupported absolute PDG code {abs_pid}")
    return CENTRAL_PT_REGION[abs_pid]


# Select the requested central pair inside its published STAR fiducial region
@cache
def central_particles(event):
    pt_min, pair_min_pt_max = _central_pt_region(event)
    records = obs.proj_central_particle_records(event)
    pair = []
    for pid_value in event.pid:
        candidates = [
            record
            for record in records
            if record["pid"] == pid_value
            and record["pt"] > pt_min
            and abs(record["eta"]) < CENTRAL_ABS_ETA_MAX
        ]
        if len(candidates) != 1:
            raise ValueError(
                "STAR central pair requires exactly one fiducial particle "
                f"with PDG {pid_value}, found {len(candidates)}"
            )
        pair.append(candidates[0])

    if pair_min_pt_max is not None and min(record["pt"] for record in pair) >= pair_min_pt_max:
        raise ValueError("STAR central pair is outside its upper transverse-momentum region")
    return [record["p4"] for record in pair]


# Compute the STAR fiducial central pair invariant mass
def invariant_mass(event):
    return sum(central_particles(event)).m


# Compute the STAR fiducial central pair rapidity
def rapidity(event):
    return sum(central_particles(event)).rapidity


# Transform the STAR fiducial central pair to the Collins-Soper frame
@cache
def collins_soper_pair(event):
    particles = central_particles(event)
    central_system = sum(particles)
    beams = [
        record["p4"]
        for record in obs.proj_event_particles(event)
        if record["status"] == pdg.INITIAL_STATE
    ]
    if len(beams) != 2:
        raise ValueError(f"STAR central pair requires two beam particles, found {len(beams)}")
    beam_1, beam_2, pair = obs.LorentFramePrepare(
        pbeam1=beams[0], pbeam2=beams[1], particles=particles, X=central_system
    )
    return obs.LorentzFrame(pb1boost=beam_1, pb2boost=beam_2, pfboost=pair, frametype="CS")


# Compute the positive central particle Collins-Soper polar angle
def collins_soper_costheta(event):
    return collins_soper_pair(event)[0].costheta


# Compute the positive central particle Collins-Soper azimuthal angle in degrees
def collins_soper_phi(event):
    return obs.rad2deg(collins_soper_pair(event)[0].phi)

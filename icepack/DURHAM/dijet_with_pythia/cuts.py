# Particle level Durham dijet selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
from core.analysis import obs
from core.io.cache import cache
from core.kinematics import fastjet

NEUTRINO_ABS_PDGS = {12, 14, 16}

cut_param = {
    "JET_RADIUS": 0.6,
    "JET_PT_MIN": 20.0,
    "JET_ABS_ETA_MAX": 2.5,
    "PARTICLE_PT_MIN": 0.1,
    "PARTICLE_ABS_ETA_MAX": 5.0,
    "DIJET_M_MIN": 40.0,
    "DIJET_M_MAX": 200.0,
}


# Compute stable visible particle arrays used as anti-kT inputs
def stable_visible_particles(event):
    particle_pt_min  = event.cut_param["PARTICLE_PT_MIN"]
    particle_eta_max = event.cut_param["PARTICLE_ABS_ETA_MAX"]
    records = [
        record
        for record in obs.proj_event_particles(event)
        if record["is_final"]
        and abs(record["pid"]) not in NEUTRINO_ABS_PDGS
        and record["pt"] > particle_pt_min
        and abs(record["eta"]) < particle_eta_max
    ]
    return {
        "px": np.asarray([record["p4"].x for record in records], dtype=np.float64),
        "py": np.asarray([record["p4"].y for record in records], dtype=np.float64),
        "pz": np.asarray([record["p4"].z for record in records], dtype=np.float64),
        "energy": np.asarray([record["p4"].t for record in records], dtype=np.float64),
    }


# Cluster particle level anti-kT jets and cache their constituent labels
@cache
def particle_level_jets(event):
    particles = stable_visible_particles(event)
    jets, labels = fastjet.cluster_antikt(
        particles["px"],
        particles["py"],
        particles["pz"],
        particles["energy"],
        event.cut_param["JET_RADIUS"],
        event.cut_param["JET_PT_MIN"],
        event.cut_param["JET_ABS_ETA_MAX"],
    )
    particles["jets"] = jets
    particles["labels"] = labels
    return particles


# Compute the leading dijet invariant mass or zero when fewer than two jets exist
def dijet_mass(jets):
    if jets.shape[0] < 2:
        return 0.0
    momentum = jets[0] + jets[1]
    return fastjet.jet_mass(momentum[0], momentum[1], momentum[2], momentum[3])


# Require two fiducial anti-kT jets inside the dijet mass window
def cut_func(event):
    jets = particle_level_jets(event)["jets"]
    mass = dijet_mass(jets)
    return (
        jets.shape[0] >= 2
        and event.cut_param["DIJET_M_MIN"] < mass < event.cut_param["DIJET_M_MAX"]
    )

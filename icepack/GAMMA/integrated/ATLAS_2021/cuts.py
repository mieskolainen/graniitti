# ATLAS 13 TeV photon-induced WW fiducial event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache

cut_param = {
    "LEPTON_ABS_ETA_MAX": 2.5,
    "LEADING_LEPTON_PT_MIN": 27.0,
    "SUBLEADING_LEPTON_PT_MIN": 20.0,
    "M_EMU_MIN": 20.0,
    "PT_EMU_MIN": 30.0,
    "TRACK_PT_MIN": 0.5,
    "TRACK_ABS_ETA_MAX": 2.5,
    "DRESSING_DR_MAX": 0.1,
}

CHARGED_ABS_PDGS = {11, 13, 15, 211, 321, 2212, 3112, 3222, 3312, 3334}


# Compute every ancestor particle in the HepMC event graph
def ancestors(particle):
    ancestors = []
    pending = []
    vertex = particle.production_vertex()
    if vertex is not None:
        pending.extend(vertex.particles_in())
    visited = set()
    while pending:
        parent = pending.pop()
        if parent.id() in visited:
            continue
        visited.add(parent.id())
        ancestors.append(parent)
        vertex = parent.production_vertex()
        if vertex is not None:
            pending.extend(vertex.particles_in())
    return ancestors


# Compute true when one ancestry contains a decaying hadron
def has_hadron_ancestor(particle):
    return any(parent.status() == 2 and 100 <= abs(parent.pid()) < 1_000_000_000 for parent in ancestors(particle))


# Compute final prompt electron and muon particle objects
@cache
def prompt_lepton_particles(event):
    candidates = [
        particle
        for particle in event.evt.particles()
        if particle.status() == 1
        and particle.end_vertex() is None
        and abs(particle.pid()) in {11, 13}
        and not has_hadron_ancestor(particle)
        and 15 not in {abs(parent.pid()) for parent in ancestors(particle)}
        and 24 in {abs(parent.pid()) for parent in ancestors(particle)}
    ]
    if len(candidates) != 2:
        return []
    selected = []
    for flavour in (11, 13):
        options = [particle for particle in candidates if abs(particle.pid()) == flavour]
        if not options:
            return []
        selected.append(options[0])
    if selected[0].pid() * selected[1].pid() >= 0:
        return []
    return selected


# Compute one lepton four momentum dressed with its QED final state photons
def dressed_momentum(event, lepton):
    records = {record["id"]: record for record in obs.proj_event_particles(event)}
    momentum = records[lepton.id()]["p4"]
    direction = momentum
    dr_max = event.cut_param["DRESSING_DR_MAX"]
    for particle in event.evt.particles():
        if particle.pid() != 22 or particle.status() != 1 or particle.end_vertex() is not None:
            continue
        if has_hadron_ancestor(particle):
            continue
        photon = records[particle.id()]["p4"]
        if direction.deltaR(photon) < dr_max:
            momentum = momentum + photon
    return momentum


# Compute the dressed prompt electron and muon generated in the selected WW channel
@cache
def selected_leptons(event):
    particles = prompt_lepton_particles(event)
    if len(particles) != 2:
        return []
    leptons = []
    for particle in particles:
        momentum = dressed_momentum(event, particle)
        leptons.append(
            {
                "id": particle.id(),
                "pid": particle.pid(),
                "eta": momentum.eta,
                "pt": momentum.pt,
                "p4": momentum,
            }
        )
    return sorted(leptons, key=lambda record: record["pt"], reverse=True)


# Compute whether an additional charged particle passes the ATLAS track selection
def has_extra_track(event):
    pt_min = event.cut_param["TRACK_PT_MIN"]
    eta_max = event.cut_param["TRACK_ABS_ETA_MAX"]
    selected_ids = {record["id"] for record in selected_leptons(event)}
    return any(
        record["is_final"]
        and abs(record["pid"]) in CHARGED_ABS_PDGS
        and record["id"] not in selected_ids
        and record["pt"] > pt_min
        and abs(record["eta"]) < eta_max
        for record in obs.proj_event_particles(event)
    )


# Apply the ATLAS particle-level e-mu fiducial selection
def cut_func(event):
    leptons = selected_leptons(event)
    if not leptons:
        return False
    leading, subleading = leptons
    eta_max = event.cut_param["LEPTON_ABS_ETA_MAX"]
    if not all(abs(record["eta"]) < eta_max for record in leptons):
        return False
    if leading["pt"] < event.cut_param["LEADING_LEPTON_PT_MIN"]:
        return False
    if subleading["pt"] <= event.cut_param["SUBLEADING_LEPTON_PT_MIN"]:
        return False
    if has_extra_track(event):
        return False
    pair = leading["p4"] + subleading["p4"]
    return bool(
        pair.m > event.cut_param["M_EMU_MIN"]
        and pair.pt > event.cut_param["PT_EMU_MIN"]
    )

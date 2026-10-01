# Construct momentum-conserving decay events with the official HepMC3 binding
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
from types import SimpleNamespace

from core.kinematics.vec4 import vec4
from pyHepMC3 import HepMC3 as h
from pyHepMC3 import std


# Construct one decay vertex from (PDG id, four-momentum, status) tuples
def decay_event(particles, *, pid=None, cut_param=None):
    event = h.GenEvent()
    vertex = h.GenVertex()
    total = vec4()
    for pdg, momentum, status in particles:
        total += momentum
        vertex.add_particle_out(h.GenParticle(
            h.FourVector(momentum.px, momentum.py, momentum.pz, momentum.e), pdg, status
        ))
    vertex.add_particle_in(h.GenParticle(
        h.FourVector(total.px, total.py, total.pz, total.e), 990, 2
    ))
    event.add_vertex(vertex)
    return SimpleNamespace(evt=event, pid=pid or [p[0] for p in particles], cut_param=cut_param or {})


# Create one HepMC particle for exclusive pair events
def make_particle(px, py, pz, energy, pid, status):
    return h.GenParticle(h.FourVector(px, py, pz, energy), pid, status)


# Construct exclusive pp -> pp + pair kinematics with on-shell particles and conserved vertices
def make_pair_event(mass=2.264, pid=2212, dphi=180.0, phi=0.0, forward_pt=0.1, costheta=0.2, beam_energy=100.0):
    proton_mass = 0.9382720813
    daughter_mass = proton_mass if pid == 2212 else 0.13957039
    beam_pz = math.sqrt(beam_energy**2 - proton_mass**2)
    angle = math.radians(dphi)
    px = -forward_pt * (1.0 + math.cos(angle))
    py = -forward_pt * math.sin(angle)
    system = vec4(px, py, 0.0, math.sqrt(mass**2 + px**2 + py**2))
    forward_energy = beam_energy - system.e / 2.0
    forward_pz = math.sqrt(forward_energy**2 - proton_mass**2 - forward_pt**2)
    beams = [vec4(0.0, 0.0, sign * beam_pz, beam_energy) for sign in (1, -1)]
    forward = [vec4(forward_pt, 0.0, forward_pz, forward_energy),
               vec4(forward_pt * math.cos(angle), forward_pt * math.sin(angle), -forward_pz, forward_energy)]
    event = h.GenEvent(h.Units.GEV, h.Units.MM)
    exchanges = []
    for beam, outgoing in zip(beams, forward, strict=True):
        transfer = beam - outgoing
        exchange = make_particle(*transfer, 990, 2)
        vertex = h.GenVertex()
        vertex.add_particle_in(make_particle(*beam, 2212, 4))
        vertex.add_particle_out(make_particle(*outgoing, 2212, 1))
        vertex.add_particle_out(exchange)
        event.add_vertex(vertex)
        exchanges.append(exchange)
    momentum = math.sqrt(mass**2 / 4.0 - daughter_mass**2)
    pt = momentum * math.sqrt(1.0 - costheta**2)
    azimuth = math.radians(phi)
    central = h.GenVertex()
    for exchange in exchanges:
        central.add_particle_in(exchange)
    for sign in (1, -1):
        daughter = vec4(sign * pt * math.cos(azimuth), sign * pt * math.sin(azimuth),
                        sign * costheta * momentum, mass / 2.0)
        daughter.boost(b=system, sign=1)
        central.add_particle_out(make_particle(*daughter, sign * pid, 1))
    event.add_vertex(central)
    return SimpleNamespace(evt=event, pid=[pid, -pid], cut_param={})


# Write conserving photon-to-muon events with physical weights and cross-section counters
def write_muon_events(path, weights=(1.0, 1.0, 1.0), *, energies=None, xsections=None,
                      attempted=None, overflow_indices=(), unit=h.Units.GEV, weight_name=None):
    writer = h.WriterAscii(str(path))
    writer.set_precision(17)
    run_info = h.GenRunInfo()
    if weight_name is not None:
        run_info.set_weight_names(std.vector_std_string([weight_name]))
    try:
        for index, weight in enumerate(weights):
            energy = 1.0 if energies is None else energies[index]
            event = h.GenEvent(h.Units.GEV, h.Units.MM)
            if weight_name is not None:
                event.set_run_info(run_info)
            event.set_event_number(index)
            data = h.GenEventData()
            event.write_data(data)
            data.weights = std.vector_double([weight])
            event.read_data(data)
            vertex = h.GenVertex()
            for sign in (-1, 1):
                vertex.add_particle_in(make_particle(0, 0, sign * energy, energy, 22, 4))
                vertex.add_particle_out(make_particle(
                    0, 0, sign * math.sqrt(energy**2 - 0.105658**2), energy, sign * 13, 1))
            event.add_vertex(vertex)
            cross_section = h.GenCrossSection()
            xs, error = (3.0, 0.3) if xsections is None else xsections[index]
            cross_section.set_cross_section(xs, error, index + 1, index + 1 if attempted is None else attempted[index])
            event.set_cross_section(cross_section)
            if index in overflow_indices:
                event.add_attribute('maximum_weight_overflow', h.DoubleAttribute(weight))
            event.set_units(unit, h.Units.MM)
            writer.write_event(event)
            assert not writer.failed()
    finally:
        writer.close()
    return str(path)


# Read complete events through the official HepMC3 binding
def read_events(path):
    events = []
    reader = h.ReaderAscii(str(path))
    try:
        while not reader.failed():
            event = h.GenEvent()
            reader.read_event(event)
            if not reader.failed():
                events.append(event)
    finally:
        reader.close()
    return events

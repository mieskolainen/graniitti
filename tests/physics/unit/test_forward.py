# Physical particle and detector acceptance tests using HepMC3
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from math import cos, hypot, pi, sin

import pytest
from core.analysis.forward import Acceptance, Plane, neutron_class
from core.kinematics.vec4 import vec4
from pyHepMC3 import HepMC3 as h


# Construct a massive on-shell forward particle
def particle(pid, px, py, pz, mass, status=1):
    p = h.GenParticle(h.FourVector(px, py, pz, hypot(px, py, pz, mass)), pid, status)
    p.set_generated_mass(mass)
    return p


# Construct two resolved ion branches and a central neutron with both beam ancestors
def event(phi=0.0, reverse=False, neutrons=(True, True)):
    record = h.GenEvent(h.Units.GEV, h.Units.MM)
    exchanges = []
    for side in (1, -1):
        sign = -side if reverse else side
        emit = neutrons[0 if side > 0 else 1]
        beam = particle(1000822080, 0, 0, sign * 500000, 193.7, 4)
        daughters = [particle(pid, pt * cos(phi), pt * sin(phi), sign * pz, mass)
                     for pid, pt, pz, mass in (
                         (2112, 0.1, 2500, 0.9396),
                         (2212, 0.2, 2400, 0.9383),
                         (22, 0.01, 20, 0),
                         (1000812060 + int(not emit) * 10, -0.31 + int(not emit) * 0.1,
                          494980 + int(not emit) * 2500, 192 + int(not emit) * 0.9396),
                     ) if pid != 2112 or emit]
        total = sum((vec4(p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().e())
                     for p in daughters), vec4())
        parent = h.GenParticle(h.FourVector(*total), 1000822080, 2)
        transfer = vec4(0, 0, sign * 500000, beam.momentum().e()) - total
        exchange = h.GenParticle(h.FourVector(*transfer), 22, 3)
        production = h.GenVertex()
        production.add_particle_in(beam)
        production.add_particle_out(parent)
        production.add_particle_out(exchange)
        record.add_vertex(production)
        decay = h.GenVertex(h.FourVector(2, 3, sign * 10, 20))
        decay.add_particle_in(parent)
        for daughter in daughters:
            decay.add_particle_out(daughter)
        record.add_vertex(decay)
        exchanges.append(exchange)
    central = h.GenVertex()
    for exchange in exchanges:
        central.add_particle_in(exchange)
    for sign in (-1, 1):
        central.add_particle_out(particle(sign * 2112, sign * 0.2, 0, sign * 3, 0.9396))
        energy = sum(p.momentum().e() for p in exchanges) / 2 - hypot(0.2, 3, 0.9396)
        central.add_particle_out(particle(22, 0, 0, sign * energy, 0))
    record.add_vertex(central)
    for vertex in record.vertices():
        incoming = sum((vec4(p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().e())
                        for p in vertex.particles_in()), vec4())
        outgoing = sum((vec4(p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().e())
                        for p in vertex.particles_out()), vec4())
        assert tuple(incoming - outgoing) == pytest.approx((0, 0, 0, 0), abs=1e-9)
    return record


# Check joint neutron, proton and gamma selections and energy sums
def test_species_and_joint_multiplicity():
    record = event()
    neutrons = Acceptance(pid=(2112,), side=1, energy=(1000, 3000))
    protons = Acceptance(pid=(2212,), side=1, pt=(0.15, 0.25))
    photons = Acceptance(pid=(22,), side=-1, energy=(10, 30))
    assert len(neutrons.particles(record)) == len(protons.particles(record)) == len(photons.particles(record)) == 1
    assert neutrons.energy_sum(record) == pytest.approx(hypot(2500, 0.1, 0.9396))
    assert not Acceptance(pid=(-2212,)).particles(record)


# Separate the two forward systems using standard vertices rather than particle attributes
def test_beam_ancestry_excludes_central_neutron():
    record = event()
    assert len(Acceptance(pid=(2112,), side=1).particles(record)) == 2
    assert len(Acceptance(pid=(2112,), beam=1).particles(record)) == 1
    assert len(Acceptance(pid=(2112,), beam=2).particles(record)) == 1


# Classify the exclusive particle sectors independently of central neutron production
@pytest.mark.parametrize("neutrons", [(False, False), (False, True), (True, False), (True, True)])
def test_neutron_classes(neutrons):
    for reverse in (False, True):
        ordered = neutrons[::-1] if reverse else neutrons
        expected = "".join("Xn" if n else "0n" for n in ordered)
        for phi in (0.0, 0.73):
            assert neutron_class(event(phi, reverse, neutrons)) == expected


# Preserve physical selections under a beam exchange and rotations about the beam
@pytest.mark.parametrize("phi", [0.0, 0.4, pi - 0.01, -pi + 0.01])
def test_beam_exchange_and_rotation(phi):
    for side in (-1, 1):
        selection = Acceptance(
            pid=(2112,), side=side, beam=1 if side == 1 else 2, abs_eta=(9, 12), theta=(0, 0.0001), pt=(0.09, 0.11)
        )
        assert len(selection.particles(event(phi, reverse=True))) == 1
    crossed = Acceptance(pid=(2112,), beam=1, phi=(3.0, -3.0))
    assert bool(crossed.particles(event(phi))) == (abs(phi) > 3.0)


# Convert MeV and cm records without mutating the input or changing acceptance
def test_units_and_serialization(tmp_path):
    record = event()
    selection = Acceptance(
        pid=(2112,), beam=1, energy=(2400, 2600), plane=Plane(140000, x=(7, 9), y=(2, 4), ct=(140000, 140020))
    )
    assert len(selection.particles(record)) == 1
    record.set_units(h.Units.MEV, h.Units.CM)
    path = str(tmp_path / "forward.hepmc3")
    writer = h.WriterAscii(path)
    writer.set_precision(17)
    writer.write_event(record)
    writer.close()
    decoded = h.GenEvent()
    reader = h.ReaderAscii(path)
    reader.read_event(decoded)
    assert not reader.failed()
    reader.close()
    assert len(selection.particles(decoded)) == 1
    assert decoded.momentum_unit() == h.Units.MEV
    assert decoded.length_unit() == h.Units.CM


# Include the production displacement and finite neutron velocity in detector arrival time
@pytest.mark.parametrize("side", [-1, 1])
def test_neutral_flight(side):
    record = event()
    neutron = Acceptance(pid=(2112,), beam=1 if side > 0 else 2).particles(record)[0]
    detector = Plane(side * 140000)
    hit = detector.position(neutron, record)
    assert hit.x() == pytest.approx(2 + 139990 * 0.1 / 2500)
    assert hit.t() == pytest.approx(20 + 139990 * neutron.momentum().e() / 2500)
    assert Plane(-side * 140000).position(neutron, record) is None
    assert not Plane(side * 140000, ct=(0, 139999)).accepts(neutron, record)


# Require a beam transport calculation for charged particles and use its positions
def test_proton_transport():
    record = event()
    with pytest.raises(ValueError, match="transport"):
        Acceptance(pid=(2212,), beam=1, plane=Plane(140000)).particles(record)

    # Propagate through a field-free drift with the exact proton velocity
    def optics(p, z, evt):
        start = p.production_vertex().position()
        momentum = p.momentum()
        distance = z - start.z()
        return h.FourVector(start.x() + distance * momentum.px() / momentum.pz(),
                           start.y() + distance * momentum.py() / momentum.pz(), z,
                           start.t() + distance * momentum.e() / momentum.pz())

    selected = Acceptance(pid=(2212,), beam=1, plane=Plane(140000, x=(13, 15), transport=optics))
    assert len(selected.particles(record)) == 1
    assert not Acceptance(pid=(2212,), beam=1, plane=Plane(140000, x=(70, 80), transport=optics)).particles(record)


# Reject unresolved nuclear states before a zero-particle veto can pass
@pytest.mark.parametrize(("pid", "status"), [(1000822080, 3), (91, 1)])
def test_unresolved_is_not_zero(pid, status):
    record = event()
    unresolved = particle(pid, 0, 0, 100, 194, status)
    record.add_particle(unresolved)
    with pytest.raises(ValueError, match="unresolved"):
        Acceptance(pid=(-2212,), energy=(1, 2)).particles(record)
    with pytest.raises(ValueError, match="unresolved"):
        neutron_class(record)


# Reject invalid intervals before processing events
@pytest.mark.parametrize(
    "params", [{"side": 2}, {"beam": 0}, {"pt": (2, 1)}, {"phi": (-4, 1)}, {"theta": (-1, 0.1)}, {"pid": (0,)}]
)
def test_invalid_selection(params):
    with pytest.raises(ValueError):
        Acceptance(**params)

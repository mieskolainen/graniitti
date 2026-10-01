# Test topology-only HepMC3 particle projectors
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from types import SimpleNamespace

import numpy as np
import pytest
from core.analysis import obs
from core.kinematics.vec4 import hepmc2vec4, vec4
from pyHepMC3 import HepMC3 as hepmc3

from icepack.UPC.EMD.ALICE_2149540 import obs as alice_emd_obs
from icepack.UPC.PHOTOPROD.ALICE_1840600.jpsi_coherent import obs as alice_coherent_obs
from icepack.UPC.PHOTOPROD.ALICE_2658375.jpsi_incoherent import cuts as alice_cuts
from icepack.UPC.PHOTOPROD.ALICE_2658375.jpsi_incoherent import obs as alice_obs


# Create one HepMC3 particle with Cartesian four-momentum
def make_particle(px, py, pz, energy, pid, status):
    return hepmc3.GenParticle(hepmc3.FourVector(px, py, pz, energy), pid, status)


# Construct pp or ion CEP with exact recoil and on-shell forward particles
def make_cep_event(particles, *, beam_pid=2212, beam_mass=0.938272, beam_energy=3500.0,
                   forward_masses=None, beam_reverse=False):
    event = hepmc3.GenEvent(hepmc3.Units.GEV, hepmc3.Units.MM)
    central = sum((momentum for _, momentum in particles), vec4())
    masses = np.broadcast_to(beam_mass, (2,))
    energies = np.broadcast_to(beam_energy, (2,))
    pids = np.broadcast_to(beam_pid, (2,))
    beams = [vec4(0.0, 0.0, sign * np.sqrt(energy**2 - mass**2), energy)
             for sign, energy, mass in zip((1, -1), energies, masses, strict=True)]
    recoil = beams[0] + beams[1] - central
    m1, m2 = forward_masses or masses
    mass = recoil.m
    q = np.sqrt((mass**2 - (m1 + m2)**2) * (mass**2 - (m1 - m2)**2)) / (2.0 * mass)
    e1 = (mass**2 + m1**2 - m2**2) / (2.0 * mass)
    pz = np.sqrt(q**2 - 0.2**2 - 0.1**2)
    forwards = [vec4(0.2, 0.1, pz, e1), vec4(-0.2, -0.1, -pz, mass - e1)]
    vertices, exchanges, transfers = [], [], []
    for beam, outgoing, pid in zip(beams, forwards, pids, strict=True):
        outgoing.boost(b=recoil, sign=1)
        transfer = beam - outgoing
        vertex = hepmc3.GenVertex()
        vertex.add_particle_in(make_particle(*beam, int(pid), 4))
        vertex.add_particle_out(make_particle(*outgoing, int(pid), 1))
        exchange = make_particle(*transfer, 99, 81)
        vertex.add_particle_out(exchange)
        vertices.append(vertex)
        exchanges.append(exchange)
        transfers.append(transfer)
    for vertex in reversed(vertices) if beam_reverse else vertices:
        event.add_vertex(vertex)
    production = hepmc3.GenVertex()
    for exchange in exchanges:
        production.add_particle_in(exchange)
    system = make_particle(*central, 90, 81)
    production.add_particle_out(system)
    event.add_vertex(production)
    decay = hepmc3.GenVertex()
    decay.add_particle_in(system)
    for pid, momentum in particles:
        decay.add_particle_out(make_particle(*momentum, pid, 1))
    event.add_vertex(decay)
    result = SimpleNamespace(evt=event, pid=[pid for pid, _ in particles], transfers=transfers, cut_param={})
    assert_conservation(result)
    return result


# Check every stored interaction or decay conserves four-momentum
def assert_conservation(record):
    for vertex in record.evt.vertices():
        incoming = sum((hepmc2vec4(p.momentum()) for p in vertex.particles_in()), vec4())
        outgoing = sum((hepmc2vec4(p.momentum()) for p in vertex.particles_out()), vec4())
        np.testing.assert_allclose(tuple(incoming), tuple(outgoing), rtol=1e-11, atol=1e-8)


# Construct an on-shell test particle from its Cartesian three-momentum
def on_shell(pid, px, py, pz):
    mass = {22: 0.0, 111: 0.134977, 211: 0.139570, 13: 0.105658, 2212: 0.938272, 2112: 0.939565}[abs(pid)]
    return vec4(px, py, pz, np.sqrt(px**2 + py**2 + pz**2 + mass**2))


# Include a central proton-antiproton pair that must not be confused with beam protons
def make_tagged_event():
    particles = [(13, on_shell(13, 12.0, 0.0, 30.0)), (-13, on_shell(-13, -11.0, 0.0, -30.0)),
                 (2212, on_shell(2212, 0.1, 0.0, 0.2)), (-2212, on_shell(-2212, -0.1, 0.0, -0.2))]
    event = make_cep_event(particles)
    event.pid = [-13, 13]
    return event


# Construct tagged dimuon photoproduction with a coherent or excited ion recoil
def make_incoherent_upc_event(target_leg, pt=0.2):
    mass, muon = 3.0969, 0.105658
    system = vec4(pt, 0.0, 0.0, np.hypot(mass, pt))
    particles = []
    for sign in (1, -1):
        momentum = vec4(sign * np.sqrt(mass**2 / 4.0 - muon**2), 0.0, 0.0, mass / 2.0)
        momentum.boost(b=system, sign=1)
        particles.append((sign * 13, momentum))
    masses = tuple(194.7 if leg == target_leg else 193.7 for leg in (1, 2))
    event = make_cep_event(particles, beam_pid=1000822080, beam_mass=193.7,
                           beam_energy=500000.0, forward_masses=masses)
    forwards = [p for p in event.evt.particles() if p.status() == 1 and p.pid() == 1000822080]
    for leg, particle in enumerate(forwards, start=1):
        if leg == target_leg:
            particle.set_pid(91)
        particle.add_attribute("graniitti_upc_leg", hepmc3.IntAttribute(leg))
        sector = "incoherent" if leg == target_leg else "coherent"
        particle.add_attribute("graniitti_upc_final_sector", hepmc3.StringAttribute(sector))
    event.cut_param = dict(alice_cuts.cut_param)
    return event


# Construct a coherent final state with two intact ion PDG identities
def make_coherent_upc_event():
    return make_incoherent_upc_event(target_leg=None)


# Include the recoil isotope in an excited Pb decay with a specified neutron multiplicity
def make_emd_event(xn_leg, multiplicity):
    event = hepmc3.GenEvent(hepmc3.Units.GEV, hepmc3.Units.MM)
    side = 1 if xn_leg == 1 else -1
    vertex = hepmc3.GenVertex()
    total = vec4()
    for _ in range(multiplicity):
        neutron = on_shell(2112, 0.0, 0.0, side * 10.0)
        total += neutron
        vertex.add_particle_out(make_particle(*neutron, 2112, 1))
    daughter_a = 208 - multiplicity
    pz = side * 10.0 * daughter_a
    daughter = vec4(0.0, 0.0, pz, np.hypot(pz, 0.9315 * daughter_a))
    total += daughter
    vertex.add_particle_out(make_particle(*daughter, 1000820000 + 10 * daughter_a, 1))
    vertex.add_particle_in(make_particle(*total, 1000822080, 2))
    event.add_vertex(vertex)
    result = SimpleNamespace(evt=event, pid=[2112], cut_param={"emd_side": side})
    assert_conservation(result)
    return result


# Construct soft CEP daughters and an optional charge-conserving forward p* -> n pi+ decay
def make_central_event(final_pids, requested_pids, forward_pid=None):
    particles = [(pid, on_shell(pid, (-1)**index * (0.2 + 0.1 * index), 0.0, 0.0))
                 for index, pid in enumerate(final_pids)]
    event = make_cep_event(particles, beam_energy=100.0,
                           forward_masses=(1.4, 0.938272) if forward_pid is not None else None)
    event.pid = requested_pids
    if forward_pid is not None:
        assert forward_pid == 211
        forward = next(p for p in event.evt.particles() if p.pid() == 2212 and p.status() == 1)
        forward.set_status(2)
        parent = hepmc2vec4(forward.momentum())
        mass, pion, neutron = parent.m, 0.139570, 0.939565
        energy = (mass**2 + pion**2 - neutron**2) / (2.0 * mass)
        p4 = vec4(np.sqrt(energy**2 - pion**2), 0.0, 0.0, energy)
        other = vec4(0.0, 0.0, 0.0, mass) - p4
        decay = hepmc3.GenVertex()
        decay.add_particle_in(forward)
        for pid, momentum in ((211, p4), (2112, other)):
            momentum.boost(b=parent, sign=1)
            decay.add_particle_out(make_particle(*momentum, pid, 1))
        event.evt.add_vertex(decay)
        assert_conservation(event)
    return event


# Verify central projection uses final-state identity without fiducial parameters
def test_central_no_cut_params():
    event = make_tagged_event()
    central = obs.proj_central_particles(event)

    assert len(central) == 2
    assert sum(central).pt == pytest.approx(1.0)
    assert obs.proj_central_system(event).pt == pytest.approx(1.0)
    assert obs.proj_1D_Pt(event) == pytest.approx(1.0)


# Check identical daughters retain finite Collins Soper angles under a common rotation
@pytest.mark.parametrize(("pid", "mass"), [(22, 0.0), (111, 0.135)])
def test_cs_identical_daughters(pid, mass):
    reference = None
    for angle in (0.0, 0.63):
        particles = []
        for px, py, pz in [(0.3, 0.1, 0.2), (-0.1, 0.2, -0.25)]:
            energy = np.sqrt(px * px + py * py + pz * pz + mass * mass)
            momentum = vec4(px * np.cos(angle) - py * np.sin(angle),
                            px * np.sin(angle) + py * np.cos(angle), pz, energy)
            particles.append((pid, momentum))
        event = make_cep_event(particles)
        frame = obs.proj_1D_CS(event)
        assert len(frame) == 2
        assert frame[0].costheta == pytest.approx(-frame[1].costheta, abs=1.0e-12)
        assert np.cos(frame[0].phi - frame[1].phi) == pytest.approx(-1.0, abs=1.0e-12)
        angles = (obs.proj_1D_costheta_CS(event), obs.proj_1D_phi_CS(event))
        assert np.isfinite(angles).all()
        if reference is not None:
            assert angles == pytest.approx(reference, abs=1.0e-10)
        reference = angles


# Retain the requirement for an unambiguous analyzer in a mixed four-body state
def test_cs_rejects_ambiguous_mixed_daughters():
    event = make_central_event([211, 211, -211, -211], [211, -211])
    with pytest.raises(ValueError, match="First requested PID"):
        obs.proj_1D_CS(event)


# Verify tagged protons are found through their beam vertices
def test_forward_projector_follows_beam_daughters():
    event = make_tagged_event()
    initial, final = obs.proj_init_final_protons(event)

    assert initial[0].z > 0.0
    assert initial[1].z < 0.0
    assert final[0].z > 0.0 and final[1].z < 0.0
    assert all(p.m == pytest.approx(0.938272, abs=1e-8) for p in [*initial, *final])


# Verify the momentum transfers use stable beam ordering and absolute histogram values
def test_forward_momentum_transfer_projectors():
    event = make_tagged_event()
    initial, final = obs.proj_init_final_protons(event)
    expected = tuple((initial[i] - final[i]).m2 for i in range(2))

    assert obs.proj_t_pair(event) == pytest.approx(expected)
    assert obs.proj_1D_t1(event) == pytest.approx(abs(expected[0]))
    assert obs.proj_1D_t2(event) == pytest.approx(abs(expected[1]))
    assert obs.proj_1D_Abs_t1t2(event) == pytest.approx(abs(sum(expected)))


# Verify tagged fragment daughters select exact target transfer in either direction
@pytest.mark.parametrize("target_leg", [1, 2])
def test_upc_incoherent_transfer_tagged_forward_leg(target_leg):
    event = make_incoherent_upc_event(target_leg)
    transfers = obs.proj_upc_forward_transfers(event)
    selected = transfers[target_leg - 1]

    assert selected["pid"] == 91
    assert selected["sector"] == "incoherent"
    assert obs.proj_upc_incoherent_abs_t(event) == pytest.approx(abs(event.transfers[target_leg - 1].m2))
    assert obs.proj_1D_Pt(event) ** 2 == pytest.approx(0.04)
    assert obs.proj_upc_incoherent_abs_t(event) != pytest.approx(
        obs.proj_1D_Pt(event) ** 2
    )


# Verify the ALICE measured bins and cuts use dimuon pT for either target direction
@pytest.mark.parametrize("target_leg", [1, 2])
@pytest.mark.parametrize(("pt", "accepted"), [(0.19, False), (0.21, True), (0.99, True), (1.01, False)])
def test_alice_incoherent_selection_measured_pt(target_leg, pt, accepted):
    event = make_incoherent_upc_event(target_leg, pt=pt)

    assert alice_obs.proj_abs_t(event) == pytest.approx(pt**2)
    assert bool(alice_cuts.cut_func(event)) is accepted


# Verify the measured coherent pT squared uses the central system rather than either ion recoil
def test_alice_coherent_pt2_central_system():
    event = make_coherent_upc_event()
    transfers = obs.proj_upc_forward_transfers(event)

    assert alice_coherent_obs.proj_pt2(event) == pytest.approx(0.04)
    assert alice_coherent_obs.proj_pt2(event) != pytest.approx(
        max(abs(record["t"]) for record in transfers)
    )


# Verify EMD multiplicity counts actual final neutrons on either detector side
@pytest.mark.parametrize(("xn_leg", "multiplicity"), [(1, 2), (2, 4)])
def test_alice_emd_multiplicity_final_particles(xn_leg, multiplicity):
    event = make_emd_event(xn_leg=xn_leg, multiplicity=multiplicity)

    assert alice_emd_obs.proj_neutron_multiplicity(event) == float(multiplicity)


# Verify the canonical observables support each soft CEP process multiplicity
@pytest.mark.parametrize(
    ("final_pids", "requested_pids"),
    [
        ([211, -211], [211, -211]),
        ([211, -211, 211, -211], [[211, -211], [211, -211]]),
        (
            [211, -211, 211, -211, 211, -211],
            [[211, -211], [211, -211], [211, -211]],
        ),
    ],
)
def test_softcep_multiplicities(final_pids, requested_pids):
    event = make_central_event(final_pids, requested_pids)
    values = (
        obs.proj_1D_M(event),
        obs.proj_1D_Rap(event),
        obs.proj_1D_Pt(event),
        obs.proj_1D_t1(event),
        obs.proj_1D_t2(event),
        obs.proj_1D_Abs_t1t2(event),
        obs.proj_1D_dPhi_pp(event),
    )

    assert all(np.isfinite(value) for value in values)


# Verify central proton identity does not select tagged beam protons
def test_central_ppbar_separation():
    event = make_central_event([2212, -2212], [2212, -2212])
    central = obs.proj_central_particle_records(event)

    assert [record["pid"] for record in central] == [2212, -2212]
    assert len(obs.proj_central_particles(event)) == 2


# Verify nested repeated identities select every central cascade leaf
def test_nested_pid_selection():
    event = make_central_event([211, -211, 211, -211], [[211, -211], [211, -211]], forward_pid=211)
    central = obs.proj_central_particle_records(event)

    assert [record["pid"] for record in central] == [211, -211, 211, -211]


# Build a two pion decay with independent HepMC record ordering and beam-axis rotation
def make_cs_event(reverse=False, rotation=0.0, requested=(211, -211), beam_reverse=False):
    pion_momenta = [(211, 0.7, -0.4, 0.8), (-211, -0.2, 0.6, -0.5)]
    particles = []
    for pid, px, py, pz in reversed(pion_momenta) if reverse else pion_momenta:
        momentum = on_shell(pid, px, py, pz)
        momentum.rotateZ(rotation)
        particles.append((pid, momentum))
    event = make_cep_event(particles, beam_reverse=beam_reverse)
    event.pid = list(requested)
    return event


# Check signed CS angles survive daughter permutation and rotations about the beams
def test_cs_angles_follow_requested_daughter():
    from core.analysis.observables.default import obs_phi_CS

    reference = make_cs_event()
    expected = [obs.proj_1D_costheta_CS(reference), obs.proj_1D_phi_CS(reference)]
    for reverse, rotation in [(True, 0.0), (False, 0.7), (True, -1.3)]:
        event = make_cs_event(reverse=reverse, rotation=rotation)
        np.testing.assert_allclose(
            [obs.proj_1D_costheta_CS(event), obs.proj_1D_phi_CS(event)], expected, atol=1.0e-10,
        )
    opposite = make_cs_event(requested=(-211, 211))
    assert obs.proj_1D_costheta_CS(opposite) == pytest.approx(-expected[0])
    angle = obs.proj_1D_phi_CS(opposite)
    assert abs(angle - expected[1]) == pytest.approx(180.0)
    assert min(angle, expected[1]) < 0.0 < max(angle, expected[1])
    entries, _ = np.histogram([angle, expected[1]], bins=obs_phi_CS["bins"])
    assert entries.sum() == 2
    exchanged = make_cs_event(beam_reverse=True)
    assert obs.proj_1D_costheta_CS(exchanged) == pytest.approx(-expected[0])


# Preserve the CS polar axis for collinear beams under rotations and beam exchange
@pytest.mark.parametrize("angle", [0.0, 0.7, 1.5707963267948966])
@pytest.mark.parametrize("reverse", [False, True])
def test_cs_collinear_polar_covariance(angle, reverse):
    from core.kinematics.vec4 import vec4

    beam1, beam2 = vec4(0.0, 0.0, 1.0, 1.0), vec4(0.0, 0.0, -1.0, 1.0)
    p = vec4(0.6, 0.0, 0.8, 1.1)
    for vector in (beam1, beam2, p):
        vector.rotateY(angle)
    if reverse:
        beam1, beam2 = beam2, beam1
    result = obs.LorentzFrame(beam1, beam2, [p])[0]
    assert result.costheta == pytest.approx(-0.8 if reverse else 0.8)
    assert result.m2 == pytest.approx(0.21)


# Normalize momentum directions without a dimensionful small-vector cutoff
@pytest.mark.parametrize("scale", [1.0e-200, 1.0, 1.0e200])
def test_unitvec_momentum_scale_invariance(scale):
    np.testing.assert_allclose(obs.unitvec(np.array([3.0, 0.0, 4.0]) * scale), [0.6, 0.0, 0.8])


# Reconstruct a vector daughter before the CS transform and preserve frame covariance
@pytest.mark.parametrize('reverse,rotation,rapidity', [(False, 0.0, 0.0), (True, 0.8, 0.6)])
def test_cs_composite_daughter(reverse, rotation, rapidity):
    from core.io.steering import select_histogram_observables

    particles = [(211, on_shell(211, 0.7, -0.4, 0.8)),
                 (111, on_shell(111, -0.1, 0.3, 0.2)),
                 (-211, on_shell(-211, -0.2, 0.1, -0.5))]
    reference = make_cep_event([(213, particles[0][1] + particles[1][1]), particles[2]])
    expected = [obs.proj_1D_costheta_CS(reference), obs.proj_1D_phi_CS(reference)]
    for _, momentum in particles:
        momentum.rotateZ(rotation)
    event = make_cep_event(list(reversed(particles)) if reverse else particles)
    for particle in event.evt.particles():
        momentum = hepmc2vec4(particle.momentum())
        momentum.boost(vec4(0.0, 0.0, np.sinh(rapidity), np.cosh(rapidity)), sign=1)
        particle.set_momentum(hepmc3.FourVector(*momentum))
    definitions = {name: {'func': function} for name, function in [
        ('cos', obs.proj_1D_costheta_CS), ('phi', obs.proj_1D_phi_CS)]}
    selected = select_histogram_observables({'hist': [
        {'obs': name, 'args': {'daughter': [211, 111]}} for name in definitions]}, definitions)
    result = [item['func'](event) for item in selected.values()]
    assert result == pytest.approx(expected, abs=1e-10)
    # Distinct daughter selections must not reuse the same cached angle
    assert obs.proj_1D_costheta_CS(event, daughter=[-211]) == pytest.approx(-result[0], abs=1e-10)
    assert obs.proj_1D_costheta_CS(event) != pytest.approx(result[0], abs=1e-5)
    assert definitions['cos']['func'] is obs.proj_1D_costheta_CS


# Reject missing and ambiguous composite daughter assignments
@pytest.mark.parametrize('daughter', [[], [211, 211], [321, 111]])
def test_cs_invalid_composite(daughter):
    event = make_cs_event()
    with pytest.raises(ValueError, match='Daughter PIDs'):
        obs.proj_1D_CS(event, daughter=daughter)

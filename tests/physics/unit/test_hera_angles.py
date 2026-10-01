# Test HERA helicity angles and fiducial selections on conserved HepMC3 events
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
import pytest
from core.analysis import obs
from core.io import readers, steering
from core.kinematics.vec4 import hepmc2vec4, vec4
from pyHepMC3 import HepMC3 as hepmc3

from icepack.PHOTOPROD._common import angles, dissociative, exclusive_jpsi
from tests.technical.data.test_hera_dissociation import photo_event
from tests.technical.support.hepdata import assert_original

ROOT = Path(__file__).resolve().parents[3]


# Build a conserved elastic ep event with two central leptons and a scattered positron
def lepton_event(pid, direction=1):
    mass = 0.000511 if pid == 11 else 0.105658
    event, _ = photo_event(direction, False, 0.938272, decay_pdg=-pid, decay_mass=mass)
    for particle in event.evt.particles():
        if particle.pid() == 90210:
            particle.set_pid(2212)
    event.pid = [-pid, pid]
    event.cut_param = {"W": [20.0, 80.0], "Q2_MAX": 2.5, "ABS_T_MAX": 1.5}
    return event


# Preserve the selected charge and both angles under rotations and arbitrary boosts
@pytest.mark.parametrize("pid", [11, 13])
@pytest.mark.parametrize("direction", [-1, 1])
def test_helicity_covariance(pid, direction):
    event = lepton_event(pid, direction)
    reference = np.array([angles.costheta(event), angles.phi(event)])
    assert exclusive_jpsi.accepted(event)
    pair = obs.proj_central_system(event)
    target = exclusive_jpsi.unique_momentum(event, 2212, False)
    daughter = next(p["p4"] for p in obs.proj_central_particle_records(event) if p["pid"] == -pid)
    a = daughter - pair * (daughter.dot4(pair) / pair.m2)
    b = target - pair * (target.dot4(pair) / pair.m2)
    assert reference[0] == pytest.approx(a.dot4(b) / np.sqrt(a.m2 * b.m2), abs=1e-10)
    for particle in event.evt.particles():
        momentum = hepmc2vec4(particle.momentum())
        momentum.rotateY(0.7)
        momentum.rotateZ(-0.4)
        momentum.boost(b=vec4(0.3, -0.2, 0.5, 1.4), sign=1)
        particle.set_momentum(hepmc3.FourVector(*momentum))
    event = readers.ice_event(event.evt, event.pid, event.cut_param)
    np.testing.assert_allclose([angles.costheta(event), angles.phi(event)], reference, atol=1e-9)
    assert exclusive_jpsi.accepted(event)
    event.pid = [pid, -pid]
    assert angles.costheta(event) == pytest.approx(-reference[0], abs=1e-9)
    assert angles.phi(event) == pytest.approx((reference[1] + np.pi) % (2 * np.pi), abs=1e-9)


# Reject each invariant cut independently without confusing central and beam leptons
@pytest.mark.parametrize("pid", [11, 13])
def test_elastic_cuts(pid):
    event = lepton_event(pid)
    values = exclusive_jpsi.kinematics(event)
    assert exclusive_jpsi.accepted(event)
    for cuts in [{"W": [values["w"] + 1, values["w"] + 2]}, {"W": [0, values["w"] - 1]},
                 {"Q2_MAX": 0.0}, {"ABS_T_MAX": values["abs_t"] / 2}]:
        original = event.cut_param
        event.cut_param = {**original, **cuts}
        assert not exclusive_jpsi.accepted(event)
        event.cut_param = original
    other = 13 if pid == 11 else 11
    event = readers.ice_event(event.evt, [-other, other], event.cut_param)
    assert not exclusive_jpsi.accepted(event)


# Reject each dissociative invariant including the reconstructed central and target masses
@pytest.mark.parametrize("fragment", [False, True])
def test_dissociative_cuts(fragment):
    event, _ = photo_event(1, fragment, 2.0)
    values = dissociative.kinematics(event)
    assert dissociative.accepted(event)
    for cuts in [{"W": [values["w"] + 1, values["w"] + 2]}, {"Q2_MAX": 0.0},
                 {"ABS_T_MAX": values["abs_t"] / 2}, {"M": [values["mass"] + 0.1, values["mass"] + 0.2]},
                 {"MY_MAX": np.sqrt(values["my2"]) / 2}]:
        original = event.cut_param
        event.cut_param = {**original, **cuts}
        assert not dissociative.accepted(event)
        event.cut_param = original


# Read the four original angular histograms through the complete icepack reader
@pytest.mark.parametrize("channel", [0, 1])
def test_angular_data(channel):
    card = ROOT / 'icepack/PHOTOPROD/ZEUS_582237/angles/dataset.json'
    dataset, path = steering.load_dataset(str(card), cdir=ROOT)
    entry = dataset['sets'][channel]
    observables, _ = steering.load_observables(entry['obs'], dataset_path=path, cdir=ROOT)
    data, _ = readers.read_hepdata(dataset=entry, datapath=dataset['datapath'], datatype=dataset['type'],
                                 all_obs=observables, cdir=str(ROOT), reader=dataset['reader'], dataset_path=path)
    for histogram in entry['hist']:
        assert_original(data[histogram['obs']], ROOT / dataset['datapath'] / histogram['file'])

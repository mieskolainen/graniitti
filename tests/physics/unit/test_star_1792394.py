# STAR 200 GeV dataset and fiducial-selection configuration tests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import math
from pathlib import Path

import numpy as np
import pytest
from core.analysis import obs
from core.io import readers, steering
from core.kinematics.vec4 import vec4

from icepack.SOFTCEP.STAR_1792394._common import observables, projectors, selections
from tests.technical.support.hepmc import make_pair_event


# Verify STAR observables select the central pair instead of tagged beam protons
def test_ppbar_obs_fiducial_central_pair():
    event = make_pair_event()
    pair = projectors.central_particles(event)

    assert len(pair) == 2
    assert projectors.invariant_mass(event) == pytest.approx(2.264)
    assert projectors.rapidity(event) == pytest.approx(0.0)
    assert projectors.collins_soper_costheta(event) == pytest.approx(0.2, rel=0.0, abs=1e-12)
    assert projectors.collins_soper_phi(event) == pytest.approx(0.0, abs=1e-12)
    for vertex in event.evt.vertices():
        incoming = sum((vec4(p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().e())
                        for p in vertex.particles_in()), vec4())
        outgoing = sum((vec4(p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().e())
                        for p in vertex.particles_out()), vec4())
        assert tuple(incoming - outgoing) == pytest.approx((0.0, 0.0, 0.0, 0.0), abs=1e-12)


# Check mass partitions using the actual invariant mass of an on-shell pion pair
def test_star_mass_selection_boundaries():
    lower = make_pair_event(mass=1.0, pid=211)
    upper = make_pair_event(mass=1.5, pid=211)
    above = make_pair_event(mass=1.500001, pid=211)
    assert selections.mass_at_most(lower, 1.0)
    assert not selections.mass_interval(lower, 1.0, 1.5)
    assert selections.mass_interval(upper, 1.0, 1.5)
    assert not selections.mass_above(upper, 1.5)
    assert selections.mass_above(above, 1.5)


# Check azimuth partitions using tagged outgoing protons in conserving events
def test_star_azimuth_selection_boundaries():
    below = make_pair_event(dphi=89.999999)
    boundary = make_pair_event(dphi=90.0)
    above = make_pair_event(dphi=90.000001)
    assert selections.azimuth_below(below, 90.0)
    assert not selections.azimuth_below(boundary, 90.0)
    assert not selections.azimuth_above(boundary, 90.0)
    assert selections.azimuth_above(above, 90.0)


# For a transverse recoil along minus x the Collins-Soper axes coincide with the pair rest axes
@pytest.mark.parametrize("phi", [-179.0, -70.0, 0.0, 70.0, 179.0])
def test_collins_soper_azimuth_sign(phi):
    event = make_pair_event(dphi=0.0, phi=phi)
    observable = observables.collins_soper_azimuthal_angle()
    assert observable["func"](event) == pytest.approx(phi, rel=0.0, abs=1e-10)
    assert projectors.collins_soper_costheta(event) == pytest.approx(0.2, rel=0.0, abs=1e-12)


# Apply each Figure 13 card selection to conserving events on both sides of its cut
@pytest.mark.parametrize("panel,above", [("left", False), ("right", True)])
def test_proton_dpt_selection(panel, above):
    root = Path(__file__).resolve().parents[3]
    dataset, path = steering.load_dataset("icepack/SOFTCEP/STAR_1792394/pipi/dataset.json", cdir=root)
    entry = next(entry for entry in dataset["sets"]
                 if any(hist["file"] == f"Figure13({panel}).json" for hist in entry["hist"]))
    selection = readers.load_cut_module(steering.resolve_python_reference(
        entry["cuts"], package="core.analysis.cuts", dataset_path=path, cdir=root,
    ))
    params = selection.cut_param
    assert params["above"] is above
    phi = params["dPhi"] / 2
    for ratio in (0.5, 1.5):
        dpt = ratio * params["dPt"]
        pt = dpt / (2 * math.sin(math.radians(phi) / 2))
        event = make_pair_event(pid=211, dphi=phi, forward_pt=pt)
        event.cut_param = params
        assert bool(selection.cut_func(event)) is (above if ratio > 1 else not above)
        _, protons = obs.proj_init_final_protons(event)
        assert (protons[0] - protons[1]).pt == pytest.approx(dpt)
        assert protons[0].pt == pytest.approx(protons[1].pt)
    for phi in (params["dPhi"], params["dPhi"] * 1.5):
        event = make_pair_event(pid=211, dphi=phi)
        event.cut_param = params
        assert not selection.cut_func(event)
    event = make_pair_event(pid=211, dphi=params["dPhi"] / 2)
    _, protons = obs.proj_init_final_protons(event)
    boundary = (protons[0] - protons[1]).pt
    for limit, expected in ((np.nextafter(boundary, -np.inf), above), (boundary, False),
                            (np.nextafter(boundary, np.inf), not above)):
        event.cut_param = {**params, "dPt": limit}
        assert bool(selection.cut_func(event)) is expected

# Physics and runtime checks for the Pandora icetune driver
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import ast
import copy
import json
import math
import subprocess
import sys
import tarfile
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyjson5
import pytest
from core.tune import push as icetune_push
from core.tune.drivers.pandora import diagnostics, plots, runtime
from core.tune.drivers.pandora import driver as pandora
from core.tune.drivers.pandora import key4hep as pandora_setup
from core.tune.drivers.pandora.tunesetup import common
from core.tune.summary import build_best_fit
from core.tune.tunesetup import load_tunesetup, save_tunesetup

from submit import campaign_source

ak = pytest.importorskip("awkward")
uproot = pytest.importorskip("uproot")
ROOT = Path(__file__).resolve().parents[3]


# Load the real parameter catalog and steering surface once for this module
@pytest.fixture(scope="module")
def tunesetup():
    return load_tunesetup(cdir=ROOT, simdriver="PANDORA", name=campaign_source("tune-pandora-v0"))


# Isolate objective variations from the production steering
@pytest.fixture
def spec(tunesetup):
    return copy.deepcopy(tunesetup.datacards[0]["objectives"]["pflow"])


# Build a truth-particle record for PF-set objective fixtures
def _truth_particle(pdg, charge, px, py, pz, mass=0.0, status=1):
    return {
        "pdg": int(pdg),
        "status": int(status),
        "charge": float(charge),
        "mass": float(mass),
        "px": float(px),
        "py": float(py),
        "pz": float(pz),
    }


# Compute the relativistic energy for a test particle record
def _particle_energy(particle):
    return float(
        np.sqrt(particle["px"] ** 2 + particle["py"] ** 2 + particle["pz"] ** 2 + particle.get("mass", 0.0) ** 2)
    )


# Build a reconstructed-particle record for PF-set objective fixtures
def _reco_particle(pdg, px, py, pz, energy=None):
    return {
        "pdg": int(pdg),
        "energy": float(np.sqrt(px * px + py * py + pz * pz) if energy is None else energy),
        "px": float(px),
        "py": float(py),
        "pz": float(pz),
    }


# Compute the baseline visible stable truth particles for PF-set tests
def _base_truth_particles():
    return [
        _truth_particle(22, 0.0, 10.0, 0.0, 0.0),
        _truth_particle(211, 1.0, 0.0, 8.0, 0.0, mass=0.139),
        _truth_particle(130, 0.0, -6.0, 0.0, 0.0, mass=0.498),
    ]


# Compute the baseline reconstructed particles for PF-set tests
def _base_reco_particles(scale=1.0):
    truth = _base_truth_particles()
    return [
        _reco_particle(22, 10.0, 0.0, 0.0, energy=scale * _particle_energy(truth[0])),
        _reco_particle(211, 0.0, 8.0, 0.0, energy=scale * _particle_energy(truth[1])),
        _reco_particle(130, -6.0, 0.0, 0.0, energy=scale * _particle_energy(truth[2])),
    ]


# Build EDM4hep arrays with hard partons in the same event as stable particles
def _pf_arrays(events, parton_costheta=None):
    angles = [0.0] * len(events) if parton_costheta is None else parton_costheta
    truth = []
    for (particles, _), angle in zip(events, angles, strict=True):
        parton = (
            [] if angle is None else [_truth_particle(2, 2 / 3, 50 * math.sqrt(1 - angle**2), 0, 50 * angle, status=23)]
        )
        truth.append(particles + parton)
    arrays = {}
    for collection, rows in (("MCParticles", truth), ("PandoraPFANewPFOs", [r for _, r in events])):
        fields = dict(pdg="PDG", px="momentum.x", py="momentum.y", pz="momentum.z")
        fields.update(
            dict(status="generatorStatus", charge="charge", mass="mass")
            if collection == "MCParticles"
            else dict(energy="energy")
        )
        for key, branch in fields.items():
            # Explicit jagged types also support an entirely empty reconstructed collection
            offsets = np.r_[0, np.cumsum([len(row) for row in rows])]
            values = np.array([p[key] for row in rows for p in row], dtype=int if key in ("pdg", "status") else float)
            arrays[f"{collection}/{collection}.{branch}"] = ak.Array(
                ak.contents.ListOffsetArray(ak.index.Index64(offsets), ak.contents.NumpyArray(values))
            )
    return arrays


# Compute a massless photon event with a controllable reconstructed energy response
def _photons(scale=1.0):
    momenta = [(10, 0, 0), (-10, 0, 0), (0, 10, 0), (0, -10, 0)]
    return (
        [_truth_particle(22, 0, *p) for p in momenta],
        [_reco_particle(22, *(scale * np.asarray(p))) for p in momenta],
    )


# Test an independently calculable shortest interval, including weighted event selection
@pytest.mark.parametrize(
    "values,weights,fraction,mean,sigma",
    [
        ([0, 1, 2, 3, 4, 5, 6, 7, 8, 1000], None, 0.9, 4, np.std(np.arange(9))),
        ([0, 1, 2, 10, 11, 12], None, 0.5, 1, np.std([0, 1, 2])),
        ([0, 1, 100], [1, 9, 0], 0.9, 1, 0),
        ([0, 1], [1, 3], 1, 0.75, math.sqrt(0.1875)),
    ],
)
def test_central_interval(values, weights, fraction, mean, sigma):
    result = diagnostics.interval_stats(values, weights=weights, fraction=fraction)
    assert result["mean"] == pytest.approx(mean)
    assert result["sigma"] == pytest.approx(sigma)


# The positive loss must reject empty events, energy loss, misclassification and extra particles
@pytest.mark.parametrize("failure", ["empty", "lost", "fake", "wrong", "response", "compensating"])
def test_loss_rejects_failures(spec, failure):
    truth, reco = _photons()
    perfect, detail = pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))
    assert perfect["pflow"] == pytest.approx(0.0, abs=1e-20)
    if failure == "empty":
        reco = []
    elif failure == "lost":
        reco = reco[1:]
    elif failure == "fake":
        reco = [*reco, _reco_particle(22, 3, 3, 0)]
    elif failure == "wrong":
        reco = [{**p, "pdg": 211} for p in reco]
    elif failure == "response":
        _, reco = _photons(0.1)
    else:
        reco = [
            _reco_particle(22, *(np.asarray([p["px"], p["py"], p["pz"]]) * scale))
            for p, scale in zip(reco, [1.9, 1.9, 0.1, 0.1], strict=True)
        ]
    bad, _ = pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))
    spread = spec["energy_bias_tolerance"] / 2
    ordinary, _ = pandora._objective_pflow(spec, _pf_arrays([_photons(1 - spread), _photons(1 + spread)]))
    assert bad["pflow"] > ordinary["pflow"] > perfect["pflow"]
    if failure == "compensating":
        assert bad["pf_four_momentum_closure_mse"] < 1e-8
        assert bad["pf_loss_response"] == pytest.approx(
            0.9**2 * (1 / 2.9 + 1 / 1.1) / spec["matched_energy_response_target"]["photon"] ** 2
        )
    assert not any("nll" in key or "fisher" in key.lower() for key in bad)
    np.testing.assert_allclose(detail["confusion"]["matrix"], detail["confusion"]["ideal"], atol=1e-14)


# Check the response scale has its declared quadratic effect on the minimized loss
def test_response_target(spec):
    truth, reco = _photons()
    reco[:2] = [_reco_particle(22, 8, 0, 0), _reco_particle(22, -12, 0, 0)]
    arrays = _pf_arrays([(truth, reco)])
    first, _ = pandora._objective_pflow(spec, arrays)
    spec["matched_energy_response_target"]["photon"] *= 2
    second, _ = pandora._objective_pflow(spec, arrays)
    assert first["pf_loss_response"] > 0
    assert second["pf_loss_response"] == pytest.approx(first["pf_loss_response"] / 4)
    assert second["pf_loss_resolution"] == pytest.approx(first["pf_loss_resolution"])
    assert second["pf_loss_momentum"] == pytest.approx(first["pf_loss_momentum"])
    assert first["pflow"] > second["pflow"]


# PID and resolution must remain active on both sides of the calibration tolerance
@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("factor", [0.5, 2.0])
def test_unified_biased_response(spec, sign, factor):
    mean = 1 + sign * factor * spec["energy_bias_tolerance"]
    truth, reco = _photons(mean)
    correct, _ = pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))
    wrong, _ = pandora._objective_pflow(spec, _pf_arrays([(truth, [{**p, "pdg": 211} for p in reco])]))
    broad, _ = pandora._objective_pflow(spec, _pf_arrays([_photons(mean - 0.1), _photons(mean + 0.1)]))
    assert correct["pf_loss_calibration"] == pytest.approx(((mean - 1) / spec["energy_bias_target"])**2)
    assert wrong["pf_loss_calibration"] == pytest.approx(correct["pf_loss_calibration"])
    assert broad["pf_loss_calibration"] == pytest.approx(correct["pf_loss_calibration"])
    assert wrong["pflow"] > correct["pflow"]
    assert broad["pflow"] > correct["pflow"]


# The calibration interval includes its boundaries and rejects shifts beyond either side
@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("factor,passed", [(0.99, 1), (1.0, 1), (1.01, 0)])
def test_calibration_boundary(spec, sign, factor, passed):
    spec["energy_bias_tolerance"] = 0.125
    metrics, _ = pandora._objective_pflow(
        spec, _pf_arrays([_photons(1 + sign * factor * spec["energy_bias_tolerance"])]))
    assert metrics["pf_calibration_passed"] == passed
    expected = (factor * spec["energy_bias_tolerance"] / spec["energy_bias_target"])**2
    assert metrics["pf_loss_calibration"] == pytest.approx(expected)
    assert metrics["pflow"] == pytest.approx(expected)


# The calibration target controls the loss while the tolerance only controls its diagnostic
def test_calibration_steering(spec):
    arrays = _pf_arrays([_photons(1 + 2 * spec["energy_bias_tolerance"])])
    first, _ = pandora._objective_pflow(spec, arrays)
    spec["energy_bias_target"] *= 2
    second, _ = pandora._objective_pflow(spec, arrays)
    assert second["pflow"] == pytest.approx(first["pflow"] / 4)
    assert first["pf_calibration_passed"] == second["pf_calibration_passed"] == 0
    spec["energy_bias_tolerance"] *= 3
    third, _ = pandora._objective_pflow(spec, arrays)
    assert third["pflow"] == pytest.approx(second["pflow"])
    assert third["pf_calibration_passed"] == 1


# Calibration stays quadratic while the summed PF risk enters through a smooth logarithm
def test_pf_components(spec):
    narrow, _ = pandora._objective_pflow(spec, _pf_arrays([_photons(0.95), _photons(1.05)]))
    broad, _ = pandora._objective_pflow(spec, _pf_arrays([_photons(0.8), _photons(1.2)]))
    truth, reco = _photons()
    wrong = [{**p, "pdg": 211} for p in reco]
    misidentified, _ = pandora._objective_pflow(spec, _pf_arrays([(truth, wrong)]))
    assert narrow["pf_loss_resolution"] == pytest.approx(0.05**2 / spec["energy_resolution_target"]**2)
    assert 0 < narrow["pflow"] < broad["pflow"]
    assert misidentified["pf_loss_resolution"] == pytest.approx(0)
    assert misidentified["pf_loss_response"] == pytest.approx(0)
    assert misidentified["pf_loss_confusion"] > 0
    for metrics in (narrow, broad, misidentified):
        combined = sum(metrics[f"pf_loss_{name}"] for name in ("calibration", "resolution", "confusion", "response", "momentum", "multiplicity"))
        calibration = metrics["pf_loss_calibration"]
        assert metrics["pflow"] == pytest.approx(calibration + math.log1p(combined - calibration))


# Calibration must use every event even when the central response interval is exactly at one
def test_calibration_full_mean(spec):
    metrics, _ = pandora._objective_pflow(spec, _pf_arrays([*[_photons()] * 9, _photons(0.5)]))
    assert metrics["pf_total_energy_response_mean90"] == pytest.approx(1)
    assert metrics["pf_total_energy_response_mean"] == pytest.approx(0.95)
    assert metrics["pf_calibration_passed"] == 0
    assert metrics["pf_loss_calibration"] == pytest.approx((0.05 / spec["energy_bias_target"])**2)


# Empty PFO collections must remain in the response mean and variance
def test_calibration_empty_reconstruction(spec):
    truth, _ = _photons()
    metrics, _ = pandora._objective_pflow(spec, _pf_arrays([(truth, []), _photons(2)]))
    assert metrics["pf_n_events"] == 2
    assert metrics["pf_total_energy_response_mean"] == pytest.approx(1)
    assert metrics["pf_total_energy_residual_variance"] == pytest.approx(1)
    assert metrics["pf_calibration_passed"] == 1


# Calibration uses physical energies at the fiducial threshold
def test_calibration_physical_energy(spec):
    truth = [_truth_particle(22, 0, 0.2, 0, 0)]
    arrays = _pf_arrays([(truth, [_reco_particle(22, energy, 0, 0)]) for energy in (0.1, 0.3)])
    metrics, detail = pandora._objective_pflow(spec, arrays)
    assert metrics["pf_total_energy_response_mean"] == pytest.approx(1)
    assert metrics["pf_total_energy_residual_variance"] == pytest.approx(0.25)
    assert metrics["pf_calibration_passed"] == 1
    np.testing.assert_allclose(detail["truth_energy_values"], [0.2, 0.2])
    np.testing.assert_allclose(detail["visible_total_energy_values"], [0.1, 0.3])


# Check physical energy fractions for soft lost, misidentified and extra particles
@pytest.mark.parametrize("failure", ["lost", "wrong", "fake"])
@pytest.mark.parametrize("scale", [1, 2])
def test_soft_confusion(spec, failure, scale):
    spec.update(gen_final_state_min_energy=0.1, reco_final_state_min_energy=0.1)
    spec["confusion_energy_fraction_target"] *= scale
    truth = [_truth_particle(22, 0, 0.9, 0, 0), _truth_particle(22, 0, 0, 0.1, 0)]
    reco = [_reco_particle(22, 0.9, 0, 0), _reco_particle(22, 0, 0.1, 0)]
    if failure == "lost":
        reco.pop()
    elif failure == "wrong":
        reco[-1]["pdg"] = 2112
    else:
        reco.append(_reco_particle(2112, -0.1, 0, 0))
    metrics, detail = pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))
    matrix = detail["confusion"]["matrix"]
    expected = np.zeros_like(matrix)
    photon, neutral = (diagnostics.PF_RESIDUAL_CLASS_NAMES.index(c) for c in ("photon", "neutral_hadron"))
    expected[photon, photon] = 1.0 if failure == "fake" else 0.9
    expected[-1 if failure == "lost" else neutral, -1 if failure == "fake" else photon] = 0.1
    np.testing.assert_allclose(matrix, expected, atol=1e-14)
    fraction, entries = (1 / 11, 1) if failure == "fake" else (0.1, 2)
    loss = entries * (fraction / spec["confusion_energy_fraction_target"])**2
    assert metrics["pf_loss_confusion"] == pytest.approx(loss)


# Equal extra energy must receive a larger count penalty when split into more accepted neutral PFOs
@pytest.mark.parametrize("count", [1, 4])
def test_neutral_multiplicity(spec, count):
    truth, reco = _photons(0.95)
    energy, mass = 2.0 / count, 0.498
    component = math.sqrt(energy**2 - mass**2) / math.sqrt(2)
    reco += [_reco_particle(130, component, component, 0, energy) for _ in range(count)]
    arrays = _pf_arrays([(truth, reco)])
    metrics, _ = pandora._objective_pflow(spec, arrays)
    assert metrics["pf_n_fake"] == count
    assert metrics["pf_total_energy_response_mean"] == pytest.approx(1)
    expected = (count / spec["multiplicity_target"])**2
    assert metrics["pf_loss_multiplicity"] == pytest.approx(expected)
    energy_loss = sum(metrics[f"pf_loss_{name}"] for name in ("resolution", "confusion", "response", "momentum"))
    assert metrics["pflow"] == pytest.approx(math.log1p(energy_loss + expected))
    spec["multiplicity_target"] *= 2
    rescaled, _ = pandora._objective_pflow(spec, arrays)
    assert rescaled["pf_loss_multiplicity"] == pytest.approx(expected / 4)
    assert rescaled["pflow"] < metrics["pflow"]


# Opposite count errors cannot cancel between classes or between events
@pytest.mark.parametrize("variation", ["classes", "events"])
def test_count_cancellation(spec, variation):
    truth, reco = _photons()
    if variation == "classes":
        reco[0]["pdg"] = 130
        events = [(truth, reco)]
        expected = 2 / spec["multiplicity_target"]**2
    else:
        extra = _reco_particle(22, 1, 1, 0)
        events = [(truth, reco + [extra]), (truth, reco[1:])]
        expected = 1 / spec["multiplicity_target"]**2
    metrics, _ = pandora._objective_pflow(spec, _pf_arrays(events))
    assert metrics["pf_loss_multiplicity"] == pytest.approx(expected)


# An unresolved soft particle is excluded, while an accepted collinear split changes PFO multiplicity
@pytest.mark.parametrize("change", ["soft", "split"])
def test_resolved_particles(spec, change):
    truth, reco = _photons()
    if change == "soft":
        energy = spec["reco_final_state_min_energy"] / 2
        reco.append(_reco_particle(22, energy, 0, 0))
    else:
        reco = [_reco_particle(22, 4, 0, 0), _reco_particle(22, 6, 0, 0), *reco[1:]]
    metrics, _ = pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))
    expected = (1 if change == "split" else 0) / spec["multiplicity_target"]**2
    assert metrics["pf_loss_multiplicity"] == pytest.approx(expected)
    assert metrics["pf_total_energy_response_mean"] == pytest.approx(1)
    if change == "soft":
        assert metrics["pflow"] == pytest.approx(0, abs=1e-20)


# Redistributing energy between soft and hard matched photons has a finite analytic cost
def test_soft_response(spec):
    spec.update(gen_final_state_min_energy=0.1, reco_final_state_min_energy=0.1)
    truth = [_truth_particle(22, 0, 0.1, 0, 0), _truth_particle(22, 0, -0.9, 0, 0)]
    reco = [_reco_particle(22, 0.2, 0, 0), _reco_particle(22, -0.8, 0, 0)]
    metrics, _ = pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))
    residual = 2 * 0.1**2 * (1 / 0.3 + 1 / 1.7)
    assert metrics["pf_total_energy_response_mean"] == pytest.approx(1)
    assert metrics["pf_loss_confusion"] == pytest.approx(0, abs=1e-14)
    assert metrics["pf_loss_response"] == pytest.approx(residual / spec["matched_energy_response_target"]["photon"]**2)
    assert metrics["pf_loss_momentum"] == pytest.approx(residual / spec["matched_momentum_response_target"]**2)


# Large PF errors must remain finite and distinguishable without bounded score saturation
def test_large_loss(spec):
    truth, reco = _photons()
    reco[:2] = [_reco_particle(22, 8, 0, 0), _reco_particle(22, -12, 0, 0)]
    arrays = _pf_arrays([(truth, reco)])
    spec["matched_energy_response_target"]["photon"] = 1e-12
    metrics, _ = pandora._objective_pflow(spec, arrays)
    assert metrics["pf_loss_response"] > 1 / np.finfo(float).eps
    assert metrics["pf_calibration_passed"] == 1
    assert math.isfinite(metrics["pflow"])
    spec["matched_energy_response_target"]["photon"] /= 2
    worse, _ = pandora._objective_pflow(spec, arrays)
    assert worse["pflow"] - metrics["pflow"] == pytest.approx(math.log(4))


# Matched momentum errors must survive cancellation of the total event momentum
@pytest.mark.parametrize("change", ["angle", "magnitude"])
def test_matched_momentum(spec, change):
    truth, reco = _photons()
    angle = spec["matching_max_angle_rad"] / 2
    scale = 0.8
    for particle in reco:
        px, py = particle["px"], particle["py"]
        particle["px"], particle["py"] = (
            (px * math.cos(angle) - py * math.sin(angle), px * math.sin(angle) + py * math.cos(angle))
            if change == "angle" else (scale * px, scale * py)
        )
    metrics, _ = pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))
    residual = 4 * math.sin(angle / 2)**2 if change == "angle" else (1 - scale)**2
    expected = residual / spec["matched_momentum_response_target"]**2
    assert metrics["pf_three_momentum_closure_mse"] == pytest.approx(0, abs=1e-20)
    assert metrics["pf_loss_response"] == pytest.approx(0, abs=1e-20)
    assert metrics["pf_loss_momentum"] == pytest.approx(expected)
    assert metrics["pflow"] == pytest.approx(math.log1p(expected))
    spec["matched_momentum_response_target"] *= 2
    rescaled, _ = pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))
    assert rescaled["pf_loss_momentum"] == pytest.approx(expected / 4)
    assert rescaled["pflow"] < metrics["pflow"]


# A substantial mean bias must remain costly even with perfect PID and zero response width
@pytest.mark.parametrize("sign", [-1, 1])
def test_calibration_priority(spec, sign):
    events = []
    for response in (0.5, 1.5):
        truth, reco = _photons(response)
        events.append((truth, [{**p, "pdg": 211} for p in reco]))
    calibrated, _ = pandora._objective_pflow(spec, _pf_arrays(events))
    biased, _ = pandora._objective_pflow(
        spec, _pf_arrays([_photons(1 + sign * 3 * spec["energy_bias_tolerance"])]))
    assert calibrated["pflow"] < biased["pflow"]
    assert biased["pflow"] == pytest.approx((3 * spec["energy_bias_tolerance"] / spec["energy_bias_target"])**2)
    assert calibrated["pf_calibration_passed"] == 1
    assert biased["pf_calibration_passed"] == 0


# Multiplicative response shrinkage can lower MSE but must not win the PF fit
def test_calibration_shrinkage(spec):
    candidates = []
    for scale in (0.8, 0.99, 1.0, 1.01, 1.2):
        metrics, _ = pandora._objective_pflow(spec, _pf_arrays([_photons(scale * r) for r in (0.5, 1.5)]))
        candidates.append(metrics)
        assert metrics["pf_loss_resolution"] == pytest.approx(0.25 / spec["energy_resolution_target"]**2)
    assert candidates[0]["pf_total_energy_response_mse"] < candidates[2]["pf_total_energy_response_mse"]
    best = min(candidates, key=lambda metrics: metrics["pflow"])
    assert best["pf_calibration_passed"] == 1
    assert best["pf_total_energy_response_mean"] == pytest.approx(1)


# Optimize a continuous gain instead of testing only a few trial energies
@pytest.mark.parametrize("responses", [(0.5, 1.5), (0.4, 0.4, 2.2), (*[0.8] * 9, 2.8), (0.0, 0.5, 2.5)])
def test_gain_minimum(spec, responses):
    from scipy.optimize import minimize_scalar

    spec["_collect_plot_details"] = False

    # Evaluate actual visible four momenta throughout the gain fit
    def loss(gain):
        metrics, _ = pandora._objective_pflow(spec, _pf_arrays([_photons(gain * r) for r in responses]))
        return metrics["pflow"]

    result = minimize_scalar(loss, bounds=(0.6, 1.4), method="bounded", options={"xatol": 1e-10})
    assert result.success
    assert result.x * np.mean(responses) == pytest.approx(1, abs=1e-7)


# Resolution, energy sharing, momentum and fake PID penalties cannot reward a uniform gain change
@pytest.mark.parametrize("fake_energy", [0.3, 2.0])
def test_gain_invariance(spec, fake_energy):
    events = []
    for response in (0.7, 1.3):
        truth, reco = _photons(response)
        reco[0] = _reco_particle(211, 8 * response, 0, 0)
        reco = [*reco[:-1], _reco_particle(22, fake_energy * response, fake_energy * response, 0)]
        events.append((truth, reco))
    results = []
    for gain in (0.6, 1, 1.4):
        scaled = [(truth, [{key: value if key == "pdg" else gain * value for key, value in p.items()}
                            for p in reco]) for truth, reco in events]
        metrics, _ = pandora._objective_pflow(spec, _pf_arrays(scaled))
        results.append(metrics)
    for term in ("resolution", "confusion", "response", "momentum", "multiplicity"):
        values = [m[f"pf_loss_{term}"] for m in results]
        assert values[1] > 0
        np.testing.assert_allclose(values, values[1], rtol=1e-10, atol=1e-12)


# The calibration tolerance is diagnostic and cannot create a jump or kink in the optimized cost
@pytest.mark.parametrize("sign", [-1, 1])
def test_loss_smoothness(spec, sign):
    tolerance = spec["energy_bias_tolerance"]
    step = tolerance * 1e-4
    center = 1 + sign * tolerance
    samples = []
    for gain in (center - step, center, center + step):
        metrics, _ = pandora._objective_pflow(spec, _pf_arrays([_photons(gain)]))
        samples.append(metrics["pflow"])
    left, middle, right = samples
    assert (right - left) / (2 * step) == pytest.approx(2 * sign * tolerance / spec["energy_bias_target"]**2, rel=1e-6)
    assert (right - 2 * middle + left) / step**2 == pytest.approx(2 / spec["energy_bias_target"]**2, rel=1e-5)


# Reconstruct invariant energies, fiducial angles and stable visible particles together
def test_visible_reconstruction(spec):
    spec.update(gen_final_state_max_energy=1.0, reco_final_state_max_energy=1.0)
    truth = [
        _truth_particle(22, 0, 0.05, 0, 0),
        _truth_particle(22, 0, 0.1, 0, 0),
        _truth_particle(130, 0, 0, 0.3, 0, mass=0.4),
        _truth_particle(22, 0, -2, 0, 0),
        _truth_particle(22, 0, 0.1, 0, 0.5),
        _truth_particle(12, 0, 0.4, 0, 0),
    ]
    reco = [_reco_particle(p["pdg"], p["px"], p["py"], p["pz"], _particle_energy(p)) for p in truth]
    metrics, detail = pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))
    assert metrics["pf_n_truth_visible"] == metrics["pf_n_reco_visible"] == 2
    np.testing.assert_allclose(detail["truth_energy_values"], [0.6])
    np.testing.assert_allclose(detail["visible_total_energy_values"], [0.6])
    assert metrics["pflow"] == pytest.approx(0, abs=1e-20)


# Check azimuthal rotations, beam reflection and particle permutations preserve the scalar objective
@pytest.mark.parametrize("angle,reflection", [(0.7, 1), (2.1, -1)])
def test_covariance_and_permutations(spec, angle, reflection):
    truth, reco = _base_truth_particles(), _base_reco_particles(1.03)
    transform = np.array(
        [[math.cos(angle), -math.sin(angle), 0], [math.sin(angle), math.cos(angle), 0], [0, 0, reflection]]
    )
    first, _ = pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))
    for collection in (truth, reco):
        for particle in collection:
            momentum = transform @ np.array([particle["px"], particle["py"], particle["pz"]])
            particle.update(zip(("px", "py", "pz"), momentum, strict=True))
        collection.reverse()
    second, _ = pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))
    for key in ("pflow", "pf_four_momentum_closure_mse", "pf_loss_response", "pf_loss_momentum", "pf_loss_confusion"):
        assert second[key] == pytest.approx(first[key], abs=1e-12)


# Maximize the number of angular matches before minimizing the sum of allowed angles
def test_matching_cardinality():
    truth_angles, reco_angles = np.array([0, 0.19]), np.array([0.1, -0.19])
    truth = np.column_stack((np.cos(truth_angles), np.sin(truth_angles), np.zeros(2)))
    reco = np.column_stack((np.cos(reco_angles), np.sin(reco_angles), np.zeros(2)))
    truth, reco = (dict(unit=unit, energy=np.ones(len(unit))) for unit in (truth, reco))
    pairs = pandora._pf_angular_matched_pairs(truth, reco, 0.2)
    assert set(map(tuple, pairs)) == {(0, 1), (1, 0)}
    empty = dict(unit=np.empty((0, 3)), energy=np.empty(0))
    assert pandora._pf_angular_matched_pairs(truth, empty, 0.2).shape == (0, 2)


# A perfect collinear particle set must retain zero loss under independent collection permutations
@pytest.mark.parametrize("truth_order,reco_order", [(1, 1), (1, -1), (-1, 1), (-1, -1)])
def test_collinear_permutations(spec, truth_order, reco_order):
    truth = [_truth_particle(22, 0, 10, 0, 0), _truth_particle(211, 1, 20, 0, 0, mass=0.139)]
    reco = [_reco_particle(p["pdg"], p["px"], p["py"], p["pz"], _particle_energy(p)) for p in truth]
    metrics, _ = pandora._objective_pflow(spec, _pf_arrays([(truth[::truth_order], reco[::reco_order])]))
    assert metrics["pf_n_matched"] == 2
    assert metrics["pf_loss_confusion"] == pytest.approx(0, abs=1e-20)
    assert metrics["pflow"] == pytest.approx(0, abs=1e-20)


# Extra collinear fragments must not take matches from particles with the correct energy
@pytest.mark.parametrize("truth_order,reco_order", [((0, 1), (0, 1, 2)), ((1, 0), (2, 1, 0)),
                                                   ((0, 1), (2, 0, 1)), ((1, 0), (1, 2, 0))])
@pytest.mark.parametrize("gain", [0.6, 1.0, 1.4])
def test_collinear_fragment(spec, truth_order, reco_order, gain):
    truth = [_truth_particle(22, 0, 10, 0, 0), _truth_particle(211, 1, 20, 0, 0, mass=0.139)]
    reco = [_reco_particle(130, 5, 0, 0),
            *[_reco_particle(p["pdg"], p["px"], p["py"], p["pz"], _particle_energy(p)) for p in truth]]
    reco = [{key: value if key == "pdg" else gain * value for key, value in p.items()} for p in reco]
    arrays = _pf_arrays([([truth[i] for i in truth_order], [reco[i] for i in reco_order])])
    particles = [pandora._pf_particles(spec, arrays, 0, side) for side in ("truth", "reco")]
    pairs = pandora._pf_angular_matched_pairs(*particles, spec["matching_max_angle_rad"])
    assert {(truth_order[t], reco_order[r]) for t, r in pairs} == {(0, 1), (1, 2)}
    for p in particles:
        p["classes"] = p["classes"][::-1]
        p["class_index"] = p["class_index"][::-1]
    np.testing.assert_array_equal(pandora._pf_angular_matched_pairs(*particles, spec["matching_max_angle_rad"]), pairs)


# Energy agreement only resolves angular ties, never a distinguishable angle difference
def test_matching_angle_priority(spec):
    truth = dict(energy=np.array([10.0]), unit=np.array([[1.0, 0.0, 0.0]]))
    angle = spec["matching_max_angle_rad"] / 2
    reco = dict(energy=np.array([5.0, 10.0]), unit=np.array([[1.0, 0.0, 0.0], [math.cos(angle), math.sin(angle), 0.0]]))
    np.testing.assert_array_equal(pandora._pf_angular_matched_pairs(truth, reco, spec["matching_max_angle_rad"]), [[0, 0]])


# Read parton acceptance from each reconstructed event even after file or event permutation
def test_reco_root_event_association(spec, tmp_path):
    arrays = _pf_arrays([_photons(), _photons(0.1), _photons(0.2), _photons(0.3)], [0, 0.8, None, -0.8])
    outputs = []
    for indices in ([0, 1, 2, 3], [3, 2, 1, 0]):
        paths = []
        for index in indices:
            path = tmp_path / f"event_{index}.root"
            with uproot.recreate(path) as output:
                output["events"] = {key: value[index : index + 1] for key, value in arrays.items()}
            paths.append(path)
        arrays_read = pandora._reco_chunks(paths, {"pflow": spec})
        metrics, pf_detail = pandora._objective_pflow(spec, arrays_read)
        detail = {"pflow": pf_detail}
        assert metrics["pf_n_events_input"] == 4
        assert metrics["pf_n_events_missing_parton"] == 1
        assert metrics["pf_n_events_parton_accepted"] == 1
        assert metrics["pf_total_energy_response_mean"] == pytest.approx(1)
        outputs.append(metrics)
    assert outputs[0] == pytest.approx(outputs[1])
    assert "ntuple" not in str(detail)


# Select the highest energy hard light quark rather than the first MC particle
def test_highest_energy_parton(spec):
    truth, reco = _photons()
    extra = _truth_particle(1, -1 / 3, 100, 0, 0, status=23)
    metrics, _ = pandora._objective_pflow(spec, _pf_arrays([(truth + [extra], reco)], [0.9]))
    assert metrics["pf_n_events_parton_accepted"] == 1


# Reject malformed particle data and unstable truth rather than silently changing visible energy
@pytest.mark.parametrize("change", ["nan", "inf", "negative", "spacelike", "unknown", "unstable"])
def test_invalid_particles(spec, change):
    truth, reco = _photons()
    if change in ("nan", "inf", "negative"):
        reco[0]["energy"] = {"nan": np.nan, "inf": np.inf, "negative": -1}[change]
    elif change == "spacelike":
        reco[0]["energy"] *= 0.9
    elif change == "unknown":
        reco[0]["pdg"] = 999999
    else:
        truth[0]["pdg"] = 111
    with pytest.raises(ValueError):
        pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))


# Neutrinos outside the visible state must not alter the objective or response
def test_neutrinos(spec):
    truth, reco = _photons()
    truth.append(_truth_particle(12, 0, 100, 0, 0))
    metrics, _ = pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))
    assert metrics["pflow"] == pytest.approx(0, abs=1e-20)
    assert metrics["pf_n_truth_visible"] == 4


# Compute genuine dataset records through the event objective API
def _record(spec, events, weight):
    metrics, detail = pandora._objective_pflow(spec, _pf_arrays(events))
    return dict(
        datacard=dict(weight=weight), metrics=metrics, details=dict(pflow=detail), dataset=dict(weight=weight)
    )


# Pool event means and variances, preserving zero weights and dataset splitting invariance
def test_weighted_combination(spec):
    records = [_record(spec, [_photons(0.8)], 1), _record(spec, [_photons(1.2)], 3), _record(spec, [_photons(0.1)], 0)]
    metrics, details = pandora._combine_datasets(records)
    assert metrics["pf_total_energy_response_mean"] == pytest.approx(1.1, abs=1e-6)
    assert metrics["pf_total_energy_residual_variance"] == pytest.approx(0.03, abs=1e-6)
    assert metrics["pf_loss_calibration"] == pytest.approx((0.1 / spec["energy_bias_target"])**2)
    assert metrics["pf_calibration_passed"] == 0
    pooled, _ = pandora._objective_pflow(spec, _pf_arrays([_photons(0.8), *[_photons(1.2)] * 3]))
    for key in (
        "pflow",
        "pf_total_energy_response_mean90",
        "pf_total_energy_response_sigma90",
        "pf_total_energy_residual_variance",
    ):
        assert metrics[key] == pytest.approx(pooled[key])
    profile = details["pflow"]["object_performance"]["classes"]["photon"]
    assert np.sum(profile["truth_weight"]) == pytest.approx(16)
    equal, _ = pandora._combine_datasets([_record(spec, [_photons(0.8)], 1), _record(spec, [_photons(1.2)], 1)])
    assert equal["pf_total_energy_residual_variance"] == pytest.approx(0.04, abs=1e-6)
    assert equal["pf_calibration_passed"] == 1
    assert equal["pflow"] < min(record["metrics"]["pflow"] for record in records[:2])


# Reject invalid weights at initialization and exclude disabled datasets exactly
@pytest.mark.parametrize("weight", [-1, np.nan, np.inf])
def test_invalid_weight(weight):
    with pytest.raises(ValueError, match="weight"):
        pandora._dataset_weight(dict(weight=weight))
    assert pandora._dataset_weight(dict(weight=0)) == 0


# Weighted efficiency and class confusion must use the same event weights as the loss
def test_weighted_efficiency_and_confusion(spec):
    truth, reco = _photons()
    wrong = [{**p, "pdg": 211} for p in reco]
    _, details = pandora._combine_datasets([_record(spec, [(truth, reco)], 1), _record(spec, [(truth, wrong)], 3)])
    detail = details["pflow"]
    profile = detail["object_performance"]["classes"]["photon"]
    result = diagnostics.binomial_profile(
        profile["truth_energy"], profile["truth_matched"], [1, 20], weights=profile["truth_weight"]
    )
    assert result["denominator"][0] == pytest.approx(16)
    assert result["rate"][0] == pytest.approx(0.25)
    matrix = diagnostics.confusion_matrices(detail["class_confusion"], low=None, high=None)
    photon, charged = (diagnostics.PF_CLASS_NAMES.index(c) for c in ("photon", "charged_hadron"))
    assert matrix["row_rate"][photon, photon] == pytest.approx(0.25)
    assert matrix["row_rate"][photon, charged] == pytest.approx(0.75)
    assert matrix["column_rate"][photon, photon] == pytest.approx(1)
    np.testing.assert_allclose(matrix["row_rate"][photon].sum(), 1)


# Lost and fake particles close the efficiency rows and purity columns under energy migration
@pytest.mark.parametrize("low,high", [(None, None), (0.1, 5.0), (5.0, 10.0), (10.0, 20.0)])
def test_confusion_migration(spec, low, high):
    truth = [_truth_particle(22, 0, 2, 0, 0), _truth_particle(130, 0, 0, 6, 0)]
    reco = [_reco_particle(211, 12, 0, 0), _reco_particle(13, -3, 0, 0), _reco_particle(11, 0, -13, 0)]
    _, detail = pandora._objective_pflow(spec, _pf_arrays([(truth, reco)]))
    matrix = diagnostics.confusion_matrices(detail["class_confusion"], low=low, high=high)
    np.testing.assert_allclose(matrix["row_numerator"].sum(axis=1), matrix["truth_denominator"])
    np.testing.assert_allclose(matrix["column_numerator"].sum(axis=0), matrix["reco_denominator"])
    np.testing.assert_allclose(matrix["row_rate"].sum(axis=1), matrix["truth_denominator"] > 0)
    np.testing.assert_allclose(matrix["column_rate"].sum(axis=0), matrix["reco_denominator"] > 0)
    photon, electron, muon, charged, neutral = (diagnostics.PF_CLASS_NAMES.index(c) for c in (
        "photon", "electron", "muon", "charged_hadron", "neutral_hadron"))
    if low is None:
        assert matrix["row_rate"][photon, charged] == pytest.approx(1)
        assert matrix["row_rate"][neutral, -1] == pytest.approx(1)
        assert matrix["column_rate"][-1, muon] == pytest.approx(1)
        assert matrix["column_rate"][-1, electron] == pytest.approx(1)
    elif low is not None and low >= 10.0:
        assert not np.any(matrix["truth_denominator"])
        assert matrix["column_rate"][photon, charged] == pytest.approx(1)
        assert matrix["column_rate"][-1, electron] == pytest.approx(1)


# Keep energy migration, empty events and weighted class counts consistent across diagnostics
def test_composition_energy_bins(spec):
    truth = [_truth_particle(22, 0, 0.5, 0, 0), _truth_particle(22, 0, -1, 0, 0),
             _truth_particle(130, 0, 0, 0.6, 0, mass=0.8), _truth_particle(211, 1, 0, -2, 0)]
    reco = [_reco_particle(22, 0.5, 0, 0, energy=1), _reco_particle(211, -1, 0, 0, energy=2),
            _reco_particle(2112, 0, 0.6, 0, energy=0.75), _reco_particle(22, 0, -1.5, 0)]
    records = [_record(spec, [(truth, reco)], 2),
               _record(spec, [([_truth_particle(22, 0, 3, 0, 0)], [])], 3)]
    _, details = pandora._combine_datasets(records)
    detail = details["pflow"]
    inclusive = diagnostics.composition_profiles(detail["events"])
    edges = (0.1, 0.5, 1, 2, 5)
    profiles = []
    for low, high in zip(edges[:-1], edges[1:], strict=True):
        profile = diagnostics.composition_profiles(detail["events"], low=low, high=high)
        profiles.append(profile)
        matrix = diagnostics.confusion_matrices(detail["class_confusion"], low=low, high=high)
        for side, prefix in (("gen", "truth"), ("pfo", "reco")):
            for name in profile["class_names"]:
                counts = profile[f"{side}_multiplicity"][name]
                index = diagnostics.PF_CLASS_NAMES.index(name)
                assert np.dot(counts, detail["weights"]) == pytest.approx(matrix[f"{prefix}_denominator"][index])
            np.testing.assert_allclose(profile[f"{side}_total_multiplicity"],
                                       np.sum(list(profile[f"{side}_multiplicity"].values()), axis=0))
        np.testing.assert_allclose(profile["pfo_total_multiplicity"][1], 0)
        np.testing.assert_allclose([v[1] for v in profile["pfo_momentum_fraction"].values()], 0)
    for side in ("gen", "pfo"):
        for name in inclusive["class_names"]:
            counts = sum(p[f"{side}_multiplicity"][name] for p in profiles)
            np.testing.assert_allclose(counts, inclusive[f"{side}_multiplicity"][name])
    np.testing.assert_allclose(profiles[1]["gen_total_multiplicity"], [1, 0])
    np.testing.assert_allclose(profiles[1]["pfo_multiplicity"]["neutral_hadron"], [1, 0])
    np.testing.assert_allclose(profiles[2]["gen_multiplicity"]["neutral_hadron"], [1, 0])
    np.testing.assert_allclose(profiles[2]["pfo_multiplicity"]["photon"], [2, 0])
    np.testing.assert_allclose(profiles[2]["gen_momentum_fraction"]["neutral_hadron"], [0.6 / 1.6, 0])
    np.testing.assert_allclose(profiles[2]["pfo_momentum_fraction"]["photon"], [1, 0])


# Check angular wrapping, wrong classes, competing candidates and absent reconstruction
def test_matching_diagnostics(spec):
    spec["matching_max_angle_rad"] = 0.1
    angles = [math.pi - 0.01, -math.pi + 0.04, 0.0, math.pi / 2]
    truth = [_truth_particle(211 if i == 3 else 22, 1 if i == 3 else 0,
                             10 * math.cos(phi), 10 * math.sin(phi), 0) for i, phi in enumerate(angles)]
    reco = [_reco_particle(22, 10 * math.cos(phi), 10 * math.sin(phi), 0)
            for phi in (-math.pi + 0.005, math.pi / 2 + 0.04)]
    arrays = _pf_arrays([(truth, reco)])
    _, detail = pandora._objective_pflow(spec, arrays)
    classes = detail["object_performance"]["classes"]
    photon = classes["photon"]
    np.testing.assert_array_equal(photon["truth_candidates"], [1, 1, 0])
    np.testing.assert_array_equal(photon["truth_assigned"], [True, False, False])
    np.testing.assert_array_equal(photon["truth_matched"], [True, False, False])
    np.testing.assert_allclose(photon["truth_angle"], [0.015, np.nan, np.nan], atol=1e-12)
    np.testing.assert_allclose(photon["truth_dr"], [0.015, np.nan, np.nan], atol=1e-12)
    charged = classes["charged_hadron"]
    np.testing.assert_array_equal(charged["truth_assigned"], [True])
    np.testing.assert_array_equal(charged["truth_matched"], [False])
    np.testing.assert_allclose(charged["truth_angle"], [0.04], atol=1e-12)
    _, empty = pandora._objective_pflow(spec, _pf_arrays([(truth, [])]))
    for name in ("photon", "charged_hadron"):
        profile = empty["object_performance"]["classes"][name]
        assert not np.any(profile["truth_assigned"])
        assert not np.any(profile["truth_candidates"])
        assert np.all(np.isnan(profile["truth_angle"]))


# Check physical energy bias and central variance with weighted, wrong and missing matches
@pytest.mark.parametrize("pdg,charge,name", [(22, 0, "photon"), (11, -1, "electron"), (13, -1, "muon"),
                                           (211, 1, "charged_hadron"), (130, 0, "neutral_hadron")])
def test_particle_resolution(spec, pdg, charge, name):
    truth = [_truth_particle(pdg, charge, 10, 0, 0)]
    records = [_record(spec, [(truth, [_reco_particle(pdg, 8, 0, 0)])], 1),
               _record(spec, [(truth, [_reco_particle(22 if pdg != 22 else 211, 12, 0, 0)])], 3),
               _record(spec, [(truth, [])], 5),
               _record(spec, [(truth, [_reco_particle(pdg, 100, 0, 0)])], 0)]
    _, details = pandora._combine_datasets(records)
    profile = details["pflow"]["object_performance"]["classes"][name]
    result = diagnostics.resolution_profile(profile, np.array([1, 10, 20, 40]))
    assert np.isnan(result["mean90"][[0, 2]]).all()
    assert result["median"][1] == pytest.approx(0.2)
    assert result["q68"][1] == pytest.approx(0.2)
    assert result["mean90"][1] == pytest.approx(0.1)
    assert result["sigma90"][1] == pytest.approx(np.sqrt(0.03))
    correct = diagnostics.resolution_profile(profile, np.array([1, 10]), correct_only=True)
    assert correct["mean90"][0] == pytest.approx(-0.2)
    assert correct["sigma90"][0] == pytest.approx(0)
    _, empty = pandora._objective_pflow(spec, _pf_arrays([(truth, [])]))
    result = diagnostics.resolution_profile(empty["object_performance"]["classes"][name], np.array([1, 20]))
    assert np.isnan(result["mean90"]).all()


# Separate central spread from bias and suppress tails in the robust resolution estimates
def test_resolution_tails(spec):
    arrays = _pf_arrays([_photons(scale) for scale in [1.1] * 9 + [1.3] * 9 + [3.0] * 2])
    _, detail = pandora._objective_pflow(spec, arrays)
    profile = detail["object_performance"]["classes"]["photon"]
    result = diagnostics.resolution_profile(profile, np.array([1, 20]))
    assert result["median"][0] == pytest.approx(0.3)
    assert result["q68"][0] == pytest.approx(0.1)
    assert result["mean90"][0] == pytest.approx(0.2)
    assert result["sigma90"][0] == pytest.approx(0.1)


# Preserve event correlations and weight scaling in statistical errors of the resolution estimators
def test_resolution_errors(spec):
    settings = json.loads((ROOT / "python/src/core/tune/drivers/pandora/settings.json").read_text())["resolution"]
    _, detail = pandora._objective_pflow(spec, _pf_arrays([_photons(s) for s in np.linspace(0.8, 1.2, 24)]))
    profile = detail["object_performance"]["classes"]["photon"]
    bins = np.array([1, 20, 40])
    full = diagnostics.resolution_profile(profile, bins, resampling=settings)
    _, selected = np.unique(profile["truth_event"], return_index=True)
    single = {k: v[selected] for k, v in profile.items() if k.startswith("truth_")}
    scaled = {**single, "truth_weight": single["truth_weight"] * 8}
    for other in (single, scaled):
        result = diagnostics.resolution_profile(other, bins, resampling=settings)
        for key, error in full["errors"].items():
            assert error[0] > 0
            assert np.isnan(error[1])
            np.testing.assert_allclose(result["errors"][key], error, rtol=1e-6, atol=1e-8)
    _, detail = pandora._objective_pflow(spec, _pf_arrays([_photons(1.1)]))
    result = diagnostics.resolution_profile(detail["object_performance"]["classes"]["photon"], bins, resampling=settings)
    assert all(np.isnan(error).all() for error in result["errors"].values())


# Close every differential contribution to the weighted objective, including wrong, lost and fake particles
def test_loss_contributions(spec):
    truth, reco = _photons(1.15)
    reco[0]["pdg"] = 211
    reco = reco[:-1] + [_reco_particle(130, 2, 3, 0)]
    _, details = pandora._combine_datasets([_record(spec, [(truth, reco)], 3),
                                          _record(spec, [_photons(0.8), _photons(1.2)], 1),
                                          _record(spec, [_photons(2)], 0)])
    detail = details["pflow"]
    fine = diagnostics.loss_data(detail, detail["object_performance"]["energy_bins_gev"])
    coarse = diagnostics.loss_data(detail, np.array([0.1, 200]))
    for name, expected in detail["loss_components"].items():
        assert fine["energy"][name].sum() == pytest.approx(expected)
        assert fine["particle"][name].sum() == pytest.approx(expected)
        np.testing.assert_allclose(fine["particle"][name], coarse["particle"][name], atol=1e-12)
    for event in detail["events"]:
        for name in ("confusion", "response", "momentum", "multiplicity"):
            assert event["loss_profile"][name].sum() == pytest.approx(event[name + "_loss"])
    classes = fine["classes"]
    assert fine["particle"]["confusion"][classes.index("neutral_hadron")] > 0
    assert fine["particle"]["multiplicity"][classes.index("charged_hadron")] > 0


# Keep compensating energy errors and covariance terms signed while preserving the total loss
def test_loss_cancellations(spec):
    truth = [_truth_particle(22, 0, 10, 0, 0), _truth_particle(211, 1, 0, 10, 0)]
    reco = [_reco_particle(22, 8, 0, 0), _reco_particle(211, 0, 14, 0)]
    perfect = [_reco_particle(22, 10, 0, 0), _reco_particle(211, 0, 10, 0)]
    _, detail = pandora._objective_pflow(spec, _pf_arrays([(truth, reco), (truth, perfect)]))
    data = diagnostics.loss_data(detail, detail["object_performance"]["energy_bins_gev"])
    photon = data["classes"].index("photon")
    for name in ("calibration", "resolution"):
        assert data["particle"][name][photon] < 0
        assert data["particle"][name].sum() == pytest.approx(detail["loss_components"][name])
    np.testing.assert_allclose(data["energy"]["multiplicity"], 0, atol=1e-12)


# Apply configured binning throughout diagnostics while conserving event weights and loss
@pytest.mark.parametrize("subdivisions", [1, 3])
def test_plot_binning(spec, tmp_path, subdivisions):
    settings = diagnostics.plot_settings()
    settings["binning"].update(energy=[2, 8], pt=[1, 4], costheta=[-1, 0, 1], eta=[-2, 0, 2],
                               extend_factor=4, subdivisions=subdivisions, matching_bins=7)
    settings["histogram"].update(min_bins=3, max_bins=3, multiplier=2)
    settings["resolution"]["resamples"] = 16
    path = tmp_path / "settings.json"
    path.write_text(json.dumps(settings))
    settings = diagnostics.plot_settings(path)
    _, detail = pandora._objective_pflow(spec, _pf_arrays([_photons(0.8), _photons(1.2)]))
    detail.update(diagnostics.plot_inputs(detail["events"], detail["weights"], settings=settings))
    prepared = diagnostics.prepare_plots(dict(pflow=detail))
    np.testing.assert_allclose(detail["object_performance"]["energy_bins_gev"], [2, 8, 32])
    np.testing.assert_allclose(prepared["bins"]["energy"][::subdivisions], [2, 8, 32])
    np.testing.assert_allclose(prepared["bins"]["pt"][::subdivisions], [1, 4, 16])
    np.testing.assert_allclose(prepared["bins"]["eta"][::subdivisions], [-2, 0, 2])
    assert len(prepared["histograms"][""]["counts"]) == 6
    assert prepared["histograms"][""]["counts"].sum() == pytest.approx(detail["weights"].sum())
    photon = prepared["particles"]["photon"]["axes"]["energy"]
    assert photon["matching"]["angle"][0].shape == (2 * subdivisions, 7)
    assert photon["rates"]["truth"]["denominator"].sum() == pytest.approx(8)
    for name, values in prepared["loss_contributions"]["energy"].items():
        assert values.sum() == pytest.approx(detail["loss_components"][name])
    compact = pandora._compact_objective_details(dict(pflow=detail))
    assert "plot_settings" not in compact["pflow"]


# Reject binning that could lose events or prevent automatic range extension
@pytest.mark.parametrize("section,key,value", [
    ("binning", "energy", [0, 1]), ("binning", "eta", [0, 0]),
    ("binning", "pt", [1, float("inf")]), ("binning", "costheta", [[-1, 1]]),
    ("binning", "extend_factor", 1), ("binning", "extend_factor", float("nan")),
    ("binning", "subdivisions", 0), ("binning", "matching_bins", 2.5),
    ("histogram", "min_bins", 0), ("histogram", "min_bins", 10000),
    ("histogram", "multiplier", True), ("resolution", "resamples", 1), ("resolution", "seed", -1)])
def test_plot_settings(section, key, value, tmp_path):
    settings = diagnostics.plot_settings()
    settings[section][key] = value
    path = tmp_path / "settings.json"
    path.write_text(json.dumps(settings))
    with pytest.raises(ValueError):
        diagnostics.plot_settings(path)


# Plot collection must leave every scalar objective diagnostic unchanged
def test_plot_collection_and_rendering(spec, tmp_path):
    arrays = _pf_arrays([_photons(0.8), _photons(1.2)])
    metrics, detail = pandora._objective_pflow(spec, arrays)
    compact, small = pandora._objective_pflow({**spec, "_collect_plot_details": False}, arrays)
    assert compact == pytest.approx(metrics)
    assert detail["cost"] == pytest.approx(metrics["pflow"])
    assert "events" not in small and "object_performance" not in small
    prepared = diagnostics.prepare_plots(dict(pflow=detail))
    photon = prepared["particles"]["photon"]
    rates = photon["axes"]["energy"]
    total = sum(result["rate"] for result in rates["outcomes"])
    np.testing.assert_allclose(total[rates["rates"]["truth"]["nonzero"]], 1)
    assert prepared["histograms"][""]["counts"].sum() == pytest.approx(detail["weights"].sum())
    paths = plots.write_plots(prepared, tmp_path)
    assert all(Path(path).stat().st_size > 0 for path in paths.values())
    assert "_events" not in pandora._compact_objective_details(dict(pflow=detail))["pflow"]


# Reject invalid objective scales before any reconstruction
@pytest.mark.parametrize(
    "key,value",
    [
        ("matched_momentum_response_target", 0),
        ("matched_momentum_response_target", np.inf),
        ("gen_final_state_min_energy", -1),
        ("reco_final_state_min_energy", np.nan),
        ("matching_max_angle_rad", 4),
        ("multiplicity_target", 0),
        ("multiplicity_target", np.nan),
        ("energy_resolution_target", 0),
        ("energy_resolution_target", np.inf),
        ("energy_bias_tolerance", 0),
        ("energy_bias_tolerance", -0.01),
        ("energy_bias_tolerance", 1),
        ("energy_bias_tolerance", np.nan),
        ("energy_bias_tolerance", np.inf),
        ("energy_bias_target", 0),
        ("energy_bias_target", -0.01),
        ("energy_bias_target", np.nan),
        ("energy_bias_target", np.inf),
    ],
)
def test_invalid_objective(spec, key, value):
    spec[key] = value
    with pytest.raises(ValueError):
        pandora._validate_pf_spec(spec)


# Resolve distinct file ranges and diagnose invalid or overlapping input selections
def test_dataset_ranges(tunesetup):
    selection = copy.deepcopy(tunesetup.param_selection)
    dataset = next(iter(selection["generator"]["datasets"].values()))
    dataset["file_ranges"] = [[0, 3], [8, 11]]
    paths = common.dataset_input_roots(tunesetup.param_paths, dataset)
    assert len(paths) == 8
    assert Path(paths[4]).name.endswith("008.root")
    assert Path(paths[-1]).name.endswith("011.root")
    for ranges in ([[0, 1], [1, 2]], [[2, 1]], [[0, 1.5]], []):
        with pytest.raises(ValueError):
            common.dataset_input_roots(tunesetup.param_paths, {**dataset, "file_ranges": ranges})


# Keep initial parameters bounded and apply changes without modifying the external XML
def test_catalog_xml_and_steering(tunesetup, tmp_path):
    catalog = tunesetup.aux_param_space["pandora_catalog"]
    source = Path(tunesetup.datacards[0]["settings_template"])
    before = source.read_bytes()
    tree = ET.parse(source)
    defaults = pandora.PandoraDriver().get_initial_param(tunesetup.param_space, tunesetup.aux_param_space, str(ROOT))
    pandora._apply_config_to_xml(tree, catalog, defaults)
    assert pandora._apply_config_to_xml(tree, catalog, defaults) == 0
    assert set(defaults) == set(tunesetup.param_space)
    active = common.active_tune_entries(catalog, tuning_tables=tunesetup.param_tuning_tables)
    effective = pandora._validate_config(catalog, defaults)
    fixed = pandora._fixed_config(catalog)
    assert not defaults.keys() & fixed.keys()
    assert effective == {**defaults, **fixed}
    for entry in active:
        lower, upper = common.default_bounds_for_entry(entry, tuning_tables=tunesetup.param_tuning_tables)
        assert lower <= effective[entry["key"]] <= upper
        if entry["tune_dtype"] == "int":
            assert float(effective[entry["key"]]).is_integer()
    entry = next(e for e in active if e["kind"] == "xml_numeric" and e["tune_dtype"] == "float")
    value = sum(entry["tune_bounds"]) / 2
    assert pandora._apply_config_to_xml(tree, catalog, {**defaults, entry["key"]: value}) == 1
    assert source.read_bytes() == before
    target, overrides = tmp_path / "reco.py", tmp_path / "overrides.json"
    pandora._write_json(overrides, {})
    pandora._copy_and_patch_steering(
        source_path=Path(tunesetup.datacards[0]["reco_steering"]),
        target_path=target,
        override_path=overrides,
        data_dir=tmp_path,
    )
    compile(target.read_text(), str(target), "exec")
    saved = tmp_path / "tunesetup.json"
    save_tunesetup(tunesetup=tunesetup, path=saved, simdriver="PANDORA", tune_default="TUNE0")
    restored = load_tunesetup(cdir=ROOT, simdriver="PANDORA", name=str(saved))
    assert restored.param_space.keys() == tunesetup.param_space.keys()
    assert pandora._validate_config(restored.aux_param_space["pandora_catalog"], defaults) == effective
    assert pandora._card_config(defaults, catalog)["parameters"] == effective


# Compare XML values independently against the algorithm scopes and bounds in the JSON card
def _assert_json_xml(path, card, bound):
    root = ET.parse(path).getroot()
    for scope, table in card["param_tuning_tables"].items():
        if scope == "wrapper":
            continue
        query = ".//" + "/.//".join(f"algorithm[@type='{name}']" for name in scope.split("/"))
        algorithms = root.findall(query) or root.findall(scope)
        assert algorithms, scope
        if scope == "ClusteringParent/ConeClustering":
            algorithms = algorithms[:1]
        for algorithm in algorithms:
            for tag, spec in table.items():
                elements = algorithm.findall(tag)
                assert len(elements) <= 1, (scope, tag)
                value = float(elements[0].text) if elements else card["param_optional_xml_numeric_defaults"][scope.split("/")[-1]][tag]
                assert value == pytest.approx(spec["bounds"][bound]), (scope, tag)


# Close the JSON-to-XML and wrapper loop through the real readers and writers at both range limits
@pytest.mark.parametrize("bound", [0, 1])
def test_json_xml_wrapper_roundtrip(tmp_path, bound):
    source = campaign_source("tune-pandora-v0")
    card = pyjson5.loads(Path(source).read_text())
    study = load_tunesetup(cdir=ROOT, simdriver="PANDORA", name=source)
    catalog = study.aux_param_space["pandora_catalog"]
    active = common.active_tune_entries(catalog, tuning_tables=study.param_tuning_tables)
    config = {entry["key"]: entry["tune_bounds"][bound] for entry in active}
    wrapper = {name: [value + 0.01 for value in values] for name, values in card["param_wrapper_defaults"].items()}
    wrapper.update({name: [spec["bounds"][bound]] for name, spec in card["param_tuning_tables"]["wrapper"].items()})
    config.update({entry["key"]: wrapper[entry["name"]][entry["index"]] for entry in catalog["wrapper_numeric"]})
    originals = {Path(study.datacards[0][key]): Path(study.datacards[0][key]).read_bytes()
                 for key in ("settings_template", "reco_steering")}
    outputs = tmp_path / "outputs"
    outputs.mkdir()
    for path, content in originals.items():
        (outputs / path.name).write_bytes(content)
    xml = outputs / Path(study.datacards[0]["settings_template"]).name
    python = outputs / Path(study.datacards[0]["reco_steering"]).name
    summary_path = tmp_path / "summary.json"
    fixed = pandora._fixed_config(catalog)
    config = {key: value for key, value in config.items() if key not in fixed}
    effective = pandora._validate_config(catalog, config)
    summary_path.write_text(json.dumps(dict(config=config, card_config=pandora._card_config(config, catalog)), allow_nan=False))
    summary = icetune_push.load_summary(summary_path)
    driver = pandora.PandoraDriver()
    rows = driver.push_parameters(summary=summary, target_path=str(outputs), cdir=str(ROOT), options={}, confirm=bool)
    assert len(rows) == len(effective)
    exported = tmp_path / "baseline.json"
    driver.push_parameters(summary=summary, target_path=str(exported), cdir=str(ROOT), options={}, confirm=bool)
    assert icetune_push.load_summary(exported)["config"] == effective
    _assert_json_xml(xml, card, bound)
    parameters = next(node.value for node in ast.walk(ast.parse(python.read_text())) if isinstance(node, ast.Assign)
                      and any(isinstance(target, ast.Attribute) and isinstance(target.value, ast.Name)
                              and target.value.id == "pandora" and target.attr == "Parameters" for target in node.targets))
    values = {key.value: ast.literal_eval(value) for key, value in zip(parameters.keys, parameters.values, strict=True)
              if key.value in wrapper}
    assert values.keys() == wrapper.keys()
    for name, expected in wrapper.items():
        assert list(map(float, values[name])) == pytest.approx(expected), name

    # Check the trial XML writer against the same JSON values and repeat the public push
    trial = ET.parse(study.datacards[0]["settings_template"])
    pandora._apply_config_to_xml(trial, catalog, config)
    trial.write(tmp_path / "trial.xml")
    _assert_json_xml(tmp_path / "trial.xml", card, bound)
    before = {path: path.read_bytes() for path in (xml, python)}
    driver.push_parameters(summary=summary, target_path=str(outputs), cdir=str(ROOT), options={}, confirm=bool)
    assert all(path.read_bytes() == content for path, content in {**originals, **before}.items())


# Discover plugin defaults in absent, empty and partially configured XML blocks
@pytest.mark.parametrize("plugin", ["", "<LCElectronId/>",
                                   "<LCElectronId><MaxProfileStart>4</MaxProfileStart></LCElectronId>"])
def test_plugin_settings(tmp_path, plugin):
    source = tmp_path / "settings.xml"
    source.write_text(f"<pandora><ElectronPlugin>LCElectronId</ElectronPlugin>{plugin}</pandora>")
    defaults = {"LCElectronId": {"MaxProfileStart": 4.5, "MaxProfileDiscrepancy": 0.6}}
    tables = {"LCElectronId": {
        "MaxProfileStart": dict(bounds=[3, 6], sampling=["between_bounds"], tune_dtype="float"),
        "MaxProfileDiscrepancy": dict(bounds=[0.4, 0.8], sampling=["between_bounds"], tune_dtype="float")}}
    catalog = common.build_catalog(source, optional_xml_numeric_defaults=defaults,
                                   wrapper_defaults={}, tuning_tables=tables)
    config = {entry["key"]: entry["tune_bounds"][1]
              for entry in common.active_tune_entries(catalog, tuning_tables=tables)}
    tree = ET.parse(source)
    pandora._apply_config_to_xml(tree, catalog, config)
    rendered, _ = pandora._render_xml_push(source, catalog, config)
    for root in (tree.getroot(), ET.fromstring(rendered)):
        for entry in catalog["xml_numeric"]:
            node = root.find(f"LCElectronId/{entry['tag']}")
            value = float(node.text) if node is not None else entry["default"]
            assert value == pytest.approx(config[entry["key"]])
    source.write_text(rendered)
    assert pandora._render_xml_push(source, catalog, config)[0] == rendered


# Optional cuts belong to each definition, while repeated named instances remain references
@pytest.mark.parametrize("fixed", [False, True])
def test_clustering_instances(tmp_path, fixed):
    source = tmp_path / "settings.xml"
    source.write_text("""<pandora>
      <algorithm type="ClusteringParent"><algorithm type="ConeClustering" instance="primary"/></algorithm>
      <algorithm type="ClusteringParent"><algorithm type="ConeClustering" instance="muon">
        <TanConeAngleFine>0.1</TanConeAngleFine></algorithm></algorithm>
      <algorithm type="ClusteringParent"><algorithm type="ConeClustering" instance="primary"/></algorithm>
    </pandora>""")
    defaults = {"ConeClustering": {"TanConeAngleFine": 0.3}}
    tables = {"ClusteringParent/ConeClustering": {
        "TanConeAngleFine": dict(bounds=[0.25, 0.25] if fixed else [0.2, 0.4],
                                sampling=["fixed"] if fixed else ["between_bounds"], tune_dtype="float")}}
    catalog = common.build_catalog(source, optional_xml_numeric_defaults=defaults,
                                   wrapper_defaults={}, tuning_tables=tables)
    entries = common.active_tune_entries(catalog, tuning_tables=tables)
    assert len(entries) == 1
    config = {} if fixed else {entries[0]["key"]: 0.25}
    tree = ET.parse(source)
    pandora._apply_config_to_xml(tree, catalog, config)
    algorithms = tree.getroot().findall("./algorithm/algorithm")
    assert float(algorithms[0].findtext("TanConeAngleFine")) == pytest.approx(0.25)
    assert float(algorithms[1].findtext("TanConeAngleFine")) == pytest.approx(0.1)
    assert algorithms[2].find("TanConeAngleFine") is None


# Reject misspelled steering parameters and empty numeric domains
def test_bad_parameter_domains(tunesetup):
    catalog = tunesetup.aux_param_space["pandora_catalog"]
    with pytest.raises(ValueError, match="unmatched"):
        common.active_tune_entries(catalog, tuning_tables={"wrapper": {"misspelled": {}}})
    entry = catalog["wrapper_numeric"][0]
    for bounds, dtype, mode in [([1, 1], "float", "between_bounds"), ([0.1, 0.2], "int", "between_bounds"),
                                ([0, np.inf], "float", "between_bounds"), ([1, 2], "float", "fixed"),
                                ([0.5, 0.5], "int", "fixed"), ([np.nan, np.nan], "float", "fixed")]:
        table = {"wrapper": {entry["name"]: dict(bounds=bounds, tune_dtype=dtype, sampling=[mode])}}
        with pytest.raises(ValueError, match="domain"):
            common.default_bounds_for_entry(entry, tuning_tables=table)


# Reject invalid values before XML mutation, including fractional hit counts
@pytest.mark.parametrize("failure", ["low", "high", "nan", "inf", "fraction", "unknown"])
@pytest.mark.parametrize("fixed", [False, True])
def test_trial_bounds(tunesetup, failure, fixed):
    catalog = tunesetup.aux_param_space["pandora_catalog"]
    config = pandora.PandoraDriver().get_initial_param(tunesetup.param_space, tunesetup.aux_param_space, ROOT)
    active = common.active_tune_entries(catalog, tuning_tables=tunesetup.param_tuning_tables)
    entry = next(e for e in active
                 if e["tune_dtype"] == "int" and (e["sampling_mode"] == "fixed") == fixed)
    low, high = entry["tune_bounds"]
    config[entry["key"]] = {"low": low - 1, "high": high + 1, "nan": np.nan,
                          "inf": np.inf, "fraction": low + 0.5, "unknown": low}[failure]
    if failure == "unknown":
        config["misspelled"] = 1
    tree = ET.parse(tunesetup.datacards[0]["settings_template"])
    before = ET.tostring(tree.getroot())
    with pytest.raises(ValueError, match="Pandora"):
        pandora._apply_config_to_xml(tree, catalog, config)
    assert ET.tostring(tree.getroot()) == before


# Detect impossible independent domains and retain ordering in saved configurations
@pytest.mark.parametrize("failure", ["overlap", "missing", "initial"])
def test_ordered_domains(tunesetup, failure):
    card = pyjson5.loads(Path(campaign_source("tune-pandora-v0")).read_text())
    left, right = card["param_ordered"][0][:2]
    left_scope, _, left_tag = left.rpartition("/")
    right_scope, _, right_tag = right.rpartition("/")
    spec = card["param_tuning_tables"][left_scope][left_tag]
    if failure == "overlap":
        spec["bounds"][1] = card["param_tuning_tables"][right_scope][right_tag]["bounds"][0] + 1
    elif failure == "initial":
        spec["initial"] = spec["bounds"][1] + 1
    else:
        card["param_ordered"][0][0] = "missing/parameter"
    with pytest.raises(ValueError, match="Pandora"):
        common.build_catalog(tunesetup.datacards[0]["settings_template"],
                             optional_xml_numeric_defaults=card["param_optional_xml_numeric_defaults"],
                             wrapper_defaults=card["param_wrapper_defaults"],
                             tuning_tables=card["param_tuning_tables"], ordered_parameters=card["param_ordered"])


# Ordered metadata must also protect manual configurations independently of domain construction
def test_manual_ordering(tunesetup):
    catalog = copy.deepcopy(tunesetup.aux_param_space["pandora_catalog"])
    config = pandora.PandoraDriver().get_initial_param(tunesetup.param_space, tunesetup.aux_param_space, ROOT)
    config = pandora._validate_config(catalog, config)
    left, right = catalog["ordered_parameters"][0][:2]
    entry = pandora._entry_by_key(catalog)[left]
    entry["tune_bounds"] = [entry["tune_bounds"][0], config[right] + 2]
    config[left] = config[right] + 1
    with pytest.raises(ValueError, match="ordering"):
        pandora._validate_config(catalog, config)
    restored = pandora._card_config(config, catalog)["catalog"]
    with pytest.raises(ValueError, match="ordering"):
        pandora._validate_config(restored, config)


# Stage files concurrently and detect mutations without copying unchecked inputs into trials
def test_atomic_input_staging(tmp_path):
    source, target = tmp_path / "input.root", tmp_path / "worker/input.root"
    with uproot.recreate(source) as output:
        output["events"] = _pf_arrays([_photons()])
    identity = runtime.file_identity(source)
    with ThreadPoolExecutor(4) as pool:
        results = list(pool.map(lambda _: runtime.stage_file(source, target, identity), range(4)))
    assert all(path == target for path in results)
    assert runtime.file_identity(target) == identity
    source.write_bytes(b"changed source")
    assert runtime.stage_file(source, target, identity) == target
    with pytest.raises(ValueError, match="checksum"):
        runtime.stage_file(source, tmp_path / "other.root", identity)
    target.chmod(0o644)
    target.write_bytes(b"corrupt cache")
    with pytest.raises(ValueError, match="[Cc]orrupt"):
        runtime.stage_file(source, target, identity)


# Verify extraction against a real XML content manifest and reject damaged cached files
def test_runtime_archive(tunesetup, tmp_path):
    source = Path(tunesetup.datacards[0]["settings_template"])
    archive, setup = tmp_path / "runtime.tar", tmp_path / "setup.sh"
    setup.write_text(pandora_setup.setup_text("/cvmfs/key4hep/setup.sh", "2026-02-26"))
    with tarfile.open(archive, "w") as output:
        output.add(source, arcname="run/settings.xml")
    manifest = dict(
        **runtime.file_identity(archive),
        archive=archive.name,
        files={"run/settings.xml": runtime.file_identity(source)},
        dependencies={str(setup): runtime.file_identity(setup)},
    )
    root = runtime.stage_runtime(manifest, tmp_path, tmp_path / "worker")
    assert (root / "run/settings.xml").read_bytes() == source.read_bytes()
    assert runtime.stage_runtime(manifest, tmp_path, tmp_path / "worker") == root
    (root / "run/settings.xml").chmod(0o644)
    (root / "run/settings.xml").write_text("changed")
    with pytest.raises(ValueError, match="[Cc]orrupt"):
        runtime.stage_runtime(manifest, tmp_path, tmp_path / "worker")


# Accept selected upstream packages and reject stale dates, hashes and absent build metadata
@pytest.mark.parametrize("change", ["none", "date", "hash", "stack", "empty", "cache", "install"])
def test_build_release_validation(tmp_path, change):
    base = "/cvmfs/sw-nightlies.hsf.org/key4hep/releases"
    platform = "x86_64-almalinux9-gcc14.2.0-opt"
    upstream = f"{base}/2026-02-25/{platform}/root/6.38.00-lovp3j"
    current = f"{base}/2026-02-26/{platform}/dd4hep/develop-fv5bd3"
    environment = dict(
        KEY4HEP_STACK=f"{base}/2026-02-26/{platform}/key4hep-stack/2026-02-26-stack/setup.sh",
        CMAKE_PREFIX_PATH=f"{upstream}:{current}",
    )
    for module in pandora_setup.MODULES:
        cache = tmp_path / module / "build/CMakeCache.txt"
        cache.parent.mkdir(parents=True)
        cache.write_text(f"ROOT_DIR:PATH={upstream}/cmake\nDD4hep_DIR:PATH={current}/cmake\n")
        (cache.parent.parent / "install").mkdir()
    if change in {"date", "hash"}:
        cache.write_text(cache.read_text().replace(
            current, current.replace("2026-02-26", "2026-02-25") if change == "date" else current + "-stale"
        ))
    elif change == "stack":
        environment["KEY4HEP_STACK"] = environment["KEY4HEP_STACK"].replace("2026-02-26", "2026-02-25")
    elif change == "empty":
        cache.write_text("# No dependency metadata\n")
    elif change == "cache":
        cache.rename(cache.with_suffix("._old"))
    elif change == "install":
        install = cache.parent.parent / "install"
        install.rename(install.with_suffix("._old"))
    if change == "none":
        pandora_setup.validate_build(tmp_path, environment, "2026-02-26")
    else:
        with pytest.raises(FileNotFoundError if change in {"cache", "install"} else ValueError):
            pandora_setup.validate_build(tmp_path, environment, "2026-02-26")


# Probe the actual pinned stack on hosts with Key4hep available
def test_pinned_key4hep_environment():
    setup = "/cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh"
    if not Path(setup).is_file():
        pytest.skip("Key4hep CVMFS is absent")
    environment = pandora_setup.key4hep_environment(setup, "2026-02-26")
    assert "/key4hep/releases/2026-02-26/" in environment["KEY4HEP_STACK"]
    assert pandora_setup.package_prefixes(environment["CMAKE_PREFIX_PATH"])
    from submit.lxplus.nodes import check_worker_runtime
    check_worker_runtime({setup: runtime.file_identity(setup)}, pandora_setup.worker_check(setup, "2026-02-26"))


# Fingerprints must distinguish objective settings and different inputs with equal basenames
def test_fingerprints(tunesetup):
    driver = pandora.PandoraDriver()
    cards = copy.deepcopy(tunesetup.datacards)
    first = driver.bootstrap_payload(runtime_sha256=None, datacards=cards)["fingerprint"]
    cards[0]["input_roots"] = ["/changed/location/" + Path(p).name for p in cards[0]["input_roots"]]
    assert first != driver.bootstrap_payload(runtime_sha256=None, datacards=cards)["fingerprint"]
    cards = copy.deepcopy(tunesetup.datacards)
    cards[0]["objectives"]["pflow"]["energy_resolution_target"] *= 2
    assert first != driver.bootstrap_payload(runtime_sha256=None, datacards=cards)["fingerprint"]


# Missing production inputs must fail preflight before saving or launching any trial
def test_preflight_missing_input(tunesetup, tmp_path):
    trial_tunesetup = SimpleNamespace(
        datacards=copy.deepcopy(tunesetup.datacards), aux_param_space=copy.deepcopy(tunesetup.aux_param_space)
    )
    trial_tunesetup.datacards[0]["input_roots"] = [str(tmp_path / "absent.root")]
    with pytest.raises(FileNotFoundError):
        pandora.PandoraDriver().prepare_tunesetup(tunesetup=trial_tunesetup, args=SimpleNamespace(cdir=str(tmp_path)))
    assert not list(tmp_path.rglob("*.root"))


# Preserve failed process diagnostics and enforce a process timeout through the real runner
def test_command_failure(tmp_path):
    with pytest.raises(RuntimeError, match="diagnostic line"):
        pandora._run_command(
            [sys.executable, "-c", "print('diagnostic line'); raise SystemExit(3)"],
            cwd=tmp_path,
            log_path=tmp_path / "run.log",
            max_t=10,
        )
    environment = runtime.command_environment()
    assert "LD_LIBRARY_PATH" not in environment
    assert environment["OMP_NUM_THREADS"] == "1"


# Retain a real input exception and its context before a long initialization log
def test_command_early_failure(tmp_path):
    command = """
import sys
import traceback
try:
    open('missing-input.root')
except OSError:
    traceback.print_exc(file=sys.stdout)
    print('Required reconstruction input could not be read')
    print('DDMarlinPandora INFO Init processor\\n' * 80, end='')
    print('EventLoopMgr ERROR Unable to initialize Algorithm: k4FWCore__Sequencer')
    print('ApplicationMgr ERROR Application Manager Terminated with error code 1')
    sys.exit(1)
"""
    log = tmp_path / "reco.log"
    with pytest.raises(RuntimeError) as caught:
        pandora._run_command([sys.executable, "-c", command], cwd=tmp_path, log_path=log, max_t=10)
    message = str(caught.value)
    assert "root_cause=FileNotFoundError:" in message
    assert "Traceback (most recent call last)" in message
    assert "Required reconstruction input could not be read" in message
    assert "Application Manager Terminated" in message
    assert "missing-input.root" in log.read_text()


# Reconstruction failures use finite trial penalties and retain the diagnostic reason
def test_failure_output(tunesetup, tmp_path):
    result = pandora.PandoraDriver().evaluate_trial_outputs(
        config={},
        param=dict(
            datacards=tunesetup.datacards,
            aux_param_space=tunesetup.aux_param_space,
            cdir=str(tmp_path),
            run_name="invalid",
            cost="pflow",
        ),
        trial_id="invalid",
        tunename="invalid",
    )
    assert result["error"] and result["metrics"]["pflow"] == pandora.PENALTY_COST
    assert result["likelihood"]["valid"] is False


# Keep retired ROOT outputs with the requested suffix
def test_reco_retirement(tmp_path):
    path = tmp_path / "reco_000.root"
    path.write_bytes(b"root")
    pandora._retire_reco_roots([path])
    assert path.with_name("reco_000.root._old").read_bytes() == b"root"


# Export parameters from the actual canonical fit summary representation
@pytest.mark.parametrize("canonical", [False, True])
def test_export(canonical, tmp_path):
    config = {"PWRAP_cal": 1.2}
    payload = (
        dict(best_fit=build_best_fit(values=config, objective_name="loss", objective_value=1, source="realized"))
        if canonical
        else dict(best_config=config)
    )
    source = tmp_path / "summary.json"
    source.write_text(json.dumps(payload))
    target = tmp_path / "baseline" / "pandora_best_config.json"
    command = [sys.executable, "-m", "core.icetune", "--cdir", str(tmp_path), "--push", source.name]
    result = subprocess.run(command + [str(target)], input="yes\n", text=True, capture_output=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(target.read_text()) == {"source": str(source), "config": config}
    assert icetune_push.load_summary(source, baseline=True)["config"] == config
    assert icetune_push.load_summary(target, baseline=True)["config"] == config

    # Apply the exported baseline through the same command to reconstruction steering
    steering = tmp_path / "run_reco_pandora.py"
    steering.write_text('pandora.Parameters = {"cal": ["1.0"]}\n')
    result = subprocess.run(command[:-1] + [str(target), str(steering)], input="yes\n",
                            text=True, capture_output=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert steering.read_text() == 'pandora.Parameters = {"cal": ["1.2"]}\n'


# Exercise Ray's actual packager so frozen runtime inputs reach workers while other scratch stays excluded
def test_ray_runtime_package(tmp_path):
    import zipfile

    from core.tune.runtime.ray import runtime_excludes
    from ray._private.runtime_env.packaging import create_package

    root = tmp_path / "checkout"
    included = ["tmp/icetune/pandora/first.tar", "tmp/icetune/pandora/second.tar", "python/input.dat"]
    files = [*included, "tmp/unrelated.txt", "tmp/icetune/pandora/cache/event.root", "python/module.py"]
    for relative in files:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(relative)
    target = tmp_path / "runtime.zip"
    create_package(str(root), target, include_gitignore=False, excludes=runtime_excludes(included))
    with zipfile.ZipFile(target) as archive:
        names = {name for name in archive.namelist() if not name.endswith("/")}
    assert names == {*included, "python/module.py"}


# Freeze actual ROOT and XML content deterministically without relying on an external build in this I/O test
def test_freeze_runtime_contents(tunesetup, tmp_path):
    import shutil

    root = tmp_path / "source"
    root.mkdir()
    shutil.copyfile(tunesetup.datacards[0]["settings_template"], root / "settings.xml")
    with uproot.recreate(root / "events.root") as output:
        output["events"] = _pf_arrays([_photons()])
    setup = root / "key4hep.sh"
    setup.write_text("# Key4hep setup input for archive I/O validation\n")
    description = {"setup": setup.read_text(), "dependencies": {str(setup): runtime.file_identity(setup)}}
    args = (root, ["settings.xml", "events.root"], tmp_path / "frozen", description)
    first, archive = runtime.freeze_runtime(*args)
    second, same = runtime.freeze_runtime(*args)
    assert first == second and archive == same
    first["archive"] = str(archive)
    staged = runtime.stage_runtime(first, tmp_path, tmp_path / "worker")
    assert (staged / "settings.xml").read_bytes() == (root / "settings.xml").read_bytes()
    (root / "settings.xml").write_text("changed XML")
    changed, other = runtime.freeze_runtime(*args)
    assert changed["sha256"] != first["sha256"] and other != archive
    setup.write_text("changed external setup")
    with pytest.raises(ValueError, match="dependency changed"):
        runtime.stage_runtime(first, tmp_path, tmp_path / "worker")


# Ignore disabled datasets before checking input existence
def test_disabled_dataset_preflight(tunesetup, tmp_path):
    code = {
        name: runtime.file_identity(Path(runtime.__file__).parent / name)
        for name in ("driver.py", "runtime.py", "tunesetup/common.py")
    }
    cards = copy.deepcopy(tunesetup.datacards) * 2
    cards = copy.deepcopy(cards)
    archive = tmp_path / "runtime.tar"
    with tarfile.open(archive, "w"):
        pass
    manifest = dict(driver_files=code, archive=str(archive), **runtime.file_identity(archive))
    cards[0] = {**cards[0], "runtime": manifest}
    cards[1] = {**cards[1], "weight": 0, "input_roots": [str(tmp_path / "absent.root")]}
    prepared = SimpleNamespace(datacards=cards)
    pandora.PandoraDriver().prepare_tunesetup(tunesetup=prepared, args=SimpleNamespace(cdir=str(tmp_path)))
    assert len(prepared.datacards) == 1


# Reject conflicting archive bytes under a compact runtime name without replacing them
def test_runtime_archive_conflict(tmp_path):
    source = tmp_path / "settings.xml"
    source.write_text("<pandora/>\n")
    args = (tmp_path, [source.name], tmp_path / "frozen", dict(setup="", dependencies={}))
    manifest, archive = runtime.freeze_runtime(*args)
    assert archive.stem == manifest["sha256"][:16]
    staged = runtime.stage_runtime(dict(manifest, archive=str(archive)), tmp_path, tmp_path / "worker")
    assert staged.name == archive.stem
    assert (staged / source.name).read_bytes() == source.read_bytes()
    archive.rename(archive.with_suffix(".tar._old"))
    archive.write_bytes(b"conflicting archive")
    with pytest.raises(ValueError, match="Conflicting Pandora runtime archive"):
        runtime.freeze_runtime(*args)
    assert archive.read_bytes() == b"conflicting archive"


# Exact energy and angular ties must not depend on collection ordering or species ordering
@pytest.mark.parametrize("mass", [0.0, 0.139])
@pytest.mark.parametrize("truth_order,reco_order", [(1, 1), (1, -1), (-1, 1), (-1, -1)])
def test_equal_energy_permutations(spec, mass, truth_order, reco_order):
    truth = [_truth_particle(22, 0, 10, 0, 0), _truth_particle(211, 1, math.sqrt(100 - mass**2), 0, 0, mass=mass)]
    reco = [_reco_particle(p["pdg"], p["px"], p["py"], p["pz"], _particle_energy(p)) for p in truth]
    metrics, _ = pandora._objective_pflow(spec, _pf_arrays([(truth[::truth_order], reco[::reco_order])]))
    assert metrics["pf_n_matched"] == 2
    assert metrics["pflow"] == pytest.approx(0, abs=1e-20)


# Compare the real assignment against exhaustive partial matchings of small particle sets
@pytest.mark.parametrize("nt,nr", [(1, 3), (3, 1), (2, 3), (3, 3)])
def test_assignment_enumeration(nt, nr):
    from itertools import product

    rng = np.random.default_rng(81)
    for _ in range(20):
        directions = [rng.normal(size=(size, 3)) for size in (nt, nr)]
        directions = [p / np.linalg.norm(p, axis=1)[:, None] for p in directions]
        truth, reco = [dict(unit=p, energy=rng.uniform(1, 10, len(p))) for p in directions]
        angles = np.arccos(np.clip(directions[0] @ directions[1].T, -1, 1))
        maximum = 1.5
        candidates = []
        for assigned in product(range(-1, nr), repeat=nt):
            pairs = [(t, r) for t, r in enumerate(assigned) if r >= 0]
            if len({r for _, r in pairs}) == len(pairs) and all(angles[t, r] <= maximum for t, r in pairs):
                candidates.append((-len(pairs), sum(angles[t, r] for t, r in pairs)))
        expected = min(candidates)
        pairs = pandora._pf_angular_matched_pairs(truth, reco, maximum)
        assert -len(pairs) == expected[0]
        assert sum(angles[t, r] for t, r in pairs) == pytest.approx(expected[1], abs=1e-14)


# Stream real ROOT data without changing event selection or combined objective statistics
@pytest.mark.parametrize("step_size", [1, 2, 10])
def test_reco_chunk_invariance(spec, tmp_path, step_size):
    arrays = _pf_arrays([_photons(s) for s in (1, 0.8, 1.2, 0.3, 1.1)], [0, 0, 0, 0.9, None])
    path = tmp_path / "events.root"
    with uproot.recreate(path) as output:
        output[spec["tree"]] = arrays
    metrics, detail = pandora._objective_pflow(spec, arrays)
    streamed, streamed_detail = pandora._objective_pflow(spec, pandora._reco_chunks(path, {"pflow": spec}, step_size=step_size))
    assert streamed == pytest.approx(metrics)
    np.testing.assert_array_equal(streamed_detail["values"], detail["values"])
    assert streamed_detail["parton_acceptance"] == detail["parton_acceptance"]


# Reject stale XML positions before either the tree or the source file can be modified
@pytest.mark.parametrize("changed", ["<Beta>2</Beta><Alpha>1</Alpha>",
    '<algorithm type="other"><Alpha>1</Alpha></algorithm>'])
def test_xml_target_identity(tmp_path, changed):
    path = tmp_path / "settings.xml"
    original = ('<algorithm type="original"><Alpha>1</Alpha></algorithm>' if "algorithm" in changed
                else "<Alpha>1</Alpha><Beta>2</Beta>")
    path.write_text(f"<pandora><Keep>3</Keep>{original}</pandora>")
    catalog = pandora.build_catalog(path, optional_xml_numeric_defaults={}, wrapper_defaults={})
    config = {entry["key"]: 7 for entry in catalog["xml_numeric"]}
    path.write_text(f"<pandora><Keep>3</Keep>{changed}</pandora>")
    tree, before = ET.parse(path), path.read_bytes()
    tree_before = ET.tostring(tree.getroot())
    with pytest.raises(ValueError, match="structure changed"):
        pandora._apply_config_to_xml(tree, catalog, config)
    assert ET.tostring(tree.getroot()) == tree_before
    with pytest.raises(ValueError, match="structure changed"):
        pandora._render_xml_push(path, catalog, config)
    assert path.read_bytes() == before


# Preserve the complete baseline and trial configuration through XML and exported catalog readers
@pytest.mark.parametrize("different", [False, True])
def test_applied_baseline(tmp_path, different):
    path = tmp_path / "settings.xml"
    path.write_text("<pandora><Alpha>1</Alpha><Beta>2</Beta></pandora>")
    catalog = pandora.build_catalog(path, optional_xml_numeric_defaults={}, wrapper_defaults={})
    alpha, beta = [entry["key"] for entry in catalog["xml_numeric"]]
    cards = [dict(weight=1, baseline={alpha: 4, beta: 5}), dict(weight=2, baseline={alpha: 4, beta: 6 if different else 5})]
    if different:
        with pytest.raises(ValueError, match="same applied"):
            pandora._trial_config(catalog, {alpha: 7}, cards)
        return
    effective = pandora._trial_config(catalog, {alpha: 7}, cards)
    assert effective == {alpha: 7, beta: 5}
    exported = json.loads(json.dumps(pandora._card_config(effective, catalog)))
    rendered, _ = pandora._render_xml_push(path, exported["catalog"], exported["parameters"])
    tree = ET.parse(path)
    pandora._apply_config_to_xml(tree, catalog, effective)
    assert [float(e.text) for e in tree.getroot()] == [7, 5]
    assert ET.tostring(ET.fromstring(rendered)) == ET.tostring(tree.getroot())


# Signed scalar literals remain safe to update without evaluating Python expressions
@pytest.mark.parametrize("literal", ["-0.5", "+0.5", "'-0.5'", '"-0.5"', "0.5"])
def test_signed_wrapper_push(tmp_path, literal):
    path = tmp_path / "steering.py"
    path.write_text(f"parameters = {{'Cut': [{literal}]}}\n")
    catalog = {"wrapper_numeric": [dict(key="PWRAP_Cut", kind="wrapper_numeric", name="Cut", index=0)]}
    rendered, old = pandora._render_python_push(path, catalog, {"PWRAP_Cut": -0.75})
    values = ast.literal_eval(ast.parse(rendered).body[0].value)
    assert float(values["Cut"][0]) == pytest.approx(-0.75)
    assert old["PWRAP_Cut"] == pytest.approx(float(ast.literal_eval(literal)))
    assert isinstance(values["Cut"][0], str) == isinstance(ast.literal_eval(literal), str)


# Refuse expressions and nonfinite or boolean wrapper values before any write
@pytest.mark.parametrize("literal", ["1 + 2", "float('nan')", "True", "'nan'", "1e999"])
def test_invalid_wrapper_literal(tmp_path, literal):
    path = tmp_path / "steering.py"
    before = f"parameters = {{'Cut': [{literal}]}}\n"
    path.write_text(before)
    catalog = {"wrapper_numeric": [dict(key="PWRAP_Cut", kind="wrapper_numeric", name="Cut", index=0)]}
    with pytest.raises(ValueError, match="not numeric"):
        pandora._render_python_push(path, catalog, {"PWRAP_Cut": 0.5})
    assert path.read_text() == before


# Repeated particles in the same event cannot increase the number of independent measurements
@pytest.mark.parametrize("copies", [1, 10])
@pytest.mark.parametrize("scale", [1.0, 3.0])
def test_event_rate_errors(spec, copies, scale):
    truth, reco = _photons()
    wrong = [{**p, "pdg": 211} for p in reco]
    _, detail = pandora._objective_pflow(spec, _pf_arrays([(truth * copies, reco * copies), (truth * copies, wrong * copies)]))
    profile = detail["object_performance"]["classes"]["photon"]
    rate = diagnostics.binomial_profile(profile["truth_energy"], profile["truth_matched"], [1, 20],
        weights=scale * profile["truth_weight"], events=profile["truth_event"])
    assert rate["rate"][0] == pytest.approx(0.5)
    np.testing.assert_allclose(rate["yerr"][:, 0], 1 / math.sqrt(12))
    matrix = diagnostics.confusion_matrices(detail["class_confusion"], low=None, high=None)
    photon = diagnostics.PF_CLASS_NAMES.index("photon")
    assert matrix["row_lower"][photon, photon] == pytest.approx(0.5 - 1 / math.sqrt(12))
    assert matrix["row_upper"][photon, photon] == pytest.approx(0.5 + 1 / math.sqrt(12))


# Four independent successful events retain the asymmetric Wilson interval at the physical boundary
def test_confusion_bounds(spec):
    _, detail = pandora._objective_pflow(spec, _pf_arrays([_photons()] * 4))
    matrix = diagnostics.confusion_matrices(detail["class_confusion"], low=None, high=None)
    photon = diagnostics.PF_CLASS_NAMES.index("photon")
    for prefix in ("row", "column"):
        assert matrix[prefix + "_rate"][photon, photon] == pytest.approx(1)
        assert matrix[prefix + "_lower"][photon, photon] == pytest.approx(0.8)
        assert matrix[prefix + "_upper"][photon, photon] == pytest.approx(1)


# Preserve binary64 parameter values even for changes smaller than a fixed decimal tolerance
@pytest.mark.parametrize("value", [np.nextafter(1.0, 2.0), 1.2345678901234567, 1e-15])
def test_parameter_precision(tmp_path, value):
    xml, python = tmp_path / "settings.xml", tmp_path / "steering.py"
    xml.write_text("<pandora><Cut>1</Cut></pandora>")
    python.write_text("parameters = {'Gain': ['1']}\n")
    catalog = pandora.build_catalog(xml, optional_xml_numeric_defaults={}, wrapper_defaults={"Gain": [1]})
    config = {entry["key"]: value for entry in pandora._entry_by_key(catalog).values()}
    tree = ET.parse(xml)
    assert pandora._apply_config_to_xml(tree, catalog, config) == 1
    rendered, _ = pandora._render_xml_push(xml, catalog, config)
    pushed, _ = pandora._render_python_push(python, catalog, config)
    applied = [float(tree.getroot()[0].text), float(ET.fromstring(rendered)[0].text),
               float(ast.literal_eval(ast.parse(pushed).body[0].value)["Gain"][0]),
               float(pandora._wrapper_overrides_from_config(catalog, config)["Gain"][0])]
    assert all(v.hex() == float(value).hex() for v in applied)


# XML boolean values are discrete and must never be silently thresholded
@pytest.mark.parametrize("value", [-1, 0.3, 2])
def test_boolean_domain(tmp_path, value):
    source = tmp_path / "settings.xml"
    source.write_text("<pandora><Enabled>true</Enabled></pandora>")
    catalog = pandora.build_catalog(source, optional_xml_numeric_defaults={}, wrapper_defaults={})
    key = catalog["xml_boolean"][0]["key"]
    with pytest.raises(ValueError, match="boolean"):
        pandora._validate_config(catalog, {key: value})


# Exported XML baselines retain target identities, including parameters from an earlier push
@pytest.mark.parametrize("reordered", [False, True])
def test_baseline_xml_identity(tmp_path, reordered):
    source, baseline = tmp_path / "settings.xml", tmp_path / "baseline.json"
    source.write_text("<pandora><Alpha>1</Alpha><Beta>2</Beta></pandora>")
    catalog = pandora.build_catalog(source, optional_xml_numeric_defaults={}, wrapper_defaults={})
    driver = pandora.PandoraDriver()
    for entry in catalog["xml_numeric"]:
        config = {entry["key"]: 7}
        driver.push_parameters(summary=dict(config=config, card_config=pandora._card_config(config, catalog)),
            target_path=str(baseline), cdir=str(tmp_path), options={}, confirm=bool)
    summary = icetune_push.load_summary(baseline, baseline=True)
    assert len(summary["config"]) == len(summary["card_config"]["catalog"]["xml_numeric"]) == 2
    if reordered:
        source.write_text("<pandora><Beta>2</Beta><Alpha>1</Alpha></pandora>")
        with pytest.raises(ValueError, match="structure changed"):
            driver.push_parameters(summary=summary, target_path=str(source), cdir=str(tmp_path), options={}, confirm=bool)
    else:
        driver.push_parameters(summary=summary, target_path=str(source), cdir=str(tmp_path), options={}, confirm=bool)
        assert [float(e.text) for e in ET.parse(source).getroot()] == [7, 7]
    with pytest.raises(ValueError, match="recorded parameter catalog"):
        driver.push_parameters(summary=dict(config=summary["config"]), target_path=str(source),
                               cdir=str(tmp_path), options={}, confirm=bool)


# Publish real Pandora diagnostics through Ray and lxplus with no reconstruction CPU available
def test_ray_plot_publication(spec, tunesetup, tmp_path):
    import time

    import ray
    from core.tune import core as icetune
    from core.tune.backends import ray as backend

    from submit.lxplus import outputs as publication

    metrics, detail = pandora._objective_pflow(spec, _pf_arrays([_photons(0.9), _photons(1.1)]))
    param = dict(cdir=str(tmp_path), run_name="plots", cost="pflow", plot=True, pickle_dump=False,
                 max_t=120, optimization=dict(optimizer="hebo"))
    outputs = icetune._normalize_trial_outputs(outputs=dict(metrics=metrics, results=dict(
        plot_details={"pflow": detail}, objective_summary=dict(metrics=metrics, plots={}))),
        config=pandora.PandoraDriver().get_initial_param(tunesetup.param_space, tunesetup.aux_param_space, str(ROOT)),
        param=param, trial_id="accepted", tunename="accepted")
    trial_dir = tmp_path / "trial"
    descriptor = icetune.maybe_dump_trial_payload(outputs=outputs, param={**param, "pickle_dump": True},
        destination=trial_dir / "results" / icetune.trial_pickle_filename(outputs))
    descriptor["relative_path"] = "results/" + descriptor["filename"]
    transfer = backend.collect_ray_trial_outputs(trial_dir=trial_dir, run_name=param["run_name"],
        descriptor=descriptor, rendered_initial=False, rendered_best=False)
    ray.shutdown()
    ray.init(num_cpus=0, resources={"icetune_head": 1}, include_dashboard=False)
    callback = None
    try:
        experiment = tmp_path / "run"
        state = backend.create_global_state(experiment_dir=str(experiment), cdir=str(tmp_path),
            run_name=param["run_name"], cost=param["cost"], render_figures=True, require_head_resource=True, keep_pickles=False)
        token = ray.get(state.publish_trial_outputs.remote(trial_id="accepted", cost=metrics["pflow"],
            transfer=transfer, descriptor=descriptor, rendered_initial=False, rendered_best=False, render_pending=True, initial=True))
        callback = backend.TrialOutputCallback(experiment_dir=str(experiment), param=param, global_state=state,
                                               simdriver=pandora.PandoraDriver())
        callback.on_trial_complete(iteration=0, trials=[], trial=SimpleNamespace(last_result=dict(
            trial_id="accepted", ray_output_token=token["output_token"])))
        deadline = time.monotonic() + 120
        while callback.rendered_trial_id != "accepted":
            callback.flush()
            assert not callback.output_jobs["render"]["failures"]
            assert time.monotonic() < deadline, "Pandora diagnostics waited for reconstruction CPU resources"
            time.sleep(0.05)
        assert callback.initial_rendered
        figure_dir = callback.figure_dir
        assert not list((experiment / "results").glob("*.pkl"))
        packed, _, revision, deferred = publication.snapshot_live_figures(figure_dir=figure_dir,
            run_name=param["run_name"], figure_revision=None)
        assert packed is not None and not deferred
        shared = tmp_path / "shared"
        fingerprint = "a" * 64
        publication.restore_outputs(payload=packed, shared_root=shared, run_name=param["run_name"], campaign_fingerprint=fingerprint)
        published = publication.campaign_figure_path(shared, param["run_name"], fingerprint)
        for relative in ("pflow/components/loss_components.png", "pflow/resolution/photon_resolution_vs_energy.png"):
            content = (published / relative).read_bytes()
            assert content.startswith(b"\x89PNG")
            assert (published / "init" / relative).read_bytes() == content
        for root in (published, published / "init"):
            summary = json.loads((root / "summary.json").read_text())
            assert summary["trial_id"] == "accepted"
            assert summary["objective_summary"]["plots"]
            assert all((root / path).is_file() and not Path(path).is_absolute()
                       for path in summary["objective_summary"]["plots"].values())
        unchanged, _, _, _ = publication.snapshot_live_figures(figure_dir=figure_dir, run_name=param["run_name"],
            figure_revision=revision)
        assert unchanged is None
    finally:
        if callback is not None and callback.output_jobs["render"]["future"] is not None:
            ray.cancel(callback.output_jobs["render"]["future"], force=True)
        ray.shutdown()

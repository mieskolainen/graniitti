# Tests for the standalone ICEBO proposal optimizer
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import builtins
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch
from core import resource

ROOT = Path(__file__).resolve().parents[3]

from core.inference.gp import ExactGaussianProcess
from core.numerics.smooth import log_fatplus
from core.tune.optimizers.icebo.config import load_icebo_config, override_icebo_config
from core.tune.optimizers.icebo.optimizer import ICEBO
from core.tune.optimizers.icebo.output import StandardizeOutput
from core.tune.optimizers.icebo.qlognei import QLogNoisyExpectedImprovement, sobol_normal
from core.tune.optimizers.icebo.space import IceboSpace
from core.tune.optimizers.icebo.turbo import TurboState, turbo_bounds
from core.tune.parameters import tools
from core.tune.search import SearchState


# Compute settings derived from the authoritative campaign JSON
def _icebo_settings(overrides: dict[str, object] | None = None) -> dict:
    settings = load_icebo_config(resource("tune/settings/icebo.json"))
    return override_icebo_config(settings, overrides or {})


# Compute fast deterministic settings for optimizer tests
def _fast_optimizer(bounds: dict, **overrides) -> ICEBO:
    runtime = {"seed": 17, "warmup": 4}
    for key in ("parameter_topology", "seed", "warmup"):
        if key in overrides:
            runtime[key] = overrides.pop(key)
    paths = {
        "direction": "search.direction",
        "dtype": "search.dtype",
        "trust_region": "turbo.enabled",
        "gradient_steps": "acquisition.gradient_steps",
    }
    config_overrides = {
        "gp.fit_steps": 12,
        "gp.fit_restarts": 1,
        "acquisition.raw_samples": 128,
        "acquisition.mc_samples": 32,
        "acquisition.restarts": 4,
        "acquisition.gradient_steps": 6,
        "acquisition.optimization_batch_size": 32,
    }
    for key, value in overrides.items():
        config_overrides[paths[key]] = value
    settings = _icebo_settings(config_overrides)
    return ICEBO(bounds, settings=settings, **runtime)


# Verify a closed Bayesian-optimization loop improves a smooth objective
def test_closed_loop_finds_quadratic_minimum():
    bounds = {"x": {"type": "uniform", "lower": -1.0, "upper": 1.0}}
    optimizer = _fast_optimizer(bounds, seed=11)
    values = []

    for _ in range(10):
        config = optimizer.suggest()[0]
        value = (config["x"] - 0.23) ** 2
        values.append(value)
        optimizer.observe(config, value)

    assert min(values) < 5.0e-3
    assert optimizer.diagnostics()["fit"]["surrogate"] == "exact_matern52_map_gp"


# Verify one call returns distinct greedily conditioned parallel proposals
def test_parallel_proposals_distinct_reserved():
    bounds = {
        "x": {"type": "uniform", "lower": 0.0, "upper": 1.0},
        "y": {"type": "uniform", "lower": 0.0, "upper": 1.0},
    }
    optimizer = _fast_optimizer(bounds, warmup=3)
    initial = optimizer.suggest(3)
    optimizer.observe(initial, [(row["x"] - 0.4) ** 2 + (row["y"] - 0.6) ** 2 for row in initial])

    batch, diagnostics = optimizer.suggest(4, return_diagnostics=True)
    tuples = {(round(row["x"], 12), round(row["y"], 12)) for row in batch}

    assert len(tuples) == 4
    assert diagnostics[0]["acquisition"] == "qlognei"
    assert all(item["joint_batch_size"] == 1 for item in diagnostics)
    assert all(item["batch_index"] == index for index, item in enumerate(diagnostics))
    assert optimizer.diagnostics()["pending"] == 4


# Verify the default acquisition remains qLogNEI in higher dimensions
def test_default_acquisition_qlognei_six_dimensions():
    bounds = {f"x{index}": {"type": "uniform", "lower": 0.0, "upper": 1.0} for index in range(6)}
    optimizer = _fast_optimizer(bounds, warmup=0)
    rows = np.random.default_rng(8).uniform(size=(7, 6))
    configs = [{f"x{index}": value for index, value in enumerate(row)} for row in rows]
    optimizer.observe(configs, np.sum((rows - 0.4) ** 2, axis=1))

    _, diagnostics = optimizer.suggest(return_diagnostics=True)

    assert diagnostics[0]["acquisition"] == "qlognei"
    assert diagnostics[0]["mc_samples"] == 32


# Verify known errors and learned noise produce separated finite uncertainties
def test_heteroscedastic_posterior_finite():
    bounds = {"x": {"type": "uniform", "lower": -1.0, "upper": 1.0}}
    optimizer = _fast_optimizer(bounds, warmup=0)
    grid = np.linspace(-1.0, 1.0, 16)
    configs = [{"x": float(value)} for value in grid]
    values = np.sin(3.0 * grid)
    errors = 0.01 + 0.15 * (grid + 1.0) / 2.0
    optimizer.observe(configs, values, errors)

    prediction = optimizer.predict([{"x": -0.5}, {"x": 0.5}])

    assert np.all(np.isfinite(prediction["mean"]))
    assert np.all(prediction["epistemic_std"] > 0.0)
    assert np.all(prediction["aleatoric_std"] > 0.0)
    assert optimizer.diagnostics()["fit"]["known_heteroscedastic"] is True


# Verify exact known variances give precise replicas the correct likelihood weight
def test_replicate_variance_weights():
    space = IceboSpace(
        {"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}},
        parameter_topology=None,
        device=torch.device("cpu"),
        dtype=torch.float64,
    )
    inputs = space.configs_to_tensor([{"x": 0.5}, {"x": 0.5}])
    surrogate = ExactGaussianProcess(
        inputs,
        torch.tensor([0.0, 10.0], dtype=torch.float64),
        Y_err=torch.tensor([1.0e-2, 1.0e2], dtype=torch.float64),
        geometry=space.geometry,
        settings=_icebo_settings()["gp"],
        device=torch.device("cpu"),
        dtype=torch.float64,
    )

    posterior = surrogate.posterior(inputs[:1])

    assert abs(float(posterior.mean[0].detach().cpu())) < 1.0e-3


# Verify all observation counts retain the exact full GP
def test_exact_full_gp_for_larger_training():
    bounds = {
        "x": {"type": "uniform", "lower": -1.0, "upper": 1.0},
        "y": {"type": "uniform", "lower": -1.0, "upper": 1.0},
    }
    optimizer = _fast_optimizer(
        bounds,
        warmup=0,
    )
    rng = np.random.default_rng(4)
    rows = rng.uniform(-1.0, 1.0, size=(48, 2))
    configs = [{"x": row[0], "y": row[1]} for row in rows]
    values = np.sum(rows**2, axis=1)
    optimizer.observe(configs, values)

    prediction = optimizer.predict([{"x": 0.0, "y": 0.0}])

    assert np.isfinite(prediction["mean"][0])
    assert (
        optimizer.diagnostics()["fit"]["surrogate"] == "exact_matern52_map_gp"
    )


# Verify phase topology identifies the two sides of a periodic branch cut
def test_space_phase_topology():
    bounds = {"phi": {"type": "uniform", "lower": -math.pi, "upper": math.pi}}
    topology = {
        "schema_version": 1,
        "groups": [
            {
                "kind": "phase",
                "base": "phi",
                "parameters": ["phi"],
                "period": 2.0 * math.pi,
            }
        ],
    }
    space = IceboSpace(
        bounds,
        parameter_topology=topology,
        device=torch.device("cpu"),
        dtype=torch.float64,
    )
    epsilon = 1.0e-6
    physical = space.configs_to_tensor([{"phi": -math.pi + epsilon}, {"phi": math.pi - epsilon}])
    unit = space.to_unit(physical)
    distance = space.minimum_topology_distance(unit[:1], unit[1:])

    assert distance[0] < 3.0e-6


# Verify ICEBO uses the coupled resonance norm, phase and projective direction
def test_space_complex_projective_topology():
    bounds = {
        "norm": {"type": "uniform", "lower": 0.0, "upper": 3.0},
        "phi": {"type": "uniform", "lower": -math.pi, "upper": math.pi},
        "theta0": {"type": "uniform", "lower": -0.5 * math.pi, "upper": 0.5 * math.pi},
        "theta1": {"type": "uniform", "lower": -0.5 * math.pi, "upper": 0.5 * math.pi},
    }
    topology = {
        "schema_version": 1,
        "groups": [
            {
                "kind": "polar_projective",
                "base": "g_ls",
                "parameters": ["norm", "phi", "theta0", "theta1"],
                "period": 2.0 * math.pi,
            }
        ],
    }
    space = IceboSpace(bounds, parameter_topology=topology, device=torch.device("cpu"), dtype=torch.float64)
    angle = 0.37
    physical = space.configs_to_tensor(
        [
            {"norm": 2.0, "phi": 0.2, "theta0": angle, "theta1": 0.5 * math.pi},
            {"norm": 2.0, "phi": 0.2 - math.pi, "theta0": -angle, "theta1": -0.5 * math.pi},
        ]
    )
    unit = space.to_unit(physical)

    distance = space.minimum_topology_distance(unit[:1], unit[1:])

    assert bool(space.periodic_mask[space.names.index("phi")])
    assert space.direction_groups == [
        {"kind": "projective", "indices": [space.names.index("theta0"), space.names.index("theta1")]}
    ]
    assert distance[0] == pytest.approx(0.0, abs=2.0e-8)


# Verify ICEBO samples a radial projective local continuum coupling
def test_space_radial_projective_topology():
    bounds = {
        "norm": {"type": "uniform", "lower": 0.0, "upper": 6.0},
        "theta0": {"type": "uniform", "lower": -0.5 * math.pi, "upper": 0.5 * math.pi},
        "theta1": {"type": "uniform", "lower": -0.5 * math.pi, "upper": 0.5 * math.pi},
    }
    topology = {
        "schema_version": 1,
        "groups": [
            {
                "kind": "radial_projective",
                "base": "helicity",
                "parameters": ["norm", "theta0", "theta1"],
                "weights": [2.0, 2.0, 1.0],
            }
        ],
    }
    space = IceboSpace(bounds, parameter_topology=topology, device=torch.device("cpu"), dtype=torch.float64)
    angle = 0.37
    physical = space.configs_to_tensor(
        [
            {"norm": 3.0, "theta0": angle, "theta1": 0.5 * math.pi},
            {"norm": 3.0, "theta0": -angle, "theta1": -0.5 * math.pi},
        ]
    )

    distance = space.minimum_topology_distance(space.to_unit(physical[:1]), space.to_unit(physical[1:]))

    assert space.direction_groups == [
        {"kind": "projective", "indices": [space.names.index("theta0"), space.names.index("theta1")]}
    ]
    assert distance[0] == pytest.approx(0.0, abs=2.0e-8)


# Verify the joint spherical resonance group wraps both circular coordinates
def test_space_joint_polar_sphere_topology():
    bounds = {
        "norm": {"type": "uniform", "lower": 0.0, "upper": 3.0},
        "phi": {"type": "uniform", "lower": -math.pi, "upper": math.pi},
        "theta0": {"type": "uniform", "lower": -0.5 * math.pi, "upper": 0.5 * math.pi},
        "theta1": {"type": "uniform", "lower": -math.pi, "upper": math.pi},
    }
    topology = {
        "schema_version": 1,
        "groups": [
            {
                "kind": "polar_sphere",
                "base": "g_ls",
                "parameters": ["norm", "phi", "theta0", "theta1"],
                "period": 2.0 * math.pi,
            }
        ],
    }
    space = IceboSpace(bounds, parameter_topology=topology, device=torch.device("cpu"), dtype=torch.float64)
    first = {"norm": 2.0, "phi": -2.0, "theta0": 0.31, "theta1": -0.47}
    direction = tools.spherical_vector_from_angles([first["theta0"], first["theta1"]])
    opposite = tools.spherical_angles_from_vector([-value for value in direction])
    second = {"norm": 2.0, "phi": -2.0 + math.pi, "theta0": opposite[0], "theta1": opposite[1]}
    physical = space.configs_to_tensor([first, second])

    distance = space.minimum_topology_distance(space.to_unit(physical[:1]), space.to_unit(physical[1:]))

    assert bool(space.periodic_mask[space.names.index("phi")])
    assert bool(space.periodic_mask[space.names.index("theta1")])
    assert space.direction_groups == [
        {"kind": "sphere", "indices": [space.names.index("theta0"), space.names.index("theta1")]}
    ]
    assert distance[0] == pytest.approx(0.0, abs=2.0e-8)


# Verify an explicit non-radian phase period identifies its branch cut
def test_space_explicit_phase_period():
    bounds = {"phase_turn": {"type": "uniform", "lower": 0.0, "upper": 1.0}}
    topology = {
        "schema_version": 1,
        "groups": [
            {
                "kind": "phase",
                "base": "phase_turn",
                "parameters": ["phase_turn"],
                "period": 1.0,
            }
        ],
    }
    space = IceboSpace(
        bounds,
        parameter_topology=topology,
        device=torch.device("cpu"),
        dtype=torch.float64,
    )
    values = space.configs_to_tensor([{"phase_turn": 0.0}, {"phase_turn": 1.0}])
    unit = space.to_unit(values)

    distance = space.minimum_topology_distance(unit[:1], unit[1:])

    assert distance[0] == pytest.approx(0.0, abs=1.0e-12)


# Verify ICEBO retains the global sign of a spherical direction
def test_space_separates_spherical_antipodes():
    bounds = {
        "theta0": {"type": "uniform", "lower": -0.5 * math.pi, "upper": 0.5 * math.pi},
        "theta1": {"type": "uniform", "lower": -math.pi, "upper": math.pi},
    }
    topology = {
        "schema_version": 1,
        "groups": [
            {
                "kind": "sphere",
                "base": "g_ls",
                "parameters": ["theta0", "theta1"],
            }
        ],
    }
    space = IceboSpace(
        bounds,
        parameter_topology=topology,
        device=torch.device("cpu"),
        dtype=torch.float64,
    )
    physical = space.configs_to_tensor(
        [{"theta0": 0.0, "theta1": 0.0}, {"theta0": 0.0, "theta1": math.pi}]
    )
    unit = space.to_unit(physical)

    distance = space.minimum_topology_distance(unit[:1], unit[1:])

    assert distance[0] == pytest.approx(2.0)


# Verify Sobol warmup directions follow the uniform spherical measure
def test_space_draws_uniform_spherical_directions():
    bounds = {
        "theta0": {"type": "uniform", "lower": -0.5 * math.pi, "upper": 0.5 * math.pi},
        "theta1": {"type": "uniform", "lower": -math.pi, "upper": math.pi},
    }
    topology = {
        "schema_version": 1,
        "groups": [
            {
                "kind": "sphere",
                "base": "g_ls",
                "parameters": ["theta0", "theta1"],
            }
        ],
    }
    space = IceboSpace(
        bounds,
        parameter_topology=topology,
        device=torch.device("cpu"),
        dtype=torch.float64,
    )
    angles = space.from_unit(space.sobol(4096, seed=41), round_integers=False)
    theta0, theta1 = angles.T
    vectors = torch.stack(
        (
            torch.cos(theta0) * torch.cos(theta1),
            torch.sin(theta0),
            torch.cos(theta0) * torch.sin(theta1),
        ),
        dim=1,
    )

    assert torch.mean(vectors.square(), dim=0).numpy() == pytest.approx(
        [1.0 / 3.0] * 3, abs=4.0e-3
    )


# Verify spherical TuRBO regions cross and canonicalize the final angle seam
def test_icebo_turbo_wraps_spherical_seam():
    bounds = {
        "theta0": {"type": "uniform", "lower": -0.5 * math.pi, "upper": 0.5 * math.pi},
        "theta1": {"type": "uniform", "lower": -math.pi, "upper": math.pi},
    }
    topology = {
        "schema_version": 1,
        "groups": [
            {
                "kind": "sphere",
                "base": "g_ls",
                "parameters": ["theta0", "theta1"],
            }
        ],
    }
    space = IceboSpace(
        bounds,
        parameter_topology=topology,
        device=torch.device("cpu"),
        dtype=torch.float64,
    )
    center = torch.tensor([0.5, 0.98], dtype=torch.float64)
    lower, upper = turbo_bounds(
        center, torch.ones(1, dtype=torch.float64), space.geometry, 0.2
    )
    wrapped = space.canonical_unit(torch.tensor([0.5, 1.05], dtype=torch.float64))

    assert lower.numpy() == pytest.approx([0.4, 0.88])
    assert upper.numpy() == pytest.approx([0.6, 1.08])
    assert wrapped.numpy() == pytest.approx([0.5, 0.05])


# Verify phase topology rejects an implicit period
def test_space_requires_explicit_phase_period():
    topology = {
        "schema_version": 1,
        "groups": [{"kind": "phase", "base": "phi", "parameters": ["phi"]}],
    }

    with pytest.raises(ValueError, match="period"):
        IceboSpace(
            {"phi": {"type": "uniform", "lower": -math.pi, "upper": math.pi}},
            parameter_topology=topology,
            device=torch.device("cpu"),
            dtype=torch.float64,
        )


# Verify pending reservations use the same phase identification as the GP
def test_pending_deduplication_phase_topology():
    bounds = {"phi": {"type": "uniform", "lower": -math.pi, "upper": math.pi}}
    topology = {
        "schema_version": 1,
        "groups": [
            {
                "kind": "phase",
                "base": "phi",
                "parameters": ["phi"],
                "period": 2.0 * math.pi,
            }
        ],
    }
    optimizer = _fast_optimizer(bounds, parameter_topology=topology)

    optimizer.reserve({"phi": -math.pi})
    optimizer.reserve({"phi": math.pi})
    assert optimizer.diagnostics()["pending"] == 1

    optimizer.observe({"phi": math.pi}, 1.0)
    assert optimizer.diagnostics()["pending"] == 0


# Verify mixed integer proposals retain typed values and maximize direction works
def test_mixed_integer_maximize_direction():
    bounds = {
        "n": {"type": "int", "lower": 0, "upper": 4},
        "x": {"type": "uniform", "lower": -1.0, "upper": 1.0},
    }
    optimizer = _fast_optimizer(
        bounds,
        direction="maximize",
        warmup=4,
    )
    initial = optimizer.suggest(4)
    values = [-((row["n"] - 2) ** 2) - (row["x"] - 0.1) ** 2 for row in initial]
    optimizer.observe(initial, values)

    proposals = optimizer.suggest(3)

    assert all(isinstance(row["n"], int) for row in proposals)
    assert all(0 <= row["n"] <= 4 and -1.0 <= row["x"] <= 1.0 for row in proposals)
    assert optimizer.diagnostics()["direction"] == "maximize"


# Verify automatic device selection reports the available accelerated backend
def test_icebo_auto_device_selection():
    bounds = {"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}
    optimizer = _fast_optimizer(bounds)
    expected = "cuda" if torch.cuda.is_available() else "cpu"

    assert optimizer.diagnostics()["device"].startswith(expected)


# Verify affine standardization exactly propagates values and variances
def test_icebo_standardization_round_trip():
    values = torch.tensor(
        [-3.0, -1.0, -0.1, 0.0, 0.2, 2.0, 20.0],
        dtype=torch.float64,
    )
    transform = StandardizeOutput()
    transform.fit(values)
    transformed = transform.forward(values)
    restored = transform.inverse(transformed)
    variance = torch.tensor([0.2, 0.7], dtype=torch.float64)

    assert torch.all(torch.diff(transformed) > 0.0)
    assert torch.allclose(restored, values, rtol=1.0e-8, atol=1.0e-8)
    assert torch.allclose(
        transform.variance_inverse(transform.variance_forward(variance)),
        variance,
        rtol=1.0e-10,
        atol=1.0e-10,
    )


# Verify interleaved optimizer instances do not share random state
def test_icebo_instances_are_rng_independent():
    bounds = {"x": {"type": "uniform", "lower": -1.0, "upper": 1.0}}
    first = _fast_optimizer(bounds, seed=31, warmup=3)
    second = _fast_optimizer(bounds, seed=31, warmup=3)

    first_initial = first.suggest(3)
    second_initial = second.suggest(3)
    values = [(row["x"] - 0.2) ** 2 for row in first_initial]
    first.observe(first_initial, values)
    second.observe(second_initial, values)
    first_next = first.suggest()[0]
    second_next = second.suggest()[0]

    assert first_initial == second_initial
    assert first_next["x"] == pytest.approx(second_next["x"], abs=1.0e-12)


# Check icetune does not import HEBO for the icebo algorithm
def test_icetune_search_state_icebo_without_hebo(monkeypatch):
    original_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name == "hebo" or name.startswith("hebo."):
            raise AssertionError("ICEBO attempted to import HEBO")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    args = SimpleNamespace(
        algorithm="icebo",
        cost="loss",
        rand_trials=2,
        rngseed=9,
        surrogate_fit_percentile=0.25,
        icebo_settings=_icebo_settings(
            {
                "gp.fit_steps": 5,
                "gp.fit_restarts": 1,
                "acquisition.raw_samples": 128,
                "acquisition.mc_samples": 32,
            }
        ),
    )
    bounds = {"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}
    search = SearchState(args=args, bounds=bounds, initial_points=None)
    records = [
        {
            "trial_id": f"trial-{index:06d}",
            "config": {"x": value},
            "metrics": {"loss": (value - 0.3) ** 2},
        }
        for index, value in enumerate([0.0, 0.4, 0.8])
    ]
    search.observe(records)

    config, payload = search.ask(3)

    assert 0.0 <= config["x"] <= 1.0
    assert payload["kind"] == "icebo"
    assert search.icebo.observation_count == 3


# Verify asynchronous results update TuRBO once per complete batch
def test_turbo_waits_for_batch():
    args = SimpleNamespace(
        algorithm="icebo",
        cost="loss",
        rand_trials=0,
        rngseed=9,
        surrogate_fit_percentile=1.0,
        icebo_settings=_icebo_settings(
            {
                "turbo.enabled": True, "gp.fit_steps": 5,
                "gp.fit_restarts": 1,
                "acquisition.raw_samples": 128,
                "acquisition.mc_samples": 32,
            }
        ),
    )
    search = SearchState(
        args=args,
        bounds={"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}},
        initial_points=None,
    )
    records = [
        {
            "trial_id": f"trial-{index:06d}",
            "config": {"x": value},
            "metrics": {"loss": loss},
            "search_payload": {
                "kind": "icebo",
                "proposal_batch_id": 0,
                "proposal_batch_size": 2,
                "proposal_index": index,
            },
        }
        for index, (value, loss) in enumerate(((0.2, 2.0), (0.8, 1.0)))
    ]

    search.observe(records[:1])
    assert search.icebo.diagnostics()["turbo"]["best_value"] is None

    search.observe(records)
    search.ask(2)
    turbo = search.icebo.diagnostics()["turbo"]
    assert turbo["best_value"] == pytest.approx(min(search.icebo.predict([row["config"] for row in records])["mean"]))
    assert turbo["batch_size"] == 2


# Verify the logarithmic positive part retains finite gradients in the far negative tail
def test_log_fatplus_stable_in_extreme_tail():
    argument = torch.tensor([-1000.0], dtype=torch.float64, requires_grad=True)
    value = log_fatplus(argument)
    value.backward()

    assert torch.isfinite(value).all()
    assert torch.isfinite(argument.grad).all()
    assert argument.grad[0] > 0.0


# Verify continuous observations reject non-finite and out-of-bounds values
@pytest.mark.parametrize("value", [-1.0e-6, 1.000001, math.inf, math.nan])
def test_invalid_continuous_observations(value):
    optimizer = _fast_optimizer(
        {"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}
    )

    with pytest.raises(ValueError):
        optimizer.observe({"x": value}, 0.0)


# Verify canonical TuRBO expansion, contraction and restart transitions
def test_icebo_turbo_state_transitions():
    state = TurboState.from_settings(2, _icebo_settings()["turbo"])
    state.best_value = 10.0
    for value in (9.0, 8.0, 7.0):
        state.update(torch.tensor([value]))
    assert state.length == pytest.approx(1.6)

    for _ in range(state.failure_tolerance):
        state.update(torch.tensor([7.0]))
    assert state.length == pytest.approx(0.8)

    state.length = 0.5 * state.length_min
    state.update(torch.tensor([7.0]))
    assert state.restart_triggered is True


# Verify QMC Normal samples are deterministic without global RNG mutation
def test_icebo_qmc_normals_are_replayable():
    first = sobol_normal(
        16,
        5,
        seed=91,
        device=torch.device("cpu"),
        dtype=torch.float64,
    )
    second = sobol_normal(
        16,
        5,
        seed=91,
        device=torch.device("cpu"),
        dtype=torch.float64,
    )

    assert torch.equal(first, second)
    assert torch.all(torch.isfinite(first))


# Verify replayable ICEBO state preserves deterministic TuRBO counters
def test_icebo_state_round_trip():
    bounds = {"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}
    first = _fast_optimizer(bounds, trust_region=True)
    first._turbo.length = 0.2
    first._turbo.failure_counter = 3
    first.set_replay_proposal_count(17)
    state = first.state_dict()
    second = _fast_optimizer(bounds, trust_region=True)

    second.load_state_dict(state)

    assert second.state_dict() == state


# Reproduce adaptive proposals after JSON restore and a different observation arrival order
def test_fitted_state_reproduces_proposals():
    bounds = {"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}
    first = _fast_optimizer(bounds, warmup=4)
    rows = first.suggest(4)
    values = [(row["x"] - 0.3)**2 for row in rows]
    first.observe(rows, values)
    pending = first.suggest(2)
    state = json.loads(json.dumps(first.state_dict()))
    second = _fast_optimizer(bounds, warmup=4)
    second.observe(rows[::-1], values[::-1], update_turbo=False)
    second.set_pending(pending)
    second.load_state_dict(state)
    assert second.suggest(2) == first.suggest(2)
    assert second.state_dict() == first.state_dict()
    values = [(row["x"] - 0.3)**2 for row in pending]
    first.observe(pending, values)
    second.observe(pending[::-1], values[::-1])
    assert second.suggest() == first.suggest()
    assert second.state_dict() == first.state_dict()


# Reevaluate the incumbent under the same noisy posterior used for a new batch
def test_turbo_reassesses_noisy_incumbent():
    optimizer = _fast_optimizer({"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}, warmup=2, trust_region=True)
    optimizer.observe([{"x": 0.2}, {"x": 0.8}], [0.0, 1.0], errors=[0.5, 0.5])
    pending = optimizer.suggest()
    optimizer._turbo.best_value = -1e6
    optimizer.observe(pending, [0.5], errors=[0.5])
    optimizer.predict(pending)
    assert optimizer.diagnostics()["turbo"]["best_value"] > -10.0


# Release every completed reservation even when its objective is invalid
@pytest.mark.parametrize("values", [[math.nan, math.inf], [1.0, math.nan]])
def test_failed_observations_release_pending(values):
    optimizer = _fast_optimizer({"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}})
    configs = optimizer.suggest(2)
    optimizer.observe(configs, values)
    assert optimizer.diagnostics()["pending"] == 0
    assert optimizer.observation_count == sum(math.isfinite(value) for value in values)


# Apply a fitted TuRBO restart before selecting the next acquisition policy
def test_restart_starts_sobol_design_immediately():
    optimizer = _fast_optimizer({"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}, warmup=2, trust_region=True)
    configs = optimizer.suggest(2)
    optimizer.observe(configs, [1.0, 1.0])
    adaptive = optimizer.suggest()
    optimizer._turbo.length = optimizer._turbo.length_min
    optimizer._turbo.failure_counter = optimizer._turbo.failure_tolerance - 1
    optimizer.observe(adaptive, [1.0])
    _, diagnostics = optimizer.suggest(return_diagnostics=True)
    assert optimizer.diagnostics()["turbo"]["restart_count"] == 1
    assert diagnostics[0]["acquisition"] == "sobol_warmup"
    assert optimizer.diagnostics()["restart_design_remaining"] == optimizer.turbo_settings["restart_min_points"] - 1


# Keep the initial design separate from adaptation and wait for its pending trials
def test_initial_design_waits_without_adapting():
    optimizer = _fast_optimizer({"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}, warmup=4, trust_region=True)
    rows = optimizer.suggest(8)
    assert len(rows) == 4
    assert optimizer.suggest() == []
    for row in rows:
        optimizer.observe(row, (row["x"] - 0.3)**2)
    optimizer.suggest()
    state = optimizer.diagnostics()["turbo"]
    assert state["success_counter"] == state["failure_counter"] == state["restart_count"] == 0


# Drain an old region before fitting an independent GP to the new initial design
def test_restart_waits_pending_fits_local_data():
    optimizer = _fast_optimizer({"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}, warmup=2, trust_region=True)
    initial = optimizer.suggest(2)
    optimizer.observe(initial, [1.0, 1.0])
    pending = optimizer.suggest(2)
    optimizer._turbo.length = optimizer._turbo.length_min
    optimizer._turbo.failure_counter = optimizer._turbo.failure_tolerance - 1
    optimizer.observe(pending[:1], [1.0])
    assert optimizer.suggest() == []
    assert optimizer.diagnostics()["turbo"]["restart_triggered"]
    optimizer.observe(pending[1:], [float("nan")])
    fresh = optimizer.suggest(8)
    assert len(fresh) == optimizer.turbo_settings["restart_min_points"]
    assert optimizer.suggest() == []
    optimizer.observe(fresh, [10.0 + row["x"] for row in fresh])
    _, diagnostics = optimizer.suggest(return_diagnostics=True)
    assert diagnostics[0]["acquisition"] == "qlognei"
    assert optimizer.diagnostics()["fit"]["training_points"] == len(fresh)
    assert optimizer._transform.location > 10.0
    assert optimizer.diagnostics()["turbo"]["success_counter"] == 0
    assert optimizer.diagnostics()["turbo"]["failure_counter"] == 0


# Resume actual Ray background fitting across failures and trust region restarts
def test_icebo_ray_restart_resume(tmp_path):
    import time

    from core.tune.backends.ray import PercentileSearch
    from ray.tune.search import Searcher

    settings = _icebo_settings({"gp.fit_steps": 12, "gp.fit_restarts": 1,
                               "acquisition.raw_samples": 64, "acquisition.mc_samples": 16,
                               "acquisition.restarts": 2, "acquisition.gradient_steps": 4,
                               "turbo.enabled": True, "turbo.length_initial": 0.01, "turbo.failure_base": 1.0,
                               "turbo.restart_min_points": 2, "turbo.restart_points_per_dimension": 1})
    args = SimpleNamespace(algorithm="icebo", cost="loss", rand_trials=2, rngseed=9,
                           surrogate_fit_percentile=1.0, icebo_settings=settings,
                           max_concurrent_trials=2, num_trials=14, proposal_batch_size=1,
                           proposal_batch_fraction=1.0)
    arguments = dict(args=args, bounds={"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}, initial_points=None)
    search = PercentileSearch(**arguments)
    checkpoint = str(tmp_path / "icebo.pkl")
    live, resumed, failed = {}, False, False
    deadline = time.monotonic() + 60
    try:
        while time.monotonic() < deadline:
            while len(live) < args.max_concurrent_trials:
                trial_id = f"worker-{search.next_index:06d}"
                config = search.suggest(trial_id)
                if config is None or config == Searcher.FINISHED:
                    break
                live[trial_id] = config
            if live:
                trial_id = next(reversed(live))
                live.pop(trial_id)
                payload = search.live_records[trial_id]["search_payload"]
                fail = not failed and payload.get("acquisition") == "qlognei"
                search.on_trial_complete(trial_id, result=None if fail else {"loss": 1.0}, error=fail)
                failed |= fail
            elif search.next_index >= args.num_trials:
                break
            else:
                time.sleep(0.01)
            if search.icebo_state is not None and search.icebo_state["turbo"]["restart_count"] > 0 and not resumed:
                search.save(checkpoint)
                if search.fit_executor is not None:
                    search.fit_executor.shutdown(wait=True)
                search = PercentileSearch(**arguments)
                search.restore(checkpoint)
                resumed = True
        else:
            pytest.fail("Ray icebo search stalled during restart")
        assert resumed and failed
        assert len(search.completed_records) + len(search.failed_records) == args.num_trials
        assert search.icebo_state["turbo"]["restart_count"] > 0
    finally:
        if search.fit_executor is not None:
            search.fit_executor.shutdown(wait=True)


# An empty optimizer needs an initial design even when warmup is zero
def test_zero_warmup_without_observations():
    optimizer = _fast_optimizer({"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}, warmup=0)
    configs, diagnostics = optimizer.suggest(return_diagnostics=True)
    assert len(configs) == 1
    assert diagnostics[0]["acquisition"] == "sobol_warmup"


# Joint raw designs must cover the full batch space without correlations between candidates
def test_icebo_raw_batches_cover_joint_space():
    from core.tune.optimizers.icebo.qlognei import QLogNEIOptimizer

    optimizer = _fast_optimizer({"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}})
    acquisition = QLogNEIOptimizer(optimizer.space, settings=optimizer._acquisition_settings(), seed=13)
    batches = acquisition._raw_batches(2, torch.zeros(1), torch.ones(1), iteration=0)
    cells = torch.floor(batches[:, :, 0] * 4).to(torch.int64)
    assert torch.unique(cells, dim=0).shape[0] == 16


# Denoise the actual completed batch when other evaluations finish between its members
def test_turbo_interleaved_batch_coords():
    optimizer = _fast_optimizer({"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}, trust_region=True)
    rows = [{"x": 0.1}, {"x": 0.5}, {"x": 0.9}]
    optimizer.observe(rows, [4.0, -10.0, 5.0], update_turbo=False)
    optimizer.update_turbo_batch([rows[0], rows[2]])
    predicted = optimizer.predict(rows)
    assert optimizer.diagnostics()["turbo"]["best_value"] == pytest.approx(min(predicted["mean"][[0, 2]]))
    assert optimizer.diagnostics()["turbo"]["best_value"] > predicted["mean"][1]


# Preserve physical bound precision when geometry is evaluated in float64
def test_geometry_double_precision_bounds():
    optimizer = _fast_optimizer({"x": {"type": "uniform", "lower": 0.1, "upper": 0.4}})
    geometry = optimizer.space.geometry
    points = torch.tensor([[0.1], [0.4]], dtype=torch.float64)
    distance = geometry.squared_distance(points, points, torch.ones(1, dtype=torch.float64))
    assert distance[0, 1].item() == pytest.approx(1.0, abs=1.0e-14)


# Evaluate shared geometry without promoting single precision surrogate inputs
@pytest.mark.parametrize("dtype", [torch.float32, torch.float64])
def test_geometry_input_dtype(dtype):
    from core.tune.parameters.kernel import ProductManifoldGeometry

    geometry = ProductManifoldGeometry(1, (), ((0.1, 0.4),))
    points = torch.tensor([[0.1], [0.4]], dtype=dtype)
    distance = geometry.squared_distance(points, points, torch.ones(1, dtype=dtype))
    assert distance.dtype == dtype
    assert distance[0, 1].item() == pytest.approx(1.0, abs=4.0 * torch.finfo(dtype).eps)


# Reject numerical controls that reverse adaptation or make its recovery ineffective
@pytest.mark.parametrize("key,value", [
    ("gp.jitter_multiplier", 1.0),
    ("acquisition.prune_baseline", 0.5),
    ("turbo.failure_length_multiplier", 1.0),
    ("turbo.success_length_multiplier", 1.0),
])
def test_icebo_rejects_inconsistent_controls(key, value):
    with pytest.raises(ValueError):
        _icebo_settings({key: value})


# Construct a real posterior with pending evaluations for acquisition checks
@pytest.fixture
def acquisition_case():
    optimizer = _fast_optimizer({"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}, warmup=0)
    optimizer.observe([{"x": x} for x in (0.1, 0.5, 0.9)], [0.16, 0.0, 0.16], [0.02, 0.01, 0.03])
    optimizer.predict({"x": 0.3})
    return optimizer


# Bound smoothing error against the exact sampled noisy batch improvement
@pytest.mark.parametrize("q", [1, 4])
def test_icebo_qlognei_smoothing_error(acquisition_case, q):
    model = acquisition_case._surrogate
    settings = acquisition_case._acquisition_settings()
    acquisition = QLogNoisyExpectedImprovement(
        model, pending=model.X[:0], q=q, settings=settings, seed=31,
    )
    points = model.X[1:2].expand(q, -1).clone().requires_grad_(True)
    samples = acquisition._conditional_samples(points.unsqueeze(0))
    exact = (acquisition.reference_best - samples.amin(-1)).clamp_min(0).mean()
    value = acquisition(points).exp().squeeze()
    bound = q ** settings.tau_max * (exact + settings.tau_relu * (math.log(2.0) + 0.1))
    assert value >= exact - 1.0e-12
    assert value <= bound + 1.0e-12
    gradient, = torch.autograd.grad(value, points)
    assert torch.isfinite(gradient).all()


# Compare cached conditional moments and gradients against the dense GP equations
def test_icebo_qlognei_cached_conditioning(acquisition_case):
    model = acquisition_case._surrogate
    acquisition = QLogNoisyExpectedImprovement(
        model, pending=model.X.new_tensor([[0.3]]), q=2,
        settings=acquisition_case._acquisition_settings(), seed=23,
    )
    points = model.X.new_tensor([[[0.25], [0.65]], [[0.35], [0.75]]], requires_grad=True)
    direct = model.posterior_cross_covariance(points, acquisition.reference)
    cached = model.kernel(points, acquisition.reference) - model.kernel(points, model.X) @ acquisition.reference_solve
    torch.testing.assert_close(cached, direct, rtol=1.0e-10, atol=1.0e-12)
    torch.testing.assert_close(torch.autograd.grad(cached.sum(), points)[0], torch.autograd.grad(direct.sum(), points)[0])
    assert torch.autograd.gradcheck(acquisition, (points,), atol=1.0e-4, rtol=1.0e-3)
    together = acquisition(points)
    separate = torch.cat([acquisition(row) for row in points])
    torch.testing.assert_close(together, separate)


# Keep an initial random design from contracting TuRBO before adaptive proposals
def test_icetune_icebo_random_design_trust_region():
    optimizer = _fast_optimizer({"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}, trust_region=True)
    args = SimpleNamespace(algorithm="icebo", cost="loss", rand_trials=40, rngseed=9,
                           surrogate_fit_percentile=1.0, icebo_settings=optimizer.settings)
    search = SearchState(args=args, bounds=optimizer.space.bounds, initial_points=None)
    search.observe([
        {"trial_id": f"trial-{index:06d}", "config": {"x": float(x)}, "metrics": {"loss": 1.0},
         "search_payload": {"kind": "random"}}
        for index, x in enumerate(np.linspace(0.0, 1.0, args.rand_trials))
    ])
    _, payload = search.ask(args.rand_trials)
    state = search.icebo.diagnostics()["turbo"]
    assert payload["acquisition"] == "qlognei"
    assert state["restart_count"] == state["failure_counter"] == state["success_counter"] == 0
    assert state["length"] == pytest.approx(optimizer.turbo_settings["length_initial"])
    assert state["best_value"] == pytest.approx(1.0)


# Keep gradient refinement invariant under its memory batch limit
def test_acquisition_refinement_batch_limit(acquisition_case):
    from dataclasses import replace

    from core.tune.optimizers.icebo.qlognei import QLogNEIOptimizer

    model = acquisition_case._surrogate
    settings = acquisition_case._acquisition_settings()
    acquisition = QLogNoisyExpectedImprovement(model, pending=model.X[:0], q=1, settings=settings, seed=23)
    seeds = model.X.new_tensor([[[0.2]], [[0.7]], [[0.8]]])
    lower, upper = model.X.new_zeros(1), model.X.new_ones(1)
    full = QLogNEIOptimizer(acquisition_case.space, settings=settings, seed=1)
    chunked = QLogNEIOptimizer(acquisition_case.space, settings=replace(settings, optimization_batch_size=1), seed=1)
    torch.testing.assert_close(full._refine(seeds, acquisition, lower, upper),
                               chunked._refine(seeds, acquisition, lower, upper), atol=1.0e-7, rtol=1.0e-6)


# Restore distributed TuRBO counters without replaying completed batches twice
def test_icetune_icebo_restores_issued_trust_region():
    optimizer = _fast_optimizer({"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}, trust_region=True)
    args = SimpleNamespace(algorithm="icebo", cost="loss", rand_trials=0, rngseed=9,
                           surrogate_fit_percentile=1.0, icebo_settings=optimizer.settings)
    first = SearchState(args=args, bounds=optimizer.space.bounds, initial_points=None)
    records = [
        {"trial_id": f"trial-{index:06d}", "config": {"x": x}, "metrics": {"loss": (x - 0.4) ** 2},
         "search_payload": {"kind": "random"}}
        for index, x in enumerate((0.0, 0.3, 0.7, 1.0))
    ]
    first.observe(records)
    config, payload = first.ask(4)
    records.append({"trial_id": "trial-000004", "config": config, "metrics": {"loss": 10.0}, "search_payload": payload})
    first.observe(records)
    config, payload = first.ask(5)
    issued = [*records, {"trial_id": "trial-000005", "config": config, "search_payload": payload}]
    state = first.icebo.state_dict()
    assert state["turbo_batches"]

    second = SearchState(args=args, bounds=optimizer.space.bounds, initial_points=None)
    second.observe(records)
    second.restore_issued(list(reversed(issued)))
    second.icebo.predict({"x": 0.4})
    assert second.icebo.state_dict() == state


# Keep Hellinger distances and derivatives finite at zero simplex probabilities
@pytest.mark.parametrize("dtype", [torch.float32, torch.float64])
def test_simplex_distance_boundary_gradients(dtype):
    from core.tune.parameters.kernel import ProductManifoldGeometry

    geometry = ProductManifoldGeometry(2, ({"kind": "simplex", "indices": (0, 1)},),
                                       ((0.0, math.pi / 2),) * 2).to(dtype=dtype)
    points = torch.tensor([[0.0, 0.0], [0.3, 0.5]], dtype=dtype, requires_grad=True)
    scale = torch.ones(1, dtype=dtype, requires_grad=True)
    distance = geometry.squared_distance(points, points, scale)
    first = torch.autograd.grad(distance.sum(), (points, scale), create_graph=True)
    second = torch.autograd.grad(sum(gradient.sum() for gradient in first), (points, scale))
    assert all(torch.isfinite(gradient).all() for gradient in (*first, *second))
    expected = 1.0 - math.cos(0.3)
    assert distance[0, 1].item() == pytest.approx(expected, abs=4 * torch.finfo(dtype).eps)
    # Preserve the inward derivative at the first-orthant boundary
    assert first[0][0, 0].item() == pytest.approx(-2 * math.sin(0.3) * math.cos(0.5),
                                                abs=4 * torch.finfo(dtype).eps)


# Honor a disabled acquisition gradient search with the actual fitted surrogate
def test_zero_gradient_steps_seeds():
    from core.tune.optimizers.icebo.qlognei import QLogNEIOptimizer, QLogNoisyExpectedImprovement

    optimizer = _fast_optimizer({"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}, gradient_steps=0)
    optimizer.observe([{"x": 0.0}, {"x": 1.0}], [1.0, 0.0])
    optimizer.predict([{"x": 0.5}])
    settings = optimizer._acquisition_settings()
    acquisition = QLogNoisyExpectedImprovement(optimizer._surrogate, pending=optimizer._pending,
                                               q=1, settings=settings, seed=11)
    maximizer = QLogNEIOptimizer(optimizer.space, settings=settings, seed=11)
    seeds = torch.tensor([[[0.3]], [[0.6]]], dtype=optimizer.dtype)
    result = maximizer._refine(seeds, acquisition, seeds.new_zeros(1), seeds.new_ones(1))
    torch.testing.assert_close(result, seeds, atol=0, rtol=0)


# Continue a fitted GP without random restarts while preserving posterior predictions
def test_exact_gp_converged_warm_start():
    x = torch.linspace(0.0, 1.0, 16, dtype=torch.float64)[:, None]
    y = torch.sin(4.0 * x[:, 0])
    settings = {**_icebo_settings()["gp"], "fit_steps": 200}
    model = ExactGaussianProcess(x, y, settings=settings, dtype=torch.float64, device="cpu", seed=7)
    initial = model.fit()
    mean, covariance = model.posterior_covariance(x)
    continued = model.fit(warm_start=True)
    new_mean, new_covariance = model.posterior_covariance(x)
    assert initial["gradient_evaluations"] < initial["evaluations"]
    assert continued["restarts"] == continued["converged_restarts"] == 1
    assert continued["best_negative_log_posterior"] <= initial["best_negative_log_posterior"] + 1e-8
    torch.testing.assert_close(new_mean, mean, atol=1e-4, rtol=1e-4)
    torch.testing.assert_close(new_covariance, covariance, atol=1e-5, rtol=1e-4)


# Exhaust configured recovery restarts if a short continuation cannot converge
def test_exact_gp_unconverged_warm_start_retries():
    x = torch.linspace(0.0, 1.0, 16, dtype=torch.float64)[:, None]
    y = torch.sin(4.0 * x[:, 0])
    settings = _icebo_settings()["gp"]
    model = ExactGaussianProcess(x, y, settings=settings, dtype=torch.float64, device="cpu", seed=7)
    result = model.fit(steps=1, warm_start=True)
    assert result["converged_restarts"] == 0
    assert result["restarts"] == settings["fit_restarts"]
    assert torch.isfinite(model.posterior(x).mean).all()

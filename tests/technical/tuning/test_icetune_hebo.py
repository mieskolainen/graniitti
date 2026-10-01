# Tests for icetune HEBO topology kernels
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import math
import random
import runpy
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch

ROOT = Path(__file__).resolve().parents[3]

pytest.importorskip("hebo")

from core.tune.optimizers.hebo import optimizer as icetune_hebo
from core.tune.optimizers.hebo.config import load_hebo_config, validate_hebo_config
from core.tune.parameters import tools, topology
from core.tune.parameters.kernel import ProductManifoldGeometry
from core.tune.search import SearchState


# Build a mixed ordinary and phase HEBO design space
def _design():
    from hebo.design_space.design_space import DesignSpace

    return DesignSpace().parse(
        [
            {"name": "scale", "type": "num", "lb": 2.0, "ub": 6.0},
            {"name": "phase", "type": "num", "lb": -math.pi, "ub": math.pi},
        ]
    )


# Compute topology metadata marking the scalar physical phase
def _topology():
    return {
        "schema_version": 1,
        "groups": [
            {
                "kind": "phase",
                "base": "amplitude",
                "parameters": ["phase"],
                "period": 2.0 * math.pi,
            },
        ],
    }


# Build the direct phase topology GP with one ordinary coordinate
def _phase_gp(scale=(2.0, 6.0), **conf):
    conf.setdefault("num_epochs", 1)
    return icetune_hebo.TopologyGP(
        2,
        0,
        1,
        topology_groups=[{"kind": "phase", "indices": [1], "period": 2.0 * math.pi}],
        numeric_bounds=[scale, [-math.pi, math.pi]],
        **conf,
    )


# Build the common optimizer arguments for one seeded HEBO toy
def _search_args(warmup, percentile=1.0):
    return SimpleNamespace(
        algorithm="hebo",
        cost="loss",
        rand_trials=warmup,
        rngseed=17,
        surrogate_fit_percentile=percentile,
    )


# Evaluate a two-dimensional quadratic with one circular coordinate
def _circular_loss(scale, phase):
    delta = (phase - (math.pi - 0.1) + math.pi) % (2.0 * math.pi) - math.pi
    return (scale - 4.0) ** 2 + delta**2


# Verify the direct kernel is continuous across the phase branch cut
def test_phase_kernel_continuous_branch_cut():
    model = _phase_gp()
    epsilon = 1.0e-5
    values = torch.tensor(
        [[4.0, -math.pi + epsilon], [4.0, math.pi - epsilon]],
        dtype=torch.float32,
    )

    covariance = model.model.conf["kern"](values).to_dense()
    correlation = covariance[0, 1] / torch.sqrt(covariance[0, 0] * covariance[1, 1])

    assert correlation.item() == pytest.approx(1.0, abs=1.0e-6)


# Verify one scalar phase does not add an internal surrogate dimension
def test_phase_model_input_dimension():
    model = _phase_gp()

    assert model.num_cont == 2
    assert model.model.num_cont == 2
    assert model.model.xscaler.__class__.__name__ == "TorchIdentityScaler"
    assert model.model.conf["kern"].base_kernel.effective_dim == 2


# Verify the exact model changes to sparse only above the configured limit
def test_topology_model_switches_sparse_above_limit():
    model = _phase_gp(
        max_exact_trials=2,
        sparse_config={"batch_size": 16, "num_epochs": 1, "num_inducing": 2},
    )
    categorical = torch.zeros((3, 0), dtype=torch.long)
    continuous = torch.tensor(
        [[3.0, -1.0], [4.0, 0.0], [5.0, 1.0]],
        dtype=torch.float32,
    )
    values = torch.tensor([[1.0], [0.0], [1.0]], dtype=torch.float32)

    model.fit(continuous[:2], categorical[:2], values[:2])
    assert model.model.__class__.__name__ == "GP"
    model.fit(continuous, categorical, values)
    assert model.model.__class__.__name__ == "SVGP"
    mean, variance = model.predict(continuous, categorical)
    assert mean.shape == values.shape
    assert variance.shape == values.shape
    assert torch.all(torch.isfinite(mean))
    assert torch.all(variance > 0.0)


# Verify the standard sparse model fits and predicts through the hybrid API
def test_standard_model_switches_sparse_above_limit():
    model = icetune_hebo.HybridGP(
        1,
        0,
        1,
        max_exact_trials=2,
        num_epochs=1,
        sparse_config={"batch_size": 16, "num_epochs": 1, "num_inducing": 2},
    )
    categorical = torch.zeros((3, 0), dtype=torch.long)
    continuous = torch.tensor([[0.0], [0.5], [1.0]], dtype=torch.float32)
    values = torch.tensor([[1.0], [0.0], [1.0]], dtype=torch.float32)

    model.fit(continuous, categorical, values)
    mean, variance = model.predict(continuous, categorical)

    assert model.model.__class__.__name__ == "SVGP"
    assert mean.shape == values.shape
    assert variance.shape == values.shape
    assert torch.all(torch.isfinite(mean))
    assert torch.all(variance > 0.0)


# Verify pending-aware synthetic rows do not trigger an early sparse switch
def test_hybrid_model_counts_only_canonical_trials():
    model = icetune_hebo.HybridGP(
        1,
        0,
        1,
        max_exact_trials=2,
        num_epochs=0,
        sparse_config={"batch_size": 16, "num_epochs": 0, "num_inducing": 2},
        synthetic_trials=1,
    )
    categorical = torch.zeros((4, 0), dtype=torch.long)
    continuous = torch.linspace(0.0, 1.0, 4).reshape(-1, 1)
    values = torch.square(continuous - 0.5)

    model.fit(continuous[:3], categorical[:3], values[:3])
    assert model.model.__class__.__name__ == "GP"
    model.fit(continuous, categorical, values)
    assert model.model.__class__.__name__ == "SVGP"


# Verify the production boundary is exact for standard and topology kernels
@pytest.mark.parametrize("topology", [False, True])
def test_hybrid_exact_sparse_limit(topology):
    limit = load_hebo_config()["gp"]["max_exact_trials"]
    sparse = {
        "batch_size": 4096,
        "num_epochs": 0,
        "num_inducing": 8,
        "verbose": False,
    }
    if topology:
        model = _phase_gp(num_epochs=0, sparse_config=sparse)
        continuous = torch.linspace(0.0, 1.0, limit + 1).reshape(-1, 1).repeat(1, 2)
        continuous[:, 0] = 2.0 + 4.0 * continuous[:, 0]
        continuous[:, 1] = -math.pi + 2.0 * math.pi * continuous[:, 1]
    else:
        model = icetune_hebo.HybridGP(
            1,
            0,
            1,
            num_epochs=0,
            sparse_config=sparse,
        )
        continuous = torch.linspace(0.0, 1.0, limit + 1).reshape(-1, 1)
    categorical = torch.zeros((limit + 1, 0), dtype=torch.long)
    values = torch.square(continuous[:, :1] - continuous[:1, :1])

    model.fit(continuous[:limit], categorical[:limit], values[:limit])
    assert model.model.__class__.__name__ == "GP"
    model.fit(continuous, categorical, values)
    assert model.model.__class__.__name__ == "SVGP"
    if topology:
        assert model.model.xscaler.__class__.__name__ == "TorchIdentityScaler"
        assert isinstance(
            model.model.conf["kern"].base_kernel,
            icetune_hebo.TopologyMaternKernel,
        )


# Verify fitted predictions identify the two representations of the phase seam
def test_phase_branch_cut():
    model = _phase_gp((0.0, 1.0), optimizer="adam", pred_likeli=False)
    Xc = torch.tensor(
        [[0.1, -3.0], [0.5, 0.0], [0.9, 3.0]],
        dtype=torch.float32,
    )
    Xe = torch.zeros((3, 0), dtype=torch.long)
    y = torch.tensor([[1.0], [0.0], [1.0]], dtype=torch.float32)
    model.fit(Xc, Xe, y)

    seam = torch.tensor(
        [[0.5, -math.pi], [0.5, math.pi]],
        dtype=torch.float32,
    )
    mean, variance = model.predict(seam, torch.zeros((2, 0), dtype=torch.long))

    assert mean[0].item() == pytest.approx(mean[1].item(), abs=1.0e-6)
    assert variance[0].item() == pytest.approx(variance[1].item(), abs=1.0e-6)


# Verify projective antipodes have unit direct-kernel correlation
def test_projective_kernel_antipodes():
    model = icetune_hebo.TopologyGP(
        3,
        0,
        1,
        topology_groups=[{"kind": "projective", "indices": [1, 2]}],
        numeric_bounds=[
            [2.0, 6.0],
            [-0.5 * math.pi, 0.5 * math.pi],
            [-0.5 * math.pi, 0.5 * math.pi],
        ],
        num_epochs=1,
    )
    angle = 0.37
    values = torch.tensor(
        [
            [4.0, angle, 0.5 * math.pi],
            [4.0, -angle, -0.5 * math.pi],
        ],
        dtype=torch.float64,
    )

    covariance = model.model.conf["kern"](values).to_dense()
    correlation = covariance[0, 1] / torch.sqrt(covariance[0, 0] * covariance[1, 1])

    assert correlation.item() == pytest.approx(1.0, abs=1.0e-12)
    assert model.model.num_cont == 3


# Verify the resonance phase compensates an antipodal signed-real LS direction
def test_polar_phase_direction_kernel():
    geometry = ProductManifoldGeometry(
        4,
        ({"kind": "polar_projective", "indices": (0, 1, 2, 3), "period": 2.0 * math.pi},),
        (
            (0.0, 3.0),
            (-math.pi, math.pi),
            (-0.5 * math.pi, 0.5 * math.pi),
            (-0.5 * math.pi, 0.5 * math.pi),
        ),
    )
    angle = 0.37
    first = torch.tensor([[2.0, 0.2, angle, 0.5 * math.pi]], dtype=torch.float64)
    same = torch.tensor([[2.0, 0.2 - math.pi, -angle, -0.5 * math.pi]], dtype=torch.float64)
    opposite = torch.tensor([[2.0, 0.2, -angle, -0.5 * math.pi]], dtype=torch.float64)
    scale = torch.ones(1, dtype=torch.float64)

    assert geometry.squared_distance(first, same, scale).item() == pytest.approx(0.0, abs=1.0e-12)
    assert geometry.squared_distance(first, opposite, scale).item() == pytest.approx(16.0 / 9.0)


# Verify zero coupling norm removes phase and LS direction coordinates
def test_polar_projective_kernel_collapses_zero_norm():
    geometry = ProductManifoldGeometry(
        3,
        ({"kind": "polar_projective", "indices": (0, 1, 2), "period": 2.0 * math.pi},),
        ((0.0, 3.0), (-math.pi, math.pi), (-0.5 * math.pi, 0.5 * math.pi)),
    )
    first = torch.tensor([[0.0, -2.7, -1.2]], dtype=torch.float64)
    second = torch.tensor([[0.0, 1.1, 0.8]], dtype=torch.float64)

    assert geometry.squared_distance(first, second, torch.ones(1, dtype=torch.float64)).item() == pytest.approx(0.0)


# Verify one Cartesian continuum coupling retains its physical complex sign
def test_cartesian_con_kernel_local_residue_sign():
    geometry = ProductManifoldGeometry(
        3,
        ({"kind": "cartesian_projective", "indices": (0, 1, 2), "weights": (2.0, 1.0)},),
        ((-3.0, 3.0), (-3.0, 3.0), (-0.5 * math.pi, 0.5 * math.pi)),
    )
    first = torch.tensor([[0.4, -0.8, 0.3]], dtype=torch.float64)
    same = -first
    scale = torch.ones(1, dtype=torch.float64)

    assert geometry.squared_distance(first, same, scale).item() > 0.0


# Verify the radial continuum chart removes only its duplicate projective orientation
def test_radial_con_antipodes_zero():
    geometry = ProductManifoldGeometry(
        3,
        ({"kind": "radial_projective", "indices": (0, 1, 2), "weights": (2.0, 2.0, 1.0)},),
        ((0.0, 6.0), (-0.5 * math.pi, 0.5 * math.pi), (-0.5 * math.pi, 0.5 * math.pi)),
    )
    angle = 0.37
    first = torch.tensor([[3.0, angle, 0.5 * math.pi]], dtype=torch.float64)
    same = torch.tensor([[3.0, -angle, -0.5 * math.pi]], dtype=torch.float64)
    zero_a = torch.tensor([[0.0, -0.8, 0.4]], dtype=torch.float64)
    zero_b = torch.tensor([[0.0, 0.7, -0.9]], dtype=torch.float64)
    scale = torch.ones(1, dtype=torch.float64)

    assert geometry.squared_distance(first, same, scale).item() == pytest.approx(0.0, abs=1.0e-12)
    assert geometry.squared_distance(zero_a, zero_b, scale).item() == pytest.approx(0.0, abs=1.0e-12)


# Verify polar and Cartesian common coefficients induce the same physical metric
def test_projective_polar_cartesian_distance():
    polar = ProductManifoldGeometry(
        3,
        ({"kind": "polar_projective", "indices": (0, 1, 2), "period": 2.0 * math.pi},),
        ((0.0, 3.0), (-math.pi, math.pi), (-0.5 * math.pi, 0.5 * math.pi)),
    )
    cartesian = ProductManifoldGeometry(
        3,
        ({"kind": "cartesian_projective", "indices": (0, 1, 2)},),
        ((-3.0, 3.0), (-3.0, 3.0), (-0.5 * math.pi, 0.5 * math.pi)),
    )
    polar_values = torch.tensor([[1.2, 0.7, -0.4], [2.1, -1.1, 0.6]], dtype=torch.float64)
    cartesian_values = polar_values.clone()
    cartesian_values[:, 0] = polar_values[:, 0] * torch.cos(polar_values[:, 1])
    cartesian_values[:, 1] = polar_values[:, 0] * torch.sin(polar_values[:, 1])
    scale = torch.ones(1, dtype=torch.float64)

    polar_distance = polar.squared_distance(polar_values[:1], polar_values[1:], scale)
    cartesian_distance = cartesian.squared_distance(cartesian_values[:1], cartesian_values[1:], scale)
    assert cartesian_distance.item() == pytest.approx(polar_distance.item(), abs=1.0e-12)


# Verify spherical antipodes remain distinct in the direct HEBO kernel
def test_spherical_kernel_antipodes():
    model = icetune_hebo.TopologyGP(
        2,
        0,
        1,
        topology_groups=[{"kind": "sphere", "indices": [0, 1]}],
        numeric_bounds=[
            [-0.5 * math.pi, 0.5 * math.pi],
            [-math.pi, math.pi],
        ],
        num_epochs=1,
    )
    values = torch.tensor([[0.0, 0.0], [0.0, math.pi]], dtype=torch.float64)

    covariance = model.model.conf["kern"](values).to_dense()
    correlation = covariance[0, 1] / torch.sqrt(covariance[0, 0] * covariance[1, 1])

    assert correlation.item() < 0.2


# Verify a GP-sized kernel is positive semidefinite and does not collapse
def test_topology_kernel_stable_gp_dimension():
    angle_counts = (2, 2, 6, 2, 6, 2, 6)
    groups = []
    bounds = []
    offset = 0
    for count in angle_counts:
        groups.append(
            {
                "kind": "polar_projective",
                "indices": tuple(range(offset, offset + count + 2)),
                "period": 2.0 * math.pi,
            }
        )
        bounds.extend([(0.0, 3.0), (-math.pi, math.pi), *((-0.5 * math.pi, 0.5 * math.pi),) * count])
        offset += count + 2
    groups.extend(
        {"kind": "phase", "indices": (offset + index,), "period": 2.0 * math.pi} for index in range(8)
    )
    bounds.extend([(-math.pi, math.pi)] * 8)
    kernel = icetune_hebo.TopologyMaternKernel(48, tuple(groups), tuple(bounds)).double()
    generator = torch.Generator().manual_seed(1234)
    values = torch.empty((80, 48), dtype=torch.float64)
    for index, (lower, upper) in enumerate(bounds):
        values[:, index] = lower + (upper - lower) * torch.rand(80, generator=generator, dtype=torch.float64)

    covariance = kernel(values).to_dense()
    off_diagonal = covariance[~torch.eye(len(values), dtype=torch.bool)]

    assert torch.linalg.eigvalsh(covariance).min().item() >= -1.0e-10
    assert 0.1 < torch.median(off_diagonal).item() < 0.9
    assert kernel.raw_group_lengthscale.numel() == 15


# Verify topology metadata selects the local HEBO surrogate adapter
def test_create_hebo_selects_phase_topology_model():
    optimizer = icetune_hebo.create_hebo(
        _design(),
        parameter_topology=_topology(),
        rand_sample=0,
        scramble_seed=17,
    )

    assert optimizer.model_name == icetune_hebo.TOPOLOGY_MODEL_NAME
    assert optimizer.model_config["topology_groups"] == [
        {
            "kind": "phase",
            "base": "amplitude",
            "indices": [1],
            "period": 2.0 * math.pi,
        },
    ]
    assert optimizer.model_config["numeric_bounds"] == [
        [2.0, 6.0],
        [-math.pi, math.pi],
    ]
    assert optimizer.model_config["max_exact_trials"] == load_hebo_config()["gp"]["max_exact_trials"]
    assert optimizer.model_config["sparse_config"]["num_inducing"] == 128


# Verify the standard hybrid retains HEBO's exact GP settings and factory API
def test_create_hebo_selects_standard_hybrid_model():
    from hebo.models.model_factory import get_model

    optimizer = icetune_hebo.create_hebo(
        _design(),
        parameter_topology=None,
        rand_sample=0,
        scramble_seed=17,
    )
    model = get_model(
        optimizer.model_name,
        optimizer.space.num_numeric,
        optimizer.space.num_categorical,
        1,
        **optimizer.model_config,
    )

    assert optimizer.model_name == icetune_hebo.HYBRID_MODEL_NAME
    assert optimizer.model_config["lr"] == pytest.approx(0.01)
    assert optimizer.model_config["noise_lb"] == pytest.approx(8.0e-4)
    assert optimizer.model_config["pred_likeli"] is False
    assert optimizer.model_config["max_exact_trials"] == load_hebo_config()["gp"]["max_exact_trials"]
    assert optimizer.model_config["sparse_config"]["num_inducing"] == 128
    assert isinstance(model, icetune_hebo.HybridGP)


# Check acquisition reproducibility despite unrelated Python, NumPy and Torch draws
def test_seeded_hebo_scopes_random_sources():
    import pandas as pd

    optimizers = [
        icetune_hebo.create_hebo(
            _design(),
            parameter_topology=None,
            rand_sample=0,
            scramble_seed=17,
        )
        for _ in range(2)
    ]

    for optimizer in optimizers:
        optimizer._model_config = {"num_epochs": 1, "verbose": False}
        optimizer.observe(pd.DataFrame({"scale": [2.5, 4.0, 5.5], "phase": [-1.0, 0.0, 1.0]}),
                          np.array([[3.25], [0.0], [3.25]]))
    random.seed(1)
    np.random.seed(2)
    torch.manual_seed(3)
    first = optimizers[0].suggest(3)
    random.seed(101)
    np.random.seed(102)
    torch.manual_seed(103)
    second = optimizers[1].suggest(3)

    assert first.equals(second)


# Verify projective metadata alone selects the topology-aware HEBO model
def test_hebo_projective_model():
    from hebo.design_space.design_space import DesignSpace

    angle_names = ["alpha_angle0@PROJECTIVE", "alpha_angle1@PROJECTIVE"]
    design = DesignSpace().parse(
        [
            {
                "name": name,
                "type": "num",
                "lb": -0.5 * math.pi,
                "ub": 0.5 * math.pi,
            }
            for name in angle_names
        ]
    )
    topology = {
        "schema_version": 1,
        "groups": [
            {"kind": "projective", "base": "alpha", "parameters": angle_names},
        ],
    }

    optimizer = icetune_hebo.create_hebo(
        design,
        parameter_topology=topology,
        rand_sample=0,
        scramble_seed=17,
    )

    assert optimizer.model_name == icetune_hebo.TOPOLOGY_MODEL_NAME
    assert optimizer.model_config["topology_groups"] == [
        {"kind": "projective", "base": "alpha", "indices": [0, 1]},
    ]


# Verify spherical metadata selects the topology aware HEBO model
def test_hebo_spherical_model():
    from hebo.design_space.design_space import DesignSpace

    angle_names = ["alpha_angle0@SPHERICAL", "alpha_angle1@SPHERICAL"]
    design = DesignSpace().parse(
        [
            {
                "name": angle_names[0],
                "type": "num",
                "lb": -0.5 * math.pi,
                "ub": 0.5 * math.pi,
            },
            {
                "name": angle_names[1],
                "type": "num",
                "lb": -math.pi,
                "ub": math.pi,
            },
        ]
    )
    topology = {
        "schema_version": 1,
        "groups": [{"kind": "sphere", "base": "alpha", "parameters": angle_names}],
    }

    optimizer = icetune_hebo.create_hebo(
        design,
        parameter_topology=topology,
        rand_sample=0,
        scramble_seed=17,
    )

    assert optimizer.model_name == icetune_hebo.TOPOLOGY_MODEL_NAME
    assert optimizer.model_config["topology_groups"] == [
        {"kind": "sphere", "base": "alpha", "indices": [0, 1]},
    ]


# Verify the complete icetune HEBO adapter lowers a convex toy objective
def test_search_state_hebo_improves_quadratic_toy(monkeypatch):
    observed = (-1.0, -0.5, 0.5, 1.0)
    records = [
        {
            "trial_id": f"trial-{index:06d}",
            "config": {"x": value},
            "metrics": {"loss": value**2},
        }
        for index, value in enumerate(observed)
    ]
    state = SearchState(
        args=_search_args(len(observed)),
        bounds={"x": {"type": "uniform", "lower": -1.0, "upper": 1.0}},
        initial_points=None,
        async_proposals=True,
    )
    state.hebo._model_config = {"num_epochs": 1, "verbose": False}
    state.observe(records)
    fitted = []
    suggest = state.hebo.suggest

    # Record the observations present at every adaptive fit
    def suggest_with_history(**kwargs):
        fitted.append(len(state.hebo.y))
        return suggest(**kwargs)

    monkeypatch.setattr(state.hebo, "suggest", suggest_with_history)

    proposals = []
    for offset in range(3):
        proposal = state.ask_many(len(records), 1)[0]
        proposals.append(proposal)
        records.append(
            {
                "trial_id": f"trial-{len(records):06d}",
                "config": proposal[0],
                "metrics": {"loss": float(proposal[0]["x"]) ** 2},
            }
        )
        if offset < 2:
            state.observe(records)
    adaptive_best = min(float(config["x"]) ** 2 for config, _ in proposals)

    assert state.hebo.model_name == icetune_hebo.HYBRID_MODEL_NAME
    assert len(proposals) == 3
    assert fitted == [4, 5, 6]
    assert all(payload["kind"] == "hebo" for _, payload in proposals)
    assert adaptive_best < 0.2 * min(value**2 for value in observed)


# Verify random proposals stop exactly at the configured trial boundary
def test_search_state_hebo_random_boundary_exact():
    state = SearchState(
        args=_search_args(5),
        bounds={"x": {"type": "uniform", "lower": -1.0, "upper": 1.0}},
        initial_points=None,
        async_proposals=True,
    )
    proposals = state.ask_many(0, 8)

    assert [payload["kind"] for _, payload in proposals] == ["random"] * 5 + ["hebo"] * 3
    assert [payload["seed"] for _, payload in proposals[:5]] == [17, 18, 19, 20, 21]
    assert all(-1.0 <= config["x"] <= 1.0 for config, _ in proposals)


# Verify failed configurations are never returned or fitted as successful trials
def test_search_state_hebo_excludes_failed_config():
    state = SearchState(
        args=_search_args(0),
        bounds={"x": {"type": "uniform", "lower": -1.0, "upper": 1.0}},
        initial_points=None,
        async_proposals=True,
    )
    state.observe(
        [{"trial_id": "trial-000000", "config": {"x": 0.0}, "metrics": {"loss": 1.0}}]
    )
    state.observe_failures(
        [{"trial_id": "trial-000001", "config": {"x": 0.5}, "search_payload": {}}]
    )
    state.hebo._model_config = {"num_epochs": 1, "verbose": False}

    proposal = state.ask_many(2, 1)

    assert not np.isclose(proposal[0][0]["x"], 0.5, rtol=0.0, atol=1e-12)
    assert not np.isclose(proposal[0][0]["x"], 0.0, rtol=0.0, atol=1e-12)
    assert state.hebo.X.to_dict(orient="records") == [{"x": 0.0}]
    assert state.hebo.y.tolist() == [[1.0]]


# Verify many live trials never inflate one HEBO acquisition request
def test_search_state_hebo_pending_request_bounded(monkeypatch):
    state = SearchState(
        args=_search_args(0),
        bounds={"x": {"type": "uniform", "lower": -1.0, "upper": 1.0}},
        initial_points=None,
        async_proposals=True,
    )
    state.observe(
        [{"trial_id": "trial-000000", "config": {"x": 0.0}, "metrics": {"loss": 1.0}}]
    )
    state.set_pending([{"x": -0.9 + 0.02 * index} for index in range(50)])
    requests = []

    state.hebo._model_config = {"num_epochs": 1, "verbose": False}
    original = state.hebo.suggest

    # Record request sizes while running the real GP and acquisition optimizer
    def suggest(**kwargs):
        requests.append(kwargs["n_suggestions"])
        return original(**kwargs)

    monkeypatch.setattr(state.hebo, "suggest", suggest)

    proposal = state.ask_many(51, 1)
    assert len(proposal) == 1
    assert -1.0 <= proposal[0][0]["x"] <= 1.0
    assert requests == [3]


# Verify canonical trial ordering removes completion-order dependence
def test_hebo_observation_order():
    records = [
        {"trial_id": "trial-000002", "config": {"x": 0.6}, "metrics": {"loss": 0.2}},
        {"trial_id": "trial-000000", "config": {"x": -0.6}, "metrics": {"loss": 0.4}},
    ]
    first = SearchState(
        args=_search_args(0),
        bounds={"x": {"type": "uniform", "lower": -1.0, "upper": 1.0}},
        initial_points=None,
        async_proposals=True,
    )
    second = SearchState(
        args=_search_args(0),
        bounds={"x": {"type": "uniform", "lower": -1.0, "upper": 1.0}},
        initial_points=None,
        async_proposals=True,
    )

    first.observe(records)
    second.observe(list(reversed(records)))

    assert first.hebo.X.equals(second.hebo.X)
    assert first.hebo.y.tolist() == second.hebo.y.tolist()


# Verify the real Pymoo acquisition is deterministic after a HEAD restart
def test_search_state_real_hebo_rebuild_seeded():
    records = [
        {"trial_id": "trial-000000", "config": {"x": -0.5}, "metrics": {"loss": 0.4}},
        {"trial_id": "trial-000001", "config": {"x": 0.5}, "metrics": {"loss": 0.2}},
    ]
    states = [
        SearchState(
            args=_search_args(0),
            bounds={"x": {"type": "uniform", "lower": -1.0, "upper": 1.0}},
            initial_points=None,
            async_proposals=True,
        )
        for _ in range(2)
    ]
    proposals = []
    for state in states:
        state.hebo._model_config = {"num_epochs": 1, "verbose": False}
        state.observe(records)
        proposals.append(state.ask_many(2, 1))

    assert proposals[0] == proposals[1]


# Verify Optuna pending replay does not duplicate native running trials
def test_search_state_optuna_pending_not_duplicated():
    optuna = pytest.importorskip("optuna")
    args = _search_args(0)
    args.algorithm = "optuna"
    state = SearchState(
        args=args,
        bounds={"x": {"type": "uniform", "lower": -1.0, "upper": 1.0}},
        initial_points=None,
        async_proposals=True,
    )
    records = [
        {"trial_id": "trial-000000", "config": {"x": 0.0}, "metrics": {"loss": 1.0}}
    ]
    state.observe(records)
    proposals = state.ask_many(1, 2)
    state.set_pending([config for config, _ in proposals])

    running = [
        trial
        for trial in state.optuna_study.trials
        if trial.state == optuna.trial.TrialState.RUNNING
    ]
    assert len(running) == 2

    state.observe(
        [
            *records,
            *(
                {
                    "trial_id": f"trial-{index:06d}",
                    "config": config,
                    "metrics": {"loss": float(config["x"]) ** 2},
                    "search_payload": payload,
                }
                for index, (config, payload) in enumerate(proposals, start=1)
            ),
        ]
    )
    assert not any(
        trial.state == optuna.trial.TrialState.RUNNING
        for trial in state.optuna_study.trials
    )


# Verify the search-state default fits every successful trial
def test_search_state_hebo_all_successes_default():
    args = _search_args(0)
    del args.surrogate_fit_percentile
    records = [
        {
            "trial_id": f"trial-{index:06d}",
            "config": {"x": value},
            "metrics": {"loss": float(index)},
        }
        for index, value in enumerate((-1.0, -0.5, 0.0, 0.5, 1.0))
    ]
    state = SearchState(
        args=args,
        bounds={"x": {"type": "uniform", "lower": -1.0, "upper": 1.0}},
        initial_points=None,
        async_proposals=True,
    )

    state.observe(records)

    assert state.surrogate_fit_percentile == 1.0
    assert state.hebo.X.to_dict(orient="records") == [record["config"] for record in records]
    assert state.hebo.y.reshape(-1).tolist() == [0.0, 1.0, 2.0, 3.0, 4.0]


# Verify topology-aware HEBO improves a nontrivial circular two-dimensional toy
def test_search_state_hebo_improves_circular_toy():
    observed = [(2.0, -2.0), (6.0, -2.0), (2.0, 2.0), (6.0, 2.0), (4.0, 0.0)]
    records = [
        {
            "trial_id": f"trial-{index:06d}",
            "config": {"scale": scale, "phase": phase},
            "metrics": {"loss": _circular_loss(scale, phase)},
        }
        for index, (scale, phase) in enumerate(observed)
    ]
    state = SearchState(
        args=_search_args(len(observed), 1.0),
        bounds={
            "scale": {"type": "uniform", "lower": 2.0, "upper": 6.0},
            "phase": {"type": "uniform", "lower": -math.pi, "upper": math.pi},
        },
        initial_points=None,
        parameter_topology=_topology(),
        async_proposals=True,
    )
    state.hebo._model_config["num_epochs"] = 1
    state.observe(records)

    proposals = state.ask_many(len(observed), 3)
    warm_best = min(record["metrics"]["loss"] for record in records)
    adaptive_best = min(
        _circular_loss(float(config["scale"]), float(config["phase"])) for config, _ in proposals
    )

    assert state.hebo.model_name == icetune_hebo.TOPOLOGY_MODEL_NAME
    assert all(payload["kind"] == "hebo" for _, payload in proposals)
    assert adaptive_best < 0.2 * warm_best


# Verify complex coupling covariance agrees with its physical vector for unequal batches
@pytest.mark.parametrize("coordinate_type", ["polar_projective", "cartesian_projective", "phase_projective"])
def test_complex_projective_rectangular_cov(coordinate_type):
    dimension = 3 if coordinate_type == "phase_projective" else 4
    group = {"kind": coordinate_type, "indices": tuple(range(dimension)), "weights": (2.0, 1.0, 2.0)}
    if coordinate_type != "cartesian_projective":
        group["period"] = 2.0 * math.pi
    geometry = ProductManifoldGeometry(dimension, (group,), ((-3.0, 3.0),) * dimension).double()
    values = torch.linspace(-1.1, 1.3, 11 * dimension, dtype=torch.float64).reshape(11, dimension)
    coefficient = (
        torch.exp(1j * values[:, 0])
        if coordinate_type == "phase_projective"
        else (
            torch.complex(values[:, 0], values[:, 1]) / 3.0
            if coordinate_type == "cartesian_projective"
            else values[:, 0] * torch.exp(1j * values[:, 1]) / 3.0
        )
    )
    a, b = values[:, -2], values[:, -1]
    direction = torch.stack((torch.cos(a) * torch.cos(b), torch.sin(a), torch.cos(a) * torch.sin(b)), dim=1)
    direction *= torch.sqrt(torch.tensor(group["weights"], dtype=torch.float64))
    direction /= torch.linalg.vector_norm(direction, dim=1, keepdim=True)
    vectors = coefficient[:, None] * direction
    expected = (vectors[:7, None] - vectors[None, 7:]).abs().square().sum(dim=-1)
    first = values[:7].clone().requires_grad_()
    distance = geometry.squared_distance(first, values[7:], torch.ones(1, dtype=torch.float64))
    torch.testing.assert_close(distance, expected)
    torch.testing.assert_close(geometry.squared_distance(values[7:], first, torch.ones(1)), distance.T)
    distance.sum().backward()
    assert torch.isfinite(first.grad).all()


# Verify both GP solvers predict complex coupling dependence across unequal batch sizes
@pytest.mark.parametrize("sparse", [False, True])
def test_complex_gp_residue_variation(sparse):
    model = icetune_hebo.TopologyGP(
        3,
        0,
        1,
        topology_groups=[{"kind": "polar_projective", "indices": [0, 1, 2], "period": 2.0 * math.pi}],
        numeric_bounds=[(0.0, 3.0), (-math.pi, math.pi), (-math.pi / 2.0, math.pi / 2.0)],
        optimizer="adam",
        num_epochs=15,
        max_exact_trials=0 if sparse else 100,
        sparse_config={"num_epochs": 30, "num_inducing": 8, "batch_size": 32},
    )
    values = torch.zeros(12, 3)
    values[:, 0] = torch.linspace(0.1, 2.9, 12)
    values[:, 1] = 0.2
    values[:, 2] = 0.3
    model.fit(values, torch.zeros((12, 0), dtype=torch.long), values[:, :1].square())
    query = values[[0, 3, 7, 11]]
    mean, variance = model.predict(query, torch.zeros((4, 0), dtype=torch.long))
    assert torch.isfinite(mean).all() and torch.isfinite(variance).all()
    assert (variance > 0.0).all()
    assert mean[-1, 0] - mean[0, 0] > 1.0


# Verify short acquisition batches keep the best predicted point and remain distinct
@pytest.mark.parametrize("count", [1, 2, 3])
def test_hebo_short_batch_best_prediction(count):
    import pandas as pd

    # Predict a known minimum with the largest variance at a different point
    def predict(continuous, categorical):
        return continuous[:, :1].square(), (continuous[:, :1] - 1.0).square() + 0.1

    recommendations = pd.DataFrame({"scale": [2.1, 3.0, 5.9], "phase": [0.0, 0.1, 0.2]})
    selected = icetune_hebo._select_recommendations(SimpleNamespace(predict=predict), _design(), recommendations, count)
    assert selected.iloc[0]["scale"] == pytest.approx(2.1)
    assert len(selected) == count
    assert not selected.duplicated().any()
    if count > 1:
        assert selected.iloc[1]["scale"] == pytest.approx(5.9)


# Reject malformed settings at input parsing instead of failing during fitting
@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("gp", "lr", 0.0),
        ("gp", "num_epochs", True),
        ("gp", "optimizer", "typo"),
        ("gp", "noise_lb", 0.1),
        ("sparse_gp", "num_inducing", 0),
        ("topology", "nu", 1.0),
        ("topology", "lengthscale_scale", math.nan),
        ("topology", "enabled", 1),
        ("acquisition", "population", 2),
        ("acquisition", "delta", 1.0),
        ("pending", "lie", "unknown"),
        ("pending", "penalty_scale", -1.0),
        ("output", "transform", "unknown"),
        ("gp", "unknown", 1),
    ],
)
def test_hebo_config_invalid_settings(section, key, value):
    settings = load_hebo_config()
    settings[section][key] = value
    with pytest.raises(ValueError, match="HEBO"):
        validate_hebo_config(settings)


# Verify a selected JSON reaches CLI parsing, history and the real surrogate factory
def test_hebo_config_reaches_cli_search_history(monkeypatch, tmp_path):
    from core.tune.core import build_optimization_metadata
    from hebo.models.model_factory import get_model

    settings = load_hebo_config()
    settings["gp"].update(lr=0.025, num_epochs=7, optimizer="adam", max_exact_trials=9, ard_kernel=False)
    settings["sparse_gp"].update(num_inducing=6, num_epochs=11)
    settings["topology"].update(nu=2.5, lengthscale_scale=0.4)
    settings["acquisition"].update(population=24, generations=12)
    settings["pending"].update(lie="mean", penalty_scale=0.0)
    settings["output"]["transform"] = "standardize"
    source = tmp_path / "hebo.json"
    source.write_text(json.dumps(settings), encoding="utf-8")
    module = runpy.run_module("core.icetune")
    monkeypatch.setattr(
        sys, "argv", ["icetune", "--algorithm", "hebo", "--cdir", str(tmp_path), "--hebo_config", "hebo.json"]
    )
    args = module["parse_arguments"]()
    bounds = {
        name: {"type": "uniform", "lower": lower, "upper": upper}
        for name, lower, upper in [("scale", 2.0, 6.0), ("phase", -math.pi, math.pi)]
    }
    state = SearchState(args=args, bounds=bounds, initial_points=None, parameter_topology=_topology())
    assert state.hebo.settings == settings
    assert build_optimization_metadata(args)["hebo"] == settings
    model = get_model(state.hebo.model_name, 2, 0, 1, **state.hebo.model_config)
    assert model.model.num_epochs == 7
    assert model.model.lr == pytest.approx(0.025)
    assert model.sparse_config["num_inducing"] == 6
    kernel = model.model.conf["kern"].base_kernel
    assert kernel.matern52
    assert kernel.raw_group_lengthscale.numel() == 1
    torch.testing.assert_close(kernel.group_lengthscale, torch.full((1, 2), 0.4 * math.sqrt(2)))
    settings["topology"]["enabled"] = False
    plain = icetune_hebo.create_hebo(
        _design(), parameter_topology=_topology(), rand_sample=0, scramble_seed=1, settings=settings
    )
    assert plain.model_name == icetune_hebo.HYBRID_MODEL_NAME
    assert state.hebo.settings["topology"]["enabled"]


# Verify adaptive acquisition consumes the configured population and generation budget
def test_hebo_acquisition_uses_json_settings(monkeypatch):
    import pandas as pd
    from hebo.acquisitions.acq import MACE

    settings = load_hebo_config()
    settings["gp"].update(num_epochs=1, optimizer="adam")
    settings["acquisition"].update(population=12, generations=3)
    optimizer = icetune_hebo.create_hebo(
        _design(), parameter_topology=_topology(), rand_sample=0, scramble_seed=17, settings=settings
    )
    optimizer.observe(
        pd.DataFrame({"scale": [2.2, 3.2, 5.5], "phase": [0.1, 0.5, -0.1]}), np.array([[3.0], [1.0], [4.0]])
    )
    calls = []
    original = icetune_hebo.SeededEvolutionOpt.optimize

    # Inspect the actual acquisition optimizer and run its evolution
    def optimize(evolution, *args, **kwargs):
        calls.append((evolution.pop, evolution.iterations))
        assert isinstance(evolution.acquisition, MACE)
        return original(evolution, *args, **kwargs)

    monkeypatch.setattr(icetune_hebo.SeededEvolutionOpt, "optimize", optimize)
    assert len(optimizer.suggest(3)) == 3
    assert calls == [(12, 3)]


# Verify real HEBO proposals improve a complex coupling objective through projective topology
def test_hebo_improves_complex_residue_objective():
    settings = load_hebo_config()
    settings["gp"].update(num_epochs=35, optimizer="adam")
    settings["acquisition"].update(population=24, generations=20)
    args = _search_args(12)
    args.rngseed = 23
    args.hebo_settings = settings
    bounds = {
        name: {"type": "uniform", "lower": lower, "upper": upper}
        for name, lower, upper in [("r", 0.0, 2.0), ("phi", -math.pi, math.pi), ("theta", -math.pi / 2, math.pi / 2)]
    }
    topology = {
        "schema_version": 1,
        "groups": [{"kind": "polar_projective", "base": "residue", "parameters": ["r", "phi", "theta"], "period": 2 * math.pi}],
    }
    state = SearchState(args=args, bounds=bounds, initial_points=None, parameter_topology=topology, async_proposals=True)
    target = 0.7 * np.exp(2.8j) * np.array([math.cos(0.4), math.sin(0.4)])
    records = []
    for index in range(18):
        config, payload = state.ask_many(index, 1)[0]
        vector = config["r"] * np.exp(1j * config["phi"]) * np.array([math.cos(config["theta"]), math.sin(config["theta"])])
        loss = float(np.sum(np.abs(vector - target) ** 2))
        records.append({"trial_id": f"trial-{index:06d}", "config": config, "metrics": {"loss": loss}, "search_payload": payload})
        state.observe(records)
    initial = min(record["metrics"]["loss"] for record in records[:12])
    best = min(record["metrics"]["loss"] for record in records[12:])
    assert best < 0.5 * initial
    assert all(record["search_payload"]["kind"] == "hebo" for record in records[12:])


# Check a single continuum helicity norm reaches the adaptive HEBO model as a scalar
def test_single_helicity_norm_hebo_fit():
    from hebo.design_space.design_space import DesignSpace

    name = "CON_GP|990:[211,211]/opposite:helicity(0,0,0)@NORM"
    names = [name]
    metadata = tools.build_parameter_topology(names)
    design = DesignSpace().parse([{"name": name, "type": "num", "lb": 0.0, "ub": 2.0}])
    optimizer = icetune_hebo.create_hebo(design, parameter_topology=metadata, rand_sample=0, scramble_seed=17)
    points = design.sample(4)
    optimizer.observe(points, (points[name].to_numpy() ** 2).reshape(-1, 1))
    optimizer._model_config.update(num_epochs=1, verbose=False)
    proposal = optimizer.suggest(n_suggestions=1)
    assert np.isfinite(proposal[name]).all()
    assert proposal[name].between(0.0, 2.0).all()
    assert not topology.indexed_parameter_topology(names, metadata)

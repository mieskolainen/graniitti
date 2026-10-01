# Tests for icescape parameter topology transformations
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
import sys
from pathlib import Path

import numpy as np
import pytest
import torch

ROOT = Path(__file__).resolve().parents[3]

from core.tune.parameters import space as parameter_space
from core.tune.parameters import tools


# Load the extensionless icescape command as a testable Python module
def _load_icescape():
    from importlib import import_module

    return import_module('core.icescape')


icescape = _load_icescape()


# Check generic scalar comparison accepts NumPy scalars without coercing strings
def test_scalar_values_equal_handles_numpy_scalars():
    assert parameter_space.scalar_values_equal(np.float64(0.5), 0.5)
    assert parameter_space.scalar_values_equal(np.int64(2), 2)
    assert not parameter_space.scalar_values_equal("2", 2)


# Select exactly one likelihood surrogate through the model option
def test_icescape_model_cli(monkeypatch, tmp_path):
    history_path = tmp_path / "runs" / "icetune" / "mesh_cli" / "history.json"
    history_path.parent.mkdir(parents=True)
    history_path.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "icescape",
            "--input",
            str(history_path),
            "--model",
            "neural",
            "--cached",
            "--mc-errors",
            "--posterior",
            "--profile",
            "true",
            "--plot-2d",
            "--profile-multistart",
            "--profile_workers",
            "6",
            "--profile_dpi",
            "320",
            "--impact-mode",
            "fixed",
            "--impact-profile-maxiter",
            "12",
            "--impact-top",
            "7",
            "--plot-brand",
            "CUSTOM",
            "--no-gp-optimize",
        ],
    )

    args = icescape.parse_args()

    assert args.model == "neural"
    assert args.cached and args.mc_errors and args.posterior
    assert args.profile and args.plot_2d and args.profile_multistart
    assert args.profile_dpi == 320
    assert args.profile_workers == 6
    assert not hasattr(args, "impacts")
    assert args.impact_top == 7
    assert args.impact_mode == "fixed" and args.impact_profile_maxiter == 12
    assert args.plot_brand == "CUSTOM"
    assert args.gp_dtype == "float64" and args.gp_validation_interval == 5
    assert args.optimizer_backend == "torch"
    assert args.optimizer_restarts == 32 and args.optimizer_scan_points == 1024
    assert args.gp_kernel == "Matern" and args.gp_target_transform == "auto"
    assert args.gp_restarts == 1
    assert not args.gp_optimize
    assert args.input_mode == "history"
    assert args.history_path == str(history_path)
    assert args.run_name == "mesh_cli"


# Check profile selection requires an explicit Boolean value
def test_icescape_profile_cli_is_explicit(monkeypatch, tmp_path):
    history_path = tmp_path / "runs" / "icetune" / "mesh_default" / "history.json"
    history_path.parent.mkdir(parents=True)
    history_path.write_text("{}", encoding="utf-8")

    monkeypatch.setattr(
        sys,
        "argv",
        ["icescape", "--input", str(history_path), "--profile", "false"],
    )
    disabled = icescape.parse_args()
    assert disabled.profile is False

    monkeypatch.setattr(
        sys,
        "argv",
        ["icescape", "--input", str(history_path), "--profile"],
    )
    with pytest.raises(SystemExit):
        icescape.parse_args()

    monkeypatch.setattr(
        sys,
        "argv",
        ["icescape", "--input", str(history_path), "--no-profile"],
    )
    with pytest.raises(SystemExit):
        icescape.parse_args()


# Check validation splitting always retains explicitly required likelihood rows
def test_training_validation_split_realized_minimum():
    training, validation = icescape.training_validation_indices(
        n_rows=20,
        validation_fraction=0.25,
        rngseed=17,
        required_training=[7],
    )

    assert 7 in training
    assert 7 not in validation
    assert len(validation) == 5
    assert len(np.intersect1d(training, validation)) == 0


# Check automatic GP target transforms preserve positive objective values
def test_gp_log_target_transform_round_trip():
    values = np.array([2.0, 5.0, 20.0], dtype=np.float64)
    transform = icescape.fit_gp_target_transform(
        values=values,
        training_indices=np.array([0, 1, 2]),
        mode="auto",
    )
    scaled, scaled_errors = icescape.transform_gp_targets(
        values=values,
        errors=np.array([0.2, 0.5, 2.0]),
        transform=transform,
    )
    restored, _ = icescape.inverse_gp_prediction(
        mean_scaled=torch.as_tensor(scaled),
        variance_scaled=None,
        transform=transform,
    )

    assert transform["kind"] == "log"
    assert transform["uncertainty"] == "delta_log_mean"
    assert np.all(scaled_errors > 0.0)
    assert np.allclose(restored.numpy().reshape(-1), values)


# Check the angular chart covers projective directions with a zero first row
def test_projective_zero_roundtrip():
    vectors = [
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 0.2, -0.3, 0.4, -0.5, 0.6, -0.7],
    ]
    for vector in vectors:
        canonical = tools.canonicalize_projective_vector(vector)
        angles = tools.projective_angles_from_vector(vector)
        decoded = tools.projective_vector_from_angles(angles)
        assert np.allclose(decoded, canonical, atol=1.0e-12)
        assert all(abs(angle) <= 0.5 * math.pi for angle in angles)


# Check the f2 chart reaches a pure non-reference LS component
def test_projective_pure_f2():
    target = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    angles = tools.projective_angles_from_vector(target)
    assert angles[-1] == pytest.approx(0.5 * math.pi)
    assert np.allclose(tools.projective_vector_from_angles(angles), target, atol=1.0e-12)


# Check unmarked Pandora-style parameters retain the exact generic unit-cube path
def test_pandora_generic_features():
    names = ["PANDORA|PhotonRecovery:MaxDistance", "PANDORA|TrackCluster:Chi2Cut"]
    bounds = np.array([[0.0, 10.0], [-2.0, 2.0]], dtype=np.float64)
    values = np.array([[2.5, -1.0], [8.0, 1.5]], dtype=np.float64)
    topology = tools.build_parameter_topology(names)

    assert topology == {}
    assert np.array_equal(
        icescape.model_features(values, bounds, names, topology),
        parameter_space.to_unit(values, bounds),
    )
    assert icescape.model_feature_names(names, topology) == names


# Check amplitude phases and constrained vectors are expanded into physical features
def test_amp_topology_groups():
    phase = tools.raw_phase_key("DECAY|225:[321,-321]:zeta.GP")
    projective = [tools.projective_angle_key("RES|f2:GP:g_ls", index) for index in range(2)]
    simplex = [tools.simplex_theta_key("RES|f2:MP:polarization.rho_population", 0)]
    names = [phase, *projective, *simplex, "REGGE|omega.GP"]
    bounds = np.array(
        [
            [-math.pi, math.pi],
            [-0.5 * math.pi, 0.5 * math.pi],
            [-0.5 * math.pi, 0.5 * math.pi],
            [0.0, 0.5 * math.pi],
            [0.0, 1.0],
        ],
        dtype=np.float64,
    )
    values = np.array([[0.5, 0.7, -0.4, 0.6, 0.25]], dtype=np.float64)
    topology = tools.build_parameter_topology(names)
    features = icescape.model_features(values, bounds, names, topology)

    expected = np.array(
        [
            [
                math.cos(0.5),
                math.sin(0.5),
                *tools.projective_embedding_from_angles([0.7, -0.4]),
                *tools.simplex_probabilities_from_angles([0.6]),
                0.25,
            ]
        ]
    )
    assert np.allclose(features, expected)
    assert features.shape == (1, 11)


# Check Torch topology features match NumPy and remain differentiable
def test_topology_autograd():
    phase = tools.raw_phase_key("DECAY|225:[321,-321]:zeta.GP")
    projective = [tools.projective_angle_key("RES|f2:GP:g_ls", index) for index in range(2)]
    simplex = [tools.simplex_theta_key("RES|f2:MP:polarization.rho_population", 0)]
    names = [phase, *projective, *simplex, "REGGE|omega.GP"]
    bounds = np.array(
        [
            [-math.pi, math.pi],
            [-0.5 * math.pi, 0.5 * math.pi],
            [-0.5 * math.pi, 0.5 * math.pi],
            [0.0, 0.5 * math.pi],
            [0.0, 1.0],
        ],
        dtype=np.float64,
    )
    values = np.array([0.5, 0.7, -0.4, 0.6, 0.25], dtype=np.float64)
    topology = tools.build_parameter_topology(names)
    torch_values = torch.tensor(values, dtype=torch.float64, requires_grad=True)

    torch_features = icescape.model_features_torch(
        torch_values,
        bounds,
        names,
        topology,
    )
    numpy_features = icescape.model_features(
        values,
        bounds,
        names,
        topology,
    )
    weights = torch.arange(
        1,
        1 + len(torch_features),
        dtype=torch.float64,
    )
    (gradient,) = torch.autograd.grad(torch.sum(weights * torch_features), torch_values)

    assert np.allclose(torch_features.detach().numpy(), numpy_features, atol=1.0e-12)
    assert torch_features.shape == (11,)
    assert torch.all(torch.isfinite(gradient))
    assert torch.count_nonzero(gradient).item() == len(values)


# Check derivative callbacks return exact physical gradients and Hessians
def test_torch_quadratic_derivatives():
    target = torch.tensor([0.2, -0.4, 0.7], dtype=torch.float64)
    scales = torch.tensor([1.0, 2.0, 3.0], dtype=torch.float64)

    # Compute one anisotropic quadratic objective
    def objective(values):
        return torch.sum(scales * (values - target).square())

    value_gradient, hessian = icescape.build_torch_derivative_callbacks(
        objective=objective,
        device="cpu",
        dtype=torch.float64,
    )
    point = np.array([0.5, -0.1, 0.1], dtype=np.float64)
    value, gradient = value_gradient(point)
    matrix = hessian(point)

    expected_gradient = 2.0 * scales.numpy() * (point - target.numpy())
    assert value == pytest.approx(np.sum(scales.numpy() * (point - target.numpy()) ** 2))
    assert np.allclose(gradient, expected_gradient, atol=1.0e-12)
    assert np.allclose(matrix, np.diag(2.0 * scales.numpy()), atol=1.0e-12)


# Check icescape identifies antipodal projective chart boundaries
def test_projective_feature_antipodes():
    projective = [tools.projective_angle_key("RES|f0:GP:g_ls", index) for index in range(2)]
    names = [*projective, "REGGE|omega.GP"]
    bounds = np.array(
        [
            [-0.5 * math.pi, 0.5 * math.pi],
            [-0.5 * math.pi, 0.5 * math.pi],
            [0.0, 1.0],
        ],
        dtype=np.float64,
    )
    angle = 0.37
    values = np.array(
        [
            [angle, 0.5 * math.pi, 0.25],
            [-angle, -0.5 * math.pi, 0.25],
        ],
        dtype=np.float64,
    )
    topology = tools.build_parameter_topology(names)

    features = icescape.model_features(values, bounds, names, topology)

    assert features.shape == (2, 7)
    assert np.allclose(features[0], features[1], atol=1.0e-12)


# Check resonance norms, phases and signed-real directions form one complex vector
def test_complex_projective_features():
    base = "RES|f0:GP:g_ls"
    phase = tools.raw_phase_key("RES|f0:GP:phi")
    norm = tools.projective_norm_key("RES|f0:GP:g_ls(0,0)")
    angles = [tools.projective_angle_key(base, index) for index in range(2)]
    names = [phase, angles[0], norm, angles[1]]
    bounds = np.array(
        [
            [-math.pi, math.pi],
            [-0.5 * math.pi, 0.5 * math.pi],
            [0.0, 3.0],
            [-0.5 * math.pi, 0.5 * math.pi],
        ],
        dtype=np.float64,
    )
    angle = 0.37
    values = np.array(
        [
            [0.2, angle, 2.0, 0.5 * math.pi],
            [0.2 - math.pi, -angle, 2.0, -0.5 * math.pi],
            [-2.5, -0.8, 0.0, 0.4],
            [1.2, 0.6, 0.0, -0.9],
        ],
        dtype=np.float64,
    )
    topology = tools.build_parameter_topology(names)
    features = icescape.model_features(values, bounds, names, topology)
    group = topology["groups"][0]

    assert group["kind"] == "polar_projective"
    assert group["parameters"] == [norm, phase, *angles]
    assert features.shape == (4, 6)
    assert np.allclose(features[0], features[1], atol=1.0e-12)
    assert np.allclose(features[2], features[3], atol=1.0e-12)


# Check the spherical resonance chart represents the physical complex LS vector
def test_polar_sphere_phase_zero():
    base = "RES|f0:GP:g_ls"
    norm = tools.projective_norm_key(f"{base}(0,0)")
    phase = tools.raw_phase_key("RES|f0:GP:phi")
    angles = [f"{base}({L},0){tools.SPHERICAL_VECTOR_SUFFIX}" for L in (2, 4)]
    names = [norm, phase, *angles]
    bounds = np.array([[0.0, 3.0], [-math.pi, math.pi], [-math.pi / 2, math.pi / 2], [-math.pi, math.pi]])
    direction = tools.normalize_spherical_vector([0.2, -0.3, 0.7])
    opposite_angles = tools.spherical_angles_from_vector([-value for value in direction])
    direct_angles = tools.spherical_angles_from_vector(direction)
    values = np.array(
        [
            [2.0, -2.0, *direct_angles],
            [2.0, -2.0 + math.pi, *opposite_angles],
            [0.0, -1.3, *direct_angles],
            [0.0, 0.8, *opposite_angles],
        ]
    )
    topology = tools.build_parameter_topology(names)

    features = icescape.model_features(values, bounds, names, topology)
    torch_features = icescape.model_features_torch(torch.tensor(values), bounds, names, topology)

    assert features.shape == (4, 6)
    assert np.allclose(features[0], features[1], atol=1.0e-12)
    assert np.allclose(features[2], features[3], atol=1.0e-12)
    assert np.allclose(torch_features.numpy(), features, atol=1.0e-12)


# Check GP continuum projective features use the normalized completed local coupling
def test_gp_con_residue_features():
    base = "CON_GP|990:[211,211]/opposite:helicity"
    norm = tools.projective_norm_key(f"{base}(0,0,-2)")
    angles = [f"{base}(0,0,{m}){tools.PROJECTIVE_VECTOR_SUFFIX}" for m in (-1, 0)]
    names = [norm, *angles]
    bounds = np.array([[0.0, 6.0], [-math.pi / 2, math.pi / 2], [-math.pi / 2, math.pi / 2]])
    values = np.array([[3.0, 0.31, -0.47], [0.0, -0.8, 0.9]])
    topology = tools.build_parameter_topology(names)

    features = icescape.model_features(values, bounds, names, topology)
    direction = np.asarray(tools.hemisphere_vector_from_angles(values[0, 1:].tolist()))
    weights = np.array([2.0, 2.0, 1.0])
    completed = np.sqrt(weights) * direction / np.sqrt(np.sum(weights * direction**2))
    expected = (values[0, 0] / 6.0) * np.asarray(tools.projective_embedding_from_vector(completed.tolist()))

    assert topology["groups"][0]["kind"] == "radial_projective"
    assert topology["groups"][0]["weights"] == weights.tolist()
    assert np.allclose(features[0], expected, atol=1.0e-12)
    assert np.allclose(features[1], 0.0, atol=1.0e-12)


# Check direct Cartesian continuum rows remain independent local couplings
def test_gp_con_cartesian_features_global_sign():
    base = "CON_GP|990:[211,211]/opposite:helicity"
    names = [
        key
        for m in (-1, 0)
        for key in (
            tools.complex_re_key(f"{base}(0,0,{m})"),
            tools.complex_im_key(f"{base}(0,0,{m})"),
        )
    ]
    bounds = np.array([[-3.0, 3.0]] * 4)
    point = np.array([0.4, -0.8, 1.2, 0.3])
    values = np.stack((point, -point, np.array([0.4, -0.8, -1.2, -0.3])))
    topology = tools.build_parameter_topology(names)

    features = icescape.model_features(values, bounds, names, topology)

    assert topology == {}
    assert not np.allclose(features[0], features[1], atol=1.0e-12)
    assert not np.allclose(features[0], features[2], atol=1.0e-12)

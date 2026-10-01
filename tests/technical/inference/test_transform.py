# Test autograd parameter maps, covariance propagation and constrained physical shifts
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
from pathlib import Path

import numpy as np
import pyjson5
import pytest
import torch
from core.inference import impact, profile
from core.stats.transform import ParameterTransform
from core.tune.drivers.graniitti.driver import GraniittiDriver
from core.tune.parameters import tools


# Decode a nonlinear map with more physical outputs than optimizer inputs
def _decode(values):
    x, y = values["x"], values["y"]
    return {"product": x * y, "sum": x + y, "fixed": 2.0}


# Check the full covariance including correlations and a rectangular Jacobian
def test_transform_full_covariance():
    transform = ParameterTransform(["x", "y"], _decode, [2.0, 3.0])
    covariance = np.array([[0.04, -0.01], [-0.01, 0.09]])
    result = transform.propagate([2.0, 3.0], covariance)
    jacobian = np.array([[3.0, 2.0], [1.0, 1.0], [0.0, 0.0]])
    np.testing.assert_allclose(result["jacobian"], jacobian)
    np.testing.assert_allclose(result["covariance"], jacobian @ covariance @ jacobian.T)
    assert result["errors"][0] == pytest.approx(math.sqrt(0.6))
    assert np.isnan(result["correlation"][2]).all()
    np.testing.assert_allclose(transform.samples([[2.0, 3.0], [4.0, 5.0]]), [[6.0, 5.0, 2.0], [20.0, 9.0, 2.0]])


# Preserve missing covariance instead of assigning artificial zero errors
def test_transform_missing_covariance():
    transform = ParameterTransform(["x", "y"], dict, [1.0, 2.0])
    result = transform.propagate([1.0, 2.0], [[0.1, np.nan], [np.nan, np.nan]])
    assert result["errors"][0] == pytest.approx(math.sqrt(0.1))
    assert math.isnan(result["errors"][1])


# Constrain the nonlinear physical value rather than a relabeled optimizer coordinate
@pytest.mark.parametrize("mode", ["fixed", "profiled"])
def test_transform_physical_impacts(mode):
    transform = ParameterTransform(["x", "y"], _decode, [2.0, 3.0])
    covariance = np.array([[0.04, 0.01], [0.01, 0.09]])
    precision = torch.linalg.inv(torch.tensor(covariance))

    # Evaluate a correlated Gaussian likelihood using autograd
    def objective(point):
        delta = torch.tensor(point - np.array([2.0, 3.0]), requires_grad=True)
        value = delta @ precision @ delta
        gradient, = torch.autograd.grad(value, delta)
        return value.item(), gradient.detach().numpy()

    shifts = impact.physical_parameter_shifts(transform=transform, best_fit=np.array([2.0, 3.0]),
        covariance=covariance, bounds=np.array([[0.5, 5.0], [0.5, 6.0]]), mode=mode,
        maxiter=100, value_gradient_func=objective)
    assert shifts["valid_minus"].tolist() == [True, True, False]
    assert shifts["valid_plus"].tolist() == [True, True, False]
    for side in ("minus", "plus"):
        actual = transform.samples(shifts[f"{side}_points"])
        np.testing.assert_allclose(np.diag(actual)[:2], shifts[f"shift_{side}"][:2], atol=1e-6)
    responses = impact.compute_direct_parameter_impacts(
        evaluate_observable_chi2=lambda points: np.sum(points**2, axis=1, keepdims=True),
        best_fit=np.array([2.0, 3.0]), bounds=np.array([[0.5, 5.0], [0.5, 6.0]]), shifts=shifts)
    assert responses["delta_plus"].shape == (1, 3)


# Cross a phase branch with periodic constraints and an autograd Jacobian
def test_transform_periodic_constraint():
    # Decode the canonical phase through a differentiable trigonometric map
    def decode(values):
        phase = values["angle"]
        return {"phase": torch.atan2(torch.sin(phase), torch.cos(phase))}

    transform = ParameterTransform(["angle"], decode, [3.1], {"phase": 2.0 * math.pi})
    result = transform.constrain([3.1], [[0.04]], [[2.0, 4.0]], 0, 0.2, maxiter=100)
    assert result["success"]
    assert result["point"][0] == pytest.approx(3.3)
    assert transform.propagate([3.1], [[0.04]])["errors"][0] == pytest.approx(0.2)


# Check the real projective steering decoder at complex amplitude level
@pytest.mark.parametrize("angle", [-0.3, 0.3])
def test_graniitti_transform_projective(angle):
    root = Path(__file__).resolve().parents[3]
    base = "RES|f0_1500:GP:g_ls"
    config = {f"{base}(0,0)@NORM": 2.0, f"{base}(2,2)@PROJECTIVE": angle}
    driver = GraniittiDriver()
    transform = driver.parameter_transform(list(config), list(config.values()), cdir=str(root),
                                           metadata={"mc_steer": {"tune_default": "TUNE0"}})
    result = transform.propagate(list(config.values()), np.eye(2))
    values = dict(zip(result["names"], result["values"], strict=True))
    assert all(":rows" not in name for name in values)
    assert f"{base}(4,4)@MAG" not in values
    expected = [2.0 * math.cos(angle), 2.0 * math.sin(angle)]
    for labels, value in zip(("0,0", "2,2"), expected, strict=True):
        amplitude = values[f"{base}({labels})@MAG"] * np.exp(1j * values[f"{base}({labels})@PHASE"])
        assert amplitude == pytest.approx(value, rel=1e-12, abs=1e-12)
    assert torch.autograd.gradcheck(transform, (torch.tensor(list(config.values()), dtype=torch.float64, requires_grad=True),))


# Check parallel physical profiles retain grid order and satisfy simultaneous constraints
@pytest.mark.parametrize("dimension", [1, 2])
def test_transformed_profile(dimension):
    # Decode two coupled physical coordinates
    def decode(values):
        return {"sum": values["x"] + values["y"], "difference": values["x"] - values["y"]}

    # Compute a Gaussian likelihood and its autograd gradient
    def objective(values):
        point = torch.tensor(values, requires_grad=True)
        value = point.square().sum()
        gradient, = torch.autograd.grad(value, point)
        return value.item(), gradient.detach().numpy()

    transform = ParameterTransform(["x", "y"], decode, [0.0, 0.0])
    grid = np.array([-1.0, 0.0, 2.0])
    values = profile.transformed_profile(transform=transform, indices=tuple(range(dimension)), grids=[grid] * dimension,
        best_fit=np.zeros(2), bounds=np.array([[-5.0, 5.0], [-5.0, 5.0]]), value_gradient_func=objective,
        maxiter=100, workers=2)
    expected = grid**2 / 2.0 if dimension == 1 else (grid[:, None]**2 + grid[None, :]**2) / 2.0
    np.testing.assert_allclose(values, expected, atol=1e-6)


# Check all supported coordinate decoders preserve autograd instead of detaching inputs
@pytest.mark.parametrize("config", [
    {tools.phase_u_key("RES|f0_1500:GP:phi"): 0.2, tools.phase_v_key("RES|f0_1500:GP:phi"): 0.4},
    {"CON_GP|990:[211,211]/opposite:helicity(0,0,0)@NORM": 2.0},
    {"RES|f0_1500:TP:[995,995]:g_tensor(0,0)@NORM": 2.0,
     "RES|f0_1500:TP:[995,995]:g_tensor(2,2)@PROJECTIVE": 0.3},
    {tools.spherical_angle_key("RES|f2_1270:MP:polarization.a_Jz", 0): 0.3},
    {tools.simplex_theta_key("RES|f2_1270:MP:polarization.rho_population", 0): 0.4,
     tools.simplex_theta_key("RES|f2_1270:MP:polarization.rho_population", 1): 0.7},
    {"SOFT|MODEL.double:EXCHANGE.P.g[0,0]@DESCENDING": 2.0,
     "SOFT|MODEL.double:EXCHANGE.P.g[1,1]@DESCENDING": 0.4},
    {"SOFT|MODEL.double:EXCHANGE.P.g[0,1]@SYMMETRIC": 0.5},
    {"SOFT|MODEL.double:FF.p.param[0,0]@ORDERED_DPOW=5": 1.0,
     "SOFT|MODEL.double:FF.p.param[0,1]@ORDERED_DPOW=5": 0.3},
])
def test_graniitti_transform_autograd(config):
    root = Path(__file__).resolve().parents[3]
    transform = GraniittiDriver().parameter_transform(list(config), list(config.values()), cdir=str(root),
        metadata={"mc_steer": {"tune_default": "TUNE0"}})
    point = torch.tensor(list(config.values()), dtype=torch.float64, requires_grad=True)
    assert torch.autograd.gradcheck(transform, (point,))
    result = transform.propagate(point.detach(), np.eye(len(config)))
    assert np.any(result["errors"] > 0.0)
    assert all(":rows" not in name for name in result["names"])


# Verify two-index helicity names and Cartesian amplitudes through the real decoder
def test_graniitti_transform_helicity(tmp_path):
    root = Path(__file__).resolve().parents[3]
    tune = tmp_path / "modeldata" / "TEST"
    (tune / "RES").mkdir(parents=True)
    card = pyjson5.loads((root / "modeldata/TUNE0/RES/f0_1500.json").read_text())
    block = card["PARAM_RES"]["MODELS"]["GP"]["[990,990]"]
    block["basis"] = "helicity"
    block["helicity"] = [[0, 0, 1.0, 0.0], [1, 1, 1.0, 0.0]]
    (tune / "RES/f0_1500.json").write_text(pyjson5.dumps(card))
    base = "RES|f0_1500:GP:helicity(0,0)"
    config = {tools.complex_re_key(base): 0.3, tools.complex_im_key(base): 0.4}
    transform = GraniittiDriver().parameter_transform(list(config), list(config.values()), cdir=str(tmp_path),
        metadata={"mc_steer": {"tune_default": "TEST"}})
    point = torch.tensor(list(config.values()), dtype=torch.float64, requires_grad=True)
    assert transform.names == (f"{base}@MAG", f"{base}@PHASE")
    magnitude, phase = transform(point).detach().numpy()
    assert magnitude * np.exp(1j * phase) == pytest.approx(0.3 + 0.4j)
    assert torch.autograd.gradcheck(transform, (point,))

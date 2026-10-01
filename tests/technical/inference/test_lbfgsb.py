# Tests for GPU capable bounded Torch L-BFGS optimization
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np
import pytest
import torch
from core.numerics import lbfgsb


# Check finite, one-sided, fixed and unbounded coordinates in the same solve
@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize("bounded", [False, True])
def test_general_bounds(batched, bounded):
    target = torch.tensor([-3.0, 4.0, 2.0, -1.0], dtype=torch.float64)
    bounds = (torch.tensor([-torch.inf, -2.0, 0.5, -torch.inf]),
              torch.tensor([-1.0, torch.inf, 0.5, -2.0])) if bounded else None
    expected = torch.tensor([-3.0, 4.0, 0.5, -2.0], dtype=torch.float64) if bounded else target
    start = torch.zeros((2, 4) if batched else (4,), dtype=torch.float64)

    # Compute independent anisotropic quadratics on the physical coordinate scale
    def objective(values):
        return ((values - target).square() * values.new_tensor([1.0, 2.0, 4.0, 8.0])).sum(dim=-1)

    minimize = lbfgsb.minimize_batch if batched else lbfgsb.minimize
    values, _, _, diagnostics = minimize(objective, start, bounds=bounds, maxiter=200)
    assert torch.as_tensor(diagnostics["success"]).all()
    torch.testing.assert_close(values, expected.expand_as(values), atol=1e-7, rtol=1e-7)


# Check optional relative objective convergence with a finite likelihood offset
@pytest.mark.parametrize("batched", [False, True])
def test_relative_objective_tolerance(batched):
    start = torch.tensor([[0.8]] if batched else [0.8], dtype=torch.float64)

    # Include a constant likelihood normalization without changing its gradient
    def objective(values):
        return 1e6 + (values - 0.4).square().sum(dim=-1)

    minimize = lbfgsb.minimize_batch if batched else lbfgsb.minimize
    _, _, gradient, diagnostics = minimize(objective, start, learning_rate=0.25, function_tolerance=1e-6)
    assert torch.as_tensor(diagnostics["success"]).all()
    assert gradient.abs().min() > 1e-7
    messages = diagnostics["messages"] if batched else [diagnostics["message"]]
    assert messages == ["relative objective tolerance reached"]
    values, _, _, diagnostics = minimize(objective, start, learning_rate=0.25)
    assert torch.as_tensor(diagnostics["success"]).all()
    torch.testing.assert_close(values, torch.full_like(values, 0.4), atol=1e-7, rtol=1e-7)


# Check coupled curvature with different active bounds against analytic KKT minima
@pytest.mark.parametrize("batched", [False, True])
def test_coupled_bound_minima(batched):
    factor = torch.tensor([[2, 0, 0, 0], [1, 3, 0, 0], [-2, 1, 2, 0], [1, -1, 2, 1]], dtype=torch.float64)
    hessian = factor @ factor.T
    expected = torch.tensor([[0.0, 0.3, 1.0, 0.8], [0.2, 0.0, 0.7, 1.0]], dtype=torch.float64)
    kkt = torch.tensor([[1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0]], dtype=torch.float64)
    target = expected - torch.linalg.solve(hessian, kkt.T).T
    if not batched:
        expected, target = expected[0], target[0]

    # Compute a quadratic whose unconstrained minimum lies outside the unit box
    def objective(values):
        residual = values - target
        return 0.5 * ((residual @ hessian) * residual).sum(dim=-1)

    minimize = lbfgsb.minimize_batch if batched else lbfgsb.minimize
    values, _, _, diagnostics = minimize(objective, torch.full_like(expected, 0.5), maxiter=200)
    assert torch.as_tensor(diagnostics["success"]).all()
    torch.testing.assert_close(values, expected, atol=1e-7, rtol=1e-7)


# Preserve interior curvature motion while rejecting outward and active-bound directions
def test_project_curvature_direction():
    values = torch.tensor([0.5, 0.0, 1.0, 0.0], dtype=torch.float64)
    gradient = torch.tensor([0.0, -1.0, -1.0, 1.0], dtype=torch.float64)
    direction = torch.tensor([0.2, -0.3, -0.4, 0.5], dtype=torch.float64)
    result = lbfgsb._project_direction(values, gradient, direction)
    torch.testing.assert_close(result, torch.tensor([0.2, 0.0, 0.0, 0.0], dtype=torch.float64))


# Check configurable step scales and sufficient decrease on an analytic quadratic
@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize(("controls", "expected"), [
    ({}, 0.6),
    ({"learning_rate": 0.25}, 0.4),
    ({"line_search_decay": 0.25}, 0.4),
    ({"line_search_armijo": 0.9}, 0.25),
])
def test_line_search_controls(batched, controls, expected):
    start = torch.tensor([[0.2]] if batched else [0.2], dtype=torch.float64)

    # Compute the same quadratic for scalar and batched optimizer calls
    def objective(values):
        return (values - 0.6).square().sum(dim=-1)

    minimize = lbfgsb.minimize_batch if batched else lbfgsb.minimize
    values, value, _, _ = minimize(objective, start, maxiter=1, **controls)
    assert values.item() == pytest.approx(expected)
    assert value.item() == pytest.approx((expected - 0.6)**2)


# Test bounded L-BFGS converges to an interior quadratic minimum
def test_minimize_quadratic_selected_device():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    target = torch.tensor([0.25, 0.75], dtype=torch.float64, device=device)

    # Evaluate one positive definite quadratic in unit coordinates
    def objective(values):
        residual = values - target
        return residual[0].square() + 4.0 * residual[1].square()

    values, objective_value, gradient, diagnostics = lbfgsb.minimize(
        objective=objective,
        start=torch.tensor([0.9, 0.1], dtype=torch.float64, device=device),
        maxiter=100,
    )

    assert diagnostics["success"]
    assert diagnostics["device"] == str(device)
    assert objective_value.item() == pytest.approx(0.0, abs=1.0e-12)
    assert torch.max(torch.abs(gradient)).item() < 1.0e-7
    np.testing.assert_allclose(values.cpu(), target.cpu(), atol=1.0e-7)


# Test projected gradients admit a valid minimum on a physical bound
def test_minimize_accepts_bound_minimum():
    # Place the unconstrained minimum below the lower bound
    def objective(values):
        return (values[0] + 1.0).square()

    values, _, gradient, diagnostics = lbfgsb.minimize(
        objective=objective,
        start=torch.tensor([0.5], dtype=torch.float64),
    )

    assert diagnostics["success"]
    assert values.item() == pytest.approx(0.0, abs=1.0e-12)
    assert gradient.item() > 0.0
    assert diagnostics["max_abs_projected_unit_gradient"] == pytest.approx(0.0)


# Test relative KKT tolerance scales from the initial projected gradient
def test_relative_projected_gradient():
    # Evaluate one steep quadratic with a large initial unit gradient
    def objective(values):
        return 1000.0 * (values[0] - 0.25).square()

    _, _, _, diagnostics = lbfgsb.minimize(
        objective=objective,
        start=torch.tensor([0.9], dtype=torch.float64),
        relative_gradient_tolerance=1.0e-6,
    )

    assert diagnostics["initial_max_abs_projected_unit_gradient"] == pytest.approx(1300.0)
    assert diagnostics["effective_gradient_tolerance"] == pytest.approx(1.3e-3)
    assert diagnostics["success"]


# Test a failed quasi-Newton search retries after clearing its curvature history
def test_line_search_history_reset(monkeypatch):
    target = torch.tensor([0.2, 0.8], dtype=torch.float64)

    # Evaluate an anisotropic quadratic that requires multiple L-BFGS steps
    def objective(values):
        residual = values - target
        return residual[0].square() + 100.0 * residual[1].square()

    original_search = lbfgsb._bounded_line_search
    calls = 0

    # Inject one failure after the first accepted step
    def intermittent_search(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            return None
        return original_search(*args, **kwargs)

    monkeypatch.setattr(lbfgsb, "_bounded_line_search", intermittent_search)
    values, _, _, diagnostics = lbfgsb.minimize(
        objective=objective,
        start=torch.tensor([0.9, 0.1], dtype=torch.float64),
        maxiter=100,
    )

    assert diagnostics["success"]
    assert diagnostics["history_resets"] == 1
    np.testing.assert_allclose(values, target, atol=1.0e-7)


# Test independent profile points converge in one batched GPU optimization
def test_minimize_batch_converges_independent_rows():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    targets = torch.tensor(
        [[0.1, 0.8], [0.3, 0.6], [0.5, 0.4], [0.7, 0.2]],
        dtype=torch.float64,
        device=device,
    )

    # Compute one independent quadratic objective per batch row
    def objective(values):
        residual = values - targets
        return residual[:, 0].square() + 4.0 * residual[:, 1].square()

    values, objective_values, _, diagnostics = lbfgsb.minimize_batch(
        objective=objective,
        start=torch.full_like(targets, 0.5),
        maxiter=100,
    )

    assert torch.all(diagnostics["success"])
    assert diagnostics["device"] == str(device)
    assert diagnostics["batch_evaluations"] < diagnostics["row_evaluations"]
    assert torch.max(torch.abs(objective_values)).item() < 1.0e-12
    np.testing.assert_allclose(values.cpu(), targets.cpu(), atol=1.0e-7)


# Test batched L-BFGS reports each completed global descent iteration
def test_minimize_batch_reports_progress():
    target = torch.tensor([[0.2], [0.8]], dtype=torch.float64)
    reports = []

    # Compute one quadratic objective per independent batch row
    def objective(values):
        return (values - target).square().reshape(-1)

    lbfgsb.minimize_batch(
        objective=objective,
        start=torch.full_like(target, 0.5),
        maxiter=20,
        progress_callback=lambda iteration, active, best: reports.append(
            (iteration, active, best)
        ),
    )

    assert reports
    assert [report[0] for report in reports] == list(range(1, len(reports) + 1))
    assert all(report[1] >= 0 for report in reports)
    assert all(np.isfinite(report[2]) for report in reports)


# Test exact autograd curvature gives the likelihood covariance convention
@pytest.mark.parametrize("vectorize", [False, True])
def test_hessian_cov_errordef_convention(vectorize):
    # Evaluate a quadratic with known physical Hessian
    def objective(values):
        return (values[0] - 0.25).square() + 4.0 * (values[1] + 0.5).square()

    point = torch.tensor([0.25, -0.5], dtype=torch.float64)
    hessian, covariance, diagnostics = lbfgsb.hessian_covariance(objective, point, vectorize=vectorize)

    assert diagnostics["positive_definite"]
    assert diagnostics["regularized_eigenvalues"] == 0
    np.testing.assert_allclose(hessian, np.diag([2.0, 8.0]), atol=1.0e-12)
    np.testing.assert_allclose(covariance, np.diag([1.0, 0.25]), atol=1.0e-12)


# Test constrained curvature excludes active-bound directions from covariance
def test_hessian_cov_free_param_subspace():
    # Use negative curvature only along a coordinate fixed at its upper bound
    def objective(values):
        return -values[0].square() + values[1].square()

    point = torch.tensor([1.0, 0.0], dtype=torch.float64)
    hessian, covariance, diagnostics = lbfgsb.hessian_covariance(
        objective,
        point,
        free_mask=torch.tensor([False, True]),
    )

    assert not diagnostics["positive_definite"]
    assert diagnostics["free_positive_definite"]
    assert diagnostics["active_dimension"] == 1
    np.testing.assert_allclose(hessian, np.diag([-2.0, 2.0]), atol=1.0e-12)
    assert np.isnan(covariance[0, 0].item())
    assert covariance[1, 1].item() == pytest.approx(1.0)


# Reject finite objective proposals whose derivatives diverge at a boundary
@pytest.mark.parametrize("batched", [False, True])
def test_line_search_rejects_infinite_gradient(batched):
    # The square root has finite values and an infinite derivative at zero
    def objective(values):
        return values.sqrt().sum(dim=-1)

    start = torch.tensor([[0.25]] if batched else [0.25], dtype=torch.float64)
    minimize = lbfgsb.minimize_batch if batched else lbfgsb.minimize
    values, objective_value, gradient, _ = minimize(objective, start, maxiter=2)
    assert torch.isfinite(gradient).all()
    assert torch.isfinite(objective_value).all()
    assert torch.all(values > 0.0)
    assert torch.all(values < start)


# A flat derivative does not make an infinite objective a converged minimum
def test_infinite_initial_objective():
    # Construct an impossible likelihood with a finite zero derivative
    def objective(values):
        return values.sum() * 0.0 + torch.inf

    _, _, _, diagnostics = lbfgsb.minimize(objective, torch.tensor([0.5], dtype=torch.float64))
    assert not diagnostics["success"]


# Recognize convergence on the final permitted scalar descent step
def test_scalar_minimize_last_step_convergence():
    # This quadratic reaches its minimum in one full Newton-scaled step
    def objective(values):
        return 0.5 * (values - 0.4).square().sum()

    _, objective_value, _, diagnostics = lbfgsb.minimize(
        objective, torch.tensor([0.8], dtype=torch.float64), maxiter=1,
    )
    assert objective_value.item() == pytest.approx(0.0)
    assert diagnostics["success"]


# Reject Armijo trials without backward passes or repeated forward evaluations
@pytest.mark.parametrize("batched", [False, True])
def test_deferred_line_search_derivatives(batched):
    counts = {"forward": 0, "backward": 0}
    start = torch.tensor([[0.9, 0.8]] if batched else [0.9, 0.8], dtype=torch.float64)

    # Count real autograd execution on an anisotropic quadratic
    def objective(values):
        counts["forward"] += 1
        loss = ((values - 0.3).square() * values.new_tensor([1.0, 20.0])).sum(dim=-1)

        # Observe backward execution without replacing the calculated derivative
        def backward(gradient):
            counts["backward"] += 1
            return gradient

        loss.register_hook(backward)
        return loss

    minimize = lbfgsb.minimize_batch if batched else lbfgsb.minimize
    values, loss, _, diagnostics = minimize(objective, start, bounds=None, maxiter=100)
    assert torch.as_tensor(diagnostics["success"]).all()
    torch.testing.assert_close(values, torch.full_like(values, 0.3), atol=1e-7, rtol=1e-7)
    assert loss.abs().max() < 1e-12
    assert counts["backward"] < counts["forward"]
    if not batched:
        assert diagnostics["function_evaluations"] == counts["forward"]
        assert diagnostics["gradient_evaluations"] == counts["backward"]


# Keep eager external derivatives and deferred autograd on the same optimization path
def test_deferred_matches_eager_evaluator():
    start = torch.tensor([0.8, 0.1], dtype=torch.float64)

    # Compute a nonlinear coupled objective with an interior minimum
    def objective(values):
        return 100.0 * (values[1] - values[0].square()).square() + (1.0 - values[0]).square()

    # Exercise the public evaluator API with actual Torch derivatives
    def evaluate(values):
        return lbfgsb.value_gradient(objective, values)

    lazy = lbfgsb.minimize(objective, start, maxiter=300)
    eager = lbfgsb.minimize(None, start, evaluator=evaluate, maxiter=300)
    for actual, expected in zip(lazy[:3], eager[:3], strict=True):
        torch.testing.assert_close(actual, expected, rtol=0.0, atol=0.0)
    assert lazy[3]["gradient_evaluations"] < eager[3]["gradient_evaluations"]


# Reject a finite objective with an infinite derivative and continue backtracking
@pytest.mark.parametrize("batched", [False, True])
def test_deferred_nonfinite_derivative(batched):
    start = torch.tensor([[0.25]] if batched else [0.25], dtype=torch.float64)

    # The lower bound has a finite value and an infinite slope
    def objective(values):
        return values.sqrt().sum(dim=-1)

    minimize = lbfgsb.minimize_batch if batched else lbfgsb.minimize
    values, loss, gradient, _ = minimize(objective, start, maxiter=1)
    torch.testing.assert_close(values, torch.full_like(values, 0.125))
    assert torch.isfinite(loss).all() and torch.isfinite(gradient).all()


# Check compact scalar and batched directions against explicit inverse BFGS updates
@pytest.mark.parametrize(("count", "scaled"), [(0, False), (1, False), (20, False), (20, True)])
def test_compact_direction(count, scaled):
    rng = torch.Generator().manual_seed(17)
    gradient = torch.randn(3, 8, generator=rng, dtype=torch.float64)
    factor = torch.randn(8, 8, generator=rng, dtype=torch.float64)
    hessian = factor @ factor.T + torch.eye(8)
    history = []
    for index in range(count):
        step = torch.randn(3, 8, generator=rng, dtype=torch.float64) * (10.0 ** -index if scaled else 1.0)
        history.append((step, step @ hessian, torch.tensor([index % 2 == 0, False, True])))
    actual = lbfgsb._batch_lbfgs_direction(gradient, history)
    for row in range(3):
        pairs = [(step[row], change[row]) for step, change, valid in history if valid[row]]
        eye = torch.eye(8, dtype=torch.float64)
        inverse = (torch.dot(*pairs[-1]) / pairs[-1][1].square().sum() if pairs else 1.0) * eye
        for step, change in pairs:
            rho = torch.reciprocal(torch.dot(step, change))
            update = eye - rho * torch.outer(step, change)
            inverse = update @ inverse @ update.T + rho * torch.outer(step, step)
        expected = -inverse @ gradient[row]
        torch.testing.assert_close(actual[row], expected, atol=1e-11, rtol=1e-11)
        scalar = lbfgsb._lbfgs_direction(gradient[row], [s for s, _ in pairs], [y for _, y in pairs])
        torch.testing.assert_close(scalar, expected, atol=1e-11, rtol=1e-11)

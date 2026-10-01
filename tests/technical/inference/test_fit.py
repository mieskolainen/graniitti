# Tests for surrogate fit optimization and diagnostics
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path

import numpy as np
import pytest
import torch
from core.inference import fit
from core.stats import objective as icecost


# Check repeated covariance eigenvalues retain finite first and second derivatives
@pytest.mark.parametrize("error", [0.0, 0.5])
@pytest.mark.parametrize("scale", [1.0e-20, 1.0, 1.0e100])
def test_correlated_chi2_degenerate_cov_derivatives(error, scale):
    counts = torch.tensor([[2.0, 3.0]], dtype=torch.float64, requires_grad=True)
    errors = torch.full_like(counts, error, requires_grad=True)

    # Evaluate the full covariance objective through the production Torch API
    def objective(counts, errors):
        return icecost.histogram_chi2_contributions_torch(
            counts=counts * scale, errors=errors * scale, data_counts=torch.ones(2, dtype=torch.float64) * scale,
            data_errors=torch.ones(2, dtype=torch.float64) * scale, mask=torch.ones(2, dtype=torch.bool),
            fit_weights=torch.ones(2, dtype=torch.float64),
            data_total_covariance=torch.eye(2, dtype=torch.float64) * scale**2,
            covariance_indices=torch.arange(2),
        ).sum()

    assert objective(counts, errors).item() == pytest.approx(5.0 / (1.0 + error**2))
    assert torch.autograd.gradcheck(objective, (counts, errors))
    assert torch.autograd.gradgradcheck(objective, (counts, errors))


# Check singular covariance accepts supported residuals and rejects discarded components
def test_correlated_chi2_singular_cov_support():
    contribution = icecost.histogram_chi2_contributions_torch(
        counts=torch.tensor([[2.0, 0.0], [2.0, 2.0]], dtype=torch.float64),
        errors=torch.zeros((2, 2), dtype=torch.float64), data_counts=torch.ones(2, dtype=torch.float64),
        data_errors=torch.ones(2, dtype=torch.float64), mask=torch.ones(2, dtype=torch.bool),
        fit_weights=torch.ones(2, dtype=torch.float64),
        data_total_covariance=torch.tensor([[1.0, -1.0], [-1.0, 1.0]], dtype=torch.float64),
        covariance_indices=torch.arange(2),
    )
    assert contribution[0].sum().item() == pytest.approx(1.0)
    assert torch.isposinf(contribution[1].sum())


# Test the Torch optimizer persists exact Hessian parameter uncertainties
def test_optimize_torch_writes_param_payload(tmp_path):
    # Evaluate one vectorized physical quadratic for sampled-point diagnostics
    def objective_numpy(values):
        values = np.asarray(values, dtype=np.float64)
        return (values[:, 0] - 0.25) ** 2 + 4.0 * (values[:, 1] + 0.5) ** 2

    # Evaluate the same quadratic without breaking the Torch graph
    def objective_torch(values):
        return (values[0] - 0.25).square() + 4.0 * (values[1] + 0.5).square()

    # Evaluate every independent quadratic restart in one Torch batch
    def objective_batch_torch(values):
        return (values[:, 0] - 0.25).square() + 4.0 * (values[:, 1] + 0.5).square()

    X = np.array([[-1.0, -1.0], [0.25, -0.5], [0.5, 0.0], [1.0, 1.0]])
    Z = objective_numpy(X)
    result, best = fit.optimize_with_torch(
        X=X,
        Z=Z,
        param_names=["x", "y"],
        func_hat=objective_numpy,
        objective_torch=objective_torch,
        objective_batch_torch=objective_batch_torch,
        bounds=np.array([[-1.0, 1.0], [-1.0, 1.0]]),
        output_root=tmp_path,
        device=torch.device("cpu"),
        dtype=torch.float64,
        realized_best=np.array([0.25, -0.5]),
        realized_Z=0.0,
        x0=np.array([-0.5, 0.5]),
    )

    assert result.fmin.is_valid
    np.testing.assert_allclose(best, [0.25, -0.5], atol=1.0e-7)
    np.testing.assert_allclose(result.covariance, np.diag([1.0, 0.25]), atol=1.0e-8)
    with open(tmp_path / "parameters.json", encoding="utf-8") as handle:
        payload = json.load(handle)
    assert payload["best_fit"]["uncertainty_method"] == "torch_autograd_hessian"
    assert payload["optimization_diagnostics"]["backend"] == "torch"
    assert payload["optimization_diagnostics"]["curvature"]["positive_definite"]
    assert payload["optimization_diagnostics"]["global_search"]["selected_restarts"] > 1


# Test batched global starts escape a deliberately poor local basin
def test_optimize_torch_escapes_local_minimum(tmp_path):
    # Evaluate a tilted double-well objective with a unique global minimum
    def objective_numpy(values):
        x = np.asarray(values, dtype=np.float64)[:, 0]
        return (x - 0.2) ** 2 * (x - 0.75) ** 2 + 0.01 * (x - 0.2) ** 2

    # Evaluate one double-well point without leaving the Torch graph
    def objective_torch(values):
        x = values[0]
        return (x - 0.2).square() * (x - 0.75).square() + 0.01 * (x - 0.2).square()

    # Evaluate all independent double-well restarts without leaving the Torch graph
    def objective_batch_torch(values):
        x = values[:, 0]
        return (x - 0.2).square() * (x - 0.75).square() + 0.01 * (x - 0.2).square()

    X = np.asarray([[0.65], [0.75], [0.85]], dtype=np.float64)
    Z = objective_numpy(X)
    result, best = fit.optimize_with_torch(
        X=X,
        Z=Z,
        param_names=["x"],
        func_hat=objective_numpy,
        objective_torch=objective_torch,
        objective_batch_torch=objective_batch_torch,
        bounds=np.asarray([[0.0, 1.0]], dtype=np.float64),
        output_root=tmp_path,
        device=torch.device("cpu"),
        dtype=torch.float64,
        x0=np.asarray([0.75]),
        restarts=8,
        scan_points=64,
        rngseed=9,
    )

    assert result.fmin.is_valid
    np.testing.assert_allclose(best, [0.2], atol=1.0e-6)


# Check real constrained fitting and covariance extraction retain the correlated free subspace
def test_torch_bound_curvature(tmp_path):
    # The free covariance is [[4, 1], [1, 9]] and x is constrained at its upper bound
    def objective(values):
        x, y, z = values[..., 0], values[..., 1], values[..., 2]
        return -x**2 + (9*y**2 - 2*y*z + 4*z**2) / 35

    points = np.array([[0., 0., 0.], [.5, -.5, .5], [1., 0., 0.]])
    result, best = fit.optimize_with_torch(
        X=points, Z=objective(points), param_names=['x', 'y', 'z'], func_hat=objective,
        objective_torch=objective, objective_batch_torch=objective,
        bounds=np.array([[0., 1.], [-1., 1.], [-1., 1.]]), output_root=tmp_path,
        device='cpu', dtype=torch.float64, restarts=4, scan_points=16,
    )
    assert result.fmin.is_valid
    np.testing.assert_allclose(best, [1., 0., 0.], atol=1e-7)
    covariance, correlation, errors, _ = fit.fitted_covariance_summary(
        m=result, param_names=['x', 'y', 'z'], surrogate_vals=best)
    np.testing.assert_allclose(covariance[1:, 1:], [[4., 1.], [1., 9.]], atol=1e-8)
    np.testing.assert_allclose(correlation[1:, 1:], [[1., 1/6], [1/6, 1.]], atol=1e-8)
    np.testing.assert_allclose(errors[1:], [2., 3.], atol=1e-8)
    assert np.all(np.isnan(correlation[0])) and np.isnan(errors[0])
    payload = json.loads((tmp_path / 'parameters.json').read_text())
    curvature = payload['optimization_diagnostics']['curvature']
    assert not curvature['positive_definite'] and curvature['free_positive_definite']


# Configure and compute on each available device through the real Torch runtime
@pytest.mark.parametrize('device', ['cpu'] + [f'cuda:{i}' for i in range(torch.cuda.device_count())])
def test_configure_torch_runtime(device):
    fit.configure_torch_runtime(device)
    values = torch.tensor([3., 4.], device=device)
    assert torch.linalg.vector_norm(values).item() == pytest.approx(5.)
    assert str(values.device) == device


# Test shared torch optimization and persisted uncertainty payload
def test_torch_fit_from_displaced_start(tmp_path):
    def objective(X):
        X = np.asarray(X, dtype=np.float64)
        return (X[:, 0] - 0.25) ** 2 + 4.0 * (X[:, 1] + 0.5) ** 2

    X = np.array(
        [
            [-1.0, -1.0],
            [0.25, -0.5],
            [0.5, 0.0],
            [1.0, 1.0],
        ],
        dtype=np.float64,
    )
    Z = objective(X)
    output_root = tmp_path / "fit"
    # Compute the quadratic objective with exact torch derivatives
    def objective_torch(values):
        return (values[..., 0] - 0.25) ** 2 + 4.0 * (values[..., 1] + 0.5) ** 2

    result, best = fit.optimize_with_torch(
        X=X,
        Z=Z,
        param_names=["x", "y"],
        func_hat=objective,
        bounds=np.array([[-1.0, 1.0], [-1.0, 1.0]], dtype=np.float64),
        output_root=output_root,
        realized_best=np.array([0.25, -0.5], dtype=np.float64),
        realized_Z=0.0,
        x0=np.array([-0.5, 0.5], dtype=np.float64),
        objective_torch=objective_torch, objective_batch_torch=objective_torch,
        device="cpu", dtype=torch.float64, restarts=2, scan_points=4,
    )

    assert np.allclose(best, [0.25, -0.5], atol=1e-4)
    assert result.covariance is not None

    with open(output_root / "parameters.json", encoding="utf-8") as handle:
        payload = json.load(handle)

    assert payload["surrogate_best"]["x"]["uncertainty"] >= 0.0
    assert payload["best_fit"]["schema_version"] == 1
    assert payload["best_fit"]["source"] == "surrogate"
    assert payload["best_fit"]["objective"]["name"] == "Z"
    assert payload["best_fit"]["parameters"]["x"]["value"] == pytest.approx(0.25)
    assert payload["best_fit"]["parameters"]["x"]["uncertainty"] >= 0.0
    assert payload["realized_best"] == {"x": 0.25, "y": -0.5}
    assert payload["optimization_diagnostics"]["realized_residual"] == 0.0
    assert payload["optimization_diagnostics"]["start"]["source"] == "explicit"
    assert payload["optimization_diagnostics"]["derivatives"] == "torch autograd"
    assert payload["optimization_diagnostics"]["hessian"] == "torch autograd"
    assert Path(payload["parameter_uncertainty_plot"]["png"]).is_file()
    assert Path(payload["parameter_uncertainty_plot"]["pdf"]).is_file()


# Test shared bound diagnostics distinguish KKT-active and free directions
def test_bound_diagnostics_reports_active_kkt_sides():
    diagnostics = fit.parameter_bound_diagnostics(
        values=np.array([0.0, 1.0, 0.5]),
        bounds=np.array([[0.0, 1.0]] * 3),
        gradient=np.array([2.0, -3.0, 0.25]),
    )

    assert diagnostics[0]["status"] == "lower"
    assert diagnostics[0]["kkt_outward"]
    assert diagnostics[0]["projected_unit_gradient"] == 0.0
    assert diagnostics[1]["status"] == "upper"
    assert diagnostics[1]["kkt_outward"]
    assert diagnostics[2]["status"] == "interior"
    assert diagnostics[2]["projected_unit_gradient"] == pytest.approx(0.25)


# Test surrogate support reports distance and optimism relative to realized trials
def test_unvalidated_surrogate_valley():
    diagnostics = fit.surrogate_support_diagnostics(
        X=np.array([[0.0, 0.0], [1.0, 1.0]]),
        Z=np.array([10.0, 12.0]),
        values=np.array([0.25, 0.25]),
        bounds=np.array([[0.0, 1.0], [0.0, 1.0]]),
        surrogate_z=4.0,
    )

    assert diagnostics["nearest_trial_index"] == 0
    assert diagnostics["nearest_trial_unit_distance"] == pytest.approx(np.sqrt(0.125))
    assert diagnostics["nearest_trial_unit_rms_distance"] == pytest.approx(0.25)
    assert diagnostics["nearest_trial_z"] == 10.0
    assert diagnostics["realized_minimum_minus_surrogate_z"] == 6.0
    assert diagnostics["surrogate_z_over_realized_minimum"] == pytest.approx(0.4)


# Test relative uncertainty is undefined only at scale-level numerical zero
def test_relative_error_near_zero():
    relative = fit.parameter_relative_uncertainties(
        values=np.asarray([1.0e-18, -0.051, 0.1]),
        errors=np.asarray([0.049, 0.06, 0.331]),
        bounds=np.asarray([[0.0, 0.8], [-np.pi, np.pi], [0.1, 5.0]]),
    )

    assert np.isnan(relative[0])
    assert relative[1] == pytest.approx(0.06 / 0.051)
    assert relative[2] == pytest.approx(3.31)


# Test the parameter plot crops empty canvas without widening its content
def test_param_plot_layout_crops_empty_canvas():
    figure, axis = fit.plt.subplots(figsize=(15.5, 4.8))
    axis.set_yticks(
        [0.0],
        labels=["CON_GP|9910:[321,321]/opposite:helicity_angle1@PROJECTIVE"],
    )
    initial_width = float(figure.get_size_inches()[0])
    left_margin, right_margin = fit._compact_parameter_plot_layout(figure, axis)
    figure.subplots_adjust(left=left_margin, right=right_margin)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    label_box = axis.get_yticklabels()[0].get_window_extent(renderer)
    final_width = float(figure.get_size_inches()[0])
    content_width = axis.get_position().width * final_width

    assert final_width < initial_width
    assert content_width == pytest.approx(0.65 * initial_width)
    assert label_box.x0 > 0.01 * figure.bbox.width
    fit.plt.close(figure)


# Test long parameter summaries are split into PNG pages and one PDF
def test_param_figure_paginates_large_param(tmp_path):
    dimension = 7
    result = fit.plot_parameter_uncertainties(
        param_names=[f"parameter_{index}" for index in range(dimension)],
        best_fit=np.linspace(0.1, 0.9, dimension),
        errors=np.full(dimension, 0.05),
        bounds=np.tile([0.0, 1.0], (dimension, 1)),
        output_root=tmp_path,
        realized_best=np.linspace(0.2, 0.8, dimension),
        parameters_per_page=3,
    )

    assert result["page_count"] == 3
    assert result["parameters_per_page"] == 3
    assert [page["first_parameter"] for page in result["pages"]] == [1, 4, 7]
    assert [page["last_parameter"] for page in result["pages"]] == [3, 6, 7]
    assert len(result["png_pages"]) == 3
    assert all(Path(path).is_file() for path in result["png_pages"])
    assert Path(result["pdf"]).is_file()


# Test optimizer initialization follows the surrogate rather than raw observations
def test_surrogate_start_minimum():
    X = np.array([[0.0], [0.5], [1.0]], dtype=np.float64)
    Z = np.array([0.0, 2.0, 3.0], dtype=np.float64)

    # Place the surrogate minimum at the final sampled point
    def surrogate(values):
        values = np.asarray(values, dtype=np.float64)
        return (values[:, 0] - 1.0) ** 2

    start, diagnostics = fit.select_surrogate_start(X, Z, surrogate)

    assert np.allclose(start, [1.0])
    assert diagnostics["observed_z"] == 3.0
    assert diagnostics["surrogate_z"] == 0.0


# Test bounded L-BFGS-B evaluates every gradient direction in one batch
def test_lbfgsb_batched_finite_diff():
    dimension = 8
    target = np.linspace(0.2, 0.8, dimension)
    batch_lengths = []

    # Record vectorized calls while evaluating a convex surrogate
    def surrogate(values):
        values = np.asarray(values, dtype=np.float64)
        batch_lengths.append(len(values))
        return np.sum((values - target) ** 2, axis=1)

    result = fit.minimize_vectorized_lbfgsb(
        func_hat=surrogate,
        start=np.full(dimension, 0.95),
        bounds=np.tile([0.0, 1.0], (dimension, 1)),
        maxiter=30,
    )

    assert result.success
    assert np.allclose(result.x, target, atol=1e-6)
    assert max(batch_lengths) == 1 + 2 * dimension


# Test bounded L-BFGS-B uses an exact physical derivative callback when supplied
def test_lbfgsb_prefers_exact_gradient_callback():
    dimension = 6
    target = np.linspace(-0.7, 0.8, dimension)
    callback_calls = []
    surrogate_calls = []

    # Record any unexpected finite-difference surrogate evaluation
    def surrogate(values):
        surrogate_calls.append(len(values))
        return np.sum((np.asarray(values) - target) ** 2, axis=1)

    # Compute the exact physical objective and gradient
    def value_gradient(values):
        callback_calls.append(1)
        delta = np.asarray(values, dtype=np.float64) - target
        return float(np.sum(delta**2)), 2.0 * delta

    result = fit.minimize_vectorized_lbfgsb(
        func_hat=surrogate,
        start=np.full(dimension, 0.95),
        bounds=np.tile([-1.0, 1.0], (dimension, 1)),
        maxiter=30,
        value_gradient_func=value_gradient,
    )

    assert result.success
    assert np.allclose(result.x, target, atol=1.0e-8)
    assert callback_calls
    assert surrogate_calls == []


# Test correlation labels and annotations scale down for large fits
def test_correlation_layout_scaling():
    compact = fit.correlation_matrix_layout(["x", "y"])
    crowded = fit.correlation_matrix_layout(
        [f"RES|f0:channel-{index}:coupling" for index in range(176)]
    )

    assert crowded["figure_extent"] > compact["figure_extent"]
    assert crowded["tick_fontsize"] < compact["tick_fontsize"]
    assert crowded["tick_fontsize"] <= 3.5
    assert crowded["annotation_fontsize"] is None
    assert crowded["dpi"] < compact["dpi"]


# Test the correlation plot masks its diagonal and resolves off-diagonal structure
def test_correlation_plot_diagonal(tmp_path, monkeypatch):
    correlation = np.array(
        [
            [1.0, 0.34, -0.08],
            [0.34, 1.0, -0.21],
            [-0.08, -0.21, 1.0],
        ]
    )
    display_values, color_limit = fit.correlation_plot_values(correlation)

    assert np.all(np.isnan(np.diag(display_values)))
    assert np.allclose(display_values[~np.eye(3, dtype=bool)], correlation[~np.eye(3, dtype=bool)])
    assert color_limit == pytest.approx(0.4)

    original_close = fit.plt.close
    monkeypatch.setattr(fit.plt, "close", lambda _figure: None)
    figure_path = fit.plot_correlation_matrix(
        correlation=correlation,
        param_names=["alpha", "beta", "gamma"],
        output_root=tmp_path,
        plot_brand="GRANIITTI",
    )
    figure = fit.plt.gcf()
    axis = figure.axes[0]
    image = axis.images[0]
    image_mask = np.ma.getmaskarray(image.get_array())

    assert np.all(np.diag(image_mask))
    assert image.get_cmap().name == "RdBu_r"
    assert np.allclose(image.get_cmap().get_bad(), [0.898, 0.906, 0.922, 1.0], atol=2.0e-3)
    assert image.get_clim() == pytest.approx((-0.4, 0.4))
    assert Path(figure_path).is_file()
    original_close(figure)


# Ignore zero-weight bins before zero-variance and nonfinite arithmetic
@pytest.mark.parametrize("excluded", [9.0, np.nan])
def test_histogram_chi2_zero_weight_bins(excluded):
    counts = torch.tensor([[1.0, excluded]], dtype=torch.float64, requires_grad=True)
    errors = torch.zeros_like(counts, requires_grad=True)
    data = torch.zeros(2, dtype=counts.dtype)
    data_errors = torch.tensor([1.0, 0.0], dtype=counts.dtype)
    weights = torch.tensor([1.0, 0.0], dtype=counts.dtype)
    value = icecost.histogram_chi2_torch(
        counts=counts, errors=errors, data_counts=data, data_errors=data_errors,
        mask=torch.ones(2, dtype=torch.bool), fit_weights=weights,
    )
    numpy_value, numpy_error = icecost.global_histogram_chi2(
        data.numpy(), data_errors.numpy(), counts.detach().numpy(), errors.detach().numpy(),
        fit_weights=weights.numpy(),
    )
    assert value.item() == pytest.approx(1.0)
    np.testing.assert_allclose(numpy_value, [1.0])
    np.testing.assert_allclose(numpy_error, [0.0])
    value.sum().backward()
    torch.testing.assert_close(counts.grad, torch.tensor([[2.0, 0.0]], dtype=counts.dtype))
    assert torch.all(torch.isfinite(errors.grad))


# Remove excluded bins and their covariance rows before computing correlated chi2
@pytest.mark.parametrize("zero_weight", [False, True])
@pytest.mark.parametrize("all_excluded", [False, True])
def test_correlated_chi2_respects_selection(zero_weight, all_excluded):
    counts = torch.tensor([[1.0, 9.0]], dtype=torch.float64, requires_grad=True)
    active = torch.tensor([not all_excluded, False])
    value = icecost.histogram_chi2_torch(
        counts=counts, errors=torch.zeros_like(counts), data_counts=torch.zeros(2, dtype=counts.dtype),
        data_errors=torch.ones(2, dtype=counts.dtype),
        mask=torch.ones(2, dtype=torch.bool) if zero_weight else active,
        fit_weights=active.to(counts.dtype) if zero_weight else torch.ones(2, dtype=counts.dtype),
        data_total_covariance=torch.tensor([[1.0, 0.5], [0.5, 1.0]], dtype=counts.dtype),
        covariance_indices=torch.arange(2),
    )
    assert value.item() == pytest.approx(0.0 if all_excluded else 1.0)
    value.sum().backward()
    torch.testing.assert_close(counts.grad, torch.tensor([[0.0 if all_excluded else 2.0, 0.0]], dtype=counts.dtype))


# Reject malformed optimizer bounds before indexing their columns
@pytest.mark.parametrize("bounds", [[0.0, 1.0], [[0.0], [1.0]]])
def test_vectorized_lbfgsb_malformed_bounds(bounds):
    with pytest.raises(ValueError, match="inconsistent dimensions"):
        fit.minimize_vectorized_lbfgsb(
            func_hat=lambda values: np.sum(values**2, axis=-1), start=np.array([0.5]), bounds=bounds,
        )

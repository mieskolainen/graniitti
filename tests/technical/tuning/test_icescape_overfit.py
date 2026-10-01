# Tests for icescape surrogate overfit guards
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np
import pytest
import torch
from core.inference import gp as gp_model
from core.inference import neural
from core.inference.gp import ExactGaussianProcess


# Test automatic GP precision uses stable float64 arithmetic on every device
def test_gp_automatic_dtype_is_float64():
    assert gp_model.resolve_gp_dtype("cpu", "auto") == torch.float64
    assert gp_model.resolve_gp_dtype("cuda", "auto") == torch.float64
    assert gp_model.resolve_gp_dtype("cuda", "float32") == torch.float32


# Test NN validation uses the holdout set when uncertainty tensors are disabled
def test_neural_holdout_no_errors():
    model = torch.nn.Linear(1, 1, bias=False)
    with torch.no_grad():
        model.weight.zero_()

    optimizer = torch.optim.SGD(model.parameters(), lr=0.0)
    scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=1.0)

    _, stats = neural.train_torch(
        model=model,
        optimizer=optimizer,
        scheduler=scheduler,
        lossfunc=neural.mae_error_loss,
        X_train=np.array([[0.0]], dtype=np.float32),
        y_train=np.array([[0.0]], dtype=np.float32),
        X_val=np.array([[0.0]], dtype=np.float32),
        y_val=np.array([[10.0]], dtype=np.float32),
        errs_train=None,
        errs_val=None,
        num_epochs=5,
        batch_size=1,
        device="cpu",
        patience=2,
        min_delta=1.0,
    )

    assert stats["train_loss"][0] == pytest.approx(0.0)
    assert stats["eval_loss"][0] == pytest.approx(10.0)
    assert stats["initial_eval_loss"] == pytest.approx(10.0)
    assert stats["best_epoch"] == 0
    assert stats["stopped_epoch"] == 2


# Test GP hyperparameter optimization restores the best validation-selected state
def test_gp_best_holdout():
    X_train = torch.tensor([[0.0], [0.5], [1.0]], dtype=torch.float64)
    Y_train = torch.sin(X_train)
    X_val = torch.tensor([[1.5]], dtype=torch.float64)
    Y_val = torch.sin(X_val)

    gp = gp_model.GaussianProcess(
        X_train=X_train,
        Y_train=Y_train,
        scale=0.25,
        noise=1e-4,
        device="cpu",
    )
    stats = gp.optimize_hyperparameters(
        num_steps=5,
        lr=0.05,
        X_val=X_val,
        Y_val=Y_val,
        patience=5,
        min_delta=0.0,
    )

    assert len(stats["validation_loss"]) == len(stats["train_loss"]) > 0
    assert stats["best_validation_loss"] is not None
    assert stats["best_step"] is not None

    with torch.no_grad():
        selected_loss = gp.validation_loss(X_val, Y_val).item()

    assert selected_loss == pytest.approx(stats["best_validation_loss"], rel=1e-8, abs=1e-8)


# Test exact GP MAP fitting reports restarts and honors validation patience
def test_exact_gp_validation_patience_restarts():
    X = torch.linspace(0.0, 1.0, 5, dtype=torch.float64).reshape(-1, 1)
    settings = gp_model.exact_gp_defaults()
    settings["fit_restarts"] = 2
    gp = ExactGaussianProcess(
        X_train=X[:3],
        Y_train=torch.sin(X[:3]),
        settings=settings,
        device="cpu",
    )

    stats = gp.optimize_hyperparameters(
        num_steps=5,
        lr=1.0e-12,
        X_val=X[3:],
        Y_val=torch.sin(X[3:]),
        patience=2,
        min_delta=1.0,
        validation_interval=2,
        log_interval=1,
        restarts=2,
    )

    assert stats["restarts"] == stats["successful_restarts"] == 2
    assert len(stats["restart_stats"]) == 2
    assert all(row["stopped_step"] == 2 for row in stats["restart_stats"])
    assert stats["validation_steps"] == [0, 2, 5, 7]


# Test the compact no-uncertainty kernel matches the explicit RBF formula
def test_gp_rbf_distances():
    X1 = torch.tensor([[0.0, 1.0], [0.5, -0.2]], dtype=torch.float64)
    X2 = torch.tensor([[0.3, 0.1], [-0.4, 0.7], [0.0, 0.0]], dtype=torch.float64)
    kernel = gp_model.MultiOutputKernel(
        input_dim=2,
        output_dim=1,
        lengthscale=0.7,
        variance=1.4,
    )

    compact = kernel.base_kernel(X1, X2)
    explicit_distance = torch.sum(
        (X1[:, None, :] - X2[None, :, :]) ** 2,
        dim=-1,
    )
    expected = 1.4 * torch.exp(-0.5 * explicit_distance / 0.7**2)

    assert torch.allclose(compact, expected, rtol=1e-12, atol=1e-12)


# Test the advertised Matern kernel follows the nu = 5/2 form
def test_gp_matern_kernel_matches_explicit_formula():
    X = torch.tensor([[0.0], [0.7]], dtype=torch.float64)
    kernel = gp_model.MultiOutputKernel(
        input_dim=1,
        output_dim=1,
        lengthscale=0.4,
        variance=1.2,
        kind="Matern",
    )

    distance = torch.cdist(X, X)
    scaled = np.sqrt(5.0) * distance / 0.4
    expected = 1.2 * (1.0 + scaled + scaled**2 / 3.0) * torch.exp(-scaled)

    assert torch.allclose(kernel.base_kernel(X, X), expected, rtol=1e-12, atol=1e-12)


# Test mean-only prediction and full-data posterior replacement stay consistent
def test_gp_mean_refit():
    X = torch.linspace(0.0, 1.0, 6, dtype=torch.float64).reshape(-1, 1)
    Y = torch.sin(X)
    gp = gp_model.GaussianProcess(
        X_train=X[:4],
        Y_train=Y[:4],
        scale=0.4,
        noise=1e-3,
        device="cpu",
    )
    gp.set_training_data(X, Y)

    with torch.no_grad():
        mean_only = gp.predict_mean(X)
        full_mean, _ = gp(X)

    assert gp.N == 6
    assert torch.allclose(mean_only, full_mean, rtol=1e-12, atol=1e-12)


# Test the frozen Matern posterior supports finite first and second input derivatives
def test_gp_frozen_posterior_supports_input_autograd():
    X = torch.linspace(0.0, 1.0, 7, dtype=torch.float64).reshape(-1, 1)
    gp = gp_model.GaussianProcess(
        X_train=X,
        Y_train=torch.sin(2.0 * X),
        scale=0.5,
        noise=1.0e-3,
        device="cpu",
        kernel="Matern",
    )
    gp.freeze_for_inference()
    query = torch.tensor([0.37], dtype=torch.float64, requires_grad=True)

    # Compute one scalar posterior mean
    def objective(values):
        return gp.predict_mean(values.reshape(1, 1))[0, 0]

    value = objective(query)
    (gradient,) = torch.autograd.grad(value, query, create_graph=True)
    hessian = torch.autograd.functional.hessian(objective, query, vectorize=True)

    assert torch.all(torch.isfinite(gradient))
    assert torch.all(torch.isfinite(hessian))
    assert not any(parameter.requires_grad for parameter in gp.parameters())


# Test sparse validation checks use the requested interval
def test_gp_validation_interval():
    X = torch.linspace(0.0, 1.0, 8, dtype=torch.float64).reshape(-1, 1)
    gp = gp_model.GaussianProcess(
        X_train=X[:6],
        Y_train=torch.sin(X[:6]),
        scale=0.3,
        noise=1e-3,
        device="cpu",
    )
    stats = gp.optimize_hyperparameters(
        num_steps=5,
        lr=0.01,
        X_val=X[6:],
        Y_val=torch.sin(X[6:]),
        validation_interval=2,
        log_interval=1,
        patience=10,
    )

    assert stats["validation_steps"] == [0, 2, 4]
    assert len(stats["validation_loss"]) == 3


# Check both GP kernels retain the analytic Matern curvature at coincidence
@pytest.mark.parametrize("scalar", [False, True])
def test_matern_coincident_curvature(scalar):
    point = torch.zeros((1, 2), dtype=torch.float64, requires_grad=True)
    reference = point.detach().clone()
    if scalar:
        gp = ExactGaussianProcess(
            X_train=torch.tensor([[0.0, 0.0], [1.0, 1.0]], dtype=torch.float64),
            Y_train=torch.tensor([0.0, 1.0], dtype=torch.float64), device="cpu",
        )
        kernel = gp.kernel
        variance = gp.outputscale
        lengthscale = gp.lengthscale
    else:
        gp = gp_model.MultiOutputKernel(input_dim=2, output_dim=1, lengthscale=0.7, variance=1.3, kind="Matern")
        kernel = gp.base_kernel
        variance = gp.variance
        lengthscale = gp.lengthscale
    hessian = torch.autograd.functional.hessian(lambda x: kernel(x, reference).sum(), point)
    expected = torch.eye(2, dtype=point.dtype) * (-5.0 * variance / (3.0 * lengthscale.square()))
    torch.testing.assert_close(hessian.reshape(2, 2), expected)
    assert torch.autograd.gradgradcheck(lambda x: kernel(x, reference), (point,))


# Preserve kernel separation and input derivatives with a large shared coordinate offset
@pytest.mark.parametrize("offset", [0.0, 1.0e8, -1.0e8])
def test_gp_distance_translation_invariance(offset):
    first = torch.tensor([[offset, offset], [offset + 2.0, offset - 1.0]], dtype=torch.float64, requires_grad=True)
    second = torch.tensor([[offset + 1.0, offset]], dtype=torch.float64)
    distance = gp_model.pairwise_squared_distance(first, second)
    torch.testing.assert_close(distance, torch.tensor([[1.0], [2.0]], dtype=first.dtype))
    gradient = torch.autograd.grad(distance.sum(), first, create_graph=True)[0]
    torch.testing.assert_close(gradient, torch.tensor([[-2.0, 0.0], [2.0, -2.0]], dtype=first.dtype))
    curvature = torch.autograd.grad(gradient.sum(), first)[0]
    torch.testing.assert_close(curvature, torch.full_like(first, 2.0))


# Learn cross-output covariance instead of keeping correlated mode fixed to independent outputs
def test_multioutput_gp_learns_output_cov():
    x = torch.linspace(0.0, 1.0, 6, dtype=torch.float64)[:, None]
    y = torch.cat((torch.sin(3.0 * x), 2.0 * torch.sin(3.0 * x)), dim=1)
    model = gp_model.GaussianProcess(x, y, noise=0.01, device="cpu")
    model.optimize_hyperparameters(num_steps=3, lr=0.01)
    covariance = model.kernel_module.L @ model.kernel_module.L.T
    assert covariance[0, 1].item() > 0.0
    mean, variance = model(x)
    assert torch.isfinite(mean).all() and torch.isfinite(variance).all()


# Preserve a better initial MAP state when an aggressive optimization step overshoots
def test_exact_gp_map_fit_best_evaluated_state():
    x = torch.linspace(0.0, 1.0, 5, dtype=torch.float64)[:, None]
    settings = gp_model.exact_gp_defaults()
    settings["learning_rate"] = 100.0
    model = ExactGaussianProcess(x, torch.sin(x), settings=settings)
    initial = float(model.negative_log_posterior().detach())
    stats = model.fit(steps=1, restarts=1)
    assert stats["best_negative_log_posterior"] <= initial
    assert float(model.negative_log_posterior().detach()) == pytest.approx(stats["best_negative_log_posterior"])


# Check shared L-BFGS reaches a stationary MAP fit with known observation errors
@pytest.mark.parametrize("kernel", ["matern52", "rbf"])
def test_exact_gp_map_gradient_convergence(kernel):
    x = torch.linspace(0.0, 1.0, 16, dtype=torch.float64)[:, None]
    settings = gp_model.exact_gp_defaults()
    settings["kernel"] = kernel
    model = ExactGaussianProcess(x, torch.sin(6.0 * x), Y_err=torch.full((16,), 0.03), settings=settings, seed=19)
    stats = model.fit()
    loss = model.negative_log_posterior()
    gradients = torch.autograd.grad(loss, tuple(model.parameters()))
    assert stats["converged_restarts"] > 0
    assert max(float(value.abs().max()) for value in gradients) < 1e-3
    assert float(loss.detach()) == pytest.approx(stats["best_negative_log_posterior"])


# Preserve the best MAP state when multistart fitting has no holdout sample
def test_gp_restart_no_holdout():
    x = torch.linspace(0.0, 1.0, 4, dtype=torch.float64)[:, None]
    model = ExactGaussianProcess(x, torch.sin(x))
    initial = float(model.negative_log_posterior().detach())
    model.optimize_hyperparameters(num_steps=1, lr=10.0, restarts=1)
    assert float(model.negative_log_posterior().detach()) <= initial


# Match mean-only values and input derivatives to the full exact posterior
def test_gp_mean_derivatives():
    x = torch.linspace(0.0, 1.0, 6, dtype=torch.float64)[:, None]
    model = ExactGaussianProcess(x, torch.sin(2.0 * x))
    model.freeze_for_inference()
    points = torch.tensor([[0.13], [0.79]], dtype=torch.float64, requires_grad=True)
    mean = model.predict_mean(points).squeeze(-1)
    posterior = model.posterior(points).mean
    torch.testing.assert_close(mean, posterior)
    grad = torch.autograd.grad(mean.sum(), points)[0]
    full_grad = torch.autograd.grad(posterior.sum(), points)[0]
    torch.testing.assert_close(grad, full_grad)

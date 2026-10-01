# Test for profiled likelihood scans
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import threading

import numpy as np
import pytest
import torch
from core.inference import profile


# Test that automatic profile jobs occupy every scheduler-visible CPU
def test_profile_jobs_use_all_available_cpus():
    available = profile.available_profile_cpus()
    count = min(4, available)
    barrier = threading.Barrier(count)
    thread_ids = set()
    lock = threading.Lock()

    # Hold every task until the full worker pool is active
    def evaluate(value):
        with lock:
            thread_ids.add(threading.get_ident())
        barrier.wait(timeout=5.0)
        return value * value

    results = list(
        profile.profile_job_results(
            range(count),
            evaluate,
            workers=0,
            description="Test profiles",
        )
    )

    assert sorted(results) == [value**2 for value in range(count)]
    assert len(thread_ids) == count
    assert profile.resolve_profile_workers(12, 20) == min(12, available)
    assert profile.resolve_profile_workers(0, 2) == min(2, available)


# Test that parallel one and two dimensional profiles render complete outputs
def test_parallel_profile_plotting_writes_all_scans(tmp_path):
    X = np.asarray([[0.1, 0.2], [0.5, 0.5], [0.8, 0.9]], dtype=np.float64)

    # Evaluate one smooth correlated likelihood batch
    def objective(values):
        values = np.asarray(values, dtype=np.float64)
        return (values[:, 0] - 0.4) ** 2 + (values[:, 1] + values[:, 0] - 0.9) ** 2

    Z = objective(X)
    profile.plot_profile_likelihoods(
        X=X,
        Z=Z,
        param_names=["a", "b"],
        func_hat=objective,
        output_root=tmp_path,
        realized_best=X[int(np.argmin(Z))],
        surrogate_best=np.asarray([0.4, 0.5]),
        realized_Z=float(np.min(Z)),
        surrogate_Z=0.0,
        bounds=np.asarray([[0.0, 1.0], [0.0, 1.0]]),
        ngrid=3,
        plot_2d=True,
        maxiter=10,
        multistart=False,
        workers=2,
        dpi=80,
    )

    assert len(list((tmp_path / "1D_profile").glob("*.png"))) == 2
    assert len(list((tmp_path / "2D_profile").glob("*.png"))) == 1


# Test that profile scans re-minimize nuisance coordinates instead of slicing
def test_profile_correlated_nuisance():
    def objective(X):
        X = np.asarray(X, dtype=np.float64)
        return (X[:, 0] + X[:, 1]) ** 2 + X[:, 0] ** 2

    gradient_calls = []

    # Compute the exact physical objective gradient
    def value_gradient(values):
        gradient_calls.append(1)
        x, y = np.asarray(values, dtype=np.float64)
        return (
            float((x + y) ** 2 + x**2),
            np.asarray([2.0 * (x + y) + 2.0 * x, 2.0 * (x + y)]),
        )

    bounds = np.array([[-2.0, 2.0], [-2.0, 2.0]], dtype=np.float64)
    x_grid = np.array([-1.0, 0.0, 1.0], dtype=np.float64)
    z, profiled = profile.scan_profile_1d(
        func_hat=objective,
        index=0,
        x_grid=x_grid,
        start=np.array([0.0, 0.0]),
        bounds=bounds,
        maxiter=100,
        value_gradient_func=value_gradient,
    )

    assert np.allclose(z, [1.0, 0.0, 1.0], atol=1e-6)
    assert np.allclose(profiled[:, 1], -x_grid, atol=1e-5)
    assert gradient_calls


# Test Torch profile scans optimize nuisance coordinates with exact autograd
def test_scan_profile_1d_uses_torch_lbfgsb(monkeypatch):
    # Evaluate one smooth correlated likelihood batch
    def objective_numpy(values):
        values = np.asarray(values, dtype=np.float64)
        return (values[:, 0] + values[:, 1]) ** 2 + values[:, 0] ** 2

    # Evaluate the same likelihood without leaving the Torch graph
    def objective_torch(values):
        return (values[0] + values[1]).square() + values[0].square()

    # Evaluate all independent likelihood rows without leaving the Torch graph
    def objective_batch_torch(values):
        return (values[:, 0] + values[:, 1]).square() + values[:, 0].square()

    # Reject an unexpected scalar recovery in this smooth batched problem
    def reject_scalar_recovery(**_kwargs):
        raise AssertionError("Batched profile unexpectedly used scalar recovery")

    batch_calls = []
    original_minimize_batch = profile.lbfgsb.minimize_batch

    # Record the profile-specific fast recovery controls
    def record_minimize_batch(**kwargs):
        batch_calls.append(kwargs)
        return original_minimize_batch(**kwargs)

    monkeypatch.setattr(profile, "profile_minimize", reject_scalar_recovery)
    monkeypatch.setattr(profile.lbfgsb, "minimize_batch", record_minimize_batch)
    bounds = np.array([[-2.0, 2.0], [-2.0, 2.0]], dtype=np.float64)
    x_grid = np.array([-1.0, 0.0, 1.0], dtype=np.float64)
    z, profiled = profile.scan_profile_1d(
        func_hat=objective_numpy,
        index=0,
        x_grid=x_grid,
        start=np.array([0.0, 0.0]),
        bounds=bounds,
        maxiter=100,
        torch_config={
            "objective": objective_torch,
            "objective_batch": objective_batch_torch,
            "device": torch.device("cpu"),
            "dtype": torch.float64,
        },
    )

    assert np.allclose(z, [1.0, 0.0, 1.0], atol=1.0e-8)
    assert np.allclose(profiled[:, 1], -x_grid, atol=1.0e-7)
    assert batch_calls[0]["retry_failed_search"] is False
    assert batch_calls[0]["relative_gradient_tolerance"] == pytest.approx(1.0e-6)


# Test a two-dimensional profile grid is optimized as one Torch batch
def test_scan_profile_2d_uses_torch_batch():
    # Evaluate one three-dimensional correlated likelihood batch
    def objective_numpy(values):
        values = np.asarray(values, dtype=np.float64)
        return (
            values[:, 0] ** 2
            + values[:, 1] ** 2
            + (values[:, 2] + values[:, 0] - values[:, 1]) ** 2
        )

    # Evaluate one scalar likelihood without leaving the Torch graph
    def objective_torch(values):
        return values[0].square() + values[1].square() + (
            values[2] + values[0] - values[1]
        ).square()

    # Evaluate all independent likelihood rows without leaving the Torch graph
    def objective_batch_torch(values):
        return values[:, 0].square() + values[:, 1].square() + (
            values[:, 2] + values[:, 0] - values[:, 1]
        ).square()

    grid = np.asarray([-0.5, 0.0, 0.5], dtype=np.float64)
    z = profile.scan_profile_2d(
        func_hat=objective_numpy,
        i=0,
        j=1,
        x_grid=grid,
        y_grid=grid,
        start=np.zeros(3, dtype=np.float64),
        bounds=np.asarray([[-1.0, 1.0]] * 3),
        maxiter=100,
        multistart=False,
        torch_config={
            "objective": objective_torch,
            "objective_batch": objective_batch_torch,
            "device": torch.device("cpu"),
            "dtype": torch.float64,
        },
    )

    expected_x, expected_y = np.meshgrid(grid, grid)
    np.testing.assert_allclose(z, expected_x**2 + expected_y**2, atol=1.0e-8)


# Test profiles reject a materially inconsistent reported surrogate optimum
def test_consistency_audit_rejects_lower_profile(tmp_path):
    # Evaluate one one-dimensional likelihood with a known minimum
    def objective(values):
        values = np.asarray(values, dtype=np.float64)
        return (values[:, 0] - 0.25) ** 2

    X = np.linspace(0.0, 1.0, 5).reshape(-1, 1)
    with pytest.raises(RuntimeError, match="Profile consistency audit rejected"):
        profile.plot_profile_likelihoods(
            X=X,
            Z=objective(X),
            param_names=["coupling"],
            func_hat=objective,
            output_root=tmp_path,
            realized_best=np.asarray([0.25]),
            surrogate_best=np.asarray([0.75]),
            realized_Z=0.0,
            surrogate_Z=1.0,
            bounds=np.asarray([[0.0, 1.0]]),
            ngrid=5,
            plot_2d=False,
            maxiter=10,
            multistart=False,
            workers=1,
            dpi=80,
        )
    assert len(list((tmp_path / "1D_profile").glob("*.png"))) == 1


# Preserve a quadratic likelihood minimum and confidence interval under parameter shifts and rescaling
@pytest.mark.parametrize("offset,scale", [(0.0, 1.0), (1.0e8, 1.0), (0.0, 1.0e-12), (0.0, 1.0e12)])
def test_quadratic_coord_invariance(offset, scale):
    local = np.linspace(-2.0, 2.0, 9)
    x = offset + scale * local
    z = (local - 0.2)**2 + 3.0
    result = profile.profile_interval(x, z, delta_Z=1.0)
    assert result is not None
    minimum, z_minimum, left, right = result
    np.testing.assert_allclose((np.array([minimum, left, right]) - offset) / scale, [0.2, -0.8, 1.2], atol=1.0e-8)
    assert z_minimum == pytest.approx(3.0)

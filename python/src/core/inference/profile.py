# Shared profile-likelihood plotting utilities for icescape and iceproxy
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import os
import pathlib
import re
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager, suppress
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import torch
from scipy.optimize import minimize
from termcolor import cprint
from tqdm import tqdm

from core.inference.fit import evaluate_batch, evaluate_scalar
from core.io.files import ensure_dir
from core.io.serialize import json_safe, write_json_file
from core.numerics import lbfgsb
from core.plot.style import add_icetune_branding


# Scan decoded physical parameters with exact autograd equality constraints
def transformed_profile(*, transform, indices, grids, best_fit, bounds, value_gradient_func, maxiter, workers=0):
    reference = transform.samples([best_fit])[0]
    mesh = np.meshgrid(*grids, indexing="ij")
    targets = np.stack([values.ravel() for values in mesh], axis=1)

    # Profile one point in the physical coordinate grid
    def scan(job):
        index, target = job
        result = transform.constrain(best_fit, np.eye(len(best_fit)), bounds, list(indices),
                                     target - reference[list(indices)], maxiter=maxiter, objective=value_gradient_func)
        return index, result["value"] if result["success"] else np.nan

    values = np.full(len(targets), np.nan)
    for index, value in profile_job_results(list(enumerate(targets)), scan, workers=workers,
                                           description="Physical likelihood profile"):
        values[index] = value
    return values.reshape(tuple(len(grid) for grid in grids))


# Plot physical likelihood scans over the sampled ranges without treating decoded axes as independent bounds
def plot_transformed_profiles(*, transform, X, best_fit, bounds, value_gradient_func, output_root,
                              ngrid, maxiter, workers=0, plot_2d=False, plot_brand="GRANIITTI", dpi=300, cmap="viridis_r"):
    if value_gradient_func is None:
        raise ValueError("Physical profiles require an autograd likelihood gradient")
    reference = transform.samples([best_fit])[0]
    samples = torch.as_tensor(transform.samples(X))
    samples = reference + transform.difference(samples, torch.as_tensor(reference)).numpy()
    lower, upper = np.min(samples, axis=0), np.max(samples, axis=0)
    active = np.flatnonzero(~np.isclose(lower, upper, rtol=0.0, atol=np.finfo(float).eps))
    grids = {i: np.unique(np.r_[np.linspace(lower[i], upper[i], ngrid), reference[i]]) for i in active}
    selected = [(i,) for i in active]
    if plot_2d:
        selected.extend(combinations(active, 2))
    baseline = float(value_gradient_func(best_fit)[0])
    for indices in selected:
        axis_grids = [grids[index] for index in indices]
        values = transformed_profile(transform=transform, indices=indices, grids=axis_grids, best_fit=best_fit,
            bounds=bounds, value_gradient_func=value_gradient_func, maxiter=maxiter, workers=workers)
        output = pathlib.Path(output_root) / f"{len(indices)}D_profile"
        ensure_dir(output)
        names = [transform.names[index] for index in indices]
        stem = "__".join(safe_name(name) for name in names)
        write_json_file(output / f"{stem}.json", json_safe({"names": names, "grid": axis_grids, "objective": values,
            "reference": reference[list(indices)], "range": "sampled physical values", "baseline": baseline}), indent=4)
        fig, axis = plt.subplots()
        if len(indices) == 1:
            axis.plot(axis_grids[0], values - baseline)
            axis.axvline(reference[indices[0]], color="black", linestyle=":")
            axis.set_ylabel(r"$\Delta Z$")
        else:
            mesh = axis.pcolormesh(*axis_grids, (values - baseline).T, shading="auto", cmap=cmap)
            fig.colorbar(mesh, ax=axis, label=r"$\Delta Z$")
            axis.scatter(*reference[list(indices)], marker="x", color="red")
            axis.set_ylabel(names[1])
        axis.set_xlabel(names[0])
        add_icetune_branding(fig, brand=plot_brand)
        fig.tight_layout()
        fig.savefig(output / f"{stem}.pdf")
        fig.savefig(output / f"{stem}.png", dpi=dpi)
        plt.close(fig)


# Convert a parameter name into a filesystem-safe plot stem
def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.@+-]+", "_", str(name)).strip("_") or "param"


# Compute the CPUs available to this process under its scheduler affinity
def available_profile_cpus() -> int:
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return max(1, os.cpu_count() or 1)


# Resolve an automatic or explicit profile worker count without oversubscription
def resolve_profile_workers(workers: int, job_count: int) -> int:
    requested = int(workers)
    jobs = max(0, int(job_count))
    if requested < 0:
        raise ValueError("Profile worker count must be non-negative")
    if jobs == 0:
        return 0
    available = available_profile_cpus()
    return min(jobs, available if requested == 0 else requested, available)


# Divide native numerical threads across concurrent profile workers
@contextmanager
def profile_cpu_thread_budget(worker_count: int):
    native_threads = max(1, available_profile_cpus() // max(1, int(worker_count)))
    previous_torch_threads = None
    try:
        import torch

        previous_torch_threads = torch.get_num_threads()
        torch.set_num_threads(native_threads)
    except (ImportError, RuntimeError):
        pass

    try:
        try:
            from threadpoolctl import threadpool_limits
        except ImportError:
            yield native_threads
        else:
            with threadpool_limits(limits=native_threads):
                yield native_threads
    finally:
        if previous_torch_threads is not None:
            with suppress(RuntimeError):
                torch.set_num_threads(previous_torch_threads)


# Evaluate independent profile jobs concurrently and yield completed results
def profile_job_results(
    jobs: Sequence, evaluate: Callable, *, workers: int = 0, description: str = "Profile jobs"
) -> Iterator:
    jobs = list(jobs)
    worker_count = resolve_profile_workers(workers, len(jobs))
    if worker_count == 0:
        return
    with profile_cpu_thread_budget(worker_count) as native_threads:
        cprint(f"{description}: {worker_count} workers, {native_threads} native threads per worker", "yellow")
        if worker_count == 1:
            for job in tqdm(jobs, desc=description):
                yield evaluate(job)
            return
        with ThreadPoolExecutor(max_workers=worker_count, thread_name_prefix="profile") as executor:
            futures = [executor.submit(evaluate, job) for job in jobs]
            for future in tqdm(as_completed(futures), total=len(futures), desc=description):
                yield future.result()


# Clip a physical parameter vector into the requested bounds
def clip_to_bounds(x: np.ndarray, bounds: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64).copy()
    bounds = np.asarray(bounds, dtype=np.float64)
    return np.minimum(np.maximum(x, bounds[:, 0]), bounds[:, 1])


# Build unique starting vectors for bounded profile minimization
def unique_starts(starts: Sequence[np.ndarray], bounds: np.ndarray) -> list[np.ndarray]:
    out: list[np.ndarray] = []
    for start in starts:
        candidate = clip_to_bounds(np.asarray(start, dtype=np.float64), bounds)
        if np.all(np.isfinite(candidate)) and not any(np.allclose(candidate, old, rtol=0.0, atol=1e-12) for old in out):
            out.append(candidate)
    return out


# Build one Torch objective over the free coordinates of a fixed profile point
def torch_profile_objective(
    objective: Callable, template: torch.Tensor, free_indices: torch.Tensor, lower: torch.Tensor, widths: torch.Tensor
) -> Callable:
    # Insert free coordinates into one fixed physical profile point
    def free_objective(unit_values):
        physical = template.clone()
        physical[free_indices] = lower + unit_values * widths
        return objective(physical)

    return free_objective


# Build one Torch objective for independent fixed profile points
def torch_profile_batch_objective(
    objective: Callable, template: torch.Tensor, free_indices: torch.Tensor, lower: torch.Tensor, widths: torch.Tensor
) -> Callable:
    # Insert each row of free coordinates into its fixed physical profile point
    def free_objective(unit_values):
        physical = template.clone()
        physical[:, free_indices] = lower[None, :] + unit_values * widths[None, :]
        return objective(physical)

    return free_objective


# Minimize all non-fixed coordinates for one profile-grid point
def profile_minimize(
    func_hat: Callable,
    fixed_indices: Sequence[int],
    fixed_values: Sequence[float],
    start: np.ndarray,
    bounds: np.ndarray,
    fallback_starts: Sequence[np.ndarray] = (),
    maxiter: int = 200,
    multistart: bool = True,
    value_gradient_func: Callable | None = None,
    torch_config: dict | None = None,
    gradient_step: float = 1.0e-4,
    return_diagnostics: bool = False,
):
    bounds = np.asarray(bounds, dtype=np.float64)
    D = bounds.shape[0]
    fixed_indices = np.asarray(fixed_indices, dtype=int)
    fixed_values = np.asarray(fixed_values, dtype=np.float64)

    fixed_mask = np.zeros(D, dtype=bool)
    fixed_mask[fixed_indices] = True
    free_indices = np.where(~fixed_mask)[0]

    base = clip_to_bounds(np.asarray(start, dtype=np.float64), bounds)
    base[fixed_indices] = fixed_values
    base = clip_to_bounds(base, bounds)

    if not np.isfinite(gradient_step) or gradient_step <= 0.0:
        raise ValueError("Profile gradient step must be positive")
    if len(free_indices) == 0:
        value = evaluate_scalar(func_hat, base)
        diagnostics = {
            "success": bool(np.isfinite(value)),
            "message": "No free parameters",
            "iterations": 0,
            "evaluations": 1,
        }
        return (base, value, diagnostics) if return_diagnostics else (base, value)

    trial_starts = [base]
    if multistart:
        for fallback in fallback_starts:
            candidate = np.asarray(fallback, dtype=np.float64).copy()
            candidate[fixed_indices] = fixed_values
            trial_starts.append(candidate)
    starts = unique_starts(trial_starts, bounds)

    best_x = base.copy()
    best_z = evaluate_scalar(func_hat, best_x)

    # Compute one constrained value and its free physical-coordinate gradient
    def objective_with_gradient(free_values: np.ndarray):
        x = base.copy()
        x[free_indices] = np.asarray(free_values, dtype=np.float64)
        x[fixed_indices] = fixed_values
        if value_gradient_func is not None:
            value, gradient = value_gradient_func(x)
            return float(value), np.asarray(gradient, dtype=np.float64)[free_indices]
        widths = bounds[free_indices, 1] - bounds[free_indices, 0]
        points = np.tile(x, (1 + 2 * len(free_indices), 1))
        denominator = np.empty(len(free_indices), dtype=np.float64)
        for local, index in enumerate(free_indices):
            lower = max(bounds[index, 0], x[index] - gradient_step * widths[local])
            upper = min(bounds[index, 1], x[index] + gradient_step * widths[local])
            points[1 + local, index] = upper
            points[1 + len(free_indices) + local, index] = lower
            denominator[local] = upper - lower
        values = evaluate_batch(func_hat, points)
        if not np.all(np.isfinite(values)):
            return float("inf"), np.zeros(len(free_indices), dtype=np.float64)
        gradient = (values[1 : 1 + len(free_indices)] - values[1 + len(free_indices) :]) / denominator
        return float(values[0]), gradient

    results = []
    for full_start in starts:
        free_start = np.asarray(full_start[free_indices], dtype=np.float64)
        z_start = evaluate_scalar(func_hat, full_start)
        if z_start < best_z:
            best_z = z_start
            best_x = full_start.copy()
            best_x[fixed_indices] = fixed_values

        if torch_config is None:
            result = minimize(
                objective_with_gradient,
                free_start,
                method="L-BFGS-B",
                jac=True,
                bounds=[tuple(bounds[i]) for i in free_indices],
                options={"maxiter": int(maxiter), "ftol": 1e-8, "gtol": 1e-6},
            )
            results.append(
                {
                    "success": bool(result.success),
                    "message": str(result.message),
                    "iterations": int(getattr(result, "nit", 0)),
                    "evaluations": int(getattr(result, "nfev", 0)),
                }
            )
            free_candidate = np.asarray(result.x, dtype=np.float64)
        else:
            dtype = torch_config["dtype"]
            device = torch_config["device"]
            objective_torch = torch_config["objective"]
            lower = torch.as_tensor(bounds[free_indices, 0], dtype=dtype, device=device)
            widths = torch.as_tensor(bounds[free_indices, 1] - bounds[free_indices, 0], dtype=dtype, device=device)
            template = torch.as_tensor(base, dtype=dtype, device=device)
            free_tensor = torch.as_tensor(free_start, dtype=dtype, device=device)
            unit_start = torch.clamp((free_tensor - lower) / widths, 0.0, 1.0)
            free_tensor_indices = torch.as_tensor(free_indices, dtype=torch.long, device=device)

            unit_candidate, z_tensor, _, diagnostics = lbfgsb.minimize(
                objective=torch_profile_objective(
                    objective=objective_torch,
                    template=template,
                    free_indices=free_tensor_indices,
                    lower=lower,
                    widths=widths,
                ),
                start=unit_start,
                maxiter=maxiter,
                gradient_tolerance=float(torch_config.get("gradient_tolerance", 1.0e-5)),
                relative_gradient_tolerance=float(torch_config.get("relative_gradient_tolerance", 1.0e-6)),
                retry_failed_search=False,
            )
            results.append(
                {
                    "success": diagnostics["success"],
                    "message": diagnostics["message"],
                    "iterations": diagnostics["iterations"],
                    "evaluations": diagnostics["function_evaluations"],
                }
            )
            free_candidate = (lower + unit_candidate * widths).detach().cpu().numpy().astype(np.float64, copy=False)
            if not np.isfinite(float(z_tensor.detach().cpu())):
                continue
        if not np.all(np.isfinite(free_candidate)):
            continue

        x_candidate = base.copy()
        x_candidate[free_indices] = free_candidate
        x_candidate[fixed_indices] = fixed_values
        z_candidate = evaluate_scalar(func_hat, x_candidate)
        if z_candidate < best_z:
            best_z = z_candidate
            best_x = x_candidate

    point = clip_to_bounds(best_x, bounds)
    diagnostics = {
        "success": any(result["success"] for result in results),
        "message": "; ".join(dict.fromkeys(result["message"] for result in results)),
        "iterations": sum(result["iterations"] for result in results),
        "evaluations": sum(result["evaluations"] for result in results),
    }
    return (point, float(best_z), diagnostics) if return_diagnostics else (point, float(best_z))


# Minimize independent fixed profile points in one Torch batch
def profile_minimize_batch(
    func_hat: Callable,
    fixed_indices: Sequence[int],
    fixed_values: np.ndarray,
    start: np.ndarray,
    bounds: np.ndarray,
    fallback_starts: Sequence[np.ndarray] = (),
    maxiter: int = 200,
    multistart: bool = True,
    torch_config: dict | None = None,
) -> tuple[np.ndarray, np.ndarray, dict]:
    if torch_config is None or not callable(torch_config.get("objective_batch")):
        raise ValueError("Batched profiles require a batched Torch objective")
    bounds = np.asarray(bounds, dtype=np.float64)
    start = np.asarray(start, dtype=np.float64)
    fixed_indices = np.asarray(fixed_indices, dtype=int).reshape(-1)
    fixed_values = np.asarray(fixed_values, dtype=np.float64)
    if fixed_values.ndim == 1:
        fixed_values = fixed_values[:, None]
    if fixed_values.ndim != 2 or fixed_values.shape[1] != len(fixed_indices):
        raise ValueError("Fixed profile values have inconsistent dimensions")
    if bounds.shape != (len(start), 2):
        raise ValueError("Profile bounds have inconsistent dimensions")

    point_count = len(fixed_values)
    dimension = len(start)
    if point_count == 0:
        return (
            np.empty((0, dimension), dtype=np.float64),
            np.empty(0, dtype=np.float64),
            {"success": np.empty(0, dtype=bool), "fallback_rows": 0},
        )

    fixed_mask = np.zeros(dimension, dtype=bool)
    fixed_mask[fixed_indices] = True
    free_indices = np.where(~fixed_mask)[0]
    starts = [start]
    if multistart:
        starts.extend(fallback_starts)
    starts = unique_starts(starts, bounds)
    if not starts:
        raise ValueError("Profile minimization requires one finite starting point")

    # Replicate every start over all fixed profile rows
    candidate_templates = []
    for candidate_start in starts:
        template = np.tile(candidate_start, (point_count, 1))
        template[:, fixed_indices] = fixed_values
        candidate_templates.append(clip_to_bounds(template, bounds))
    templates = np.concatenate(candidate_templates, axis=0)
    start_count = len(starts)

    if len(free_indices) == 0:
        candidate_values = evaluate_batch(func_hat, templates).reshape(start_count, point_count)
        selected = np.argmin(candidate_values, axis=0)
        columns = np.arange(point_count)
        points = templates.reshape(start_count, point_count, dimension)[selected, columns]
        values = candidate_values[selected, columns]
        success = np.isfinite(values)
        return points, values, {"success": success, "fallback_rows": 0}

    dtype = torch_config["dtype"]
    device = torch_config["device"]
    lower = torch.as_tensor(bounds[free_indices, 0], dtype=dtype, device=device)
    widths = torch.as_tensor(bounds[free_indices, 1] - bounds[free_indices, 0], dtype=dtype, device=device)
    template_tensor = torch.as_tensor(templates, dtype=dtype, device=device)
    free_tensor = template_tensor[:, free_indices]
    unit_start = torch.clamp((free_tensor - lower[None, :]) / widths[None, :], 0.0, 1.0)
    free_tensor_indices = torch.as_tensor(free_indices, dtype=torch.long, device=device)
    unit_points, objective_values, _, batch_diagnostics = lbfgsb.minimize_batch(
        objective=torch_profile_batch_objective(
            objective=torch_config["objective_batch"],
            template=template_tensor,
            free_indices=free_tensor_indices,
            lower=lower,
            widths=widths,
        ),
        start=unit_start,
        maxiter=maxiter,
        gradient_tolerance=float(torch_config.get("gradient_tolerance", 1.0e-5)),
        relative_gradient_tolerance=float(torch_config.get("relative_gradient_tolerance", 1.0e-6)),
        retry_failed_search=False,
    )
    physical_points = template_tensor.clone()
    physical_points[:, free_tensor_indices] = lower[None, :] + unit_points * widths[None, :]
    candidate_points = (
        physical_points.detach()
        .cpu()
        .numpy()
        .astype(np.float64, copy=False)
        .reshape(start_count, point_count, dimension)
    )
    candidate_values = (
        objective_values.detach().cpu().numpy().astype(np.float64, copy=False).reshape(start_count, point_count)
    )
    candidate_success = batch_diagnostics["success"].detach().cpu().numpy().reshape(start_count, point_count)
    candidate_iterations = batch_diagnostics["iterations"].detach().cpu().numpy().reshape(start_count, point_count)
    finite = np.isfinite(candidate_values)
    converged_values = np.where(finite & candidate_success, candidate_values, np.inf)
    finite_values = np.where(finite, candidate_values, np.inf)
    converged_selected = np.argmin(converged_values, axis=0)
    finite_selected = np.argmin(finite_values, axis=0)
    has_converged = np.any(finite & candidate_success, axis=0)
    selected = np.where(has_converged, converged_selected, finite_selected)
    columns = np.arange(point_count)
    points = candidate_points[selected, columns].copy()
    values = candidate_values[selected, columns].copy()
    success = candidate_success[selected, columns].copy()
    selected_iterations = candidate_iterations[selected, columns]

    # Recover only rows whose selected batched solve did not converge
    recovery_rows = np.where(~success | ~np.isfinite(values))[0]
    for row in recovery_rows:
        recovery_start = points[row] if np.all(np.isfinite(points[row])) else start
        point, value, diagnostics = profile_minimize(
            func_hat=func_hat,
            fixed_indices=fixed_indices,
            fixed_values=fixed_values[row],
            start=recovery_start,
            bounds=bounds,
            fallback_starts=(),
            maxiter=max(5, int(maxiter) - int(selected_iterations[row])),
            multistart=False,
            torch_config=torch_config,
            return_diagnostics=True,
        )
        points[row] = point
        values[row] = value
        success[row] = diagnostics["success"] and np.isfinite(value)

    diagnostics = {
        "success": success,
        "fallback_rows": len(recovery_rows),
        "batch_evaluations": batch_diagnostics["batch_evaluations"],
        "row_evaluations": batch_diagnostics["row_evaluations"],
    }
    return points, values, diagnostics


# Scan a true 1D profile by re-minimizing all other coordinates
def scan_profile_1d(
    func_hat: Callable,
    index: int,
    x_grid: np.ndarray,
    start: np.ndarray,
    bounds: np.ndarray,
    fallback_starts: Sequence[np.ndarray] = (),
    maxiter: int = 200,
    multistart: bool = True,
    value_gradient_func: Callable | None = None,
    torch_config: dict | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    x_grid = np.asarray(x_grid, dtype=np.float64)
    if torch_config is not None and callable(torch_config.get("objective_batch")):
        profiled, values, _ = profile_minimize_batch(
            func_hat=func_hat,
            fixed_indices=[index],
            fixed_values=x_grid[:, None],
            start=start,
            bounds=bounds,
            fallback_starts=fallback_starts,
            maxiter=maxiter,
            multistart=multistart,
            torch_config=torch_config,
        )
        return values, profiled

    order = np.argsort(x_grid)
    x_sorted = x_grid[order]
    start = np.asarray(start, dtype=np.float64)
    pivot = int(np.argmin(np.abs(x_sorted - start[index])))
    sorted_profiled = np.full((len(x_sorted), len(start)), np.nan, dtype=np.float64)
    sorted_z = np.full(len(x_sorted), np.nan, dtype=np.float64)

    # Fit one grid point and retain its nuisance coordinates as a warm start
    def fit_grid_point(position: int, current: np.ndarray, use_multistart: bool):
        return profile_minimize(
            func_hat=func_hat,
            fixed_indices=[index],
            fixed_values=[float(x_sorted[position])],
            start=current,
            bounds=bounds,
            fallback_starts=fallback_starts,
            maxiter=maxiter,
            multistart=use_multistart,
            value_gradient_func=value_gradient_func,
            torch_config=torch_config,
        )

    center, center_z = fit_grid_point(pivot, start, multistart)
    sorted_profiled[pivot] = center
    sorted_z[pivot] = center_z
    current = center
    for position in range(pivot - 1, -1, -1):
        current, value = fit_grid_point(position, current, False)
        sorted_profiled[position] = current
        sorted_z[position] = value
    current = center
    for position in range(pivot + 1, len(x_sorted)):
        current, value = fit_grid_point(position, current, False)
        sorted_profiled[position] = current
        sorted_z[position] = value

    inverse_order = np.argsort(order)
    return sorted_z[inverse_order], sorted_profiled[inverse_order]


# Scan a true 2D profile by re-minimizing all other coordinates
def scan_profile_2d(
    func_hat: Callable,
    i: int,
    j: int,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    start: np.ndarray,
    bounds: np.ndarray,
    fallback_starts: Sequence[np.ndarray] = (),
    maxiter: int = 200,
    multistart: bool = True,
    value_gradient_func: Callable | None = None,
    torch_config: dict | None = None,
) -> np.ndarray:
    x_grid = np.asarray(x_grid, dtype=np.float64)
    y_grid = np.asarray(y_grid, dtype=np.float64)
    if torch_config is not None and callable(torch_config.get("objective_batch")):
        x_mesh, y_mesh = np.meshgrid(x_grid, y_grid)
        fixed_values = np.column_stack((x_mesh.reshape(-1), y_mesh.reshape(-1)))
        _, values, _ = profile_minimize_batch(
            func_hat=func_hat,
            fixed_indices=[i, j],
            fixed_values=fixed_values,
            start=start,
            bounds=bounds,
            fallback_starts=fallback_starts,
            maxiter=maxiter,
            multistart=multistart,
            torch_config=torch_config,
        )
        return values.reshape(len(y_grid), len(x_grid))

    z_grid = np.full((len(y_grid), len(x_grid)), np.nan, dtype=np.float64)
    current = np.asarray(start, dtype=np.float64)

    for iy, y_value in enumerate(y_grid):
        row_start = current.copy()
        for ix, x_value in enumerate(x_grid):
            row_start, z = profile_minimize(
                func_hat=func_hat,
                fixed_indices=[i, j],
                fixed_values=[float(x_value), float(y_value)],
                start=row_start,
                bounds=bounds,
                fallback_starts=fallback_starts,
                maxiter=maxiter,
                multistart=multistart,
                value_gradient_func=value_gradient_func,
                torch_config=torch_config,
            )
            z_grid[iy, ix] = z
        current = row_start.copy()

    return z_grid


# Estimate a smooth local profile minimum from neighboring grid points
def profile_quadratic_minimum(x: np.ndarray, z: np.ndarray, i0: int) -> tuple[float, float]:
    if i0 <= 0 or i0 >= len(x) - 1:
        return float(x[i0]), float(z[i0])

    xs = x[i0 - 1 : i0 + 2]
    zs = z[i0 - 1 : i0 + 2]
    if len(np.unique(xs)) < 3:
        return float(x[i0]), float(z[i0])

    origin = float(x[i0])
    scale = float(np.max(np.abs(xs - origin)))
    xs = (xs - origin) / scale
    baseline = float(z[i0])
    try:
        a, b, c = np.polyfit(xs, zs - baseline, 2)
    except Exception:
        return float(x[i0]), float(z[i0])

    if not np.all(np.isfinite([a, b, c])) or a <= 0.0:
        return float(x[i0]), float(z[i0])

    x_hat = -b / (2.0 * a)
    if x_hat < xs[0] or x_hat > xs[-1]:
        return float(x[i0]), float(z[i0])

    z_hat = a * x_hat * x_hat + b * x_hat + c + baseline
    if not np.isfinite(z_hat):
        return float(x[i0]), float(z[i0])

    return float(origin + scale * x_hat), float(z_hat)


# Estimate a threshold crossing using a bracketed quadratic fit
def profile_quadratic_crossing(x: np.ndarray, z: np.ndarray, z_level: float, i_inside: int, i_outside: int) -> float:
    za = z[i_inside]
    zb = z[i_outside]
    if zb == za:
        return float(x[i_inside])

    t = np.clip((z_level - za) / (zb - za), 0.0, 1.0)
    x_linear = float(x[i_inside] + t * (x[i_outside] - x[i_inside]))

    n = len(x)
    if i_inside < i_outside:
        indices = [max(0, i_inside - 1), i_inside, i_outside]
        if len(set(indices)) < 3:
            indices = [i_inside, i_outside, min(n - 1, i_outside + 1)]
    else:
        indices = [i_outside, i_inside, min(n - 1, i_inside + 1)]
        if len(set(indices)) < 3:
            indices = [max(0, i_outside - 1), i_outside, i_inside]

    indices = sorted(set(indices))
    if len(indices) < 3:
        return x_linear

    xs = x[indices]
    zs = z[indices]
    if len(np.unique(xs)) < 3:
        return x_linear

    origin = float(x[i_inside])
    scale = float(np.max(np.abs(xs - origin)))
    xs = (xs - origin) / scale
    zs = zs - z_level
    z_scale = float(np.max(np.abs(zs)))
    if not z_scale > 0.0:
        return x_linear
    try:
        a, b, c = np.polyfit(xs, zs / z_scale, 2)
    except Exception:
        return x_linear

    if not np.all(np.isfinite([a, b, c])):
        return x_linear

    if abs(a) < 1e-14:
        if abs(b) < 1e-14:
            return x_linear
        root = origin - scale * c / b
        lo = min(x[i_inside], x[i_outside])
        hi = max(x[i_inside], x[i_outside])
        return float(root) if lo <= root <= hi else x_linear

    disc = b * b - 4.0 * a * c
    if disc < 0.0:
        return x_linear

    roots = [origin + scale * r for r in ((-b - np.sqrt(disc)) / (2.0 * a), (-b + np.sqrt(disc)) / (2.0 * a))]
    lo = min(x[i_inside], x[i_outside])
    hi = max(x[i_inside], x[i_outside])
    bracket_roots = [r for r in roots if np.isfinite(r) and lo <= r <= hi]
    if not bracket_roots:
        return x_linear

    return float(min(bracket_roots, key=lambda r: abs(r - x_linear)))


# Compute a likelihood interval from a scanned profile curve
def profile_interval(x: np.ndarray, z: np.ndarray, delta_Z: float) -> tuple[float, float, float, float] | None:
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    z = np.asarray(z, dtype=np.float64).reshape(-1)
    finite = np.isfinite(x) & np.isfinite(z)
    x = x[finite]
    z = z[finite]
    if len(x) < 2:
        return None
    order = np.argsort(x)
    x = x[order]
    z = z[order]

    i0 = int(np.argmin(z))
    x_hat, z_hat = profile_quadratic_minimum(x, z, i0)
    z_level = float(z_hat + delta_Z)
    inside = z <= z_level

    left_i = i0
    while left_i > 0 and inside[left_i - 1]:
        left_i -= 1
    right_i = i0
    while right_i < len(z) - 1 and inside[right_i + 1]:
        right_i += 1

    x_left = float(x[0]) if left_i == 0 else profile_quadratic_crossing(x, z, z_level, left_i, left_i - 1)
    x_right = float(x[-1]) if right_i == len(z) - 1 else profile_quadratic_crossing(x, z, z_level, right_i, right_i + 1)
    return float(x_hat), float(z_hat), min(x_left, x_right), max(x_left, x_right)


# Format finite values with the requested number of significant figures
def format_sigfig(value: float, sigfig: int = 2) -> str:
    if not np.isfinite(value):
        return "nan"
    if value == 0.0:
        return "0"
    return f"{value:.{sigfig}g}"


# Compute decimal places implied by a significant-figure uncertainty
def uncertainty_decimals(value: float, sigfig: int = 2) -> int | None:
    if not np.isfinite(value) or value <= 0.0:
        return None
    exponent = int(np.floor(np.log10(abs(value))))
    return max(0, sigfig - 1 - exponent)


# Format a central value with precision set by asymmetric uncertainties
def format_value_with_errors(value: float, error_pos: float, error_neg: float, sigfig: int = 2) -> str:
    decimals = [
        d for d in [uncertainty_decimals(error_pos, sigfig), uncertainty_decimals(error_neg, sigfig)] if d is not None
    ]
    if not np.isfinite(value):
        return "nan"
    if not decimals:
        return f"{value:.4g}"
    return f"{value:.{max(decimals)}f}"


# Build a compact legend entry for one profile interval
def profile_interval_label(name: str, interval: tuple[float, float, float, float] | None) -> str:
    if interval is None:
        return f"{name}: unavailable"
    x0, _, x_left, x_right = interval
    error_neg = x0 - x_left
    error_pos = x_right - x0
    rel_neg = 100.0 * error_neg / abs(x0) if x0 != 0.0 else np.nan
    rel_pos = 100.0 * error_pos / abs(x0) if x0 != 0.0 else np.nan

    value_text = format_value_with_errors(x0, error_pos, error_neg, sigfig=2)
    error_pos_text = format_sigfig(error_pos, sigfig=2)
    error_neg_text = format_sigfig(error_neg, sigfig=2)
    rel_pos_text = format_sigfig(rel_pos, sigfig=2)
    rel_neg_text = format_sigfig(rel_neg, sigfig=2)

    return f"{name}: {value_text} $\\pm$ (+{error_pos_text}, -{error_neg_text}), (+{rel_pos_text}, -{rel_neg_text}) %"


# Draw the realized and surrogate optima on one two-dimensional profile axis
def plot_best_markers(ax, i: int, j: int, markers: Sequence[tuple[np.ndarray, str, str, float]]) -> None:
    for point, color, name, objective in markers:
        ax.scatter(point[i], point[j], color=color, s=60, edgecolors="k", label=f"{name} best (Z: {objective:.2f})")


# Plot one- and two-dimensional profile-likelihood scans
def plot_profile_likelihoods(
    X: np.ndarray,
    Z: np.ndarray,
    param_names: Sequence[str],
    func_hat: Callable,
    output_root: pathlib.Path,
    realized_best: np.ndarray,
    surrogate_best: np.ndarray,
    realized_Z: float,
    surrogate_Z: float,
    bounds: np.ndarray,
    ngrid: int = 80,
    colorbar_title: str = "$Z = -2 \\log L(\\theta)$",
    cmap: str = "viridis_r",
    vmax: float | None = None,
    plot_2d: bool = False,
    maxiter: int = 200,
    multistart: bool = True,
    value_gradient_func: Callable | None = None,
    torch_config: dict | None = None,
    workers: int = 0,
    dpi: int = 300,
    plot_brand: str = "GRANIITTI",
) -> None:
    X = np.asarray(X, dtype=np.float64)
    Z = np.asarray(Z, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)
    realized_best = clip_to_bounds(np.asarray(realized_best, dtype=np.float64), bounds)
    surrogate_best = clip_to_bounds(np.asarray(surrogate_best, dtype=np.float64), bounds)

    vmin = float(np.nanmin(Z))
    if vmax is None:
        vmax = vmin + 2.0
    if vmax < vmin:
        raise ValueError(f"plot_profile_likelihoods: vmax = {vmax} < min(Z) = {vmin}")

    D = X.shape[1]
    output_root = pathlib.Path(output_root)
    fallback_starts = [surrogate_best, realized_best]
    if torch_config is not None and str(torch_config["device"]).startswith("cuda"):
        workers = 1
    if torch_config is not None and callable(torch_config.get("objective_batch")):
        cprint("Torch profile batching enabled for independent grid points and starts", "yellow")
    best_markers = ((realized_best, "g", "Realized", realized_Z), (surrogate_best, "r", "Surrogate", surrogate_Z))

    cprint("Plotting 1D profile-likelihood scans ...", "yellow")
    fullpath_1d = output_root / "1D_profile"
    ensure_dir(fullpath_1d)
    print(f"Saving 1D profile figures to {fullpath_1d}")

    # Compute one independent profile scan for each physical parameter
    def scan_1d_job(index: int):
        x_range = bounds[index]
        x_fine = np.linspace(float(x_range[0]), float(x_range[1]), int(ngrid))
        x_fine = np.unique(np.sort(np.concatenate((x_fine, [surrogate_best[index], realized_best[index]]))))

        z_profile, _ = scan_profile_1d(
            func_hat=func_hat,
            index=index,
            x_grid=x_fine,
            start=surrogate_best,
            bounds=bounds,
            fallback_starts=fallback_starts,
            maxiter=maxiter,
            multistart=multistart,
            value_gradient_func=value_gradient_func,
            torch_config=torch_config,
        )
        return index, x_fine, z_profile

    scans_1d = profile_job_results(range(D), scan_1d_job, workers=workers, description="1D profile scans")
    profile_failures = []
    for i, x_fine, z_profile in scans_1d:
        x_range = bounds[i]
        finite_indices = np.flatnonzero(np.isfinite(z_profile))
        if finite_indices.size:
            minimum_index = int(finite_indices[np.argmin(z_profile[finite_indices])])
            profile_minimum = float(z_profile[minimum_index])
            profile_parameter_value = float(x_fine[minimum_index])
        else:
            profile_minimum = np.inf
            profile_parameter_value = np.nan
        consistency_tolerance = max(1.0e-5, 1.0e-6 * abs(float(surrogate_Z)))
        if float(surrogate_Z) > profile_minimum + consistency_tolerance:
            profile_failures.append(
                {
                    "parameter": str(param_names[i]),
                    "parameter_value": profile_parameter_value,
                    "profile_z": profile_minimum,
                    "surrogate_z": float(surrogate_Z),
                    "difference": float(surrogate_Z - profile_minimum),
                }
            )
            cprint(
                "ERROR: profile scan found a lower surrogate minimum for "
                f"{param_names[i]} at {param_names[i]} = {profile_parameter_value:.9g}: "
                f"profile Z = {profile_minimum:.9g}, "
                f"reported Z = {float(surrogate_Z):.9g}, "
                f"Delta Z = {surrogate_Z - profile_minimum:.9g}",
                "red",
                attrs=["bold"],
            )
        interval_68 = profile_interval(x_fine, z_profile, delta_Z=1.0)
        interval_95 = profile_interval(x_fine, z_profile, delta_Z=3.84)
        z_profile_min = float(interval_68[1]) if interval_68 is not None else float(np.nanmin(z_profile))

        delta_profile = z_profile - z_profile_min
        fig, ax = plt.subplots(figsize=(6.4, 5.6))
        if interval_95 is not None:
            ax.axvspan(
                interval_95[2],
                interval_95[3],
                color="#eeeeee",
                label=profile_interval_label("95% CL", interval_95),
                zorder=0,
            )
        if interval_68 is not None:
            ax.axvspan(
                interval_68[2],
                interval_68[3],
                color="#d3d3d3",
                label=profile_interval_label("68% CL", interval_68),
                zorder=1,
            )
        ax.plot(x_fine, delta_profile, color="black", linewidth=1.6, label="Profile likelihood", zorder=4)
        ax.axvline(
            surrogate_best[i],
            color="#e42536",
            linestyle="--",
            linewidth=1.0,
            label=f"Surrogate optimum: {surrogate_best[i]:.3g}",
            zorder=3,
        )
        ax.axvline(
            realized_best[i],
            color="#2e8b57",
            linestyle=":",
            linewidth=1.0,
            label=f"Best trial: {realized_best[i]:.3g}",
            zorder=3,
        )
        ax.set_xlabel(f"{param_names[i]}")
        ax.set_ylabel(r"$\Delta Z = Z - Z_{\min}$")
        ax.set_xlim(float(x_range[0]), float(x_range[1]))
        finite_delta = delta_profile[np.isfinite(delta_profile)]
        upper = max(4.5, 1.05 * float(np.max(finite_delta))) if finite_delta.size else 4.5
        ax.set_ylim(0.0, upper)
        ax.margins(x=0)
        ax.grid(axis="y", color="#e5e5e5", linewidth=0.6)
        ax.legend(
            loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, borderaxespad=0.0, frameon=False, fontsize=7.5
        )
        fig.subplots_adjust(top=0.91, bottom=0.27, left=0.15, right=0.98)
        add_icetune_branding(fig, brand=plot_brand, x=0.015, y=0.985)
        fig.savefig(fullpath_1d / f"hmesh__{safe_name(param_names[i])}.png", bbox_inches="tight", dpi=int(dpi))
        plt.close(fig)

    if profile_failures:
        failures = ", ".join(
            f"{item['parameter']} (value = {item['parameter_value']:.9g}, "
            f"profile Z = {item['profile_z']:.9g}, "
            f"reported Z = {item['surrogate_z']:.9g}, "
            f"Delta Z = {item['difference']:.9g})"
            for item in profile_failures
        )
        raise RuntimeError("Profile consistency audit rejected the reported surrogate optimum: " + failures)

    if not plot_2d:
        cprint("Skipping 2D profile-likelihood scans; use --plot-2d to enable them", "yellow")
        return
    if D < 2:
        cprint("Skipping 2D profile-likelihood scans because D < 2", "yellow")
        return

    cprint("Plotting 2D profile-likelihood scans ...", "yellow")
    fullpath_2d = output_root / "2D_profile"
    ensure_dir(fullpath_2d)
    print(f"Saving 2D profile figures to {fullpath_2d}")

    pairs = [(i, j) for i in range(D - 1) for j in range(i + 1, D)]

    # Compute one independent profile surface for each parameter pair
    def scan_2d_job(pair: tuple[int, int]):
        i, j = pair
        x_range = bounds[i]
        y_range = bounds[j]
        x_fine = np.linspace(float(x_range[0]), float(x_range[1]), int(ngrid))
        y_fine = np.linspace(float(y_range[0]), float(y_range[1]), int(ngrid))
        z_grid = scan_profile_2d(
            func_hat=func_hat,
            i=i,
            j=j,
            x_grid=x_fine,
            y_grid=y_fine,
            start=surrogate_best,
            bounds=bounds,
            fallback_starts=fallback_starts,
            maxiter=maxiter,
            multistart=multistart,
            value_gradient_func=value_gradient_func,
            torch_config=torch_config,
        )
        return i, j, x_fine, y_fine, z_grid

    scans_2d = profile_job_results(pairs, scan_2d_job, workers=workers, description="2D profile scans")
    for i, j, x_fine, y_fine, z_grid in scans_2d:
        p = X[:, [i, j]]
        x_range = bounds[i]
        y_range = bounds[j]
        x_grid, y_grid = np.meshgrid(x_fine, y_fine)

        fig, ax = plt.subplots(1, 2, figsize=(12, 6))
        sc = ax[0].scatter(p[:, 0], p[:, 1], c=Z, s=20, ec="k", cmap=cmap, vmin=vmin, vmax=vmax)
        fig.colorbar(sc, ax=ax[0], label=colorbar_title)
        ax[0].set_facecolor(np.zeros(3))
        plot_best_markers(ax[0], i, j, best_markers)
        ax[0].set_title("All trial points", fontsize=7)
        ax[0].set_xlabel(f"{param_names[i]}")
        ax[0].set_ylabel(f"{param_names[j]}")
        ax[0].legend(loc="upper center", fontsize=10)

        z_min = float(np.nanmin(z_grid))
        levels = [z_min + 2.30, z_min + 5.99]
        if np.nanmax(z_grid) > levels[0]:
            CS = ax[1].contour(x_grid, y_grid, z_grid, levels=levels, colors="k", linestyles=["-", "--"])
            ax[1].clabel(CS, inline=True, fontsize=8, fmt={levels[0]: "CL68", levels[1]: "CL95"})
        dy = y_range[1] - y_range[0]
        aspect_ratio = (x_range[1] - x_range[0]) / dy if dy != 0 else 1.0
        ax[1].set_aspect(aspect_ratio)
        ax[1].set_title("Profile likelihood contours", fontsize=7)
        ax[1].set_xlabel(f"{param_names[i]}")
        ax[1].set_ylabel(f"{param_names[j]}")
        plot_best_markers(ax[1], i, j, best_markers)
        ax[1].legend(loc="upper center", fontsize=10)
        add_icetune_branding(fig, brand=plot_brand, x=0.01, y=0.995)

        name_i = safe_name(param_names[i])
        name_j = safe_name(param_names[j])
        fig.savefig(fullpath_2d / f"hmesh__{name_i}--{name_j}.png", bbox_inches="tight", dpi=int(dpi))
        plt.close(fig)

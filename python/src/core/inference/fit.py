# Shared torch surrogate optimization and uncertainty reporting
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import pathlib
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import torch
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea
from scipy.optimize import minimize
from tabulate import tabulate
from termcolor import colored, cprint
from tqdm import tqdm

from core.io.files import ensure_dir
from core.io.serialize import json_safe, write_json_file
from core.numerics import lbfgsb
from core.plot.style import add_icetune_branding
from core.stats.transform import covariance_summary
from core.tune import summary as fit_summary
from core.tune.parameters.space import from_unit, to_unit


@dataclass(frozen=True)
class StyledCell:
    """Store one table value with optional terminal styling."""

    value: Any
    color: str | None = None
    attrs: tuple[str, ...] = ()


# Compute one printable table cell
def _table_cell(value: Any) -> str:
    if isinstance(value, StyledCell):
        value = value.value
    if value is None:
        return "N/A"
    if isinstance(value, (float, np.floating)):
        if np.isposinf(value):
            return "inf"
        if np.isneginf(value):
            return "-inf"
        if not np.isfinite(value):
            return "N/A"
        return f"{value:.6g}"
    return str(value)


# Render one cell with optional terminal styling
def _render_table_cell(value: Any) -> str:
    text = _table_cell(value)
    if isinstance(value, StyledCell) and value.color is not None:
        return colored(text, value.color, attrs=list(value.attrs))
    return text


# Compute one ASCII table as a string
def format_table(headers: Sequence[str], rows: Sequence[Sequence[Any]], align: Sequence[str] | None = None) -> str:
    header_cells = [_render_table_cell(header) for header in headers]
    if not header_cells:
        return ""
    ncols = len(header_cells)
    table_rows = [
        [_render_table_cell(cell) for cell in list(row)[:ncols]] + [""] * max(0, ncols - len(row)) for row in rows
    ]
    requested_align = list(align or ())
    column_align = tuple(requested_align[index] if index < len(requested_align) else "left" for index in range(ncols))
    return tabulate(table_rows, headers=header_cells, tablefmt="presto", colalign=column_align, disable_numparse=True)


# Print one titled ASCII table with optional colored title
def print_table(
    title: str | None,
    headers: Sequence[str],
    rows: Sequence[Sequence[Any]],
    align: Sequence[str] | None = None,
    color: str = "yellow",
) -> None:
    if title:
        cprint(title, color)
    table = format_table(headers=headers, rows=rows, align=align)
    if table:
        print(table)
    print()


# Print one shared GP validation model-selection summary
def print_gp_validation_summary(stats: dict, title: str = "GP validation summary:") -> None:
    rows = [["optimized steps", len(stats.get("train_loss", []))]]
    if "restarts" in stats:
        rows.append(["restarts", stats.get("restarts")])
    if "successful_restarts" in stats:
        rows.append(["successful restarts", stats.get("successful_restarts")])
    if "validation_loss" in stats:
        rows.append(["validation evaluations", len(stats.get("validation_loss", []))])
    rows.extend(
        [
            ["last train loss", stats.get("train_loss", [None])[-1] if stats.get("train_loss") else None],
            ["initial validation NLPD", stats.get("initial_validation_loss")],
            ["best validation NLPD", stats.get("best_validation_loss")],
            ["best step", "initial" if stats.get("best_step") == -1 else stats.get("best_step")],
            ["best restart", stats.get("best_restart")],
            ["stopped step", stats.get("stopped_step")],
            ["stopped restart", stats.get("stopped_restart")],
        ]
    )
    for key, label in (("elapsed_seconds", "elapsed seconds"), ("device", "device"), ("dtype", "dtype")):
        if key in stats:
            rows.append([label, stats.get(key)])
    print_table(title, ["Quantity", "Value"], rows)


# Print one shared neural validation model-selection summary
def print_neural_validation_summary(stats: dict, title: str = "Neural validation summary:") -> None:
    print_table(
        title,
        ["Quantity", "Value"],
        [
            ["trained epochs", len(stats.get("train_loss", []))],
            ["last train loss", stats.get("train_loss", [None])[-1] if stats.get("train_loss") else None],
            ["initial validation loss", stats.get("initial_eval_loss")],
            ["last validation loss", stats.get("eval_loss", [None])[-1] if stats.get("eval_loss") else None],
            ["best validation loss", stats.get("best_val_loss")],
            ["best epoch", "initial" if stats.get("best_epoch") == 0 else stats.get("best_epoch")],
            ["stopped epoch", stats.get("stopped_epoch")],
        ],
    )


# Configure CUDA arithmetic for high-throughput surrogate work
def configure_torch_runtime(device: str) -> None:
    if not str(device).startswith("cuda"):
        return
    cuda_device = torch.device(device)
    if cuda_device.index is not None:
        torch.cuda.set_device(cuda_device)
    torch.set_float32_matmul_precision("high")
    torch.backends.cuda.matmul.allow_tf32 = True


# Split rows while retaining required likelihood points in training
def training_validation_indices(
    n_rows: int, validation_fraction: float, rngseed: int, required_training: Sequence[int]
) -> tuple[np.ndarray, np.ndarray]:
    if n_rows < 3:
        raise ValueError("At least three rows are required for a validation split")
    if not 0.0 < float(validation_fraction) < 1.0:
        raise ValueError("Validation fraction must be in (0, 1)")
    required = np.unique(np.asarray(required_training, dtype=np.int64))
    if np.any(required < 0) or np.any(required >= n_rows):
        raise ValueError("Required training row is outside the available range")
    available = np.setdiff1d(np.arange(n_rows, dtype=np.int64), required)
    n_validation = max(1, int(round(float(validation_fraction) * n_rows)))
    n_validation = min(n_validation, len(available), n_rows - 2)
    rng = np.random.default_rng(int(rngseed))
    validation = np.sort(rng.permutation(available)[:n_validation])
    training_mask = np.ones(n_rows, dtype=bool)
    training_mask[validation] = False
    return np.flatnonzero(training_mask), validation


# Standardize columns using all rows or one selected training subset
def standardize(
    values: np.ndarray, rows: np.ndarray | Sequence[int] | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=np.float64)
    reference = values if rows is None else values[np.asarray(rows, dtype=np.int64)]
    location = np.mean(reference, axis=0)
    scale = np.std(reference, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 0.0), scale, 1.0)
    return (values - location) / scale, location, scale


# Build physical parameter bounds from metadata or sampled coordinates
def parameter_bounds(X: np.ndarray, parameter_space: Sequence[dict] | None = None) -> np.ndarray:
    values = np.asarray(X, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] == 0:
        raise ValueError("Parameter coordinates must be a non-empty matrix")
    if parameter_space:
        bounds = np.asarray(
            [[float(item["lower"]), float(item["upper"])] for item in parameter_space], dtype=np.float64
        )
    else:
        bounds = np.vstack((np.min(values, axis=0), np.max(values, axis=0))).T
    if bounds.shape != (values.shape[1], 2):
        raise ValueError("Parameter bounds have inconsistent dimensions")
    if np.any(~np.isfinite(bounds)) or np.any(bounds[:, 1] <= bounds[:, 0]):
        raise ValueError("Parameter bounds must be finite and non-degenerate")
    return bounds


# Apply one shared likelihood-quality cut and retain the realized minimum
def quality_mask(
    Z: np.ndarray, *, mode: str, max_delta_Z: float, quality_percentile: float, max_Z: float
) -> tuple[np.ndarray, float]:
    values = np.asarray(Z, dtype=np.float64)
    if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("Quality cut requires one finite non-empty objective vector")
    z_best = float(np.min(values))
    if mode == "none":
        threshold = np.inf
        keep = np.ones(len(values), dtype=bool)
    elif mode == "delta":
        threshold = z_best + float(max_delta_Z)
        keep = values <= threshold
    elif mode == "percentile":
        threshold = float(np.percentile(values, float(quality_percentile)))
        keep = values <= threshold
    elif mode == "absolute":
        threshold = float(max_Z)
        keep = values <= threshold
    else:
        raise ValueError(f'Unknown quality mode "{mode}"')
    keep[int(np.argmin(values))] = True
    return keep, threshold


# Compute one finite float formatted for user-facing fit tables
def _format_float(value: Any, digits: int = 3) -> str:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return "N/A"
    if not np.isfinite(out):
        return "N/A"
    return f"{out:.{digits}f}"


# Format one physical parameter with at most three decimal places
def _format_parameter_value(value: Any) -> str:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return "N/A"
    if not np.isfinite(out):
        return "N/A"
    text = f"{out:.3f}".rstrip("0").rstrip(".")
    return text if "." in text else f"{text}.0"


# Format one relative uncertainty percentage without redundant trailing zeroes
def _format_relative_percent(relative: float) -> str:
    text = f"{float(relative) * 100:.1f}".rstrip("0").rstrip(".")
    return f"{text}%"


# Compute relative uncertainties while treating scale-level numerical zero as undefined
def parameter_relative_uncertainties(
    values: np.ndarray, errors: np.ndarray, bounds: np.ndarray | None = None
) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    errors = np.asarray(errors, dtype=np.float64)
    if values.shape != errors.shape:
        raise ValueError("Parameter values and errors must have matching shapes")
    scale = np.abs(errors)
    if bounds is not None:
        limits = np.asarray(bounds, dtype=np.float64)
        if limits.shape != values.shape + (2,):
            raise ValueError("Parameter bounds have incompatible shape")
        scale = np.maximum(scale, np.abs(limits[..., 1] - limits[..., 0]))
    tolerance = np.sqrt(np.finfo(np.float64).eps) * np.maximum(scale, np.finfo(np.float64).eps)
    denominator = np.abs(values)
    valid = np.isfinite(denominator) & np.isfinite(errors) & (errors >= 0.0) & (denominator > tolerance)
    return np.divide(errors, denominator, out=np.full(values.shape, np.nan, dtype=np.float64), where=valid)


# Compute one user-facing parameter bound cell
def _format_bounds(bounds: np.ndarray | None, index: int, diagnostic: dict | None = None) -> str:
    if bounds is None:
        return "N/A"
    low: Any = _format_parameter_value(bounds[index, 0])
    high: Any = _format_parameter_value(bounds[index, 1])
    if diagnostic is not None and diagnostic["near_bound"]:
        color = "red" if diagnostic["at_bound"] else "yellow"
        attrs = ("bold",) if diagnostic["at_bound"] else ()
        if diagnostic["side"] == "lower":
            low = colored(low, color, attrs=list(attrs))
        else:
            high = colored(high, color, attrs=list(attrs))
    return f"[{low}, {high}]"


# Compute the common parameter-summary column labels for terminal or figure output
def parameter_summary_headers(*, include_parameter: bool, latex: bool = False, fit_label: str = "Surrogate") -> list[str]:
    surrogate = fit_label + (r" $\pm$ Hessian" if latex else " +/- Hessian")
    headers = ["Bounds", "Best trial", surrogate, "Relative uncertainty"]
    return ["Parameter", *headers] if include_parameter else headers


# Format one surrogate value and its absolute Hessian uncertainty
def format_surrogate_hessian(value: float, error: float, *, latex: bool = False) -> str:
    text = _format_parameter_value(value)
    if not np.isfinite(error):
        return text
    separator = r" $\pm$ " if latex else " +/- "
    return text + separator + _format_parameter_value(error)


# Build the common parameter-summary records used by terminal and figure outputs
def parameter_summary_records(
    *,
    param_names: Sequence[str],
    surrogate_values: np.ndarray,
    errors: np.ndarray,
    bounds: np.ndarray | None,
    realized_best: np.ndarray | None,
    gradient: np.ndarray | None = None,
) -> list[dict]:
    names = [str(name) for name in param_names]
    values = np.asarray(surrogate_values, dtype=np.float64).reshape(-1)
    errors = np.asarray(errors, dtype=np.float64).reshape(-1)
    limits = None if bounds is None else np.asarray(bounds, dtype=np.float64)
    realized = (
        np.full(len(names), np.nan, dtype=np.float64)
        if realized_best is None
        else np.asarray(realized_best, dtype=np.float64).reshape(-1)
    )
    if values.shape != (len(names),) or errors.shape != (len(names),):
        raise ValueError("Parameter summary vectors have inconsistent dimensions")
    if realized.shape != (len(names),):
        raise ValueError("Best-trial parameter vector has inconsistent dimensions")
    if limits is not None and limits.shape != (len(names), 2):
        raise ValueError("Parameter summary bounds have inconsistent dimensions")
    relative = parameter_relative_uncertainties(values, errors, limits)
    diagnostics = (
        parameter_bound_diagnostics(values, limits, gradient=gradient) if limits is not None else [None] * len(names)
    )
    return [
        {
            "name": name,
            "bounds": _format_bounds(limits, index, diagnostics[index]),
            "best_trial": _format_parameter_value(realized[index]),
            "surrogate": float(values[index]),
            "error": float(errors[index]),
            "relative_uncertainty": float(relative[index]),
            "relative_uncertainty_text": _plot_relative_uncertainty(relative[index]),
            "bound_diagnostic": diagnostics[index],
        }
        for index, name in enumerate(names)
    ]


# Evaluate a vectorized surrogate at one physical parameter point
def evaluate_scalar(func_hat: Callable, x: np.ndarray) -> float:
    try:
        z = np.asarray(func_hat(np.asarray(x, dtype=np.float64)[None, :]), dtype=np.float64).reshape(-1)
    except Exception:
        return float("inf")
    if len(z) == 0 or not np.isfinite(z[0]):
        return float("inf")
    return float(z[0])


# Evaluate a vectorized surrogate in bounded-size batches
def evaluate_batch(func_hat: Callable, X: np.ndarray, batch_size: int = 2048) -> np.ndarray:
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2:
        raise ValueError("Surrogate batch input must be a matrix")
    if batch_size < 1:
        raise ValueError("Surrogate batch size must be positive")
    values = []
    for start in range(0, len(X), int(batch_size)):
        stop = min(len(X), start + int(batch_size))
        prediction = np.asarray(func_hat(X[start:stop]), dtype=np.float64).reshape(-1)
        if len(prediction) != stop - start:
            raise ValueError("Surrogate returned an inconsistent batch length")
        values.append(prediction)
    return np.concatenate(values) if values else np.empty(0, dtype=np.float64)


# Build exact Torch value-gradient and Hessian callbacks for numerical optimizers
def build_torch_derivative_callbacks(
    objective: Callable, device: str | torch.device, dtype: torch.dtype
) -> tuple[Callable, Callable]:
    # Build one leaf tensor in the selected model arithmetic
    def differentiable_tensor(values: np.ndarray) -> torch.Tensor:
        return (
            torch.as_tensor(np.asarray(values, dtype=np.float64), dtype=dtype, device=device)
            .clone()
            .requires_grad_(True)
        )

    # Evaluate one scalar objective and its physical-coordinate gradient
    def value_gradient(values: np.ndarray) -> tuple[float, np.ndarray]:
        with torch.enable_grad():
            x = differentiable_tensor(values)
            value = objective(x)
            (gradient,) = torch.autograd.grad(value, x)
        return (float(value.detach().cpu()), gradient.detach().cpu().numpy().astype(np.float64, copy=False))

    # Evaluate the physical-coordinate Hessian with batched reverse AD
    def hessian(values: np.ndarray) -> np.ndarray:
        with torch.enable_grad():
            x = differentiable_tensor(values)
            try:
                matrix = torch.autograd.functional.hessian(objective, x, vectorize=True)
            except (NotImplementedError, RuntimeError):
                matrix = torch.autograd.functional.hessian(objective, x, vectorize=False)
        matrix = 0.5 * (matrix + matrix.T)
        return matrix.detach().cpu().numpy().astype(np.float64, copy=False)

    return value_gradient, hessian


# Evaluate one bounded central finite-difference gradient in a single batch
def vectorized_surrogate_gradient(
    func_hat: Callable, x: np.ndarray, bounds: np.ndarray, gradient_step: float = 1.0e-4
) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)
    widths = bounds[:, 1] - bounds[:, 0]
    unit = np.clip(to_unit(x[None, :], bounds)[0], 0.0, 1.0)
    dimension = len(unit)
    points = np.tile(unit, (2 * dimension, 1))
    denominator = np.empty(dimension, dtype=np.float64)
    for index in range(dimension):
        lower = max(0.0, unit[index] - gradient_step)
        upper = min(1.0, unit[index] + gradient_step)
        points[index, index] = upper
        points[dimension + index, index] = lower
        denominator[index] = (upper - lower) * widths[index]
    physical = from_unit(points, bounds)
    values = evaluate_batch(func_hat, physical)
    if not np.all(np.isfinite(values)):
        return np.zeros(dimension, dtype=np.float64)
    return (values[:dimension] - values[dimension:]) / denominator


# Select the sampled point with the lowest surrogate prediction
def select_surrogate_start(
    X: np.ndarray,
    Z: np.ndarray,
    func_hat: Callable,
    x0: np.ndarray | None = None,
    predictions: np.ndarray | None = None,
) -> tuple[np.ndarray, dict]:
    X = np.asarray(X, dtype=np.float64)
    Z = np.asarray(Z, dtype=np.float64)
    predictions = (
        evaluate_batch(func_hat, X) if predictions is None else np.asarray(predictions, dtype=np.float64).reshape(-1)
    )
    if len(predictions) != len(X):
        raise ValueError("Surrogate start predictions have inconsistent dimensions")
    finite = np.isfinite(predictions)
    if x0 is not None:
        selected = np.asarray(x0, dtype=np.float64).reshape(-1)
        source = "explicit"
        selected_prediction = evaluate_scalar(func_hat, selected)
        selected_observed = None
    elif np.any(finite):
        index = int(np.argmin(np.where(finite, predictions, np.inf)))
        selected = X[index].copy()
        source = "lowest surrogate prediction over sampled points"
        selected_prediction = float(predictions[index])
        selected_observed = float(Z[index])
    else:
        index = int(np.argmin(Z))
        selected = X[index].copy()
        source = "lowest observed point fallback"
        selected_prediction = evaluate_scalar(func_hat, selected)
        selected_observed = float(Z[index])
    return selected, {"source": source, "surrogate_z": selected_prediction, "observed_z": selected_observed}


# Build diverse deterministic starts for a bounded Torch global search
def torch_multistart_candidates(
    *,
    X: np.ndarray,
    Z: np.ndarray,
    func_hat: Callable,
    bounds: np.ndarray,
    x0: np.ndarray | None,
    restart_count: int,
    scan_points: int,
    rngseed: int,
) -> tuple[np.ndarray, dict, dict]:
    X = np.asarray(X, dtype=np.float64)
    Z = np.asarray(Z, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)
    if restart_count < 1 or scan_points < 0:
        raise ValueError("Torch global-search controls have inconsistent values")
    dimension = bounds.shape[0]
    sample_predictions = evaluate_batch(func_hat, X)
    primary, primary_diagnostics = select_surrogate_start(
        X=X, Z=Z, func_hat=func_hat, x0=x0, predictions=sample_predictions
    )

    midpoint = 0.5 * (bounds[:, 0] + bounds[:, 1])
    primary_unit = np.clip(to_unit(primary[None, :], bounds)[0], 0.0, 1.0)
    coordinate_units = []
    for index in range(dimension):
        for displacement in (-0.25, -0.10, 0.10, 0.25):
            point = primary_unit.copy()
            point[index] = np.clip(point[index] + displacement, 0.0, 1.0)
            coordinate_units.append(point)
    coordinate_points = (
        from_unit(np.asarray(coordinate_units, dtype=np.float64), bounds)
        if coordinate_units
        else np.empty((0, dimension), dtype=np.float64)
    )
    sobol_points = np.empty((0, dimension), dtype=np.float64)
    if scan_points > 0:
        engine = torch.quasirandom.SobolEngine(dimension=dimension, scramble=True, seed=int(rngseed))
        sobol_unit = engine.draw(int(scan_points)).to(dtype=torch.float64).cpu().numpy()
        sobol_points = from_unit(sobol_unit, bounds)

    extra_points = np.vstack((midpoint[None, :], coordinate_points, sobol_points))
    extra_predictions = evaluate_batch(func_hat, extra_points)
    pool = np.vstack((primary[None, :], X, extra_points))
    pool_values = np.concatenate(([float(primary_diagnostics["surrogate_z"])], sample_predictions, extra_predictions))
    pool_sources = (
        [f"primary: {primary_diagnostics['source']}"]
        + ["training sample"] * len(X)
        + ["box midpoint"]
        + ["coordinate perturbation"] * len(coordinate_points)
        + ["Sobol scan"] * len(sobol_points)
    )
    finite = np.isfinite(pool_values) & np.all(np.isfinite(pool), axis=1)
    if not np.any(finite):
        raise RuntimeError("Torch global search found no finite start candidates")

    order = np.argsort(np.where(finite, pool_values, np.inf))
    primary_index = 0
    selected_indices = [int(order[0])]
    if finite[primary_index] and primary_index != selected_indices[0]:
        selected_indices.append(primary_index)
    pool_unit = np.clip(to_unit(pool, bounds), 0.0, 1.0)

    # Favor distinct low-surrogate basins before filling any remaining slots
    for candidate_index in order:
        candidate_index = int(candidate_index)
        if len(selected_indices) >= restart_count or not finite[candidate_index]:
            break
        if candidate_index in selected_indices:
            continue
        separation = np.linalg.norm(pool_unit[selected_indices] - pool_unit[candidate_index], axis=1)
        if np.all(separation >= 0.05):
            selected_indices.append(candidate_index)
    for candidate_index in order:
        candidate_index = int(candidate_index)
        if len(selected_indices) >= restart_count or not finite[candidate_index]:
            break
        if candidate_index not in selected_indices:
            selected_indices.append(candidate_index)

    selected_indices = selected_indices[:restart_count]
    diagnostics = {
        "candidate_count": int(np.count_nonzero(finite)),
        "requested_restarts": int(restart_count),
        "selected_restarts": len(selected_indices),
        "sobol_points": int(scan_points),
        "coordinate_points": len(coordinate_points),
        "best_candidate_z": float(pool_values[selected_indices[0]]),
        "selected_candidate_z": [float(pool_values[index]) for index in selected_indices],
        "selected_sources": [pool_sources[index] for index in selected_indices],
    }
    return pool[selected_indices], primary_diagnostics, diagnostics


# Minimize a bounded surrogate with one vectorized finite-difference batch
def minimize_vectorized_lbfgsb(
    func_hat: Callable,
    start: np.ndarray,
    bounds: np.ndarray,
    maxiter: int = 200,
    gradient_step: float = 1.0e-4,
    value_gradient_func: Callable | None = None,
):
    start = np.asarray(start, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)
    if bounds.shape != (len(start), 2):
        raise ValueError("Optimizer bounds have inconsistent dimensions")
    widths = bounds[:, 1] - bounds[:, 0]
    if np.any(~np.isfinite(bounds)) or np.any(widths <= 0.0):
        raise ValueError("Vectorized L-BFGS-B requires finite non-degenerate bounds")
    if maxiter < 1 or not np.isfinite(gradient_step) or gradient_step <= 0.0:
        raise ValueError("Vectorized L-BFGS-B controls must be positive")

    # Map unit-cube optimizer rows into physical parameter coordinates
    def physical_points(unit_points):
        unit_points = np.asarray(unit_points, dtype=np.float64)
        return from_unit(unit_points, bounds)

    # Compute one objective and physical autograd gradient in unit coordinates
    def autograd_objective_with_gradient(unit_point):
        physical = physical_points(np.asarray(unit_point, dtype=np.float64))
        value, physical_gradient = value_gradient_func(physical)
        gradient = np.asarray(physical_gradient, dtype=np.float64) * widths
        if not np.isfinite(value) or not np.all(np.isfinite(gradient)):
            return float("inf"), np.zeros(len(unit_point), dtype=np.float64)
        return float(value), gradient

    # Compute one objective and central finite-difference fallback gradient
    def objective_with_gradient(unit_point):
        unit_point = np.asarray(unit_point, dtype=np.float64)
        dimension = len(unit_point)
        points = np.tile(unit_point, (1 + 2 * dimension, 1))
        denominator = np.empty(dimension, dtype=np.float64)
        for index in range(dimension):
            lower = max(0.0, unit_point[index] - gradient_step)
            upper = min(1.0, unit_point[index] + gradient_step)
            points[1 + index, index] = upper
            points[1 + dimension + index, index] = lower
            denominator[index] = upper - lower
        values = evaluate_batch(func_hat, physical_points(points))
        if not np.all(np.isfinite(values)):
            return float("inf"), np.zeros(dimension, dtype=np.float64)
        gradient = (values[1 : 1 + dimension] - values[1 + dimension :]) / denominator
        return float(values[0]), gradient

    unit_start = np.clip(to_unit(start[None, :], bounds)[0], 0.0, 1.0)
    result = minimize(
        (autograd_objective_with_gradient if value_gradient_func is not None else objective_with_gradient),
        unit_start,
        method="L-BFGS-B",
        jac=True,
        bounds=[(0.0, 1.0)] * len(start),
        options={"maxiter": int(maxiter), "ftol": 1.0e-10, "gtol": 1.0e-6, "maxls": 30},
    )
    result.x_unit = np.asarray(result.x, dtype=np.float64)
    result.x = physical_points(result.x_unit)
    return result


# Extract fitted covariance diagnostics with NaN placeholders when unavailable
def fitted_covariance_summary(
    m, param_names: Sequence[str], surrogate_vals: np.ndarray, bounds: np.ndarray | None = None
):
    D = len(param_names)
    cov_matrix = np.full((D, D), np.nan, dtype=np.float64)

    try:
        raw_cov = getattr(m, "covariance", None)
        if raw_cov is None:
            raise ValueError("missing covariance")

        candidate = np.asarray(raw_cov, dtype=np.float64)
        candidate = candidate.reshape(1, 1) if candidate.shape == () and D == 1 else np.atleast_2d(candidate)

        if candidate.shape != (D, D):
            raise ValueError(f"unexpected covariance shape {candidate.shape}")

        cov_matrix = candidate
    except Exception as exc:
        cprint(f"Fitted covariance unavailable ({exc}); storing NaN uncertainties", "yellow")

    surrogate_err, corr_matrix = covariance_summary(cov_matrix)

    surrogate_rel_err = parameter_relative_uncertainties(
        values=np.asarray(surrogate_vals, dtype=np.float64), errors=surrogate_err, bounds=bounds
    )

    return cov_matrix, corr_matrix, surrogate_err, surrogate_rel_err


# Extract fitted covariance and errors with a robust unavailable fallback
def fitted_covariance(m, param_names: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    dimension = len(param_names)
    if m is None:
        return (np.zeros((dimension, dimension), dtype=np.float64), np.full(dimension, np.nan, dtype=np.float64))
    values = np.asarray([m.values[name] for name in param_names], dtype=np.float64)
    covariance, _, errors, _ = fitted_covariance_summary(m, param_names, values)
    return covariance, errors


# Classify fitted parameters by distance and gradient at physical bounds
def parameter_bound_diagnostics(
    values: np.ndarray,
    bounds: np.ndarray,
    gradient: np.ndarray | None = None,
    *,
    bound_tolerance_fraction: float = 1.0e-6,
    near_bound_fraction: float = 0.01,
) -> list[dict]:
    fitted = np.asarray(values, dtype=np.float64)
    limits = np.asarray(bounds, dtype=np.float64)
    if limits.shape != (len(fitted), 2):
        raise ValueError("Bound diagnostics have inconsistent dimensions")
    widths = limits[:, 1] - limits[:, 0]
    if np.any(~np.isfinite(limits)) or np.any(widths <= 0.0):
        raise ValueError("Bound diagnostics require finite non-degenerate bounds")
    unit = (fitted - limits[:, 0]) / widths
    gradients = (
        np.full(len(fitted), np.nan, dtype=np.float64) if gradient is None else np.asarray(gradient, dtype=np.float64)
    )
    if gradients.shape != fitted.shape:
        raise ValueError("Bound diagnostics gradient has inconsistent dimensions")
    unit_gradients = gradients * widths

    output = []
    for index in range(len(fitted)):
        lower_distance = float(unit[index])
        upper_distance = float(1.0 - unit[index])
        side = "lower" if lower_distance <= upper_distance else "upper"
        distance_fraction = min(lower_distance, upper_distance)
        at_bound = bool(distance_fraction <= bound_tolerance_fraction)
        near_bound = bool(distance_fraction <= near_bound_fraction)
        projected_gradient = float(unit_gradients[index])
        outward_sign = 1.0 if side == "lower" else -1.0
        kkt_outward = bool(at_bound and outward_sign * projected_gradient >= 0.0)
        if kkt_outward:
            projected_gradient = 0.0
        output.append(
            {
                "status": side if at_bound else (f"near_{side}" if near_bound else "interior"),
                "side": side,
                "at_bound": at_bound,
                "near_bound": near_bound,
                "distance": float(distance_fraction * widths[index]),
                "distance_fraction": float(distance_fraction),
                "physical_gradient": (float(gradients[index]) if np.isfinite(gradients[index]) else None),
                "unit_gradient": (float(unit_gradients[index]) if np.isfinite(unit_gradients[index]) else None),
                "projected_unit_gradient": (projected_gradient if np.isfinite(projected_gradient) else None),
                "kkt_outward": kkt_outward,
            }
        )
    return output


# Attach parameter names to diagnostics for coordinates fixed at a bound
def named_active_bounds(param_names: Sequence[str], diagnostics: Sequence[dict]) -> list[dict]:
    return [
        {"parameter": str(name), **item}
        for name, item in zip(param_names, diagnostics, strict=True)
        if item["at_bound"]
    ]


# Quantify realized-trial support around one surrogate optimum
def surrogate_support_diagnostics(
    X: np.ndarray, Z: np.ndarray, values: np.ndarray, bounds: np.ndarray, surrogate_z: float
) -> dict:
    coordinates = np.asarray(X, dtype=np.float64)
    objective = np.asarray(Z, dtype=np.float64).reshape(-1)
    optimum = np.asarray(values, dtype=np.float64).reshape(-1)
    limits = np.asarray(bounds, dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape[0] != len(objective):
        raise ValueError("Surrogate support inputs have inconsistent rows")
    if coordinates.shape[1] != len(optimum) or limits.shape != (len(optimum), 2):
        raise ValueError("Surrogate support inputs have inconsistent dimensions")
    widths = limits[:, 1] - limits[:, 0]
    if np.any(~np.isfinite(widths)) or np.any(widths <= 0.0):
        raise ValueError("Surrogate support requires finite non-degenerate bounds")

    unit_coordinates = (coordinates - limits[:, 0]) / widths
    unit_optimum = (optimum - limits[:, 0]) / widths
    distances = np.linalg.norm(unit_coordinates - unit_optimum, axis=1)
    nearest = int(np.argmin(distances))
    realized_minimum = float(np.min(objective))
    optimism_gap = realized_minimum - float(surrogate_z)
    return {
        "nearest_trial_index": nearest,
        "nearest_trial_unit_distance": float(distances[nearest]),
        "nearest_trial_unit_rms_distance": float(distances[nearest] / np.sqrt(len(optimum))),
        "nearest_trial_z": float(objective[nearest]),
        "realized_minimum_z": realized_minimum,
        "surrogate_z": float(surrogate_z),
        "realized_minimum_minus_surrogate_z": float(optimism_gap),
        "surrogate_z_over_realized_minimum": (
            float(surrogate_z) / realized_minimum if realized_minimum != 0.0 else None
        ),
    }


# Compute one compact physical value for the parameter figure
def _plot_value(value: float) -> str:
    return _format_parameter_value(value)


# Compute one compact relative Hessian uncertainty for the parameter figure
def _plot_relative_uncertainty(relative: float) -> str:
    if not np.isfinite(relative):
        return "N/A"
    return _format_relative_percent(relative)


# Draw one interval with only its active boundary highlighted
def _add_plot_bound_interval(
    ax, *, row: int, bounds: np.ndarray, diagnostic: dict, x_position: float
) -> AnnotationBbox:
    low_color = "black"
    high_color = "black"
    low_weight = "normal"
    high_weight = "normal"
    if diagnostic["near_bound"]:
        color = "#e42536" if diagnostic["at_bound"] else "#f89c20"
        weight = "bold" if diagnostic["at_bound"] else "normal"
        if diagnostic["side"] == "lower":
            low_color, low_weight = color, weight
        else:
            high_color, high_weight = color, weight
    parts = (
        ("[", "black", "normal"),
        (_plot_value(bounds[row, 0]), low_color, low_weight),
        (", ", "black", "normal"),
        (_plot_value(bounds[row, 1]), high_color, high_weight),
        ("]", "black", "normal"),
    )
    packed = HPacker(
        children=[
            TextArea(text, textprops={"fontsize": 7, "color": color, "fontweight": weight})
            for text, color, weight in parts
        ],
        align="center",
        pad=0,
        sep=0,
    )
    annotation = AnnotationBbox(
        packed,
        (x_position, row),
        xycoords=ax.get_yaxis_transform(),
        box_alignment=(0.0, 0.5),
        frameon=False,
        pad=0.0,
        annotation_clip=False,
    )
    ax.add_artist(annotation)
    return annotation


# Shrink the parameter figure canvas while preserving the plotting width
def _compact_parameter_plot_layout(fig, ax) -> tuple[float, float]:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    label_width = max((label.get_window_extent(renderer).width for label in ax.get_yticklabels()), default=0.0)
    initial_width = float(fig.get_size_inches()[0])
    content_width = 0.65 * initial_width
    label_width_inches = label_width / max(float(fig.dpi), 1.0)
    left_padding = 0.30
    right_padding = 0.15
    figure_width = label_width_inches + left_padding + content_width + right_padding
    fig.set_size_inches(figure_width, fig.get_size_inches()[1], forward=True)
    left_margin = (label_width_inches + left_padding) / figure_width
    right_margin = 1.0 - right_padding / figure_width
    return float(left_margin), float(right_margin)


# Render one page of normalized parameter positions and physical values
def _parameter_uncertainty_figure(
    *,
    param_names: Sequence[str],
    best_fit: np.ndarray,
    errors: np.ndarray,
    bounds: np.ndarray,
    realized_best: np.ndarray | None,
    plot_brand: str,
    page_number: int,
    page_count: int,
    first_parameter: int,
    fit_label: str,
):
    best_fit = np.asarray(best_fit, dtype=np.float64)
    errors = np.asarray(errors, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)
    midpoint = 0.5 * (bounds[:, 0] + bounds[:, 1])
    half_width = 0.5 * (bounds[:, 1] - bounds[:, 0])
    normalized = (best_fit - midpoint) / half_width
    normalized_error = errors / half_width
    realized = None if realized_best is None else (np.asarray(realized_best, dtype=np.float64) - midpoint) / half_width
    summary_records = parameter_summary_records(
        param_names=param_names, surrogate_values=best_fit, errors=errors, bounds=bounds, realized_best=realized_best
    )

    y = np.arange(len(param_names), dtype=np.float64)
    height = max(4.0, 0.26 * len(param_names) + 1.5)
    fig = plt.figure(figsize=(15.5, height))
    grid = fig.add_gridspec(1, 2, width_ratios=(1.05, 1.35), wspace=0.025)
    ax = fig.add_subplot(grid[0, 0])
    ax_values = fig.add_subplot(grid[0, 1], sharey=ax)
    for row in range(len(param_names)):
        if (first_parameter + row) % 2 != 0:
            continue
        for axis in (ax, ax_values):
            axis.axhspan(row - 0.5, row + 0.5, color="#f3f3f3", zorder=0)
    ax.axvspan(-1.0, 1.0, color="#d9d9d9", alpha=0.5, zorder=0)
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.axvline(-1.0, color="#e42536", linewidth=0.8, linestyle=":")
    ax.axvline(1.0, color="#e42536", linewidth=0.8, linestyle=":")
    finite_best = np.isfinite(normalized)
    valid_error = finite_best & np.isfinite(normalized_error)
    ax.scatter(
        normalized[finite_best], y[finite_best], marker="o", color="black", s=18,
        label="Surrogate optimum" if fit_label == "Surrogate" else fit_label, zorder=3
    )
    ax.errorbar(
        normalized[valid_error],
        y[valid_error],
        xerr=normalized_error[valid_error],
        fmt="none",
        color="black",
        capsize=2.0,
        linewidth=1.0,
        label="Hessian uncertainty",
        zorder=3,
    )
    if realized is not None:
        finite_realized = np.isfinite(realized)
        ax.scatter(
            realized[finite_realized],
            y[finite_realized],
            marker="x",
            color="#e42536",
            s=25,
            label="Best realized trial",
            zorder=4,
        )
    active_bound = np.asarray([record["bound_diagnostic"]["at_bound"] for record in summary_records], dtype=bool)
    if np.any(active_bound):
        ax.scatter(
            normalized[active_bound],
            y[active_bound],
            marker="o",
            facecolors="none",
            edgecolors="#e42536",
            linewidths=1.3,
            s=55,
            label="At parameter bound",
            zorder=5,
        )
    ax.set_yticks(y, labels=[str(name) for name in param_names])
    ax.tick_params(axis="y", labelsize=7)
    ax.set_xlim(-1.3, 1.3)
    ax.set_xlabel(r"$(\theta-\theta_{\rm mid})/(\Delta\theta/2)$")
    ax.grid(axis="x", color="white", linewidth=0.8)
    ax.set_ylim(len(param_names) - 0.75, -0.25)

    ax_values.set_xlim(0.0, 1.0)
    ax_values.tick_params(axis="both", which="both", left=False, labelleft=False, bottom=False, labelbottom=False)
    for spine in ax_values.spines.values():
        spine.set_visible(False)
    column_positions = {"bounds": 0.02, "best_trial": 0.21, "surrogate": 0.42, "relative_uncertainty": 0.73}
    header_positions = (
        column_positions["bounds"],
        column_positions["best_trial"],
        column_positions["surrogate"],
        column_positions["relative_uncertainty"],
    )
    headers = zip(parameter_summary_headers(include_parameter=False, latex=True, fit_label=fit_label), header_positions, strict=True)
    for label, x_position in headers:
        ax_values.text(
            x_position, 1.002, label, transform=ax_values.transAxes, fontsize=8, fontweight="bold", va="bottom"
        )
    for row, record in enumerate(summary_records):
        diagnostic = record["bound_diagnostic"]
        realized_text = record["best_trial"]
        surrogate_text = format_surrogate_hessian(record["surrogate"], record["error"], latex=True)
        relative_text = record["relative_uncertainty_text"]
        ax_values.text(column_positions["best_trial"], row, realized_text, fontsize=7, va="center")
        ax_values.text(column_positions["surrogate"], row, surrogate_text, fontsize=7, va="center")
        ax_values.text(column_positions["relative_uncertainty"], row, relative_text, fontsize=7, va="center")
        _add_plot_bound_interval(
            ax_values, row=row, bounds=bounds, diagnostic=diagnostic, x_position=column_positions["bounds"]
        )

    left_margin, right_margin = _compact_parameter_plot_layout(fig, ax)
    figure_height = float(fig.get_size_inches()[1])
    fig.subplots_adjust(
        left=left_margin, right=right_margin, top=1.0 - 0.42 / figure_height, bottom=1.0 / figure_height
    )
    plot_box = ax.get_position()
    legend_x = plot_box.x0 + 0.5 * plot_box.width
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(legend_x, 0.012),
        ncol=2,
        fontsize=7,
        frameon=False,
        labelspacing=0.25,
        handletextpad=0.5,
        borderaxespad=0.0,
    )
    add_icetune_branding(fig, brand=plot_brand, x=0.025, y=0.995)
    if page_count > 1:
        final_parameter = first_parameter + len(param_names)
        fig.text(
            0.985,
            0.995,
            f"Parameters {first_parameter + 1}-{final_parameter}",
            ha="right",
            va="top",
            fontsize=7,
            color="#555555",
        )
        fig.text(0.985, 0.015, f"Page {page_number}/{page_count}", ha="right", va="bottom", fontsize=7, color="#555555")
    return fig


# Plot normalized parameter positions together with their physical values
def plot_parameter_uncertainties(
    *,
    param_names: Sequence[str],
    best_fit: np.ndarray,
    errors: np.ndarray,
    bounds: np.ndarray,
    output_root: pathlib.Path,
    realized_best: np.ndarray | None = None,
    plot_brand: str = "GRANIITTI",
    parameters_per_page: int = 30,
    fit_label: str = "Surrogate",
) -> dict:
    param_names = [str(name) for name in param_names]
    best_fit = np.asarray(best_fit, dtype=np.float64).reshape(-1)
    errors = np.asarray(errors, dtype=np.float64).reshape(-1)
    bounds = np.asarray(bounds, dtype=np.float64)
    realized = None if realized_best is None else np.asarray(realized_best, dtype=np.float64).reshape(-1)
    dimension = len(param_names)
    if dimension < 1:
        raise ValueError("Parameter uncertainty plot requires at least one parameter")
    if best_fit.shape != (dimension,) or errors.shape != (dimension,):
        raise ValueError("Parameter uncertainty vectors have inconsistent dimensions")
    if bounds.shape != (dimension, 2):
        raise ValueError("Parameter uncertainty bounds have inconsistent dimensions")
    if realized is not None and realized.shape != (dimension,):
        raise ValueError("Realized parameter vector has inconsistent dimensions")
    if int(parameters_per_page) < 1:
        raise ValueError("Parameters per page must be positive")

    output_root = pathlib.Path(output_root)
    ensure_dir(output_root)
    page_size = int(parameters_per_page)
    page_count = int(np.ceil(dimension / page_size))
    pdf_path = output_root / "parameter_uncertainties.pdf"
    png_paths = []
    pages = []
    with PdfPages(pdf_path) as pdf:
        for page_index, start in enumerate(range(0, dimension, page_size)):
            stop = min(start + page_size, dimension)
            page_number = page_index + 1
            fig = _parameter_uncertainty_figure(
                param_names=param_names[start:stop],
                best_fit=best_fit[start:stop],
                errors=errors[start:stop],
                bounds=bounds[start:stop],
                realized_best=None if realized is None else realized[start:stop],
                plot_brand=plot_brand,
                page_number=page_number,
                page_count=page_count,
                first_parameter=start,
                fit_label=fit_label,
            )
            png_name = (
                "parameter_uncertainties.png"
                if page_count == 1
                else f"parameter_uncertainties_page_{page_number:03d}.png"
            )
            png_path = output_root / png_name
            fig.savefig(png_path, dpi=220)
            pdf.savefig(fig)
            plt.close(fig)
            png_paths.append(str(png_path))
            pages.append(
                {"page": page_number, "first_parameter": start + 1, "last_parameter": stop, "png": str(png_path)}
            )
    return {
        "png": png_paths[0],
        "png_pages": png_paths,
        "pdf": str(pdf_path),
        "page_count": page_count,
        "parameters_per_page": page_size,
        "pages": pages,
    }


# Select a readable correlation-matrix layout for one parameter list
def correlation_matrix_layout(param_names: Sequence[str]) -> dict:
    dimension = len(param_names)
    longest_label = max((len(str(name)) for name in param_names), default=1)
    label_margin = min(4.0, 0.07 * longest_label)
    figure_extent = float(np.clip(6.0 + 0.18 * dimension + label_margin, 8.0, 24.0))
    tick_fontsize = float(np.clip(7.5 - 0.025 * max(0, dimension - 10), 2.5, 7.5))
    annotation_fontsize = float(np.clip(7.0 - 0.12 * dimension, 3.0, 6.0)) if dimension <= 30 else None
    return {
        "figure_extent": figure_extent,
        "tick_fontsize": tick_fontsize,
        "annotation_fontsize": annotation_fontsize,
        "title_fontsize": float(np.clip(10.0 + 0.03 * dimension, 10.0, 13.0)),
        "dpi": 300 if dimension <= 80 else (220 if dimension <= 140 else 180),
    }


# Mask the uninformative diagonal and select a symmetric off-diagonal color scale
def correlation_plot_values(correlation: np.ndarray) -> tuple[np.ndarray, float]:
    values = np.asarray(correlation, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        raise ValueError("Correlation matrix must be square")
    display_values = np.array(values, copy=True)
    np.fill_diagonal(display_values, np.nan)
    finite = np.abs(display_values[np.isfinite(display_values)])
    maximum = float(np.max(finite)) if finite.size else 0.0
    raw_limit = float(np.clip(max(0.1, maximum), 0.1, 1.0))
    color_limit = min(1.0, np.ceil(raw_limit * 10.0 - 1.0e-12) / 10.0)
    return display_values, float(color_limit)


# Plot one zero-centered parameter correlation matrix without its unit diagonal
def plot_correlation_matrix(
    *, correlation: np.ndarray, param_names: Sequence[str], output_root: pathlib.Path, plot_brand: str
) -> pathlib.Path:
    layout = correlation_matrix_layout(param_names)
    display_values, color_limit = correlation_plot_values(correlation)
    figure_extent = layout["figure_extent"]
    fig, ax = plt.subplots(figsize=(figure_extent, figure_extent))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e5e7eb")
    image = ax.imshow(
        np.ma.masked_invalid(display_values), interpolation="nearest", cmap=cmap, vmin=-color_limit, vmax=color_limit
    )
    ax.set_title("Parameter correlations", loc="left", fontsize=layout["title_fontsize"], fontweight="semibold", pad=10)
    colorbar = fig.colorbar(
        image,
        ax=ax,
        fraction=0.04,
        pad=0.025,
        shrink=0.86,
        aspect=28,
        ticks=np.linspace(-color_limit, color_limit, 5),
        format="%.2f",
    )
    colorbar.set_label(
        r"Correlation coefficient $\rho_{ij}$", fontsize=max(7.0, layout["tick_fontsize"] + 1.0), labelpad=8
    )
    colorbar.ax.tick_params(labelsize=layout["tick_fontsize"], length=3, width=0.6, colors="#374151")
    colorbar.outline.set_edgecolor("#9ca3af")
    colorbar.outline.set_linewidth(0.6)
    ticks = np.arange(len(param_names))
    ax.set_xticks(ticks, labels=param_names, rotation=90, ha="center", fontsize=layout["tick_fontsize"])
    ax.set_yticks(ticks, labels=param_names, fontsize=layout["tick_fontsize"])
    ax.tick_params(axis="both", which="both", length=0, pad=2, colors="#374151")
    for spine in ax.spines.values():
        spine.set_color("#9ca3af")
        spine.set_linewidth(0.6)
    if len(param_names) <= 30:
        boundaries = np.arange(-0.5, len(param_names), 1.0)
        ax.set_xticks(boundaries, minor=True)
        ax.set_yticks(boundaries, minor=True)
        ax.grid(which="minor", color="white", linewidth=0.45, alpha=0.35)
    annotation_fontsize = layout["annotation_fontsize"]
    if annotation_fontsize is not None:
        for row in range(len(param_names)):
            for column in range(len(param_names)):
                if row == column:
                    continue
                value = display_values[row, column]
                if not np.isfinite(value):
                    continue
                rgba = cmap(image.norm(value))
                luminance = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
                ax.text(
                    column,
                    row,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    color="#111827" if luminance > 0.55 else "white",
                    fontsize=annotation_fontsize,
                )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.985), pad=0.8)
    add_icetune_branding(fig, brand=plot_brand, x=0.01, y=0.995)
    output_root = pathlib.Path(output_root)
    ensure_dir(output_root)
    figure_path = output_root / "correlation_matrix.png"
    fig.savefig(figure_path, dpi=layout["dpi"], bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return figure_path


# Print and persist surrogate parameter diagnostics
def save_fit_parameters(
    m,
    param_names: Sequence[str],
    output_root: pathlib.Path,
    bounds: np.ndarray | None = None,
    realized_best: np.ndarray | None = None,
    realized_Z: float | None = None,
    optimization_diagnostics: dict | None = None,
    final_gradient: np.ndarray | None = None,
    plot_brand: str = "GRANIITTI",
    optimizer_label: str = "Torch L-BFGS-B + autograd Hessian",
    uncertainty_method: str = "torch_autograd_hessian",
) -> dict:
    output_root = pathlib.Path(output_root)
    ensure_dir(output_root)
    bounds_arr = None if bounds is None else np.asarray(bounds, dtype=np.float64)

    surrogate_vals = np.array([m.values[p] for p in param_names], dtype=np.float64)
    cov_matrix, corr_matrix, surrogate_err, _ = fitted_covariance_summary(
        m=m, param_names=param_names, surrogate_vals=surrogate_vals, bounds=bounds_arr
    )
    summary_records = parameter_summary_records(
        param_names=param_names,
        surrogate_values=surrogate_vals,
        errors=surrogate_err,
        bounds=bounds_arr,
        realized_best=realized_best,
        gradient=final_gradient,
    )
    bound_rows = [record["bound_diagnostic"] for record in summary_records]

    title_parts = [f"Surrogate parameter summary after {optimizer_label} (Z = {float(m.fmin.fval):.3f})"]
    if realized_Z is not None:
        title_parts.append(f"realized Z = {float(realized_Z):.3f}")
    rows = []
    for record in summary_records:
        surrogate_value = format_surrogate_hessian(record["surrogate"], record["error"])
        boundary = record["bound_diagnostic"]
        if boundary is not None and boundary["at_bound"]:
            surrogate_value = StyledCell(surrogate_value, color="red", attrs=("bold",))
        elif boundary is not None and boundary["near_bound"]:
            surrogate_value = StyledCell(surrogate_value, color="yellow", attrs=("bold",))
        rows.append(
            [
                record["name"],
                record["bounds"],
                record["best_trial"],
                surrogate_value,
                record["relative_uncertainty_text"],
            ]
        )
    print_table(
        " | ".join(title_parts),
        parameter_summary_headers(include_parameter=True),
        rows,
        align=["left", "right", "right", "right", "right"],
    )

    valid_correlation = np.isfinite(np.diag(corr_matrix))
    if np.any(valid_correlation):
        fig_path = plot_correlation_matrix(
            correlation=corr_matrix, param_names=param_names, output_root=output_root / "optimizer", plot_brand=plot_brand
        )
        offdiag = corr_matrix[~np.eye(len(param_names), dtype=bool)]
        finite_offdiag = offdiag[np.isfinite(offdiag)]
        correlation_rows = [
            ["dimension", len(param_names)],
            ["available parameters", int(np.count_nonzero(valid_correlation))],
            ["unavailable parameters", int(np.count_nonzero(~valid_correlation))],
            ["finite entries", int(np.count_nonzero(np.isfinite(corr_matrix)))],
            [
                "max |rho_ij| off diagonal",
                (_format_float(float(np.max(np.abs(finite_offdiag))), digits=3) if finite_offdiag.size else "N/A"),
            ],
            [
                "mean |rho_ij| off diagonal",
                (_format_float(float(np.mean(np.abs(finite_offdiag))), digits=3) if finite_offdiag.size else "N/A"),
            ],
            ["figure", fig_path],
        ]
        print_table("Correlation matrix summary:", ["Quantity", "Value"], correlation_rows)
    else:
        print_table(
            "Correlation matrix summary:",
            ["Quantity", "Value"],
            [
                ["dimension", len(param_names)],
                ["finite entries", int(np.count_nonzero(np.isfinite(corr_matrix)))],
                ["figure", "skipped because fitted covariance is unavailable"],
            ],
        )

    parameter_plot = None
    if bounds_arr is not None:
        parameter_plot = plot_parameter_uncertainties(
            param_names=param_names,
            best_fit=surrogate_vals,
            errors=surrogate_err,
            bounds=bounds_arr,
            output_root=output_root / "optimizer",
            realized_best=realized_best,
            plot_brand=plot_brand,
        )

    surrogate_best = {
        param_names[i]: {
            "value": float(surrogate_vals[i]),
            "uncertainty": float(surrogate_err[i]),
            "relative_error_percent": float(summary_records[i]["relative_uncertainty"] * 100),
        }
        for i in range(len(param_names))
    }
    data_to_save = {
        "basis": "optimizer",
        "best_fit": fit_summary.build_best_fit(
            values={param_names[i]: surrogate_vals[i] for i in range(len(param_names))},
            uncertainties={param_names[i]: surrogate_err[i] for i in range(len(param_names))},
            objective_name="Z",
            objective_value=m.fmin.fval,
            source="surrogate",
            uncertainty_method=uncertainty_method,
        ),
        "surrogate_best": surrogate_best,
        "surrogate_Z": float(m.fmin.fval),
        "covariance_matrix": cov_matrix.tolist(),
        "correlation_matrix": corr_matrix.tolist(),
        "bound_diagnostics": {param_names[i]: bound_rows[i] for i in range(len(param_names))}
        if bounds_arr is not None
        else None,
        "realized_best": None,
        "realized_Z": None,
        "optimization_diagnostics": optimization_diagnostics,
        "parameter_uncertainty_plot": parameter_plot,
        "plot_brand": str(plot_brand),
    }

    if realized_best is not None and realized_Z is not None:
        data_to_save["realized_best"] = {param_names[i]: float(realized_best[i]) for i in range(len(param_names))}
        data_to_save["realized_Z"] = float(realized_Z)

    json_filename = output_root / "parameters.json"
    write_json_file(json_filename, json_safe(data_to_save), indent=4)
    outputs = [["parameter JSON", json_filename]]
    if parameter_plot is not None:
        for page in parameter_plot["pages"]:
            label = (
                "parameter uncertainty PNG"
                if parameter_plot["page_count"] == 1
                else f"parameter uncertainty PNG page {page['page']}"
            )
            outputs.append([label, page["png"]])
        outputs.append(["parameter uncertainty PDF", parameter_plot["pdf"]])
    print_table("Saved outputs:", ["Output", "Path"], outputs)

    return data_to_save


# Save decoded parameter uncertainties without assigning independent bounds to coupled outputs
def save_physical_parameters(*, transform, best_fit, covariance, output_root, realized_best=None,
                             plot_brand="GRANIITTI", fit_label="Surrogate optimum"):
    output_root = pathlib.Path(output_root) / "physical"
    ensure_dir(output_root)
    result = transform.propagate(best_fit, covariance)
    result["basis"] = "physical"
    realized = None if realized_best is None else transform.samples([realized_best])[0]
    if realized is not None:
        realized = result["values"] + transform.difference(torch.as_tensor(realized), torch.as_tensor(result["values"])).numpy()
    result["realized_best"] = realized
    write_json_file(output_root / "parameters.json", json_safe(result), indent=4)
    if np.any(np.isfinite(result["correlation"])):
        plot_correlation_matrix(correlation=result["correlation"], param_names=result["names"],
                                output_root=output_root, plot_brand=plot_brand)
    with PdfPages(output_root / "parameter_uncertainties.pdf") as pdf:
        for start in range(0, len(result["names"]), 8):
            count = min(8, len(result["names"]) - start)
            columns = min(2, count)
            rows = (count + columns - 1) // columns
            fig, axes = plt.subplots(rows, columns, figsize=(11.7, 2.1 * rows + 0.8), squeeze=False)
            for index, axis in enumerate(axes.flat, start):
                if index >= len(result["names"]):
                    axis.set_visible(False)
                    continue
                value, error = result["values"][index], result["errors"][index]
                axis.set_title(result["names"][index], fontsize=8)
                axis.errorbar(value, 0.0, xerr=error if np.isfinite(error) else None,
                              fmt="o", color="black", capsize=3, label=fit_label)
                if realized is not None:
                    axis.scatter(realized[index], 0.0, marker="x", color="#e42536", label="Best realized trial")
                axis.set_yticks([])
                axis.set_xlabel(format_surrogate_hessian(value, error, latex=True), fontsize=8)
                axis.grid(axis="x", alpha=0.3)
            axes.flat[0].legend(fontsize=7)
            fig.suptitle("Physical parameters (local covariance approximation)", fontsize=10)
            add_icetune_branding(fig, brand=plot_brand)
            fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
            filename = "parameter_uncertainties.png" if len(result["names"]) <= 8 else f"parameter_uncertainties_page_{start // 8 + 1:03d}.png"
            fig.savefig(output_root / filename, dpi=220)
            pdf.savefig(fig)
            plt.close(fig)
    return result


# Optimize a surrogate with GPU capable bounded Torch L-BFGS
def optimize_with_torch(
    X: np.ndarray,
    Z: np.ndarray,
    param_names: Sequence[str],
    func_hat: Callable,
    objective_torch: Callable,
    objective_batch_torch: Callable,
    bounds: np.ndarray,
    output_root: pathlib.Path,
    device: str | torch.device,
    dtype: torch.dtype,
    realized_best: np.ndarray | None = None,
    realized_Z: float | None = None,
    x0: np.ndarray | None = None,
    maxiter: int = 200,
    restarts: int = 32,
    scan_points: int = 1024,
    rngseed: int = 1234,
    gradient_tolerance: float = 1.0e-6,
    relative_gradient_tolerance: float = 1.0e-6,
    plot_brand: str = "GRANIITTI",
) -> tuple[lbfgsb.Result, np.ndarray]:
    X = np.asarray(X, dtype=np.float64)
    Z = np.asarray(Z, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)
    if bounds.shape != (len(param_names), 2):
        raise ValueError("Surrogate optimizer bounds have inconsistent dimensions")
    if np.any(~np.isfinite(bounds)) or np.any(bounds[:, 1] <= bounds[:, 0]):
        raise ValueError("Surrogate optimizer requires finite non-degenerate bounds")
    if dtype != torch.float64:
        raise ValueError("Torch L-BFGS surrogate optimization requires float64 arithmetic")

    starts, start_diagnostics, global_search = torch_multistart_candidates(
        X=X,
        Z=Z,
        func_hat=func_hat,
        bounds=bounds,
        x0=x0,
        restart_count=int(restarts),
        scan_points=int(scan_points),
        rngseed=int(rngseed),
    )
    lower = torch.as_tensor(bounds[:, 0], dtype=dtype, device=device)
    widths = torch.as_tensor(bounds[:, 1] - bounds[:, 0], dtype=dtype, device=device)
    start_tensor = torch.as_tensor(starts, dtype=dtype, device=device)
    unit_starts = torch.clamp((start_tensor - lower[None, :]) / widths[None, :], 0.0, 1.0)

    # Map normalized optimizer coordinates into physical model parameters
    def unit_objective(unit_values):
        return objective_torch(lower + unit_values * widths)

    # Map a batch of normalized rows into physical model parameters
    def unit_objective_batch(unit_values):
        return objective_batch_torch(lower[None, :] + unit_values * widths[None, :])

    realized_surrogate_z = evaluate_scalar(func_hat, realized_best) if realized_best is not None else None
    residual = (
        float(realized_Z - realized_surrogate_z)
        if realized_Z is not None and realized_surrogate_z is not None
        else None
    )
    cprint("Torch batched multistart L-BFGS-B step:", "yellow")
    with tqdm(total=maxiter, desc="Torch L-BFGS-B multistart", unit="iter") as progress:
        # Update the batched descent display after every global L-BFGS iteration
        def update_batch_progress(iteration: int, active: int, best: float) -> None:
            progress.update(iteration - progress.n)
            progress.set_postfix(best=f"{best:.6g}", active=active)

        batch_units, batch_objectives, _, batch_descent = lbfgsb.minimize_batch(
            objective=unit_objective_batch,
            start=unit_starts,
            maxiter=maxiter,
            gradient_tolerance=gradient_tolerance,
            relative_gradient_tolerance=relative_gradient_tolerance,
            progress_callback=update_batch_progress,
        )
    finite_batch = torch.isfinite(batch_objectives)
    if not bool(torch.any(finite_batch)):
        raise RuntimeError("Torch multistart optimization returned no finite minimum")
    best_restart = int(torch.argmin(torch.where(finite_batch, batch_objectives, torch.inf)).item())

    # Polish the best batched basin with one strict scalar L-BFGS solve
    with tqdm(total=maxiter, desc="Torch L-BFGS-B polish", unit="iter") as progress:
        # Update the selected-basin descent display after every L-BFGS iteration
        def update_polish_progress(iteration: int, objective: float) -> None:
            progress.update(iteration - progress.n)
            progress.set_postfix(objective=f"{objective:.6g}")

        unit_values, objective_value, unit_gradient, polish = lbfgsb.minimize(
            objective=unit_objective,
            start=batch_units[best_restart],
            maxiter=maxiter,
            gradient_tolerance=gradient_tolerance,
            relative_gradient_tolerance=relative_gradient_tolerance,
            progress_callback=update_polish_progress,
        )
    batch_best_value = batch_objectives[best_restart]
    if not torch.isfinite(objective_value) or bool(objective_value > batch_best_value):
        unit_values = batch_units[best_restart]
        objective_value = batch_best_value
        _, unit_gradient = lbfgsb.value_gradient(unit_objective, unit_values)
    projected_unit_gradient = lbfgsb.projected_gradient(unit_values, unit_gradient)
    batch_success = batch_descent["success"].detach().cpu().numpy().astype(bool)
    batch_values = batch_objectives.detach().cpu().numpy().astype(np.float64, copy=False)
    agreement_tolerance = max(1.0e-7, 1.0e-7 * abs(float(batch_best_value.detach().cpu())))
    agreement_count = int(
        np.count_nonzero(np.isfinite(batch_values) & (batch_values <= batch_values[best_restart] + agreement_tolerance))
    )
    selected_batch_success = bool(batch_success[best_restart])
    polish_success = bool(polish["success"])
    descent = {
        "success": bool(polish_success or selected_batch_success),
        "message": (
            str(polish["message"])
            if polish_success or not selected_batch_success
            else str(batch_descent["messages"][best_restart])
        ),
        "iterations": int(polish["iterations"]),
        "function_evaluations": int(polish["function_evaluations"]),
        "gradient_evaluations": int(polish["gradient_evaluations"]),
        "max_abs_unit_gradient": float(torch.max(torch.abs(unit_gradient)).item()),
        "max_abs_projected_unit_gradient": float(torch.max(torch.abs(projected_unit_gradient)).item()),
        "initial_max_abs_projected_unit_gradient": float(polish["initial_max_abs_projected_unit_gradient"]),
        "effective_gradient_tolerance": float(polish["effective_gradient_tolerance"]),
        "relative_gradient_tolerance": float(relative_gradient_tolerance),
        "device": str(unit_values.device),
        "dtype": str(unit_values.dtype),
        "history_size": int(polish["history_size"]),
        "history_resets": int(polish["history_resets"]),
    }
    scan_tolerance = max(1.0e-7, 1.0e-7 * abs(float(global_search["best_candidate_z"])))
    if float(objective_value.detach().cpu()) > global_search["best_candidate_z"] + scan_tolerance:
        raise RuntimeError("Torch optimizer failed its global candidate objective audit")

    physical_values = lower + unit_values * widths
    physical_objective_value, physical_gradient = lbfgsb.value_gradient(objective_torch, physical_values)
    surrogate_best = physical_values.detach().cpu().numpy().astype(np.float64, copy=False)
    final_gradient = physical_gradient.detach().cpu().numpy().astype(np.float64, copy=False)
    final_bounds = parameter_bound_diagnostics(values=surrogate_best, bounds=bounds, gradient=final_gradient)
    free_mask = torch.as_tensor(
        [not item["at_bound"] for item in final_bounds], dtype=torch.bool, device=physical_values.device
    )
    hessian, covariance, curvature = lbfgsb.hessian_covariance(
        objective=objective_torch, values=physical_values, errordef=1.0, free_mask=free_mask
    )
    if not curvature["free_positive_definite"]:
        descent["success"] = False
        descent["message"] += ": free parameter Hessian is not positive definite"
    result = lbfgsb.build_result(
        param_names=param_names,
        values=physical_values,
        objective_value=physical_objective_value,
        covariance=covariance,
        success=descent["success"],
        message=descent["message"],
    )
    active_bounds = named_active_bounds(param_names, final_bounds)
    sample_support = surrogate_support_diagnostics(
        X=X, Z=Z, values=surrogate_best, bounds=bounds, surrogate_z=result.fmin.fval
    )
    optimization_diagnostics = {
        "backend": "torch",
        "start": start_diagnostics,
        "global_search": {
            **global_search,
            "selected_restart": best_restart,
            "successful_restarts": int(np.count_nonzero(batch_success)),
            "best_basin_agreement": agreement_count,
            "agreement_tolerance": agreement_tolerance,
            "batch_evaluations": int(batch_descent["batch_evaluations"]),
            "row_evaluations": int(batch_descent["row_evaluations"]),
        },
        "derivatives": "torch autograd",
        "hessian": "torch autograd",
        "device": str(device),
        "dtype": str(dtype),
        "realized_observed_z": realized_Z,
        "realized_surrogate_z": realized_surrogate_z,
        "realized_residual": residual,
        "lbfgsb": {**descent, "function_value": float(objective_value.detach().cpu()), "active_bounds": active_bounds},
        "curvature": {**curvature, "hessian_frobenius_norm": float(torch.linalg.matrix_norm(hessian).cpu())},
        "final": {
            "surrogate_z": result.fmin.fval,
            "observed_minimum_minus_surrogate_optimum": (
                None if realized_Z is None else float(realized_Z - result.fmin.fval)
            ),
            "max_abs_physical_gradient": float(np.max(np.abs(final_gradient))),
            "active_bound_count": len(active_bounds),
            "active_bounds": active_bounds,
            "sample_support": sample_support,
        },
    }
    print_table(
        None,
        ["Field", "Value"],
        [
            ["success", descent["success"]],
            ["message", descent["message"]],
            ["function value", result.fmin.fval],
            ["restart candidates", global_search["selected_restarts"]],
            ["successful restarts", int(np.count_nonzero(batch_success))],
            ["restarts agreeing at best basin", agreement_count],
            ["global scan candidates", global_search["candidate_count"]],
            ["best scanned objective", global_search["best_candidate_z"]],
            ["batched function evaluations", batch_descent["batch_evaluations"]],
            ["iterations", descent["iterations"]],
            ["function evaluations", descent["function_evaluations"]],
            ["max |unit gradient|", descent["max_abs_unit_gradient"]],
            ["max |projected unit gradient|", descent["max_abs_projected_unit_gradient"]],
            ["effective gradient tolerance", descent["effective_gradient_tolerance"]],
            ["L-BFGS history resets", descent["history_resets"]],
            ["batched line-search retries", batch_descent["line_search_retries"]],
            ["device", descent["device"]],
            ["dtype", descent["dtype"]],
            ["full Hessian positive definite", curvature["positive_definite"]],
            ["free Hessian positive definite", curvature["free_positive_definite"]],
            ["minimum free Hessian eigenvalue", curvature["free_minimum_eigenvalue"]],
            ["active bound parameters", curvature["active_dimension"]],
            ["regularized free eigenvalues", curvature["regularized_eigenvalues"]],
        ],
    )
    print()
    if agreement_count < 2:
        cprint(
            "WARNING: only one restart reached the best surrogate basin; the profile consistency audit will verify this minimum",
            "red",
            attrs=["bold"],
        )
        print()
    if active_bounds:
        cprint(
            "WARNING: surrogate optimum is constrained by parameter bounds: "
            + ", ".join(f"{item['parameter']} ({item['side']})" for item in active_bounds),
            "red",
            attrs=["bold"],
        )
        print()
    save_fit_parameters(
        m=result,
        param_names=param_names,
        output_root=output_root,
        bounds=bounds,
        realized_best=realized_best,
        realized_Z=realized_Z,
        optimization_diagnostics=optimization_diagnostics,
        final_gradient=final_gradient,
        plot_brand=plot_brand,
        optimizer_label="Torch L-BFGS-B + autograd Hessian",
        uncertainty_method="torch_autograd_hessian",
    )
    return result, surrogate_best

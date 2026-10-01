# Per-observable parameter impact surrogates and plots
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import pathlib
import re
from collections.abc import Callable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import torch
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

from core.inference import fit, profile
from core.io.files import ensure_dir
from core.io.serialize import finite_or_none, json_safe, write_json_file
from core.plot import plot
from core.plot.style import add_icetune_branding


# Build one random Fourier feature map for a smooth RBF surrogate
def build_rff_map(input_dim: int, n_features: int, scale: float, rng: np.random.Generator) -> dict:
    if input_dim < 1 or n_features < 1:
        raise ValueError("RFF dimensions must be positive")
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("RFF scale must be positive")
    omega = rng.normal(loc=0.0, scale=1.0 / float(scale), size=(int(input_dim), int(n_features)))
    phase = rng.uniform(0.0, 2.0 * np.pi, size=int(n_features))
    return {
        "omega": omega.astype(np.float64),
        "phase": phase.astype(np.float64),
        "scale": float(scale),
        "n_features": int(n_features),
    }


# Estimate a stable RBF scale from random pairwise feature distances
def median_pairwise_scale(X_model: np.ndarray, rng: np.random.Generator, max_pairs: int = 4096) -> float:
    X_model = np.asarray(X_model, dtype=np.float64)
    if X_model.ndim != 2 or X_model.shape[0] < 2:
        raise ValueError("At least two feature rows are required for automatic RBF scaling")
    n_pairs = min(int(max_pairs), max(32, 4 * X_model.shape[0]))
    left = rng.integers(0, X_model.shape[0], size=n_pairs)
    right = rng.integers(0, X_model.shape[0], size=n_pairs)
    distance = np.linalg.norm(X_model[left] - X_model[right], axis=1)
    distance = distance[np.isfinite(distance) & (distance > 0.0)]
    if distance.size == 0:
        return 1.0
    return float(np.median(distance))


# Evaluate bias, linear and random Fourier features
def rff_feature_matrix(X_model: np.ndarray, rff: dict) -> np.ndarray:
    X_model = np.asarray(X_model, dtype=np.float64)
    if X_model.ndim != 2:
        raise ValueError("X_model must be a two-dimensional array")
    omega = np.asarray(rff["omega"], dtype=np.float64)
    phase = np.asarray(rff["phase"], dtype=np.float64)
    if omega.shape[0] != X_model.shape[1]:
        raise ValueError("RFF input dimension does not match X_model")
    n_features = int(rff["n_features"])
    nonlinear = np.sqrt(2.0 / n_features) * np.cos(X_model @ omega + phase)
    return np.hstack((np.ones((X_model.shape[0], 1), dtype=np.float64), X_model, nonlinear))


# Fit ridge coefficients in the existing bias-first matrix convention
def fit_ridge_coefficients(Phi: np.ndarray, target: np.ndarray, alpha: float) -> np.ndarray:
    fitted = Ridge(alpha=float(alpha)).fit(Phi[:, 1:], target)
    return np.vstack((np.atleast_1d(fitted.intercept_), np.atleast_2d(fitted.coef_).T))


# Solve ridge regression by SVD while leaving the intercept unpenalized
def fit_ridge_tensor(Phi: torch.Tensor, target: torch.Tensor, alpha: float) -> torch.Tensor:
    if not np.isfinite(alpha) or alpha < 0.0:
        raise ValueError("Ridge regularization must be finite and nonnegative")
    features = Phi[:, 1:].to(dtype=torch.float64)
    target = target.to(device=features.device, dtype=features.dtype)
    if target.ndim == 1:
        target = target[:, None]
    feature_mean, target_mean = features.mean(0), target.mean(0)
    centered = features - feature_mean
    left, singular, right = torch.linalg.svd(centered, full_matrices=False)
    if alpha > 0.0:
        gain = singular / (singular.square() + float(alpha))
    else:
        threshold = torch.finfo(singular.dtype).eps * max(centered.shape) * singular.max()
        retained = singular > threshold
        gain = torch.where(retained, 1.0 / torch.where(retained, singular, 1.0), 0.0)
    slopes = right.T @ (gain[:, None] * (left.T @ (target - target_mean)))
    intercept = target_mean - feature_mean @ slopes
    return torch.cat((intercept[None, :], slopes), dim=0)


# Group observable columns that share the same finite trial mask
def finite_mask_groups(Y: np.ndarray, min_trials: int) -> list[tuple[np.ndarray, list[int]]]:
    groups = {}
    for column in range(Y.shape[1]):
        mask = np.isfinite(Y[:, column])
        if int(np.count_nonzero(mask)) < int(min_trials):
            continue
        groups.setdefault(mask.tobytes(), {"mask": mask, "columns": []})["columns"].append(column)
    return [(entry["mask"], entry["columns"]) for entry in groups.values()]


# Fit a smooth multi-output surrogate to observable-level chi2 values
def fit_observable_surrogate(
    X_model: np.ndarray,
    observable_chi2: np.ndarray,
    *,
    n_features: int = 256,
    scale: float | None = None,
    alpha: float = 1.0e-3,
    validation_fraction: float = 0.2,
    rngseed: int = 1234,
    min_trials: int = 30,
) -> dict:
    X_model = np.asarray(X_model, dtype=np.float64)
    Y = np.asarray(observable_chi2, dtype=np.float64)
    if X_model.ndim != 2 or Y.ndim != 2 or X_model.shape[0] != Y.shape[0]:
        raise ValueError("Impact surrogate inputs must be aligned two-dimensional arrays")
    if not (0.0 < validation_fraction < 1.0):
        raise ValueError("Impact validation fraction must be in (0, 1)")
    if min_trials < 3:
        raise ValueError("Impact minimum trial count must be at least 3")
    # Correlated precision partitions can have signed observable contributions

    rng = np.random.default_rng(int(rngseed) + 15485863)
    selected_scale = median_pairwise_scale(X_model, rng) if scale is None else float(scale)
    rff = build_rff_map(input_dim=X_model.shape[1], n_features=int(n_features), scale=selected_scale, rng=rng)
    Phi = rff_feature_matrix(X_model, rff)
    coefficients = np.full((Phi.shape[1], Y.shape[1]), np.nan, dtype=np.float64)
    response_mean = np.full(Y.shape[1], np.nan, dtype=np.float64)
    response_std = np.full(Y.shape[1], np.nan, dtype=np.float64)
    validation_r2 = np.full(Y.shape[1], np.nan, dtype=np.float64)
    validation_mae = np.full(Y.shape[1], np.nan, dtype=np.float64)
    trial_count = np.zeros(Y.shape[1], dtype=np.int64)

    for group_index, (finite_mask, columns) in enumerate(finite_mask_groups(Y, min_trials)):
        row_indices = np.flatnonzero(finite_mask)
        group_rng = np.random.default_rng(int(rngseed) + 32452843 + group_index)
        shuffled = group_rng.permutation(row_indices)
        n_validation = max(1, int(round(validation_fraction * len(shuffled))))
        n_validation = min(n_validation, len(shuffled) - 2)
        validation_rows = shuffled[:n_validation]
        training_rows = shuffled[n_validation:]
        group_columns = np.asarray(columns, dtype=np.int64)

        training_target = Y[np.ix_(training_rows, group_columns)]
        train_scaled, train_mean, train_std = fit.standardize(training_target)
        validation_coef = fit_ridge_coefficients(Phi[training_rows], train_scaled, alpha)
        validation_prediction = (Phi[validation_rows] @ validation_coef) * train_std + train_mean
        validation_target = Y[np.ix_(validation_rows, group_columns)]
        validation_r2[group_columns] = r2_score(validation_target, validation_prediction, multioutput="raw_values")
        validation_mae[group_columns] = np.mean(np.abs(validation_target - validation_prediction), axis=0)

        full_target = Y[np.ix_(row_indices, group_columns)]
        full_scaled, full_mean, full_std = fit.standardize(full_target)
        coefficients[:, group_columns] = fit_ridge_coefficients(Phi[row_indices], full_scaled, alpha)
        response_mean[group_columns] = full_mean
        response_std[group_columns] = full_std
        trial_count[group_columns] = len(row_indices)

    return {
        "rff": rff,
        "coefficients": coefficients,
        "response_mean": response_mean,
        "response_std": response_std,
        "validation_r2": validation_r2,
        "validation_mae": validation_mae,
        "trial_count": trial_count,
        "trained": np.all(np.isfinite(coefficients), axis=0),
        "target_transform": "standard",
        "reference_X_model": X_model.copy(),
        "reference_chi2": Y.copy(),
        "alpha": float(alpha),
        "min_trials": int(min_trials),
    }


# Predict all trained observable chi2 values from topology-aware features
def predict_observable_chi2(model: dict, X_model: np.ndarray) -> np.ndarray:
    Phi = rff_feature_matrix(X_model, model["rff"])
    prediction = Phi @ np.asarray(model["coefficients"], dtype=np.float64)
    prediction = prediction * np.asarray(model["response_std"], dtype=np.float64) + np.asarray(
        model["response_mean"], dtype=np.float64
    )
    prediction[:, ~np.asarray(model["trained"], dtype=bool)] = np.nan
    return prediction


# Select a realized chi2 contribution nearest to one model point
def nearest_realized_chi2(model: dict, X_model: np.ndarray) -> dict:
    query = np.asarray(X_model, dtype=np.float64).reshape(1, -1)
    reference_X = np.asarray(model["reference_X_model"], dtype=np.float64)
    reference_Y = np.asarray(model["reference_chi2"], dtype=np.float64)
    if reference_X.ndim != 2 or reference_X.shape[1] != query.shape[1]:
        raise ValueError("Impact reference coordinates have inconsistent dimensions")
    if reference_Y.ndim != 2 or reference_Y.shape[0] != reference_X.shape[0]:
        raise ValueError("Impact reference chi2 values have inconsistent rows")
    distance = np.linalg.norm(reference_X - query, axis=1)
    values = np.full(reference_Y.shape[1], np.nan, dtype=np.float64)
    indices = np.full(reference_Y.shape[1], -1, dtype=np.int64)
    distances = np.full(reference_Y.shape[1], np.nan, dtype=np.float64)
    for column in range(reference_Y.shape[1]):
        valid = np.isfinite(reference_Y[:, column])
        if not np.any(valid):
            continue
        index = int(np.argmin(np.where(valid, distance, np.inf)))
        values[column] = reference_Y[index, column]
        indices[column] = index
        distances[column] = distance[index]
    return {"kind": "nearest_realized_trial", "chi2": values, "trial_index": indices, "feature_distance": distances}


# Build fixed one-at-a-time Hessian shift points
def fixed_parameter_shifts(best_fit: np.ndarray, errors: np.ndarray, bounds: np.ndarray) -> dict:
    best_fit = np.asarray(best_fit, dtype=np.float64)
    errors = np.asarray(errors, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)
    if bounds.shape != (len(best_fit), 2) or errors.shape != best_fit.shape:
        raise ValueError("Impact best fit, errors and bounds have inconsistent shapes")

    lower_points = np.tile(best_fit, (len(best_fit), 1))
    upper_points = np.tile(best_fit, (len(best_fit), 1))
    valid_error = np.isfinite(errors) & (errors > 0.0)
    diagonal = np.arange(len(best_fit))
    lower_points[diagonal, diagonal] = np.where(valid_error, np.maximum(bounds[:, 0], best_fit - errors), best_fit)
    upper_points[diagonal, diagonal] = np.where(valid_error, np.minimum(bounds[:, 1], best_fit + errors), best_fit)
    shift_minus = np.diag(lower_points)
    shift_plus = np.diag(upper_points)
    valid_minus = valid_error & (shift_minus < best_fit)
    valid_plus = valid_error & (shift_plus > best_fit)
    valid_error = valid_minus | valid_plus
    return {
        "mode": "fixed",
        "minus_points": lower_points,
        "plus_points": upper_points,
        "shift_minus": shift_minus,
        "shift_plus": shift_plus,
        "valid_minus": valid_minus,
        "valid_plus": valid_plus,
        "valid_error": valid_error,
        "profile": None,
    }


# Compute physical Hessian shifts through the full decoder and covariance
def physical_parameter_shifts(*, transform, best_fit, covariance, bounds, mode, maxiter, value_gradient_func=None):
    propagated = transform.propagate(best_fit, covariance)
    values, errors = propagated["values"], propagated["errors"]
    valid = np.isfinite(errors) & (errors > 0.0)
    if mode == "profiled" and value_gradient_func is None:
        raise ValueError("Physical profiled impacts require an autograd likelihood gradient")
    shifts = {"mode": "physical_profiled" if mode == "profiled" else "physical_covariance",
              "physical_values": values, "physical_errors": errors, "valid_error": valid,
              "profile": {"minus": [None] * len(values), "plus": [None] * len(values)}}
    baseline = None if value_gradient_func is None else float(value_gradient_func(best_fit)[0])
    for side, sign in (("minus", -1.0), ("plus", 1.0)):
        shifts[f"{side}_points"] = np.tile(best_fit, (len(values), 1))
        shifts[f"shift_{side}"] = values.copy()
        shifts[f"valid_{side}"] = np.zeros(len(values), dtype=bool)
        for index in np.flatnonzero(valid):
            result = transform.constrain(best_fit, covariance, bounds, index, sign * errors[index], maxiter=maxiter,
                                         objective=value_gradient_func if mode == "profiled" else None)
            result["value"] = None if baseline is None else float(value_gradient_func(result["point"])[0])
            result["delta_z"] = None if baseline is None else result["value"] - baseline
            shifts["profile"][side][index] = result
            if not result["success"]:
                continue
            shifts[f"{side}_points"][index] = result["point"]
            shifts[f"shift_{side}"][index] += sign * errors[index]
            shifts[f"valid_{side}"][index] = True
    shifts["valid_error"] = shifts["valid_minus"] | shifts["valid_plus"]
    return shifts


# Compute one covariance-guided starting point for a constrained profile
def covariance_profile_start(
    *, best_fit: np.ndarray, covariance: np.ndarray, bounds: np.ndarray, fixed_index: int, fixed_value: float
) -> np.ndarray:
    start = np.asarray(best_fit, dtype=np.float64).copy()
    covariance = np.asarray(covariance, dtype=np.float64)
    variance = covariance[fixed_index, fixed_index] if covariance.shape == (len(best_fit), len(best_fit)) else np.nan
    if np.isfinite(variance) and variance > 0.0 and np.all(np.isfinite(covariance[:, fixed_index])):
        displacement = float(fixed_value - best_fit[fixed_index])
        start += covariance[:, fixed_index] * (displacement / variance)
    start = np.clip(start, bounds[:, 0], bounds[:, 1])
    start[fixed_index] = fixed_value
    return start


# Minimize the global surrogate with one physical parameter fixed
def profile_constrained_point(
    *,
    func_hat: Callable[[np.ndarray], np.ndarray],
    best_fit: np.ndarray,
    bounds: np.ndarray,
    fixed_index: int,
    fixed_value: float,
    initial: np.ndarray,
    maxiter: int,
    gradient_step: float = 1.0e-4,
    value_gradient_func: Callable | None = None,
    torch_config: dict | None = None,
) -> dict:
    point, value, diagnostics = profile.profile_minimize(
        func_hat=func_hat,
        fixed_indices=[fixed_index],
        fixed_values=[fixed_value],
        start=initial,
        bounds=bounds,
        fallback_starts=[best_fit],
        maxiter=maxiter,
        multistart=True,
        value_gradient_func=value_gradient_func,
        torch_config=torch_config,
        gradient_step=gradient_step,
        return_diagnostics=True,
    )
    return {"point": point, "value": value, **diagnostics}


# Profile all remaining parameters at every one-sided Hessian shift
def profiled_parameter_shifts(
    *,
    func_hat: Callable[[np.ndarray], np.ndarray],
    best_fit: np.ndarray,
    errors: np.ndarray,
    bounds: np.ndarray,
    covariance: np.ndarray,
    maxiter: int = 30,
    value_gradient_func: Callable | None = None,
    torch_config: dict | None = None,
    workers: int = 0,
) -> dict:
    if maxiter < 1:
        raise ValueError("Profiled impact maxiter must be at least 1")
    shifts = fixed_parameter_shifts(best_fit, errors, bounds)
    best_fit = np.asarray(best_fit, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)
    covariance = np.asarray(covariance, dtype=np.float64)
    baseline_z = float(np.asarray(func_hat(best_fit[None, :])).reshape(-1)[0])
    profiles = {"baseline_z": baseline_z, "minus": [None] * len(best_fit), "plus": [None] * len(best_fit)}

    jobs = []
    for index in range(len(best_fit)):
        for side, valid_key, point_key, shift_key in (
            ("minus", "valid_minus", "minus_points", "shift_minus"),
            ("plus", "valid_plus", "plus_points", "shift_plus"),
        ):
            if not bool(shifts[valid_key][index]):
                continue
            jobs.append((index, side, point_key, float(shifts[shift_key][index])))

    # Profile one independent parameter side in each CPU worker
    def evaluate_shift(job):
        index, side, point_key, fixed_value = job
        initial = covariance_profile_start(
            best_fit=best_fit, covariance=covariance, bounds=bounds, fixed_index=index, fixed_value=fixed_value
        )
        result = profile_constrained_point(
            func_hat=func_hat,
            best_fit=best_fit,
            bounds=bounds,
            fixed_index=index,
            fixed_value=fixed_value,
            initial=initial,
            maxiter=maxiter,
            value_gradient_func=value_gradient_func,
            torch_config=torch_config,
        )
        return index, side, point_key, result

    if torch_config is not None and str(torch_config["device"]).startswith("cuda"):
        workers = 1
    results = profile.profile_job_results(
        jobs, evaluate_shift, workers=workers, description="Profiled impact shifts"
    )
    for index, side, point_key, result in results:
        shifts[point_key][index] = result["point"]
        profiles[side][index] = {
            **result,
            "point": result["point"].tolist(),
            "delta_z": float(result["value"] - baseline_z),
        }
    shifts["mode"] = "profiled"
    shifts["profile"] = profiles
    return shifts


# Evaluate per-observable chi2 responses at prepared parameter shift points
def compute_parameter_impacts(
    *,
    model: dict,
    best_fit: np.ndarray,
    bounds: np.ndarray,
    feature_transform: Callable[[np.ndarray], np.ndarray],
    shifts: dict,
) -> dict:
    best_fit = np.asarray(best_fit, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)
    lower_points = np.asarray(shifts["minus_points"], dtype=np.float64)
    upper_points = np.asarray(shifts["plus_points"], dtype=np.float64)

    physical_points = np.vstack((best_fit[None, :], lower_points, upper_points))
    model_points = np.asarray(feature_transform(physical_points), dtype=np.float64)
    prediction = predict_observable_chi2(model, model_points)
    baseline_reference = nearest_realized_chi2(model, model_points[0])
    return assemble_parameter_impacts(
        prediction=prediction, best_fit=best_fit, bounds=bounds, shifts=shifts, baseline_reference=baseline_reference
    )


# Assemble normalized pulls and observable responses from evaluated shift points
def assemble_parameter_impacts(
    *,
    prediction: np.ndarray,
    best_fit: np.ndarray,
    bounds: np.ndarray,
    shifts: dict,
    baseline_reference: dict | None = None,
) -> dict:
    best_fit = np.asarray(shifts.get("physical_values", best_fit), dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)
    if "physical_values" in shifts:
        errors = shifts["physical_errors"]
        bounds = np.column_stack((best_fit - errors, best_fit + errors))
    prediction = np.asarray(prediction, dtype=np.float64)
    valid_minus = np.asarray(shifts["valid_minus"], dtype=bool)
    valid_plus = np.asarray(shifts["valid_plus"], dtype=bool)
    valid_error = np.asarray(shifts["valid_error"], dtype=bool)
    response_baseline = prediction[0]
    n_parameters = len(best_fit)
    delta_minus = (prediction[1 : 1 + n_parameters] - response_baseline).T
    delta_plus = (prediction[1 + n_parameters :] - response_baseline).T
    delta_minus[:, ~valid_minus] = np.nan
    delta_plus[:, ~valid_plus] = np.nan

    half_width = 0.5 * (bounds[:, 1] - bounds[:, 0])
    midpoint = 0.5 * (bounds[:, 1] + bounds[:, 0])
    normalized_best = np.divide(
        best_fit - midpoint, half_width, out=np.full_like(best_fit, np.nan), where=half_width > 0.0
    )
    normalized_error_minus = np.divide(
        best_fit - np.asarray(shifts["shift_minus"], dtype=np.float64),
        half_width,
        out=np.full_like(best_fit, np.nan),
        where=half_width > 0.0,
    )
    normalized_error_plus = np.divide(
        np.asarray(shifts["shift_plus"], dtype=np.float64) - best_fit,
        half_width,
        out=np.full_like(best_fit, np.nan),
        where=half_width > 0.0,
    )

    reference = baseline_reference or {
        "kind": "histogram_surrogate_optimum",
        "chi2": response_baseline,
        "trial_index": np.full(len(response_baseline), -1, dtype=np.int64),
        "feature_distance": np.full(len(response_baseline), np.nan, dtype=np.float64),
    }
    return {
        "baseline": np.asarray(reference["chi2"], dtype=np.float64),
        "baseline_reference": reference,
        "best_fit": best_fit,
        "mode": str(shifts["mode"]),
        "profile": shifts.get("profile"),
        "delta_minus": delta_minus,
        "delta_plus": delta_plus,
        "normalized_best": normalized_best,
        "normalized_error_minus": normalized_error_minus,
        "normalized_error_plus": normalized_error_plus,
        "shift_minus": np.asarray(shifts["shift_minus"], dtype=np.float64),
        "shift_plus": np.asarray(shifts["shift_plus"], dtype=np.float64),
        "valid_minus": valid_minus,
        "valid_plus": valid_plus,
        "valid_error": valid_error,
    }


# Evaluate direct histogram-surrogate impacts without fitting a second surrogate
def compute_direct_parameter_impacts(
    *,
    evaluate_observable_chi2: Callable[[np.ndarray], np.ndarray],
    best_fit: np.ndarray,
    bounds: np.ndarray,
    shifts: dict,
) -> dict:
    best_fit = np.asarray(best_fit, dtype=np.float64)
    physical_points = np.vstack((best_fit[None, :], shifts["minus_points"], shifts["plus_points"]))
    prediction = np.asarray(evaluate_observable_chi2(physical_points), dtype=np.float64)
    return assemble_parameter_impacts(prediction=prediction, best_fit=best_fit, bounds=bounds, shifts=shifts)


# Compute a filesystem-safe ASCII slug
def safe_slug(value: str) -> str:
    text = str(value).encode("ascii", errors="ignore").decode("ascii")
    text = re.sub(r"[^A-Za-z0-9._-]+", "-", text).strip("-")
    return text or "observable"


# Rank one observable by the largest absolute one-sided impact
def ranked_parameter_indices(impacts: dict, observable_index: int, top_n: int) -> np.ndarray:
    minus = np.asarray(impacts["delta_minus"][observable_index], dtype=np.float64)
    plus = np.asarray(impacts["delta_plus"][observable_index], dtype=np.float64)
    candidates = np.vstack((np.abs(minus), np.abs(plus)))
    magnitude = np.max(np.where(np.isfinite(candidates), candidates, -np.inf), axis=0)
    magnitude[~np.isfinite(magnitude)] = np.nan
    finite = np.flatnonzero(np.isfinite(magnitude))
    if finite.size == 0:
        return np.empty(0, dtype=np.int64)
    order = finite[np.argsort(magnitude[finite])[::-1]]
    return order[: max(1, int(top_n))]


# Compute the machine-readable definition for one impact mode
def impact_definition(mode: str) -> str:
    if mode == "physical_covariance":
        return "Histogram chi2 response at decoded Hessian shifts with minimum input covariance distance"
    if mode == "physical_profiled":
        return "Histogram chi2 response after likelihood refits constrained to decoded Hessian shifts"
    if mode == "fixed":
        return "Histogram chi2 response at one-at-a-time global Hessian shifts with other parameters fixed"
    if mode == "profiled":
        return "Histogram chi2 response at global Hessian shifts after constrained refits of all other parameters"
    raise ValueError(f"Unknown impact mode: {mode}")


# Compute the compact plot subtitle for one impact mode
def impact_subtitle(mode: str) -> str:
    if mode == "physical_covariance":
        return "Physical Hessian shifts with minimum covariance distance"
    if mode == "physical_profiled":
        return "Physical Hessian shifts with constrained likelihood refits"
    if mode == "fixed":
        return "One-at-a-time global Hessian shifts with other parameters fixed"
    if mode == "profiled":
        return "Global Hessian shifts with all remaining parameters reprofiled"
    raise ValueError(f"Unknown impact mode: {mode}")


# Append a math-typeset validation coefficient to one impact subtitle
def impact_validation_subtitle(mode: str, validation_r2: float | None) -> str:
    subtitle = impact_subtitle(mode)
    if validation_r2 is None:
        return subtitle
    return subtitle + f"; validation $R^2$ = {validation_r2:.3f}"


# Move impact axes left until rendered parameter labels reach a compact margin
def compact_impact_left_margin(fig, ax, *, target_fraction: float = 0.015) -> float:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    labels = [label for label in ax.get_yticklabels() if label.get_visible()]
    if not labels:
        return float(fig.subplotpars.left)
    minimum_x = min(label.get_window_extent(renderer).x0 for label in labels)
    target_x = float(target_fraction) * float(fig.bbox.width)
    correction = (minimum_x - target_x) / float(fig.bbox.width)
    new_left = float(np.clip(fig.subplotpars.left - correction, 0.08, 0.42))
    fig.subplots_adjust(left=new_left)
    return new_left


# Compute compact profile diagnostics without duplicating full parameter vectors
def compact_profile_record(profile: dict | None, side: str, index: int) -> dict | None:
    if profile is None or profile[side][index] is None:
        return None
    record = profile[side][index]
    return {
        "global_surrogate_z": finite_or_none(record["value"]),
        "delta_global_surrogate_z": finite_or_none(record["delta_z"]),
        "success": bool(record["success"]),
        "message": str(record["message"]),
        "iterations": int(record["iterations"]),
        "evaluations": int(record["evaluations"]),
    }


# Draw and persist one CMS Combine style observable impact plot
def plot_observable_impact(
    *,
    output_dir: pathlib.Path,
    observable: dict,
    observable_index: int,
    param_names: Sequence[str],
    impacts: dict,
    model: dict,
    top_n: int,
    plot_brand: str = "GRANIITTI",
) -> dict | None:
    indices = ranked_parameter_indices(impacts=impacts, observable_index=observable_index, top_n=top_n)
    if indices.size == 0:
        return None

    output_dir = pathlib.Path(output_dir)
    ensure_dir(output_dir)
    dataset = int(observable["dataset"])
    subset = int(observable["subset"])
    observable_name = str(observable["observable"])
    stem = f"dataset-{dataset:03d}_subset-{subset:03d}_{safe_slug(observable_name)}"
    y = np.arange(len(indices), dtype=np.float64)
    minus = np.asarray(impacts["delta_minus"][observable_index, indices], dtype=np.float64)
    plus = np.asarray(impacts["delta_plus"][observable_index, indices], dtype=np.float64)
    pull = np.asarray(impacts["normalized_best"][indices], dtype=np.float64)
    pull_minus = np.asarray(impacts["normalized_error_minus"][indices], dtype=np.float64)
    pull_plus = np.asarray(impacts["normalized_error_plus"][indices], dtype=np.float64)

    height = max(4.8, 0.36 * len(indices) + 1.8)
    fig, (ax_pull, ax_impact) = plt.subplots(
        1, 2, figsize=(13.0, height), sharey=True, gridspec_kw={"width_ratios": [1.0, 1.45], "wspace": 0.04}
    )
    for row in range(0, len(indices), 2):
        for axis in (ax_pull, ax_impact):
            axis.axhspan(row - 0.5, row + 0.5, color="#f3f3f3", zorder=0)

    ax_pull.axvspan(-1.0, 1.0, color="#d9d9d9", alpha=0.55, zorder=0)
    ax_pull.axvline(0.0, color="black", linewidth=0.8)
    ax_pull.errorbar(
        pull,
        y,
        xerr=np.vstack((pull_minus, pull_plus)),
        fmt="o",
        color="black",
        markersize=4.0,
        capsize=2.0,
        linewidth=1.0,
        zorder=4,
    )
    ax_pull.set_xlim(-1.55, 1.55)
    ax_pull.set_xlabel(r"$(p-\hat{p})/\sigma_p$" if impacts["mode"].startswith("physical_")
                       else r"$(\hat{\theta}-\theta_{\rm mid})/(\Delta\theta/2)$")
    ax_pull.set_yticks(y, labels=[param_names[index] for index in indices])
    ax_pull.tick_params(axis="y", labelsize=7)
    ax_pull.grid(axis="x", color="white", linewidth=0.8)

    bar_height = 0.32
    ax_impact.barh(
        y - 0.17, minus, height=bar_height, color=plot.colors(1), label=r"$\theta-\sigma_{\rm Hessian}$", zorder=3
    )
    ax_impact.barh(
        y + 0.17, plus, height=bar_height, color=plot.colors(0), label=r"$\theta+\sigma_{\rm Hessian}$", zorder=3
    )
    impact_limit = float(np.nanmax(np.abs(np.concatenate((minus, plus)))))
    impact_limit = 1.0 if not np.isfinite(impact_limit) or impact_limit == 0.0 else 1.18 * impact_limit
    ax_impact.set_xlim(-impact_limit, impact_limit)
    ax_impact.axvline(0.0, color="black", linewidth=0.8)
    ax_impact.set_xlabel(r"Impact on histogram $\chi^2$")
    ax_impact.grid(axis="x", color="#d0d0d0", linewidth=0.6, alpha=0.8)
    ax_impact.legend(loc="best", fontsize=8, frameon=False)
    ax_pull.invert_yaxis()

    title = f"dataset {dataset}, subset {subset}: {observable_name}"
    add_icetune_branding(fig, brand=plot_brand, x=0.015, y=0.985)
    fig.text(0.98, 0.985, title, fontsize=10, ha="right", va="top")
    validation_r2 = finite_or_none(model["validation_r2"][observable_index])
    fig.text(0.98, 0.958, impact_validation_subtitle(impacts["mode"], validation_r2), fontsize=7, ha="right", va="top")
    fig.subplots_adjust(left=0.30, right=0.98, top=0.94, bottom=0.09)
    compact_impact_left_margin(fig, ax_pull)

    png_path = output_dir / f"{stem}.png"
    pdf_path = output_dir / f"{stem}.pdf"
    fig.savefig(png_path, dpi=220)
    fig.savefig(pdf_path)
    plt.close(fig)

    parameter_records = []
    for index in indices:
        parameter_records.append(
            {
                "name": str(param_names[index]),
                "best_fit": finite_or_none(impacts["best_fit"][index]),
                "shift_minus": finite_or_none(impacts["shift_minus"][index]),
                "shift_plus": finite_or_none(impacts["shift_plus"][index]),
                "delta_chi2_minus": finite_or_none(impacts["delta_minus"][observable_index, index]),
                "delta_chi2_plus": finite_or_none(impacts["delta_plus"][observable_index, index]),
                "profile_minus": compact_profile_record(impacts["profile"], "minus", int(index)),
                "profile_plus": compact_profile_record(impacts["profile"], "plus", int(index)),
            }
        )
    json_path = output_dir / f"{stem}.json"
    payload = {
        "impact_schema_version": 1,
        "mode": impacts["mode"],
        "plot_brand": str(plot_brand),
        "definition": impact_definition(impacts["mode"]),
        "dataset": dataset,
        "subset": subset,
        "observable": observable_name,
        "ndf": finite_or_none(observable.get("ndf")),
        "fitw": finite_or_none(observable.get("fitw")),
        "baseline_chi2": finite_or_none(impacts["baseline"][observable_index]),
        "baseline_reference": {
            "kind": str(impacts["baseline_reference"]["kind"]),
            "trial_index": (
                int(impacts["baseline_reference"]["trial_index"][observable_index])
                if int(impacts["baseline_reference"]["trial_index"][observable_index]) >= 0
                else None
            ),
            "feature_distance": finite_or_none(impacts["baseline_reference"]["feature_distance"][observable_index]),
        },
        "validation_r2": validation_r2,
        "validation_mae": finite_or_none(model["validation_mae"][observable_index]),
        "trial_count": int(model["trial_count"][observable_index]),
        "parameters": parameter_records,
        "png": str(png_path),
        "pdf": str(pdf_path),
    }
    write_json_file(json_path, json_safe(payload), indent=4)
    return {
        "dataset": dataset,
        "subset": subset,
        "observable": observable_name,
        "png": str(png_path),
        "pdf": str(pdf_path),
        "json": str(json_path),
    }


# Render all trained observable impacts and write their common manifest
def _render_impact_outputs(
    *,
    impacts: dict,
    model: dict,
    observables: Sequence[dict],
    param_names: Sequence[str],
    output_dir: pathlib.Path,
    top_n: int,
    plot_brand: str,
    surrogate: dict,
    physical: dict | None = None,
) -> dict:
    physical_output = None
    if physical is not None:
        options = {key: value for key, value in physical.items() if key not in {"evaluate", "output_dir"}}
        shifts = physical_parameter_shifts(**options)
        transformed = compute_direct_parameter_impacts(
            evaluate_observable_chi2=physical["evaluate"], best_fit=physical["best_fit"],
            bounds=physical["bounds"], shifts=shifts,
        )
        physical_output = _render_impact_outputs(
            impacts=transformed, model=model, observables=observables, param_names=physical["transform"].names,
            output_dir=physical["output_dir"], top_n=top_n, plot_brand=plot_brand, surrogate=surrogate,
        )
    trained = np.asarray(model["trained"], dtype=bool)
    outputs = []
    for observable_index, observable in enumerate(observables):
        if not bool(trained[observable_index]):
            continue
        output = plot_observable_impact(
            output_dir=output_dir,
            observable=observable,
            observable_index=observable_index,
            param_names=param_names,
            impacts=impacts,
            model=model,
            top_n=top_n,
            plot_brand=plot_brand,
        )
        if output is not None:
            outputs.append(output)

    output_dir = pathlib.Path(output_dir)
    ensure_dir(output_dir)
    profile_records = (
        [record for side in ("minus", "plus") for record in impacts["profile"][side] if record is not None]
        if impacts["profile"] is not None
        else []
    )
    manifest_path = output_dir / "manifest.json"
    manifest = {
        "impact_schema_version": 1,
        "mode": impacts["mode"],
        "plot_brand": str(plot_brand),
        "definition": impact_definition(impacts["mode"]),
        "observable_count": len(observables),
        "trained_observable_count": int(np.count_nonzero(trained)),
        "published_observable_count": len(outputs),
        "parameter_count": len(param_names),
        "valid_parameter_error_count": int(np.count_nonzero(impacts["valid_error"])),
        "profiled_shift_count": len(profile_records),
        "successful_profiled_shift_count": int(sum(bool(record["success"]) for record in profile_records)),
        "profile": impacts["profile"],
        "surrogate": surrogate,
        "outputs": outputs,
    }
    write_json_file(manifest_path, json_safe(manifest), indent=4)
    return {**manifest, "manifest": str(manifest_path), **({"physical": physical_output} if physical is not None else {})}


# Fit observable surrogates and render every available histogram impact plot
def render_parameter_impacts(
    *,
    X_model: np.ndarray,
    observable_chi2: np.ndarray,
    observables: Sequence[dict],
    param_names: Sequence[str],
    best_fit: np.ndarray,
    errors: np.ndarray,
    bounds: np.ndarray,
    feature_transform: Callable[[np.ndarray], np.ndarray],
    output_dir: pathlib.Path,
    shifts: dict | None = None,
    n_features: int = 256,
    scale: float | None = None,
    alpha: float = 1.0e-3,
    validation_fraction: float = 0.2,
    rngseed: int = 1234,
    min_trials: int = 30,
    top_n: int = 20,
    plot_brand: str = "GRANIITTI",
    physical: dict | None = None,
) -> dict:
    model = fit_observable_surrogate(
        X_model=X_model,
        observable_chi2=observable_chi2,
        n_features=n_features,
        scale=scale,
        alpha=alpha,
        validation_fraction=validation_fraction,
        rngseed=rngseed,
        min_trials=min_trials,
    )
    if shifts is None:
        shifts = fixed_parameter_shifts(best_fit=best_fit, errors=errors, bounds=bounds)
    impacts = compute_parameter_impacts(
        model=model, best_fit=best_fit, bounds=bounds, feature_transform=feature_transform, shifts=shifts
    )

    # Evaluate the same observable surrogate at decoded physical shift points
    def evaluate(points):
        return predict_observable_chi2(model, feature_transform(points))

    return _render_impact_outputs(
        impacts=impacts,
        model=model,
        observables=observables,
        param_names=param_names,
        output_dir=output_dir,
        top_n=top_n,
        plot_brand=plot_brand,
        physical=None if physical is None else {**physical, "evaluate": evaluate},
        surrogate={
            "kind": "rff_ridge",
            "target_transform": model["target_transform"],
            "impact_response": "surrogate_differences",
            "baseline_reference": "nearest_realized_trial",
            "features": int(n_features),
            "scale": float(model["rff"]["scale"]),
            "alpha": float(alpha),
            "validation_fraction": float(validation_fraction),
            "min_trials": int(min_trials),
        },
    )


# Render direct impacts from a histogram surrogate already used by the fit
def render_direct_parameter_impacts(
    *,
    evaluate_observable_chi2: Callable[[np.ndarray], np.ndarray],
    observables: Sequence[dict],
    param_names: Sequence[str],
    best_fit: np.ndarray,
    bounds: np.ndarray,
    shifts: dict,
    output_dir: pathlib.Path,
    validation_r2: np.ndarray | None = None,
    validation_mae: np.ndarray | None = None,
    trial_count: int = 0,
    top_n: int = 20,
    plot_brand: str = "GRANIITTI",
    physical: dict | None = None,
) -> dict:
    impacts = compute_direct_parameter_impacts(
        evaluate_observable_chi2=evaluate_observable_chi2, best_fit=best_fit, bounds=bounds, shifts=shifts
    )
    n_observables = len(observables)
    model = {
        "trained": np.ones(n_observables, dtype=bool),
        "validation_r2": (
            np.full(n_observables, np.nan, dtype=np.float64)
            if validation_r2 is None
            else np.asarray(validation_r2, dtype=np.float64)
        ),
        "validation_mae": (
            np.full(n_observables, np.nan, dtype=np.float64)
            if validation_mae is None
            else np.asarray(validation_mae, dtype=np.float64)
        ),
        "trial_count": np.full(n_observables, int(trial_count), dtype=np.int64),
    }
    return _render_impact_outputs(
        impacts=impacts,
        model=model,
        observables=observables,
        param_names=param_names,
        output_dir=output_dir,
        top_n=top_n,
        plot_brand=plot_brand,
        surrogate={"kind": "direct_histogram_replica"},
        physical=None if physical is None else {**physical, "evaluate": evaluate_observable_chi2},
    )

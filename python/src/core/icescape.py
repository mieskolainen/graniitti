# Fit, optimize and plot an icetune likelihood surrogate
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import json
import pathlib
import pickle
import sys
from functools import partial

import numpy as np
import torch
from termcolor import cprint

from core import __AUTHOR__, __RELEASE__, __version__
from core.inference import fit, gp, impact, neural, surrogate_cli, surrogate_data
from core.io.files import ensure_dir
from core.io.serialize import finite_or_none, json_safe
from core.stats import cov, objective
from core.tune import history as icetune_history
from core.tune import likelihood as icetune_likelihood
from core.tune.parameters import topology as parameter_topology
from core.tune.parameters.kernel import ProductManifoldGeometry


# Argument Parsing and Data Reading
# Parse command-line arguments for icescape
def parse_args():
    return surrogate_cli.parse_args(
        description=f"%(prog)s {__version__} {__RELEASE__} [{__AUTHOR__}]",
        replica=False,
    )


# Re-export shared topology and derivative machinery for the CLI module
_indexed_parameter_topology = parameter_topology.indexed_parameter_topology
model_feature_names = parameter_topology.model_feature_names
model_features = parameter_topology.model_features
model_features_torch = parameter_topology.model_features_torch
build_torch_derivative_callbacks = fit.build_torch_derivative_callbacks


def columnize(x: np.ndarray) -> np.ndarray:
    """Return a one-column array for single-output surrogate training"""

    return np.asarray(x, dtype=np.float64).reshape(-1, 1)


# Re-export the shared validation splitter
training_validation_indices = fit.training_validation_indices


# Fit a positivity-preserving GP target transform from training rows
def fit_gp_target_transform(values: np.ndarray, training_indices: np.ndarray, mode: str) -> dict:
    values = columnize(values)
    requested = str(mode)
    kind = (
        "log"
        if requested == "log" or (requested == "auto" and np.all(values > 0.0))
        else "standard"
    )
    offset = 0.0
    transformed = values.copy()
    if kind == "log":
        minimum = float(np.min(values))
        offset = max(0.0, 1.0 - minimum)
        transformed = np.log(values + offset)
    location = float(np.mean(transformed[training_indices]))
    scale = float(np.std(transformed[training_indices]))
    if not np.isfinite(scale) or scale <= 0.0:
        scale = 1.0
    return {
        "requested": requested,
        "kind": kind,
        "offset": offset,
        "location": location,
        "scale": scale,
        "uncertainty": (
            "delta_log_mean" if kind == "log" else "gaussian"
        ),
    }


# Transform objective values and standard deviations into GP target units
def transform_gp_targets(
    values: np.ndarray, errors: np.ndarray, transform: dict
) -> tuple[np.ndarray, np.ndarray]:
    values = columnize(values)
    errors = columnize(errors)
    if transform["kind"] == "log":
        shifted = values + float(transform["offset"])
        transformed = np.log(shifted)
        transformed_errors = np.divide(
            errors,
            shifted,
            out=np.zeros_like(errors),
            where=shifted > 0.0,
        )
    else:
        transformed = values
        transformed_errors = errors
    scale = float(transform["scale"])
    centered = (transformed - float(transform["location"])) / scale
    return centered, transformed_errors / scale


# Map standardized GP latent moments back into objective units
def inverse_gp_prediction(
    mean_scaled: torch.Tensor, variance_scaled: torch.Tensor | None, transform: dict
) -> tuple[torch.Tensor, torch.Tensor | None]:
    scale = float(transform["scale"])
    latent_mean = mean_scaled * scale + float(transform["location"])
    latent_variance = variance_scaled * scale**2 if variance_scaled is not None else None
    if transform["kind"] == "log":
        max_log = float(np.log(torch.finfo(latent_mean.dtype).max)) - 2.0
        shifted_median = torch.exp(torch.clamp(latent_mean, max=max_log))
        objective_mean = shifted_median - float(transform["offset"])
        objective_variance = (
            shifted_median.square() * latent_variance if latent_variance is not None else None
        )
        return objective_mean, objective_variance
    return latent_mean, latent_variance


# Compute objective-space GP holdout residual diagnostics
def gp_holdout_metrics(
    model, X_validation: np.ndarray, Z_validation: np.ndarray, transform: dict
) -> dict:
    X_tensor = torch.as_tensor(
        np.asarray(X_validation),
        dtype=model.dtype,
        device=model.device,
    )
    with torch.no_grad():
        mean_scaled = model.predict_mean(X_tensor)
        prediction, _ = inverse_gp_prediction(
            mean_scaled=mean_scaled,
            variance_scaled=None,
            transform=transform,
        )
    target = np.asarray(Z_validation, dtype=np.float64).reshape(-1)
    predicted = prediction.detach().cpu().numpy().reshape(-1)
    residual = predicted - target
    low_threshold = float(np.percentile(target, 20.0))
    low_mask = target <= low_threshold
    relative = np.divide(
        np.abs(residual),
        np.abs(target),
        out=np.full_like(target, np.nan),
        where=target != 0.0,
    )
    return {
        "count": len(target),
        "mae": float(np.mean(np.abs(residual))),
        "median_absolute_relative_error": float(np.nanmedian(relative)),
        "low_20_percent_threshold": low_threshold,
        "low_20_percent_count": int(np.count_nonzero(low_mask)),
        "low_20_percent_mae": float(np.mean(np.abs(residual[low_mask]))),
        "minimum_prediction": float(np.min(predicted)),
    }


# Re-export the shared Torch runtime configuration
configure_torch_runtime = fit.configure_torch_runtime


# Build NumPy evaluation and exact derivatives around one Torch model-space predictor
def physical_surrogate_callbacks(
    *,
    predict_model_torch,
    physical_features,
    physical_features_torch,
    device,
    dtype,
):
    # Evaluate physical Torch coordinates without breaking the autograd graph
    def predict_physical_torch(values):
        model_values = physical_features_torch(values)
        if model_values.ndim == 1:
            model_values = model_values.unsqueeze(0)
        return predict_model_torch(model_values).reshape(-1)

    # Evaluate physical NumPy coordinates without retaining an autograd graph
    def predict_physical(values):
        model_values = physical_features(np.atleast_2d(np.asarray(values, dtype=np.float64)))
        with torch.inference_mode():
            prediction = predict_model_torch(
                torch.as_tensor(model_values, dtype=dtype, device=device)
            )
        return prediction.detach().cpu().numpy().reshape(-1)

    # Compute the scalar objective required by the shared derivative builder
    def objective(values):
        return predict_physical_torch(values)[0]

    # Compute one objective value per row for batched profile minimization
    def objective_batch(values):
        return predict_physical_torch(values)

    gradient, hessian = build_torch_derivative_callbacks(
        objective=objective,
        device=device,
        dtype=dtype,
    )
    return predict_physical, gradient, hessian, objective, objective_batch


# Update a JSON object on disk while preserving existing keys
def update_json_file(path: pathlib.Path, updates: dict) -> None:
    path = pathlib.Path(path)
    payload = {}
    if path.is_file():
        with open(path, encoding="utf-8") as handle:
            payload = json.load(handle)
    payload.update(json_safe(updates))
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(json_safe(payload), handle, indent=4)


# Print and persist GP posterior uncertainty for the optimized surrogate Z value
def report_surrogate_z_uncertainty(
    *,
    output_root: pathlib.Path,
    surrogate_best: np.ndarray,
    surrogate_Z: float,
    z_uncertainty_func: callable | None,
) -> None:
    if z_uncertainty_func is None:
        return
    mean_z, sigma_z = z_uncertainty_func(np.asarray(surrogate_best, dtype=np.float64))
    mean_z = finite_or_none(mean_z)
    sigma_z = finite_or_none(sigma_z)
    if mean_z is None or sigma_z is None:
        return
    fit.print_table(
        "GP predictive uncertainty at surrogate optimum:",
        ["Quantity", "Value"],
        [
            ["surrogate optimum Z", f"{float(surrogate_Z):0.6g}"],
            ["GP predictive central Z_hat", f"{mean_z:0.6g}"],
            ["GP predictive sigma_Z", f"{sigma_z:0.6g}"],
            ["Z_hat +- sigma_Z", f"{mean_z:0.6g} +- {sigma_z:0.6g}"],
        ],
    )
    update_json_file(
        pathlib.Path(output_root) / "parameters.json",
        {
            "surrogate_Z_predictive": {
                "central": mean_z,
                "sigma": sigma_z,
                "source": "gp_posterior_at_minimum",
            }
        },
    )


# Load the global and per-observable history surfaces with aligned trial rows
def load_history_surface(args) -> dict:
    """Load the posterior-ready likelihood surface from schema-v1 history.json"""

    source = pathlib.Path(args.history_path or args.run_name)
    history_path = icetune_history.resolve_history_path(source, pathlib.Path(args.cdir))
    if history_path.is_dir():
        raise ValueError("icescape requires a compact history.json input")
    cprint(f"Reading likelihood history: {history_path}", "yellow")
    payload = icetune_likelihood.load_likelihood_history(history_path)
    X = np.asarray(payload["X"], dtype=np.float64)
    X_unit = np.asarray(payload["X_unit"], dtype=np.float64)
    logL = np.asarray(payload["logL"], dtype=np.float64)
    two_nll = np.asarray(payload["two_nll"], dtype=np.float64)
    valid = np.asarray(payload["valid"], dtype=bool)

    valid = (
        valid
        & np.all(np.isfinite(X), axis=1)
        & np.all(np.isfinite(X_unit), axis=1)
        & np.isfinite(logL)
        & np.isfinite(two_nll)
    )
    X = X[valid]
    X_unit = X_unit[valid]
    logL = logL[valid]
    two_nll = two_nll[valid]
    mc_uncertainty = [
        item for item, keep in zip(payload["mc_uncertainty"], valid, strict=False) if keep
    ]
    theta_hash = [item for item, keep in zip(payload["theta_hash"], valid, strict=False) if keep]
    observable_chi2 = np.asarray(payload["observable_chi2"], dtype=np.float64)[valid]

    if args.max_trials is not None:
        X = X[: args.max_trials]
        X_unit = X_unit[: args.max_trials]
        logL = logL[: args.max_trials]
        two_nll = two_nll[: args.max_trials]
        mc_uncertainty = mc_uncertainty[: args.max_trials]
        theta_hash = theta_hash[: args.max_trials]
        observable_chi2 = observable_chi2[: args.max_trials]

    sigma_two_nll = np.array(
        [
            np.nan if item is None else finite_or_none(item.get("sigma_two_nll"))
            for item in mc_uncertainty
        ],
        dtype=np.float64,
    )
    Z_err = np.where(np.isfinite(sigma_two_nll), sigma_two_nll, 0.0)
    Z = two_nll
    bounds = fit.parameter_bounds(X, payload["parameter_space"])
    topology = payload.get("parameter_topology", {})
    X_model = model_features(X, bounds, payload["param_names"], topology)

    return {
        "input_mode": "history",
        "history_path": history_path,
        "X": X,
        "X_model": X_model,
        "Z": Z,
        "Z_err": Z_err,
        "logL": logL,
        "param_names": payload["param_names"],
        "bounds": bounds,
        "parameter_space": payload["parameter_space"],
        "optimization": payload["optimization"],
        "simdriver": payload["simdriver"],
        "mc_steer": payload["mc_steer"],
        "plot_brand": payload["plot_brand"],
        "covariance_mode": payload["covariance_mode"],
        "parameter_topology": topology,
        "model_feature_names": model_feature_names(payload["param_names"], topology),
        "theta_hash": theta_hash,
        "observable_chi2": observable_chi2,
        "observables": payload["observables"],
    }


# Build per-observable chi-square surfaces from full histogram tensors
def full_observable_chi2(
    *,
    data: np.ndarray,
    data_error: np.ndarray,
    mc: np.ndarray,
    mc_error: np.ndarray,
    manifest: list[dict],
    covariance_payload: dict | None,
) -> tuple[np.ndarray, list[dict]]:
    """Return aligned fit-weighted marginal histogram chi-square values"""
    data = np.asarray(data, dtype=np.float64)
    data_error = np.asarray(data_error, dtype=np.float64)
    mc = np.asarray(mc, dtype=np.float64)
    mc_error = np.asarray(mc_error, dtype=np.float64)
    manifest_by_key = {
        (int(item["dataset"]), int(item["subset"]), str(item["observable"])): item
        for item in manifest
    }
    values = []
    observables = []
    records = covariance_payload["layout"] if covariance_payload is not None else manifest
    for record in records:
        key = (
            int(record["dataset"]),
            int(record["subset"]),
            str(record["observable"]),
        )
        item = manifest_by_key[key]
        fitw = float(item.get("fitw", 1.0))
        if covariance_payload is None:
            indices = np.arange(int(item["start"]), int(item["stop"]), dtype=np.int64)
            valid = np.asarray(item["valid"], dtype=bool)
            active = valid[None, :] & ((data[indices] != 0.0) | (mc[:, indices] != 0.0))
            chi2, _ = objective.global_histogram_chi2(
                Y_data=data[indices],
                E_data=data_error[indices],
                Y_hat=mc[:, indices],
                E_hat=mc_error[:, indices],
                mask=active,
                fit_weights=np.full(len(indices), fitw, dtype=np.float64),
            )
            ndf = int(np.count_nonzero(valid))
        else:
            bins, covariance = cov.observable_covariance_selection(
                payload=covariance_payload,
                dataset=key[0],
                subset=key[1],
                observable=key[2],
                bin_count=int(item["stop"]) - int(item["start"]),
            )
            indices = int(item["start"]) + bins
            chi2, _ = objective.correlated_chi2_batch(
                data_values=data[indices],
                mc_prediction=mc[:, indices],
                mc_stat_uncertainty=mc_error[:, indices],
                fit_weights=np.full(len(indices), fitw, dtype=np.float64),
                data_total_covariance=covariance,
            )
            ndf = len(indices)
        values.append(np.asarray(chi2, dtype=np.float64))
        observables.append(
            {
                "dataset": key[0],
                "subset": key[1],
                "observable": key[2],
                "ndf": ndf,
                "fitw": fitw,
            }
        )
    matrix = np.column_stack(values) if values else np.empty((len(mc), 0), dtype=np.float64)
    return matrix, observables


# Load the explicit histogram payload surface without compact impact records
def load_full_surface(args) -> dict:
    """Load full histogram payloads and recompute global 2NLL"""

    cprint("Using full histogram input path", "yellow")
    surface = surrogate_data.load_histogram_surface(args, active_nonzero_bins=False)
    observable_chi2, observables = full_observable_chi2(
        data=surface["Y_data"],
        data_error=surface["E_data"],
        mc=surface["Y_hat"],
        mc_error=surface["E_hat"],
        manifest=surface["histogram_manifest"],
        covariance_payload=surface["covariance_payload"],
    )
    surface.update(
        {
            "input_mode": "full",
            "history_path": None,
            "logL": -0.5 * surface["Z"],
            "parameter_space": [],
            "optimization": {},
            "theta_hash": [],
            "observable_chi2": observable_chi2,
            "observables": observables,
        }
    )
    return surface


def load_surface(args) -> dict:
    """Load either schema-v1 history input or the explicit full histogram input"""

    if args.input_mode == "full":
        return load_full_surface(args)
    return load_history_surface(args)


# Render observable-level chi2 impacts using fixed or profiled Hessian shifts
def render_surrogate_impacts(
    *,
    surface: dict,
    X_model: np.ndarray,
    observable_chi2: np.ndarray,
    param_names: list[str],
    bounds: np.ndarray,
    parameter_topology: dict,
    fit_result,
    surrogate_best: np.ndarray,
    func_hat: callable,
    value_gradient_func: callable | None,
    torch_config: dict | None,
    model_name: str,
    args,
) -> dict:
    if observable_chi2.ndim != 2 or observable_chi2.shape[1] == 0:
        raise RuntimeError(
            "Impact rendering requires per-observable chi2 records from the selected input"
        )

    covariance, _, errors, _ = fit.fitted_covariance_summary(
        m=fit_result,
        param_names=param_names,
        surrogate_vals=np.asarray(surrogate_best, dtype=np.float64),
    )
    if not np.any(np.isfinite(errors) & (errors > 0.0)):
        raise RuntimeError(
            "Cannot construct impacts because the global Hessian covariance is unavailable"
        )

    # Preserve phase and constrained-amplitude topology in every parameter shift
    def feature_transform(X_physical):
        return model_features(
            np.asarray(X_physical, dtype=np.float64),
            bounds,
            param_names,
            parameter_topology,
        )

    if args.impact_mode == "fixed":
        shifts = impact.fixed_parameter_shifts(
            best_fit=surrogate_best,
            errors=errors,
            bounds=bounds,
        )
    else:
        cprint(
            "Profiling all remaining parameters at every Hessian impact shift",
            "yellow",
        )
        shifts = impact.profiled_parameter_shifts(
            func_hat=func_hat,
            best_fit=surrogate_best,
            errors=errors,
            bounds=bounds,
            covariance=covariance,
            maxiter=args.impact_profile_maxiter,
            value_gradient_func=value_gradient_func,
            torch_config=torch_config,
            workers=args.profile_workers,
        )

    output_dir = pathlib.Path(args.output_root) / model_name / "optimizer" / "impacts" / args.impact_mode
    cprint(
        f"Rendering per-histogram parameter impacts under {output_dir}",
        "yellow",
    )
    result = impact.render_parameter_impacts(
        X_model=X_model,
        observable_chi2=observable_chi2,
        observables=surface["observables"],
        param_names=param_names,
        best_fit=surrogate_best,
        errors=errors,
        bounds=bounds,
        feature_transform=feature_transform,
        output_dir=output_dir,
        shifts=shifts,
        n_features=args.impact_features,
        physical={"transform": surrogate_cli.parameter_transform(args, surface), "best_fit": surrogate_best,
                  "covariance": covariance, "bounds": bounds, "mode": args.impact_mode,
                  "output_dir": pathlib.Path(args.output_root) / model_name / "physical" / "impacts" / args.impact_mode,
                  "maxiter": args.impact_profile_maxiter, "value_gradient_func": value_gradient_func},
        scale=args.impact_scale,
        alpha=args.impact_ridge,
        validation_fraction=args.validation_fraction,
        rngseed=args.rngseed,
        min_trials=args.impact_min_trials,
        top_n=args.impact_top,
        plot_brand=args.plot_brand,
    )
    fit.print_table(
        "Histogram impact summary:",
        ["Quantity", "Value"],
        [
            ["mode", result["mode"]],
            ["history observables", result["observable_count"]],
            ["trained observables", result["trained_observable_count"]],
            ["published plots", result["published_observable_count"]],
            ["valid Hessian parameter errors", result["valid_parameter_error_count"]],
            ["profiled shifts", result["profiled_shift_count"]],
            ["successful profiled shifts", result["successful_profiled_shift_count"]],
            ["manifest", result["manifest"]],
        ],
        color="green",
    )
    return result


# Fit one likelihood surrogate and render every requested post-fit diagnostic
def fit_likelihood_surrogate(
    *,
    surface: dict,
    model_name: str,
    func_hat: callable,
    realized_best: np.ndarray,
    realized_Z: float,
    value_gradient_func: callable | None,
    hessian_func: callable | None,
    objective_torch: callable,
    objective_batch_torch: callable,
    torch_device,
    torch_dtype,
    args,
    z_uncertainty_func: callable | None = None,
):
    output_root = pathlib.Path(args.output_root) / model_name
    torch_config = surrogate_cli.torch_minimizer_config(
        args=args,
        objective_torch=objective_torch,
        objective_batch_torch=objective_batch_torch,
        device=torch_device,
        dtype=torch_dtype,
    )
    m, surrogate_best = surrogate_cli.optimize_surrogate(
        args=args,
        X=surface["X"],
        Z=surface["Z"],
        param_names=surface["param_names"],
        func_hat=func_hat,
        bounds=surface["bounds"],
        output_root=output_root,
        realized_best=realized_best,
        realized_Z=realized_Z,
        objective_torch=objective_torch,
        objective_batch_torch=objective_batch_torch,
        device=torch_device,
        dtype=torch_dtype,
        transform=surrogate_cli.parameter_transform(args, surface),
    )
    surrogate_Z = float(m.fmin.fval)
    report_surrogate_z_uncertainty(
        output_root=output_root,
        surrogate_best=surrogate_best,
        surrogate_Z=surrogate_Z,
        z_uncertainty_func=z_uncertainty_func,
    )
    surrogate_cli.render_profiles(
        args=args,
        X=surface["X"],
        Z=surface["Z"],
        param_names=surface["param_names"],
        func_hat=func_hat,
        output_root=output_root,
        realized_best=realized_best,
        surrogate_best=surrogate_best,
        realized_Z=realized_Z,
        surrogate_Z=surrogate_Z,
        bounds=surface["bounds"],
        value_gradient_func=value_gradient_func,
        torch_config=torch_config,
        transform=surrogate_cli.parameter_transform(args, surface),
    )
    render_surrogate_impacts(
        surface=surface,
        X_model=surface["X_model"],
        observable_chi2=np.asarray(surface["observable_chi2"], dtype=np.float64),
        param_names=surface["param_names"],
        bounds=surface["bounds"],
        parameter_topology=surface["parameter_topology"],
        fit_result=m,
        surrogate_best=surrogate_best,
        func_hat=func_hat,
        value_gradient_func=value_gradient_func,
        torch_config=torch_config,
        model_name=model_name,
        args=args,
    )
    surrogate_cli.run_posterior(
        args=args,
        objective_torch=objective_torch,
        device=torch_device,
        dtype=torch_dtype,
        bounds=surface["bounds"],
        param_names=surface["param_names"],
        output_root=output_root,
        start=surrogate_best,
        transform=surrogate_cli.parameter_transform(args, surface),
    )
    return m, surrogate_best


# Re-export the shared histogram chi2 convention
global_2NLL = objective.global_histogram_chi2
# Main Routine
# Run the icescape command-line workflow
def main():
    args = parse_args()
    start_time = surrogate_cli.start_run(
        args,
        tool="icescape",
        version=__version__,
        command=sys.argv,
    )

    surface = load_surface(args)
    surrogate_cli.record_surface(
        args,
        surface,
        inputs={
            "mode": surface["input_mode"],
            "history_path": surface["history_path"],
            "simdriver": surface.get("simdriver"),
            "covariance_mode": surface.get("covariance_mode"),
        },
    )
    surface, _ = surrogate_cli.quality_filter_surface(
        surface,
        finite_keys=("X", "X_model", "Z"),
        row_keys=("X", "X_model", "Z", "Z_err", "observable_chi2"),
        args=args,
    )
    param_names = surface["param_names"]
    X = surface["X"]
    X_model = surface["X_model"]
    Z = surface["Z"]
    Z_err = surface["Z_err"]
    bounds = surface["bounds"]
    parameter_topology = surface["parameter_topology"]
    observable_chi2 = np.asarray(surface["observable_chi2"], dtype=np.float64)
    physical_features = partial(
        model_features,
        bounds=bounds,
        param_names=param_names,
        topology=parameter_topology,
    )
    physical_features_torch = partial(
        model_features_torch,
        bounds=bounds,
        param_names=param_names,
        topology=parameter_topology,
    )

    Z_err = np.where(np.isfinite(Z_err), Z_err, 0.0)

    min_idx = np.argmin(Z)
    realized_best = X[min_idx].copy()
    realized_Z = Z[min_idx].copy()

    # Save the trial data as a pickle file
    pickle_filename = pathlib.Path(args.output_root) / "icescape.pkl"
    ensure_dir(pickle_filename.parent)
    with open(pickle_filename, "wb") as f:
        pickle.dump(
            {
                "X": X,
                "X_model": X_model,
                "Z": Z,
                "Z_err": Z_err,
                "param_names": param_names,
                "bounds": bounds,
                "parameter_space": surface["parameter_space"],
                "optimization": surface["optimization"],
                "parameter_topology": parameter_topology,
                "model_feature_names": surface["model_feature_names"],
                "input_mode": surface["input_mode"],
                "history_path": surface["history_path"],
                "simdriver": surface.get("simdriver"),
                "plot_brand": args.plot_brand,
                "observable_chi2": observable_chi2,
                "observables": surface["observables"],
            },
            f,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
        fit.print_table(
            "Saved outputs:", ["Output", "Path"], [["icescape pickle", pickle_filename]]
        )
    # Split data into training and validation sets
    Z_fit = columnize(Z)
    training_indices, validation_indices = training_validation_indices(
        n_rows=len(Z),
        validation_fraction=args.validation_fraction,
        rngseed=args.rngseed,
        required_training=[int(np.argmin(Z))],
    )
    Z_scaled, Z_mu, Z_std = fit.standardize(Z_fit, training_indices)
    Z_mu_scalar = Z_mu.item()
    Z_std_scalar = Z_std.item()
    Z_err_scaled = columnize(Z_err) / Z_std
    X_train = X_model[training_indices]
    X_val = X_model[validation_indices]
    y_train = Z_scaled[training_indices]
    y_val = Z_scaled[validation_indices]
    err_train = Z_err_scaled[training_indices]
    err_val = Z_err_scaled[validation_indices]

    fit.print_table(
        "Training split:",
        ["Array", "Shape"],
        [
            ["X_train", X_train.shape],
            ["X_val", X_val.shape],
            ["y_train", y_train.shape],
            ["y_val", y_val.shape],
            ["realized minimum in training", bool(np.argmin(Z) in training_indices)],
            ["target mean from training", Z_mu_scalar],
            ["target std from training", Z_std_scalar],
        ],
    )
    in_dim = X_train.shape[1]
    out_dim = y_train.shape[1]
    # Fit Gaussian Process surrogate

    if args.model == "gp":
        # Convert one numerical array into the model-selected GP precision
        def gp_tensor(values):
            return torch.as_tensor(np.asarray(values))

        gp_target_transform = fit_gp_target_transform(
            values=Z_fit,
            training_indices=training_indices,
            mode=args.gp_target_transform,
        )
        gp_Z_scaled, gp_Z_err_scaled = transform_gp_targets(
            values=Z_fit,
            errors=columnize(Z_err),
            transform=gp_target_transform,
        )
        fit.print_table(
            "GP target transform:",
            ["Quantity", "Value"],
            [
                ["requested", gp_target_transform["requested"]],
                ["selected", gp_target_transform["kind"]],
                ["offset", gp_target_transform["offset"]],
                ["transformed location", gp_target_transform["location"]],
                ["transformed scale", gp_target_transform["scale"]],
                ["uncertainty propagation", gp_target_transform["uncertainty"]],
            ],
        )
        run_gp_optimization = args.gp_optimize and args.gp_opt_steps > 0
        gp_X = X[training_indices] if run_gp_optimization else X
        gp_y = gp_Z_scaled[training_indices] if run_gp_optimization else gp_Z_scaled
        gp_err = gp_Z_err_scaled[training_indices] if run_gp_optimization else gp_Z_err_scaled
        gp_holdout = None
        gp_settings = gp.exact_gp_defaults()
        gp_settings.update(
            {
                "fit_steps": int(args.gp_opt_steps),
                "fit_restarts": int(args.gp_restarts),
                "learning_rate": float(args.gp_lr),
                "initial_lengthscale": float(args.gp_scale),
                "initial_noise": float(args.gp_noise),
            }
        )
        gp_geometry = ProductManifoldGeometry.from_metadata(
            param_names=param_names,
            numeric_bounds=tuple(tuple(row) for row in np.asarray(bounds, dtype=np.float64)),
            topology=parameter_topology,
        )
        torch_model_gp = gp.ExactGaussianProcess(
            X_train=gp_tensor(gp_X),
            Y_train=gp_tensor(gp_y),
            Y_err=gp_tensor(gp_err) if args.mc_errors else None,
            geometry=gp_geometry,
            settings=gp_settings,
            device=args.device,
            dtype=gp.resolve_gp_dtype(args.device, args.gp_dtype),
            seed=args.rngseed,
        )
        if run_gp_optimization:
            gp_stats = torch_model_gp.optimize_hyperparameters(
                num_steps=args.gp_opt_steps,
                lr=args.gp_lr,
                X_val=gp_tensor(X[validation_indices]),
                Y_val=gp_tensor(gp_Z_scaled[validation_indices]),
                Y_val_err=gp_tensor(gp_Z_err_scaled[validation_indices]) if args.mc_errors else None,
                patience=args.gp_patience,
                min_delta=args.gp_min_delta,
                validation_interval=args.gp_validation_interval,
                log_interval=args.gp_log_interval,
                restarts=args.gp_restarts,
            )
            fit.print_gp_validation_summary(gp_stats)
            gp_holdout = gp_holdout_metrics(
                model=torch_model_gp,
                X_validation=X[validation_indices],
                Z_validation=Z_fit[validation_indices],
                transform=gp_target_transform,
            )
            fit.print_table(
                "GP objective-space holdout summary:",
                ["Quantity", "Value"],
                [
                    ["validation points", gp_holdout["count"]],
                    ["MAE in Z", gp_holdout["mae"]],
                    [
                        "median absolute relative error",
                        gp_holdout["median_absolute_relative_error"],
                    ],
                    ["lowest 20 percent Z threshold", gp_holdout["low_20_percent_threshold"]],
                    ["lowest 20 percent points", gp_holdout["low_20_percent_count"]],
                    ["lowest 20 percent MAE in Z", gp_holdout["low_20_percent_mae"]],
                    ["minimum predicted validation Z", gp_holdout["minimum_prediction"]],
                ],
            )
            cprint(
                f"{gp.training_timestamp()} Refitting GP posterior with all "
                f"{len(X)} retained trials",
                "yellow",
            )
            torch_model_gp.set_training_data(
                X_train=gp_tensor(X),
                Y_train=gp_tensor(gp_Z_scaled),
                Y_err=gp_tensor(gp_Z_err_scaled) if args.mc_errors else None,
            )

        torch_model_gp.freeze_for_inference()

        # Map the GP latent mean back into physical objective units
        def gp_predict_model(X_new_model):
            mean_scaled = torch_model_gp.predict_mean(X_new_model)
            mean_pred, _ = inverse_gp_prediction(
                mean_scaled=mean_scaled,
                variance_scaled=None,
                transform=gp_target_transform,
            )
            return mean_pred.reshape(-1)

        (
            func_hat_gp,
            gp_value_gradient,
            gp_hessian,
            gp_objective,
            gp_objective_batch,
        ) = physical_surrogate_callbacks(
            predict_model_torch=gp_predict_model,
            physical_features=lambda values: np.asarray(values, dtype=np.float64),
            physical_features_torch=lambda values: values,
            device=torch_model_gp.device,
            dtype=torch_model_gp.dtype,
        )

        # Compute GP predictive mean and standard deviation at one point
        def gp_z_uncertainty(X_new):
            X_new_model = np.asarray(X_new, dtype=np.float64).reshape(1, -1)
            X_new_tensor = torch.as_tensor(
                X_new_model,
                dtype=torch_model_gp.dtype,
                device=torch_model_gp.device,
            )
            with torch.inference_mode():
                mean_scaled, variance_scaled = torch_model_gp(X_new_tensor)
                mean_pred, var_pred = inverse_gp_prediction(
                    mean_scaled, variance_scaled, gp_target_transform
                )
            sigma_pred = torch.sqrt(torch.clamp(var_pred, min=0.0))
            return float(mean_pred.item()), float(sigma_pred.item())

        fit_likelihood_surrogate(
            surface=surface,
            model_name="gp",
            func_hat=func_hat_gp,
            realized_best=realized_best,
            realized_Z=realized_Z,
            value_gradient_func=gp_value_gradient,
            hessian_func=gp_hessian,
            objective_torch=gp_objective,
            objective_batch_torch=gp_objective_batch,
            torch_device=torch_model_gp.device,
            torch_dtype=torch_model_gp.dtype,
            z_uncertainty_func=gp_z_uncertainty,
            args=args,
        )
        update_json_file(
            pathlib.Path(args.output_root) / "gp" / "parameters.json",
            {
                "gp_target_transform": gp_target_transform,
                "gp_holdout": gp_holdout,
            },
        )
    # Fit Neural Network surrogate

    if args.model == "neural":
        model_param = neural.lzmlp_parameters(
            in_dim=in_dim,
            out_dim=out_dim,
            hidden_dim=in_dim,
            hidden_layers=4,
        )
        _, torch_model_nn, stats = neural.train_validated_lzmlp(
            model_param=model_param,
            training=(X_train, y_train, err_train),
            validation=(X_val, y_val, err_val),
            full=(X_model, Z_scaled, Z_err_scaled),
            options=surrogate_cli.neural_training_options(args),
            use_errors=args.mc_errors,
            print_model=True,
        )
        fit.print_neural_validation_summary(stats)

        full_refit_epochs = int(stats.get("best_epoch") or 0)
        cprint(
            f"Neural full refit used {full_refit_epochs} selected epochs "
            f"and all {len(X_model)} retained trials",
            "yellow",
        )
        nn_parameter = next(torch_model_nn.parameters())

        # Map the standardized neural prediction back into physical objective units
        def nn_predict_model(X_new_model):
            prediction = torch_model_nn(X_new_model)
            return (prediction * Z_std_scalar + Z_mu_scalar).reshape(-1)

        (
            func_hat_nn,
            nn_value_gradient,
            nn_hessian,
            nn_objective,
            nn_objective_batch,
        ) = physical_surrogate_callbacks(
            predict_model_torch=nn_predict_model,
            physical_features=physical_features,
            physical_features_torch=physical_features_torch,
            device=nn_parameter.device,
            dtype=nn_parameter.dtype,
        )

        fit_likelihood_surrogate(
            surface=surface,
            model_name="neural",
            func_hat=func_hat_nn,
            realized_best=realized_best,
            realized_Z=realized_Z,
            value_gradient_func=nn_value_gradient,
            hessian_func=nn_hessian,
            objective_torch=nn_objective,
            objective_batch_torch=nn_objective_batch,
            torch_device=nn_parameter.device,
            torch_dtype=nn_parameter.dtype,
            args=args,
        )

        # Save the neural network surrogate model for later inference
        surrogate_save_path = pathlib.Path(args.output_root) / "neural" / "neural_likelihood.pt"
        ensure_dir(surrogate_save_path.parent)

        prior_boundaries = bounds.T.astype(np.float32)

        # Build a checkpoint dictionary. You can also add any extra info you need
        checkpoint = {
            "model_param": model_param,
            "model_state_dict": torch_model_nn.state_dict(),
            "prior_boundaries": torch.tensor(prior_boundaries, dtype=torch.float32),
            "config_keys": param_names,
            "parameter_space": surface["parameter_space"],
            "parameter_topology": parameter_topology,
            "model_feature_names": surface["model_feature_names"],
            "optimization": surface["optimization"],
            "input_mode": surface["input_mode"],
            "history_path": surface["history_path"],
            "Z_mu": Z_mu,
            "Z_std": Z_std,
            "training_stats": json_safe(stats),
        }

        torch.save(checkpoint, surrogate_save_path)
        fit.print_table(
            "Saved outputs:",
            ["Output", "Path"],
            [["neural likelihood checkpoint", surrogate_save_path]],
        )
    surrogate_cli.complete_run(args, tool="icescape", started=start_time)


if __name__ == "__main__":
    main()

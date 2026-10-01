# Shared command-line interface for likelihood and histogram surrogate tools
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import argparse
import os
import pathlib
import time

import numpy as np
import torch
from termcolor import cprint

from core.inference import fit, mcmc, profile
from core.io import cli
from core.io.files import ensure_dir
from core.io.serialize import json_safe, write_json_file
from core.plot.style import resolve_plot_brand
from core.stats.transform import ParameterTransform
from core.tune import manifest as run_manifest
from core.tune.drivers.registry import create_driver
from core.tune.runtime import process as iceruntime


# Add input, quality-cut and runtime arguments
def _add_input_arguments(parser: argparse.ArgumentParser, replica: bool) -> None:
    cli.add_value(parser, "cdir", os.getcwd())
    cli.add_value(parser, "rngseed", 1234, help="Random seed")
    cli.add_value(parser, "device", "auto", help="Torch device <auto|cpu|cuda>")
    parser.add_argument(
        "--input",
        required=True,
        help=("Input icetune run directory" if replica else "Input history JSON or icetune run directory"),
    )
    cli.add_value(parser, "output_name", None, help="Optional output label")
    cli.add_value(parser, "max_trials", 1000000, help="Maximum trials to read")
    cli.add_value(parser, "max_delta_Z", 300.0, help="Maximum delta Z")
    cli.add_value(parser, "quality_percentile", 90.0, help="Quality percentile")
    cli.add_value(parser, "max_Z", float("inf"), help="Maximum absolute Z")
    cli.add_value(
        parser,
        "quality_mode",
        "percentile",
        choices=["delta", "percentile", "absolute", "none"],
        help="Quality cut mode",
    )
    cli.add_flag(parser, "cached", help="Use cached input")
    cli.add_flag(parser, "mc-errors", help="Use MC statistical uncertainties during surrogate training")
    cli.add_value(parser, "validation_fraction", 0.2, help="Validation split fraction")


# Add Gaussian-process arguments
def _add_gp_arguments(parser: argparse.ArgumentParser, replica: bool) -> None:
    cli.add_value(parser, "gp_kernel", "Matern", choices=["Matern"], help="GP kernel")
    cli.add_value(parser, "gp_scale", 0.25, help="GP length scale")
    cli.add_value(parser, "gp_noise", 1e-4, help="GP diagonal noise")
    cli.add_value(parser, "gp_dtype", "float64", choices=["auto", "float32", "float64"], help="GP precision")
    cli.add_toggle(parser, "gp-optimize", help="Optimize GP hyperparameters")
    cli.add_value(parser, "gp_opt_steps", 150, help="GP optimization steps")
    if not replica:
        cli.add_value(parser, "gp_restarts", 1, help="GP optimization restarts")
    cli.add_value(parser, "gp_lr", 0.02, help="GP learning rate")
    cli.add_value(parser, "gp_patience", 40, help="GP validation patience")
    cli.add_value(parser, "gp_min_delta", 1e-4, help="GP minimum improvement")
    cli.add_value(parser, "gp_validation_interval", 5, help="GP validation interval")
    cli.add_value(parser, "gp_log_interval", 20, help="GP logging interval")
    if not replica:
        cli.add_value(
            parser, "gp_target_transform", "auto", choices=["auto", "standard", "log"], help="GP objective transform"
        )


# Add neural-surrogate arguments
def _add_neural_arguments(parser: argparse.ArgumentParser, replica: bool) -> None:
    cli.add_value(
        parser,
        "nn_loss",
        "MSE" if replica else "MAE",
        choices=["MAE", "MSE"] if replica else None,
        help="Neural loss function",
    )
    cli.add_value(parser, "nn_num_epochs", 3000 if replica else 5000, help="Training epochs")
    cli.add_value(parser, "nn_batch_size", 128, help="Training batch size")
    cli.add_value(parser, "nn_lr", 1e-2, help="Learning rate")
    cli.add_value(parser, "nn_weight_decay", 1e-3, help="Weight decay")
    cli.add_value(parser, "nn_lipschitz", 1e-4, help="Lipschitz penalty")
    cli.add_value(parser, "nn_gamma", 1e-4, help="Scheduler gamma")
    cli.add_value(parser, "nn_patience", 250, help="Validation patience")
    cli.add_value(parser, "nn_min_delta", 1e-4, help="Minimum improvement")
    if replica:
        cli.add_value(parser, "nn_hidden_dim", 0, help="Hidden layer width")


# Add minimization, profile and posterior arguments
def _add_fit_arguments(parser: argparse.ArgumentParser, replica: bool) -> None:
    if replica:
        cli.add_value(parser, "random_samples", 100000, help="Random scan samples")
        cli.add_value(parser, "random_batch_size", 512, help="Random batch size")
        cli.add_toggle(parser, "optimize", help="Run surrogate refinement")
    cli.add_value(
        parser, "optimizer_backend", "torch", choices=["torch"], help="Bounded surrogate optimizer"
    )
    cli.add_value(parser, "optimizer_lbfgs_maxiter", 200, help="Maximum L-BFGS-B iterations")
    cli.add_value(parser, "optimizer_restarts", 32, help="Batched Torch optimizer starts")
    cli.add_value(parser, "optimizer_scan_points", 1024, help="Sobol global scan points")
    parser.add_argument(
        "--profile",
        type=cli.parse_bool,
        default=True,
        metavar="<true|false>",
        help="Render profile scans, default true",
    )
    cli.add_flag(parser, "plot-2d", help="Render two-dimensional profiles")
    cli.add_value(parser, "profile_maxiter", 40, help="Profile iterations")
    cli.add_value(parser, "profile_workers", 0, help="Parallel profile workers, 0 uses all available CPUs")
    cli.add_flag(parser, "profile-multistart", help="Use multiple profile starts")
    cli.add_flag(parser, "posterior", help="Sample the posterior with Torch HMC")
    cli.add_value(parser, "posterior_samples", 20000, help="Posterior samples")
    cli.add_value(parser, "posterior_burnin", 2000, help="Posterior burn-in")
    cli.add_value(parser, "posterior_thin", 1, help="Posterior thinning")
    cli.add_value(parser, "posterior_step", 0.05, help="Initial HMC step in unit coordinates")
    cli.add_value(parser, "posterior_leapfrog", 10, help="Leapfrog steps per HMC trajectory")
    cli.add_value(parser, "posterior_target_accept", 0.8, help="HMC warmup acceptance target")
    cli.add_value(parser, "posterior_jitter", 0.1, help="Relative HMC step jitter")


# Add plotting and impact arguments
def _add_plot_arguments(parser: argparse.ArgumentParser, replica: bool) -> None:
    cli.add_value(parser, "vmax", None, value_type=float, help="Visualization maximum")
    cli.add_value(parser, "ngrid", 20, help="Grid points per dimension")
    cli.add_value(parser, "profile_dpi", 300, help="Profile PNG resolution in DPI")
    cli.add_value(parser, "ncontour", 8, help="Contour levels")
    cli.add_value(parser, "cmap", "viridis_r", help="Colormap")
    cli.add_value(parser, "plot-brand", None, help="Plot producer name")
    cli.add_value(parser, "impact-mode", "profiled", choices=["fixed", "profiled"], help="Impact shift mode")
    cli.add_value(parser, "impact-profile-maxiter", 30, help="Profiled impact iterations")
    cli.add_value(parser, "impact-top", 20, help="Maximum impact parameters")
    if replica:
        cli.add_toggle(parser, "predictive-plots", help="Render predictive tune figures")
    else:
        cli.add_value(parser, "impact-features", 256, help="RFF impact features")
        cli.add_value(parser, "impact-scale", None, value_type=float, help="RFF impact scale")
        cli.add_value(parser, "impact-ridge", 1e-3, help="Impact ridge penalty")
        cli.add_value(parser, "impact-min-trials", 30, help="Minimum impact trials")


# Add histogram-replica model arguments
def _add_replica_arguments(parser: argparse.ArgumentParser) -> None:
    cli.add_value(parser, "pca_variance", 0.995, help="Retained PCA variance fraction")
    cli.add_value(parser, "pca_max_modes", 24, help="Maximum PCA modes")
    cli.add_value(parser, "rff_features", 1024, help="RFF features")
    cli.add_value(parser, "rff_scale", 0.25, help="RFF length scale")
    cli.add_value(parser, "rff_ensemble", 1, help="RFF ensemble members")
    cli.add_value(parser, "ridge_alpha", 1e-5, help="Ridge penalty")


# Resolve one standard icetune run directory into its base path and run name
def _standard_run_identity(run_directory: pathlib.Path) -> tuple[str, str]:
    resolved = pathlib.Path(run_directory).resolve()
    parts = resolved.parts
    marker = next(
        (index for index in range(len(parts) - 1) if parts[index] == "runs" and parts[index + 1] == "icetune"), None
    )
    if marker is None or marker + 2 >= len(parts):
        raise ValueError("Full histogram input must be under <base>/runs/icetune/<run name>")
    cdir = pathlib.Path(*parts[:marker])
    run_name = pathlib.Path(*parts[marker + 2 :])
    return str(cdir), str(run_name)


# Infer the surrogate input reader and internal run identity from one path
def normalize_input_args(args: argparse.Namespace, replica: bool) -> argparse.Namespace:
    input_path = pathlib.Path(args.input).expanduser().resolve()
    history_path = None
    run_directory = None
    input_mode = None

    if input_path.is_file() or input_path.suffix.lower() in {".json", ".pkl"}:
        suffix = input_path.suffix.lower()
        if suffix == ".json":
            input_mode = "history"
            history_path = input_path
            run_directory = input_path.parent
        elif suffix == ".pkl" and input_path.parent.name == "results":
            input_mode = "full"
            run_directory = input_path.parent.parent
        else:
            raise ValueError("--input must be a history JSON or an icetune run input")
    elif input_path.is_dir():
        if input_path.name == "results":
            history_candidates = []
            full_directory = input_path
            candidate_run = input_path.parent
        else:
            history_candidates = [input_path / "history.json"]
            full_directory = input_path / "results"
            candidate_run = input_path
        available_history = next((path for path in history_candidates if path.is_file()), None)
        if not replica and available_history is not None:
            input_mode = "history"
            history_path = available_history
            run_directory = candidate_run
        elif full_directory is not None and full_directory.is_dir():
            input_mode = "full"
            run_directory = candidate_run
        else:
            raise FileNotFoundError(f"Could not infer input type from {input_path}")
    else:
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

    if replica and input_mode != "full":
        raise ValueError("iceproxy requires an icetune run with full trial pickle payloads")

    try:
        cdir, run_name = _standard_run_identity(run_directory)
    except ValueError:
        if input_mode == "full":
            raise
        cdir = str(pathlib.Path(args.cdir).expanduser().resolve())
        run_name = run_directory.name

    args.input = str(input_path)
    args.input_mode = input_mode
    args.history_path = None if history_path is None else str(history_path)
    args.cdir = cdir
    args.run_name = run_name
    return args


# Normalize shared runtime selections and validate their ranges
def _normalize_args(args: argparse.Namespace, replica: bool) -> argparse.Namespace:
    args = normalize_input_args(args, replica)
    if args.device == "auto":
        args.device = "cuda" if torch.cuda.is_available() else "cpu"

    rules = [
        (args.posterior_thin >= 1, "--posterior_thin must be >= 1"),
        (args.posterior_samples >= 1, "--posterior_samples must be >= 1"),
        (args.posterior_burnin >= 0, "--posterior_burnin must be >= 0"),
        (args.posterior_leapfrog >= 1, "--posterior_leapfrog must be >= 1"),
        (np.isfinite(args.posterior_step) and args.posterior_step > 0, "--posterior_step must be finite and positive"),
        (0 < args.posterior_target_accept < 1, "--posterior_target_accept must be in (0, 1)"),
        (0 <= args.posterior_jitter < 1, "--posterior_jitter must be in [0, 1)"),
        (0.0 < args.validation_fraction < 1.0, "--validation_fraction must be in (0, 1)"),
        (args.gp_patience >= 0, "--gp_patience must be >= 0"),
        (args.gp_opt_steps >= 0, "--gp_opt_steps must be >= 0"),
        (args.gp_min_delta >= 0, "--gp_min_delta must be >= 0"),
        (args.gp_validation_interval >= 1, "--gp_validation_interval must be >= 1"),
        (args.gp_log_interval >= 1, "--gp_log_interval must be >= 1"),
        (args.profile_dpi >= 72, "--profile_dpi must be at least 72"),
        (args.profile_workers >= 0, "--profile_workers must be non-negative"),
        (args.nn_patience >= 0, "--nn_patience must be >= 0"),
        (args.nn_min_delta >= 0, "--nn_min_delta must be >= 0"),
        (args.impact_profile_maxiter >= 1, "--impact-profile-maxiter must be >= 1"),
        (args.impact_top >= 1, "--impact-top must be >= 1"),
        (args.optimizer_lbfgs_maxiter >= 1, "L-BFGS-B iterations must be positive"),
        (args.optimizer_restarts >= 1, "Torch optimizer restarts must be positive"),
        (args.optimizer_scan_points >= 0, "Torch optimizer scan points must be non-negative"),
    ]
    if replica:
        rules.extend(
            [
                (args.nn_num_epochs >= 1, "--nn_num_epochs must be >= 1"),
                (args.nn_batch_size >= 1, "--nn_batch_size must be >= 1"),
                (args.random_samples >= 1, "--random_samples must be >= 1"),
                (args.random_batch_size >= 1, "--random_batch_size must be >= 1"),
                (0.0 < args.pca_variance <= 1.0, "--pca_variance must be in (0, 1]"),
                (args.pca_max_modes >= 1, "--pca_max_modes must be >= 1"),
                (args.rff_features >= 1, "--rff_features must be >= 1"),
                (args.rff_scale > 0.0, "--rff_scale must be > 0"),
                (args.rff_ensemble >= 1, "--rff_ensemble must be >= 1"),
                (args.ridge_alpha >= 0.0, "--ridge_alpha must be >= 0"),
            ]
        )
    else:
        rules.extend(
            [
                (args.gp_lr > 0.0, "--gp_lr must be > 0"),
                (args.gp_restarts >= 1, "--gp_restarts must be >= 1"),
                (args.impact_features >= 1, "--impact-features must be >= 1"),
                (args.impact_scale is None or args.impact_scale > 0.0, "--impact-scale must be > 0"),
                (args.impact_ridge >= 0.0, "--impact-ridge must be >= 0"),
                (args.impact_min_trials >= 3, "--impact-min-trials must be >= 3"),
            ]
        )
    cli.validate(rules)
    return args


# Parse the shared surrogate interface for one workflow
def parse_args(*, description: str, replica: bool) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    _add_input_arguments(parser, replica)
    cli.add_value(
        parser, "model", "gp", choices=["rff", "gp", "neural"] if replica else ["gp", "neural"], help="Surrogate model"
    )
    if replica:
        _add_replica_arguments(parser)
    _add_gp_arguments(parser, replica)
    _add_neural_arguments(parser, replica)
    _add_fit_arguments(parser, replica)
    _add_plot_arguments(parser, replica)
    return _normalize_args(parser.parse_args(), replica)


# Project the shared neural command-line options onto the trainer interface
def neural_training_options(args: argparse.Namespace) -> dict:
    names = ("lr", "weight_decay", "gamma", "lipschitz", "num_epochs", "batch_size", "patience", "min_delta")
    options = {name: getattr(args, f"nn_{name}") for name in names}
    return options | {"loss_name": args.nn_loss, "device": args.device, "dtype": torch.float64}


# Build the shared Torch configuration for profile and impact minimization
def torch_minimizer_config(
    *, args: argparse.Namespace, objective_torch, objective_batch_torch, device: str | torch.device, dtype: torch.dtype
) -> dict | None:
    if args.optimizer_backend != "torch":
        return None
    return {
        "objective": objective_torch,
        "objective_batch": objective_batch_torch,
        "device": device,
        "dtype": dtype,
        "gradient_tolerance": 1.0e-5,
        "relative_gradient_tolerance": 1.0e-6,
    }


# Dispatch one surrogate fit through the shared Torch backend
def optimize_surrogate(
    *,
    args: argparse.Namespace,
    X: np.ndarray,
    Z: np.ndarray,
    param_names,
    func_hat,
    bounds: np.ndarray,
    output_root: pathlib.Path,
    realized_best: np.ndarray,
    realized_Z: float,
    objective_torch,
    objective_batch_torch,
    device: str | torch.device,
    dtype: torch.dtype,
    x0: np.ndarray | None = None,
    transform=None,
):
    common = {
        "X": X,
        "Z": Z,
        "param_names": param_names,
        "func_hat": func_hat,
        "bounds": bounds,
        "output_root": output_root,
        "realized_best": realized_best,
        "realized_Z": realized_Z,
        "x0": x0,
        "plot_brand": args.plot_brand,
    }
    if args.optimizer_backend == "torch":
        result = fit.optimize_with_torch(
            **common,
            objective_torch=objective_torch,
            objective_batch_torch=objective_batch_torch,
            device=device,
            dtype=dtype,
            maxiter=int(args.optimizer_lbfgs_maxiter),
            restarts=int(args.optimizer_restarts),
            scan_points=int(args.optimizer_scan_points),
            rngseed=int(args.rngseed),
        )
    else:
        raise ValueError(f"Unknown optimizer backend = {args.optimizer_backend}")
    if transform is not None:
        covariance, _ = fit.fitted_covariance(result[0], param_names)
        fit.save_physical_parameters(transform=transform, best_fit=result[1], covariance=covariance,
                                     output_root=output_root, realized_best=realized_best, plot_brand=args.plot_brand)
    return result


# Bind the simulator decoder once for every physical parameter diagnostic
def parameter_transform(args, surface):
    if "parameter_transform" not in surface:
        names, reference = surface["param_names"], surface["X"][0]
        driver = surface.get("simdriver")
        surface["parameter_transform"] = (
            create_driver(driver).parameter_transform(names, reference, cdir=args.cdir, metadata=surface)
            if driver else ParameterTransform(names, dict, reference)
        )
    return surface["parameter_transform"]


# Compute the covariance method associated with the selected shared backend
def optimizer_uncertainty_method(args: argparse.Namespace) -> str:
    methods = {"torch": "torch_autograd_hessian"}
    try:
        return methods[args.optimizer_backend]
    except KeyError as exc:
        raise ValueError(f"Unknown optimizer backend = {args.optimizer_backend}") from exc


# Render one shared profile-likelihood configuration
def render_profiles(
    *,
    args: argparse.Namespace,
    X: np.ndarray,
    Z: np.ndarray,
    param_names,
    func_hat,
    output_root: pathlib.Path,
    realized_best: np.ndarray,
    surrogate_best: np.ndarray,
    realized_Z: float,
    surrogate_Z: float,
    bounds: np.ndarray,
    value_gradient_func=None,
    torch_config: dict | None = None,
    transform=None,
) -> bool:
    if not args.profile:
        cprint("Profile-likelihood scans disabled by --profile false", "yellow")
        return False
    profile.plot_profile_likelihoods(
        X=X,
        Z=Z,
        param_names=param_names,
        func_hat=func_hat,
        output_root=pathlib.Path(output_root) / "optimizer",
        realized_best=realized_best,
        surrogate_best=surrogate_best,
        realized_Z=realized_Z,
        surrogate_Z=surrogate_Z,
        bounds=bounds,
        ngrid=int(args.ngrid),
        colorbar_title="$Z = -2 \\log L(\\theta) = 2\\,\\mathrm{NLL}$",
        cmap=args.cmap,
        vmax=args.vmax,
        plot_2d=bool(args.plot_2d),
        maxiter=int(args.profile_maxiter),
        multistart=bool(args.profile_multistart),
        value_gradient_func=value_gradient_func,
        torch_config=torch_config,
        workers=int(args.profile_workers),
        dpi=int(args.profile_dpi),
        plot_brand=args.plot_brand,
    )
    if transform is not None:
        profile.plot_transformed_profiles(transform=transform, X=X, best_fit=surrogate_best, bounds=bounds,
            value_gradient_func=value_gradient_func, output_root=pathlib.Path(output_root) / "physical",
            ngrid=int(args.ngrid), maxiter=int(args.profile_maxiter), workers=int(args.profile_workers),
            plot_2d=bool(args.plot_2d), plot_brand=args.plot_brand, dpi=int(args.profile_dpi), cmap=args.cmap)
    return True


# Run the shared optional posterior sampler from one fitted surrogate point
def run_posterior(
    *, args: argparse.Namespace, objective_torch, device, dtype, bounds: np.ndarray, param_names,
    output_root: pathlib.Path, start: np.ndarray,
    transform=None,
) -> None:
    result = mcmc.run_configured_sampling(
        args=args,
        objective_torch=objective_torch,
        device=device,
        dtype=dtype,
        bounds=bounds,
        param_names=param_names,
        output_dir=pathlib.Path(output_root) / "optimizer" / "posterior",
        start=start,
    )
    if result is not None and transform is not None:
        output = pathlib.Path(output_root) / "physical" / "posterior"
        ensure_dir(output)
        samples = transform.samples(result[0])
        reference = transform.samples([start])[0]
        unwrapped = reference + transform.difference(torch.as_tensor(samples), torch.as_tensor(reference)).numpy()
        np.savez(output / "posterior_samples.npz", samples=samples, samples_unwrapped=unwrapped,
                 param_names=np.asarray(transform.names), periods=transform.periods)
        write_json_file(output / "posterior_summary.json", json_safe({
            **mcmc.summarize_posterior_samples(unwrapped, transform.names), "reference": reference,
            "periods": transform.periods, "method": "decoded posterior samples",
        }), indent=4)
        mcmc.plot_posterior_samples(unwrapped, transform.names, output)


# Remove invalid and poor-quality rows from one aligned surrogate surface
def quality_filter_surface(
    surface: dict, *, finite_keys: tuple[str, ...], row_keys: tuple[str, ...], args, minimum_rows: int = 3
) -> tuple[dict, float]:
    rows_before = len(surface["Z"])
    finite = np.ones(rows_before, dtype=bool)
    for key in finite_keys:
        values = np.asarray(surface[key])
        finite &= np.all(np.isfinite(values.reshape(rows_before, -1)), axis=1)
    for key in row_keys:
        surface[key] = np.asarray(surface[key])[finite]
    rows_finite = len(surface["Z"])

    keep, threshold = fit.quality_mask(
        surface["Z"],
        mode=args.quality_mode,
        max_delta_Z=args.max_delta_Z,
        quality_percentile=args.quality_percentile,
        max_Z=args.max_Z,
    )
    for key in row_keys:
        surface[key] = surface[key][keep]
    if len(surface["Z"]) < minimum_rows:
        raise RuntimeError(f"Need at least {minimum_rows} trials after quality cut, got {len(surface['Z'])}")
    fit.print_table(
        "Surface quality cut:",
        ["Quantity", "Value"],
        [
            ["input rows", rows_before],
            ["finite rows", rows_finite],
            ["quality mode", args.quality_mode],
            ["Z threshold", threshold],
            ["retained rows", len(surface["Z"])],
            ["Z range", f"{np.min(surface['Z']):0.3f} to {np.max(surface['Z']):0.3f}"],
        ],
        color="green",
    )
    return surface, float(threshold)


# Initialize one surrogate run with reproducible outputs and runtime state
def start_run(args, *, tool: str, version, command: list[str]) -> float:
    started = time.perf_counter()
    args.output_root = run_manifest.create_run_output(
        cdir=args.cdir,
        tool=tool,
        tool_version=version,
        input_run_name=args.run_name,
        output_name=args.output_name,
        arguments=vars(args),
        command=command,
    )
    fit.configure_torch_runtime(args.device)
    fit.print_table(f"{tool} arguments:", ["Argument", "Value"], sorted(vars(args).items()))
    iceruntime.set_random_seeds(args.rngseed, use_torch=True)
    return started


# Record the shared surface identity, dimensions and plotting presentation
def record_surface(args, surface: dict, *, inputs: dict, **dimensions: int) -> None:
    args.plot_brand = resolve_plot_brand(
        override=args.plot_brand,
        saved_brand=surface.get("plot_brand"),
        simdriver=surface.get("simdriver"),
        param_names=surface["param_names"],
    )
    run_manifest.update_run_manifest(
        args.output_root,
        {
            "inputs": inputs,
            "presentation": {"plot_brand": args.plot_brand},
            "surface": {
                "trial_count": int(len(surface["Z"])),
                "parameter_count": int(len(surface["param_names"])),
                **dimensions,
            },
        },
    )


# Complete one surrogate run and publish its elapsed runtime
def complete_run(args, *, tool: str, started: float, model: str | None = None) -> None:
    elapsed = time.perf_counter() - started
    updates = {"status": "completed", "elapsed_seconds": elapsed}
    if model is not None:
        updates["surrogate_model"] = model
    run_manifest.update_run_manifest(args.output_root, updates)
    cprint(iceruntime.format_done_message(tool, elapsed), "yellow")

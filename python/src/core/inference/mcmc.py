# Markov Chain Monte Carlo algorithms and Bayesian analysis
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import glob
import math
import os
import pathlib
import pickle
from collections.abc import Callable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import torch
from scipy.signal import correlate
from termcolor import cprint
from tqdm import tqdm

from core.inference.profile import safe_name
from core.io.files import ensure_dir
from core.io.serialize import json_safe, write_json_file
from core.tune.parameters.space import from_unit, to_unit


# [REFERENCE: Neal, MCMC Using Hamiltonian Dynamics, 2011, Section 5.5.1.5]
# Reflect one free Hamiltonian position update exactly through finite box boundaries
def reflective_position_update(pos, mom, dt, lower_bounds, upper_bounds):
    width = upper_bounds - lower_bounds
    if torch.any(width <= 0.0):
        raise ValueError("Reflective boundaries require lower < upper in every dimension")
    unfolded = pos + dt * mom - lower_bounds
    phase = torch.remainder(unfolded, 2.0 * width)
    ascending = phase <= width
    reflected = torch.where(ascending, lower_bounds + phase, upper_bounds - (phase - width))
    quotient = unfolded / width
    nearest_segment = torch.round(quotient)
    roundoff = 8.0 * torch.finfo(quotient.dtype).eps * torch.maximum(torch.abs(quotient), torch.ones_like(quotient))
    final_boundary = torch.abs(quotient - nearest_segment) <= roundoff
    boundary_segment = nearest_segment.to(torch.int64)
    segment = torch.where(final_boundary, boundary_segment, torch.floor(quotient).to(torch.int64))
    boundary_position = torch.where(torch.remainder(boundary_segment, 2) == 0, lower_bounds, upper_bounds)
    reflected = torch.where(final_boundary, boundary_position, reflected)
    segment -= (final_boundary & (dt * mom > 0.0)).to(segment.dtype)
    same_direction = torch.remainder(segment, 2) == 0
    reflected_momentum = torch.where(same_direction, mom, -mom)
    return torch.clamp(reflected, min=lower_bounds, max=upper_bounds), reflected_momentum


# Evaluate and detach the log density and its autograd gradient at one position
def _logp_grad(position, log_posterior_fn):
    if not torch.isfinite(position).all():
        return position.new_tensor(math.nan), None
    with torch.enable_grad():
        position = position.detach().requires_grad_(True)
        logp, valid = log_posterior_fn(position)
        if logp.ndim != 0:
            raise ValueError("HMC log posterior must be a scalar Torch tensor")
        if not valid or not torch.isfinite(logp):
            return logp.detach(), None
        grad = torch.autograd.grad(logp, position)[0]
    return logp.detach(), grad.detach() if torch.isfinite(grad).all() else None


# Integrate reversible kick-drift-kick steps with unit mass and optional box reflections
def leapfrog(position, momentum, grad, step, n_steps, log_posterior_fn, lower=None, upper=None):
    for _ in range(n_steps):
        momentum = momentum + 0.5 * step * grad
        if lower is None:
            position = position + step * momentum
        else:
            position, momentum = reflective_position_update(position, momentum, step, lower, upper)
        logp, grad = _logp_grad(position, log_posterior_fn)
        if grad is None:
            return position, momentum, logp, None
        momentum = momentum + 0.5 * step * grad
    return position, momentum, logp, grad


# Draw a Hamiltonian chain including the initial adaptation_window warmup transitions
def hmc_sampler(
    initial_position: torch.Tensor,
    n_samples: int,
    step_size: float,
    n_leapfrog_steps: int,
    log_posterior_fn: Callable[[torch.Tensor], tuple[torch.Tensor, bool]],
    adapt_step_size: bool = True,
    target_accept: float = 0.8,
    adaptation_window: int = 100,
    lower_bounds: torch.Tensor | None = None,
    upper_bounds: torch.Tensor | None = None,
    device: str = "cpu",
    generator: torch.Generator | None = None,
    step_jitter: float = 0.0,
    max_energy_error: float = 1000.0,
    step_size_bounds: tuple[float, float] = (1e-10, 1e10),
) -> tuple[torch.Tensor, dict]:
    if n_samples < 1 or n_leapfrog_steps < 1:
        raise ValueError("HMC sample and leapfrog counts must be positive")
    if not np.isfinite(step_size) or step_size <= 0.0:
        raise ValueError("HMC step size must be finite and positive")
    if not 0.0 <= step_jitter < 1.0 or not math.isfinite(max_energy_error) or max_energy_error <= 0.0:
        raise ValueError("HMC requires jitter in [0, 1) and a finite positive energy error threshold")
    if not 0.0 < step_size_bounds[0] <= step_size <= step_size_bounds[1] < math.inf:
        raise ValueError("HMC step size must lie within finite positive adaptation bounds")
    if (lower_bounds is None) != (upper_bounds is None):
        raise ValueError("HMC requires both lower and upper bounds")
    if initial_position.ndim != 1 or initial_position.numel() == 0 or not torch.isfinite(initial_position).all():
        raise ValueError("HMC initial position must be a finite non-empty vector")
    if not initial_position.is_floating_point():
        raise ValueError("HMC initial position must have floating point dtype")
    if adapt_step_size and (not 0.0 < target_accept < 1.0 or adaptation_window < 0):
        raise ValueError("HMC adaptation requires an acceptance target in (0, 1) and a nonnegative window")
    print(f"Running HMC sampler on device: {device}")
    samples = []
    current_position = initial_position.clone().detach().to(device)
    current_step_size = step_size
    bounded = lower_bounds is not None and upper_bounds is not None
    if bounded:
        lower_bounds = lower_bounds.to(current_position)
        upper_bounds = upper_bounds.to(current_position)
        if (
            lower_bounds.shape != current_position.shape
            or upper_bounds.shape != current_position.shape
            or not torch.isfinite(lower_bounds).all()
            or not torch.isfinite(upper_bounds).all()
            or not torch.all(lower_bounds < upper_bounds)
        ):
            raise ValueError("HMC bounds must be finite aligned vectors with lower < upper")
        if torch.any(current_position < lower_bounds) or torch.any(current_position > upper_bounds):
            raise ValueError("HMC initial position must lie inside its bounds")

    current_logp, current_grad = _logp_grad(current_position, log_posterior_fn)
    if current_grad is None:
        raise ValueError("HMC initial position has an invalid log posterior or non-finite posterior gradient")
    diagnostics = {key: [] for key in (
        "accept_history", "accept_prob_history", "energy_history", "energy_error_history",
        "divergent_history", "log_posterior_history", "step_size_history",
    )}
    warmup = min(adaptation_window, n_samples) if adapt_step_size else 0
    log_center, log_average, error_average = math.log(10.0 * step_size), math.log(step_size), 0.0
    random_options = {"dtype": current_position.dtype, "device": current_position.device, "generator": generator}
    for i in tqdm(range(n_samples), desc="HMC sampling"):
        momentum = torch.randn(current_position.shape, **random_options)
        step = current_step_size
        if step_jitter > 0.0:
            step *= 1.0 + step_jitter * (2.0 * torch.rand((), **random_options).item() - 1.0)
        current_H = -current_logp + 0.5 * momentum.square().sum()
        position, momentum, logp, grad = leapfrog(
            current_position, momentum, current_grad, step, n_leapfrog_steps,
            log_posterior_fn, lower_bounds, upper_bounds,
        )
        # Momentum reversal completes the involution without changing kinetic energy
        new_H = -logp + 0.5 * momentum.square().sum() if grad is not None else current_H.new_tensor(math.inf)
        delta = (new_H - current_H).item()
        finite = math.isfinite(delta)
        probability = math.exp(min(0.0, -delta)) if finite else 0.0
        accept = finite and torch.rand((), **random_options).item() < probability
        if accept:
            current_position, current_logp, current_grad = position, logp, grad
        samples.append(current_position.cpu())
        for key, value in zip(diagnostics, (
            float(accept), probability, new_H.item(), delta, not finite or delta > max_energy_error,
            current_logp.item(), step,
        ), strict=True):
            diagnostics[key].append(value)

        # Adapt only during warmup, then freeze the averaged log step
        # [REFERENCE: Hoffman and Gelman, JMLR 15 (2014) 1593-1623, Algorithm 5]
        if i < warmup:
            iteration = i + 1
            weight = 1.0 / (iteration + 10.0)
            error_average += weight * (target_accept - probability - error_average)
            log_step = log_center - math.sqrt(iteration) * error_average / 0.05
            log_step = min(math.log(step_size_bounds[1]), max(math.log(step_size_bounds[0]), log_step))
            log_average += iteration**(-0.75) * (log_step - log_average)
            current_step_size = math.exp(log_average if iteration == warmup else log_step)

    diagnostics["acceptance_rate"] = float(np.mean(diagnostics["accept_history"]))
    diagnostics["final_step_size"] = current_step_size
    diagnostics["warmup"] = warmup
    diagnostics["target_accept"] = target_accept
    print(f"Final acceptance rate: {diagnostics['acceptance_rate']:.4f}")
    return torch.stack(samples), diagnostics


# Compute marginal posterior summaries and correlation diagnostics
def summarize_posterior_samples(samples: np.ndarray, param_names: Sequence[str]) -> dict:
    samples = np.asarray(samples, dtype=np.float64)
    summary = {}
    for index, name in enumerate(param_names):
        values = samples[:, index]
        summary[name] = {
            "mean": float(np.mean(values)),
            "std": float(np.std(values)),
            "median": float(np.percentile(values, 50)),
            "lower68": float(np.percentile(values, 16)),
            "upper68": float(np.percentile(values, 84)),
            "lower95": float(np.percentile(values, 2.5)),
            "upper95": float(np.percentile(values, 97.5)),
        }
    covariance = np.atleast_2d(np.cov(samples, rowvar=False)) if len(samples) > 1 else np.full(
        (len(param_names), len(param_names)), np.nan,
    )
    scale = np.sqrt(np.maximum(0.0, np.diag(covariance)))
    denominator = np.outer(scale, scale)
    correlation = np.divide(covariance, denominator, out=np.full_like(covariance, np.nan), where=denominator > 0.0)
    return {
        "parameters": summary,
        "covariance_matrix": np.atleast_2d(covariance).tolist(),
        "correlation_matrix": np.atleast_2d(correlation).tolist(),
    }


# Write trace and marginal histogram plots for posterior samples
def plot_posterior_samples(samples: np.ndarray, param_names: Sequence[str], output_dir: pathlib.Path) -> None:
    output_dir = pathlib.Path(output_dir)
    ensure_dir(output_dir)
    for index, name in enumerate(param_names):
        fig, axes = plt.subplots(2, 1, figsize=(7, 5))
        axes[0].plot(samples[:, index], lw=0.6)
        axes[0].set_ylabel(name)
        axes[1].hist(samples[:, index], bins=50, density=True, alpha=0.75)
        axes[1].set_xlabel(name)
        axes[1].set_ylabel("density")
        fig.tight_layout()
        fig.savefig(output_dir / f"posterior__{safe_name(name)}.png", bbox_inches="tight")
        plt.close(fig)


# Sample exp(-0.5 Z_surrogate) with uniform priors over parameter bounds
def run_posterior_sampling(
    objective_torch: Callable,
    bounds: np.ndarray,
    param_names: Sequence[str],
    output_dir: pathlib.Path,
    start: np.ndarray,
    rngseed: int = 1234,
    n_samples: int = 20000,
    burnin: int = 2000,
    thin: int = 1,
    step: float = 0.05,
    n_leapfrog_steps: int = 10,
    target_accept: float = 0.8,
    step_jitter: float = 0.1,
    device: str = "cpu",
    dtype: torch.dtype = torch.float64,
) -> tuple[np.ndarray, dict]:
    if thin < 1:
        raise ValueError("thin must be >= 1")
    if n_samples < 1:
        raise ValueError("n_samples must be >= 1")
    if burnin < 0:
        raise ValueError("burnin must be >= 0")
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("Posterior HMC step must be finite and positive")

    bounds = np.asarray(bounds, dtype=np.float64)
    dimension = len(param_names)
    if (
        dimension < 1
        or bounds.shape != (dimension, 2)
        or not np.isfinite(bounds).all()
        or np.any(bounds[:, 0] >= bounds[:, 1])
    ):
        raise ValueError("Posterior bounds must be finite aligned intervals with lower < upper")
    start = np.asarray(start, dtype=np.float64)
    if start.shape != (dimension,) or not np.isfinite(start).all():
        raise ValueError("Posterior start must be a finite vector aligned with the parameters")
    current = torch.as_tensor(to_unit(start[None, :], bounds)[0], device=device, dtype=dtype)
    lower, width = torch.as_tensor(bounds[:, 0], device=device, dtype=dtype), torch.as_tensor(
        bounds[:, 1] - bounds[:, 0], device=device, dtype=dtype,
    )

    # The affine unit transform has a constant Jacobian, with uniform priors in optimizer coordinates
    def log_posterior_unit(unit):
        return -0.5 * objective_torch(lower + width * unit), True

    chain, diagnostics = hmc_sampler(
        initial_position=current, n_samples=burnin + n_samples * thin,
        step_size=step, n_leapfrog_steps=n_leapfrog_steps, log_posterior_fn=log_posterior_unit,
        adapt_step_size=burnin > 0, adaptation_window=burnin, target_accept=target_accept,
        lower_bounds=torch.zeros_like(current), upper_bounds=torch.ones_like(current), device=device,
        generator=torch.Generator(device=device).manual_seed(int(rngseed)), step_jitter=step_jitter,
        step_size_bounds=(1e-10, max(1.0, step)),
    )
    samples_unit = chain[burnin::thin].numpy()
    samples = from_unit(samples_unit, bounds)
    logp_values = np.asarray(diagnostics["log_posterior_history"])[burnin::thin]

    output_dir = pathlib.Path(output_dir)
    ensure_dir(output_dir)
    np.savez(
        output_dir / "posterior_samples.npz",
        samples=samples,
        samples_unit=samples_unit,
        log_posterior=logp_values,
        param_names=np.asarray(param_names),
    )
    summary = {
        "sampler": "reflective_hmc",
        "target": "exp(-0.5 * Z_surrogate)",
        "prior": "uniform within optimizer bounds",
        "metric": "unit mass in unit box coordinates",
        "initial_step_unit": float(step),
        "final_step_unit": diagnostics["final_step_size"],
        "n_leapfrog_steps": n_leapfrog_steps,
        "target_accept": target_accept,
        "step_jitter": step_jitter,
        "rngseed": int(rngseed),
        "device": str(device),
        "dtype": str(dtype),
        "burnin": int(burnin),
        "thin": int(thin),
        "n_samples": int(len(samples)),
        "acceptance_rate": float(np.mean(diagnostics["accept_history"][burnin:])),
        "mean_accept_probability": float(np.mean(diagnostics["accept_prob_history"][burnin:])),
        "divergences": int(np.sum(diagnostics["divergent_history"][burnin:])),
        "warmup_divergences": int(np.sum(diagnostics["divergent_history"][:burnin])),
        **summarize_posterior_samples(samples, param_names),
    }
    write_json_file(output_dir / "posterior_summary.json", json_safe(summary), indent=4)
    np.savez(output_dir / "posterior_diagnostics.npz", **diagnostics)
    cprint(f"HMC after warmup: acceptance {summary['acceptance_rate']:.4f}, divergences {summary['divergences']}", "yellow")
    plot_posterior_samples(samples, param_names, output_dir)
    cprint(f"Saved posterior samples into: {output_dir}", "yellow")
    return samples, summary


# Run posterior sampling when enabled by shared surrogate CLI arguments
def run_configured_sampling(
    *,
    args,
    objective_torch: Callable,
    device: str,
    dtype: torch.dtype,
    bounds: np.ndarray,
    param_names: Sequence[str],
    output_dir: pathlib.Path,
    start: np.ndarray,
):
    if not args.posterior:
        return None
    return run_posterior_sampling(
        objective_torch=objective_torch,
        device=device,
        dtype=dtype,
        bounds=bounds,
        param_names=param_names,
        output_dir=output_dir,
        start=start,
        rngseed=args.rngseed,
        n_samples=args.posterior_samples,
        burnin=args.posterior_burnin,
        thin=args.posterior_thin,
        step=args.posterior_step,
        n_leapfrog_steps=args.posterior_leapfrog,
        target_accept=args.posterior_target_accept,
        step_jitter=args.posterior_jitter,
    )


# Summarize posterior samples with marginal intervals and joint covariance
def analyze_samples(samples: torch.Tensor, config_keys: list = None) -> dict:
    samples_np = samples.numpy()
    n_params = samples_np.shape[1]

    if config_keys is None:
        config_keys = ["$\\theta_{i}$" for _ in range(n_params)]

    quantiles = np.percentile(samples_np, [50.0, 16.0, 84.0, 2.5, 97.5], axis=0)
    medians = quantiles[0]
    ci68 = quantiles[[1, 2]].T
    ci95 = quantiles[[3, 4]].T

    cov_matrix = np.cov(samples_np, rowvar=False)
    corr_matrix = np.corrcoef(samples_np, rowvar=False)

    analysis = {"median": medians, "ci68": ci68, "ci95": ci95, "covariance": cov_matrix, "correlation": corr_matrix}

    print("Parameter Analysis:")
    for i in range(n_params):
        print(
            f"{config_keys[i]:30s}: median = {medians[i]:.3f}, 68% CI = [{ci68[i, 0]:.3f}, {ci68[i, 1]:.3f}], 95% CI = [{ci95[i, 0]:.3f}, {ci95[i, 1]:.3f}]"
        )
    print("Covariance Matrix:\n", cov_matrix)
    print("Correlation Matrix:\n", corr_matrix)

    return analysis


# Compute the finite-sample autocorrelation through one vectorized correlation
def compute_autocorrelation(x, max_lag: int = 50):
    x = np.asarray(x)
    if x.ndim != 1 or x.size == 0:
        raise ValueError("Autocorrelation requires a non-empty vector")
    if max_lag < 0:
        raise ValueError("Maximum autocorrelation lag must be non-negative")
    centered = x - np.mean(x)
    n = len(centered)
    max_lag = min(max_lag, n - 1)
    lags = np.arange(max_lag + 1)
    products = correlate(centered, centered, mode="full", method="auto")[n - 1 : n + max_lag]
    if products[0] <= 0.0:
        return lags, np.full(len(lags), np.nan)
    autocorr = products * n / ((n - lags) * products[0])
    return lags, autocorr


# Save one completed MCMC figure and retain the interactive study behavior
def _save_figure(fig, filename: str, description: str):
    ensure_dir(pathlib.Path(filename).parent)
    fig.tight_layout()
    fig.savefig(filename)
    print(f"Saving {description} to: {filename}")
    plt.show()


# Plot per-parameter trace, marginal density, and autocorrelation diagnostics
def plot_parameter_traces(samples: torch.Tensor, config_keys: list = None, save_dir: str = "."):
    # Convert samples to numpy array if necessary
    samples_np = samples.numpy() if not isinstance(samples, np.ndarray) else samples

    _, n_params = samples_np.shape

    if config_keys is None:
        config_keys = [f"$\\theta_{{{i}}}$" for i in range(n_params)]

    for i in range(n_params):
        fig, axs = plt.subplots(1, 3, figsize=(18, 4))

        # Trace plot
        axs[0].plot(samples_np[:, i], lw=0.8)
        axs[0].set(title=f"Trace for {config_keys[i]}", xlabel="Iteration", ylabel="Value")

        # Histogram with median
        axs[1].hist(samples_np[:, i], bins=50, density=True, alpha=0.7)
        median = np.median(samples_np[:, i])
        axs[1].axvline(median, color="r", linestyle="--", label=f"Median: {median:.2f}")
        axs[1].set_title(f"Histogram for {config_keys[i]}")
        axs[1].legend()

        # Autocorrelation plot using matplotlib
        max_lag = 50
        lags, acorr = compute_autocorrelation(samples_np[:, i], max_lag=max_lag)
        axs[2].stem(lags, acorr)
        axs[2].set(title=f"Autocorrelation for {config_keys[i]}", xlabel="Lag", ylabel="Autocorrelation")
        axs[2].set_ylim(-1, 1)

        filename = os.path.join(save_dir, f"parameter_trace_{i}.pdf")
        _save_figure(fig, filename, "parameter trace plot")


# Plot joint posterior summaries and Hamiltonian diagnostics
def visualize_samples_and_diagnostics(samples: torch.Tensor, diagnostics: dict, savedir: str, config_keys: list = None):
    samples_np = samples.numpy()
    n_params = samples_np.shape[1]

    if config_keys is None:
        config_keys = ["$\\theta_{i}$" for _ in range(n_params)]
    # Diagnostics Plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    for i in range(min(4, n_params)):
        axes[0, 0].plot(samples_np[:, i], label=config_keys[i])
    axes[0, 0].legend()
    accept_hist = np.asarray(diagnostics["accept_history"])
    ends = np.arange(1, len(accept_hist) + 1)
    starts = np.maximum(0, ends - 100)
    cumulative = np.pad(np.cumsum(accept_hist), (1, 0))
    running_accept = (cumulative[ends] - cumulative[starts]) / (ends - starts)
    series = (
        (None, "Parameter Traces", "Sample", "Parameter Value"),
        (diagnostics["energy_history"], "Energy Trace", "Sample", "Energy"),
        (diagnostics["step_size_history"], "Step Size Adaptation", "Iteration", "Step Size"),
        (running_accept, "Running Acceptance Rate", "Sample", "Acceptance Rate"),
    )
    for axis, (values, title, xlabel, ylabel) in zip(axes.flat, series, strict=True):
        if values is not None:
            axis.plot(values)
        axis.set(title=title, xlabel=xlabel, ylabel=ylabel)
    axes[1, 1].axhline(y=diagnostics["target_accept"], color="r", linestyle="--", label="Target")
    axes[1, 1].legend()
    filename = os.path.join(savedir, "hmc_diagnostics.pdf")
    _save_figure(fig, filename, "diagnostic plots")
    # Analyze samples to compute central estimates and credibility intervals
    analysis = analyze_samples(samples=samples, config_keys=config_keys)
    # Parameter Histograms with Legends
    n_cols = min(4, n_params)
    n_rows = (n_params + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 3 * n_rows), squeeze=False)
    axes = axes.flatten()

    for i in range(n_params):
        axes[i].hist(samples_np[:, i], bins=50, density=True, alpha=0.7)
        median = analysis["median"][i]
        ci68_lower, ci68_upper = analysis["ci68"][i]
        ci95_lower, ci95_upper = analysis["ci95"][i]
        axes[i].set_title(config_keys[i])

        # Add a legend with the numerical summary
        axes[i].legend(
            [
                f"Median: {median:.2f}\n68% CI: [{ci68_lower:.2f}, {ci68_upper:.2f}]\n95% CI: [{ci95_lower:.2f}, {ci95_upper:.2f}]"
            ]
        )
    for i in range(n_params, len(axes)):
        axes[i].axis("off")

    filename = os.path.join(savedir, "parameter_distributions.pdf")
    _save_figure(fig, filename, "parameter distribution plots")
    # Matrix Plot: Diagonal (histograms) and lower triangle (scatter with correlations)
    fig, axes = plt.subplots(n_params, n_params, figsize=(3 * n_params, 3 * n_params), squeeze=False)
    for i in range(n_params):
        for j in range(n_params):
            if i == j:
                # Diagonal: Histogram with median line
                axes[i, j].hist(samples_np[:, i], bins=30, density=True, alpha=0.7)
                axes[i, j].axvline(analysis["median"][i], color="r", linestyle="--")
                axes[i, j].set_title(config_keys[i], fontsize=10)
            elif i > j:
                # Lower triangle: Scatter plot of parameter j vs parameter i
                axes[i, j].scatter(samples_np[:, j], samples_np[:, i], s=1, alpha=0.5)
                # Annotate with correlation coefficient
                corr_coef = analysis["correlation"][i, j]
                axes[i, j].annotate(
                    f"r={corr_coef:.2f}", xy=(0.05, 0.9), xycoords="axes fraction", fontsize=8, color="red"
                )
            else:
                axes[i, j].set_visible(False)

        axes[i, 0].set_ylabel(config_keys[i], fontsize=10)

    for j in range(n_params):
        axes[n_params - 1, j].set_xlabel(config_keys[j], fontsize=10)

    filename = os.path.join(savedir, "matrix_plot.pdf")
    _save_figure(fig, filename, "matrix plot")


# Compute sorted trial pickle files for an icetune run
def collect_trial_pickles(run_name: str, cdir: str) -> list[str]:
    dirname = os.path.join(cdir, "runs/icetune", run_name)
    files = glob.glob(os.path.join(dirname, "results", "TUNE_icetune_*.pkl"))
    files_sorted = sorted(files, key=os.path.getctime)
    if len(files_sorted) == 0:
        raise Exception(__name__ + f".collect_trial_pickles: Did not find any trials under {dirname}")
    return files_sorted


# Extract one histogram count and uncertainty pair from a hdata object
def _hist_counts_errs(hist_record: dict) -> tuple[np.ndarray, np.ndarray]:
    hdata = hist_record["hdata"]
    counts = np.asarray(hdata.counts_scaled, dtype=float)
    errs = np.asarray(hdata.errs_scaled, dtype=float)
    if counts.ndim != 1 or errs.shape != counts.shape:
        raise ValueError("Trial histogram counts and errors must be aligned vectors")
    return counts, errs


# Compute stable histogram keys present in both MC and data payloads
def _histogram_keys_from_results(results: dict) -> list[tuple[int, int, str]]:
    mc = results.get("mc")
    data = results.get("data")
    if mc is None or data is None:
        raise ValueError("results payload must contain 'mc' and 'data'")

    keys = []
    for i in range(len(mc)):
        if mc[i] is None or data[i] is None:
            continue
        for j in range(len(mc[i])):
            if mc[i][j] is None or data[i][j] is None:
                continue
            obs_keys = sorted(set(mc[i][j].keys()) & set(data[i][j].keys()))
            keys.extend((i, j, obs) for obs in obs_keys)
    if len(keys) == 0:
        raise ValueError("results payload does not contain any common MC/data histograms")
    return keys


# Flatten nested histogram payloads using a stable key manifest
def _flatten_histograms(source: list, keys: list[tuple[int, int, str]]) -> tuple[np.ndarray, np.ndarray]:
    counts_all = []
    errs_all = []
    for i, j, obs in keys:
        counts, errs = _hist_counts_errs(source[i][j][obs])
        counts_all.extend(counts)
        errs_all.extend(errs)
    return np.asarray(counts_all, dtype=float), np.asarray(errs_all, dtype=float)


# Build stable bin slices and plotting metadata for every histogram
def _histogram_manifest(results: dict) -> list[dict]:
    manifest = []
    offset = 0
    for dataset, subset, observable in _histogram_keys_from_results(results):
        mc_record = results["mc"][dataset][subset][observable]
        data_record = results["data"][dataset][subset][observable]
        mc_hdata = mc_record["hdata"]
        data_hdata = data_record["hdata"]
        counts, _ = _hist_counts_errs(mc_record)
        data_counts, _ = _hist_counts_errs(data_record)
        mc_valid = np.asarray(
            getattr(mc_hdata, "valid", np.ones_like(counts, dtype=np.bool_)), dtype=np.bool_
        )
        data_valid = np.asarray(getattr(data_hdata, "valid", np.ones_like(data_counts, dtype=np.bool_)), dtype=np.bool_)
        if data_counts.shape != counts.shape or mc_valid.shape != counts.shape or data_valid.shape != counts.shape:
            raise ValueError("Trial MC and data histogram bins and validity masks must align")
        valid = mc_valid & data_valid
        bins = getattr(data_hdata, "bins", None)
        mc_bins = getattr(mc_hdata, "bins", None)
        if bins is not None or mc_bins is not None:
            bins, mc_bins = np.asarray(bins, dtype=float), np.asarray(mc_bins, dtype=float)
            if (
                bins.shape != (len(counts) + 1,) or mc_bins.shape != bins.shape
                or not np.all(np.isfinite(bins)) or not np.all(np.diff(bins) > 0.0)
                or not np.allclose(mc_bins, bins, rtol=1.0e-12, atol=0.0)
            ):
                raise ValueError("Trial MC and data histogram edges must match")
        stop = offset + len(counts)
        manifest.append(
            {
                "dataset": int(dataset),
                "subset": int(subset),
                "observable": str(observable),
                "start": int(offset),
                "stop": int(stop),
                "valid": valid.tolist(),
                "fitw": float(data_record.get("fitw", 1.0)),
                "bins": None if bins is None else bins.tolist(),
                "ndf": int(np.count_nonzero(valid & ((counts != 0.0) | (data_counts != 0.0)))),
            }
        )
        offset = stop
    return manifest


# Check that every training trial represents the same physical bin definitions
def _matching_histogram_manifest(first: list[dict], second: list[dict]) -> bool:
    if len(first) != len(second):
        return False
    for left, right in zip(first, second, strict=True):
        if any(left[key] != right[key] for key in ("dataset", "subset", "observable", "start", "stop", "valid")):
            return False
        if not np.isclose(left["fitw"], right["fitw"], rtol=1.0e-12, atol=0.0):
            return False
        if left["bins"] is None or right["bins"] is None:
            if left["bins"] is not right["bins"]:
                return False
        elif not np.allclose(left["bins"], right["bins"], rtol=1.0e-12, atol=0.0):
            return False
    return True


# Extract simulator metadata and source steering from one full trial payload
def _trial_producer_metadata(payload: dict) -> dict:
    param = payload.get("param", {}) or {}
    if param.get("cost") == "gaussian":
        raise ValueError("Gaussian ampfit requires icescape history input. Histogram surrogates do not preserve "
                         "its source MC statistics and full MC covariance")
    return {
        "simdriver": copy.deepcopy(param.get("simdriver")),
        "mc_steer": copy.deepcopy(param.get("mc_steer", {})),
        "plot_brand": copy.deepcopy(param.get("plot_brand")),
        "data_covariance_mode": copy.deepcopy(param.get("data_covariance_mode", "diagonal")),
    }


# Load one current-format icetune trial pickle with full histogram results
def _load_trial_payload(filename: str) -> dict:
    with open(filename, "rb") as f:
        payload = pickle.load(f)
    if payload.get("replica_schema_version") != 1:
        raise ValueError(f"trial pickle has unsupported schema: {filename}")
    results = payload.get("results")
    if results is None:
        raise ValueError(f"trial pickle has no full 'results' payload: {filename}")
    _histogram_keys_from_results(results)
    return payload


# Load a deep-copy plotting template and its stable histogram manifest
def collect_replica_template(run_name: str, cdir: str) -> dict:
    files_sorted = collect_trial_pickles(run_name=run_name, cdir=cdir)
    for filename in files_sorted:
        try:
            payload = _load_trial_payload(filename)
        except Exception as exc:
            cprint(__name__ + f".collect_replica_template: Skipping invalid trial pickle {filename}: {exc}", "red")
            continue
        results = copy.deepcopy(payload["results"])
        return {
            "results": results,
            "histogram_manifest": _histogram_manifest(results),
            "source": filename,
            **_trial_producer_metadata(payload),
        }
    raise RuntimeError("No valid full-result trial is available for replica plotting")


def collect_simu(run_name: str, cdir: str, max_trials: int = int(1e6), use_cached: bool = False):
    """
    Collect observed data and MC predictions from pickle files.

    Parameters:
        run_name (str): Name of the simulation run.
        cdir (str): Base directory.
        max_trials (int): Maximum number of trials to load.
        use_cached (bool): If True and a cached combined file exists, load and return it.

    Returns:
        dict: A dictionary containing 'config_keys', 'X', 'Y', and 'E'.
    """

    dirname = os.path.join(cdir, "runs/icetune", run_name)
    if max_trials < 1:
        raise ValueError("Maximum trial count must be positive")
    cprint(__name__ + f".collect_simu: dirname = {dirname}", "magenta")

    # Path for the combined pickle file
    combined_pickle_path = os.path.join(dirname, "results", "combined.pkl")

    # If use_cached is True and the combined file exists, load and return it
    if use_cached and os.path.exists(combined_pickle_path):
        cprint(__name__ + f".collect_simu: Loading cached output from {combined_pickle_path}", "yellow")
        with open(combined_pickle_path, "rb") as f:
            output = pickle.load(f)
        if output.get("replica_input_schema_version") == 1:
            for name in ("X", "Y", "E"):
                output[name] = output[name][:max_trials]
            return output
        cprint(__name__ + ".collect_simu: Ignoring stale combined.pkl replica schema", "yellow")

    files_sorted = collect_trial_pickles(run_name=run_name, cdir=cdir)
    cprint(__name__ + f".collect_simu: Found {len(files_sorted)} trials under {dirname}", "yellow")
    # Gather MC predictions from all files
    config_keys = None
    key_manifest = None
    histogram_manifest = None
    producer_metadata = None
    observed_counts = None
    observed_errs = None
    samples = {name: [] for name in ("X", "Y", "E")}

    for filename in files_sorted:
        try:
            data = _load_trial_payload(filename)
        except Exception as exc:
            cprint(__name__ + f".collect_simu: Skipping invalid trial pickle {filename}: {exc}", "red")
            continue
        print(filename)

        trial_producer = _trial_producer_metadata(data)
        if producer_metadata is None:
            producer_metadata = trial_producer
        elif trial_producer != producer_metadata:
            raise ValueError("Full trial pickles contain inconsistent simulator or plot-brand metadata")

        config = data["config"]
        trial_config_keys = sorted(config.keys())
        if config_keys is None:
            config_keys = trial_config_keys
        elif trial_config_keys != config_keys:
            raise ValueError("Full trial pickles contain inconsistent parameter keys")
        x_vec = np.array([config[key] for key in config_keys], dtype=float)
        results = data["results"]
        keys = _histogram_keys_from_results(results)
        trial_manifest = _histogram_manifest(results)
        if key_manifest is None:
            key_manifest = keys
            histogram_manifest = trial_manifest
        elif keys != key_manifest:
            cprint(__name__ + f".collect_simu: Skipping trial with non-matching histogram manifest {filename}", "red")
            continue
        elif not _matching_histogram_manifest(histogram_manifest, trial_manifest):
            raise ValueError("Full trial pickles contain inconsistent histogram bins, masks or fit weights")
        trial_counts, trial_errs = _flatten_histograms(results["data"], keys)
        if observed_counts is None:
            observed_counts = trial_counts
            observed_errs = trial_errs
        elif not (
            np.allclose(trial_counts, observed_counts, rtol=1.0e-12, atol=0.0)
            and np.allclose(trial_errs, observed_errs, rtol=1.0e-12, atol=0.0)
        ):
            raise ValueError("Full trial pickles contain inconsistent observed histograms")
        y_vec, e_vec = _flatten_histograms(results["mc"], keys)

        for name, value in zip(samples, (x_vec, y_vec, e_vec), strict=True):
            samples[name].append(value)

        if len(samples["X"]) >= max_trials:
            cprint(__name__ + f".collect_simu: max_trials = {max_trials} reached, break", "red")
            break

    if len(samples["X"]) == 0:
        raise Exception(__name__ + f".collect_simu: Did not find any valid full-result trials under {dirname}")

    samples = {name: np.asarray(values) for name, values in samples.items()}
    for name, values in samples.items():
        print(f"{name} shape:", values.shape)

    output = {
        "replica_input_schema_version": 1,
        "config_keys": config_keys,
        "Y_data": observed_counts,
        "E_data": observed_errs,
        **samples,
        "histogram_manifest": histogram_manifest,
        **(producer_metadata or {"simdriver": None, "plot_brand": None, "data_covariance_mode": "diagonal"}),
    }

    # Dump the output dict into the combined pickle file
    ensure_dir(pathlib.Path(combined_pickle_path).parent)
    with open(combined_pickle_path, "wb") as f:
        cprint(__name__ + f".collect_simu: Saving simulation pickle to: {combined_pickle_path}", "yellow")
        pickle.dump(output, f)

    return output

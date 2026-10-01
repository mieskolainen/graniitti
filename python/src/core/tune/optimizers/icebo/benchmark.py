# Reproducible toy and scaling benchmarks for ICEBO
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import json
import math
import time
import warnings
from collections.abc import Callable
from contextlib import redirect_stdout
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

import numpy as np
import torch

from core import resource
from core.io.files import ensure_dir
from core.tune.optimizers.icebo.config import load_icebo_config, override_icebo_config
from core.tune.optimizers.icebo.optimizer import ICEBO

DEFAULT_ICEBO_CONFIG = resource("tune/settings/icebo.json")


# Compute benchmark settings derived from the authoritative steering JSON
def benchmark_icebo_settings(overrides: dict[str, object]) -> dict:
    settings = load_icebo_config(DEFAULT_ICEBO_CONFIG)
    return override_icebo_config(settings, overrides)


@dataclass(frozen=True)
class ToyProblem:
    """Bounded scalar minimization problem with an optional noise model"""

    name: str
    bounds: tuple[tuple[float, float], ...]
    optimum: float
    function: Callable[[np.ndarray], float]
    noise_std: float = 0.0
    noise_radial_scale: float = 0.0

    # Validate the benchmark domain and observation-noise model
    def __post_init__(self) -> None:
        if not self.bounds:
            raise ValueError("Benchmark problem bounds cannot be empty")
        if any(not lower < upper for lower, upper in self.bounds):
            raise ValueError("Benchmark problem bounds must have positive width")
        if self.noise_std < 0.0 or self.noise_radial_scale < 0.0:
            raise ValueError("Benchmark noise settings must be non-negative")

    # Compute the normalized dictionary space used by proposal optimizers
    def parameter_space(self) -> dict:
        return {
            f"x{index}": {"type": "uniform", "lower": lower, "upper": upper}
            for index, (lower, upper) in enumerate(self.bounds)
        }

    # Evaluate the latent objective at one configuration
    def latent(self, config: dict) -> float:
        values = np.asarray([config[f"x{index}"] for index in range(len(self.bounds))])
        return float(self.function(values))

    # Compute the configured input-dependent observation standard deviation
    def observation_std(self, config: dict) -> float:
        if self.noise_std == 0.0:
            return 0.0
        values = np.asarray([config[f"x{index}"] for index in range(len(self.bounds))])
        lower = np.asarray([item[0] for item in self.bounds])
        upper = np.asarray([item[1] for item in self.bounds])
        centered = 2.0 * (values - lower) / (upper - lower) - 1.0
        radius = float(np.sqrt(np.mean(centered**2)))
        return self.noise_std * (1.0 + self.noise_radial_scale * radius)


@dataclass(frozen=True)
class ObjectiveEvaluation:
    """One latent and observed benchmark objective evaluation"""

    latent: float
    observed: float
    standard_deviation: float


# Compute the standard Branin test function
def branin(values: np.ndarray) -> float:
    first, second = values
    term = second - 5.1 * first**2 / (4.0 * math.pi**2) + 5.0 * first / math.pi - 6.0
    return term**2 + 10.0 * (1.0 - 1.0 / (8.0 * math.pi)) * math.cos(first) + 10.0


# Compute the six-dimensional Hartmann test function
def hartmann6(values: np.ndarray) -> float:
    alpha = np.asarray([1.0, 1.2, 3.0, 3.2])
    matrix = np.asarray(
        [
            [10.0, 3.0, 17.0, 3.5, 1.7, 8.0],
            [0.05, 10.0, 17.0, 0.1, 8.0, 14.0],
            [3.0, 3.5, 1.7, 10.0, 17.0, 8.0],
            [17.0, 8.0, 0.05, 10.0, 0.1, 14.0],
        ]
    )
    centers = 1.0e-4 * np.asarray(
        [
            [1312, 1696, 5569, 124, 8283, 5886],
            [2329, 4135, 8307, 3736, 1004, 9991],
            [2348, 1451, 3522, 2883, 3047, 6650],
            [4047, 8828, 8732, 5743, 1091, 381],
        ]
    )
    return -float(np.sum(alpha * np.exp(-np.sum(matrix * (values[None, :] - centers) ** 2, axis=1))))


# Compute the d-dimensional Levy test function
def levy(values: np.ndarray) -> float:
    transformed = 1.0 + (values - 1.0) / 4.0
    first = math.sin(math.pi * transformed[0]) ** 2
    middle = np.sum((transformed[:-1] - 1.0) ** 2 * (1.0 + 10.0 * np.sin(math.pi * transformed[:-1] + 1.0) ** 2))
    last = (transformed[-1] - 1.0) ** 2 * (1.0 + math.sin(2.0 * math.pi * transformed[-1]) ** 2)
    return float(first + middle + last)


# Compute the d-dimensional Ackley test function
def ackley(values: np.ndarray) -> float:
    squared = np.mean(values**2)
    cosine = np.mean(np.cos(2.0 * math.pi * values))
    return float(-20.0 * np.exp(-0.2 * np.sqrt(squared)) - np.exp(cosine) + 20.0 + math.e)


# Compute the d-dimensional Rastrigin test function
def rastrigin(values: np.ndarray) -> float:
    return float(10.0 * len(values) + np.sum(values**2 - 10.0 * np.cos(2.0 * math.pi * values)))


# Compute the fixed toy suite used for paired optimizer comparisons
def toy_problems() -> list[ToyProblem]:
    return [
        ToyProblem("branin2", ((-5.0, 10.0), (0.0, 15.0)), 0.397887, branin),
        ToyProblem("hartmann6", ((0.0, 1.0),) * 6, -3.322368, hartmann6),
        ToyProblem("levy8", ((-10.0, 10.0),) * 8, 0.0, levy),
        ToyProblem("hartmann6_noisy", ((0.0, 1.0),) * 6, -3.322368, hartmann6, noise_std=0.05, noise_radial_scale=1.0),
        ToyProblem("ackley10_noisy", ((-5.0, 5.0),) * 10, 0.0, ackley, noise_std=0.2, noise_radial_scale=1.0),
        ToyProblem("rastrigin12_noisy", ((-5.12, 5.12),) * 12, 0.0, rastrigin, noise_std=1.0, noise_radial_scale=1.0),
        ToyProblem("levy16_noisy", ((-10.0, 10.0),) * 16, 0.0, levy, noise_std=0.5, noise_radial_scale=1.0),
    ]


# Select benchmark problems by their stable names
def select_toy_problems(names: list[str] | None) -> list[ToyProblem]:
    available = {problem.name: problem for problem in toy_problems()}
    if names is None:
        return list(available.values())
    unknown = sorted(set(names) - set(available))
    if unknown:
        raise ValueError(f"Unknown benchmark problems: {unknown}")
    if not names:
        raise ValueError("At least one benchmark problem is required")
    return [available[name] for name in names]


# Compute a stable integer tag for one named benchmark problem
def _problem_seed_tag(name: str) -> int:
    return sum((index + 1) * ord(character) for index, character in enumerate(name))


# Evaluate one paired noisy objective observation
def evaluate_observation(problem: ToyProblem, config: dict, *, seed: int, evaluation_index: int) -> ObjectiveEvaluation:
    if evaluation_index < 0:
        raise ValueError("Benchmark evaluation index must be non-negative")
    latent = problem.latent(config)
    standard_deviation = problem.observation_std(config)
    if standard_deviation == 0.0:
        return ObjectiveEvaluation(latent, latent, 0.0)
    sequence = np.random.SeedSequence([int(seed) & 0xFFFFFFFF, int(evaluation_index), _problem_seed_tag(problem.name)])
    fluctuation = float(np.random.default_rng(sequence).normal())
    return ObjectiveEvaluation(latent, latent + standard_deviation * fluctuation, standard_deviation)


# Evaluate the common initial design with a paired noise stream
def evaluate_initial_design(problem: ToyProblem, initial: list[dict], seed: int) -> list[ObjectiveEvaluation]:
    return [
        evaluate_observation(problem, config, seed=seed, evaluation_index=index) for index, config in enumerate(initial)
    ]


# Compute one common scrambled Sobol initial design
def initial_design(problem: ToyProblem, count: int, seed: int) -> list[dict]:
    engine = torch.quasirandom.SobolEngine(len(problem.bounds), scramble=True, seed=int(seed))
    unit = engine.draw(count).double().numpy()
    configs = []
    for row in unit:
        configs.append(
            {f"x{index}": lower + row[index] * (upper - lower) for index, (lower, upper) in enumerate(problem.bounds)}
        )
    return configs


# Compute the declared ICEBO numerical profile for one benchmark mode
def icebo_benchmark_profile(fast: bool) -> dict:
    if not fast:
        return {}
    return {
        "exact_gp_steps": 30,
        "exact_gp_restarts": 1,
        "raw_samples": 512,
        "acquisition_restarts": 8,
        "acquisition_steps": 20,
        "mc_samples": 128,
    }


# Run ICEBO after observing a common initial design
def run_icebo(problem: ToyProblem, initial: list[dict], budget: int, seed: int, *, fast: bool, batch_size: int = 1) -> dict:
    evaluations = evaluate_initial_design(problem, initial, seed)
    profile = icebo_benchmark_profile(fast)
    overrides = {}
    if fast:
        overrides = {
            "gp.fit_steps": profile["exact_gp_steps"],
            "gp.fit_restarts": profile["exact_gp_restarts"],
            "acquisition.raw_samples": profile["raw_samples"],
            "acquisition.restarts": profile["acquisition_restarts"],
            "acquisition.gradient_steps": profile["acquisition_steps"],
            "acquisition.mc_samples": profile["mc_samples"],
        }
    settings = benchmark_icebo_settings(overrides)
    optimizer = ICEBO(problem.parameter_space(), settings=settings, seed=seed, warmup=0)
    optimizer.observe(initial, [item.observed for item in evaluations], update_turbo=False)
    proposal_seconds = 0.0
    while len(evaluations) < budget:
        start = time.perf_counter()
        configs = optimizer.suggest(min(batch_size, budget - len(evaluations)))
        proposal_seconds += time.perf_counter() - start
        batch = [evaluate_observation(problem, config, seed=seed, evaluation_index=len(evaluations) + index)
                 for index, config in enumerate(configs)]
        optimizer.observe(configs, [item.observed for item in batch])
        evaluations.extend(batch)
    return _run_record("icebo", problem, seed, evaluations, proposal_seconds, len(initial), batch_size)


# Run native HEBO after observing the common initial design
def run_hebo(problem: ToyProblem, initial: list[dict], budget: int, seed: int, *, batch_size: int = 1) -> dict:
    import pandas as pd
    from hebo.design_space.design_space import DesignSpace
    from hebo.optimizers.hebo import HEBO

    captured = StringIO()
    with redirect_stdout(captured):
        np.random.seed(int(seed))
        torch.manual_seed(int(seed))
        design = DesignSpace().parse(
            [
                {"name": name, "type": "num", "lb": spec["lower"], "ub": spec["upper"]}
                for name, spec in problem.parameter_space().items()
            ]
        )
        optimizer = HEBO(design, rand_sample=0, scramble_seed=int(seed))
        evaluations = evaluate_initial_design(problem, initial, seed)
        optimizer.observe(pd.DataFrame(initial), np.asarray([item.observed for item in evaluations])[:, None])
        proposal_seconds = 0.0
        while len(evaluations) < budget:
            start = time.perf_counter()
            proposals = optimizer.suggest(min(batch_size, budget - len(evaluations)))
            proposal_seconds += time.perf_counter() - start
            batch = [evaluate_observation(problem, config, seed=seed, evaluation_index=len(evaluations) + index)
                     for index, config in enumerate(proposals.to_dict("records"))]
            optimizer.observe(proposals, np.asarray([item.observed for item in batch])[:, None])
            evaluations.extend(batch)
    record = _run_record("hebo", problem, seed, evaluations, proposal_seconds, len(initial), batch_size)
    messages = captured.getvalue()
    record["numerics"] = {
        "gp_fit_failures": messages.count("jitter is too large, give up fitting GP"),
        "random_prediction_fallbacks": messages.count("jitter is too large, output random predictions"),
    }
    return record


# Run Optuna TPE after installing the common initial design
def run_optuna(problem: ToyProblem, initial: list[dict], budget: int, seed: int, *, batch_size: int = 1) -> dict:
    import optuna
    from optuna.distributions import FloatDistribution

    optuna.logging.set_verbosity(optuna.logging.ERROR)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", optuna.exceptions.ExperimentalWarning)
        sampler = optuna.samplers.TPESampler(seed=int(seed), n_startup_trials=0, multivariate=True, group=True,
                                            constant_liar=batch_size > 1)
    study = optuna.create_study(direction="minimize", sampler=sampler)
    distributions = {
        f"x{index}": FloatDistribution(lower, upper) for index, (lower, upper) in enumerate(problem.bounds)
    }
    evaluations = evaluate_initial_design(problem, initial, seed)
    for config, evaluation in zip(initial, evaluations, strict=True):
        study.add_trial(
            optuna.trial.create_trial(params=config, distributions=distributions, value=evaluation.observed)
        )
    proposal_seconds = 0.0
    while len(evaluations) < budget:
        start = time.perf_counter()
        trials = [study.ask(fixed_distributions=distributions) for _ in range(min(batch_size, budget - len(evaluations)))]
        proposal_seconds += time.perf_counter() - start
        for trial in trials:
            evaluation = evaluate_observation(problem, trial.params, seed=seed, evaluation_index=len(evaluations))
            study.tell(trial, evaluation.observed)
            evaluations.append(evaluation)
    return _run_record("optuna", problem, seed, evaluations, proposal_seconds, len(initial), batch_size)


# Build one JSON-safe run result
def _run_record(
    optimizer: str,
    problem: ToyProblem,
    seed: int,
    evaluations: list[ObjectiveEvaluation],
    proposal_seconds: float,
    initial_count: int,
    batch_size: int = 1,
) -> dict:
    latent_values = np.asarray([item.latent for item in evaluations])
    observed_values = np.asarray([item.observed for item in evaluations])
    latent_trace = np.minimum.accumulate(latent_values)
    observed_trace = np.minimum.accumulate(observed_values)
    regret = max(float(latent_trace[-1] - problem.optimum), 0.0)
    initial_regret = max(float(np.min(latent_values[:initial_count]) - problem.optimum), 1.0e-12)
    return {
        "optimizer": optimizer,
        "problem": problem.name,
        "dimension": len(problem.bounds),
        "noise_std": problem.noise_std,
        "noise_radial_scale": problem.noise_radial_scale,
        "seed": int(seed),
        "budget": len(evaluations),
        "initial_points": int(initial_count),
        "batch_size": int(batch_size),
        "best": float(latent_trace[-1]),
        "best_observed": float(observed_trace[-1]),
        "optimum": problem.optimum,
        "regret": regret,
        "normalized_regret": regret / initial_regret,
        "proposal_seconds": float(proposal_seconds),
        "latent_trace": [float(value) for value in latent_trace],
        "observed_trace": [float(value) for value in observed_trace],
        "observation_std": [float(item.standard_deviation) for item in evaluations],
    }


# Run paired toy comparisons with common initial points
def benchmark_optimizers(
    *,
    optimizers: list[str],
    seeds: list[int],
    budget: int,
    initial_points: int | None,
    fast: bool,
    problem_names: list[str] | None = None,
    batch_size: int = 1,
) -> dict:
    runners = {"icebo": run_icebo, "hebo": run_hebo, "optuna": run_optuna}
    if not optimizers:
        raise ValueError("At least one benchmark optimizer is required")
    if not seeds:
        raise ValueError("At least one benchmark seed is required")
    if not isinstance(batch_size, int) or isinstance(batch_size, bool) or batch_size < 1:
        raise ValueError("Benchmark batch size must be a positive integer")
    unknown = sorted(set(optimizers) - set(runners))
    if unknown:
        raise ValueError(f"Unknown benchmark optimizers: {unknown}")
    runs = []
    problems = select_toy_problems(problem_names)
    for problem in problems:
        initial_count = (
            min(budget - 1, max(4, len(problem.bounds) + 1)) if initial_points is None else int(initial_points)
        )
        if not 1 <= initial_count < budget:
            raise ValueError("Benchmark initial point count must be below the budget")
        for seed in seeds:
            initial = initial_design(problem, initial_count, seed)
            for name in optimizers:
                runner = runners[name]
                if name == "icebo":
                    runs.append(runner(problem, initial, budget, seed, fast=fast, batch_size=batch_size))
                else:
                    runs.append(runner(problem, initial, budget, seed, batch_size=batch_size))
    return {
        "schema_version": 1,
        "settings": {
            "optimizers": list(optimizers),
            "problems": [problem.name for problem in problems],
            "seeds": [int(seed) for seed in seeds],
            "budget": int(budget),
            "initial_points": initial_points,
            "batch_size": batch_size,
            "fast": bool(fast),
            "icebo_config": str(DEFAULT_ICEBO_CONFIG),
            "icebo_numerics": icebo_benchmark_profile(fast),
            "noise_protocol": {
                "draws": "paired_standard_normal_by_evaluation_index",
                "icebo_receives_known_standard_deviation": False,
                "metric": "best_latent_value_at_evaluated_configurations",
            },
        },
        "runs": runs,
        "summary": summarize_runs(runs),
        "problem_summary": summarize_problems(runs),
    }


# Aggregate paired normalized regret, runtime, and win rates
def summarize_runs(runs: list[dict]) -> dict:
    names = sorted({run["optimizer"] for run in runs})
    summary = {}
    for name in names:
        selected = [run for run in runs if run["optimizer"] == name]
        regrets = np.asarray([run["normalized_regret"] for run in selected])
        times = np.asarray([run["proposal_seconds"] for run in selected])
        summary[name] = {
            "runs": len(selected),
            "median_regret": float(np.median([run["regret"] for run in selected])),
            "median_normalized_regret": float(np.median(regrets)),
            "geometric_mean_normalized_regret": float(np.exp(np.mean(np.log(np.maximum(regrets, 1.0e-12))))),
            "median_proposal_seconds": float(np.median(times)),
            "wins": 0,
            "win_rate": 0.0,
            "gp_fit_failures": sum(run.get("numerics", {}).get("gp_fit_failures", 0) for run in selected),
            "random_prediction_fallbacks": sum(
                run.get("numerics", {}).get("random_prediction_fallbacks", 0) for run in selected
            ),
        }
    groups = {}
    for run in runs:
        groups.setdefault((run["problem"], run["seed"]), []).append(run)
    for group in groups.values():
        best = min(item["normalized_regret"] for item in group)
        for item in group:
            if math.isclose(item["normalized_regret"], best, rel_tol=1.0e-12, abs_tol=1.0e-12):
                summary[item["optimizer"]]["wins"] += 1
    comparison_count = len(groups)
    for values in summary.values():
        values["win_rate"] = values["wins"] / comparison_count
    return summary


# Aggregate paired benchmark metrics separately for every problem
def summarize_problems(runs: list[dict]) -> dict:
    output = {}
    problem_names = sorted({run["problem"] for run in runs})
    optimizer_names = sorted({run["optimizer"] for run in runs})
    for problem_name in problem_names:
        problem_runs = [run for run in runs if run["problem"] == problem_name]
        seed_groups = {}
        for run in problem_runs:
            seed_groups.setdefault(run["seed"], []).append(run)
        problem_output = {}
        for optimizer_name in optimizer_names:
            selected = [run for run in problem_runs if run["optimizer"] == optimizer_name]
            if not selected:
                continue
            normalized = np.asarray([run["normalized_regret"] for run in selected])
            wins = 0
            for group in seed_groups.values():
                candidate = next((run for run in group if run["optimizer"] == optimizer_name), None)
                if candidate is None:
                    continue
                best = min(run["normalized_regret"] for run in group)
                if math.isclose(candidate["normalized_regret"], best, rel_tol=1.0e-12, abs_tol=1.0e-12):
                    wins += 1
            problem_output[optimizer_name] = {
                "runs": len(selected),
                "median_regret": float(np.median([run["regret"] for run in selected])),
                "median_normalized_regret": float(np.median(normalized)),
                "normalized_regret_q25": float(np.quantile(normalized, 0.25)),
                "normalized_regret_q75": float(np.quantile(normalized, 0.75)),
                "median_proposal_seconds": float(np.median([run["proposal_seconds"] for run in selected])),
                "wins": wins,
                "gp_fit_failures": sum(run.get("numerics", {}).get("gp_fit_failures", 0) for run in selected),
                "random_prediction_fallbacks": sum(
                    run.get("numerics", {}).get("random_prediction_fallbacks", 0) for run in selected
                ),
            }
        output[problem_name] = problem_output
    return output


# Benchmark exact GP fit and prediction throughput as N grows
def benchmark_scaling(*, sizes: list[int], dimension: int, seed: int, device: str) -> list[dict]:
    rng = np.random.default_rng(seed)
    bounds = {f"x{index}": {"type": "uniform", "lower": -1.0, "upper": 1.0} for index in range(dimension)}
    output = []
    for size in sizes:
        values = rng.uniform(-1.0, 1.0, size=(size, dimension))
        target = np.sum(values**2, axis=1) + 0.1 * np.sin(5.0 * values[:, 0])
        configs = [{f"x{index}": row[index] for index in range(dimension)} for row in values]
        settings = benchmark_icebo_settings(
            {
                "search.device": device,
                "gp.fit_steps": 30,
                "gp.fit_restarts": 1,
                "acquisition.raw_samples": 512,
                "acquisition.restarts": 8,
                "acquisition.gradient_steps": 10,
                "acquisition.mc_samples": 128,
            }
        )
        optimizer = ICEBO(bounds, settings=settings, seed=seed, warmup=0)
        optimizer.observe(configs, target)
        start = time.perf_counter()
        optimizer.predict(configs[: min(1024, size)])
        fit_predict_seconds = time.perf_counter() - start
        proposal_start = time.perf_counter()
        optimizer.suggest()
        proposal_seconds = time.perf_counter() - proposal_start
        diagnostics = optimizer.diagnostics()
        output.append(
            {
                "N": size,
                "dimension": dimension,
                "device": diagnostics["device"],
                "surrogate": diagnostics["fit"]["surrogate"],
                "fit_and_predict_seconds": fit_predict_seconds,
                "proposal_seconds": proposal_seconds,
            }
        )
    return output


# Save benchmark results as indented deterministic JSON
def save_benchmark(result: dict, path: str | Path) -> None:
    output = Path(path)
    ensure_dir(output.parent)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

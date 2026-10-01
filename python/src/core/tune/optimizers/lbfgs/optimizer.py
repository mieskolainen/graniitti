# L-BFGS backend for icetune
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import copy
import os
import pathlib
import tempfile
import time
import uuid
from dataclasses import dataclass
from functools import partial

import numpy as np
import ray
import torch

from core.io import logger as log
from core.io.files import ensure_dir
from core.io.serialize import json_safe, write_json_file
from core.numerics import lbfgsb
from core.tune import core as icetune_main
from core.tune import history as icetune_history
from core.tune import summary as fit_summary
from core.tune.io import time_fields, trial_times
from core.tune.parameters.space import normalize_continuous_param_space as normalize_param_space
from core.tune.runtime import process as iceruntime
from core.tune.runtime.ray import worker_runtime

logger = log.get_logger(__name__)


# Convert a parameter vector into a config dictionary
def unpack_config(values, names: list[str]) -> dict[str, float]:
    arr = np.asarray(values, dtype=float)
    if arr.size != len(names):
        raise ValueError("Parameter vector length does not match parameter names")
    return {name: float(arr[i]) for i, name in enumerate(names)}


# Choose the LBFGS start config from initial points or bounds midpoints
def initial_config(
    bounds: dict[str, dict[str, float]], initial_points: dict | None, no_initial_point: bool
) -> dict[str, float]:
    names = sorted(bounds)
    if no_initial_point or initial_points is None:
        return {name: 0.5 * (bounds[name]["lower"] + bounds[name]["upper"]) for name in names}
    config = {name: float(initial_points[name]) for name in names}
    for name, value in config.items():
        lower = bounds[name]["lower"]
        upper = bounds[name]["upper"]
        if value < lower or value > upper:
            raise ValueError(f'Initial value for "{name}" = {value} is outside lbfgs bounds [{lower}, {upper}]')
    return config


# Solve finite-difference weights for the first derivative at x0
def finite_difference_weights(x0: float, points: list[float]) -> np.ndarray:
    offsets = np.array([float(p) - float(x0) for p in points], dtype=float)
    matrix = np.vstack([offsets**k for k in range(len(points))])
    rhs = np.zeros(len(points), dtype=float)
    rhs[1] = 1.0
    return np.linalg.solve(matrix, rhs)


@dataclass
class GradientStencil:
    """Finite-difference stencil for one parameter derivative"""

    index: int
    points: list[float]
    weights: list[float]


# Build bounded finite-difference stencils for all parameters
def build_gradient_stencils(
    values, bounds: dict[str, dict[str, float]], names: list[str], rel_step: float, abs_step: float
) -> list[GradientStencil]:
    x = np.asarray(values, dtype=float)
    if x.size != len(names):
        raise ValueError("Gradient vector length does not match parameter names")

    rel_step = float(rel_step)
    abs_step = float(abs_step)
    if rel_step <= 0.0 or abs_step <= 0.0:
        raise ValueError("lbfgs finite-difference steps must be positive")

    stencils = []
    for i, name in enumerate(names):
        lower = float(bounds[name]["lower"])
        upper = float(bounds[name]["upper"])
        span = upper - lower
        xi = min(max(float(x[i]), lower), upper)
        h = max(abs_step, rel_step * max(abs(xi), span))
        h = min(h, 0.5 * span)

        if lower <= xi - h and xi + h <= upper:
            points = [xi - h, xi + h]
        elif xi + 2.0 * h <= upper:
            points = [xi, xi + h, xi + 2.0 * h]
        elif lower <= xi - 2.0 * h:
            points = [xi - 2.0 * h, xi - h, xi]
        elif xi < upper:
            right = min(upper, xi + h)
            points = [xi, right]
        elif lower < xi:
            left = max(lower, xi - h)
            points = [left, xi]
        else:
            raise ValueError(f'Cannot build finite-difference stencil for "{name}"')

        if len(set(points)) != len(points):
            raise ValueError(f'Degenerate finite-difference stencil for "{name}"')
        weights = finite_difference_weights(xi, points)
        stencils.append(
            GradientStencil(index=i, points=[float(p) for p in points], weights=[float(w) for w in weights])
        )
    return stencils


# Evaluate a finite-difference gradient from precomputed function values
def finite_difference_gradient_from_values(
    stencils: list[GradientStencil], values_by_point: dict[tuple[int, float], float]
) -> np.ndarray:
    grad = np.zeros(len(stencils), dtype=float)
    for stencil in stencils:
        total = 0.0
        for point, weight in zip(stencil.points, stencil.weights, strict=False):
            total += weight * float(values_by_point[(stencil.index, point)])
        grad[stencil.index] = total
    return grad


# Compute a stable cache key for a config and random seed
def _cache_key(config: dict, seed: int) -> tuple:
    return (int(seed), tuple((key, float(config[key]).hex()) for key in sorted(config)))


# Rerun and render the selected fit entirely in a Ray worker
@ray.remote
def _render_best_remote(config: dict, param: dict, simdriver, metrics: dict, seed: int) -> dict | None:
    from core.tune.backends import ray as icetune_ray

    param = worker_runtime(param, simdriver)
    iceruntime.set_random_seeds(seed)
    trial_id = "lbfgs-final-best"
    with icetune_main.trial_evaluation(config=config, param=param, simdriver=simdriver, trial_id=trial_id) as outputs:
        if outputs.get("results") is None:
            return None
        summary = icetune_main.build_trial_summary_payload(outputs=outputs, param=param, metrics=metrics)
        tmp_root = pathlib.Path(param["cdir"]) / "tmp"
        ensure_dir(tmp_root)
        with tempfile.TemporaryDirectory(prefix="lbfgs-figures-", dir=tmp_root) as stage:
            published = icetune_main._render_figure_candidate(
                tmp_dir=stage,
                run_name=param["run_name"],
                render_fn=partial(
                    simdriver.render_trial_figures_to_dir, outputs=outputs, param=param, summary_payload=summary
                ),
            )
            return {"summary": published, "figures": icetune_ray.pack_ray_figure_tree(pathlib.Path(stage))}


# Transfer selected worker figures to the head without invoking the simulator
def _publish_best(*, args, config: dict, param: dict, simdriver, metrics: dict) -> dict | None:
    from core.tune.backends import ray as icetune_ray

    result = ray.get(
        _render_best_remote.options(
            num_cpus=max(1, int(getattr(args, "cpu_per_trial", 1))),
            num_gpus=max(0, int(getattr(args, "gpu_per_trial", 0))),
        ).remote(config, param, simdriver, metrics, int(args.rngseed))
    )
    if result is None:
        return None
    tmp_root = pathlib.Path(param["cdir"]) / "tmp"
    ensure_dir(tmp_root)
    with tempfile.TemporaryDirectory(prefix="lbfgs-transfer-", dir=tmp_root) as stage:
        source = icetune_ray.extract_ray_figure_tree(payload=result["figures"], stage=pathlib.Path(stage))
        target = pathlib.Path(icetune_main.icetune_figure_dir(cdir=param["cdir"], run_name=param["run_name"]))
        icetune_ray._publish_ray_figure_tree(source=source, target=target, param=param)
    return result["summary"]


# Evaluate one simulator config in a Ray worker
@ray.remote
def _evaluate_config_remote(config: dict, param: dict, simdriver, trial_id: str, seed: int) -> dict:
    started = time.time()
    outputs = None
    error = None
    likelihood = None
    theta_hash = None
    pickle_payload = None
    try:
        with iceruntime.trial_limit(time.monotonic() + float(param["max_t"])):
            param = worker_runtime(param, simdriver)
            iceruntime.set_random_seeds(int(seed))
            with icetune_main.trial_evaluation(
                config=config, param=param, simdriver=simdriver, trial_id=trial_id
            ) as outputs:
                metrics = copy.deepcopy(outputs["metrics"])
                cost = icetune_main.require_finite_trial_cost(outputs=outputs, cost_key=param["cost"])
                card_config = copy.deepcopy(outputs.get("card_config"))
                result_config = outputs.get("config", config)
                theta = copy.deepcopy(outputs.get("theta"))
                likelihood = copy.deepcopy(outputs.get("likelihood"))
                theta_hash = copy.deepcopy(outputs.get("theta_hash"))
                tunename = outputs.get("tunename", trial_id)
                if param.get("pickle_dump"):
                    pickle_payload = icetune_main._trial_pickle_payload(outputs=outputs, param=param)
                success = True
    except Exception as exc:
        error = str(exc)
        cost = float(np.inf)
        card_config = None
        metrics = {param.get("cost", "loss"): cost}
        result_config = config
        theta = None
        tunename = trial_id
        success = False
    completed = time.time()
    return {
        "card_config": card_config,
        "config": copy.deepcopy(result_config),
        "cost": cost,
        **trial_times(started, completed),
        "error": error,
        "metrics": metrics,
        "likelihood": likelihood,
        **({"pickle_payload": pickle_payload} if pickle_payload is not None else {}),
        "node_id": str(param.get("node_id") or "ray"),
        "search_payload": {"kind": "lbfgs"},
        "seed": int(seed),
        "success": success,
        "theta": theta,
        "theta_hash": theta_hash,
        "trial_id": trial_id,
        "tunename": tunename,
    }


# Stop before launching trials beyond the configured evaluation budget
class EvaluationLimit(RuntimeError):
    pass


class DistributedLBFGSObjective:
    """LBFGS objective wrapper with Ray-distributed finite-difference gradients"""

    # Initialize the objective state, cache, and diagnostics files
    def __init__(
        self, *, args, bounds: dict[str, dict[str, float]], names: list[str], param: dict, simdriver, output_dir: str
    ):
        self.args = args
        self.bounds = bounds
        self.names = names
        self.param = copy.deepcopy(param)
        self.param["plot"] = False
        self.simdriver = simdriver
        self.cache: dict[tuple, dict] = {}
        self.records: list[dict] = []
        self.best_record: dict | None = None
        self.eval_counter = 0
        self.history_path = os.path.join(output_dir, "history.json")
        ensure_dir(output_dir)

    # Convert a vector into a bounded config
    def _config_from_values(self, values) -> dict[str, float]:
        config = unpack_config(values, self.names)
        return {
            name: min(max(value, self.bounds[name]["lower"]), self.bounds[name]["upper"])
            for name, value in config.items()
        }

    # Record one request and update the best successful result
    def _observe_result(self, record: dict) -> None:
        self.records.append(copy.deepcopy(record))
        if not record.get("success"):
            return
        if not np.isfinite(float(record["cost"])):
            return
        if self.best_record is None or float(record["cost"]) < float(self.best_record["cost"]):
            self.best_record = copy.deepcopy(record)

    # Compute the current best cost as a compact string for progress logs
    def _best_cost_text(self) -> str:
        if self.best_record is None:
            return "none"
        return f"{float(self.best_record['cost']):.6g}"

    # Evaluate several configs with bounded Ray concurrency and cache reuse
    def _evaluate_many(self, requests: list[dict]) -> list[dict]:
        results: list[dict | None] = [None] * len(requests)
        fresh: dict[tuple, list[tuple[int, dict]]] = {}

        for i, request in enumerate(requests):
            key = _cache_key(request["config"], request["seed"])
            if key in self.cache:
                record = copy.deepcopy(self.cache[key])
                record["cached"] = True
                record["label"] = request["label"]
                results[i] = record
            else:
                fresh.setdefault(key, []).append((i, request))

        batch_size = max(1, int(getattr(self.args, "max_concurrent_trials", 1)))
        unique = list(fresh.items())
        for start in range(0, len(unique), batch_size):
            remaining = int(self.args.num_trials) - len(self.cache)
            if remaining <= 0:
                for record in results:
                    if record is not None:
                        self._observe_result(record)
                self.publish_history(plot=False)
                raise EvaluationLimit
            chunk = unique[start : start + min(batch_size, remaining)]
            refs = []
            for _, group in chunk:
                request = group[0][1]
                refs.append(
                    _evaluate_config_remote.options(
                        num_cpus=max(1, int(getattr(self.args, "cpu_per_trial", 1))),
                        num_gpus=max(0, int(getattr(self.args, "gpu_per_trial", 0))),
                    ).remote(request["config"], self.param, self.simdriver, request["trial_id"], int(request["seed"]))
                )

            for (key, group), record in zip(chunk, ray.get(refs), strict=False):
                payload = record.pop("pickle_payload", None)
                if payload is not None:
                    # Persist worker histograms on the head before adding the compact history record
                    payload["param"] = copy.deepcopy(self.param)
                    destination = icetune_main.trial_pickle_dir(self.param) / icetune_main.trial_pickle_filename(
                        payload
                    )
                    icetune_main._publish_trial_pickle(
                        payload=payload, param=self.param, destination=destination, verify=True
                    )
                request = group[0][1]
                record = {**record, "cached": False, "label": request["label"]}
                record.setdefault("node_id", str(self.param.get("node_id") or "ray"))
                record.setdefault("search_payload", {"kind": "lbfgs"})
                self.cache[key] = copy.deepcopy(record)
                for duplicate, (index, item) in enumerate(group):
                    observed = copy.deepcopy(record)
                    observed["cached"] = bool(duplicate)
                    observed["label"] = item["label"]
                    results[index] = observed

        observed = [record for record in results if record is not None]
        if any(record is None for record in results):
            for record in observed:
                self._observe_result(record)
            self.publish_history(plot=False)
            raise EvaluationLimit
        for record in observed:
            self._observe_result(record)
        self.publish_history(plot=False)
        return observed

    # Write the standard completed trial history and its cost evolution plot
    def publish_history(self, plot=True) -> str | None:
        updated = time.time()
        optimization = copy.deepcopy(self.param.get("optimization") or {})
        optimization.update({"backend": "ray", "cost": self.param["cost"], "optimizer": "lbfgs"})
        trials = [
            record
            for record in self.records
            if not record.get("cached") and record.get("success") and np.isfinite(float(record["cost"]))
        ]
        history = {
            "history_schema_version": icetune_history.HISTORY_SCHEMA_VERSION,
            "optimization": optimization,
            "parameter_space": self.param.get("parameter_space") or icetune_main.build_parameter_space(self.bounds),
            "parameter_topology": copy.deepcopy(self.param.get("parameter_topology") or {}),
            "mc_steer": copy.deepcopy(self.param.get("mc_steer", {})),
            "plot_brand": self.param.get("plot_brand"),
            "requests": json_safe(self.records),
            "simdriver": self.param.get("simdriver"),
            "trials": json_safe(trials),
            **time_fields("updated_at", updated),
        }
        write_json_file(self.history_path, json_safe(history), indent=2, sort_keys=True, newline=True)
        if not plot:
            return None
        return str(
            icetune_history.create_history_plot(
                source=pathlib.Path(self.history_path),
                cdir=pathlib.Path(self.param["cdir"]),
                cost=self.param["cost"],
                run_name=self.param["run_name"],
            )
        )

    # Evaluate the scalar objective for LBFGS
    def value(self, *values) -> float:
        config = self._config_from_values(values)
        eval_id = self.eval_counter
        trial_id = f"lbfgs-f-{eval_id:06d}-{uuid.uuid4().hex[:8]}"
        self.eval_counter += 1
        record = self._evaluate_many(
            [{"config": config, "label": "objective", "seed": int(self.args.rngseed), "trial_id": trial_id}]
        )[0]
        cost = float(record["cost"])
        logger.info(
            ".icetune_lbfgs: objective eval=%06d cost=%.6g best=%s cached=%s trial_id=%s",
            eval_id,
            cost,
            self._best_cost_text(),
            bool(record.get("cached")),
            trial_id,
        )
        return cost

    # Evaluate the finite-difference gradient for LBFGS
    def gradient(self, *values) -> np.ndarray:
        x = np.asarray(values, dtype=float)
        stencils = build_gradient_stencils(
            x,
            self.bounds,
            self.names,
            rel_step=float(self.args.lbfgs_settings["grad_step_rel"]),
            abs_step=float(self.args.lbfgs_settings["grad_step_abs"]),
        )
        eval_id = self.eval_counter
        batch_id = f"lbfgs-g-{eval_id:06d}-{uuid.uuid4().hex[:8]}"
        self.eval_counter += 1
        requests = []
        point_keys = []
        for stencil in stencils:
            for point in stencil.points:
                shifted = np.array(x, copy=True)
                shifted[stencil.index] = point
                config = self._config_from_values(shifted)
                requests.append(
                    {
                        "config": config,
                        "label": "gradient",
                        "seed": int(self.args.rngseed),
                        "trial_id": f"{batch_id}-p{stencil.index}-{len(requests)}",
                    }
                )
                point_keys.append((stencil.index, point))

        logger.info(
            ".icetune_lbfgs: gradient eval=%06d launching probes=%d ndim=%d max_concurrent=%d",
            eval_id,
            len(requests),
            len(self.names),
            max(1, int(getattr(self.args, "max_concurrent_trials", 1))),
        )
        records = self._evaluate_many(requests)
        values_by_point = {key: float(record["cost"]) for key, record in zip(point_keys, records, strict=False)}
        gradient = finite_difference_gradient_from_values(stencils, values_by_point)
        costs = np.array([float(record["cost"]) for record in records], dtype=float)
        cached = sum(1 for record in records if record.get("cached"))
        logger.info(
            ".icetune_lbfgs: gradient eval=%06d done probes=%d cached=%d cost[min,mean,max]=[%.6g, %.6g, %.6g] |grad|=%.6g best=%s",
            eval_id,
            len(records),
            cached,
            float(np.min(costs)),
            float(np.mean(costs)),
            float(np.max(costs)),
            float(np.linalg.norm(gradient)),
            self._best_cost_text(),
        )
        return gradient


# Run bounded torch L-BFGS with Ray-distributed finite-difference gradients
def run_lbfgs_backend(*, args, tunesetup, simdriver, initial_points, param, start_time, storage_path=None):
    bounds = normalize_param_space(tunesetup.param_space)
    names = sorted(bounds)
    start = initial_config(bounds, initial_points, bool(args.no_initial_point))
    output_dir = pathlib.Path(storage_path or pathlib.Path(param["cdir"]) / "runs/icetune") / param["run_name"]
    param = copy.deepcopy(param)
    param["pickle_dir"] = param.get("pickle_dir") or str(output_dir / "results")
    objective = DistributedLBFGSObjective(args=args, bounds=bounds, names=names, param=param,
                                         simdriver=simdriver, output_dir=str(output_dir))
    lower = torch.tensor([bounds[name]["lower"] for name in names], dtype=torch.float64)
    width = torch.tensor([bounds[name]["upper"] - bounds[name]["lower"] for name in names], dtype=torch.float64)
    unit = (lower.new_tensor([start[name] for name in names]) - lower) / width
    settings = {key: value for key, value in args.lbfgs_settings.items() if not key.startswith("grad_step_")}

    # Defer distributed gradients until a line search accepts the objective value
    def evaluate(point):
        values = (lower + point * width).numpy()
        value = point.new_tensor(objective.value(*values))
        return value, lambda: point.new_tensor(objective.gradient(*values)) * width

    try:
        _, _, _, diagnostics = lbfgsb.minimize(None, unit, evaluator=evaluate, **settings)
    except EvaluationLimit:
        diagnostics = {"success": False, "message": "trial evaluation budget reached"}
    best = objective.best_record
    if best is None:
        raise RuntimeError("L-BFGS completed without a successful finite trial")
    config, metrics = best["config"], best["metrics"]
    published = _publish_best(args=args, config=config, param=param, simdriver=simdriver, metrics=metrics) if param.get("plot") else None
    history_plot = objective.publish_history()
    summary = {"algorithm": "lbfgs", "best_config": config, "best_metrics": metrics, "bounds": bounds,
        "best_fit": fit_summary.build_best_fit(values=config, uncertainties=None, objective_name=param["cost"],
            objective_value=metrics[param["cost"]], source="realized"),
        "diagnostics": diagnostics, "history_path": objective.history_path, "history_plot_path": history_plot,
        "published_summary": published, "num_unique_points": len(objective.cache), "start_config": start,
        **time_fields("completed_at", time.time()), "elapsed_seconds": time.time() - start_time}
    write_json_file(output_dir / "summary.json", json_safe(summary), indent=2, sort_keys=True, newline=True)
    return summary

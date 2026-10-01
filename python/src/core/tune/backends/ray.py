# Ray-distributed backend for icetune
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import copy
import hashlib
import importlib.util
import io
import json
import multiprocessing
import os
import pathlib
import pickle
import re
import shutil
import socket
import sys
import tarfile
import tempfile
import time
import uuid
from concurrent.futures import Future, ProcessPoolExecutor, ThreadPoolExecutor
from functools import wraps

import numpy as np
import ray
from ray import tune
from ray.tune.search import Searcher

from core.analysis import obs
from core.io import logger as log
from core.io.files import check_dependencies, ensure_dir, is_retryable_io_error
from core.io.serialize import finite_or_none, json_safe, selected_fields
from core.tune import core as icetune_main
from core.tune import events as event_output
from core.tune import history as icetune_history
from core.tune import short_id
from core.tune.cache import json_fingerprint, publish_file_immutable, sha256_file
from core.tune.io import (
    atomic_write_json,
    proposal_batch_size,
    proposal_slots,
    read_json,
    time_fields,
    trial_times,
    write_once_json,
)
from core.tune.parameters.space import config_contains, normalize_param_space, sample_config, typed_config
from core.tune.runtime import process as iceruntime
from core.tune.runtime import ray as ray_runtime

logger = log.get_logger(__name__)

RAY_TUNE_PENDING_ENV = "TUNE_MAX_PENDING_TRIALS_PG"


# Route Ray memory kills through the native actor cleanup and Tune retry callbacks
def install_ray_oom_handler() -> None:
    from ray.air.execution._internal.actor_manager import RayActorManager
    from ray.exceptions import OutOfMemoryError, RayActorError, RayTaskError

    original = RayActorManager._actor_task_failed
    if getattr(original, "icetune_oom", False):
        return

    # Keep this handler stateless so concurrent actor managers share no trial state
    @wraps(original)
    def actor_failed(self, tracked_actor_task, exception):
        if isinstance(exception, OutOfMemoryError) and not isinstance(exception, RayTaskError):
            error = RayActorError(str(exception))
            error.__cause__ = exception
            exception = error
        return original(self, tracked_actor_task, exception)

    actor_failed.icetune_oom = True
    RayActorManager._actor_task_failed = actor_failed


# Identify infrastructure exceptions without retrying invalid inputs or simulator errors
class RayClusterError(RuntimeError):
    """No usable trial allocation remains after the worker recovery interval."""


# Identify infrastructure exceptions without retrying invalid inputs or simulator errors
def ray_infrastructure_error(exc: BaseException) -> bool:
    from ray.exceptions import ActorUnschedulableError, OutOfMemoryError, RayActorError, RaySystemError

    seen = set()
    while exc is not None and id(exc) not in seen:
        seen.add(id(exc))
        if (isinstance(exc, OSError) and not isinstance(exc, (FileNotFoundError, PermissionError))
                and is_retryable_io_error(exc)):
            return True
        if isinstance(exc, RayActorError) and getattr(exc, "_actor_init_failed", False):
            return False
        if isinstance(exc, (RayClusterError, ActorUnschedulableError, OutOfMemoryError, RayActorError, RaySystemError)):
            return True
        exc = exc.__cause__
    return False


# Restore unfinished trials after bounded infrastructure failures in Tune itself
def fit_ray_tuner(*, tuner, restore, limit: int = 2, delay_s: float = 10):
    for attempt in range(limit + 1):
        try:
            return tuner.fit()
        except Exception as exc:
            if attempt >= limit or not ray_infrastructure_error(exc):
                raise
            logger.exception(".icetune_ray: Restoring interrupted campaign, attempt %d/%d", attempt + 1, limit)
            time.sleep(delay_s * (attempt + 1))
            tuner = restore()


# Start a fresh Python actor in the trial runtime without evaluating a physics point
class RayProbe:
    # Check that the uploaded runtime and its imports are available on this worker
    def ready(self, upload: bool, dependencies: dict | None = None) -> str:
        check_dependencies(dependencies or {})
        if upload:
            ray_runtime.ray_worker_root()
        return ray.get_runtime_context().get_node_id()


# Track allocations by Condor marker and Ray node ID, independently of physical hosts
class RayWorkers:
    # Keep health checks and retry counters local to this campaign
    def __init__(self, upload: bool, timeout_s: float = 600, poll_s: float = 60,
                 failure_limit: int = 2, recovery_s: float = 600, dependencies: dict | None = None):
        self.dependencies = dependencies or {}
        self.upload = upload
        self.timeout_s, self.poll_s = timeout_s, poll_s
        self.failure_limit, self.recovery_s = failure_limit, recovery_s
        self.node_ids, self.pending, self.failed = {}, {}, {}
        self.checked_at = 0.0
        self.empty_since = None

    # Record infrastructure failures with enough node information for allocation replacement
    def failure(self, error: str) -> None:
        for marker, node_id in self.node_ids.items():
            if node_id in error:
                previous = self.failed.get(marker, {})
                self.failed[marker] = {"node_id": node_id, "error": error[-16000:],
                    "count": int(previous.get("count", 0)) + 1}

    # Probe new nodes and observe lost nodes without blocking the Tune control thread
    def poll(self) -> dict:
        if not ray.is_initialized():
            return self.failed
        now = time.monotonic()
        for marker, (actor, future, started) in tuple(self.pending.items()):
            ready, _ = ray.wait([future], timeout=0)
            if not ready and now - started < self.timeout_s:
                continue
            self.pending.pop(marker)
            try:
                if not ready:
                    raise TimeoutError("Ray worker admission timed out")
                ray.get(future)
                self.failed.pop(marker, None)
            except Exception as exc:
                self.failure(f"{self.node_ids[marker]}: {exc}")
            finally:
                ray.kill(actor, no_restart=True)
        if now - self.checked_at < self.poll_s:
            return self.failed
        self.checked_at = now
        from ray.util.scheduling_strategies import NodeAffinitySchedulingStrategy

        live = {str(marker): str(node["NodeID"]) for node in ray.nodes() if node.get("Alive")
            for marker in node.get("Resources", {}) if marker.startswith("icetune_worker_")}
        for marker, node_id in self.node_ids.items():
            if marker not in live:
                self.failure(f"{node_id}: Ray worker disappeared")
        for marker, node_id in live.items():
            changed = self.node_ids.get(marker) != node_id
            if changed:
                self.failed.pop(marker, None)
                old = self.pending.pop(marker, None)
                if old is not None:
                    ray.kill(old[0], no_restart=True)
            self.node_ids[marker] = node_id
            if marker in self.pending or (not changed and marker not in self.failed):
                continue
            actor = (ray.remote(num_cpus=0, max_restarts=0)(RayProbe)
                .options(scheduling_strategy=NodeAffinitySchedulingStrategy(node_id=node_id, soft=False)) .remote())
            self.pending[marker] = (actor, actor.ready.remote(self.upload, self.dependencies), now)
        return self.failed

    # Stop bounded recovery if every known allocation remains unusable past its timeout
    def require_capacity(self) -> None:
        if not self.node_ids or any(
            self.failed.get(marker, {}).get("count", 0) < self.failure_limit for marker in self.node_ids):
            self.empty_since = None
            return
        now = time.monotonic()
        if self.empty_since is None:
            self.empty_since = now
        if now - self.empty_since >= self.recovery_s:
            raise RayClusterError(f"No usable Ray worker allocations after {self.recovery_s:g} seconds of allocation recovery")

    # Release admission actors when the campaign completes or its driver fails
    def close(self) -> None:
        for actor, _, _ in self.pending.values():
            try:
                ray.kill(actor, no_restart=True)
            except Exception as exc:
                logger.warning(".icetune_ray: Worker probe cleanup deferred: %s", exc)
        self.pending.clear()


# Fit one immutable optimizer snapshot outside the Ray Tune control thread
def _fit_search_state(request: dict) -> dict:
    from core.tune.search import SearchState

    arguments = {name: request[name] for name in ('args', 'bounds', 'initial_points', 'parameter_topology')}
    result = {'finished': False, 'diagnostics': []}
    if request['args'].algorithm not in {'hebo', 'icebo', 'hyperopt', 'optuna', 'ampfit'}:
        search = PercentileSearch(**arguments)
        search.completed_records = request['completed_records']
        search.failed_records = request['failed_records']
        search.live_records = {f'pending-{i}': {'config': config, 'trial_id': f'pending-{i}'}
                               for i, config in enumerate(request['pending_configs'])}
        search.param_space, search.next_index = request['param_space'], request['next_index']
        return {**result, 'proposals': search._ask_native_batch(request['count'])}
    search = SearchState(**arguments, async_proposals=True)
    search.observe(request['completed_records'])
    search.observe_failures(request['failed_records'])
    search.restore_issued(request['issued_records'])
    if search.icebo is not None and request.get('icebo_state') is not None:
        search.icebo.load_state_dict(request['icebo_state'])
    search.set_pending(request['pending_configs'])
    result['proposals'] = search.ask_many(request['next_index'], request['count'])
    if search.amplitude is not None:
        result.update(finished=search.amplitude.finished, diagnostics=search.amplitude.diagnostics)
    if search.icebo is not None:
        result.update(diagnostics=search.icebo.diagnostics(), icebo_state=search.icebo.state_dict())
    return result


# Format the complete retained simulator diagnostics for one Ray failure
def ray_failure_text(*, exc: BaseException, outputs: dict | None) -> str:
    lines = [str(exc)]
    if not isinstance(outputs, dict):
        return lines[0]
    lines.extend(f"{key}: {outputs[key]}" for key in ("error_type", "stage", "error") if outputs.get(key) is not None)
    failure = outputs.get("failure")
    if isinstance(failure, dict):
        lines.extend(("failure:", json.dumps(json_safe(failure), indent=2, sort_keys=True)))
    traceback_text = outputs.get("traceback")
    if traceback_text:
        lines.extend(("driver traceback:", str(traceback_text).rstrip()))
    return "\n".join(lines)


# Read the complete Ray traceback available to a driver callback
def ray_error_text(trial) -> str:
    getter = getattr(trial, "get_error", None)
    if callable(getter):
        try:
            error = getter()
            if error:
                return str(error).rstrip()
        except Exception as exc:
            logger.warning(".icetune_ray: Ray stored error read deferred: %s", exc)
    error_file = getattr(trial, "error_file", None)
    if error_file:
        try:
            return pathlib.Path(error_file).read_text(encoding="utf-8").rstrip()
        except OSError as exc:
            logger.warning(".icetune_ray: Ray local error read deferred: %s", exc)
    return "Ray trial failed without a persisted traceback"


# Compute one JSON safe snapshot of current live Ray node resources
def ray_resource_state() -> dict:
    if not ray.is_initialized():
        return {"available": {}, "cluster": {}, "nodes": {}}
    nodes = {str(node["NodeID"]): {"node_ip": str(node.get("NodeManagerAddress") or ""),
             "resources": {str(key): float(value) for key, value in (node.get("Resources") or {}).items()}}
             for node in ray.nodes() if node.get("Alive") and node.get("NodeID")}
    return {"available": {str(key): float(value) for key, value in ray.available_resources().items()},
        "cluster": {str(key): float(value) for key, value in ray.cluster_resources().items()}, "nodes": nodes}


# Publish live Ray state and rate limit automatic history plots
class HistoryPlotCallback(tune.Callback):
    # Initialize one completion driven plot timer
    def __init__(self, *, experiment_dir: str, cdir: str, cost: str, run_name: str, interval_s: float, search_alg,
        param: dict, num_trials: int, max_concurrent_trials: int,):
        self.experiment_dir, self.cdir = pathlib.Path(experiment_dir), pathlib.Path(cdir)
        self.cost, self.run_name = str(cost), str(run_name)
        self.interval_s = max(1.0, float(interval_s))
        self.search_alg = search_alg
        self.param = copy.deepcopy(param)
        self.num_trials, self.max_concurrent_trials = int(num_trials), int(max_concurrent_trials)
        self.history_path = self.experiment_dir / "history.json"
        self.status_path = self.experiment_dir / "status.json"
        self.failure_dir = self.experiment_dir / "failures"
        self.last_plot_at = self.last_history_at = self.last_resource_at = self.resource_observed_at = self.last_status_at = 0.0
        self.last_history = None
        self.plot_pool = self.plot_future = self.plot_history = None
        self.resource_interval_s = float(param.get("ray_resource_status_interval_s", 5))
        self.status_interval_s = float(param.get("ray_status_interval_s", 60))
        self.resource_state = {"available": {}, "cluster": {}, "nodes": {}}
        self.last_status = ()
        self.trial_status = {}
        self.workers = RayWorkers(
            bool(param.get("ray_upload_runtime")), timeout_s=param.get("ray_worker_admission_timeout_s", 600),
            poll_s=param.get("ray_worker_poll_interval_s", 60), failure_limit=param.get("ray_worker_failure_limit", 2),
            recovery_s=param.get("ray_worker_recovery_timeout_s", 600),
            dependencies=param.get("aux_param_space", {}).get("runtime_dependencies", {}))

    # Write one full failure record outside compact status and history
    def _write_failure(self, trial, *, retry: bool = False) -> pathlib.Path:
        ray_trial_id = str(getattr(trial, "trial_id", "unknown"))
        record = next((item for item in reversed(self.search_alg.failed_records)
                if str(item.get("ray_trial_id")) == ray_trial_id), {})
        trial_id = str(record.get("trial_id") or ray_trial_id)
        filename = re.sub(r"[^A-Za-z0-9_.-]+", "-", trial_id).strip("-") or "unknown"
        error = (getattr(trial, "last_result", None) or {}).get("trial_error") or ray_error_text(trial)
        payload = {"failure_schema_version": 1,
            "config": copy.deepcopy(record.get("config", getattr(trial, "config", None))), "error": error,
            "num_failures": int(getattr(trial, "num_failures", 1)), "ray_trial_id": ray_trial_id,
            "search_payload": copy.deepcopy(record.get("search_payload")), "trial_id": trial_id,
            **time_fields("failed_at", time.time())}
        suffix = f"-retry-{payload['num_failures']}" if retry else ""
        path = self.failure_dir / f"{filename}{suffix}.json"
        if path.is_file():
            return path
        atomic_write_json(str(path), json_safe(payload))
        logger.error('.icetune_ray: Trial failure written to "%s"\n%s', path, error)
        return path

    # Keep only the Ray trial identifiers and scheduler states
    def _capture(self, trials) -> None:
        self.trial_status = {str(trial.trial_id): str(trial.status).upper() for trial in trials
            if getattr(trial, "trial_id", None) is not None}

    # Write compact status and schema-1 likelihood history from canonical feedback
    def _write_state(self, *, force_history: bool = False) -> bool:
        now = time.time()
        done = list(self.search_alg.completed_records)
        failed = list(self.search_alg.failed_records)
        live = dict(self.search_alg.live_records)
        proposal_count = len(self.search_alg.proposal_queue)
        running = sum(self.trial_status.get(trial_id) == "RUNNING" for trial_id in live)
        queued = len(live) - running + proposal_count
        diagnostics = copy.deepcopy(getattr(self.search_alg, "optimizer_diagnostics", []))
        optimizer_finished = bool(getattr(self.search_alg, "optimizer_finished", False))
        history_signature = (tuple(record["trial_id"] for record in done), tuple(record["trial_id"] for record in failed),
                             json.dumps(diagnostics, sort_keys=True), optimizer_finished)
        history_due = history_signature != self.last_history and (
            force_history or self.last_history is None or now - self.last_history_at >= self.interval_s)
        if history_due:
            history = {**{key: self.param[key] for key in ("optimization", "parameter_space", "parameter_topology")},
                       **selected_fields(self.param, ("plot_brand", "simdriver"), include_missing=True),
                       "history_schema_version": icetune_history.HISTORY_SCHEMA_VERSION, "likelihood_default": "objective",
                       "optimizer_diagnostics": diagnostics, "mc_steer": self.param.get("mc_steer", {}),
                       "failed": failed, "num_trials": self.num_trials, "trials": done, **time_fields("updated_at", now)}
            atomic_write_json(str(self.history_path), json_safe(history))
            self.last_history = history_signature
            self.last_history_at = now

        best = min(done, key=lambda record: float(record["metrics"][self.cost]), default=None)
        finished = (self.search_alg.next_index >= self.num_trials or optimizer_finished) and not live
        if self.last_resource_at == 0.0 or now - self.last_resource_at >= self.resource_interval_s:
            self.last_resource_at = now
            try:
                self.resource_state = ray_resource_state()
                self.resource_state["unhealthy_workers"] = copy.deepcopy(self.workers.poll())
                self.resource_observed_at = now
            except Exception as exc:
                logger.warning(".icetune_ray: Ray resource status deferred: %s", exc)
        resources = self.resource_state
        status = {"backend": "ray", "best_result": None if best is None else selected_fields(best, (
                    "completed_at_datetime", "completed_at_unix", "config", "cost_definition", "metrics",
                    "search_payload", "theta_hash", "trial_id", "tunename"), include_missing=True),
            "completed": len(done), "done": finished,
            "optimizer_finished": optimizer_finished, "optimizer_diagnostics": diagnostics,
            "failed": len(failed), "max_concurrent_trials": self.max_concurrent_trials, "num_trials": self.num_trials,
            "pending": len(live) + proposal_count, "phase": "done" if finished else "running", "queued": queued,
            "ray_resources": resources, "ray_resources_observed_at_unix": self.resource_observed_at,
            "remaining_trials": max(0, self.num_trials - self.search_alg.next_index), "running": running,
            **time_fields("updated_at", now)}
        status_signature = (history_signature, len(failed), tuple(sorted(live)), queued, finished, tuple(sorted(
                    (key, value) for key, value in resources["cluster"].items()
                    if key == "icetune_head" or key.startswith("icetune_worker_"))))
        if status_signature != self.last_status or now - self.last_status_at >= self.status_interval_s:
            atomic_write_json(str(self.status_path), json_safe(status))
            self.last_status = status_signature
            self.last_status_at = now
        return bool(done) and history_due

    # Collect a finished plot without waiting in the Tune control thread
    def _collect_plot(self, *, wait: bool = False) -> pathlib.Path | None:
        if self.plot_future is None or (not wait and not self.plot_future.done()):
            return None
        try:
            path = self.plot_future.result()
        except Exception as exc:
            self.plot_history = None
            self.plot_pool.shutdown(wait=False, cancel_futures=True)
            self.plot_pool = None
            logger.warning(".icetune_ray: automatic history plot deferred: %s", exc)
            return None
        finally:
            self.plot_future = None
            self.last_plot_at = time.time()
        logger.info('.icetune_ray: History plot updated: "%s"', path)
        return path

    # Keep one head process rendering the latest history while Tune accepts results
    def plot(self, *, force: bool = False) -> pathlib.Path | None:
        self._write_state(force_history=force)
        path = self._collect_plot(wait=force)
        if self.plot_future is not None or not self.search_alg.completed_records:
            return path
        if not force and (self.last_history == self.plot_history or time.time() - self.last_plot_at < self.interval_s):
            return path
        try:
            if self.plot_pool is None:
                self.plot_pool = ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context("spawn"))
            self.plot_future = self.plot_pool.submit(icetune_history.publish_history_plot,
                history=read_json(self.history_path), source=self.history_path.resolve(), cdir=self.cdir.resolve(),
                cost=self.cost, run_name=self.run_name)
            self.plot_history = self.last_history
        except Exception as exc:
            self.last_plot_at = time.time()
            if self.plot_pool is not None:
                self.plot_pool.shutdown(wait=False, cancel_futures=True)
                self.plot_pool = None
            logger.warning(".icetune_ray: automatic history plot deferred: %s", exc)
        return self._collect_plot(wait=True) if force else path

    # Preserve final accepted feedback and release the head plotting process
    def close(self) -> None:
        try:
            self._write_state(force_history=True)
            self._collect_plot(wait=True)
        finally:
            if self.plot_pool is not None:
                self.plot_pool.shutdown(wait=True)
                self.plot_pool = None

    # Exclude active process handles from Ray callback checkpoints
    def __getstate__(self) -> dict:
        return dict(self.__dict__, plot_pool=None, plot_future=None, plot_history=None)

    # Preserve nested trial results before Ray flattens the completion feedback
    def on_trial_result(self, *, iteration, trials, trial, result, **info) -> None:
        live = self.search_alg.live_records.get(trial.trial_id)
        if live is not None:
            live["result"] = copy.deepcopy(result)

    # Refresh after trial completion using the configured rate limit
    def on_trial_complete(self, *, iteration, trials, trial, **info) -> None:
        self._capture(trials)
        if (getattr(trial, "last_result", None) or {}).get("trial_failed"):
            self._write_failure(trial)
        self.plot()

    # Publish new Ray trial states immediately
    def on_trial_start(self, *, iteration, trials, trial, **info) -> None:
        self._capture(trials)
        live = self.search_alg.live_records.get(trial.trial_id)
        if live is not None:
            started = finite_or_none(getattr(trial, "start_time", None)) or time.time()
            live.update(time_fields("started_at", started))
        self._write_state()

    # Refresh live counts while Ray waits for trial outcomes
    def on_step_end(self, *, iteration, trials, **info) -> None:
        self._capture(trials)
        self.plot()
        if self.search_alg.live_records:
            self.workers.require_capacity()

    # Publish failed trial counts immediately
    def on_trial_error(self, *, iteration, trials, trial, retry=False, **info) -> None:
        self._capture(trials)
        self._write_failure(trial, retry=retry)
        self._worker_failure(trial)
        self._write_state()

    # Preserve diagnostics for retries without adding failed points to optimizer feedback
    def on_trial_recover(self, *, iteration, trials, trial, **info) -> None:
        self.on_trial_error(iteration=iteration, trials=trials, trial=trial, retry=True, **info)

    # Request replacement only for infrastructure failures that identify a Ray node
    def _worker_failure(self, trial) -> None:
        error = (getattr(trial, "last_result", None) or {}).get("trial_error") or ray_error_text(trial)
        if any(text in error for text in ("worker startup repeatedly failed",
                "icetune infrastructure failure on Ray node", "OutOfMemoryError", "node running low on memory",
                "node has died")):
            self.workers.failure(error)
            self.last_resource_at = 0.0


# Copy one worker-produced figure tree through the canonical directory swap
def _publish_ray_figure_tree(*, source: pathlib.Path, target: pathlib.Path, param: dict) -> None:
    token = icetune_main.acquire_publish_lock(cdir=param["cdir"], run_name=param["run_name"])
    stage = icetune_main.default_tmp_stage_path(target)
    try:
        shutil.copytree(source, stage, dirs_exist_ok=True)
        names = ("covariance", icetune_main.ICETUNE_INITIAL_FIGURE_SUBDIR, *icetune_main.ICETUNE_HISTORY_FIGURES)
        for name in names:
            old = target / name
            new = stage / name
            if new.exists() or not old.exists():
                continue
            (shutil.copytree if old.is_dir() else shutil.copy2)(old, new)
        icetune_main._commit_figure_dir(target_dir=str(target), tmp_dir=str(stage))
    finally:
        if stage.exists():
            shutil.rmtree(stage)
        icetune_main.release_publish_lock(cdir=param["cdir"], run_name=param["run_name"], token=token)


# Validate figures and full pickles transferred to the dedicated Ray head
class TrialOutputCallback(tune.Callback):
    # Initialize output targets and recover only the selected cost identity
    def __init__(self, *, experiment_dir: str, param: dict, global_state, simdriver=None, runtime_param=None):
        self.experiment_dir = pathlib.Path(experiment_dir)
        self.param = copy.deepcopy(param)
        self.global_state = global_state
        self.figure_dir = pathlib.Path(icetune_main.icetune_figure_dir(cdir=param["cdir"], run_name=param["run_name"]))
        self.pending, self.submitted, self.restored = {}, set(), False
        self.simdriver = simdriver
        self.runtime_param = copy.deepcopy(param if runtime_param is None else runtime_param)
        self.render_records = {}
        self.output_jobs = {name: dict(record=None, future=None, publishing=False, failures={})
                            for name in ("render", "covariance")}
        self.covariance_worker = None
        self.covariance_next_at = 0.0
        self.recovery_refs = set()
        summary = read_json(self.figure_dir / "summary.json", default={}) or {}
        self.rendered_trial_id = summary.get("trial_id")
        covariance = read_json(self.figure_dir / "covariance" / "covariance.json", default={}) or {}
        self.covariance_trial_id = covariance.get("trial_id")
        self.rendered_cost = finite_or_none((summary.get("metrics") or {}).get(param["cost"]))
        self.initial_rendered = (self.figure_dir / icetune_main.ICETUNE_INITIAL_FIGURE_SUBDIR / "summary.json"
        ).is_file()

    # Submit one idempotent output operation without waiting for head storage
    def _submit(self, method: str, trial_id: str, token=None) -> None:
        key = (method, str(trial_id), token if token is not None else uuid.uuid4().hex)
        if key in self.submitted:
            return
        arguments = (trial_id,) if token is None else (trial_id, token)
        future = getattr(self.global_state, method).remote(*arguments)
        self.pending[future] = key
        self.submitted.add(key)

    # Collect completed transfers while keeping the Tune control thread free
    def flush(self, *, wait: bool = False) -> None:
        while self.pending:
            ready, _ = ray.wait(list(self.pending), num_returns=1, timeout=None if wait else 0)
            if not ready:
                break
            future = ready[0]
            key = self.pending.pop(future)
            self.recovery_refs.discard(future)
            try:
                receipt = ray.get(future)
                if isinstance(receipt, dict) and receipt.get("render_token"):
                    self.render_records[receipt["render_token"]] = receipt["render_record"]
            except Exception:
                self.submitted.discard(key)
                raise
        for name in self.output_jobs:
            self._poll_output(name, wait=wait)

    # Keep only the initial point and the best accepted render candidate
    def _next_render(self) -> dict | None:
        candidates = list(self.render_records.values())
        best = min(candidates, key=lambda record: record["cost"], default=None)
        initial = next((record for record in candidates if record["initial"] and not self.initial_rendered), None)
        rendered = next((record for record in candidates if record["trial_id"] == self.rendered_trial_id), None)
        active = [job["record"] for job in self.output_jobs.values()]
        keep = {record["token"] for record in (best, initial, rendered, *active) if record is not None}
        for token in list(self.render_records):
            if token not in keep:
                self.global_state.discard_render.remote(token)
                self.render_records.pop(token)
        selected = next((record for record in (initial, best) if record is not None
                         and self.output_jobs["render"]["failures"].get(record["token"], 0) < 3), None)
        if selected is None:
            return None
        improving = self.rendered_cost is None or selected["cost"] < self.rendered_cost
        if not improving and selected is not initial:
            return None
        return {**selected, "render_initial": selected is initial, "render_best": improving}

    # Select the latest rendered fit once the covariance interval has elapsed
    def _next_covariance(self, *, wait: bool) -> dict | None:
        if not wait and time.time() < self.covariance_next_at:
            return None
        return next((record for record in self.render_records.values() if record["trial_id"] == self.rendered_trial_id
                     and record["trial_id"] != self.covariance_trial_id
                     and self.output_jobs["covariance"]["failures"].get(record["trial_id"], 0) < 3), None)

    # Poll independent figure and covariance jobs through computation and publication
    def _poll_output(self, name: str, *, wait: bool = False) -> None:
        optimization = self.param.get("optimization", {})
        covariance = name == "covariance"
        if self.simdriver is None or (not covariance and self.recovery_refs):
            return
        if covariance and (optimization.get("optimizer") != "ampfit" or not optimization["ampfit"]["covariance"]):
            return
        job = self.output_jobs[name]
        while True:
            if job["future"] is None:
                record = self._next_covariance(wait=wait) if covariance else self._next_render()
                if record is None:
                    return
                payload = self.global_state.render_payload.remote(record["token"])
                options = {"num_cpus": 1}
                if not covariance and ray.cluster_resources().get("icetune_head", 0) > 0:
                    options = {"num_cpus": 0, "resources": {"icetune_head": 0.001}}
                if covariance:
                    if self.covariance_worker is None:
                        self.covariance_worker = CovarianceWorker.options(**options).remote(self.runtime_param, self.simdriver)
                    future = self.covariance_worker.compute.remote(payload)
                else:
                    future = render_completed_trial.options(**options).remote(
                        payload, self.runtime_param, self.simdriver, record["render_initial"], record["render_best"])
                job.update(record=record, future=future, publishing=False)
            ready, _ = ray.wait([job["future"]], timeout=None if wait else 0)
            if not ready:
                return
            record = job["record"]
            try:
                if not job["publishing"]:
                    future = getattr(self.global_state, f"publish_{name}").remote(record, job["future"])
                    job.update(future=future, publishing=True)
                    continue
                ray.get(job["future"])
                if covariance:
                    self.covariance_trial_id = record["trial_id"]
                else:
                    if record["render_best"]:
                        self.rendered_cost, self.rendered_trial_id = record["cost"], record["trial_id"]
                    self.initial_rendered |= record["render_initial"]
            except Exception as exc:
                key = record["trial_id" if covariance else "token"]
                job["failures"][key] = job["failures"].get(key, 0) + 1
                logger.warning(".icetune_ray: %s failed for %s: %s", name, record["trial_id"], exc)
                if covariance:
                    ray.kill(self.covariance_worker)
                    self.covariance_worker = None
            job.update(future=None, record=None)
            if covariance:
                self.covariance_next_at = time.time() + optimization["ampfit"]["covariance_interval"]
            elif not wait:
                return

    # Recover accepted output tokens before polling asynchronous publication
    def on_step_end(self, *, iteration, trials, **info) -> None:
        if not self.restored:
            for trial in trials:
                if trial.status == "TERMINATED":
                    metrics = trial.last_result or {}
                    if metrics.get("trial_failed"):
                        self.on_trial_error(iteration=iteration, trials=trials, trial=trial)
                    elif metrics.get("ray_output_token"):
                        future = self.global_state.recover_trial_outputs.remote(metrics)
                        self.pending[future] = ("recover_trial_outputs", trial.trial_id, metrics["ray_output_token"])
                elif trial.status == "ERROR":
                    self.on_trial_error(iteration=iteration, trials=trials, trial=trial)
            self.restored = True
            self.recovery_refs = set(self.pending)
        self.flush()

    # Reconstruct pending publication from terminal Tune results after restore
    def __getstate__(self) -> dict:
        jobs = {name: dict(job, future=None, record=None, publishing=False) for name, job in self.output_jobs.items()}
        return dict(self.__dict__, pending={}, submitted=set(), restored=False, output_jobs=jobs,
                    covariance_worker=None, recovery_refs=set())

    # Commit head-staged outputs only after Tune accepts a successful result
    def on_trial_complete(self, *, iteration, trials, trial, **info) -> None:
        metrics = getattr(trial, "last_result", None) or {}
        if metrics.get("trial_failed"):
            self.on_trial_error(iteration=0, trials=[], trial=trial)
            return
        output_token = metrics.get("ray_output_token")
        if output_token is not None:
            self._submit("commit_trial_outputs", metrics.get("trial_id"), output_token)

    # Release staged outputs after a failed Ray trial
    def on_trial_error(self, *, iteration, trials, trial, **info) -> None:
        metrics = getattr(trial, "last_result", None) or {}
        trial_id = metrics.get("trial_id") or getattr(trial, "trial_id", None)
        if trial_id is not None:
            self._submit("discard_trial_outputs", trial_id)

    # Release failed outputs before Tune repeats the same trial
    def on_trial_recover(self, **info) -> None:
        self.on_trial_error(**info)

    # Validate required canonical outputs against Ray's exact final best result
    def validate(self, *, best_result, require_initial: bool) -> dict:
        self.flush(wait=True)
        metrics = copy.deepcopy(best_result.metrics or {})
        trial_id = metrics.get("trial_id")
        if self.param.get("plot", False):
            summary = read_json(self.figure_dir / "summary.json", default={}) or {}
            if summary.get("trial_id") != trial_id:
                raise RuntimeError("Ray final best figure output is incomplete")
            initial = self.figure_dir / icetune_main.ICETUNE_INITIAL_FIGURE_SUBDIR / "summary.json"
            if require_initial and not initial.is_file():
                raise RuntimeError("Ray initial figure output is incomplete")
        descriptor = metrics.get("trial_pickle")
        if self.param.get("pickle_dump", False):
            if not isinstance(descriptor, dict):
                raise RuntimeError("Ray final best full trial pickle descriptor is missing")
            ray_trial_pickle_path(trial_dir=self.experiment_dir, descriptor=descriptor)
        return {"best_trial_id": trial_id, "figure_dir": str(self.figure_dir)}


# Validate and return one trial pickle relative path
def ray_trial_pickle_relative(descriptor: dict) -> pathlib.PurePosixPath:
    relative = pathlib.PurePosixPath(str(descriptor.get("relative_path") or ""))
    if (relative.is_absolute() or ".." in relative.parts or len(relative.parts) != 2 or relative.parts[0] != "results"
        or not re.fullmatch(r"TUNE_icetune_[^/]+\.pkl", relative.name)
        or str(descriptor.get("filename") or "") != relative.name):
        raise ValueError("Ray trial pickle has an invalid relative path")
    return relative


# Validate and return one trial pickle path below its worker-local directory
def ray_trial_pickle_path(*, trial_dir: pathlib.Path, descriptor: dict) -> pathlib.Path:
    relative = ray_trial_pickle_relative(descriptor)
    source = trial_dir.joinpath(*relative.parts)
    if (not source.is_file() or source.stat().st_size != int(descriptor["size"])
        or sha256_file(source) != descriptor["sha256"]):
        raise RuntimeError(f'Ray trial pickle is incomplete: "{source}"')
    return source


# Package one worker-local figure tree for transfer through the Ray object store
def pack_ray_figure_tree(source: pathlib.Path) -> bytes:
    if not source.is_dir():
        raise RuntimeError(f'Ray trial figure tree is missing: "{source}"')
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w:gz") as archive:
        archive.add(source, arcname="figures", recursive=True)
    return payload.getvalue()


# Collect only durable selected outputs from one worker-local Ray trial
def collect_ray_trial_outputs(
    *, trial_dir: pathlib.Path, run_name: str, descriptor: dict | None, rendered_initial: bool, rendered_best: bool
) -> dict[str, bytes | None]:
    pickle_payload = (ray_trial_pickle_path(trial_dir=trial_dir, descriptor=descriptor).read_bytes()
                      if descriptor is not None else None)

    figure_payload = None
    if rendered_initial or rendered_best:
        source = trial_dir / "figs" / "icetune" / run_name
        for enabled, subdir, label in ((rendered_initial, icetune_main.ICETUNE_INITIAL_FIGURE_SUBDIR, "initial"),
                                      (rendered_best, "", "best")):
            if enabled and not (source / subdir / "summary.json").is_file():
                raise RuntimeError(f'Ray {label} figure tree is incomplete: "{source / subdir}"')
        figure_payload = pack_ray_figure_tree(source)
    return {"figures": figure_payload, "pickle": pickle_payload}


# Extract one trusted worker figure payload into a new head-local stage directory
def extract_ray_figure_tree(*, payload: bytes, stage: pathlib.Path) -> pathlib.Path:
    if not isinstance(payload, bytes) or not payload:
        raise ValueError("Ray figure transfer payload must be non-empty bytes")
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
        members = archive.getmembers()
        for member in members:
            relative = pathlib.PurePosixPath(member.name)
            if (not relative.parts or relative.parts[0] != "figures" or relative.is_absolute() or ".." in relative.parts
                or not (member.isdir() or member.isfile())):
                raise ValueError(f"Invalid Ray figure transfer path: {member.name}")
        archive.extractall(stage, members=members, filter="data")
    source = stage / "figures"
    if not source.is_dir():
        raise RuntimeError("Ray figure transfer has no figure directory")
    return source


# Compute whether optimizer metadata marks one trial as the initial point
def _search_payload_is_initial(payload: dict | None) -> bool:
    return isinstance(payload, dict) and str(payload.get("kind")) == "initial"


# Build Ray initialization arguments for a local runtime or existing cluster
def ray_init_arguments(*, args, runtime_env) -> dict:
    address = (args.address or "local").strip()
    if address.lower() != "local":
        return {"address": address, "runtime_env": runtime_env}

    # Keep a local fit independent of workstation network and VPN changes
    options = {"include_dashboard": False, "_node_ip_address": "127.0.0.1", "num_gpus": int(getattr(args, "ray_worker_gpu", 0))}
    if getattr(args, "ray_worker_cpu", None) is not None:
        options["num_cpus"] = int(args.ray_worker_cpu)
    if getattr(args, "ray_temp_dir", None):
        ensure_dir(pathlib.Path(args.ray_temp_dir))
        options["_temp_dir"] = str(pathlib.Path(args.ray_temp_dir).resolve())
    return options


# Match Ray Tune placement group staging to the effective trial concurrency
def configure_tune_pending_trials(args) -> int:
    target = min(max(1, int(args.max_concurrent_trials)), max(1, int(args.num_trials)))
    raw = os.environ.get(RAY_TUNE_PENDING_ENV)
    if raw is None or raw.strip().lower() == "auto":
        os.environ[RAY_TUNE_PENDING_ENV] = str(target)
        logger.info(".icetune_ray: Ray Tune pending placement group limit=%d from trial concurrency", target)
        return target
    try:
        configured = int(raw)
    except ValueError as exc:
        raise ValueError(f"{RAY_TUNE_PENDING_ENV} must be a positive integer or auto") from exc
    if configured <= 0:
        raise ValueError(f"{RAY_TUNE_PENDING_ENV} must be a positive integer or auto")
    logger.info(".icetune_ray: Ray Tune pending placement group limit=%d from environment", configured)
    return configured


# Replace checkout paths with paths relative to the uploaded Ray runtime
def portable_ray_environment(variables: dict[str, str], *, cdir: str) -> dict[str, str]:
    root = pathlib.Path(cdir).resolve()
    output = copy.deepcopy(variables)
    paths = []
    for value in output.get("PYTHONPATH", "").split(os.pathsep):
        if not value:
            continue
        try:
            relative = pathlib.Path(value).resolve().relative_to(root)
        except (OSError, ValueError):
            paths.append(value)
            continue
        paths.append(relative.as_posix())
    output["PYTHONPATH"] = os.pathsep.join(paths)
    output["RAY_CHDIR_TO_TRIAL_DIR"] = "0"
    return output


# Inspect one Ray node without consuming a scheduled CPU
def _ray_cluster_node_probe(*, cdir: str, libdir: str, storage_path: str) -> dict:
    runtime_root = pathlib.Path(cdir)
    return {"cdir": runtime_root.is_dir(), "entrypoint": importlib.util.find_spec("core.icetune") is not None,
        "hostname": socket.gethostname(), "libdir": pathlib.Path(libdir).is_dir(),
        "python_version": f"{sys.version_info.major}.{sys.version_info.minor}", "ray_version": ray.__version__,
        "storage": pathlib.Path(storage_path).is_dir(), "storage_writable": os.access(storage_path, os.R_OK | os.W_OK)}


# Validate each worker node runtime configuration
def validate_ray_node_contracts(results: list[dict], *, expected_python: str, expected_ray: str) -> None:
    problems = []
    for result in results:
        missing = [key for key in ("cdir", "entrypoint", "libdir", "storage", "storage_writable") if not result[key]]
        for key, expected in (("python_version", expected_python), ("ray_version", expected_ray)):
            if result[key] != expected:
                missing.append(f"{key}={result[key]} expected={expected}")
        if missing:
            problems.append(f"{result['hostname']}: {', '.join(missing)}")
    if problems:
        raise RuntimeError("Ray cluster nodes do not share the requested runtime contract: " + "; ".join(problems))


# Verify that every connected node sees the same runtime and Tune storage paths
def validate_ray_cluster_paths(*, cdir: str, libdir: str, storage_path: str) -> list[dict]:
    from ray.util.scheduling_strategies import NodeAffinitySchedulingStrategy

    live_nodes = [node for node in ray.nodes() if bool(node.get("Alive", False))]
    if not live_nodes:
        raise RuntimeError("Ray reports no live nodes after initialization")

    probe = ray.remote(num_cpus=0)(_ray_cluster_node_probe)
    references = [
        probe.options(scheduling_strategy=NodeAffinitySchedulingStrategy(node_id=node["NodeID"], soft=False)).remote(
            cdir=cdir, libdir=libdir, storage_path=storage_path) for node in live_nodes]
    results = ray.get(references)
    expected_python = f"{sys.version_info.major}.{sys.version_info.minor}"
    validate_ray_node_contracts(results, expected_python=expected_python, expected_ray=ray.__version__)
    return results


# Wait for worker allocations and a complete trial bundle within one allocation
def wait_ray_trial_capacity(args) -> dict[str, dict[str, float]]:
    required_cpus = float(args.cpu_per_trial)
    required_gpus = float(args.gpu_per_trial)
    wait_workers = getattr(args, "wait_workers", None)
    timeout = float(getattr(args, "wait_workers_timeout_s", 300.0))
    deadline = time.monotonic() + timeout
    while True:
        live = [node for node in ray.nodes() if bool(node.get("Alive", False))]
        node_resources = [
            {str(key): float(value) for key, value in dict(node.get("Resources") or {}).items()} for node in live]
        worker_resources = [
            item for item in node_resources if item.get("icetune_head", 0.0) <= 0.0 and item.get("CPU", 0.0) > 0.0]
        resources = ray_resource_state()
        cluster = resources["cluster"]
        workers_ready = wait_workers is None or len(worker_resources) >= int(wait_workers)
        bundle_ready = any(
            item.get("CPU", 0.0) >= required_cpus and (required_gpus <= 0.0 or item.get("GPU", 0.0) >= required_gpus)
            for item in worker_resources)
        if workers_ready and bundle_ready:
            return resources
        if wait_workers is None or time.monotonic() >= deadline:
            raise RuntimeError(f"Ray has no complete trial resource bundle after {timeout:.0f} seconds, "
                f"workers={len(worker_resources)}, required_cpu={required_cpus:g}, "
                f"required_gpu={required_gpus:g}, node_resources={node_resources}, " f"cluster_resources={cluster}")
        time.sleep(1.0)


class PercentileSearch(Searcher):
    """Ray Tune searcher backed by percentile-filtered optimizer state."""

    CHECKPOINT_SCHEMA_VERSION = 1

    # Initialize a Ray searcher that replays only the selected surrogate-fit records
    def __init__(self, *, args, bounds: dict, initial_points: dict | None, parameter_topology: dict | None = None):
        super().__init__(metric=getattr(args, "cost", None), mode="min")
        self.args = copy.deepcopy(args)
        self.bounds = copy.deepcopy(bounds)
        self.param_space = None
        self.initial_points = copy.deepcopy(initial_points)
        self.parameter_topology = copy.deepcopy(parameter_topology or {})
        self.completed_records, self.failed_records, self.live_records = [], [], {}
        self.next_index, self.initial_completed = 0, False
        self.optimizer_finished, self.optimizer_diagnostics = False, []
        self.proposal_queue, self.reserve = [], []
        self.fit_feedback = None
        self.fit_request = self.fit_future = self.fit_started_at = self.fit_executor = None
        self.icebo_state, self.icebo_waiting = None, False

    # Compute whether one Ray proposal belongs to the independent random regime
    def _independent_record(self, record: dict) -> bool:
        payload = record.get("search_payload")
        kind = payload.get("kind") if isinstance(payload, dict) else None
        return str(kind) in icetune_history.RANDOM_PROPOSAL_TYPES

    # Compute the independent proposal target before adaptive Ray fitting
    def _independent_target(self, total: int) -> int:
        if str(self.args.algorithm) == "basic":
            return 0
        requested = max(0, int(getattr(self.args, "rand_trials", 0)))
        return min(int(total), max(1, requested))

    # Compute a config coerced to normalized parameter types when bounds are known
    def _typed_config(self, config: dict) -> dict:
        return typed_config(config, self.bounds) if self.bounds else copy.deepcopy(config)

    # Compute completed records selected for the current surrogate fit
    def _selected_records(self) -> list[dict]:
        return icetune_main.select_surrogate_fit_records(
            self.completed_records, self.metric, getattr(self.args, "surrogate_fit_percentile", 0.95))

    # Compute one hashable optimizer configuration identity
    def _config_key(self, config: dict) -> tuple:
        typed = self._typed_config(config)
        return tuple((key, typed[key]) for key in sorted(self.bounds))

    # Compute every config already live or reserved by the Ray searcher
    def _pending_configs(self) -> list[dict]:
        configs = [record["config"] for record in self.live_records.values() if isinstance(record.get("config"), dict)]
        configs.extend(config for config, _ in self.proposal_queue)
        configs.extend(config for config, _ in self.reserve)
        return [copy.deepcopy(config) for config in configs if isinstance(config, dict)]

    # Copy only optimizer feedback into a fit snapshot
    def _fit_records(self, records: list[dict]) -> list[dict]:
        return [{**{key: record[key] for key in ('trial_id', 'config', 'metrics', 'search_payload', 'gradient')
                    if key in record},
                 'likelihood': {key: value for key, value in (record.get('likelihood') or {}).items()
                                if key in ('objective', 'mc_uncertainty')}} for record in records]

    # Build one immutable fit request for the background worker
    def _request_fit(self, *, index: int, count: int) -> dict:
        completed, failed = self._fit_records(self.completed_records), self._fit_records(self.failed_records)
        issued = sorted([*completed, *failed, *self._fit_records(list(self.live_records.values()))],
                        key=lambda record: str(record['trial_id']))
        return copy.deepcopy(dict(args=self.args, bounds=self.bounds, initial_points=self.initial_points,
                                  parameter_topology=self.parameter_topology, param_space=self.param_space,
                                  completed_records=completed, failed_records=failed, issued_records=issued,
                                  next_index=int(index), count=int(count), pending_configs=self._pending_configs(),
                                  icebo_state=self.icebo_state))

    # Submit one fit without waiting in Ray Tune's control thread
    def _submit_fit(self, request: dict) -> Future:
        if self.fit_executor is None:
            self.fit_executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="icetune-fit")
        return self.fit_executor.submit(_fit_search_state, request)

    # Start a refit with pending trials before the rolling reserve runs empty
    def _start_fit(self) -> None:
        if self.args.algorithm == "basic" or self.fit_future is not None or self.optimizer_finished:
            return
        feedback = len(self.completed_records) + len(self.failed_records)
        if self.args.algorithm == "icebo" and self.icebo_waiting and self.fit_feedback == feedback:
            return
        if self.fit_request is None:
            total = int(getattr(self.args, "num_trials", self.next_index + 1))
            reserved = len(self.proposal_queue) + len(self.reserve)
            if self.args.algorithm == "ampfit":
                if self.fit_feedback == feedback:
                    return
                count = int(self.args.max_concurrent_trials)
            else:
                if self.next_index < self._independent_target(total) or not self.completed_records:
                    return
                batch = proposal_batch_size(self.args)
                target = min(max(1, int(getattr(self.args, "max_concurrent_trials", 1))), 2 * batch)
                due = self.fit_feedback is None or feedback - int(self.fit_feedback) >= batch or reserved == 0
                if not due or reserved > max(0, target - batch):
                    return
                count = max(batch, target - reserved)
            count = min(count, total - self.next_index - reserved)
            if count <= 0:
                return
            self.fit_request = self._request_fit(index=self.next_index + reserved, count=count)
            self.fit_feedback = feedback
        request = self.fit_request
        selected = len(icetune_main.select_surrogate_fit_records(
                request["completed_records"], self.metric, getattr(self.args, "surrogate_fit_percentile", 0.95)))
        logger.info(
            ".icetune_ray: Starting background %s fit with completed=%d selected=%d failed=%d pending=%d requested=%d",
            self.args.algorithm, len(request["completed_records"]), selected, len(request["failed_records"]),
            len(request["pending_configs"]), request["count"])
        self.fit_started_at = time.monotonic()
        self.fit_future = self._submit_fit(request)

    # Move a completed background fit into the rolling proposal reserve
    def _collect_fit(self) -> None:
        future = self.fit_future
        if future is None or not future.done():
            return
        started = self.fit_started_at or time.monotonic()
        request = self.fit_request or {}
        self.fit_future = self.fit_request = self.fit_started_at = None
        try:
            proposals = future.result()
        except Exception:
            logger.exception(".icetune_ray: Background %s fit failed after %.1f seconds", self.args.algorithm,
                time.monotonic() - started)
            raise
        self.optimizer_finished = proposals['finished']
        self.optimizer_diagnostics = proposals['diagnostics']
        if 'icebo_state' in proposals:
            self.icebo_state = proposals['icebo_state']
            self.icebo_waiting = not proposals['proposals']
        proposals = proposals['proposals']
        if self.args.algorithm not in {"ampfit", "icebo"} and len(proposals) != int(request.get("count", len(proposals))):
            raise RuntimeError(
                f"{self.args.algorithm} returned {len(proposals)} of {request.get('count')} Ray proposals")
        total = int(getattr(self.args, "num_trials", self.next_index + len(proposals)))
        available = max(0, total - self.next_index - len(self.proposal_queue) - len(self.reserve))
        blocked = {self._config_key(config) for config in self._pending_configs() if isinstance(config, dict)}
        accepted = []
        for config, payload in proposals:
            if len(accepted) >= available:
                break
            typed = self._typed_config(config)
            key = self._config_key(typed)
            if key in blocked:
                continue
            accepted.append((typed, copy.deepcopy(payload)))
            blocked.add(key)
        self.reserve.extend(accepted)
        logger.info(
            ".icetune_ray: Completed background %s fit with proposals=%d accepted=%d reserve=%d in %.1f seconds",
            self.args.algorithm, len(proposals), len(accepted), len(self.reserve), time.monotonic() - started)

    # Compute a deterministic uniform proposal for random-start phases
    def _uniform_config(self, index: int) -> dict:
        rng = np.random.RandomState(int(getattr(self.args, "rngseed", 0)) + int(index))
        return sample_config(bounds=self.bounds, uniform=lambda lower, upper: float(rng.uniform(lower, upper)),
            integer=lambda lower, upper: int(rng.randint(lower, upper + 1)))

    # Ask BayesianOptimization for one pending aware proposal group
    def _ask_bayesopt_batch(self, count: int) -> list[tuple[dict, dict]]:
        selected = self._selected_records()
        proposals = []

        from bayes_opt import BayesianOptimization

        pbounds = {key: (float(spec["lower"]), float(spec["upper"])) for key, spec in sorted(self.bounds.items())}
        optimizer_kwargs = {"f": None, "pbounds": pbounds, "random_state": int(getattr(self.args, "rngseed", 0)),
            "verbose": 0}
        try:
            from bayes_opt.acquisition import ConstantLiar, UpperConfidenceBound

            acquisition = UpperConfidenceBound(kappa=2.576, random_state=self.args.rngseed)
            optimizer_kwargs["acquisition_function"] = (
                ConstantLiar(acquisition, strategy="max", random_state=self.args.rngseed) if count > 1 else acquisition)
            optimizer_kwargs["allow_duplicate_points"] = True
        except Exception:
            pass

        optimizer = BayesianOptimization(**optimizer_kwargs)
        for record in selected:
            optimizer.register(params={key: float(record["config"][key]) for key in sorted(self.bounds)},
                target=-float(record["metrics"][self.metric]))
        pending_configs = [record["config"] for record in self.live_records.values()]
        if selected and pending_configs:
            targets = np.array([-float(record["metrics"][self.metric]) for record in selected], dtype=float)
            worst = float(np.min(targets))
            penalty = max(float(np.std(targets)), abs(worst) * 1.0e-6, 1.0e-6)
            for config in pending_configs:
                optimizer.register(
                    params={key: float(config[key]) for key in sorted(self.bounds)}, target=worst - penalty)

        for _ in range(count):
            try:
                suggestion = optimizer.suggest()
            except TypeError:
                from bayes_opt import UtilityFunction

                suggestion = optimizer.suggest(UtilityFunction(kind="ucb", kappa=2.576, xi=0.0))
            proposals.append((self._typed_config({key: suggestion[key] for key in sorted(self.bounds)}),
                    {"kind": self.args.algorithm}))
        return proposals

    # Ask Ax directly with completed, failed and live trial state
    def _ask_ax_batch(self, count: int) -> list[tuple[dict, dict]]:
        from ax.service.ax_client import AxClient
        from ax.service.utils.instantiation import ObjectiveProperties
        from ray.tune.search.ax import AxSearch

        selected_trial_ids = {record["trial_id"] for record in self._selected_records()}
        ax_client = AxClient(verbose_logging=False, enforce_sequential_optimization=False, random_seed=self.args.rngseed)
        ax_client.create_experiment(parameters=AxSearch.convert_search_space(self.param_space),
            objectives={self.metric: ObjectiveProperties(minimize=self.mode != "max")},
            choose_generation_strategy_kwargs={"num_initialization_trials": 0}, is_test=True)
        for record in self.completed_records:
            _, trial_index = ax_client.attach_trial(record["config"])
            if record["trial_id"] in selected_trial_ids:
                ax_client.complete_trial(
                    trial_index=trial_index, raw_data={self.metric: float(record["metrics"][self.metric])})
            else:
                ax_client.abandon_trial(trial_index=trial_index, reason="excluded from percentile surrogate fit")
        for record in self.failed_records:
            _, trial_index = ax_client.attach_trial(record["config"])
            ax_client.abandon_trial(trial_index=trial_index, reason="Ray trial failed")
        for record in self.live_records.values():
            ax_client.attach_trial(record["config"])
        proposals = []
        for _ in range(count):
            config, _ = ax_client.get_next_trial()
            proposals.append((self._typed_config(config), {"kind": self.args.algorithm}))
        return proposals

    # Fit native optimizers only from the isolated background snapshot
    def _ask_native_batch(self, count: int) -> list[tuple[dict, dict]]:
        if self.args.algorithm == "bayesopt":
            return self._ask_bayesopt_batch(count)
        if self.args.algorithm == "ax":
            return self._ask_ax_batch(count)
        if self.args.algorithm == "scikit":
            from skopt import Optimizer
            from skopt.space import Integer, Real

            names = sorted(self.bounds)
            dimensions = [(Integer if spec.get("type") == "int" else Real)(spec["lower"], spec["upper"])
                for name in names for spec in [self.bounds[name]]]
            optimizer = Optimizer(dimensions, n_initial_points=0, random_state=self.args.rngseed)
            records = self._selected_records()
            costs = [float(record["metrics"][self.metric]) for record in records]
            pending = self._pending_configs()
            configs = [record["config"] for record in records] + pending
            optimizer.tell([[config[name] for name in names] for config in configs], costs + [max(costs)] * len(pending)
            )
            return [(self._typed_config(dict(zip(names, point, strict=True))), {"kind": "scikit"})
                for point in optimizer.ask(n_points=count, strategy="cl_max")]
        raise ValueError(f"Unsupported native optimizer: {self.args.algorithm}")

    # Accept Ray Tune metric and mode propagation without changing the fixed bounds
    def set_search_properties(self, metric, mode, config, **spec) -> bool:
        if isinstance(metric, str):
            self._metric = metric
            self.args.cost = metric
        if isinstance(mode, str):
            self._mode = mode
        if self.param_space is None:
            self.param_space = copy.deepcopy(config)
        return True

    # Refill and consume the rolling proposal reserve
    def _take_proposals(self, count: int) -> list[tuple[dict, dict]]:
        self._collect_fit()
        self._start_fit()
        self._collect_fit()
        take = min(count, len(self.reserve))
        proposals, self.reserve = self.reserve[:take], self.reserve[take:]
        return proposals

    # Build proposals for every currently free Ray trial slot
    def _new_proposal_batch(self) -> list[tuple[dict, dict]]:
        capacity = max(1, int(getattr(self.args, "max_concurrent_trials", 1)))
        if self.args.algorithm == "ampfit":
            return self._take_proposals(max(0, capacity - len(self.live_records)))
        total = int(getattr(self.args, "num_trials", self.next_index + capacity))
        independent_target = self._independent_target(total)
        independent_terminal = sum(
            self._independent_record(record) for record in [*self.completed_records, *self.failed_records])
        independent_live = sum(self._independent_record(record) for record in self.live_records.values())
        independent_phase = independent_terminal + independent_live < independent_target
        if independent_phase:
            independent_needed = max(0, independent_target - independent_terminal - independent_live)
            count = min(capacity - len(self.live_records), max(0, total - self.next_index), independent_needed)
        else:
            count = min(proposal_slots(self.args, completed=len(self.completed_records), pending=len(self.live_records),
                    next_number=self.next_index), max(0, total - self.next_index))
        if count == 0:
            return []
        algorithm = self.args.algorithm
        initial_live = any(
            _search_payload_is_initial(record.get("search_payload")) for record in self.live_records.values())
        initial_completed = bool(self.initial_completed
            or any(_search_payload_is_initial(record.get("search_payload")) for record in self.completed_records))
        wait_for_initial = (self.initial_points is not None and not initial_completed and self.args.algorithm != "basic"
            and int(getattr(self.args, "rand_trials", 0)) == 0)
        if wait_for_initial and initial_live:
            return []
        if wait_for_initial:
            return [(copy.deepcopy(self.initial_points), {"kind": "initial"})]

        adaptive_cold = algorithm != "basic" and not self.completed_records and not independent_phase
        if adaptive_cold and self.live_records:
            return []
        if adaptive_cold:
            return [(self._uniform_config(self.next_index), {"kind": "cold"})]

        if independent_phase or algorithm == "basic":
            kind = "random" if independent_phase else "basic"
            proposals = [(self._uniform_config(self.next_index + offset), {"kind": kind}) for offset in range(count)]
        else:
            proposals = self._take_proposals(count)
        if self.initial_points is not None and not initial_completed and not initial_live:
            initial = (copy.deepcopy(self.initial_points), {"kind": "initial"})
            if proposals:
                proposals[0] = initial
            else:
                proposals = [initial]
        return proposals

    # Serve cached proposals while adaptive fitting continues in the background
    def suggest(self, trial_id: str):
        capacity = max(1, int(getattr(self.args, "max_concurrent_trials", 1)))
        if len(self.live_records) >= capacity:
            return None
        if not self.proposal_queue:
            self.proposal_queue = self._new_proposal_batch()
        if not self.proposal_queue:
            total = int(getattr(self.args, "num_trials", self.next_index + 1))
            finished = self.next_index >= total or self.optimizer_finished
            return Searcher.FINISHED if finished else None
        config, payload = self.proposal_queue.pop(0)
        if isinstance(config, dict):
            config = self._typed_config(config)
        self.live_records[trial_id] = {"config": copy.deepcopy(config), "ray_trial_id": str(trial_id),
            "search_payload": copy.deepcopy(payload), "trial_id": f"trial-{self.next_index:06d}"}
        self.next_index += 1
        self._start_fit()
        return config

    # Store finite completed metrics for future percentile-filtered replay
    def on_trial_complete(self, trial_id: str, result: dict | None = None, error: bool = False) -> None:
        live = self.live_records.pop(trial_id, None)
        if live is None:
            return
        result = live.pop("result", result)
        value = finite_or_none(result.get(self.metric)) if isinstance(result, dict) else None
        result = result if isinstance(result, dict) else {}
        timing = trial_times(finite_or_none(result.get("started_at_unix", live.get("started_at_unix"))),
                             finite_or_none(result.get("completed_at_unix")) or time.time())
        if error or value is None or result.get("trial_failed"):
            self.failed_records.append({**live, **timing, "error": result.get("trial_error") or "Ray trial failed"})
            self._start_fit()
            return
        if _search_payload_is_initial(live.get("search_payload")):
            self.initial_completed = True
        metrics = copy.deepcopy(result.get("optimizer_metrics"))
        if not isinstance(metrics, dict):
            metrics = selected_fields(result, (self.metric, f"{self.metric}_error", f"{self.metric}_uncertainty", "objective_error"))
        metrics[self.metric] = value
        record = {**live, **timing, **selected_fields(result, (
            "card_config", "cost_definition", "likelihood", "gradient", "theta", "theta_hash", "tunename"),
            include_missing=True), **selected_fields(live, ("config", "search_payload")),
            "metrics": metrics, "node_id": result.get("node_id", "ray"), "status": "success"}
        self.completed_records.append(record)
        self._start_fit()

    # Save the searcher state to a checkpoint file
    def save(self, checkpoint_path: str) -> None:
        ensure_dir(pathlib.Path(checkpoint_path).parent)
        with open(checkpoint_path, "wb") as handle:
            pickle.dump(self.__getstate__(), handle)

    # Restore the searcher state from a checkpoint file
    def restore(self, checkpoint_path: str) -> None:
        with open(checkpoint_path, "rb") as handle:
            self.__setstate__(pickle.load(handle))

    # Exclude active thread objects while preserving a restartable fit request
    def __getstate__(self) -> dict:
        return dict(self.__dict__, _checkpoint_schema_version=self.CHECKPOINT_SCHEMA_VERSION,
                    fit_executor=None, fit_future=None, fit_started_at=None)

    # Restore one validated Ray searcher checkpoint into an idle runtime state
    def __setstate__(self, state: dict) -> None:
        if not isinstance(state, dict):
            raise TypeError("Ray searcher checkpoint state must be a dictionary")
        state = dict(state)
        version = state.pop("_checkpoint_schema_version", None)
        if version != self.CHECKPOINT_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported Ray searcher checkpoint schema {version}, expected {self.CHECKPOINT_SCHEMA_VERSION}")
        missing = sorted({"reserve", "fit_request"} - state.keys())
        if missing:
            raise ValueError(f"Ray searcher checkpoint is missing {', '.join(missing)}")
        self.__dict__.update(dict(optimizer_finished=False, optimizer_diagnostics=[], icebo_waiting=False), **state)
        for record in [*self.completed_records, *self.failed_records]:
            record.setdefault("elapsed_seconds",
                trial_times(record.get("started_at_unix"), record.get("completed_at_unix") or time.time())[
                    "elapsed_seconds"])
        self.fit_future = self.fit_started_at = self.fit_executor = None


class BestState:
    """Track the best trial and receive selected outputs on the Ray head."""

    # Initialize the best completed objective and optional head-local output paths
    def __init__(self, *, experiment_dir: str | None = None, cdir: str | None = None, run_name: str | None = None,
        cost: str | None = None, render_figures: bool = False, keep_pickles: bool = True,):
        paths = (experiment_dir, cdir, run_name, cost)
        configured = all(value is not None for value in paths)
        if not configured and any(value is not None for value in paths):
            raise ValueError("Ray head output state requires all output path settings")
        self.experiment_dir = self.param = self.figure_dir = None
        summary = {}
        if configured:
            self.experiment_dir = pathlib.Path(str(experiment_dir))
            self.param = dict(cdir=str(pathlib.Path(str(cdir)).resolve()), cost=str(cost), run_name=str(run_name))
            self.figure_dir = pathlib.Path(icetune_main.icetune_figure_dir(cdir=self.param["cdir"], run_name=str(run_name)))
            summary = read_json(self.figure_dir / "summary.json", default={}) or {}
        self.best_cost = self.published_cost = finite_or_none((summary.get("metrics") or {}).get(str(cost)))
        trial_id = summary.get("trial_id")
        self.best_trial_id = self.published_trial_id = (
            str(trial_id) if self.best_cost is not None and trial_id is not None else None)
        self.pending_outputs, self.committed_outputs = {}, {}
        self.render_figures = bool(render_figures)
        self.keep_pickles = bool(keep_pickles)
        self.output_stage_root = self._prepare_output_stage_root() if configured else None
        if self.output_stage_root is not None:
            self._restore_output_stages()

    # Create one private campaign-specific output staging directory on head-local scratch
    def _prepare_output_stage_root(self) -> pathlib.Path:
        cdir = pathlib.Path(self.param["cdir"]).resolve()
        temporary_root = cdir / "tmp"
        ensure_dir(temporary_root)
        try:
            temporary_root.resolve().relative_to(cdir)
        except ValueError as exc:
            raise RuntimeError("Ray head output staging escaped the local runtime") from exc
        identity = read_json(self.experiment_dir / "icetune_campaign.json", default={}) or {}
        key = "\0".join(
            (str(self.experiment_dir.resolve()), str(identity.get("fingerprint") or "unfenced"), self.param["run_name"])
        )
        token = hashlib.sha256(key.encode("utf-8")).hexdigest()[:16]
        previous = temporary_root / f"ray-output-{token}"
        ensure_dir(self.experiment_dir)
        root = self.experiment_dir / ".outputs"
        if not root.exists() and previous.is_dir():
            previous.rename(root)
        if root.is_symlink() or (root.exists() and not root.is_dir()):
            raise RuntimeError(f'Ray head output stage is not a private directory: "{root}"')
        ensure_dir(root, mode=0o700)
        if root.stat().st_uid != os.getuid():
            raise RuntimeError(f'Ray head output stage is not owned by this process: "{root}"')
        os.chmod(root, 0o700)
        receipt_root = root / "committed"
        if receipt_root.is_symlink() or (receipt_root.exists() and not receipt_root.is_dir()):
            raise RuntimeError(f'Ray output receipt stage is not a directory: "{receipt_root}"')
        ensure_dir(receipt_root, mode=0o700)
        os.chmod(receipt_root, 0o700)
        ensure_dir(root / "renders", mode=0o700)
        return root

    # Compute one validated staged-output record from its private directory
    def _load_output_stage(self, stage: pathlib.Path) -> dict:
        if (self.output_stage_root is None or stage.is_symlink() or stage.parent != self.output_stage_root
            or re.fullmatch(r"pending-[0-9a-f]{32}", stage.name) is None):
            raise RuntimeError(f'Invalid Ray head output stage: "{stage}"')
        manifest = read_json(stage / "manifest.json", default={}) or {}
        token = str(manifest.get("token") or "")
        trial_id = str(manifest.get("trial_id") or "")
        if (manifest.get("schema_version") != 1 or stage.name != f"pending-{token}"
            or re.fullmatch(r"[0-9a-f]{32}", token) is None or not trial_id):
            raise RuntimeError(f'Invalid Ray staged output manifest: "{stage}"')
        descriptor = manifest.get("descriptor")
        pickle_path = stage / "trial.pkl"
        if descriptor is not None:
            if not isinstance(descriptor, dict):
                raise RuntimeError(f'Invalid Ray staged pickle manifest: "{stage}"')
            ray_trial_pickle_relative(descriptor)
            if (pickle_path.is_symlink() or not pickle_path.is_file()
                or pickle_path.stat().st_size != int(descriptor["size"])
                or sha256_file(pickle_path) != descriptor["sha256"]):
                raise RuntimeError(f'Incomplete Ray staged pickle: "{stage}"')
        elif pickle_path.exists():
            raise RuntimeError(f'Unexpected Ray staged pickle: "{stage}"')
        figure_path = stage / "figure-tree" / "figures"
        has_figures = bool(manifest.get("has_figures"))
        if has_figures != figure_path.is_dir() or (has_figures and figure_path.is_symlink()):
            raise RuntimeError(f'Incomplete Ray staged figure tree: "{stage}"')
        cost = float(manifest["cost"])
        if not np.isfinite(cost):
            raise RuntimeError(f'Invalid Ray staged output cost: "{stage}"')
        if has_figures and (not isinstance(manifest.get("archive_size"), int) or manifest["archive_size"] <= 0
                            or re.fullmatch(r"[0-9a-f]{64}", str(manifest.get("archive_sha256") or "")) is None):
            raise RuntimeError(f'Invalid Ray staged figure manifest: "{stage}"')
        if not has_figures and any(manifest.get(key) is not None for key in ("archive_size", "archive_sha256")):
            raise RuntimeError(f'Unexpected Ray staged figure metadata: "{stage}"')
        return {**selected_fields(manifest, ("archive_sha256", "archive_size", "descriptor"), include_missing=True),
                **{key: bool(manifest.get(key)) for key in ("rendered_best", "rendered_initial", "render_pending", "initial")},
                "cost": cost, "created_at_ns": int(manifest["created_at_ns"]),
                "figure_path": str(figure_path) if has_figures else None,
                "pickle_path": str(pickle_path) if descriptor is not None else None,
                "stage_dir": str(stage), "token": token, "trial_id": trial_id}

    # Restore pending manifests and completed receipts after a Ray actor restart
    def _restore_output_stages(self) -> None:
        receipt_root = self.output_stage_root / "committed"
        for path in sorted(receipt_root.glob("*.json")):
            if path.is_symlink() or path.parent != receipt_root:
                logger.warning('.BestState: Ignoring invalid Ray output receipt "%s"', path)
                continue
            receipt = read_json(path, default={}) or {}
            token = str(receipt.get("token") or "")
            if (receipt.get("schema_version") == 1 and path.name == f"{token}.json"
                and re.fullmatch(r"[0-9a-f]{32}", token) is not None and receipt.get("trial_id") is not None
                and isinstance(receipt.get("published"), dict)):
                self.committed_outputs[token] = receipt
        for stage in sorted(self.output_stage_root.glob("pending-*")):
            try:
                record = self._load_output_stage(stage)
            except Exception as exc:
                logger.warning('.BestState: Ignoring invalid staged Ray output "%s": %s', stage, exc)
                continue
            if record["token"] in self.committed_outputs:
                self._remove_output_stage(record)
                continue
            previous = self.pending_outputs.get(record["trial_id"])
            if previous is not None and previous["created_at_ns"] >= record["created_at_ns"]:
                self._remove_output_stage(record)
                continue
            if previous is not None:
                self._remove_output_stage(previous)
            self.pending_outputs[record["trial_id"]] = record

    # Select histograms for later rendering without making completed workers wait
    def finalize_trial_result(self, trial_id, cost):
        finite = finite_or_none(cost)
        improving = finite is not None and (self.best_cost is None or finite < self.best_cost)
        if improving:
            self.best_cost = finite
            self.best_trial_id = str(trial_id)
        store = (
            self.render_figures and finite is not None and (self.published_cost is None or finite < self.published_cost)
        )
        return {"locally_improving": improving, "store_render": store}

    # Remove one validated private staged output tree
    def _remove_output_stage(self, record: dict) -> None:
        stage = pathlib.Path(str(record.get("stage_dir") or ""))
        if (self.output_stage_root is None or stage.is_symlink() or stage.parent != self.output_stage_root
            or re.fullmatch(r"pending-[0-9a-f]{32}", stage.name) is None):
            raise RuntimeError(f'Refusing to remove invalid Ray output stage: "{stage}"')
        if stage.exists():
            shutil.rmtree(stage)

    # Write one worker payload into a private head-local staging file
    def _write_output_payload(self, path: pathlib.Path, payload: bytes) -> None:
        if not isinstance(payload, bytes):
            raise TypeError("Ray trial output transfer payload must be bytes")
        with open(path, "xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(path, 0o600)

    # Validate staged figure summaries before accepting one worker transfer
    def _validate_staged_figures(self, *, source: pathlib.Path, trial_id: str,
                                rendered_initial: bool, rendered_best: bool) -> None:
        for enabled, subdir, label in ((rendered_initial, icetune_main.ICETUNE_INITIAL_FIGURE_SUBDIR, 'initial'),
                                      (rendered_best, '', 'best')):
            if enabled and str((read_json(source / subdir / 'summary.json', default={}) or {}).get('trial_id')) != trial_id:
                raise RuntimeError(f"Ray {label} figure transfer has the wrong trial identity")

    # Persist one validated worker transfer below private head-local scratch
    def _stage_trial_outputs(self, record: dict, transfer: dict) -> dict:
        if self.output_stage_root is None:
            raise RuntimeError("Ray head output state has no configured output paths")
        token = uuid.uuid4().hex
        temporary, stage = (self.output_stage_root / f"{prefix}-{token}" for prefix in ('.stage', 'pending'))
        ensure_dir(temporary, exist_ok=False, mode=0o700)
        try:
            if transfer.get('pickle') is not None:
                self._write_output_payload(temporary / 'trial.pkl', transfer['pickle'])
            if transfer.get('figures') is not None:
                figure_stage = temporary / 'figure-tree'
                ensure_dir(figure_stage, exist_ok=False, mode=0o700)
                source = extract_ray_figure_tree(payload=transfer['figures'], stage=figure_stage)
                self._validate_staged_figures(source=source, **selected_fields(
                    record, ('trial_id', 'rendered_initial', 'rendered_best')))
            manifest = dict(record, created_at_ns=time.time_ns(), has_figures=transfer.get('figures') is not None,
                            schema_version=1, token=token)
            atomic_write_json(str(temporary / 'manifest.json'), json_safe(manifest))
            os.replace(temporary, stage)
            return self._load_output_stage(stage)
        finally:
            if temporary.exists():
                shutil.rmtree(temporary)

    # Publish one immutable staged worker pickle into head-local Tune storage
    def _publish_pickle(self, *, source: pathlib.Path, descriptor: dict) -> None:
        if self.experiment_dir is None or self.param is None:
            raise RuntimeError("Ray head output state has no configured output paths")
        relative = ray_trial_pickle_relative(descriptor)
        if (not source.is_file() or source.stat().st_size != int(descriptor["size"])
            or sha256_file(source) != descriptor["sha256"]):
            raise RuntimeError(f'Ray staged trial pickle is incomplete: "{source}"')
        target = self.experiment_dir / "results" / relative.name
        installed = publish_file_immutable(source, str(target))
        if not installed and sha256_file(target) != descriptor["sha256"]:
            raise RuntimeError(f'Conflicting Ray trial pickle: "{target}"')

    # Publish selected worker figures into the canonical head-local tree
    def _publish_figures(
        self, *, source: pathlib.Path, trial_id: str, cost: float, rendered_initial: bool, rendered_best: bool
    ) -> dict[str, bool]:
        if self.figure_dir is None or self.param is None:
            raise RuntimeError("Ray head output state has no configured output paths")
        published = {"published_best": False, "published_initial": False}
        self._validate_staged_figures(source=source, trial_id=trial_id,
                                      rendered_initial=rendered_initial, rendered_best=rendered_best)
        if rendered_initial:
            initial = source / icetune_main.ICETUNE_INITIAL_FIGURE_SUBDIR
            _publish_ray_figure_tree(
                source=initial, target=self.figure_dir / icetune_main.ICETUNE_INITIAL_FIGURE_SUBDIR, param=self.param)
            published["published_initial"] = True
        already_published = bool(self.published_trial_id == trial_id and self.published_cost is not None
            and np.isclose(float(cost), self.published_cost, rtol=1.0e-15, atol=0.0))
        if rendered_best and (already_published or self.published_cost is None or float(cost) < self.published_cost):
            if not already_published:
                _publish_ray_figure_tree(source=source, target=self.figure_dir, param=self.param)
                self.published_cost = float(cost)
                self.published_trial_id = trial_id
            published["published_best"] = True
        return published

    # Stage selected worker outputs on the head until Tune commits the result
    def publish_trial_outputs(self, *, trial_id, cost, transfer, descriptor, rendered_initial, rendered_best,
                              render_pending=False, initial=False):
        if not isinstance(transfer, dict):
            raise TypeError("Ray trial output transfer must be a dictionary")
        record = dict(trial_id=str(trial_id), cost=float(cost), descriptor=copy.deepcopy(descriptor),
                      rendered_initial=bool(rendered_initial), rendered_best=bool(rendered_best),
                      render_pending=bool(render_pending), initial=bool(initial))
        pickle_payload, figure_payload = transfer.get('pickle'), transfer.get('figures')
        if render_pending and pickle_payload is None:
            raise RuntimeError("Ray rendering requires saved trial histograms")
        if (pickle_payload is None) != (descriptor is None):
            raise RuntimeError("Ray trial pickle transfer and descriptor must be provided together")
        if (figure_payload is None) != (not (rendered_initial or rendered_best)):
            raise RuntimeError("Ray rendered trial output and figure transfer must be provided together")
        if pickle_payload is None and figure_payload is None:
            raise RuntimeError("Ray trial output transfer is empty")
        if not np.isfinite(record['cost']):
            raise ValueError("Ray trial output cost must be finite")
        for label, payload in (('trial pickle', pickle_payload), ('figure', figure_payload)):
            if payload is not None and not isinstance(payload, bytes):
                raise TypeError(f"Ray {label} transfer payload must be bytes")
        if pickle_payload is not None:
            if not isinstance(descriptor, dict):
                raise RuntimeError("Ray trial pickle transfer has no descriptor")
            ray_trial_pickle_relative(descriptor)
            if len(pickle_payload) != int(descriptor['size']):
                raise RuntimeError("Ray trial pickle transfer has the wrong size")
            if hashlib.sha256(pickle_payload).hexdigest() != descriptor['sha256']:
                raise RuntimeError("Ray trial pickle transfer has the wrong checksum")
        record.update(archive_size=None if figure_payload is None else len(figure_payload),
                      archive_sha256=None if figure_payload is None else hashlib.sha256(figure_payload).hexdigest())
        pending = self.pending_outputs.get(record['trial_id'])
        if (pending is not None and np.isclose(record['cost'], pending['cost'], rtol=1.0e-15, atol=0.0)
                and all(value == pending[key] for key, value in record.items() if key != 'cost')):
            return {'output_token': pending['token']}
        staged = self._stage_trial_outputs(record, transfer)
        self.pending_outputs[record['trial_id']] = staged
        if pending is not None:
            self._remove_output_stage(pending)
        return {'output_token': staged['token']}

    # Promote staged outputs only after Tune accepts one successful trial result
    def commit_trial_outputs(self, trial_id, output_token):
        trial_id, output_token = str(trial_id), str(output_token)
        receipt = self.committed_outputs.get(output_token)
        if receipt is not None and str(receipt["trial_id"]) == trial_id:
            return copy.deepcopy(receipt["published"])
        pending = self.pending_outputs.get(trial_id)
        if pending is None or pending["token"] != output_token:
            raise RuntimeError("Ray completed trial has no matching staged head output")
        published = {"published_best": False, "published_initial": False}
        if pending["pickle_path"] is not None and self.keep_pickles:
            self._publish_pickle(source=pathlib.Path(pending["pickle_path"]), descriptor=pending["descriptor"])
        if pending["render_pending"]:
            source = pathlib.Path(pending["pickle_path"])
            destination = self.output_stage_root / "renders" / f"{output_token}.pkl"
            publish_file_immutable(source, str(destination))
            published.update({"render_token": output_token, "render_record": {"token": output_token,
                        "trial_id": trial_id, "cost": pending["cost"], "initial": pending["initial"]}})
        if pending["figure_path"] is not None:
            published = self._publish_figures(source=pathlib.Path(pending["figure_path"]),
                **selected_fields(pending, ("trial_id", "cost", "rendered_initial", "rendered_best")))
        receipt = json_safe({"published": published, "schema_version": 1, "token": output_token, "trial_id": trial_id,
                **time_fields("committed_at", time.time())})
        receipt_path = self.output_stage_root / "committed" / f"{output_token}.json"
        atomic_write_json(str(receipt_path), receipt)
        self.committed_outputs[output_token] = receipt
        self.pending_outputs.pop(trial_id, None)
        try:
            self._remove_output_stage(pending)
        except OSError as exc:
            logger.warning('.BestState: Unable to clean committed Ray output "%s": %s', output_token, exc)
        return published

    # Recover accepted transfers or verify outputs already returned with older Tune state
    def recover_trial_outputs(self, metrics: dict) -> dict:
        token, trial_id = metrics["ray_output_token"], metrics["trial_id"]
        if token in self.committed_outputs or trial_id in self.pending_outputs:
            return self.commit_trial_outputs(trial_id, token)
        descriptor = metrics.get("trial_pickle")
        if descriptor is not None:
            ray_trial_pickle_path(trial_dir=self.experiment_dir, descriptor=descriptor)
        elif not metrics.get("rendered_best") and not metrics.get("rendered_initial"):
            raise RuntimeError("Ray restored trial has no accepted output")
        return {}

    # Read accepted histograms for a separate rendering worker
    def render_payload(self, token: str) -> bytes:
        if token not in self.committed_outputs:
            raise RuntimeError("Ray render payload has no accepted trial")
        return (self.output_stage_root / "renders" / f"{token}.pkl").read_bytes()

    # Release superseded histograms after a better accepted point is available
    def discard_render(self, token: str) -> None:
        if token not in self.committed_outputs:
            raise RuntimeError("Ray render payload has no accepted trial")
        (self.output_stage_root / "renders" / f"{token}.pkl").unlink(missing_ok=True)

    # Publish separately rendered figures using the accepted trial identity
    def publish_render(self, record: dict, transfer: dict) -> dict:
        staged = self.publish_trial_outputs(trial_id=record["trial_id"], cost=record["cost"], transfer=transfer,
            descriptor=None, rendered_initial=record["render_initial"], rendered_best=record["render_best"])
        return self.commit_trial_outputs(record["trial_id"], staged["output_token"])

    # Retain the latest completed covariance independently of newer best fit histograms
    def publish_covariance(self, record: dict, payload: bytes) -> None:
        with tempfile.TemporaryDirectory(prefix="covariance-", dir=self.output_stage_root) as stage:
            source = extract_ray_figure_tree(payload=payload, stage=pathlib.Path(stage))
            covariance = read_json(source / "covariance.json", default={}) or {}
            if covariance.get("trial_id") != record["trial_id"]:
                raise RuntimeError("Ray covariance transfer has the wrong trial identity")
            if any(path.name not in {"covariance.json", "optimizer", "physical"} for path in source.iterdir()):
                raise RuntimeError("Ray covariance transfer contains unrelated figure outputs")
            target = self.figure_dir / "covariance"
            _publish_ray_figure_tree(source=source, target=target, param=self.param)
            logger.info(".icetune_ray: Saved covariance for %s to %s", record["trial_id"], target)

    # Discard worker outputs when Tune rejects or loses the corresponding trial
    def discard_trial_outputs(self, trial_id):
        trial_id = str(trial_id)
        pending = self.pending_outputs.pop(trial_id, None)
        if pending is not None:
            self._remove_output_stage(pending)

    # Compute the current best identity for diagnostics
    def snapshot(self):
        return {"best_cost": self.best_cost, "best_trial_id": self.best_trial_id,
            "pending_outputs": sorted(self.pending_outputs), "published_cost": self.published_cost,
            "published_trial_id": self.published_trial_id}


GlobalState = ray.remote(num_cpus=0, max_restarts=-1, max_task_retries=-1)(BestState)


# Create the shared result state on the dedicated zero CPU Ray head
def create_global_state(*, experiment_dir: str, cdir: str, run_name: str, cost: str, render_figures: bool,
    require_head_resource: bool, keep_pickles: bool = True,):
    head_capacity = ray.cluster_resources().get("icetune_head", 0.0)
    heads = [node["Resources"] for node in ray.nodes()
             if node.get("Alive") and (node.get("Resources") or {}).get("icetune_head", 0.0) > 0.0]
    if require_head_resource:
        if head_capacity <= 0.0 or len(heads) != 1:
            raise RuntimeError("Distributed Ray requires exactly one dedicated icetune_head node")
        if any(heads[0].get(name, 0.0) > 0.0 for name in ("CPU", "GPU")):
            raise RuntimeError("The dedicated icetune_head must advertise zero CPUs and GPUs")
    options = {"resources": {"icetune_head": 0.001}} if head_capacity > 0.0 else {}
    return GlobalState.options(**options).remote(experiment_dir=experiment_dir, cdir=cdir, run_name=run_name, cost=cost,
        render_figures=render_figures, keep_pickles=keep_pickles)


def set_search_algo(args, tunesetup, initial_points, parameter_topology=None):
    """Set the Ray Tune search algorithm."""

    logger.info(".icetune_ray: Optimization algorithm: %s", args.algorithm)

    bounds = normalize_param_space(tunesetup.param_space)
    if args.algorithm == "bayesopt" and any(spec.get("type") == "int" for spec in bounds.values()):
        raise ValueError(
            "Ray bayesopt does not support integer parameter domains; use basic, icebo, hebo, hyperopt, optuna, or ax")
    search_alg = PercentileSearch(args=args, bounds=bounds,
        initial_points=None if args.no_initial_point else initial_points, parameter_topology=parameter_topology)
    search_alg.param_space = copy.deepcopy(tunesetup.param_space)
    return search_alg


# Render one selected trial directly into its persistent Tune trial directory
def _render_ray_trial_figures(*, outputs: dict, param: dict, simdriver, initial: bool, improving: bool) -> None:
    if initial:
        outputs = {**outputs, "search_payload": {**(outputs.get("search_payload") or {}), "kind": "initial"}}
    summary = icetune_main.build_trial_summary_payload(outputs=outputs, param=param)
    if initial or improving:
        publish = icetune_main.publish_trial_figures if improving else icetune_main.publish_initial_trial_figures
        publish(outputs=outputs, param=param, simdriver=simdriver, summary_payload=summary)


# Render saved histograms without waiting for amplitude covariance
@ray.remote(max_retries=0)
def render_completed_trial(payload: bytes, param: dict, simdriver, initial: bool, improving: bool) -> dict:
    param = ray_runtime.worker_runtime(param, simdriver)
    iceruntime.configure_numerical_threads(1)
    outputs = pickle.loads(payload)
    scratch = pathlib.Path(param["cdir"]) / "tmp"
    ensure_dir(scratch)
    with (iceruntime.trial_limit(time.monotonic() + float(param["max_t"])),
        tempfile.TemporaryDirectory(prefix="icetune-render-", dir=scratch) as stage,):
        output_param = {**param, "cdir": stage}
        _render_ray_trial_figures(
            outputs=outputs, param=output_param, simdriver=simdriver, initial=initial, improving=improving)
        return collect_ray_trial_outputs(trial_dir=pathlib.Path(stage), run_name=str(param["run_name"]),
            descriptor=None, rendered_initial=initial, rendered_best=improving)


# Reuse prepared amplitude banks in one sequential covariance worker
@ray.remote(max_restarts=0, max_task_retries=0)
class CovarianceWorker:
    # Resolve the runtime once and retain driver caches between accepted best fits
    def __init__(self, param: dict, simdriver):
        self.param = ray_runtime.worker_runtime(param, simdriver)
        self.simdriver = simdriver
        iceruntime.configure_numerical_threads(1)

    # Compute exact amplitude curvature without generating new events
    def compute(self, payload: bytes) -> bytes:
        from core.tune.optimizers.ampfit.report import trial_covariance

        param, simdriver = self.param, self.simdriver
        outputs = pickle.loads(payload)
        scratch = pathlib.Path(param["cdir"]) / "tmp"
        ensure_dir(scratch)
        started = time.monotonic()
        logger.info(".icetune_ray: Computing covariance for %s", outputs["trial_id"])
        with tempfile.TemporaryDirectory(prefix="icetune-covariance-", dir=scratch) as stage:
            try:
                with iceruntime.trial_limit(started + float(param["max_t"])):
                    trial_covariance(driver=simdriver, record=outputs, param=param, output_root=pathlib.Path(stage))
            except Exception as exc:
                atomic_write_json(str(pathlib.Path(stage) / "covariance.json"), {
                    "trial_id": outputs["trial_id"], "status": "unavailable", "error": str(exc)})
                logger.warning(".icetune_ray: ampfit covariance unavailable for %s: %s", outputs["trial_id"], exc)
            logger.info(".icetune_ray: Covariance for %s completed in %.4f s", outputs["trial_id"], time.monotonic() - started)
            return pack_ray_figure_tree(pathlib.Path(stage))


# Evaluate one physics point using its Ray-assigned trial identity
class CFunc(tune.Trainable):
    # Initialize one Ray Tune trainable instance
    def setup(self, config: dict, simdriver, param: dict, global_state, initial_points=None):
        self.started_at = time.time()
        self.deadline = time.monotonic() + float(param["max_t"])
        self.setup_error = None
        self.param = param
        self.config = config
        self.simdriver = simdriver
        try:
            with iceruntime.trial_limit(self.deadline):
                self.param = ray_runtime.worker_runtime(param, simdriver)
                if self.param.get("optimization", {}).get("optimizer") == "ampfit":
                    iceruntime.configure_numerical_threads(self.param["processes"])
                self.global_state = global_state
                self.initial_points = copy.deepcopy(initial_points)
                # Initialize Numba before the Ray trial forks histogram workers
                obs.preflight_numba_runtime()
                logger.info(".icetune_ray: Numba runtime preflight completed for trial worker")
        except iceruntime.TrialTimeout as exc:
            self.setup_error = str(exc)

    # Reuse initialized amplitude banks while resetting one trial's configuration and time budget
    def reset_config(self, new_config: dict) -> bool:
        if self.setup_error is not None:
            return False
        self.config = new_config
        self.started_at = time.time()
        self.deadline = time.monotonic() + float(self.param["max_t"])
        return True

    # Report computed histograms before separately rendering accepted trial figures
    def step(self):
        config, param, trial_id = self.config, self.param, self.trial_id
        outputs = None
        try:
            with iceruntime.trial_limit(self.deadline):
                if self.setup_error:
                    raise iceruntime.TrialTimeout(self.setup_error)
                is_initial = config_contains(config, self.initial_points)
                outputs = icetune_main.evaluate_trial_outputs(
                    config=config, param=param, simdriver=self.simdriver, trial_id=trial_id)
                metrics = copy.deepcopy(outputs["metrics"])
                icetune_main.require_finite_trial_cost(outputs=outputs, cost_key=param["cost"])
                saved_events = None
                if param.get("save_events"):
                    if not outputs.get("event_samples"):
                        raise event_output.EventTransferError("Successful trial returned no event samples")
                    saved_events = event_output.publish_samples(
                        source=pathlib.Path(outputs["event_samples"]), destination=pathlib.Path(param["events_dir"]),
                        trial_id=trial_id, metadata={"config": outputs["config"], "card_config": outputs["card_config"],
                            "aux_param_space": param["aux_param_space"], "mc_steer": param["mc_steer"],
                            "rngseed": param.get("rngseed"), "run_name": param["run_name"],
                            "campaign_fingerprint": param["event_campaign"], "metrics": metrics})
                render_pending = False
                if param["plot"]:
                    decision = ray.get(self.global_state.finalize_trial_result.remote(trial_id, metrics[param["cost"]]))
                    render_pending = bool(is_initial or decision["store_render"])
                descriptor = output_token = None
                if param["pickle_dump"] or render_pending:
                    relative = pathlib.Path("results", icetune_main.trial_pickle_filename(outputs))
                    descriptor = icetune_main.maybe_dump_trial_payload(
                        outputs=outputs, param={**param, "pickle_dump": True}, destination=pathlib.Path(self.logdir, relative))
                    descriptor["relative_path"] = relative.as_posix()
                    transfer = collect_ray_trial_outputs(trial_dir=pathlib.Path(self.logdir),
                        run_name=str(param["run_name"]), descriptor=descriptor, rendered_initial=False,
                        rendered_best=False)
                    staged = ray.get(self.global_state.publish_trial_outputs.remote(trial_id=trial_id,
                            cost=metrics[param["cost"]], transfer=transfer, descriptor=descriptor,
                            rendered_initial=False, rendered_best=False, render_pending=render_pending,
                            initial=is_initial))
                    output_token = staged["output_token"]
                    pathlib.Path(self.logdir, relative).unlink(missing_ok=True)
                result = {**metrics, **selected_fields(outputs, (
                    "card_config", "cost_definition", "likelihood", "gradient", "theta", "theta_hash"),
                    include_missing=True), "is_initial": is_initial, "node_id": socket.getfqdn(),
                    "optimizer_metrics": metrics, "ray_output_token": output_token,
                    "rendered_best": False, "rendered_initial": False, "trial_id": trial_id,
                    "trial_pickle": descriptor if param["pickle_dump"] else None,
                    "event_samples": saved_events, "tunename": outputs["tunename"]}
        except Exception as exc:
            if ray_infrastructure_error(exc):
                if ray.is_initialized():
                    exc.add_note(f"icetune infrastructure failure on Ray node {ray.get_runtime_context().get_node_id()}")
                raise
            if isinstance(exc, event_output.EventTransferError):
                raise
            result = {param["cost"]: float("nan"), "trial_failed": True,
                "trial_error": ray_failure_text(exc=exc, outputs=outputs), "trial_id": trial_id}
        finally:
            if isinstance(outputs, dict):
                if outputs.get("event_samples"):
                    shutil.rmtree(outputs["event_samples"], ignore_errors=True)
                if outputs.get("tunename"):
                    icetune_main.cleanup_trial_outputs(simdriver=self.simdriver, param=param, tunename=outputs["tunename"])
        return {**result, **trial_times(self.started_at, time.time())}

    # Save the trial configuration for native Tune recovery
    def save_checkpoint(self, checkpoint_dir):
        ensure_dir(checkpoint_dir)
        with open(os.path.join(checkpoint_dir, "checkpoint.pkl"), "wb") as f:
            pickle.dump(self.config, f, protocol=pickle.HIGHEST_PROTOCOL)
        return checkpoint_dir

    # Restore the trial configuration before a new attempt
    def load_checkpoint(self, checkpoint_dir):
        with open(os.path.join(checkpoint_dir, "checkpoint.pkl"), "rb") as f:
            self.config = pickle.load(f)


# Compute whether one Ray experiment contains restorable Tuner state
def ray_tuner_state_exists(experiment_dir: str) -> bool:
    return os.path.isfile(os.path.join(experiment_dir, "tuner.pkl"))


# Compute whether Ray's own newest checkpoint contains only unfinished trials
def ray_native_tune_state_is_empty(root: pathlib.Path) -> bool:
    state_paths = sorted(root.glob("experiment_state-*.json"))
    if not state_paths:
        return False
    try:
        state = read_json(state_paths[-1], default={}) or {}
        trial_data = state["trial_data"]
        if not isinstance(trial_data, list):
            return False
        for item in trial_data:
            if not isinstance(item, list | tuple) or len(item) != 2:
                return False
            trial = json.loads(item[0])
            runtime = json.loads(item[1])
            if trial.get("status") in {"ERROR", "TERMINATED"}:
                return False
            result = runtime.get("last_result") or {}
            if not isinstance(result, dict):
                return False
            evaluated = int(result.get("training_iteration", 0) or 0) > 0 or any(key in result
                for key in ("completed_at_unix", "optimizer_metrics", "theta_hash", "trial_pickle", "tunename"))
            if evaluated or int(runtime.get("num_failures", 0)) > 0:
                return False
        for search_path in sorted(root.glob("searcher-state-*.pkl")):
            with open(search_path, "rb") as handle:
                search_state = pickle.load(handle)
            if not isinstance(search_state, dict):
                return False
            if search_state.get("completed_records") or search_state.get("failed_records"):
                return False
    except (AttributeError, EOFError, ImportError, KeyError, OSError, TypeError, ValueError, json.JSONDecodeError,
        pickle.PickleError):
        return False
    return True


# Compute whether one saved Ray campaign has no completed or failed evaluations
def ray_tune_state_is_empty(*, experiment_dir: str, param: dict) -> bool:
    root = pathlib.Path(experiment_dir)
    if not ray_tuner_state_exists(experiment_dir) or not ray_native_tune_state_is_empty(root):
        return False
    status = read_json(root / "status.json", default={}) or {}
    if (status.get("backend") != "ray" or status.get("completed") != 0 or status.get("failed") != 0
        or status.get("best_result") is not None):
        return False
    history_path = root / "history.json"
    if history_path.is_file():
        history = read_json(history_path, default={}) or {}
        if history.get("trials") or history.get("failed"):
            return False
    if any(root.rglob("TUNE_icetune_*.pkl")):
        return False
    figure_dir = pathlib.Path(icetune_main.icetune_figure_dir(cdir=param["cdir"], run_name=param["run_name"]))
    return not ((figure_dir / "summary.json").is_file()
        or (figure_dir / icetune_main.ICETUNE_INITIAL_FIGURE_SUBDIR / "summary.json").is_file())


# Move one incompatible empty Ray state aside and retain only current INIT data
def reset_empty_ray_tune_state(*, experiment_dir: str) -> pathlib.Path:
    root = pathlib.Path(experiment_dir)
    backup = root.with_name(f"{root.name}._old-{time.time_ns()}-{os.getpid()}")
    root.replace(backup)
    ensure_dir(root, exist_ok=False)
    for name in ("ray_init.json", "results", "submit", "ray", "init", "init_logs", "campaigns"):
        source = backup / name
        if source.exists():
            source.replace(root / name)
    return backup


# Rebind restored Tune state to current head-local callbacks and storage controls
def rebind_restored_ray_tuner(*, tuner, callbacks: list, history_plot, retry_limit: int = 3) -> Searcher:
    internal = getattr(tuner, "_local_tuner", None)
    if internal is None:
        raise RuntimeError("Ray Tune restore requires a direct Ray connection for head-local storage")
    restored_search = internal._tune_config.search_alg
    if not isinstance(restored_search, Searcher):
        raise RuntimeError("Restored Ray Tune state has no compatible search algorithm")
    history_plot.search_alg = restored_search
    internal._run_config.callbacks = callbacks
    internal._run_config.sync_config = tune.SyncConfig(sync_artifacts=False)
    internal._run_config.checkpoint_config = tune.CheckpointConfig(checkpoint_at_end=False, checkpoint_frequency=0)
    internal._run_config.failure_config = tune.FailureConfig(max_failures=retry_limit)
    return restored_search


# Install or validate the immutable campaign identity beside Ray Tuner state
def ensure_ray_campaign_identity(
    *, args, tunesetup, initial_points, param: dict, experiment_dir: str, physics_fingerprint=None) -> dict:
    path = os.path.join(experiment_dir, "icetune_campaign.json")
    missing_identity = ray_tuner_state_exists(experiment_dir) and not os.path.isfile(path)
    physics_keys = ("aux_param_space", "cost", "cost_avg", "cost_definition", "cost_rho", "data_covariance_mode",
        "datacards", "max_t", "mc_steer", "obs_module", "pickle_dump", "save_events", "plot", "plot_brand", "simdriver")
    identity = json_safe({"algorithm": args.algorithm,
            "initial_points": None if args.no_initial_point else initial_points,
            "optimization": param.get("optimization") or icetune_main.build_optimization_metadata(args),
            "parameter_space": param.get("parameter_space")
            or icetune_main.build_parameter_space(normalize_param_space(tunesetup.param_space)),
            "parameter_topology": param.get("parameter_topology", {}), "physics_fingerprint": physics_fingerprint,
            "physics": {key: param.get(key) for key in physics_keys}, "run_name": args.run_name,
            "tunesetup": getattr(args, "tunesetup", "")})
    fingerprint = json_fingerprint(identity)
    stored = read_json(path, default={}, retry_missing=True) or {}
    expected = {"schema_version": 1, "fingerprint": fingerprint, "identity": identity}
    compatible = expected.items() <= stored.items()
    if not compatible and ray_tuner_state_exists(experiment_dir):
        if not ray_tune_state_is_empty(experiment_dir=experiment_dir, param=param):
            if missing_identity:
                raise RuntimeError("Existing Ray Tuner state has no icetune campaign identity; use a new run name")
            raise RuntimeError(
                "Ray Tune storage belongs to a different physics or optimizer campaign; use a new run name")
        backup = reset_empty_ray_tune_state(experiment_dir=experiment_dir)
        logger.warning('.icetune_ray: Replaced incompatible empty Tune state; old local state is "%s"', backup)
        stored = {}
    elif not compatible and stored:
        raise RuntimeError("Ray Tune storage belongs to a different physics or optimizer campaign; use a new run name")
    write_once_json(path, expected)
    stored = read_json(path, default={}, retry_missing=True) or {}
    if not expected.items() <= stored.items():
        raise RuntimeError("Ray campaign identity could not be installed atomically")
    return stored


# Run the Ray backend through final publication validation
def run_ray_backend(*, args, tunesetup, simdriver, initial_points, param: dict, storage_path: str, start_time: float
) -> None:
    runtime_variables = simdriver.runtime_environment(
        cdir=args.cdir, libdir=args.libdir, python_version=param["PYTHON_VERSION"])

    upload_runtime = bool(getattr(args, "ray_upload_runtime", False))
    if upload_runtime:
        runtime_variables = portable_ray_environment(runtime_variables, cdir=args.cdir)
        os.environ["RAY_RUNTIME_ENV_IGNORE_GITIGNORE"] = "1"
        runtime_env = ray.runtime_env.RuntimeEnv(
            env_vars=runtime_variables, excludes=ray_runtime.runtime_excludes(simdriver.runtime_files(param)),
            working_dir=getattr(args, "ray_runtime_uri", None) or args.cdir)
    else:
        runtime_env = ray.runtime_env.RuntimeEnv(env_vars=runtime_variables)

    local = (args.address or "local").strip().lower() == "local"
    storage_path = str(pathlib.Path(storage_path).resolve())
    ensure_dir(pathlib.Path(storage_path))

    if local:
        simdriver.runtime_environment(
            cdir=args.cdir, libdir=args.libdir, python_version=param["PYTHON_VERSION"], apply=True)
    ray.init(**ray_init_arguments(args=args, runtime_env=runtime_env))

    try:
        resources = wait_ray_trial_capacity(args)
        logger.info(".icetune_ray: Trial capacity ready with resources=%s", resources["cluster"])

        if upload_runtime:
            logger.info(".icetune_ray: Ray owns uploaded runtime installation on admitted nodes")
        else:
            node_contracts = validate_ray_cluster_paths(cdir=args.cdir, libdir=args.libdir, storage_path=storage_path)
            logger.info(".icetune_ray: Validated shared runtime on %s Ray node(s)", len(node_contracts))

        if args.algorithm == "lbfgs":
            from core.tune.optimizers.lbfgs import optimizer as icetune_lbfgs

            summary = icetune_lbfgs.run_lbfgs_backend(args=args, tunesetup=tunesetup, simdriver=simdriver,
                initial_points=initial_points, param=param, start_time=start_time, storage_path=storage_path)
            logger.info(".icetune_ray: lbfgs best config:\n%s", summary.get("best_config"))
            logger.info(".icetune_ray: lbfgs best metrics:\n%s", summary.get("best_metrics"))
            return

        # Ray Tune otherwise stages one pending placement group for a custom searcher
        configure_tune_pending_trials(args)
        install_ray_oom_handler()

        experiment_dir = os.path.join(storage_path, args.run_name)
        physics_fingerprint = getattr(simdriver, "physics_fingerprint", lambda **_: None)(param=param, tunesetup=tunesetup)
        identity = ensure_ray_campaign_identity(args=args, tunesetup=tunesetup, initial_points=initial_points,
            param=param, experiment_dir=experiment_dir, physics_fingerprint=physics_fingerprint)
        if param.get("save_events"):
            param["event_campaign"] = identity["fingerprint"]
            param["events_dir"] = str(pathlib.Path(param.get("events_dir") or (
                pathlib.Path(experiment_dir) / "campaigns" / short_id(identity["fingerprint"]) / "events")).resolve())
            ensure_dir(pathlib.Path(param["events_dir"]))
        if ray_tuner_state_exists(experiment_dir) and not args.restore:
            raise RuntimeError("Ray Tuner state already exists; use --restore or choose a new run name")

        if args.algorithm == "ampfit":
            from core.tune.drivers.graniitti.ampfit.distributed import start_workers

            param["ampfit_workers"] = start_workers(simdriver, param)
            param["processes"] = 1
        callback_cdir = "." if upload_runtime else args.cdir
        callback_experiment = str(pathlib.Path("runs", "icetune", args.run_name)) if upload_runtime else experiment_dir
        callback_param = copy.deepcopy(param)
        callback_param["cdir"] = callback_cdir
        global_state = create_global_state(experiment_dir=experiment_dir,
            cdir=str(pathlib.Path(callback_cdir).resolve()), run_name=args.run_name, cost=args.cost,
            render_figures=bool(param["plot"]), require_head_resource=not local,
            keep_pickles=bool(param["pickle_dump"]))

        logger.info(".icetune_ray: Initial parameter values=%s", None if args.no_initial_point else initial_points)

        search_alg = set_search_algo(args=args, tunesetup=tunesetup, initial_points=initial_points,
            parameter_topology=param.get("parameter_topology", {}))
        trainable = tune.with_resources(tune.with_parameters(CFunc, simdriver=simdriver, param=param,
                global_state=global_state, initial_points=None if args.no_initial_point else initial_points),
            tune.PlacementGroupFactory([{"CPU": 1 if param.get("ampfit_workers") else args.cpu_per_trial,
                                         "GPU": args.gpu_per_trial}]))

        history_plot = HistoryPlotCallback(experiment_dir=callback_experiment, cdir=callback_cdir, cost=args.cost,
            run_name=args.run_name, interval_s=getattr(args, "history_interval_s", 60.0), search_alg=search_alg,
            param=callback_param, num_trials=args.num_trials, max_concurrent_trials=args.max_concurrent_trials)
        trial_outputs = TrialOutputCallback(experiment_dir=callback_experiment, param=callback_param,
            global_state=global_state, simdriver=simdriver, runtime_param=param)
        callbacks = [trial_outputs, history_plot]

        # Build one fresh tuner with head process history plotting enabled
        def create_tuner():
            return tune.Tuner(trainable=trainable, param_space=tunesetup.param_space, tune_config=tune.TuneConfig(
                    reuse_actors=args.algorithm == "ampfit", metric=args.cost, mode="min", search_alg=search_alg,
                    num_samples=args.num_trials, max_concurrent_trials=args.max_concurrent_trials),
                run_config=tune.RunConfig(callbacks=callbacks,
                    failure_config=tune.FailureConfig(max_failures=getattr(args, "ray_trial_retry_limit", 3)),
                    checkpoint_config=tune.CheckpointConfig(checkpoint_at_end=False, checkpoint_frequency=0),
                    name=args.run_name, sync_config=tune.SyncConfig(sync_artifacts=False), verbose=args.verbose,
                    stop={"training_iteration": 1}, storage_path=storage_path))

        # Resume unfinished trials while retaining terminal trial outcomes and their identities
        def restore_tuner():
            if not ray_tuner_state_exists(experiment_dir):
                raise RuntimeError("Ray campaign recovery has no saved Tune state")
            tuner = tune.Tuner.restore(experiment_dir, trainable=trainable, param_space=tunesetup.param_space)
            rebind_restored_ray_tuner(tuner=tuner, callbacks=callbacks, history_plot=history_plot,
                retry_limit=getattr(args, "ray_trial_retry_limit", 3))
            logger.info(".icetune_ray: Rebound restored Tune state to head-local output transport")
            return tuner

        tuner = restore_tuner() if args.restore and ray_tuner_state_exists(experiment_dir) else create_tuner()

        try:
            results = fit_ray_tuner(tuner=tuner, restore=restore_tuner,
                limit=getattr(args, "ray_restore_retry_limit", 2),
                delay_s=getattr(args, "ray_restore_retry_delay_s", 10))
            best_result = results.get_best_result(metric=args.cost, mode="min")
            output_state = trial_outputs.validate(
                best_result=best_result, require_initial=bool(param["plot"] and not args.no_initial_point
                                                             and isinstance(initial_points, dict)))
            output_state.update(ray.get(global_state.snapshot.remote()))
            logger.info(".icetune_ray: Final worker output state:\n%s", output_state)
            if history_plot.plot(force=True) is None:
                raise RuntimeError("Ray final history plot is incomplete")
        finally:
            try:
                history_plot.close()
            finally:
                history_plot.workers.close()
                if history_plot.search_alg.fit_executor is not None:
                    history_plot.search_alg.fit_executor.shutdown(wait=True)

        for label, value in (("config", best_result.config), ("card config", best_result.metrics.get("card_config")),
                             ("metrics", best_result.metrics)):
            logger.info(".icetune_ray: Best result %s:\n%s", label, value)
    finally:
        for actor in {entry[0] for entries in param.get("ampfit_workers", {}).values() for entry in entries}:
            ray.kill(actor)
        logger.info(".icetune_ray: Disconnecting from the Ray cluster")
        ray.shutdown()
        logger.info(".icetune_ray: Ray client shutdown completed")

# icetune generic steering utilities and driver interfaces
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import contextlib
import copy
import hashlib
import json
import math
import os
import pathlib
import pickle
import re
import shutil
import socket
import threading
import time
import uuid
from collections.abc import Callable
from functools import partial
from glob import glob

from core.io import logger as log
from core.io.files import ensure_dir, is_retryable_io_error
from core.io.serialize import finite_or_none, selected_fields
from core.tune import likelihood as icetune_likelihood
from core.tune import summary as fit_summary
from core.tune.cache import PermanentConfigurationError, publish_file_immutable, sha256_file
from core.tune.io import atomic_write_json as _atomic_write_json
from core.tune.io import read_json as _read_json
from core.tune.io import time_fields
from core.tune.runtime import process as iceruntime

logger = log.get_logger(__name__)

FIGURE_PUBLISH_REQUEST_WAIT_S = 5.0
FIGURE_CLAIM_TIMEOUT_S = 600.0
PLOT_PUBLISH_LOCK_WAIT_S = 120.0
PLOT_PUBLISH_LOCK_STALE_S = FIGURE_CLAIM_TIMEOUT_S
ICETUNE_INITIAL_FIGURE_SUBDIR = "init"
ICETUNE_HISTORY_FIGURES = (
    "cost_evolution.png",
    "cost_evolution_logy.png",
    "cost_sorted.png",
    "cost_sorted_logy.png",
    "parameter_evolution_physical.pdf",
    "parameter_evolution_optimizer.pdf",
)

# Fields preserved in an immutable full trial payload
TRIAL_PICKLE_FIELDS = (  # noqa: SIM905
    "trial_id tunename card_config config results search_payload metrics gradient cost_definition theta theta_hash likelihood cost_arr ndf_arr weight_arr valid_arr error"
).split()


# Compute the finite scalar cost value for optimizer fitting decisions
def finite_cost_or_none(record: dict, cost_key: str):
    value = ((record or {}).get("metrics") or {}).get(cost_key)
    return finite_or_none(value)


# Compute one finite optimizer cost or raise the driver-reported trial failure
def require_finite_trial_cost(*, outputs: dict, cost_key: str) -> float:
    detail = outputs.get("error")
    failure = outputs.get("failure")
    if not detail and isinstance(failure, dict):
        detail = failure.get("summary") or failure.get("root_cause") or "driver reported failure"
    if detail:
        raise RuntimeError(f"Trial evaluation failed: {detail}")
    value = finite_cost_or_none(outputs, cost_key)
    if value is not None:
        return value
    raise RuntimeError(f'configured cost metric "{cost_key}" is not finite')


# Compute completed records in the lowest-cost percentile for surrogate fitting
def select_surrogate_fit_records(records: list[dict], cost_key: str, percentile: float | None) -> list[dict]:
    try:
        cut = float(percentile)
    except (TypeError, ValueError):
        cut = 1.0
    cut = min(1.0, max(0.0, cut))

    finite_records = [(finite_cost_or_none(record, cost_key), index, record) for index, record in enumerate(records)]
    finite_records = [item for item in finite_records if item[0] is not None]
    if not finite_records:
        return []
    if cut >= 1.0:
        return [record for _, _, record in finite_records]

    keep_count = max(1, int(math.ceil(cut * len(finite_records))))
    keep_indices = {index for _, index, _ in sorted(finite_records, key=lambda item: (item[0], item[1]))[:keep_count]}
    return [record for _, index, record in finite_records if index in keep_indices]


# Compute ordered parameter-space metadata for one normalized bound map
def build_parameter_space(bounds: dict) -> list[dict]:
    out = []
    for name, spec in sorted(bounds.items()):
        if spec.get("type") == "int":
            lower = int(spec["lower"])
            upper = int(spec["upper"])
            prior = {"type": "int", "lower": lower, "upper": upper}
        else:
            lower = float(spec["lower"])
            upper = float(spec["upper"])
            prior = {"type": "uniform", "lower": lower, "upper": upper}
        out.append({"name": str(name), "lower": lower, "upper": upper, "prior": prior})
    return out


def theta_from_config(config: dict, parameter_space: list[dict]) -> list[float]:
    """Return a parameter vector ordered exactly as parameter_space"""

    return [float(config[item["name"]]) for item in parameter_space]


def theta_hash(theta: list[float]) -> str:
    """Return a stable fingerprint for an ordered parameter vector"""

    payload = json.dumps([float(x) for x in theta], allow_nan=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def build_theta_record(config: dict, parameter_space: list[dict]) -> tuple[list[float] | None, str | None]:
    """Return theta and theta_hash when the parameter metadata is available"""

    if not parameter_space:
        return None, None
    theta = theta_from_config(config, parameter_space)
    return theta, theta_hash(theta)


# Compute the stable mathematical definition of one optimizer cost
def optimizer_cost_definition(cost: object) -> str | None:
    if cost is None:
        return None
    name = str(cost)
    if name == "ratio2":
        return "symmetric_asinh_ratio_chi2_v3"
    if name == "wasserstein":
        return "piecewise_constant_density_cdf_l1_v2"
    return name


def build_optimization_metadata(args) -> dict:
    """Return common optimizer settings for posterior-ready history files"""

    names = (
        "backend",
        "cost",
        "cost_avg",
        "cost_rho",
        "max_concurrent_trials",
        "mc_correlation_events",
        "num_trials",
        "proposal_batch_fraction",
        "proposal_batch_size",
        "rand_trials",
        "rngseed",
    )
    metadata = {name: copy.deepcopy(getattr(args, name, None)) for name in names}
    metadata.update(
        cost_definition=optimizer_cost_definition(metadata["cost"]),
        data_covariance_mode=copy.deepcopy(getattr(args, "data_covariance_mode", "diagonal")),
        icebo=copy.deepcopy(getattr(args, "icebo_settings", None)),
        hebo=copy.deepcopy(getattr(args, "hebo_settings", None)),
        ampfit=copy.deepcopy(getattr(args, "ampfit_settings", None)),
        lbfgs=copy.deepcopy(getattr(args, "lbfgs_settings", None)),
        mc_correlation_weighting=copy.deepcopy(getattr(args, "mc_correlation_weighting", "card")),
        optimizer=copy.deepcopy(getattr(args, "algorithm", None)),
        surrogate_fit_percentile=copy.deepcopy(getattr(args, "surrogate_fit_percentile", 1.0)),
    )
    return metadata


def _sanitize_tunename_component(value) -> str:
    value = str(value)
    value = value.replace(".", "-")
    value = re.sub(r"[^A-Za-z0-9_-]+", "-", value)
    value = value.strip("-")
    return value or "unknown"


def unique_trial_identifier(
    *, prefix: str, trial_id: str | None = None, node_id: str | None = None, pid: int | None = None
) -> str:
    """Return a collision-resistant driver-owned trial identifier."""

    try:
        suffix = iceruntime.get_ip()
    except Exception:
        suffix = socket.gethostname().split(".")[0]

    parts = [_sanitize_tunename_component(prefix), _sanitize_tunename_component(suffix)]
    if node_id is not None:
        parts.append(_sanitize_tunename_component(node_id))
    if trial_id is not None:
        parts.append(_sanitize_tunename_component(trial_id))
    if pid is not None:
        parts.append(f"pid{int(pid)}")

    return "_".join(parts)


def default_tmp_token() -> str:
    """Return a cross-node collision-resistant token for temporary filesystem outputs."""

    host = _sanitize_tunename_component(socket.gethostname().split(".")[0])
    return ".".join([host, f"p{os.getpid():x}", f"t{time.time_ns():x}", uuid.uuid4().hex[:12]])


def default_tmp_path(path: str) -> str:
    """Return a collision-resistant sibling temporary path for atomic file replacement."""

    return f"{path}.tmp.{default_tmp_token()}"


def default_tmp_stage_path(path: str | os.PathLike[str]) -> pathlib.Path:
    """Return a hidden collision-resistant sibling temporary directory path."""

    target = pathlib.Path(path)
    return target.parent / f".{target.name}.tmp.{default_tmp_token()}"


def icetune_figure_dir(*, cdir: str, run_name: str) -> str:
    """Return the published figure directory for an icetune run."""

    return os.path.join(cdir, "figs", "icetune", run_name)


def icetune_initial_figure_dir(*, cdir: str, run_name: str) -> str:
    """Return the preserved initial-trial figure directory for an icetune run."""

    return os.path.join(icetune_figure_dir(cdir=cdir, run_name=run_name), ICETUNE_INITIAL_FIGURE_SUBDIR)


def icetune_figure_publish_lock_dir(*, cdir: str, run_name: str) -> str:
    """Return the per-run publish lock directory for figure publication."""

    return os.path.join(cdir, "runs", "icetune", run_name, "figure_publish.lock")


def icetune_figure_publish_tmp_dir(*, cdir: str, run_name: str) -> str:
    """Return the temporary directory used while rendering a granted publish."""

    return os.path.join(cdir, "runs", "icetune", run_name, "figure_publish.tmp")


def figure_summary_timestamped(payload: dict, timestamp: float | int | None = None) -> dict:
    """Return a figure-summary payload with local and unix creation timestamps"""

    output = copy.deepcopy(payload)
    fields = time_fields("created_at", time.time() if timestamp is None else timestamp)
    for key, value in fields.items():
        output.setdefault(key, value)
    return output


def _publish_lock_owner_path(lock_dir: str) -> str:
    return os.path.join(lock_dir, "owner.json")


def _publish_lock_token() -> dict:
    return {
        "host": socket.gethostname().split(".")[0],
        "lock_id": uuid.uuid4().hex,
        "pid": os.getpid(),
        "thread_id": threading.get_ident(),
        **time_fields("acquired_at", time.time()),
    }


def _lock_token_matches(owner: dict | None, token: dict) -> bool:
    if not isinstance(owner, dict):
        return False
    return (
        owner.get("host") == token.get("host")
        and owner.get("lock_id") == token.get("lock_id")
        and int(owner.get("pid", -1)) == int(token.get("pid", -2))
        and int(owner.get("thread_id", -1)) == int(token.get("thread_id", -2))
    )


# Renew one owned figure lock while rendering or transfer remains active
def _renew_publish_lock(stop: threading.Event, *, owner_path: str, token: dict, interval_s: float) -> None:
    while not stop.wait(interval_s):
        try:
            owner = _read_json(owner_path, default={}) or {}
            if not _lock_token_matches(owner, token):
                return
            os.utime(owner_path, None)
        except OSError as exc:
            if not is_retryable_io_error(exc):
                return


def acquire_publish_lock(
    *, cdir: str, run_name: str, wait_s: float = PLOT_PUBLISH_LOCK_WAIT_S, stale_s: float = PLOT_PUBLISH_LOCK_STALE_S
) -> dict:
    """Acquire the per-run plot publication lock."""

    lock_dir = icetune_figure_publish_lock_dir(cdir=cdir, run_name=run_name)
    token = _publish_lock_token()
    owner_path = _publish_lock_owner_path(lock_dir)
    parent = os.path.dirname(lock_dir)
    if parent:
        ensure_dir(parent)

    start = time.time()
    while True:
        try:
            ensure_dir(lock_dir, exist_ok=False)
            _atomic_write_json(owner_path, token)
            stop = threading.Event()
            thread = threading.Thread(
                target=_renew_publish_lock,
                kwargs={
                    "stop": stop,
                    "owner_path": owner_path,
                    "token": token,
                    "interval_s": max(0.01, min(60.0, stale_s / 3.0)),
                },
                daemon=True,
            )
            thread.start()
            token["_renewal"] = (stop, thread)
            return token
        except FileExistsError:
            owner = _read_json(owner_path, default={}) or {}
            owner_ts = float(owner.get("acquired_at_unix", 0.0) or 0.0)
            with contextlib.suppress(OSError):
                owner_ts = max(owner_ts, os.path.getmtime(owner_path))
            if owner_ts <= 0.0:
                with contextlib.suppress(OSError):
                    owner_ts = os.path.getmtime(lock_dir)
            if owner_ts > 0.0 and (time.time() - owner_ts) > stale_s:
                stale_name = f"{lock_dir}.stale.{int(time.time())}.{os.getpid()}"
                try:
                    os.rename(lock_dir, stale_name)
                    shutil.rmtree(stale_name, ignore_errors=True)
                    continue
                except OSError:
                    pass
            if (time.time() - start) > wait_s:
                raise TimeoutError(f'Could not acquire figure publish lock "{lock_dir}" within {wait_s} s') from None
            time.sleep(1.0)


def release_publish_lock(*, cdir: str, run_name: str, token: dict) -> None:
    """Release the per-run plot publication lock when owned by this process."""

    renewal = token.get("_renewal")
    if isinstance(renewal, tuple) and len(renewal) == 2:
        renewal[0].set()
        renewal[1].join(timeout=2.0)
    lock_dir = icetune_figure_publish_lock_dir(cdir=cdir, run_name=run_name)
    owner_path = _publish_lock_owner_path(lock_dir)
    owner = _read_json(owner_path, default={}) or {}
    if not _lock_token_matches(owner, token):
        return
    with contextlib.suppress(FileNotFoundError):
        os.remove(owner_path)
    with contextlib.suppress(OSError):
        os.rmdir(lock_dir)


def _remove_tree(path: str) -> None:
    if not os.path.lexists(path):
        return
    if os.path.islink(path) or os.path.isfile(path):
        with contextlib.suppress(FileNotFoundError):
            os.remove(path)
        return
    shutil.rmtree(path, ignore_errors=True)


# Count rendered files in a candidate figure directory, excluding the summary file
def figure_output_counts(path: str) -> dict[str, int]:
    root = pathlib.Path(path)
    counts = {"published_output_count": 0, "published_png_count": 0, "published_pdf_count": 0}
    if not root.is_dir():
        return counts
    for item in root.rglob("*"):
        if not item.is_file():
            continue
        try:
            if item.resolve() == (root / "summary.json").resolve():
                continue
        except OSError:
            if item.name == "summary.json" and item.parent == root:
                continue
        counts["published_output_count"] += 1
        suffix = item.suffix.lower()
        if suffix == ".png":
            counts["published_png_count"] += 1
        elif suffix == ".pdf":
            counts["published_pdf_count"] += 1
    return counts


def cleanup_figure_publish_outputs(*, cdir: str, run_name: str) -> None:
    """Remove stale temporary and backup directories left by interrupted publishes."""

    _remove_tree(icetune_figure_publish_tmp_dir(cdir=cdir, run_name=run_name))
    target = icetune_figure_dir(cdir=cdir, run_name=run_name)
    backups = glob(os.path.join(os.path.dirname(target), f"{run_name}.backup.*"))
    if backups and not os.path.lexists(target):
        newest = max(backups, key=os.path.getmtime)
        os.rename(newest, target)
        backups.remove(newest)
    for path in backups:
        _remove_tree(path)


# Atomically replace one published figure directory with rollback
def _commit_figure_dir(*, target_dir: str, tmp_dir: str) -> None:
    parent_dir = os.path.dirname(target_dir)
    name = os.path.basename(target_dir)
    backup_dir = os.path.join(parent_dir, f"{name}.backup.{time.time_ns()}")
    had_backup = False

    ensure_dir(parent_dir)
    if os.path.abspath(tmp_dir) == os.path.abspath(target_dir):
        return

    try:
        if os.path.lexists(target_dir):
            if os.path.islink(target_dir) or os.path.isfile(target_dir):
                os.remove(target_dir)
            else:
                os.rename(target_dir, backup_dir)
                had_backup = True
        os.rename(tmp_dir, target_dir)
    except Exception:
        if had_backup and not os.path.lexists(target_dir):
            with contextlib.suppress(OSError):
                os.rename(backup_dir, target_dir)
        raise
    else:
        if had_backup:
            _remove_tree(backup_dir)


def _initial_trial_publish_payload(payload: dict | None) -> bool:
    """Return true when a published payload represents the initial trial."""

    if not isinstance(payload, dict):
        return False
    search_payload = payload.get("search_payload")
    if isinstance(search_payload, dict) and "kind" in search_payload:
        return str(search_payload.get("kind")) == "initial"
    return False


# Copy one figure tree without recursively copying the staged initial subfolder
def _copy_figure_tree_contents(*, source_dir: str, target_dir: str, exclude_names: set[str] | None = None) -> None:
    exclude = set(exclude_names or set())
    _remove_tree(target_dir)
    ensure_dir(target_dir)
    for name in os.listdir(source_dir):
        if name in exclude:
            continue
        source = os.path.join(source_dir, name)
        target = os.path.join(target_dir, name)
        if os.path.isdir(source) and not os.path.islink(source):
            shutil.copytree(source, target, symlinks=True)
        else:
            shutil.copy2(source, target)


# Stage the initial figures inside the next canonical published tree
def _stage_initial_figure_snapshot(*, cdir: str, run_name: str, tmp_dir: str, payload: dict | None) -> None:
    target_dir = os.path.join(tmp_dir, ICETUNE_INITIAL_FIGURE_SUBDIR)
    stage_dir = os.path.join(os.path.dirname(tmp_dir), f".{os.path.basename(tmp_dir)}.init.{default_tmp_token()}")
    if _initial_trial_publish_payload(payload):
        source_dir = tmp_dir
        exclude_names = {ICETUNE_INITIAL_FIGURE_SUBDIR}
    else:
        source_dir = icetune_initial_figure_dir(cdir=cdir, run_name=run_name)
        exclude_names = set()
        if not os.path.isdir(source_dir):
            return

    try:
        _copy_figure_tree_contents(source_dir=source_dir, target_dir=stage_dir, exclude_names=exclude_names)
        _remove_tree(target_dir)
        os.rename(stage_dir, target_dir)
    finally:
        _remove_tree(stage_dir)


# Preserve the latest history plots across a selected trial figure swap
def _stage_history_plots(*, target_dir: str, tmp_dir: str) -> None:
    for name in ICETUNE_HISTORY_FIGURES:
        source = os.path.join(target_dir, name)
        target = os.path.join(tmp_dir, name)
        if os.path.isfile(source) and not os.path.exists(target):
            shutil.copy2(source, target)


# Render, validate, and timestamp one candidate figure directory
def _render_figure_candidate(*, tmp_dir: str, render_fn, run_name: str) -> dict:
    summary_file = os.path.join(tmp_dir, "summary.json")
    payload = render_fn(output_dir=tmp_dir, summary_file=summary_file)
    output_counts = figure_output_counts(tmp_dir)
    if output_counts["published_output_count"] <= 0:
        raise RuntimeError(f'Figure publication for "{run_name}" produced no figure outputs')
    if isinstance(payload, dict):
        payload.update(output_counts)
        payload = figure_summary_timestamped(payload)
        _atomic_write_json(summary_file, payload)
    return payload


# Publish one rendered figure directory under the shared per-run lock
def _publish_directory_with_lock(
    *,
    cdir: str,
    run_name: str,
    target_dir: str,
    tmp_dir: str,
    render_fn,
    before_render=None,
    before_commit=None,
    after_commit=None,
    preserve_initial: bool = False,
):
    token = acquire_publish_lock(cdir=cdir, run_name=run_name)
    published = False
    try:
        if callable(before_render):
            before_render()
        if preserve_initial:
            cleanup_figure_publish_outputs(cdir=cdir, run_name=run_name)
        else:
            _remove_tree(tmp_dir)
        ensure_dir(tmp_dir, exist_ok=False)
        payload = _render_figure_candidate(tmp_dir=tmp_dir, render_fn=render_fn, run_name=run_name)
        if preserve_initial:
            _stage_initial_figure_snapshot(cdir=cdir, run_name=run_name, tmp_dir=tmp_dir, payload=payload)
            _stage_history_plots(target_dir=target_dir, tmp_dir=tmp_dir)
        if callable(before_commit):
            before_commit(payload)
        _commit_figure_dir(target_dir=target_dir, tmp_dir=tmp_dir)
        published = True
        if callable(after_commit):
            after_commit(payload)
        return payload
    finally:
        if not published:
            _remove_tree(tmp_dir)
        release_publish_lock(cdir=cdir, run_name=run_name, token=token)


# Run one canonical figure publication transaction
def _publish_with_lock(
    *, cdir: str, run_name: str, render_fn, before_render=None, before_commit=None, after_commit=None
):
    return _publish_directory_with_lock(
        cdir=cdir,
        run_name=run_name,
        target_dir=icetune_figure_dir(cdir=cdir, run_name=run_name),
        tmp_dir=icetune_figure_publish_tmp_dir(cdir=cdir, run_name=run_name),
        render_fn=render_fn,
        before_render=before_render,
        before_commit=before_commit,
        after_commit=after_commit,
        preserve_initial=True,
    )


# Build common bootstrap initialization callbacks for one simulator driver
def bootstrap_callbacks(
    *,
    args,
    tunesetup,
    mc_steer: dict,
    simdriver,
    fingerprint_getter: Callable[[], str],
    bootstrap_builder: Callable[[], dict],
    bootstrap_stager: Callable[[dict], None] | None = None,
) -> dict:
    # Initialize one driver data bundle once unless full trial payloads are requested
    def initialize_data(pickle_dump: bool) -> None:
        if getattr(simdriver, "initialized", False) and not pickle_dump:
            return
        simdriver.init_data(
            run_name=args.run_name,
            datacards=copy.deepcopy(tunesetup.datacards),
            obs_module=args.obs_module,
            cdir=args.cdir,
            pickle_dump=pickle_dump,
        )

    # Initialize local driver data before computing the expected fingerprint
    def expected_fingerprint() -> str:
        initialize_data(False)
        return fingerprint_getter()

    # Build the immutable initialization record and stage driver outputs if needed
    def bootstrapper() -> dict:
        initialize_data(bool(args.pickle_dump))
        bootstrap = bootstrap_builder()
        if bootstrap_stager is not None:
            bootstrap_stager(bootstrap)
        return bootstrap

    # Compute the simulator baseline point published by the initialization node
    def initial_point_getter() -> dict:
        return simdriver.get_initial_param(
            param_space=tunesetup.param_space,
            aux_param_space=tunesetup.aux_param_space,
            cdir=args.cdir,
            tune_default=mc_steer["tune_default"],
        )

    return {
        "bootstrap_fingerprint_getter": expected_fingerprint,
        "bootstrapper": bootstrapper,
        "initial_point_getter": initial_point_getter,
        "initial_points": None,
    }


# Load and validate one completed Ray initialization state
def load_ray_init_state(
    *, path: str | os.PathLike[str], fingerprint: str, fingerprint_inputs: dict | None = None
) -> tuple[dict, dict | None]:
    from core.tune import init as icetune_init

    state_path = pathlib.Path(path).resolve()
    state = _read_json(state_path, default={}) or {}
    bootstrap = state.get("bootstrap") if isinstance(state, dict) else None
    if (
        not isinstance(state, dict)
        or state.get("backend") != "ray"
        or state.get("protocol") != icetune_init.PROTOCOL_NAME
        or state.get("protocol_version") != icetune_init.PROTOCOL_VERSION
        or state.get("status") != "completed"
        or not isinstance(bootstrap, dict)
    ):
        raise PermanentConfigurationError(f"Ray INIT state is incomplete: {state_path}")
    if bootstrap.get("fingerprint") != fingerprint:
        detail = f"INIT={bootstrap.get('fingerprint')}, runtime={fingerprint}"
        if fingerprint_inputs is not None:
            comparison = state_path.parent / "failures" / f"ray_init_{uuid.uuid4().hex}.json"
            try:
                _atomic_write_json(comparison, {
                    "init_fingerprint": bootstrap.get("fingerprint"), "runtime_fingerprint": fingerprint,
                    "init_inputs": bootstrap.get("fingerprint_inputs"), "runtime_inputs": fingerprint_inputs,
                })
                detail += f", compare inputs in {comparison}"
            except OSError as exc:
                detail += f", could not save {comparison}: {exc}"
        raise PermanentConfigurationError(f"Ray INIT bootstrap fingerprint differs from the staged runtime: {detail}")
    points = state.get("initial_points")
    if points is not None and not isinstance(points, dict):
        raise PermanentConfigurationError("Ray INIT initial_points is not a parameter mapping")
    return bootstrap, copy.deepcopy(points)


def cleanup_trial_outputs(*, simdriver, param: dict, tunename: str) -> None:
    """Delegate temporary trial-output cleanup to the active simulator driver."""

    cleanup = getattr(simdriver, "cleanup_trial_outputs", None)
    try:
        if callable(cleanup):
            cleanup(tunename=tunename, datacards=param["datacards"], mc_steer=param["mc_steer"], cdir=param["cdir"])
    except Exception as cleanup_error:
        logger.warning('.cleanup_trial_outputs: cleanup failed for "%s" (%s)', tunename, cleanup_error)


def _normalize_trial_outputs(*, outputs: dict, config: dict, param: dict, trial_id: str, tunename: str) -> dict:
    if not isinstance(outputs, dict):
        raise TypeError("Driver evaluate_trial_outputs() must return a dictionary")

    metrics = outputs.get("metrics")
    if not isinstance(metrics, dict):
        raise TypeError('Trial outputs are missing a "metrics" dictionary')
    if param["cost"] not in metrics:
        raise KeyError(f'Trial outputs are missing the configured cost metric "{param["cost"]}"')

    normalized_config = copy.deepcopy(outputs.get("config", config))
    card_config = copy.deepcopy(outputs.get("card_config"))
    if not isinstance(card_config, dict):
        card_config = {"schema_version": 1, "parameters": copy.deepcopy(normalized_config), "tables": {}}
    theta = copy.deepcopy(outputs.get("theta"))
    theta_hash_value = copy.deepcopy(outputs.get("theta_hash"))
    if theta is None or theta_hash_value is None:
        try:
            theta, theta_hash_value = build_theta_record(normalized_config, param.get("parameter_space", []))
        except Exception:
            theta, theta_hash_value = None, None

    likelihood = copy.deepcopy(outputs.get("likelihood"))
    if not isinstance(likelihood, dict):
        likelihood = icetune_likelihood.build_metric_likelihood_payload(
            metrics=metrics, cost_key=param["cost"], error=outputs.get("error")
        )

    return {
        **outputs,
        "card_config": card_config,
        "config": normalized_config,
        "cost_definition": optimizer_cost_definition(param["cost"]),
        "error": outputs.get("error"),
        "likelihood": likelihood,
        "results": outputs.get("results"),
        **selected_fields(
            outputs,
            ("cost_arr", "metrics", "ndf_arr", "search_payload", "valid_arr", "weight_arr"),
            include_missing=True,
        ),
        "theta": theta,
        "theta_hash": theta_hash_value,
        "trial_id": outputs.get("trial_id", trial_id),
        "tunename": outputs.get("tunename", tunename),
    }


def evaluate_trial_outputs(*, config: dict, param: dict, simdriver, trial_id: str, tunename: str | None = None) -> dict:
    """Evaluate one trial through the active driver and return a normalized output bundle."""

    logger.info("[Trial ID: %s] Starting trial step with config = %s", trial_id, config)
    tunename = tunename or simdriver.trial_tunename(trial_id=trial_id, node_id=param.get("node_id"))

    evaluator = getattr(simdriver, "evaluate_trial_outputs", None)
    if not callable(evaluator):
        raise NotImplementedError(f"{simdriver.__class__.__name__} does not implement evaluate_trial_outputs()")

    outputs = evaluator(config=config, param=param, trial_id=trial_id, tunename=tunename)
    return _normalize_trial_outputs(outputs=outputs, config=config, param=param, trial_id=trial_id, tunename=tunename)


# Evaluate one trial and guarantee cleanup of its driver-owned temporary outputs
@contextlib.contextmanager
def trial_evaluation(*, config: dict, param: dict, simdriver, trial_id: str, tunename: str | None = None):
    outputs = evaluate_trial_outputs(
        config=config, param=param, simdriver=simdriver, trial_id=trial_id, tunename=tunename
    )
    try:
        yield outputs
    finally:
        cleanup_trial_outputs(simdriver=simdriver, param=param, tunename=outputs["tunename"])


# Compute the conventional or explicitly configured full-pickle directory
def trial_pickle_dir(param: dict) -> pathlib.Path:
    configured = param.get("pickle_dir")
    if configured:
        return pathlib.Path(configured)
    return pathlib.Path(param["cdir"], "runs", "icetune", param["run_name"], "results")


# Compute the deterministic iceproxy-compatible filename for one trial
def trial_pickle_filename(outputs: dict) -> str:
    trial_id = outputs.get("trial_id")
    if not trial_id:
        raise ValueError("Full trial pickle requires a non-empty trial_id")
    return f"TUNE_icetune_{_sanitize_tunename_component(trial_id)}.pkl"


# Build the full observable payload consumed by iceproxy and icescape
def _trial_pickle_payload(*, outputs: dict, param: dict) -> dict:
    return {
        "replica_schema_version": 1,
        "param": copy.deepcopy(param),
        **selected_fields(outputs, TRIAL_PICKLE_FIELDS, include_missing=True),
        **time_fields("created_at", time.time()),
    }


# Require an existing immutable trial pickle to represent the same parameter point
def _validate_trial_pickle_identity(path: pathlib.Path, payload: dict) -> None:
    with open(path, "rb") as source:
        existing = pickle.load(source)
    if not isinstance(existing, dict):
        raise RuntimeError(f'Invalid full trial pickle payload at "{path}"')
    for key in ("replica_schema_version", "trial_id", "card_config", "config", "theta_hash"):
        if existing.get(key) != payload.get(key):
            raise RuntimeError(f'Conflicting full trial pickle at "{path}" for identity field "{key}"')


# Serialize locally and atomically publish one immutable full trial pickle
def _publish_trial_pickle(*, payload: dict, param: dict, destination: pathlib.Path, verify: bool) -> dict:
    stage_dir = pathlib.Path(param["cdir"], "tmp")
    ensure_dir(stage_dir)
    ensure_dir(destination.parent)
    stage_path = stage_dir / f".{destination.name}.stage.{uuid.uuid4().hex}.tmp"
    try:
        with open(stage_path, "wb") as output:
            pickle.dump(payload, output, protocol=pickle.HIGHEST_PROTOCOL)
            output.flush()
            os.fsync(output.fileno())
        stage_checksum = sha256_file(stage_path)
        stage_size = stage_path.stat().st_size
        installed = publish_file_immutable(stage_path, str(destination))
        if not installed:
            _validate_trial_pickle_identity(destination, payload)
        destination_checksum = sha256_file(destination) if verify or not installed else stage_checksum
        if destination.stat().st_size != stage_size or destination_checksum != stage_checksum:
            raise RuntimeError(f'Full trial pickle checksum mismatch at "{destination}"')
        return {"filename": destination.name, "schema_version": 1, "sha256": destination_checksum, "size": stage_size}
    finally:
        stage_path.unlink(missing_ok=True)


# Persist a full trial payload when pickle dumping is enabled
def maybe_dump_trial_payload(
    *, outputs: dict, param: dict, destination: str | os.PathLike[str] | None = None
) -> dict | None:
    if not param["pickle_dump"]:
        return None
    if outputs.get("results") is None:
        raise RuntimeError(f'Full trial pickle requires observable results for "{outputs.get("trial_id")}"')
    payload = _trial_pickle_payload(outputs=outputs, param=param)
    target = (
        trial_pickle_dir(param) / trial_pickle_filename(outputs) if destination is None else pathlib.Path(destination)
    )
    return _publish_trial_pickle(payload=payload, param=param, destination=target, verify=destination is None)


def build_trial_summary_payload(
    *,
    outputs: dict,
    param: dict,
    extra: dict | None = None,
    trial_id: str | None = None,
    tunename: str | None = None,
    metrics: dict | None = None,
) -> dict:
    """Build the metadata payload written into a published figure summary."""

    selected_metrics = copy.deepcopy(metrics if metrics is not None else outputs["metrics"])
    config = copy.deepcopy(outputs["config"])
    payload = {
        "best_fit": fit_summary.build_best_fit(
            values=config,
            objective_name=param["cost"],
            objective_value=selected_metrics.get(param["cost"]),
            source="realized",
        ),
        "card_config": copy.deepcopy(outputs.get("card_config")),
        "config": config,
        "cost": param["cost"],
        "cost_definition": optimizer_cost_definition(param["cost"]),
        **time_fields("created_at", time.time()),
        "likelihood": copy.deepcopy(outputs.get("likelihood")),
        "metrics": selected_metrics,
        "search_payload": copy.deepcopy(outputs.get("search_payload")),
        "theta": copy.deepcopy(outputs.get("theta")),
        "theta_hash": copy.deepcopy(outputs.get("theta_hash")),
        "trial_id": trial_id if trial_id is not None else outputs["trial_id"],
        "tunename": tunename if tunename is not None else outputs["tunename"],
    }
    if param.get("simdriver") is not None:
        payload["simdriver"] = copy.deepcopy(param["simdriver"])
    if param.get("plot_brand") is not None:
        payload["plot_brand"] = copy.deepcopy(param["plot_brand"])
    if extra:
        payload.update(copy.deepcopy(extra))
    return payload


def publish_trial_figures(
    *,
    outputs: dict,
    param: dict,
    simdriver,
    summary_payload: dict,
    before_render=None,
    before_commit=None,
    after_commit=None,
) -> dict:
    """Render a trial's figures and publish them via the canonical directory swap."""

    return _publish_with_lock(
        cdir=param["cdir"],
        run_name=param["run_name"],
        before_render=before_render,
        before_commit=before_commit,
        after_commit=after_commit,
        # Render one trial's figure tree into a caller-provided directory
        render_fn=partial(
            simdriver.render_trial_figures_to_dir, outputs=outputs, param=param, summary_payload=summary_payload
        ),
    )


# Publish an initial trial directly into figs/icetune/<run>/init
def publish_initial_trial_figures(*, outputs: dict, param: dict, simdriver, summary_payload: dict) -> dict:
    target_dir = icetune_initial_figure_dir(cdir=param["cdir"], run_name=param["run_name"])
    return _publish_directory_with_lock(
        cdir=param["cdir"],
        run_name=param["run_name"],
        target_dir=target_dir,
        tmp_dir=str(default_tmp_stage_path(target_dir)),
        render_fn=partial(
            simdriver.render_trial_figures_to_dir, outputs=outputs, param=param, summary_payload=summary_payload
        ),
    )


# Evaluate, publish, and clean one plot-only rerun through a selected publisher
def _publish_trial_rerun(
    *,
    config: dict,
    param: dict,
    simdriver,
    trial_id: str,
    publisher: Callable[..., dict],
    metrics: dict | None = None,
    summary_tunename: str | None = None,
    search_payload: dict | Callable[[dict], dict] | None = None,
    extra: dict | Callable[[dict], dict] | None = None,
    after_publish: Callable[[dict, dict], None] | None = None,
) -> dict | None:
    tunename = simdriver.trial_tunename(
        trial_id=trial_id, node_id=param.get("node_id"), pid=os.getpid(), purpose="publish"
    )
    with trial_evaluation(
        config=config, param=param, simdriver=simdriver, trial_id=trial_id, tunename=tunename
    ) as outputs:
        if outputs.get("results") is None:
            return None
        if search_payload is not None:
            selected = search_payload(outputs) if callable(search_payload) else search_payload
            outputs["search_payload"] = copy.deepcopy(selected)
        summary_extra = extra(outputs) if callable(extra) else extra
        summary = build_trial_summary_payload(
            outputs=outputs,
            param=param,
            trial_id=trial_id,
            tunename=summary_tunename,
            metrics=metrics,
            extra=summary_extra,
        )
        published = publisher(outputs=outputs, param=param, simdriver=simdriver, summary_payload=summary)
        if callable(after_publish):
            after_publish(summary, published)
        return published


# Evaluate and publish one canonical plot-only rerun
def publish_trial_rerun(
    *,
    config: dict,
    param: dict,
    simdriver,
    trial_id: str,
    metrics: dict | None = None,
    summary_tunename: str | None = None,
    search_payload: dict | Callable[[dict], dict] | None = None,
    extra: dict | Callable[[dict], dict] | None = None,
    after_publish: Callable[[dict, dict], None] | None = None,
) -> dict | None:
    return _publish_trial_rerun(
        config=config,
        param=param,
        simdriver=simdriver,
        trial_id=trial_id,
        publisher=publish_trial_figures,
        metrics=metrics,
        summary_tunename=summary_tunename,
        search_payload=search_payload,
        extra=extra,
        after_publish=after_publish,
    )

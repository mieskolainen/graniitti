# Read compact and Ray Tune icetune histories
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

"""Print and plot completed trial costs from an icetune run."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import re
import tempfile
import uuid
from dataclasses import dataclass
from datetime import UTC, datetime

from core.io.files import ensure_dir
from core.io.serialize import load_json_file
from core.tune.cache import json_fingerprint

HISTORY_SCHEMA_VERSION = 1
RANDOM_PROPOSAL_TYPES = frozenset({"basic", "cold", "initial", "random"})
RANDOM_BUDGET_PROPOSAL_TYPES = frozenset({"bayesopt", "hebo", "icebo"})
SURROGATE_PROPOSAL_TYPES = frozenset({"ax", "bayesopt", "hebo", "hyperopt", "icebo", "optuna", "scikit"})
RAY_CAMPAIGN_KEYS = ("ALGORITHM", "COST", "PLOT_BRAND", "RAND_TRIALS", "SIMDRIVER")


# Store one completed trial used by the table and cost plot
@dataclass(frozen=True)
class CostPoint:
    trial_id: str
    cost: float
    completed_at: str
    completed_unix: float
    node_id: str
    proposal_type: str | None = None
    acquisition: str | None = None
    started_unix: float | None = None
    config: dict | None = None
    card_config: dict | None = None
    start: int | None = None


# Create the standalone parser used by the icetune history command
def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="icetune history", description=__doc__)
    parser.add_argument("source", type=pathlib.Path, help="run name, run directory, or history.json")
    parser.add_argument("--cdir", type=pathlib.Path, default=pathlib.Path.cwd(), help="main run directory")
    parser.add_argument("--cost", help="Ray cost field override")
    parser.add_argument("--sort", action="store_true", help="put the best trial last")
    parser.add_argument("--run-name", help="figure directory name, inferred from the history path by default")
    return parser


# Compute Ray trial result files below one run directory
def _ray_result_paths(path: pathlib.Path) -> list[pathlib.Path]:
    if path.is_file() and path.name == "result.json":
        return [path]
    if not path.is_dir():
        return []
    return sorted(candidate for candidate in path.rglob("result.json") if candidate.is_file())


# Select the most recently updated compact campaign history below one run directory
def _nested_history_path(path: pathlib.Path) -> pathlib.Path | None:
    if not path.is_dir():
        return None
    histories = [candidate for candidate in (path / "campaigns").glob("*/history.json") if candidate.is_file()]
    if not histories:
        return None
    return max(histories, key=lambda candidate: (candidate.stat().st_mtime_ns, str(candidate)))


# Resolve a run name or local path to one compact or Ray Tune history source
def resolve_history_path(source: pathlib.Path, cdir: pathlib.Path) -> pathlib.Path:
    root = cdir.expanduser().resolve()
    requested = source.expanduser()
    if not requested.is_absolute():
        requested = root / requested
    requested = requested.resolve()

    candidates = []
    if requested.is_dir():
        candidates.append(requested / "history.json")
    candidates.append(requested)
    if not source.is_absolute() and len(source.parts) == 1:
        run_dir = root / "runs" / "icetune" / source
        candidates.extend((run_dir / "history.json", run_dir))

    candidates = list(dict.fromkeys(candidates))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
        nested = _nested_history_path(candidate)
        if nested is not None:
            return nested
        if _ray_result_paths(candidate):
            return candidate

    paths = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"No compact or Ray Tune history found. Checked: {paths}")


# Build a natural ordering key for trial identifiers
def trial_id_key(trial_id: str) -> tuple[tuple[int, int | str], ...]:
    return tuple(
        (0, int(part)) if part.isascii() and part.isdigit() else (1, part) for part in re.split(r"([0-9]+)", trial_id)
    )


# Require one nonempty string field from a trial record
def _trial_text(trial: dict, trial_id: str, field: str) -> str:
    value = trial.get(field)
    if not isinstance(value, str) or not value.strip() or any(ord(char) < 32 for char in value):
        raise ValueError(f'trial "{trial_id}" is missing the string field "{field}"')
    return value


# Convert one required numeric field to a finite float
def _finite_number(value, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} is not numeric")
    try:
        number = float(value)
    except OverflowError:
        raise ValueError(f"{label} is nonfinite") from None
    if not math.isfinite(number):
        raise ValueError(f"{label} is nonfinite")
    return number


# Read the completion time used for chronological cost evolution
def _completion_time(trial: dict, trial_id: str, completed_at: str) -> float:
    value = trial.get("completed_at_unix")
    if value is not None:
        return _finite_number(value, f'trial "{trial_id}" "completed_at_unix"')

    normalized = completed_at.strip()
    if normalized.endswith("Z"):
        normalized = normalized[:-1] + "+00:00"
    try:
        parsed = datetime.fromisoformat(normalized)
    except ValueError:
        match = re.fullmatch(r"(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) [A-Z]{2,6}", normalized)
        if match is None:
            raise ValueError(f'trial "{trial_id}" has an invalid completion timestamp') from None
        parsed = datetime.strptime(match.group(1), "%Y-%m-%d %H:%M:%S").replace(tzinfo=UTC)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=UTC)
    return parsed.timestamp()


# Read an optional finite trial start time
def _started_time(trial: dict, trial_id: str) -> float | None:
    value = trial.get("started_at_unix")
    if value is None:
        return None
    return _finite_number(value, f'trial "{trial_id}" "started_at_unix"')


# Read optional proposal metadata used to locate the surrogate transition
def _proposal_metadata(trial: dict) -> tuple[str | None, str | None, int | None]:
    payload = trial.get("search_payload")
    if not isinstance(payload, dict):
        return None, None, None
    proposal = payload.get("kind")
    acquisition = payload.get("acquisition")
    proposal_type = proposal.strip().lower() if isinstance(proposal, str) else None
    acquisition_type = acquisition.strip().lower() if isinstance(acquisition, str) else None
    start = payload.get("start")
    if start is not None and (type(start) is not int or start < 0):
        raise ValueError("Trial search_payload.start must be a nonnegative integer")
    return proposal_type or None, acquisition_type or None, start


# Load and validate completed trial costs from compact schema version 1
def _load_compact_history(history_path: pathlib.Path) -> tuple[dict, str, list[CostPoint]]:
    with history_path.open(encoding="utf-8") as handle:
        history = load_json_file(handle.name)
    if not isinstance(history, dict):
        raise ValueError("history root must be an object")

    schema_version = history.get("history_schema_version")
    if type(schema_version) is not int or schema_version != HISTORY_SCHEMA_VERSION:
        raise ValueError(f'"history_schema_version" must be {HISTORY_SCHEMA_VERSION}, got {schema_version!r}')

    optimization = history.get("optimization")
    cost_key = optimization.get("cost") if isinstance(optimization, dict) else None
    if not isinstance(cost_key, str) or not cost_key.strip() or any(ord(character) < 32 for character in cost_key):
        raise ValueError('history is missing the string field "optimization.cost"')

    trials = history.get("trials")
    if not isinstance(trials, list):
        raise ValueError('history is missing the array field "trials"')
    failed = history.get("failed", [])
    if not isinstance(failed, list):
        raise ValueError('history field "failed" is not an array')

    points = []
    trial_ids = set()
    for index, trial in enumerate(trials):
        if not isinstance(trial, dict):
            raise ValueError(f'"trials[{index}]" must be an object')
        trial_id = _trial_text(trial, f"trials[{index}]", "trial_id")
        if trial_id in trial_ids:
            raise ValueError(f'duplicate trial identifier "{trial_id}"')
        trial_ids.add(trial_id)

        metrics = trial.get("metrics")
        value = metrics.get(cost_key) if isinstance(metrics, dict) else None
        cost = _finite_number(value, f'trial "{trial_id}" "{cost_key}" cost')

        completed_at = _trial_text(trial, trial_id, "completed_at_datetime")
        completed_unix = _completion_time(trial, trial_id, completed_at)
        started_unix = _started_time(trial, trial_id)
        node_id = _trial_text(trial, trial_id, "node_id")
        proposal_type, acquisition, start = _proposal_metadata(trial)
        points.append(
            CostPoint(
                trial_id=trial_id,
                cost=cost,
                completed_at=completed_at,
                completed_unix=completed_unix,
                node_id=node_id,
                proposal_type=proposal_type,
                acquisition=acquisition,
                start=start,
                started_unix=started_unix,
                config=trial.get("config") if isinstance(trial.get("config"), dict) else None,
                card_config=trial.get("card_config") if isinstance(trial.get("card_config"), dict) else None,
            )
        )

    return history, cost_key, sorted(points, key=lambda point: trial_id_key(point.trial_id))


# Read complete JSON objects from one active or finished Ray result stream
def _ray_result_rows(path: pathlib.Path) -> list[dict]:
    lines = path.read_text(encoding="utf-8").splitlines()
    rows = []
    for index, line in enumerate(lines):
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            if index == len(lines) - 1:
                continue
            raise ValueError(f'Ray result file "{path}" contains invalid JSON') from None
        if not isinstance(row, dict):
            raise ValueError(f'Ray result file "{path}" contains a non-object row')
        rows.append(row)
    return rows


# Read campaign metadata from the immutable Ray campaign identity
def _ray_identity_environment(path: pathlib.Path) -> dict:
    payload = load_json_file(path)
    identity = payload.get("identity") if isinstance(payload, dict) else None
    if not isinstance(identity, dict):
        raise ValueError(f'Ray campaign identity "{path}" has no identity object')
    if payload.get("schema_version") != 1 or payload.get("fingerprint") != json_fingerprint(identity):
        raise ValueError(f'Ray campaign identity "{path}" has an invalid fingerprint')
    optimization = identity.get("optimization")
    physics = identity.get("physics")
    if not isinstance(optimization, dict) or not isinstance(physics, dict):
        raise ValueError(f'Ray campaign identity "{path}" is incomplete')
    environment = {
        "ALGORITHM": identity.get("algorithm"),
        "COST": optimization.get("cost"),
        "PLOT_BRAND": physics.get("plot_brand"),
        "RAND_TRIALS": optimization.get("rand_trials"),
        "SIMDRIVER": physics.get("simdriver"),
    }
    missing = [key for key in RAY_CAMPAIGN_KEYS if key != "PLOT_BRAND" and environment[key] is None]
    if missing:
        raise ValueError(f'Ray campaign identity "{path}" is missing {", ".join(missing)}')
    return environment


# Read immutable campaign metadata when it accompanies a Ray experiment
def _ray_campaign_environment(run_dir: pathlib.Path) -> dict:
    identity_path = run_dir / "icetune_campaign.json"
    if identity_path.is_file():
        return _ray_identity_environment(identity_path)
    return {}


# Infer the Ray Tune cost field without guessing among unrelated metrics
def _ray_cost_key(rows: list[dict], environment: dict, override: str | None) -> str:
    if override is not None:
        cost = str(override).strip()
        if not cost:
            raise ValueError("Ray cost override cannot be empty")
        return cost
    configured = environment.get("COST")
    if isinstance(configured, str) and configured.strip():
        return configured.strip()
    definitions = {
        row.get("cost_definition")
        for row in rows
        if isinstance(row.get("cost_definition"), str) and row.get("cost_definition").strip()
    }
    if len(definitions) == 1:
        definition = definitions.pop()
        if definition == "piecewise_constant_density_cdf_l1_v2":
            return "wasserstein"
        ratio2 = {"symmetric_log_ratio_chi2_v2", "symmetric_asinh_ratio_chi2_v3"}
        return "ratio2" if definition in ratio2 else definition
    known = [key for key in ("gaussian", "ratio2", "chi2", "wasserstein", "loss", "objective") if any(key in row for row in rows)]
    if len(known) == 1:
        return known[0]
    raise ValueError("Ray cost field is ambiguous. Pass --cost COST")


# Compute one finite Ray cost or none for a nonterminal row
def _ray_cost(row: dict, cost_key: str) -> float | None:
    value = row.get(cost_key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    try:
        cost = float(value)
    except OverflowError:
        return None
    return cost if math.isfinite(cost) else None


# Infer a Ray trial start time from its elapsed runtime when available
def _ray_started_time(row: dict, completed_unix: float) -> float | None:
    elapsed = _ray_cost(row, "time_total_s")
    return completed_unix - elapsed if elapsed is not None and elapsed >= 0.0 else None


# Load completed costs from Ray Tune line delimited trial results
def _load_ray_history(run_dir: pathlib.Path, cost_override: str | None) -> tuple[dict, str, list[CostPoint]]:
    streams = [(path, _ray_result_rows(path)) for path in _ray_result_paths(run_dir)]
    rows = [row for _, stream in streams for row in stream]
    if not rows:
        raise ValueError(f'Ray run "{run_dir}" has no complete result rows')
    metadata_dir = run_dir if run_dir.is_dir() else run_dir.parent.parent
    environment = _ray_campaign_environment(metadata_dir)
    cost_key = _ray_cost_key(rows, environment, cost_override)
    if not any(_ray_cost(row, cost_key) is not None for row in rows):
        raise ValueError(f'Ray results contain no finite "{cost_key}" cost')
    points = []
    trial_ids = set()
    parameter_names = set()
    for path, stream in streams:
        completed = [(row, cost) for row in stream if (cost := _ray_cost(row, cost_key)) is not None]
        if not completed:
            continue
        row, cost = completed[-1]
        trial_id = row.get("trial_id") or row.get("trial_name") or path.parent.name
        if not isinstance(trial_id, str) or not trial_id.strip():
            raise ValueError(f'Ray result file "{path}" has no trial identifier')
        if trial_id in trial_ids:
            raise ValueError(f'duplicate Ray trial identifier "{trial_id}"')
        trial_ids.add(trial_id)
        completed_unix = _finite_number(row.get("timestamp"), f'Ray trial "{trial_id}" completion timestamp')
        completed = datetime.fromtimestamp(completed_unix, tz=UTC).astimezone()
        completed_at = completed.strftime("%Y-%m-%d %H:%M:%S %Z")
        node_id = row.get("hostname") or row.get("node_ip") or "ray"
        config = row.get("config")
        if isinstance(config, dict):
            parameter_names.update(str(name) for name in config)
        card_config = row.get("card_config")
        proposal_type, acquisition, start = _proposal_metadata(row)
        points.append(
            CostPoint(
                trial_id=trial_id,
                cost=cost,
                completed_at=completed_at,
                completed_unix=completed_unix,
                node_id=str(node_id),
                proposal_type=proposal_type,
                acquisition=acquisition,
                start=start,
                started_unix=_ray_started_time(row, completed_unix),
                config=config if isinstance(config, dict) else None,
                card_config=card_config if isinstance(card_config, dict) else None,
            )
        )
    optimization = {"backend": "ray", "cost": cost_key}
    optimizer = environment.get("ALGORITHM")
    if isinstance(optimizer, str) and optimizer.strip():
        optimization["optimizer"] = optimizer.strip().lower()
    raw_rand_trials = environment.get("RAND_TRIALS")
    try:
        rand_trials = int(raw_rand_trials)
    except (TypeError, ValueError, OverflowError):
        rand_trials = -1
    if not isinstance(raw_rand_trials, bool) and rand_trials >= 0:
        optimization["rand_trials"] = rand_trials
    history = {
        "history_schema_version": HISTORY_SCHEMA_VERSION,
        "optimization": optimization,
        "parameter_space": {name: {} for name in sorted(parameter_names)},
        "plot_brand": environment.get("PLOT_BRAND"),
        "simdriver": environment.get("SIMDRIVER"),
    }
    return history, cost_key, sorted(points, key=lambda point: trial_id_key(point.trial_id))


# Load one resolved compact or Ray Tune history source
def load_history(history_path: pathlib.Path, cost_override: str | None = None) -> tuple[dict, str, list[CostPoint]]:
    if history_path.is_dir() or history_path.name == "result.json":
        return _load_ray_history(history_path, cost_override)
    return _load_compact_history(history_path)


# Format validated history rows as one aligned table
def format_cost_values(cost_key: str, points: list[CostPoint]) -> str:
    headers = ("trial_id", cost_key, "completed_at_datetime", "node_id")
    values = [(point.trial_id, f"{point.cost:.4f}", point.completed_at, point.node_id) for point in points]
    widths = [max([len(headers[index]), *(len(row[index]) for row in values)]) for index in range(len(headers))]
    lines = [
        " | ".join(value.ljust(width) for value, width in zip(headers, widths, strict=True)),
        "-+-".join("-" * width for width in widths),
    ]
    for row in values:
        columns = [row[0].ljust(widths[0]), row[1].rjust(widths[1]), row[2].ljust(widths[2]), row[3].ljust(widths[3])]
        lines.append(" | ".join(columns))
    return "\n".join(lines)


# Compute chronological costs with running mean, spread and minimum
def cost_evolution(points: list[CostPoint]):
    ordered = sorted(points, key=lambda point: (point.completed_unix, trial_id_key(point.trial_id)))
    if not ordered:
        return [], [], [], [], [], []
    origin = ordered[0].completed_unix
    elapsed = [point.completed_unix - origin for point in ordered]
    costs = [point.cost for point in ordered]
    average = []
    sigma = []
    best = []
    mean = 0.0
    squared_deviation = 0.0
    minimum = math.inf
    for count, cost in enumerate(costs, start=1):
        delta = cost - mean
        mean += delta / count
        squared_deviation += delta * (cost - mean)
        average.append(mean)
        sigma.append(math.sqrt(max(0.0, squared_deviation / count)))
        minimum = min(minimum, cost)
        best.append(minimum)
    return ordered, elapsed, costs, average, sigma, best


# Sort trial costs from highest to lowest with the minimum last
def sorted_costs(points: list[CostPoint]) -> list[CostPoint]:
    return sorted(points, key=lambda point: (-point.cost, trial_id_key(point.trial_id)))


# Keep finite numeric parameter values from one trial mapping
def _numeric_parameters(values: dict | None) -> dict[str, float]:
    if not isinstance(values, dict):
        return {}
    output = {}
    for name, value in values.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            continue
        number = float(value)
        if math.isfinite(number):
            output[str(name)] = number
    return output


# Locate the initial trial independently of completion order
def _initial_index(ordered: list[CostPoint]) -> int | None:
    initial = next((index for index, point in enumerate(ordered) if point.proposal_type == "initial"), None)
    if initial is not None:
        return initial
    return next((index for index, point in enumerate(ordered) if re.fullmatch(r"(?:trial-)?0+", point.trial_id)), None)


# Compute chronological raw or physical parameter series from completed trials
def parameter_evolution(
    history: dict, points: list[CostPoint], *, physical: bool
) -> tuple[list[CostPoint], list[str], dict[str, list[float]]]:
    ordered = sorted(points, key=lambda point: (point.completed_unix, trial_id_key(point.trial_id)))
    if physical:
        mappings = [
            _numeric_parameters(point.card_config.get("parameters"))
            if point.card_config is not None else {}
            for point in ordered
        ]
        names = sorted({name for values in mappings for name in values})
    else:
        mappings = [_numeric_parameters(point.config) for point in ordered]
        parameter_space = history.get("parameter_space")
        declared = []
        if isinstance(parameter_space, list):
            declared = [str(item["name"]) for item in parameter_space if isinstance(item, dict) and "name" in item]
        elif isinstance(parameter_space, dict):
            declared = [str(name) for name in parameter_space]
        available = {name for values in mappings for name in values}
        names = [name for name in declared if name in available]
        names.extend(sorted(available - set(names)))
    series = {name: [values.get(name, math.nan) for values in mappings] for name in names}
    return ordered, names, series


# Compute raw optimizer bounds keyed by parameter name
def _parameter_bounds(history: dict) -> dict[str, tuple[float, float]]:
    parameter_space = history.get("parameter_space")
    if not isinstance(parameter_space, list):
        return {}
    bounds = {}
    for item in parameter_space:
        if not isinstance(item, dict) or "name" not in item:
            continue
        try:
            lower = float(item["lower"])
            upper = float(item["upper"])
        except (KeyError, TypeError, ValueError, OverflowError):
            continue
        if math.isfinite(lower) and math.isfinite(upper):
            bounds[str(item["name"])] = (lower, upper)
    return bounds


# Identify local optimization from the run or recorded proposals
def _local_optimizer(history: dict, ordered: list[CostPoint]) -> bool:
    optimization = history.get("optimization") or {}
    return optimization.get("optimizer") in {"ampfit", "lbfgs"} or any(
        point.proposal_type in {"ampfit", "lbfgs"} or point.start is not None for point in ordered
    )


# Group ampfit evaluations by their recorded start with stable colors
def start_trajectories(ordered: list[CostPoint], history: dict | None = None) -> list[tuple[str, object, list[int]]]:
    import matplotlib as mpl

    starts = sorted({point.start for point in ordered if point.start is not None})
    if not starts:
        return []
    count = _start_count(history or {}, ordered)
    groups = [
        (f"Start {start}", mpl.colormaps["turbo"](start / (count - 1)) if count > 10 else
         mpl.colormaps["tab20"](2 * start),
         [index for index, point in enumerate(ordered) if point.start == start])
        for start in starts
    ]
    missing = [index for index, point in enumerate(ordered) if point.start is None]
    if missing:
        groups.append(("Unassigned", "#777777", missing))
    return groups


# Preserve the start color scale across partial histories of the same campaign
def _start_count(history, ordered):
    settings = (history.get("optimization") or {}).get("ampfit") or {}
    return max(settings.get("starts", 0), max((point.start + 1 for point in ordered if point.start is not None), default=0))


# Use a start color scale when individual legend entries would crowd the plot
def _start_colorbar(ax, history, ordered):
    import matplotlib as mpl

    count = _start_count(history, ordered)
    if count > 10:
        scale = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, count - 1), cmap="turbo")
        bar = ax.figure.colorbar(scale, ax=ax, pad=0.02, fraction=0.025)
        bar.set_label("Start", fontsize=7)
        bar.ax.tick_params(labelsize=7)
        return bar


# Draw evaluations within each start without connecting independent trajectories
def _draw_trajectories(ax, groups, x_values, values, *, connect=True, legend=True):
    for label, color, indices in groups:
        ax.plot(
            [x_values[index] for index in indices], [values[index] for index in indices],
            color=color, label=label if legend else "_nolegend_", linewidth=0.8 if connect and label != "Unassigned" else 0,
            marker="o", markersize=3, markeredgewidth=0.3, markeredgecolor="white", zorder=2,
        )


# Classify completed trials into exploration and optimizer proposal stages
def trial_stages(history: dict, ordered: list[CostPoint]) -> list[bool]:
    if _local_optimizer(history, ordered):
        return [True] * len(ordered)
    optimization = history.get("optimization")
    optimization = optimization if isinstance(optimization, dict) else {}
    rand_trials = optimization.get("rand_trials")
    has_budget = not isinstance(rand_trials, bool) and isinstance(rand_trials, int) and rand_trials >= 0
    stages = []
    for index, point in enumerate(ordered):
        surrogate = _uses_surrogate(point)
        if surrogate is not None:
            stages.append(surrogate)
        elif point.proposal_type is not None:
            stages.append(True)
        elif has_budget:
            stages.append(index >= rand_trials)
        else:
            stages.append(True)
    return stages


# Select readable elapsed time units for the horizontal axis
def _time_axis(elapsed: list[float]) -> tuple[float, str]:
    duration = max(elapsed, default=0.0)
    if duration >= 2.0 * 86400.0:
        return 86400.0, "d"
    if duration >= 2.0 * 3600.0:
        return 3600.0, "h"
    if duration >= 120.0:
        return 60.0, "min"
    return 1.0, "s"


# Compute whether one proposal was made before or after surrogate fitting
def _uses_surrogate(point: CostPoint) -> bool | None:
    if point.proposal_type in RANDOM_PROPOSAL_TYPES:
        return False
    if point.proposal_type == "icebo" and point.acquisition == "sobol_warmup":
        return False
    if point.proposal_type in SURROGATE_PROPOSAL_TYPES:
        return True
    return None


# Locate the last random completion on the elapsed completion time axis
def random_regime_end(history: dict, ordered: list[CostPoint]) -> float | None:
    if not ordered or _local_optimizer(history, ordered):
        return None
    origin = ordered[0].completed_unix
    states = [_uses_surrogate(point) for point in ordered]
    random_seen = any(state is False for state in states)
    if random_seen and any(state is True for state in states):
        # Random and surrogate trials can finish in either order with parallel workers
        transition = max(point.completed_unix for point, state in zip(ordered, states, strict=True) if state is False)
        return max(0.0, transition - origin)
    if random_seen and all(state is not True for state in states):
        return math.inf
    if any(state is not None for state in states):
        return None

    optimization = history.get("optimization")
    if not isinstance(optimization, dict):
        return None
    optimizer = optimization.get("optimizer")
    optimizer = optimizer.strip().lower() if isinstance(optimizer, str) else None
    rand_trials = optimization.get("rand_trials")
    if (
        optimizer not in RANDOM_BUDGET_PROPOSAL_TYPES
        or isinstance(rand_trials, bool)
        or not isinstance(rand_trials, int)
        or rand_trials <= 0
    ):
        return None
    if len(ordered) <= rand_trials:
        return math.inf
    transition = ordered[rand_trials - 1].completed_unix
    return max(0.0, transition - origin)


# Extract the run name and campaign identifier from one nested history path
def campaign_parts(history_path: pathlib.Path) -> tuple[str, str] | None:
    parts = history_path.expanduser().resolve().parts
    for index in range(len(parts) - 4):
        if parts[index : index + 2] != ("runs", "icetune"):
            continue
        if parts[index + 3] == "campaigns":
            return parts[index + 2], parts[index + 4]
    return None


# Resolve and validate the directory name used for history figures
def resolve_run_name(history_path: pathlib.Path, override: str | None = None) -> str:
    name = str(override).strip() if override is not None else ""
    if not name:
        campaign = campaign_parts(history_path)
        if campaign is not None:
            name = campaign[0]
        elif history_path.is_dir():
            name = history_path.name
        elif history_path.name == "result.json":
            name = history_path.parent.parent.name
        else:
            name = history_path.parent.name
    if (
        not name
        or name in {".", ".."}
        or pathlib.Path(name).name != name
        or any(ord(character) < 32 for character in name)
    ):
        raise ValueError(f'Invalid history run name "{name}"')
    return name


# Resolve the canonical figure directory for one flat or nested history source
def resolve_plot_dir(*, history_path: pathlib.Path, cdir: pathlib.Path, run_name: str) -> pathlib.Path:
    plot_dir = cdir.expanduser().resolve() / "figs" / "icetune" / run_name
    campaign = campaign_parts(history_path)
    if campaign is not None:
        plot_dir /= pathlib.Path("campaigns") / campaign[1]
    return plot_dir


# Draw one chronological trial and running minimum cost figure
def plot_cost_evolution(
    *,
    history: dict,
    cost_key: str,
    points: list[CostPoint],
    output_path: pathlib.Path,
    run_name: str,
    log_y: bool = False,
) -> pathlib.Path:
    import matplotlib.pyplot as plt

    from core.plot.style import add_icetune_branding, resolve_plot_brand

    ordered, elapsed, costs, average, sigma, _ = cost_evolution(points)
    trajectories = start_trajectories(ordered, history)
    scale, unit = _time_axis(elapsed)
    x_values = [value / scale for value in elapsed]
    parameter_space = history.get("parameter_space")
    parameter_names = parameter_space if isinstance(parameter_space, dict) else ()
    plot_brand = resolve_plot_brand(
        saved_brand=history.get("plot_brand"), simdriver=history.get("simdriver"), param_names=parameter_names
    )

    ensure_dir(output_path.parent)
    with plt.rc_context({"axes.labelsize": 10, "font.size": 9, "xtick.labelsize": 8, "ytick.labelsize": 8}):
        fig, ax = plt.subplots(figsize=(7.2, 4.4))
        if log_y:
            ax.set_yscale("log")
        log_valid = not log_y or all(cost > 0.0 for cost in costs)
        if ordered and log_valid:
            x_limit = max(1.0, max(x_values) * 1.03)
            regime_end = random_regime_end(history, ordered)
            if regime_end is not None:
                shade_end = x_limit if math.isinf(regime_end) else min(x_limit, regime_end / scale)
                if shade_end > 0.0:
                    ax.axvspan(0.0, shade_end, color="#eeeeee", alpha=0.85, linewidth=0.0, zorder=0)
            if trajectories:
                _draw_trajectories(ax, trajectories, x_values, costs, legend=_start_count(history, ordered) <= 10)
                _start_colorbar(ax, history, ordered)
            else:
                ax.scatter(x_values, costs, s=24, color="#777777", edgecolor="white", linewidth=0.5, zorder=3)
            initial_index = _initial_index(ordered)
            if initial_index is not None:
                initial_xy = (x_values[initial_index], costs[initial_index])
                ax.scatter(
                    [initial_xy[0]],
                    [initial_xy[1]],
                    marker="o",
                    s=50,
                    facecolor="#228833",
                    edgecolor="#228833",
                    linewidth=1.2,
                    zorder=5,
                    clip_on=False,
                )
            if not trajectories:
                lower = [mean - spread for mean, spread in zip(average, sigma, strict=True)]
                upper = [mean + spread for mean, spread in zip(average, sigma, strict=True)]
                ax.fill_between(x_values, lower, upper, color="#4477AA", alpha=0.18, linewidth=0.0, zorder=1)
                ax.plot(x_values, average, color="#4477AA", linewidth=1.7, zorder=2)
            if trajectories:
                ax.legend(frameon=False, fontsize=7, ncol=3)
            best_index = min(range(len(costs)), key=costs.__getitem__)
            ax.scatter(
                [x_values[best_index]],
                [costs[best_index]],
                marker="o",
                s=150,
                facecolor="none",
                edgecolor="#CC3311",
                linewidth=1.6,
                zorder=6,
                clip_on=False,
            )
            ax.autoscale(enable=True, axis="y", tight=True)
            if log_y:
                lower, upper = ax.get_ylim()
                ax.set_ylim(bottom=lower * (lower / upper)**plt.rcParams["axes.ymargin"])
            elif all(cost > 0.0 for cost in costs):
                ax.set_ylim(bottom=0.0)
            ax.set_xlim(0.0, x_limit)
            failed = history.get("failed")
            failed_count = len(failed) if isinstance(failed, list) else 0
            target = history.get("num_trials")
            target_count = int(target) if type(target) is int and target > 0 else None
            count = f"successful {len(costs)}"
            if target_count is not None:
                count += f"/{target_count}"
            if failed_count:
                count += f", failed {failed_count}"
            ax.set_title(f"best cost {costs[best_index]:.4g}, {count} ({ordered[-1].completed_at})")
        elif ordered:
            ax.text(
                0.5,
                0.5,
                "Log scale requires positive costs",
                color="#666666",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_xlim(0.0, 1.0)
            ax.set_title("log cost evolution unavailable")
        else:
            ax.text(0.5, 0.5, "No completed trials", color="#666666", ha="center", va="center", transform=ax.transAxes)
            ax.set_xlim(0.0, 1.0)
            ax.set_title("best cost n/a, trials 0 (n/a)")

        ax.set_xlabel(f"Time [{unit}]")
        ax.set_ylabel(f"Cost ({cost_key})")
        ax.grid(axis="y", color="#dddddd", linewidth=0.7)
        if not log_y:
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4), useMathText=True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.subplots_adjust(left=0.13, right=0.97, bottom=0.15, top=0.87)
        add_icetune_branding(fig, brand=plot_brand, x=0.02, y=0.98)
        fig.text(0.97, 0.98, run_name, ha="right", va="top", fontsize=9, color="#555555")
        temporary = output_path.with_name(f".{output_path.stem}.tmp.{uuid.uuid4().hex}{output_path.suffix}")
        try:
            fig.savefig(temporary, dpi=220, bbox_inches="tight", facecolor="white")
            temporary.replace(output_path)
        finally:
            plt.close(fig)
            temporary.unlink(missing_ok=True)
    return output_path


# Draw one descending cost curve with the minimum at the final trial index
def plot_sorted_costs(
    *,
    history: dict,
    cost_key: str,
    points: list[CostPoint],
    output_path: pathlib.Path,
    run_name: str,
    log_y: bool = False,
) -> pathlib.Path:
    import matplotlib.pyplot as plt

    from core.plot.style import add_icetune_branding, resolve_plot_brand

    ordered = sorted_costs(points)
    trajectories = start_trajectories(ordered, history)
    costs = [point.cost for point in ordered]
    x_values = list(range(1, len(ordered) + 1))
    parameter_space = history.get("parameter_space")
    parameter_names = parameter_space if isinstance(parameter_space, dict) else ()
    plot_brand = resolve_plot_brand(
        saved_brand=history.get("plot_brand"), simdriver=history.get("simdriver"), param_names=parameter_names
    )

    ensure_dir(output_path.parent)
    with plt.rc_context({"axes.labelsize": 10, "font.size": 9, "xtick.labelsize": 8, "ytick.labelsize": 8}):
        fig, ax = plt.subplots(figsize=(7.2, 4.4))
        if log_y:
            ax.set_yscale("log")
        log_valid = not log_y or all(cost > 0.0 for cost in costs)
        if ordered and log_valid:
            ax.plot(x_values, costs, color="#4477AA", linewidth=1.7, zorder=2)
            if trajectories:
                _draw_trajectories(ax, trajectories, x_values, costs, connect=False, legend=_start_count(history, ordered) <= 10)
                _start_colorbar(ax, history, ordered)
                if _start_count(history, ordered) <= 10:
                    ax.legend(frameon=False, fontsize=7, ncol=3)
            ax.autoscale(enable=True, axis="y", tight=True)
            if not log_y and all(cost > 0.0 for cost in costs):
                ax.set_ylim(bottom=0.0)
            if len(ordered) == 1:
                ax.set_xlim(0.5, 1.5)
            else:
                ax.set_xlim(1.0, len(ordered))
            ax.set_title(f"minimum cost {costs[-1]:.4g}, successful {len(costs)}")
        elif ordered:
            ax.text(
                0.5,
                0.5,
                "Log scale requires positive costs",
                color="#666666",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_xlim(0.0, 1.0)
            ax.set_title("log sorted costs unavailable")
        else:
            ax.text(0.5, 0.5, "No completed trials", color="#666666", ha="center", va="center", transform=ax.transAxes)
            ax.set_xlim(0.0, 1.0)
            ax.set_title("minimum cost n/a, trials 0")

        ax.set_xlabel("Sorted trials")
        ax.set_ylabel(f"Cost ({cost_key})")
        ax.grid(axis="y", color="#dddddd", linewidth=0.7)
        if not log_y:
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4), useMathText=True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.subplots_adjust(left=0.13, right=0.97, bottom=0.15, top=0.87)
        add_icetune_branding(fig, brand=plot_brand, x=0.02, y=0.98)
        fig.text(0.97, 0.98, run_name, ha="right", va="top", fontsize=9, color="#555555")
        temporary = output_path.with_name(f".{output_path.stem}.tmp.{uuid.uuid4().hex}{output_path.suffix}")
        try:
            fig.savefig(temporary, dpi=220, bbox_inches="tight", facecolor="white")
            temporary.replace(output_path)
        finally:
            plt.close(fig)
            temporary.unlink(missing_ok=True)
    return output_path


# Draw paginated parameter values against elapsed runtime into one PDF
def plot_parameter_evolution(
    *, history: dict, points: list[CostPoint], output_path: pathlib.Path, run_name: str, physical: bool
) -> pathlib.Path:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    from core.plot.style import add_icetune_branding, resolve_plot_brand

    ordered, names, series = parameter_evolution(history, points, physical=physical)
    trajectories = start_trajectories(ordered, history)
    initial_index = _initial_index(ordered)
    best_index = min(range(len(ordered)), key=lambda index: ordered[index].cost) if ordered else None
    origin = ordered[0].completed_unix if ordered else 0.0
    elapsed = [point.completed_unix - origin for point in ordered]
    scale, unit = _time_axis(elapsed)
    x_values = [value / scale for value in elapsed]
    x_limit = max(1.0, max(x_values, default=0.0) * 1.03)
    regime_end = random_regime_end(history, ordered)
    bounds = {} if physical else _parameter_bounds(history)
    parameter_space = history.get("parameter_space")
    parameter_names = parameter_space if isinstance(parameter_space, dict) else ()
    plot_brand = resolve_plot_brand(
        saved_brand=history.get("plot_brand"), simdriver=history.get("simdriver"), param_names=parameter_names
    )
    basis = "Physical" if physical else "Raw optimizer"
    optimizer_stage = trial_stages(history, ordered)
    rows, columns = 4, 2
    per_page = rows
    pages = max(1, math.ceil(len(names) / per_page))

    ensure_dir(output_path.parent)
    temporary = output_path.with_name(f".{output_path.stem}.tmp.{uuid.uuid4().hex}{output_path.suffix}")
    try:
        with PdfPages(temporary) as pdf:
            for page in range(pages):
                with plt.rc_context({"axes.labelsize": 9, "font.size": 8, "xtick.labelsize": 7, "ytick.labelsize": 7}):
                    fig, axes = plt.subplots(rows, columns, figsize=(11.7, 8.3), squeeze=False)
                    colorbar_legends = []
                    selected = names[page * per_page : (page + 1) * per_page]
                    if not selected:
                        axes[0][0].text(
                            0.5,
                            0.5,
                            f"No {basis.lower()} parameter values available",
                            color="#666666",
                            ha="center",
                            va="center",
                            transform=axes[0][0].transAxes,
                        )
                        axes[0][0].set_xlim(0.0, 1.0)
                    for row, name in enumerate(selected):
                        trace_axis = axes[row][0]
                        density_axis = axes[row][1]
                        if regime_end is not None:
                            shade_end = x_limit if math.isinf(regime_end) else min(x_limit, regime_end / scale)
                            if shade_end > 0.0:
                                trace_axis.axvspan(0.0, shade_end, color="#eeeeee", alpha=0.85, linewidth=0.0, zorder=0)
                        values = series[name]
                        for index, title, color, size, facecolor, linewidth in (
                            (initial_index, "Initial", "#228833", 50, "#228833", 1.2),
                            (best_index, "Best", "#CC3311", 150, "none", 1.6),
                        ):
                            if index is None or not math.isfinite(values[index]):
                                continue
                            value = values[index]
                            label = f"{title}: {value:.4g}"
                            for axis, x, y in ((trace_axis, x_values[index], value), (density_axis, value, 0.0)):
                                axis.scatter(
                                    [x],
                                    [y],
                                    marker="o",
                                    s=size,
                                    facecolor=facecolor,
                                    edgecolor=color,
                                    linewidth=linewidth,
                                    zorder=4,
                                    clip_on=False,
                                )
                            density_axis.axvline(
                                value, color=color, linewidth=2.5, linestyle="-", alpha=0.5, label=label, zorder=4
                            )
                        if trajectories:
                            _draw_trajectories(trace_axis, trajectories, x_values, values)
                        else:
                            trace_axis.scatter(
                                x_values, values, s=15, color="#777777", edgecolor="white", linewidth=0.4, zorder=2
                            )
                        if name in bounds:
                            lower, upper = bounds[name]
                            trace_axis.axhline(lower, color="#999999", linewidth=0.7, zorder=1)
                            trace_axis.axhline(upper, color="#999999", linewidth=0.7, zorder=1)
                            density_axis.axvline(lower, color="#999999", linewidth=0.7, zorder=1)
                            density_axis.axvline(upper, color="#999999", linewidth=0.7, zorder=1)
                        random_values = [
                            value
                            for value, is_optimizer in zip(values, optimizer_stage, strict=True)
                            if not is_optimizer and math.isfinite(value)
                        ]
                        optimized_values = [
                            value
                            for value, is_optimizer in zip(values, optimizer_stage, strict=True)
                            if is_optimizer and math.isfinite(value)
                        ]
                        combined = random_values + optimized_values
                        if combined:
                            lower = min(combined)
                            upper = max(combined)
                            if math.isclose(lower, upper, rel_tol=1.0e-12, abs_tol=1.0e-12):
                                padding = max(0.5, abs(lower) * 0.05)
                                lower -= padding
                                upper += padding
                            bins = min(50, max(5, math.ceil(math.sqrt(len(combined)))))
                            distributions = (
                                [(label, color, [values[index] for index in indices if math.isfinite(values[index])])
                                 for label, color, indices in trajectories]
                                if trajectories else
                                [("Exploration", "#222222", random_values), ("Optimizer", "#CC3311", optimized_values)]
                            )
                            for label, color, samples in distributions:
                                if not samples:
                                    continue
                                density_axis.hist(
                                    samples,
                                    bins=bins,
                                    range=(lower, upper),
                                    density=True,
                                    histtype="step",
                                    color=color,
                                    linewidth=1.4,
                                    label=label if not trajectories or _start_count(history, ordered) <= 10 else "_nolegend_",
                                    zorder=2,
                                )
                            bar = _start_colorbar(density_axis, history, ordered) if trajectories else None
                            legend = density_axis.legend(
                                frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0
                            )
                            if bar is not None:
                                colorbar_legends.append((bar, legend))
                            if best_index is not None:
                                trial = ordered[best_index].trial_id.removeprefix("trial-")
                                density_axis.annotate(
                                    f"Best trial: {int(trial) if trial.isascii() and trial.isdigit() else trial}",
                                    xy=(0, 0), xycoords=legend, xytext=(0, -5), textcoords="offset points",
                                    ha="left", va="top", fontsize=7, annotation_clip=False,
                                )
                        else:
                            density_axis.text(
                                0.5,
                                0.5,
                                "No parameter values",
                                color="#666666",
                                ha="center",
                                va="center",
                                transform=density_axis.transAxes,
                            )
                        trace_axis.set_title(name, loc="left", fontsize=7.5, pad=3.0)
                        trace_axis.set_xlim(0.0, x_limit)
                        trace_axis.set_xlabel(f"Time [{unit}]")
                        density_axis.set_xlabel("Parameter value")
                        density_axis.set_ylabel("Density", fontsize=8)
                        density_axis.set_ylim(bottom=0.0)
                        trace_axis.grid(axis="y", color="#dddddd", linewidth=0.6)
                        density_axis.grid(axis="y", color="#dddddd", linewidth=0.6)
                        trace_axis.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4), useMathText=True)
                        density_axis.ticklabel_format(axis="both", style="sci", scilimits=(-3, 4), useMathText=True)
                        for axis in (trace_axis, density_axis):
                            axis.spines["top"].set_visible(False)
                            axis.spines["right"].set_visible(False)
                    first_unused_row = len(selected) if selected else 1
                    for row in range(first_unused_row, rows):
                        axes[row][0].set_visible(False)
                        axes[row][1].set_visible(False)
                    fig.subplots_adjust(
                        left=0.05, right=0.82 if colorbar_legends else 0.90,
                        bottom=0.07, top=0.92, hspace=0.62, wspace=0.12,
                    )
                    if colorbar_legends:
                        fig.canvas.draw()
                        renderer = fig.canvas.get_renderer()
                        for bar, legend in colorbar_legends:
                            # Place the summary beyond the colorbar ticks and label
                            box = bar.ax.get_tightbbox(renderer).transformed(fig.transFigure.inverted())
                            legend.set_bbox_to_anchor(
                                (box.x1 + 0.01, legend.axes.get_position().y1), transform=fig.transFigure
                            )
                    add_icetune_branding(fig, brand=plot_brand, x=0.02, y=0.985)
                    fig.text(0.98, 0.985, run_name, ha="right", va="top", fontsize=9, color="#555555")
                    fig.text(0.98, 0.02, f"{page + 1}/{pages}", ha="right", va="bottom", fontsize=8, color="#777777")
                    pdf.savefig(fig, facecolor="white")
                    plt.close(fig)
        temporary.replace(output_path)
    finally:
        temporary.unlink(missing_ok=True)
    return output_path


# Resolve, load, and draw one history source without command line output
def create_history_plot(
    *, source: pathlib.Path, cdir: pathlib.Path, cost: str | None = None, run_name: str | None = None
) -> pathlib.Path:
    history_path = resolve_history_path(source, cdir)
    history, cost_key, points = load_history(history_path, cost)
    resolved_name = resolve_run_name(history_path, run_name)
    plot_dir = resolve_plot_dir(history_path=history_path, cdir=cdir, run_name=resolved_name)
    evolution_path = plot_cost_evolution(
        history=history,
        cost_key=cost_key,
        points=points,
        output_path=plot_dir / "cost_evolution.png",
        run_name=resolved_name,
    )
    plot_cost_evolution(
        history=history,
        cost_key=cost_key,
        points=points,
        output_path=plot_dir / "cost_evolution_logy.png",
        run_name=resolved_name,
        log_y=True,
    )
    plot_sorted_costs(
        history=history,
        cost_key=cost_key,
        points=points,
        output_path=plot_dir / "cost_sorted.png",
        run_name=resolved_name,
    )
    plot_sorted_costs(
        history=history,
        cost_key=cost_key,
        points=points,
        output_path=plot_dir / "cost_sorted_logy.png",
        run_name=resolved_name,
        log_y=True,
    )
    plot_parameter_evolution(
        history=history,
        points=points,
        output_path=plot_dir / "parameter_evolution_physical.pdf",
        run_name=resolved_name,
        physical=True,
    )
    plot_parameter_evolution(
        history=history,
        points=points,
        output_path=plot_dir / "parameter_evolution_optimizer.pdf",
        run_name=resolved_name,
        physical=False,
    )
    return evolution_path


# Render a completed history snapshot before briefly locking the published figures
def publish_history_plot(
    *, history: dict, source: pathlib.Path, cdir: pathlib.Path, cost: str, run_name: str
) -> pathlib.Path:
    from core.tune import core as icetune_main
    from core.tune.io import atomic_write_json

    target = resolve_plot_dir(history_path=source, cdir=cdir, run_name=run_name)
    scratch = cdir / "tmp"
    ensure_dir(scratch)
    with tempfile.TemporaryDirectory(prefix="icetune-history-", dir=scratch) as temporary:
        stage = pathlib.Path(temporary)
        snapshot = stage / "history.json"
        atomic_write_json(str(snapshot), history)
        path = create_history_plot(source=snapshot, cdir=stage, cost=cost, run_name=run_name)
        token = icetune_main.acquire_publish_lock(cdir=str(cdir), run_name=run_name)
        try:
            ensure_dir(target)
            for name in icetune_main.ICETUNE_HISTORY_FIGURES:
                (path.parent / name).replace(target / name)
        finally:
            icetune_main.release_publish_lock(cdir=str(cdir), run_name=run_name, token=token)
    return target / path.name


# Run the history reader and convert malformed input into a command line error
def main(argv: list[str] | None = None) -> None:
    parser = create_parser()
    args = parser.parse_args(argv)
    try:
        history_path = resolve_history_path(args.source, args.cdir)
        history, cost_key, points = load_history(history_path, args.cost)
        plot_path = create_history_plot(source=history_path, cdir=args.cdir, cost=args.cost, run_name=args.run_name)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    if args.sort:
        points.sort(key=lambda point: point.cost, reverse=True)
    print(format_cost_values(cost_key, points))
    print(f"\nCost evolution: {plot_path}")
    print(f"Log cost evolution: {plot_path.with_name('cost_evolution_logy.png')}")
    print(f"Sorted costs: {plot_path.with_name('cost_sorted.png')}")
    print(f"Log sorted costs: {plot_path.with_name('cost_sorted_logy.png')}")
    print(f"Physical parameters: {plot_path.with_name('parameter_evolution_physical.pdf')}")
    print(f"Optimizer parameters: {plot_path.with_name('parameter_evolution_optimizer.pdf')}")


if __name__ == "__main__":
    main()

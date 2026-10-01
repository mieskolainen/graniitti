# Shared tuning I/O and proposal primitives
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import contextlib
import copy
import json
import math
import os
import time
import uuid
from datetime import UTC, datetime

from core.io.files import IO_RETRIES, ensure_dir, is_retryable_io_error, retry_delay
from core.io.files import retry_io as _retry
from core.io.serialize import load_json_file

AUTO_BATCH_MAX = 16


# Compute the bounded joint acquisition block within the concurrency limit
def proposal_batch_size(args) -> int:
    target = max(1, int(getattr(args, "max_concurrent_trials", 1)))
    configured = max(0, int(getattr(args, "proposal_batch_size", 0)))
    fraction = float(getattr(args, "proposal_batch_fraction", 0.1))
    automatic = min(AUTO_BATCH_MAX, max(2, math.ceil(fraction * target)))
    return min(target, configured or automatic)


# Compute whether proposals at one trial index need no pending result
def proposals_independent(args, index: int) -> bool:
    if str(getattr(args, "algorithm", "basic")) == "basic":
        return True
    warmup = max(0, int(getattr(args, "rand_trials", 0)))
    return int(index) < warmup


# Compute free proposal slots using the shared asynchronous phase rules
def proposal_slots(args, *, completed: int, pending: int, next_number: int) -> int:
    target = max(1, int(getattr(args, "max_concurrent_trials", 1)))
    total = int(getattr(args, "num_trials", int(next_number) + target))
    available = max(0, min(total - int(completed) - int(pending), target - int(pending)))
    if available == 0:
        return 0
    algorithm = str(getattr(args, "algorithm", "basic"))
    if algorithm == "basic":
        return available
    warmup = max(0, int(getattr(args, "rand_trials", 0)))
    if proposals_independent(args, next_number):
        return min(available, warmup - int(next_number))
    return available


# Build paired unix and display timestamp fields
def time_fields(name: str, timestamp: float | int | None) -> dict:
    value = None if timestamp is None else float(timestamp)
    return {
        f"{name}_unix": value,
        f"{name}_datetime": None
        if value is None
        else datetime.fromtimestamp(value, tz=UTC).astimezone().strftime("%Y-%m-%d %H:%M:%S %Z"),
    }


# Record worker wall time with consistent history timestamps
def trial_times(started: float | None, completed: float) -> dict:
    return {
        **time_fields("started_at", started),
        **time_fields("completed_at", completed),
        "elapsed_seconds": None if started is None else max(0.0, completed - started),
    }


# Compute whether a path exists
def path_exists(path: str, *, retries: int = IO_RETRIES) -> bool:
    return _retry(lambda: os.stat(path), retries=retries, missing=False) is not False


# Unlink one private temporary path
def unlink_path(path: str, *, retries: int = IO_RETRIES) -> bool:
    return _retry(lambda: os.unlink(path) or True, retries=retries, missing=False)


# Compute a collision-resistant temporary sibling path
def tmp_path(path: str) -> str:
    return f"{path}.tmp.{os.getpid()}.{time.time_ns()}.{uuid.uuid4().hex}"


# Serialize one JSON payload and prepare its parent
def _json_data(path: str, payload: dict, ensure_parent: bool) -> str:
    if ensure_parent:
        ensure_dir(os.path.dirname(path))
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


# Write one complete JSON payload through a private temporary file
def atomic_write_json(path: str, payload: dict, *, ensure_parent: bool = True, retries: int = IO_RETRIES) -> None:
    data = _json_data(path, payload, ensure_parent)
    attempts = max(1, int(retries))
    for attempt in range(1, attempts + 1):
        temporary = tmp_path(path)
        try:
            with open(temporary, "x", encoding="utf-8") as handle:
                handle.write(data)
            os.replace(temporary, path)
            return
        except OSError as exc:
            unlink_path(temporary)
            if attempt >= attempts or not is_retryable_io_error(exc):
                raise
            time.sleep(retry_delay(attempt, base=0.05, cap=2.0))


# Write one immutable JSON payload with a hard-link commit
def write_once_json(path: str, payload: dict, *, ensure_parent: bool = True) -> bool:
    data = _json_data(path, payload, ensure_parent)
    for attempt in range(1, IO_RETRIES + 1):
        temporary = tmp_path(path)
        try:
            with open(temporary, "x", encoding="utf-8") as handle:
                handle.write(data)
            os.link(temporary, path)
            return True
        except FileExistsError:
            return False
        except OSError as exc:
            with contextlib.suppress(OSError):
                if os.path.samefile(temporary, path):
                    return True
            if path_exists(path) is True:
                return False
            if attempt >= IO_RETRIES or not is_retryable_io_error(exc):
                raise
            time.sleep(retry_delay(attempt, base=0.05, cap=2.0))
        finally:
            with contextlib.suppress(OSError):
                os.unlink(temporary)
    return False


# Read JSON with bounded retries for transient and partial visibility
def read_json(path: str, default=None, *, retries: int = IO_RETRIES, retry_missing: bool = False):
    for attempt in range(1, max(1, int(retries)) + 1):
        try:
            with open(path, encoding="utf-8") as handle:
                return load_json_file(handle.name)
        except FileNotFoundError:
            if not retry_missing or attempt >= retries:
                return copy.deepcopy(default)
        except (ValueError, OSError) as exc:
            # Inspect reference-reader causes without retrying invalid card structure
            cause = exc
            while cause.__cause__ is not None:
                cause = cause.__cause__
            retryable = isinstance(cause, json.JSONDecodeError) or (
                isinstance(cause, OSError) and is_retryable_io_error(cause)
            )
            if attempt >= retries or not retryable:
                raise
        time.sleep(retry_delay(attempt, base=0.05, cap=2.0))
    return copy.deepcopy(default)

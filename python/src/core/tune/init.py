# Ray initialization state and bootstrap execution
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import os
import time
import uuid
from collections.abc import Callable

from core.io.files import ensure_dir
from core.io.serialize import json_safe
from core.tune.cache import PermanentConfigurationError
from core.tune.io import atomic_write_json, read_json, time_fields

PROTOCOL_NAME = "core.tune.ray.bootstrap"
PROTOCOL_VERSION = 1


# Run one serialized Ray bootstrap initialization job
def run_initialization(
    *,
    args,
    bootstrapper: Callable[[], dict],
    fingerprint_getter: Callable[[], str] | None = None,
    initial_point_getter: Callable[[], dict] | None = None,
) -> dict:
    if str(getattr(args, "backend", "")) != "ray":
        raise ValueError("Ray initialization requires --backend ray")
    if getattr(args, "phase", None) != "init":
        raise ValueError("Ray initialization requires --phase init")
    if not getattr(args, "coord_dir", None):
        raise ValueError("--coord_dir is required with --phase init")

    coord_dir = os.path.abspath(args.coord_dir)
    ensure_dir(coord_dir)
    init_path = os.path.join(coord_dir, "init.json")
    expected = fingerprint_getter() if fingerprint_getter is not None else None
    previous = read_json(init_path, default={}) or {}
    if previous.get("status") == "completed":
        fingerprint = (previous.get("bootstrap") or {}).get("fingerprint")
        if (
            previous.get("protocol") != PROTOCOL_NAME
            or previous.get("protocol_version") != PROTOCOL_VERSION
            or (expected is not None and fingerprint != expected)
        ):
            raise PermanentConfigurationError("Existing initialization directory has incompatible state")
        return previous

    attempt_id = uuid.uuid4().hex
    started = time.time()
    running = {
        "attempt_id": attempt_id,
        "backend": "ray",
        "protocol": PROTOCOL_NAME,
        "protocol_version": PROTOCOL_VERSION,
        "status": "running",
        **time_fields("started_at", started),
        **time_fields("updated_at", started),
    }
    atomic_write_json(init_path, running, ensure_parent=False)
    try:
        bootstrap = bootstrapper()
        if not isinstance(bootstrap, dict) or not bootstrap.get("fingerprint"):
            raise PermanentConfigurationError("Ray init produced no bootstrap fingerprint")
        if expected is not None and bootstrap["fingerprint"] != expected:
            raise PermanentConfigurationError("Ray init bootstrap fingerprint differs from its expected identity")
        points = None
        if not getattr(args, "no_initial_point", False) and initial_point_getter is not None:
            points = initial_point_getter()
        completed = time.time()
        state = {
            "attempt_id": attempt_id,
            "backend": "ray",
            "bootstrap": json_safe(bootstrap),
            "initial_points": json_safe(points),
            "protocol": PROTOCOL_NAME,
            "protocol_version": PROTOCOL_VERSION,
            "status": "completed",
            **time_fields("completed_at", completed),
            **time_fields("started_at", started),
            **time_fields("updated_at", completed),
        }
        atomic_write_json(init_path, state, ensure_parent=False)
        return state
    except Exception as exc:
        failed = {
            **running,
            "error": str(exc),
            "error_type": exc.__class__.__name__,
            "status": "failed",
            **time_fields("updated_at", time.time()),
        }
        try:
            atomic_write_json(init_path, failed, ensure_parent=False)
        except OSError as error:
            exc.add_note(f"Could not write initialization failure to {init_path}: {error}")
        raise

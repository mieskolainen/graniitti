# Shared timestamped output directories and reproducibility manifests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import pathlib
import re
from datetime import UTC, datetime

from core.io.files import ensure_dir
from core.io.serialize import json_safe, load_json_file
from core.tune.io import atomic_write_json

RUN_MANIFEST_SCHEMA_VERSION = 1


# Convert run arguments into lossless JSON-compatible values
def _manifest_value(value):
    return json_safe(value, nonfinite="string")


# Compute one filesystem-safe run label
def _run_label(value: str) -> str:
    candidate = pathlib.Path(str(value)).name
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", candidate).strip("._")
    return label or "run"


# Write one JSON object atomically inside its output directory
def _atomic_write_json(path: pathlib.Path, payload: dict) -> None:
    atomic_write_json(str(path), _manifest_value(payload))


# Allocate one collision-safe timestamped run output and initial manifest
def create_run_output(
    *,
    cdir: str,
    tool: str,
    tool_version,
    input_run_name: str,
    arguments: dict,
    command: list[str],
    output_name: str | None = None,
    timestamp: datetime | None = None,
) -> pathlib.Path:
    created = timestamp or datetime.now(UTC)
    if created.tzinfo is None:
        created = created.replace(tzinfo=UTC)
    created = created.astimezone(UTC)
    timestamp_label = created.strftime("%Y%m%dT%H%M%S.%fZ")
    base_label = _run_label(output_name or input_run_name)
    parent = pathlib.Path(cdir) / "figs" / str(tool)
    ensure_dir(parent)

    output_root = parent / f"{base_label}__{timestamp_label}"
    collision = 0
    while True:
        candidate = output_root if collision == 0 else output_root.with_name(f"{output_root.name}__{collision:02d}")
        try:
            ensure_dir(candidate, exist_ok=False)
            output_root = candidate
            break
        except FileExistsError:
            collision += 1

    manifest = {
        "run_manifest_schema_version": RUN_MANIFEST_SCHEMA_VERSION,
        "tool": str(tool),
        "tool_version": list(tool_version),
        "run_id": output_root.name,
        "input_run_name": str(input_run_name),
        "created_at_utc": created.isoformat(),
        "status": "running",
        "command": [str(item) for item in command],
        "arguments": _manifest_value(arguments),
        "output_root": str(output_root.resolve()),
    }
    _atomic_write_json(output_root / "manifest.json", manifest)
    return output_root


# Merge reproducibility or completion fields into one run manifest
def update_run_manifest(output_root: pathlib.Path, updates: dict) -> dict:
    manifest_path = pathlib.Path(output_root) / "manifest.json"
    with manifest_path.open("r", encoding="utf-8") as handle:
        manifest = load_json_file(handle.name)
    manifest.update(_manifest_value(updates))
    _atomic_write_json(manifest_path, manifest)
    return manifest

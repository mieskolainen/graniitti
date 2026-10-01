# Tests for timestamped tuning run manifests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from datetime import UTC, datetime

from core.tune import manifest as run_manifest


# Check timestamped run outputs coexist and preserve exact input arguments
def test_timestamped_run_manifests_collision_safe(tmp_path):
    timestamp = datetime(2026, 7, 28, 12, 34, 56, tzinfo=UTC)
    arguments = {
        "run_name": "GP/input",
        "max_Z": float("inf"),
        "models": ["gp", "neural"],
    }
    first = run_manifest.create_run_output(
        cdir=str(tmp_path),
        tool="icescape",
        tool_version=(0, 0, 7),
        input_run_name="GP/input",
        arguments=arguments,
        command=["python", "python/src/core/icescape.py"],
        timestamp=timestamp,
    )
    second = run_manifest.create_run_output(
        cdir=str(tmp_path),
        tool="icescape",
        tool_version=(0, 0, 7),
        input_run_name="GP/input",
        arguments=arguments,
        command=["python", "python/src/core/icescape.py"],
        timestamp=timestamp,
    )

    assert first.name == "input__20260728T123456.000000Z"
    assert second.name == "input__20260728T123456.000000Z__01"
    first_payload = json.loads((first / "manifest.json").read_text(encoding="utf-8"))
    assert first_payload["arguments"]["max_Z"] == "inf"
    assert first_payload["status"] == "running"

    completed = run_manifest.update_run_manifest(
        first,
        {"status": "completed", "elapsed_seconds": 12.5},
    )

    assert completed["status"] == "completed"
    assert completed["elapsed_seconds"] == 12.5

# HTCondor scratch paths, process identity and submission logging
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import hashlib
import os
import pathlib
import re
from datetime import datetime

from core.io.files import ensure_dir


# Print one timestamped steering message immediately
def write_log(message: str, *, file=None) -> None:
    timestamp = datetime.now().astimezone().isoformat(timespec="seconds")
    print(f"[{timestamp}] {message}", file=file, flush=True)


# Compute the DAGMan exit code for one steering failure
def steer_exit_code(exc: BaseException) -> int:
    return 64 if isinstance(exc, (KeyError, ValueError)) else 1


# Compute the stable campaign scratch root shared by head worker tasks
def runtime_work_root(environment: dict[str, str]) -> pathlib.Path:
    run_name = str(environment["RUN_NAME"])
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", run_name) is None:
        raise ValueError(f"Invalid Ray lxplus run name: {run_name!r}")
    return condor_scratch_dir() / f"icetune-ray-{run_name}"


# Compute fast local scratch for the current HTCondor process
def condor_scratch_dir() -> pathlib.Path:
    scratch = pathlib.Path(os.environ.get("_CONDOR_SCRATCH_DIR", "/tmp"))
    if not scratch.is_absolute():
        raise ValueError(f"HTCondor scratch path must be absolute: {scratch}")
    ensure_dir(scratch)
    return scratch


# Compute cache and temporary environment paths below local Condor scratch
def condor_scratch_environment() -> dict[str, str]:
    scratch = condor_scratch_dir()
    cache = scratch / "cache"
    paths = {"MPLCONFIGDIR": cache / "matplotlib", "NUMBA_CACHE_DIR": cache / "numba",
        "PYTHONPYCACHEPREFIX": cache / "pycache", }
    for path in paths.values():
        ensure_dir(path)
    return {"TMPDIR": str(scratch), "TMP": str(scratch), "TEMP": str(scratch), "XDG_CACHE_HOME": str(cache),
        **{key: str(path) for key, path in paths.items()}, }


# Compute the head local process record used for graceful driver shutdown
def runtime_pid_path(environment: dict[str, str]) -> pathlib.Path:
    return runtime_work_root(environment) / "icetune-driver.json"


# Compute one Linux process birth marker so reused PIDs cannot inherit leases
def process_birth(pid: int) -> str | None:
    try:
        payload = pathlib.Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
    except OSError:
        return None
    fields = payload.rsplit(")", 1)
    if len(fields) != 2:
        return None
    values = fields[1].split()
    return values[19] if len(values) > 19 else None


# Compute the SHA256 identity of one full trial pickle
def file_sha256(path: pathlib.Path) -> str:
    with open(path, "rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()

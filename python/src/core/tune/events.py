# Publish complete MC event samples from trial scratch to shared storage
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import json
import os
import pathlib
import re
import shutil
import tempfile

from core.io.files import ensure_dir
from core.io.serialize import json_safe
from core.tune import short_id
from core.tune.cache import json_fingerprint, sha256_file


# Signal a transfer failure which Ray must retry instead of accepting a trial penalty
class EventTransferError(RuntimeError):
    pass


# Verify every file named by a completed event manifest
def verify_samples(directory: pathlib.Path, manifest: dict) -> None:
    for name, info in manifest["files"].items():
        relative = pathlib.PurePosixPath(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe event sample path: {name}")
        path = directory / relative
        if path.is_symlink() or not path.is_file() or path.stat().st_size != info["size"] or sha256_file(path) != info["sha256"]:
            raise EventTransferError(f"Incomplete event sample: {path}")


# Stream and verify a trial sample tree before atomically publishing its manifest
def publish_samples(*, source: pathlib.Path, destination: pathlib.Path, trial_id: str, metadata: dict) -> dict:
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", trial_id) is None:
        raise ValueError("Event output trial_id must be one safe path component")
    files = {}
    for path in sorted(source.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Event sample contains a symbolic link: {path}")
        if path.is_file():
            files[path.relative_to(source).as_posix()] = {"size": path.stat().st_size, "sha256": sha256_file(path)}
    if not any(name.endswith(".hepmc3") for name in files):
        raise EventTransferError("Trial produced no event samples to save")
    manifest = json_safe({"schema_version": 1, "status": "complete", "trial_id": trial_id,
                          "metadata": metadata, "files": files})
    attempt = json_fingerprint(manifest)
    target = destination / trial_id / short_id(attempt)
    try:
        if target.exists():
            if json.loads((target / "manifest.json").read_text()) != manifest:
                raise EventTransferError(f"Conflicting event manifest: {target}")
            verify_samples(target, manifest)
        else:
            ensure_dir(target.parent)
            with tempfile.TemporaryDirectory(prefix=f".{short_id(attempt)}.", dir=target.parent) as temporary:
                stage = pathlib.Path(temporary)
                for name in files:
                    output = stage / name
                    ensure_dir(output.parent)
                    with (source / name).open("rb") as incoming, output.open("xb") as outgoing:
                        shutil.copyfileobj(incoming, outgoing)
                        outgoing.flush()
                        os.fsync(outgoing.fileno())
                verify_samples(stage, manifest)
                with (stage / "manifest.json").open("x") as handle:
                    json.dump(manifest, handle, sort_keys=True, indent=2)
                    handle.flush()
                    os.fsync(handle.fileno())
                try:
                    stage.rename(target)
                except OSError:
                    if not target.is_dir() or json.loads((target / "manifest.json").read_text()) != manifest:
                        raise
                    verify_samples(target, manifest)
        return {"manifest": str(target / "manifest.json"), "attempt": attempt, "trial_id": trial_id}
    except Exception as exc:
        raise EventTransferError(f"Event output transfer failed for {trial_id}: {exc}") from exc

# Publish and restore icetune outputs from the dedicated Ray head
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import hashlib
import io
import json
import math
import os
import pathlib
import re
import shutil
import sys
import tarfile
import tempfile

from core.io.files import ensure_dir

from submit.lxplus import common

# Keep the submit process independent of the icetune Python environment
HISTORY_FIGURES = (
    "cost_evolution.png",
    "cost_evolution_logy.png",
    "cost_sorted.png",
    "cost_sorted_logy.png",
    "parameter_evolution_physical.pdf",
    "parameter_evolution_optimizer.pdf",
)


# Publish a complete output only after its bytes reach the filesystem
def atomic_write(target: pathlib.Path, payload: bytes) -> None:
    ensure_dir(target.parent)
    descriptor, name = tempfile.mkstemp(prefix=f".{target.stem}-", suffix=".tmp", dir=target.parent)
    stage = pathlib.Path(name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        stage.replace(target)
    finally:
        stage.unlink(missing_ok=True)


# Compute whether one Tune-relative path is a transported worker output
def is_ray_trial_output(path: pathlib.PurePath) -> bool:
    parts = path.parts
    for index, part in enumerate(parts):
        if part == "figs" and index + 1 < len(parts) and parts[index + 1] == "icetune":
            return True
    return bool(len(parts) >= 2 and parts[-2] == "results" and re.fullmatch(r"TUNE_icetune_[^/]+\.pkl", parts[-1]))


# Compute a content revision for one compact output tree
def output_tree_revision(
    path: pathlib.Path, *, exclude_top: set[str] | None = None, exclude_files: set[str] | None = None) -> str | None:
    if not path.is_dir():
        return None
    excluded = set(exclude_top or set())
    excluded_files = set(exclude_files or set())
    digest = hashlib.sha256()
    found = False
    for item in sorted(path.rglob("*")):
        if not item.is_file():
            continue
        relative = item.relative_to(path)
        if ((relative.parts and relative.parts[0] in excluded) or relative.as_posix() in excluded_files
            or is_ray_trial_output(relative)):
            continue
        found = True
        stat = item.stat()
        digest.update(relative.as_posix().encode("utf-8"))
        digest.update(f"\0{stat.st_size}\0{stat.st_mtime_ns}\0".encode("ascii"))
    return digest.hexdigest() if found else None


# Compress one local output tree with a fixed archive root
def pack_output_tree(source: pathlib.Path, archive_root: pathlib.Path, *, exclude_top: set[str] | None = None) -> bytes:
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as archive:
        excluded = set(exclude_top or set())

        # Drop worker outputs because the incremental immutable stream owns them
        def include(info: tarfile.TarInfo) -> tarfile.TarInfo | None:
            relative = pathlib.PurePosixPath(info.name).relative_to(archive_root)
            if relative.parts and relative.parts[0] in excluded:
                return None
            return None if is_ray_trial_output(relative) else info

        archive.add(source, arcname=archive_root.as_posix(), recursive=True, filter=include)
    return output.getvalue()


# Compute one compact revision for a readable live state file
def state_revision(path: pathlib.Path) -> tuple[int, int] | None:
    try:
        if not path.is_file():
            return None
        stat = path.stat()
    except OSError:
        return None
    return stat.st_size, stat.st_mtime_ns


# Read one file only when its size and modification time remain stable
def read_stable_file(path: pathlib.Path, revision: tuple[int, int]) -> bytes:
    payload = path.read_bytes()
    if state_revision(path) != revision:
        raise RuntimeError(f"Ray lxplus output changed while reading: {path}")
    return payload


# Package changed figures, live state and new full trial pickles
def pack_runtime_outputs(*, figure_dir: pathlib.Path | None, result_files: list[pathlib.Path], run_name: str,
    state_files: list[pathlib.Path] | None = None,
    state_payloads: dict[pathlib.Path, tuple[bytes, tuple[int, int]]] | None = None,
    state_root: pathlib.Path | None = None, exclude_figures: tuple[str, ...] = (), ) -> bytes | None:
    state_files = list(state_files or [])
    state_payloads = dict(state_payloads or {})
    if (figure_dir is None or not figure_dir.is_dir()) and not result_files and not state_files and not state_payloads:
        return None
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as archive:
        if figure_dir is not None and figure_dir.is_dir():
            root = pathlib.PurePosixPath("figs", "icetune", run_name)

            # Exclude independently published history files only at the figure root
            def include(info: tarfile.TarInfo) -> tarfile.TarInfo | None:
                return None if pathlib.PurePosixPath(info.name).relative_to(root).as_posix() in exclude_figures else info

            archive.add(figure_dir, arcname=root.as_posix(), recursive=True, filter=include)
        for path in result_files:
            archive.add(path, arcname=(pathlib.Path("runs") / "icetune" / run_name / "results" / path.name).as_posix(),
                recursive=False, )
        for path in state_files:
            relative = path.relative_to(state_root) if state_root is not None else pathlib.Path(path.name)
            archive.add(
                path, arcname=(pathlib.Path("runs") / "icetune" / run_name / relative).as_posix(), recursive=False)
        for relative, (payload, revision) in sorted(state_payloads.items(), key=lambda item: item[0].as_posix()):
            member = pathlib.Path("runs") / "icetune" / run_name / relative
            info = tarfile.TarInfo(member.as_posix())
            info.mode = 0o644
            info.mtime = revision[1] / 1_000_000_000
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
    return output.getvalue()


# Compute completed full-trial descriptors only when figures match Ray's checkpoint
def recovery_result_descriptors(
    *, checkpoint: bytes, state_dir: pathlib.Path, figure_dir: pathlib.Path, run_name: str, cost: str, plot: bool
) -> list[dict] | None:
    try:
        with tarfile.open(fileobj=io.BytesIO(checkpoint), mode="r:gz") as archive:
            state_names = sorted(name for name in archive.getnames()
                if re.fullmatch(rf"{re.escape(run_name)}/experiment_state-[^/]+\.json", name))
            if not state_names:
                return None
            state = json.loads(archive.extractfile(state_names[-1]).read())
        completed_ids = set()
        completed = []
        descriptors = {}
        for item in state["trial_data"]:
            if not isinstance(item, list | tuple) or len(item) != 2:
                return None
            trial = json.loads(item[0])
            runtime = json.loads(item[1])
            if trial.get("status") != "TERMINATED":
                continue
            trial_id = str(trial.get("trial_id") or "")
            if not trial_id:
                return None
            completed_ids.add(trial_id)
            result = runtime.get("last_result") or {}
            result_cost = result.get(cost)
            if isinstance(result_cost, int | float) and math.isfinite(result_cost):
                completed.append({"cost": float(result_cost), "initial": bool(result.get("is_initial")),
                        "theta_hash": result.get("theta_hash"), "trial_id": trial_id, })
            descriptor = result.get("trial_pickle")
            if descriptor is None:
                continue
            if not isinstance(descriptor, dict):
                return None
            name = str(descriptor.get("filename") or "")
            checksum = str(descriptor.get("sha256") or "")
            size = descriptor.get("size")
            path = state_dir / "results" / name
            if (re.fullmatch(r"TUNE_icetune_[^/]+\.pkl", name) is None
                or re.fullmatch(r"[0-9a-f]{64}", checksum) is None or not isinstance(size, int) or size < 0
                or not path.is_file() or path.stat().st_size != size or common.file_sha256(path) != checksum):
                return None
            selected = {"filename": name, "sha256": checksum, "size": size}
            if name in descriptors and descriptors[name] != selected:
                return None
            descriptors[name] = selected
        if plot and completed:
            minimum = min(item["cost"] for item in completed)
            summary_path = figure_dir / "summary.json"
            if not summary_path.is_file():
                return None
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            summary_cost = (summary.get("metrics") or {}).get(cost)
            summary_id = str(summary.get("trial_id") or "")
            best = next((item for item in completed if item["trial_id"] == summary_id), None)
            if (best is None or not isinstance(summary_cost, int | float)
                or not math.isclose(float(summary_cost), best["cost"], rel_tol=1e-12, abs_tol=1e-12)
                or not math.isclose(best["cost"], minimum, rel_tol=1e-12, abs_tol=1e-12)
                or (best["theta_hash"] is not None and summary.get("theta_hash") != best["theta_hash"])):
                return None
            initial = next((item for item in completed if item["initial"]), None)
            if initial is not None:
                initial_path = figure_dir / "init" / "summary.json"
                if not initial_path.is_file():
                    return None
                initial_summary = json.loads(initial_path.read_text(encoding="utf-8"))
                if str(initial_summary.get("trial_id") or "") != initial["trial_id"]:
                    return None
        elif any(str((json.loads(path.read_text(encoding="utf-8")) or {}).get("trial_id") or "") not in completed_ids
            for path in (figure_dir / "summary.json", figure_dir / "init" / "summary.json") if path.is_file()):
            return None
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return None
    return [descriptors[name] for name in sorted(descriptors)]


# Read changed live JSON files without coupling their publication to figures
def snapshot_live_state(*, state_dir: pathlib.Path, run_name: str, live_revision: dict[str, tuple[int, int]] | None
) -> tuple[bytes | None, dict[str, tuple[int, int]], dict, list[str]]:
    paths = [state_dir / name for name in ("status.json", "history.json", "summary.json")]
    paths.extend(sorted((state_dir / "failures").glob("*.json")))
    previous = dict(live_revision or {})
    accepted = dict(previous)
    current = {}
    payloads = {}
    deferred = []
    status_payload = None
    for path in paths:
        relative = path.relative_to(state_dir).as_posix()
        revision = state_revision(path)
        if revision is None:
            continue
        current[relative] = revision
        payload = None
        if revision != previous.get(relative):
            try:
                payload = read_stable_file(path, revision)
                state = json.loads(payload)
                if not isinstance(state, dict):
                    raise ValueError(f"Ray lxplus JSON is not a dictionary: {path}")
            except (OSError, RuntimeError, TypeError, ValueError, json.JSONDecodeError):
                deferred.append(relative)
                continue
            payloads[pathlib.Path(relative)] = (payload, revision)
            accepted[relative] = revision
        if relative == "status.json":
            status_payload = payload
    accepted = {name: revision for name, revision in accepted.items() if name in current}
    status_state = {}
    status_path = state_dir / "status.json"
    status_revision = current.get("status.json")
    if status_payload is None and status_revision is not None:
        try:
            status_payload = read_stable_file(status_path, status_revision)
        except (OSError, RuntimeError):
            status_payload = None
    if status_payload is not None:
        try:
            status_state = json.loads(status_payload)
        except (TypeError, ValueError, json.JSONDecodeError):
            status_state = {}
        if not isinstance(status_state, dict):
            status_state = {}
    output = pack_runtime_outputs(figure_dir=None, result_files=[], run_name=run_name, state_payloads=payloads)
    return output, accepted, status_state, deferred


# Package only changed figure data and isolate concurrent writes from live JSON
def snapshot_live_figures(*, figure_dir: pathlib.Path, run_name: str, figure_revision: dict | None
) -> tuple[bytes | None, dict | None, dict, list[str]]:
    previous = dict(figure_revision or {})
    accepted = dict(previous)
    deferred = []
    output = None
    history = {}
    try:
        best_revision = output_tree_revision(figure_dir, exclude_files=set(HISTORY_FIGURES))
    except OSError:
        best_revision = None
        deferred.append("figures")
    if best_revision is None:
        accepted.pop("best", None)
    elif best_revision != previous.get("best"):
        try:
            output = pack_runtime_outputs(figure_dir=figure_dir, result_files=[], run_name=run_name,
                                          exclude_figures=HISTORY_FIGURES)
            if output_tree_revision(figure_dir, exclude_files=set(HISTORY_FIGURES)) != best_revision:
                raise RuntimeError(f"Ray lxplus figures changed while snapshotting: {figure_dir}")
        except (OSError, RuntimeError):
            output = None
            deferred.append("figures")
        else:
            accepted["best"] = best_revision
    previous_history = previous.get("history", {})
    accepted["history"] = dict(previous_history)
    for name in HISTORY_FIGURES:
        path = figure_dir / name
        revision = state_revision(path)
        if revision is None:
            accepted["history"].pop(name, None)
        elif revision != previous_history.get(name):
            try:
                history[name] = read_stable_file(path, revision)
            except (OSError, RuntimeError):
                deferred.append(name)
            else:
                accepted["history"][name] = revision
    return output, history or None, accepted, deferred


# Package changed live outputs without reading mutable Tune checkpoint state
def snapshot_live_outputs(*, environment: dict[str, str], figure_revision: dict | None,
    live_revision: dict[str, tuple[int, int]] | None = None, ) -> dict:
    runtime = common.runtime_work_root(environment) / "runtime"
    run_name = str(environment["RUN_NAME"])
    figure_dir = runtime / "figs" / "icetune" / run_name
    state_dir = runtime / "runs" / "icetune" / run_name
    outputs, live_revision, status, live_deferred = snapshot_live_state(
        state_dir=state_dir, run_name=run_name, live_revision=live_revision)
    figure_outputs, history_figures, figure_revision, figure_deferred = snapshot_live_figures(
        figure_dir=figure_dir, run_name=run_name, figure_revision=figure_revision)
    return {"deferred": live_deferred + figure_deferred, "figure_outputs": figure_outputs,
        "figure_revision": figure_revision, "history_figures": history_figures, "live_revision": live_revision,
        "outputs": outputs, "ray_resources": status.get("ray_resources"),
        "ray_resources_observed_at_unix": status.get("ray_resources_observed_at_unix", status.get("updated_at_unix")), }


# Package one final checkpoint only after the icetune driver has exited
def snapshot_final_outputs(*, environment: dict[str, str]) -> dict:
    pid_path = common.runtime_pid_path(environment)
    if pid_path.exists():
        raise RuntimeError(f"Refusing to checkpoint while the icetune driver record exists: {pid_path}")
    runtime = common.runtime_work_root(environment) / "runtime"
    run_name = str(environment["RUN_NAME"])
    figure_dir = runtime / "figs" / "icetune" / run_name
    state_dir = runtime / "runs" / "icetune" / run_name
    live_files = [state_dir / name for name in ("status.json", "history.json", "summary.json")]
    live_files.extend(sorted((state_dir / "failures").glob("*.json")))
    live_files = [path for path in live_files if path.is_file()]
    checkpoint = None
    recovery = None
    results = []
    if (state_dir / "tuner.pkl").is_file():
        state_revision_before = output_tree_revision(state_dir, exclude_top={"results"})
        checkpoint = pack_output_tree(state_dir, pathlib.Path(run_name), exclude_top={"results"})
        if output_tree_revision(state_dir, exclude_top={"results"}) != state_revision_before:
            raise RuntimeError(f"Ray lxplus final state changed while packing: {state_dir}")
        descriptors = recovery_result_descriptors(checkpoint=checkpoint, state_dir=state_dir,
            figure_dir=figure_dir, run_name=run_name, cost=str(environment.get("COST", "loss")),
            plot=str(environment.get("PLOT", "0")) == "1", )
        if descriptors is not None:
            names = {str(item["filename"]) for item in descriptors}
            result_dir = state_dir / "results"
            results = [path for path in sorted(result_dir.glob("TUNE_icetune_*.pkl")) if path.name in names]
            recovery_outputs = pack_runtime_outputs(figure_dir=figure_dir if figure_dir.is_dir() else None,
                result_files=[], run_name=run_name, state_files=live_files, state_root=state_dir, )
            if recovery_outputs is not None:
                recovery = pack_recovery_generation(
                    checkpoint=checkpoint, outputs=recovery_outputs, result_descriptors=descriptors, run_name=run_name)
    if environment.get("ALGORITHM") == "lbfgs":
        results = sorted((state_dir / "results").glob("TUNE_icetune_*.pkl"))
    outputs = pack_runtime_outputs(figure_dir=figure_dir if figure_dir.is_dir() else None, result_files=results,
        run_name=run_name, state_files=live_files, state_root=state_dir, )
    return {"checkpoint": checkpoint, "history_figures": None, "outputs": outputs, "recovery": recovery,
        "result_names": [path.name for path in results], }


# Compute the external run directory for one immutable tuning campaign
def campaign_run_path(shared_root: pathlib.Path, run_name: str, campaign_fingerprint: str) -> pathlib.Path:
    if re.fullmatch(r"[0-9a-f]{64}", campaign_fingerprint) is None:
        raise ValueError("Ray campaign fingerprint must be lowercase hexadecimal")
    return shared_root / "runs" / "icetune" / run_name / "campaigns" / campaign_fingerprint[:16]


# Compute the external figure directory for one immutable tuning campaign
def campaign_figure_path(shared_root: pathlib.Path, run_name: str, campaign_fingerprint: str) -> pathlib.Path:
    if re.fullmatch(r"[0-9a-f]{64}", campaign_fingerprint) is None:
        raise ValueError("Ray campaign fingerprint must be lowercase hexadecimal")
    return shared_root / "figs" / "icetune" / run_name / "campaigns" / campaign_fingerprint[:16]


# Compute the single persistent Ray Tune checkpoint path
def checkpoint_path(shared_root: pathlib.Path, run_name: str, campaign_fingerprint: str) -> pathlib.Path:
    return campaign_run_path(shared_root, run_name, campaign_fingerprint) / "checkpoint.tar.gz"


# Compute the atomic restart generation containing state and selected outputs
def recovery_path(shared_root: pathlib.Path, run_name: str, campaign_fingerprint: str) -> pathlib.Path:
    return campaign_run_path(shared_root, run_name, campaign_fingerprint) / "recovery.tar.gz"


# Append one in-memory regular file to a tar stream
def add_tar_bytes(archive: tarfile.TarFile, name: str, payload: bytes) -> None:
    info = tarfile.TarInfo(name=name)
    info.size = len(payload)
    info.mode = 0o600
    archive.addfile(info, io.BytesIO(payload))


# Package one checkpoint and its matching compact outputs as one generation
def pack_recovery_generation(*, checkpoint: bytes, outputs: bytes, result_descriptors: list[dict], run_name: str
) -> bytes:
    manifest = {"checkpoint_sha256": hashlib.sha256(checkpoint).hexdigest(), "checkpoint_size": len(checkpoint),
        "outputs_sha256": hashlib.sha256(outputs).hexdigest(), "outputs_size": len(outputs),
        "results": result_descriptors, "run_name": run_name, "schema_version": 1, }
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w:gz") as archive:
        add_tar_bytes(archive, "manifest.json", json.dumps(manifest, indent=2, sort_keys=True).encode("utf-8"))
        add_tar_bytes(archive, "checkpoint.tar.gz", checkpoint)
        add_tar_bytes(archive, "outputs.tar.gz", outputs)
    return payload.getvalue()


# Validate and unpack one complete recovery generation
def unpack_recovery_generation(payload: bytes, *, run_name: str) -> tuple[bytes, bytes, list[dict]]:
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
        members = archive.getmembers()
        if {member.name for member in members} != {"checkpoint.tar.gz", "manifest.json", "outputs.tar.gz"} or any(
            not member.isfile() for member in members):
            raise ValueError("Ray lxplus recovery generation has unknown contents")
        manifest = json.loads(archive.extractfile("manifest.json").read())
        checkpoint = archive.extractfile("checkpoint.tar.gz").read()
        outputs = archive.extractfile("outputs.tar.gz").read()
    if (manifest.get("schema_version") != 1 or manifest.get("run_name") != run_name
        or manifest.get("checkpoint_size") != len(checkpoint)
        or manifest.get("checkpoint_sha256") != hashlib.sha256(checkpoint).hexdigest()
        or manifest.get("outputs_size") != len(outputs)
        or manifest.get("outputs_sha256") != hashlib.sha256(outputs).hexdigest()
        or not isinstance(manifest.get("results"), list)):
        raise ValueError("Ray lxplus recovery generation manifest is invalid")
    return checkpoint, outputs, manifest["results"]


# Validate and extract one compressed Ray Tune checkpoint into local scratch
def extract_checkpoint(*, checkpoint: pathlib.Path, runtime: pathlib.Path, run_name: str) -> None:
    destination = runtime / "runs" / "icetune"
    ensure_dir(destination)
    with tarfile.open(checkpoint, mode="r:gz") as archive:
        members = archive.getmembers()
        for member in members:
            parts = pathlib.PurePosixPath(member.name).parts
            if not parts or parts[0] != run_name:
                raise ValueError(f"Unknown Ray lxplus checkpoint path: {member.name}")
        if f"{run_name}/tuner.pkl" not in {member.name for member in members}:
            raise ValueError(f"Ray lxplus checkpoint has no tuner.pkl: {checkpoint}")
        archive.extractall(destination, members=members, filter="data")


# Validate one immutable full-trial file referenced by a recovery generation
def recovery_result_path(*, shared_root: pathlib.Path, run_name: str, campaign_fingerprint: str, descriptor: dict
) -> pathlib.Path:
    if not isinstance(descriptor, dict):
        raise ValueError("Ray recovery result descriptor must be a dictionary")
    name = str(descriptor.get("filename") or "")
    checksum = str(descriptor.get("sha256") or "")
    size = descriptor.get("size")
    if (re.fullmatch(r"TUNE_icetune_[^/]+\.pkl", name) is None or re.fullmatch(r"[0-9a-f]{64}", checksum) is None
        or not isinstance(size, int) or size < 0):
        raise ValueError("Ray recovery result descriptor is invalid")
    source = campaign_run_path(shared_root, run_name, campaign_fingerprint) / "results" / name
    if not source.is_file() or source.stat().st_size != size or common.file_sha256(source) != checksum:
        raise RuntimeError(f"Ray recovery result is missing or incomplete: {source}")
    return source


# Restore one complete atomic recovery generation into head-local scratch
def stage_recovery_generation(*, recovery: pathlib.Path, shared_root: pathlib.Path, runtime: pathlib.Path,
    run_name: str, campaign_fingerprint: str, ) -> None:
    checkpoint, outputs, descriptors = unpack_recovery_generation(recovery.read_bytes(), run_name=run_name)
    temporary_root = runtime / "tmp"
    ensure_dir(temporary_root)
    descriptor, checkpoint_name = tempfile.mkstemp(prefix="recovery-checkpoint-", suffix=".tar.gz", dir=temporary_root)
    checkpoint_file = pathlib.Path(checkpoint_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(checkpoint)
        extract_checkpoint(checkpoint=checkpoint_file, runtime=runtime, run_name=run_name)
    finally:
        checkpoint_file.unlink(missing_ok=True)
    extract_runtime_outputs(payload=outputs, runtime=runtime, run_name=run_name)
    result_target = runtime / "runs" / "icetune" / run_name / "results"
    ensure_dir(result_target)
    for item in descriptors:
        source = recovery_result_path(
            shared_root=shared_root, run_name=run_name, campaign_fingerprint=campaign_fingerprint, descriptor=item)
        shutil.copy2(source, result_target / source.name)
    restore_outputs(
        payload=outputs, shared_root=shared_root, run_name=run_name, campaign_fingerprint=campaign_fingerprint)


# Restore prior Tune state from one checkpoint file into local scratch
def stage_previous_outputs(
    *, shared_root: pathlib.Path, runtime: pathlib.Path, run_name: str, campaign_fingerprint: str
) -> list[pathlib.Path]:
    recovery = recovery_path(shared_root, run_name, campaign_fingerprint)
    if recovery.is_file():
        stage_recovery_generation(recovery=recovery, shared_root=shared_root, runtime=runtime, run_name=run_name,
            campaign_fingerprint=campaign_fingerprint, )
        return [recovery]
    checkpoint = checkpoint_path(shared_root, run_name, campaign_fingerprint)
    if checkpoint.is_file():
        extract_checkpoint(checkpoint=checkpoint, runtime=runtime, run_name=run_name)
        common.write_log("Restoring checkpoint without loose EOS figures or trial files because "
            "no matching recovery generation exists", file=sys.stderr, )
        return [checkpoint]

    return []


# Atomically publish one returned output directory with immediate rollback
def publish_output_dir(*, source: pathlib.Path, target: pathlib.Path) -> None:
    ensure_dir(target.parent)
    with tempfile.TemporaryDirectory(prefix=".figures-", dir=target.parent) as directory:
        stage = pathlib.Path(directory)
        prepared = stage / "new"
        backup = stage / "previous"
        shutil.copytree(source, prepared)
        for name in HISTORY_FIGURES:
            if (target / name).is_file() and not (prepared / name).exists():
                shutil.copy2(target / name, prepared / name)
        if target.exists():
            target.replace(backup)
        try:
            prepared.replace(target)
        except Exception:
            if backup.exists() and not target.exists():
                backup.replace(target)
            raise


# Atomically publish each independently changing Ray history figure
def publish_history_figures(*, payload: dict[str, bytes], shared_root: pathlib.Path, run_name: str,
    campaign_fingerprint: str) -> None:
    """Publish changed history plots without replacing the best trial figures"""
    for name, data in payload.items():
        signature = b"%PDF-" if name.endswith(".pdf") else b"\x89PNG\r\n\x1a\n"
        if name not in HISTORY_FIGURES or not isinstance(data, bytes) or not data.startswith(signature):
            raise ValueError(f"Invalid Ray lxplus history figure: {name}")
    target = campaign_figure_path(shared_root, run_name, campaign_fingerprint)
    ensure_dir(target)
    for name, data in payload.items():
        atomic_write(target / name, data)


# Publish immutable trial histograms and data covariance for iceproxy
def publish_result_files(*, source: pathlib.Path, target: pathlib.Path) -> None:
    if not source.is_dir():
        return
    ensure_dir(target)
    paths = sorted(source.glob("TUNE_icetune_*.pkl"))
    paths.extend(path for name in ("data_covariance.npz", "data_covariance.json") if (path := source / name).is_file())
    for path in paths:
        destination = target / path.name
        if destination.exists():
            if path.stat().st_size != destination.stat().st_size or common.file_sha256(path) != common.file_sha256(
                destination):
                raise RuntimeError(f"Conflicting Ray lxplus histogram or covariance file: {destination}") from None
        else:
            atomic_write(destination, path.read_bytes())


# Publish full Ray trial failure logs into the canonical run directory
def publish_failure_files(*, source: pathlib.Path, target: pathlib.Path) -> None:
    if not source.is_dir():
        return
    ensure_dir(target)
    for path in sorted(source.glob("*.json")):
        destination = target / path.name
        if destination.is_file():
            if path.read_bytes() != destination.read_bytes():
                raise RuntimeError(f"Conflicting Ray lxplus failure log: {destination}")
            continue
        atomic_write(destination, path.read_bytes())


# Validate compact returned-output archive members for one campaign
def validate_output_members(members: list[tarfile.TarInfo], *, run_name: str) -> None:
    for member in members:
        parts = pathlib.PurePosixPath(member.name).parts
        valid_figure = (parts and parts[0] == "figs" and (len(parts) < 2 or parts[1] == "icetune")
            and (len(parts) < 3 or parts[2] == run_name))
        valid_result = (parts and parts[0] == "runs" and (len(parts) < 2 or parts[1] == "icetune")
            and (len(parts) < 3 or parts[2] == run_name) and (len(parts) < 4 or parts[3] == "results")
            and len(parts) <= 5 and (len(parts) < 5 or re.fullmatch(r"TUNE_icetune_[^/]+\.pkl", parts[4])
                or parts[4] in {"data_covariance.npz", "data_covariance.json"}))
        valid_state = (len(parts) == 4 and parts[:3] == ("runs", "icetune", run_name)
            and parts[3] in {"history.json", "status.json", "summary.json"})
        valid_failure = (len(parts) == 5 and parts[:4] == ("runs", "icetune", run_name, "failures")
            and re.fullmatch(r"[A-Za-z0-9_.-]+\.json", parts[4]) is not None)
        if not valid_figure and not valid_result and not valid_state and not valid_failure:
            raise ValueError(f"Unknown Ray lxplus output path: {member.name}")


# Extract one validated compact output archive below a local runtime
def extract_runtime_outputs(*, payload: bytes, runtime: pathlib.Path, run_name: str) -> None:
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
        members = archive.getmembers()
        validate_output_members(members, run_name=run_name)
        archive.extractall(runtime, members=members, filter="data")


# Restore Tune results and figures below the shared GRANIITTI output trees
def restore_outputs(*, payload: bytes, shared_root: pathlib.Path, run_name: str, campaign_fingerprint: str) -> None:
    run_target = campaign_run_path(shared_root, run_name, campaign_fingerprint)
    ensure_dir(run_target)
    with tempfile.TemporaryDirectory(prefix="return-", dir=run_target) as stage_name:
        stage = pathlib.Path(stage_name)
        extract_runtime_outputs(payload=payload, runtime=stage, run_name=run_name)
        source = stage / "figs" / "icetune" / run_name
        if source.is_dir():
            publish_output_dir(
                source=source, target=campaign_figure_path(shared_root, run_name, campaign_fingerprint))
        publish_result_files(source=stage / "runs" / "icetune" / run_name / "results", target=run_target / "results")
        publish_failure_files(source=stage / "runs" / "icetune" / run_name / "failures", target=run_target / "failures")
        for name in ("status.json", "history.json", "summary.json"):
            source = stage / "runs" / "icetune" / run_name / name
            if source.is_file():
                target = run_target / name
                ensure_dir(target.parent)
                atomic_write(target, source.read_bytes())


# Atomically publish one compressed Ray Tune checkpoint file
def publish_checkpoint(*, payload: bytes, shared_root: pathlib.Path, run_name: str, campaign_fingerprint: str) -> None:
    target = checkpoint_path(shared_root, run_name, campaign_fingerprint)
    ensure_dir(target.parent)
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
        members = archive.getmembers()
        names = {member.name for member in members}
        if f"{run_name}/tuner.pkl" not in names:
            raise ValueError("Ray lxplus checkpoint payload has no tuner.pkl")
        for member in members:
            parts = pathlib.PurePosixPath(member.name).parts
            if not parts or parts[0] != run_name:
                raise ValueError(f"Unknown Ray lxplus checkpoint path: {member.name}")
    atomic_write(target, payload)


# Atomically publish one validated restart generation after its result files
def publish_recovery_generation(
    *, payload: bytes, shared_root: pathlib.Path, run_name: str, campaign_fingerprint: str) -> None:
    _, _, descriptors = unpack_recovery_generation(payload, run_name=run_name)
    for item in descriptors:
        recovery_result_path(
            shared_root=shared_root, run_name=run_name, campaign_fingerprint=campaign_fingerprint, descriptor=item)
    target = recovery_path(shared_root, run_name, campaign_fingerprint)
    ensure_dir(target.parent)
    atomic_write(target, payload)


# Fetch and restore changed live outputs from the Ray head worker
def return_live_outputs(*, client, head_key: str, environment: dict[str, str], shared_root: pathlib.Path,
    run_name: str, campaign_fingerprint: str, figure_revision: dict | None,
    live_revision: dict[str, tuple[int, int]] | None, timeout_s: float | None = None,
) -> tuple[dict | None, dict[str, tuple[int, int]], set[str], dict | None]:
    future = client.submit(snapshot_live_outputs, environment=environment, figure_revision=figure_revision,
        live_revision=live_revision, workers=[head_key], allow_other_workers=False, pure=False, )
    try:
        snapshot = future.result() if timeout_s is None else future.result(timeout=timeout_s)
    except BaseException:
        future.cancel()
        raise
    if not isinstance(snapshot, dict):
        raise RuntimeError("Ray head returned no live output snapshot")
    deferred = sorted(set(snapshot.get("deferred") or []))
    if deferred:
        common.write_log(
            f"Ray lxplus live output changed during reading and will retry: {', '.join(deferred)}", file=sys.stderr)
    for key, publish in (("outputs", restore_outputs), ("figure_outputs", restore_outputs),
                         ("history_figures", publish_history_figures)):
        if snapshot.get(key) is not None:
            publish(payload=snapshot[key], shared_root=shared_root, run_name=run_name,
                    campaign_fingerprint=campaign_fingerprint)
    previous_figure_revision = figure_revision if isinstance(figure_revision, dict) else {}
    current_figure_revision = snapshot.get("figure_revision")
    figure_updates = set()
    if isinstance(current_figure_revision, dict):
        for name in ("best", "history"):
            if current_figure_revision.get(name) is not None and current_figure_revision.get(name
            ) != previous_figure_revision.get(name):
                figure_updates.add(name)
    return (current_figure_revision, snapshot.get("live_revision"), figure_updates,
        {"observed_at": snapshot.get("ray_resources_observed_at_unix"), "resources": snapshot.get("ray_resources")}
        if snapshot.get("ray_resources") is not None else None, )


# Publish the complete result returned by the head local icetune driver
def publish_runtime_result(*, result: dict, shared_root: pathlib.Path, run_name: str, campaign_fingerprint: str) -> int:
    if not isinstance(result, dict):
        raise RuntimeError("Ray head returned no final runtime result")
    for key, publish in (("outputs", restore_outputs), ("history_figures", publish_history_figures),
                         ("recovery", publish_recovery_generation), ("checkpoint", publish_checkpoint)):
        if result.get(key) is not None:
            publish(payload=result[key], shared_root=shared_root, run_name=run_name,
                    campaign_fingerprint=campaign_fingerprint)
    return int(result["returncode"])

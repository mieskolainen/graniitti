# Immutable Pandora reconstruction inputs and worker staging
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import fcntl
import hashlib
import json
import os
import shutil
import tarfile
import uuid
import xml.etree.ElementTree as ET
from pathlib import Path

from core.io.files import check_dependencies, ensure_dir, file_identity
from core.tune import short_id


# Copy and verify one file atomically under a process shared lock
def stage_file(source, target, identity):
    target = Path(target)
    ensure_dir(target.parent)
    with target.with_suffix(target.suffix + ".lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not target.exists():
            pending = target.with_name(target.name + ".pending." + uuid.uuid4().hex)
            shutil.copyfile(source, pending)
            if file_identity(pending) != identity:
                raise ValueError(f"Pandora input checksum mismatch: {source}")
            pending.chmod(0o444)
            pending.replace(target)
        status = target.stat()
        stamp = json.dumps([identity, status.st_size, status.st_mtime_ns, status.st_ctime_ns, status.st_ino])
        lock.seek(0)
        if lock.read() != stamp:
            if file_identity(target) != identity:
                raise ValueError(f"Pandora checksum mismatch or corrupt cached file: {target}")
            lock.seek(0)
            lock.truncate()
            lock.write(stamp)
            lock.flush()
    return target


# Isolate reconstruction from inherited runtime paths and bound numerical threads
def command_environment(*, inherit=True):
    environment = {key: os.environ[key] for key in ("HOME", "USER", "LOGNAME", "TMPDIR") if inherit and key in os.environ}
    return dict(environment, PATH=os.defpath, **{key: "1" for key in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")})


# Freeze external files into a deterministic reconstruction archive
def freeze_runtime(root, paths, output_dir, description):
    import io

    root, output_dir = Path(root).resolve(), Path(output_dir)
    check_dependencies(description["dependencies"])
    files = set()
    for pattern in paths:
        matches = list(root.glob(pattern))
        if not matches:
            raise FileNotFoundError(f"Pandora runtime input not found: {pattern}")
        for path in matches:
            files.update(p for p in (path.rglob("*") if path.is_dir() else [path]) if p.is_file())
    ensure_dir(output_dir)
    pending = output_dir / (uuid.uuid4().hex + ".pending.tar")
    manifest = {}
    with tarfile.open(pending, "w", dereference=True) as archive:
        for path in sorted(files):
            relative = path.relative_to(root).as_posix()
            identity = file_identity(path)
            info = archive.gettarinfo(str(path), arcname=relative)
            info.uid = info.gid = info.mtime = 0
            info.uname = info.gname = ""
            info.mode = 0o555 if os.access(path, os.X_OK) else 0o444
            with path.open("rb") as stream:
                archive.addfile(info, stream)
            if file_identity(path) != identity:
                raise ValueError(f"Pandora runtime changed while freezing: {path}")
            manifest[relative] = identity
        content = description["setup"].encode()
        info = tarfile.TarInfo("setup.sh")
        info.size, info.mode = len(content), 0o444
        archive.addfile(info, io.BytesIO(content))
        manifest["setup.sh"] = dict(sha256=hashlib.sha256(content).hexdigest(), size=len(content))
    identity = file_identity(pending)
    target = output_dir / (short_id(identity["sha256"]) + ".tar")
    with target.with_suffix(".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if target.exists() and file_identity(target) != identity:
            raise ValueError(f"Conflicting Pandora runtime archive: {target}")
        pending.replace(target)
    target.chmod(0o444)
    return dict(**identity, files=manifest, dependencies=description["dependencies"]), target


# Expand one verified archive once per worker and check every extracted input
def stage_runtime(manifest, cdir, cache):
    identity = {key: manifest[key] for key in ("sha256", "size")}
    archive = stage_file(Path(cdir) / manifest["archive"], cache / (short_id(manifest["sha256"]) + ".tar"), identity)
    root = cache / short_id(manifest["sha256"])
    with root.with_suffix(".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not root.exists():
            pending = root.with_name(root.name + ".pending." + uuid.uuid4().hex)
            with tarfile.open(archive) as source:
                source.extractall(pending, filter="data")
            for relative, expected in manifest["files"].items():
                if file_identity(pending / relative) != expected:
                    raise ValueError(f"Corrupt Pandora runtime input: {relative}")
            pending.replace(root)
        else:
            for relative, expected in manifest["files"].items():
                stage_file(root / relative, root / relative, expected)
    check_dependencies(manifest["dependencies"])
    return root


# Freeze the selected runtime and input identities before saving a campaign
def prepare_tunesetup(tunesetup, cdir):
    from core.tune import push as icetune_push
    from core.tune.drivers.pandora.driver import _push_catalog, _xml_targets

    cards = tunesetup.datacards
    code = {
        name: file_identity(Path(__file__).parent / name) for name in ("driver.py", "runtime.py", "tunesetup/common.py")
    }
    if all("runtime" in card for card in cards):
        if any(card["runtime"]["driver_files"] != code for card in cards):
            raise ValueError("Pandora driver changed after preflight, start a new campaign")
        for card in cards:
            manifest = card["runtime"]
            if file_identity(Path(cdir) / manifest["archive"]) != {key: manifest[key] for key in ("sha256", "size")}:
                raise ValueError("Frozen Pandora runtime checksum mismatch")
        return
    if any("runtime" in card for card in cards):
        raise ValueError("Pandora campaign contains a mixture of frozen and unfrozen runtimes")
    for card in cards:
        card["inputs"] = [dict(path=str(Path(path).resolve()), **file_identity(path)) for path in card["input_roots"]]
    root = Path(cards[0]["pandora_dir"]).resolve()
    if any(Path(card["pandora_dir"]).resolve() != root for card in cards):
        raise ValueError("Pandora datasets must share one reconstruction installation")
    catalog = tunesetup.aux_param_space["pandora_catalog"]
    expected_xml = catalog["settings_xml_sha256"]
    if any(file_identity(card["settings_template"])["sha256"] != expected_xml for card in cards):
        raise ValueError("Pandora XML changed after parameter catalog construction")
    from core.tune.drivers.pandora import key4hep

    key4hep.prepare(tunesetup=tunesetup, environment=os.environ)
    description = tunesetup.aux_param_space.pop("pandora_runtime")
    manifest, archive = freeze_runtime(root, tunesetup.aux_param_space["pandora_paths"]["runtime_paths"],
                                       Path(cdir) / "tmp/icetune/pandora", description)
    for card in cards:
        relative_xml = Path(card["settings_template"]).relative_to(root).as_posix()
        if manifest["files"][relative_xml]["sha256"] != expected_xml:
            raise ValueError("Frozen Pandora XML disagrees with its parameter catalog")
    relative = archive.relative_to(Path(cdir).resolve()).as_posix()
    manifest.update(archive=relative, driver_files=code)
    tunesetup.aux_param_space["runtime_files"] = {relative: str(archive)}
    tunesetup.aux_param_space["runtime_dependencies"] = manifest["dependencies"]
    for card in cards:
        card["runtime"] = manifest
        for key in ("settings_template", "reco_steering", "run_dir"):
            card[key] = Path(card[key]).relative_to(root).as_posix()
        if card.get("baseline_config"):
            summary = icetune_push.load_summary(card["baseline_config"], baseline=True)
            card["baseline"] = summary["config"]
            if any(key.startswith(("PXML_", "PBOOL_")) for key in card["baseline"]):
                source = root / card["settings_template"]
                list(_xml_targets(ET.parse(source).getroot(), _push_catalog(summary, source), card["baseline"]))
        card.pop("baseline_config", None)
        card.pop("input_roots")
        card.pop("pandora_dir")


# Require Pandora initialization to contain metadata without generator caches
def stage_bootstrap(*, cdir, bootstrap):
    if (bootstrap.get("archive_sha256") is not None or bootstrap.get("archive_url") is not None
            or bootstrap.get("archive_size") != 0 or bootstrap.get("files") != []):
        raise RuntimeError("Ray PANDORA INIT bootstrap must not contain staged files")

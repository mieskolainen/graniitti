# Scratch local icetune cache construction and shared storage transport
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import errno
import hashlib
import json
import logging
import os
import pathlib
import shutil
import subprocess
import tempfile
import uuid
from collections.abc import Callable, Iterable
from urllib.parse import urlsplit

from core.io.files import ensure_dir
from core.io.serialize import load_json_file, sha256_file, write_json_file
from core.tune import short_id

COMMAND_TIMEOUT_S = max(1.0, float(os.environ.get("ICETUNE_CACHE_COMMAND_TIMEOUT_S", "600")))


class PermanentConfigurationError(RuntimeError):
    """Error which must abort a run instead of consuming transient retries."""


class ContentArchiveError(RuntimeError):
    """Cached bytes do not match the bootstrap manifest."""


# Compute true for an XRootD URL
def is_xrootd_url(value: str | os.PathLike[str]) -> bool:
    return os.fspath(value).startswith("root://")


# Compute a stable SHA-256 digest for one JSON-compatible payload
def json_fingerprint(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


# Compute local regular files below one path
def files_below(root: pathlib.Path, relative: str) -> list[pathlib.Path]:
    path = root / relative
    if path.is_file():
        return [path]
    return sorted(item for item in path.rglob("*") if item.is_file()) if path.is_dir() else []


# Compute checksummed records for a collection of local files
def file_records(
    paths: Iterable[str | os.PathLike[str]],
    *,
    root: pathlib.Path | None = None,
    allowed: Callable[[pathlib.PurePosixPath], bool] | None = None,
) -> list[dict]:
    from core import resource

    records = []
    for path in sorted({pathlib.Path(value).resolve() for value in paths}):
        try:
            relative = path.relative_to(root).as_posix() if root is not None else str(path)
        except ValueError:
            if allowed is not None:
                raise PermanentConfigurationError(f"File lies outside cache root: {path}") from None
            relative = str(path)
        if path.is_relative_to(resource("")):
            relative = "package:" + path.relative_to(resource("")).as_posix()
        if allowed is not None and not allowed(pathlib.PurePosixPath(relative)):
            raise PermanentConfigurationError(f"Refusing to cache unsupported output: {path}")
        records.append({"path": relative, "sha256": sha256_file(path), "size": path.stat().st_size})
    return sorted(records, key=lambda record: record["path"])


# Compute one validated relative path from a content archive manifest
def safe_manifest_path(value: object, *, allowed: Callable[[pathlib.PurePosixPath], bool]) -> pathlib.PurePosixPath:
    relative = pathlib.PurePosixPath(str(value))
    if relative.is_absolute() or not relative.parts or ".." in relative.parts or not allowed(relative):
        raise PermanentConfigurationError(f"Unsafe bootstrap manifest path: {relative}")
    return relative


# Compute the XRootD server and absolute namespace path for one URL
def _split_xrootd_url(url: str) -> tuple[str, str]:
    parsed = urlsplit(url)
    if parsed.scheme != "root" or not parsed.netloc:
        raise PermanentConfigurationError(f"Invalid XRootD URL: {url}")
    return f"root://{parsed.netloc}", "/" + parsed.path.lstrip("/")


# Run one external command and raise with its captured output on failure
def _run_command(cmd: list[str], *, context: str, env: dict | None = None) -> str:
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
            check=False,
            timeout=COMMAND_TIMEOUT_S,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"{context} timed out after {COMMAND_TIMEOUT_S:.1f} s") from exc
    if result.returncode != 0:
        output = (result.stdout or "").strip()
        raise RuntimeError(f"{context} ({result.returncode}): {' '.join(cmd)}\n{output}")
    return result.stdout or ""


# Run one external transport command and raise with captured output on failure
def _run_transport(cmd: list[str]) -> str:
    return _run_command(cmd, context="Transport command failed")


# Compute whether one local or XRootD object exists
def location_exists(location: str) -> bool:
    if not is_xrootd_url(location):
        return pathlib.Path(location.removeprefix("file://")).exists()
    server, path = _split_xrootd_url(location)
    try:
        result = subprocess.run(
            ["xrdfs", server, "stat", path],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=COMMAND_TIMEOUT_S,
        )
    except subprocess.TimeoutExpired as exc:
        raise OSError(errno.ETIMEDOUT, f"Timed out checking output location {location}") from exc
    if result.returncode == 0:
        return True
    output = (result.stdout or "").strip()
    lowered = output.lower()
    if "no such file" in lowered or "not found" in lowered:
        return False
    raise OSError(errno.EIO, f"Unable to stat output location {location}: {output}")


# Copy one local or XRootD object to a local file
def copy_to_local(source: str, destination: str | os.PathLike[str]) -> None:
    target = pathlib.Path(destination)
    ensure_dir(target.parent)
    if is_xrootd_url(source):
        _run_transport(["xrdcp", "--force", source, str(target)])
        return
    local = pathlib.Path(source.removeprefix("file://"))
    if target.exists() and os.path.samefile(local, target):
        return
    shutil.copy2(local, target)


# Publish one immutable local or remote object
def _publish_immutable(source: pathlib.Path, destination: str) -> bool:
    if not is_xrootd_url(destination):
        target = pathlib.Path(destination.removeprefix("file://"))
        ensure_dir(target.parent)
        if target.exists():
            return False
        temporary = target.with_name(f".{target.name}.tmp.{uuid.uuid4().hex}")
        try:
            # Commit an independent copy before exposing it on shared storage
            with source.open("rb") as incoming, temporary.open("xb") as outgoing:
                shutil.copyfileobj(incoming, outgoing)
                outgoing.flush()
                os.fsync(outgoing.fileno())
            os.link(temporary, target)
        except FileExistsError:
            return False
        finally:
            temporary.unlink(missing_ok=True)
        return True

    server, final_path = _split_xrootd_url(destination)
    parent = str(pathlib.PurePosixPath(final_path).parent)
    _run_transport(["xrdfs", server, "mkdir", "-p", parent])
    result = subprocess.run(
        ["xrdcp", "--posc", str(source), destination],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=COMMAND_TIMEOUT_S,
    )
    if result.returncode == 0:
        return True
    if not location_exists(destination):
        raise RuntimeError(f"Unable to publish immutable cache object {destination}: {(result.stdout or '').strip()}")
    return False


# Publish one local file without replacing an existing immutable object
def publish_file_immutable(source: str | os.PathLike[str], destination: str) -> bool:
    return _publish_immutable(pathlib.Path(source), destination)


# Preserve one unusable cache file before rebuilding its content
def _retire_cache_file(location: str) -> None:
    previous = f"{location}.{uuid.uuid4().hex}._old"
    if is_xrootd_url(location):
        server, path = _split_xrootd_url(location)
        _, old_path = _split_xrootd_url(previous)
        _run_transport(["xrdfs", server, "mv", path, old_path])
    else:
        try:
            pathlib.Path(location.removeprefix("file://")).rename(previous.removeprefix("file://"))
        except FileNotFoundError:
            return
    logging.getLogger(__name__).warning("Preserved unusable bootstrap cache file at %s", previous)


# Run one archive command and return its captured output
def _run_archive(command: list[str], *, label: str, env: dict | None = None) -> str:
    return _run_command(command, context=f"{label} bootstrap command failed", env=env)


# Compute one validated SHA-256 digest
def _cache_digest(value: object, *, label: str) -> str:
    digest = str(value)
    if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
        raise PermanentConfigurationError(f"Invalid {label} SHA-256 digest: {digest}")
    return digest


# Compute the immutable input manifest location for one bootstrap fingerprint
def _content_manifest_location(cache_base_url: str, kind: str, fingerprint: str) -> str:
    digest = _cache_digest(fingerprint, label="bootstrap fingerprint")
    return f"{cache_base_url.rstrip('/')}/{kind}/{short_id(digest)}.json"


# Compute the immutable archive location for one archive checksum
def _content_archive_location(cache_base_url: str, kind: str, archive_sha256: str) -> str:
    digest = _cache_digest(archive_sha256, label="bootstrap archive")
    return f"{cache_base_url.rstrip('/')}/{kind}/objects/{short_id(digest)}.tar.zst"


# Download and validate one content-addressed manifest
def _load_content_manifest(
    location: str, temporary_dir: pathlib.Path, *, schema_version: int, label: str
) -> dict | None:
    if not location_exists(location):
        return None
    local = temporary_dir / "bootstrap-manifest.json"
    copy_to_local(location, local)
    manifest = load_json_file(local)
    if manifest.get("schema_version") != schema_version:
        raise PermanentConfigurationError(f"Unsupported {label} bootstrap manifest at {location}")
    return manifest


# Download and checksum one content-addressed archive
def _verified_content_archive(archive_url: str, manifest: dict, temporary_dir: pathlib.Path) -> pathlib.Path:
    archive = temporary_dir / "bootstrap.tar.zst"
    copy_to_local(archive_url, archive)
    actual_size = archive.stat().st_size
    expected_size = int(manifest["archive_size"])
    if actual_size != expected_size:
        raise ContentArchiveError(
            f"Bootstrap archive size mismatch at {archive_url}: expected {expected_size}, got {actual_size}"
        )
    actual_sha256 = sha256_file(archive)
    if actual_sha256 != manifest["archive_sha256"]:
        raise ContentArchiveError(
            f"Bootstrap archive checksum mismatch at {archive_url}: "
            f"expected {manifest['archive_sha256']}, got {actual_sha256}"
        )
    return archive


# Compute an existing verified content-addressed cache entry
def reuse_content_archive(
    *, cache_base_url: str, kind: str, fingerprint: str, temporary_dir: pathlib.Path, schema_version: int, label: str
) -> dict | None:
    manifest_url = _content_manifest_location(cache_base_url, kind, fingerprint)
    manifest = _load_content_manifest(manifest_url, temporary_dir, schema_version=schema_version, label=label)
    if manifest is None:
        return None
    archive_url = _content_archive_location(
        cache_base_url, kind, _cache_digest(manifest.get("archive_sha256"), label="bootstrap archive")
    )
    if (
        manifest.get("fingerprint") != fingerprint
        or manifest.get("kind") != kind
        or manifest.get("archive_url") != archive_url
    ):
        raise PermanentConfigurationError(f"Bootstrap manifest identity mismatch at {manifest_url}")
    if not location_exists(archive_url):
        _retire_cache_file(manifest_url)
        return None
    try:
        _verified_content_archive(archive_url, manifest, temporary_dir)
    except ContentArchiveError as exc:
        logging.getLogger(__name__).warning("Rebuilding bootstrap cache: %s", exc)
        _retire_cache_file(manifest_url)
        _retire_cache_file(archive_url)
        return None
    return manifest


# Create one deterministic zstd-compressed tar archive from manifest records
def _create_content_archive(*, root: pathlib.Path, records: list[dict], output: pathlib.Path, label: str) -> None:
    relative_paths = [record["path"] for record in records]
    if not relative_paths:
        raise PermanentConfigurationError(f"{label} initialization produced no reusable files")
    ensure_dir(output.parent)
    environment = os.environ.copy()
    environment["ZSTD_CLEVEL"] = "1"
    _run_archive(
        [
            "tar",
            "--zstd",
            "--sort=name",
            "--mtime=@0",
            "--owner=0",
            "--group=0",
            "--numeric-owner",
            "--mode=u+rwX,go+rX,go-w",
            "-cf",
            str(output),
            "-C",
            str(root),
            "--",
            *relative_paths,
        ],
        label=label,
        env=environment,
    )


# Publish and verify one immutable content-addressed archive and manifest
def publish_content_archive(
    *,
    cache_base_url: str,
    kind: str,
    fingerprint: str,
    root: pathlib.Path,
    records: list[dict],
    temporary_dir: pathlib.Path,
    schema_version: int,
    label: str,
    manifest_extra: dict | None = None,
) -> dict:
    existing = reuse_content_archive(
        cache_base_url=cache_base_url,
        kind=kind,
        fingerprint=fingerprint,
        temporary_dir=temporary_dir,
        schema_version=schema_version,
        label=label,
    )
    if existing is not None:
        return existing
    archive = temporary_dir / "bootstrap.tar.zst"
    _create_content_archive(root=root, records=records, output=archive, label=label)
    archive_sha256 = sha256_file(archive)
    archive_url = _content_archive_location(cache_base_url, kind, archive_sha256)
    manifest_url = _content_manifest_location(cache_base_url, kind, fingerprint)
    manifest = {
        "archive_sha256": archive_sha256,
        "archive_size": archive.stat().st_size,
        "archive_url": archive_url,
        "files": records,
        "fingerprint": fingerprint,
        "kind": kind,
        "schema_version": schema_version,
        **(manifest_extra or {}),
    }
    manifest_file = temporary_dir / "manifest.json"
    write_json_file(manifest_file, manifest, indent=2, sort_keys=True, newline=True)
    publish_file_immutable(archive, archive_url)
    verification_dir = temporary_dir / "published"
    ensure_dir(verification_dir)
    try:
        _verified_content_archive(archive_url, manifest, verification_dir)
    except ContentArchiveError as exc:
        logging.getLogger(__name__).warning("Republishing bootstrap cache: %s", exc)
        _retire_cache_file(archive_url)
        publish_file_immutable(archive, archive_url)
        _verified_content_archive(archive_url, manifest, verification_dir)
    publish_file_immutable(manifest_file, manifest_url)
    verified = reuse_content_archive(
        cache_base_url=cache_base_url,
        kind=kind,
        fingerprint=fingerprint,
        temporary_dir=temporary_dir,
        schema_version=schema_version,
        label=label,
    )
    if verified is None:
        raise RuntimeError(f"Committed bootstrap cache is incomplete at {manifest_url}")
    return verified


# Verify every staged file against the immutable manifest
def _verify_staged_records(stage_root: pathlib.Path, records: list[dict], allowed: Callable) -> None:
    for record in records:
        relative = safe_manifest_path(record["path"], allowed=allowed)
        path = stage_root.joinpath(*relative.parts)
        if not path.is_file() or path.stat().st_size != int(record["size"]):
            raise RuntimeError(f"Missing or truncated staged bootstrap file: {path}")
        if sha256_file(path) != record["sha256"]:
            raise RuntimeError(f"Checksum mismatch for staged bootstrap file: {path}")


# Install one verified local cache path
def _install_staged_path(staged: pathlib.Path, target: pathlib.Path) -> None:
    ensure_dir(target.parent)
    if (target.is_file() and not target.is_symlink() and staged.is_file()
            and target.stat().st_size == staged.stat().st_size and sha256_file(target) == sha256_file(staged)):
        return
    if os.path.lexists(target):
        target.rename(target.with_name(f"{target.name}.{uuid.uuid4().hex}._old"))
    os.rename(staged, target)


# Install complete cache directories and then isolated manifest files
def _install_staged_cache(
    *, root: pathlib.Path, stage_root: pathlib.Path, records: list[dict], cache_dirs: tuple[str, ...]
) -> None:
    cache_roots = tuple(pathlib.PurePosixPath(value) for value in cache_dirs)
    for relative in cache_roots:
        staged = stage_root.joinpath(*relative.parts)
        ensure_dir(staged)
        ignore = root.joinpath(*relative.parts, ".gitignore")
        if ignore.is_file() and not (staged / ".gitignore").exists():
            shutil.copy2(ignore, staged / ".gitignore")
        _install_staged_path(staged, root.joinpath(*relative.parts))
    for record in records:
        relative = pathlib.PurePosixPath(record["path"])
        if any(relative == cache_root or cache_root in relative.parents for cache_root in cache_roots):
            continue
        _install_staged_path(stage_root.joinpath(*relative.parts), root.joinpath(*relative.parts))


# Extract, verify, and install one immutable content-addressed archive
def stage_content_archive(
    *,
    root: pathlib.Path,
    bootstrap: dict,
    kind: str,
    cache_dirs: Iterable[str],
    allowed: Callable[[pathlib.PurePosixPath], bool],
    label: str,
) -> None:
    if bootstrap.get("kind") != kind:
        raise PermanentConfigurationError(f"{label} driver received a non-{label} bootstrap")
    root = root.resolve()
    cache_dirs = tuple(str(directory) for directory in cache_dirs)
    ensure_dir(root / "tmp")
    records = bootstrap.get("files", [])
    if not isinstance(records, list):
        raise PermanentConfigurationError("Bootstrap manifest files field is not a list")
    expected_paths = {
        safe_manifest_path(record.get("path"), allowed=allowed).as_posix()
        for record in records
        if isinstance(record, dict)
    }
    if len(expected_paths) != len(records):
        raise PermanentConfigurationError("Bootstrap manifest has malformed or duplicate file records")
    with tempfile.TemporaryDirectory(prefix="icetune-stage-", dir=root / "tmp") as temporary:
        stage_root = pathlib.Path(temporary) / "tree"
        ensure_dir(stage_root, exist_ok=False)
        archive = _verified_content_archive(str(bootstrap["archive_url"]), bootstrap, pathlib.Path(temporary))
        archive_paths = {
            safe_manifest_path(line.strip().rstrip("/"), allowed=allowed).as_posix()
            for line in _run_archive(["tar", "--zstd", "-tf", str(archive)], label=label).splitlines()
            if line.strip().rstrip("/")
        }
        if archive_paths != expected_paths:
            raise PermanentConfigurationError("Bootstrap archive contents do not match its file manifest")
        # Reuse verified immutable outputs already present in a shared checkout
        cache_roots = tuple(pathlib.PurePosixPath(value) for value in cache_dirs)
        selected = []
        for record in records:
            relative = safe_manifest_path(record["path"], allowed=allowed)
            target = root.joinpath(*relative.parts)
            in_cache = any(relative == prefix or prefix in relative.parents for prefix in cache_roots)
            if (in_cache or not target.is_file() or target.is_symlink() or target.stat().st_size != record["size"]
                    or sha256_file(target) != record["sha256"]):
                selected.append(record)
        selection = pathlib.Path(temporary) / "selected.bin"
        selection.write_bytes(b"".join(record["path"].encode("utf-8") + b"\0" for record in selected))
        # Match literal process names even when TAR_OPTIONS enables wildcards
        _run_archive(["tar", "--zstd", "-xf", str(archive), "-C", str(stage_root),
                      "--no-wildcards", "--null", "-T", str(selection)], label=label)
        _verify_staged_records(stage_root, selected, allowed)
        _install_staged_cache(root=root, stage_root=stage_root, records=selected, cache_dirs=cache_dirs)

# Shared file identities, directory creation and filesystem retries
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import datetime
import errno
import hashlib
import os
import random
import time
from pathlib import Path

IO_RETRIES = int(os.environ.get("ICETUNE_IO_RETRIES", "5"))
DIR_RETRIES = max(1, int(os.environ.get("ICETUNE_IO_DIR_RETRIES", "60")))


# Create a UTC dated output directory without replacing an existing run
def dated_directory(parent, prefix=""):
    stamp = prefix + datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%d_%H-%M-%S_UTC")
    index = 1
    while True:
        path = Path(parent) / (stamp if index == 1 else f"{stamp}_{index}")
        try:
            path.mkdir(parents=True)
            return path
        except FileExistsError:
            index += 1


# Compute bounded exponential retry delay with jitter
def retry_delay(attempt: int, *, base: float = 0.02, cap: float = 1.0) -> float:
    return min(cap, base * 2 ** max(0, int(attempt) - 1)) * (0.5 + random.random())


# Compute whether a storage I/O error is normally transient
def is_retryable_io_error(exc: BaseException) -> bool:
    retryable = {
        errno.EAGAIN,
        errno.EBUSY,
        errno.EIO,
        errno.ENOENT,
        errno.EPERM,
        errno.ETIMEDOUT,
        *(
            getattr(errno, name)
            for name in ("ENETDOWN", "ENETRESET", "ENETUNREACH", "ENOTCONN", "ESTALE")
            if hasattr(errno, name)
        ),
    }
    return getattr(exc, "errno", None) in retryable


# Ensure a directory becomes visible despite shared filesystem races
def ensure_dir(path: str | os.PathLike, *, exist_ok: bool = True, mode: int = 0o777) -> None:
    """Create missing parents with shared filesystem retries"""
    if not path:
        return
    for attempt in range(1, DIR_RETRIES + 1):
        try:
            os.makedirs(path, mode=mode, exist_ok=exist_ok)
            if not exist_ok or os.path.isdir(path):
                return
        except FileExistsError:
            if not exist_ok:
                raise
            if os.path.isdir(path):
                return
            raise NotADirectoryError(errno.ENOTDIR, "Directory path is occupied by a file", path) from None
        except OSError as exc:
            if not is_retryable_io_error(exc):
                raise
        if attempt == DIR_RETRIES:
            raise FileNotFoundError(errno.ENOENT, "Directory was not visible after creation", path)
        time.sleep(retry_delay(attempt, base=0.05, cap=2.0))


# Run one filesystem operation with bounded transient retries
def retry_io(action, *, retries: int = IO_RETRIES, missing=None, base: float = 0.05, cap: float = 2.0):
    for attempt in range(1, max(1, int(retries)) + 1):
        try:
            return action()
        except FileNotFoundError:
            if missing is not None:
                return missing
            if attempt >= retries:
                raise
        except OSError as exc:
            if attempt >= retries or not is_retryable_io_error(exc):
                raise
        time.sleep(retry_delay(attempt, base=base, cap=cap))
    raise RuntimeError("Filesystem retry loop exhausted")


# Hash a stable input and allow CVMFS cold reads time to recover
def file_identity(path):
    path = Path(path)

    # Restart the complete read after transient storage failures
    def read():
        before = path.stat()
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError(f"Input changed while hashing: {path}")
        return dict(sha256=digest, size=after.st_size)

    try:
        if path.is_relative_to("/cvmfs"):
            return retry_io(read, retries=max(IO_RETRIES, 9), base=1.0, cap=15.0)
        return retry_io(read)
    except OSError as exc:
        if exc.filename is None:
            exc.filename = str(path)
        raise


# Read and validate external runtime inputs before admitting a worker
def check_dependencies(dependencies):
    for path, expected in dependencies.items():
        if file_identity(path) != expected:
            raise ValueError(f"Runtime dependency changed: {path}")

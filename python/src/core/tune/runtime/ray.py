# Explicit Ray runtime configuration for restricted multi-user systems
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import copy
import os
import pathlib


# Configure Ray before importing the tuning backend
def configure_ray_runtime() -> None:
    os.environ.setdefault("RAY_ENABLE_UV_RUN_RUNTIME_ENV", "0")

    try:
        import ray._private.ray_constants as ray_constants

        ray_constants.RAY_ENABLE_UV_RUN_RUNTIME_ENV = False
    except Exception:
        pass

    try:
        import ray._private.node as ray_node

        original = ray_node.Node._get_system_processes_for_resource_isolation

        # Fall back to known Ray processes when procfs access is restricted
        def safe_system_processes(node):
            try:
                return original(node)
            except Exception as exc:
                if exc.__class__.__name__ != "AccessDenied":
                    raise
                pids = [str(process[0].process.pid) for process in node.all_processes.values()]
                return ",".join(pids)

        ray_node.Node._get_system_processes_for_resource_isolation = safe_system_processes
    except Exception:
        pass


RAY_RUNTIME_EXCLUDES = (".git", "build", "figs", "output", "runs", "tmp", "**/__pycache__", "**/*.pyc", "**/*._old")


# Compute the uploaded worker local simulator runtime root
def ray_worker_root() -> pathlib.Path:
    return pathlib.Path.cwd().resolve()


# Include declared driver runtime files while keeping other scratch excluded
def runtime_excludes(files):
    excludes = list(RAY_RUNTIME_EXCLUDES)
    paths = [pathlib.PurePosixPath(path) for path in files]
    if any(path.is_absolute() or ".." in path.parts for path in paths):
        raise ValueError("Driver runtime files must be relative to the checkout")
    parents = {parent for path in paths for parent in path.parents
               if str(parent) != "." and parent.parts[0] in RAY_RUNTIME_EXCLUDES}
    for parent in sorted(parents, key=lambda path: (len(path.parts), str(path))):
        excludes.extend([f"!/{parent}", f"/{parent}/*"])
    return excludes + [f"!/{path}" for path in paths]


# Rebase checkout paths from the Ray head runtime into the uploaded worker runtime
def rebase_ray_paths(value, *, source_root: pathlib.Path, worker_root: pathlib.Path, shared_paths=frozenset()):
    if isinstance(value, str):
        if value in shared_paths:
            return value
        path = pathlib.Path(value).expanduser()
        if not path.is_absolute():
            return value
        try:
            relative = path.relative_to(source_root)
        except ValueError:
            return value
        return str(worker_root / relative)
    if isinstance(value, dict):
        return {
            key: rebase_ray_paths(item, source_root=source_root, worker_root=worker_root, shared_paths=shared_paths) for key, item in value.items()
        }
    if isinstance(value, list):
        return [rebase_ray_paths(item, source_root=source_root, worker_root=worker_root, shared_paths=shared_paths) for item in value]
    if isinstance(value, tuple):
        return tuple(rebase_ray_paths(item, source_root=source_root, worker_root=worker_root, shared_paths=shared_paths) for item in value)
    return value


# Prepare simulator paths and libraries in the assigned Ray worker
def worker_runtime(param: dict, simdriver) -> dict:
    param = copy.deepcopy(param)
    if param.get("ray_upload_runtime"):
        source_root = pathlib.Path(param["cdir"]).expanduser()
        worker_root = ray_worker_root()
        events_dir = param.get("events_dir")
        param = rebase_ray_paths(param, source_root=source_root, worker_root=worker_root,
                                 shared_paths=frozenset(simdriver.shared_runtime_files(param)))
        # Event samples go to shared storage even when the generator checkout is uploaded
        if events_dir:
            param["events_dir"] = events_dir
        for name in ("dataset_paths", "cuts"):
            if hasattr(simdriver, name):
                setattr(
                    simdriver,
                    name,
                    rebase_ray_paths(getattr(simdriver, name), source_root=source_root, worker_root=worker_root),
                )
        param["cdir"] = str(worker_root)
    if {"cdir", "libdir", "PYTHON_VERSION"}.issubset(param):
        simdriver.runtime_environment(
            cdir=param["cdir"], libdir=param["libdir"], python_version=param["PYTHON_VERSION"], apply=True
        )
    return param

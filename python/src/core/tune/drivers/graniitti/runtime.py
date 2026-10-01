# GRANIITTI runtime path and version helpers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>

from __future__ import annotations

import os
import pathlib
import sys

from core.io.serialize import load_json_file


# Accept only generator caches and declared icetune initialization outputs
def allowed_bootstrap_path(relative: pathlib.PurePosixPath) -> bool:
    if relative.parts and relative.parts[0] in {"eikonal", "sudakov", "vgrid"}:
        return True
    if len(relative.parts) < 5 or relative.parts[:2] != ("runs", "icetune") or relative.parts[3] != "results":
        return False
    return (len(relative.parts) > 5 and relative.parts[4] == "amplitude") or (
        len(relative.parts) == 5 and relative.name in {"data_covariance.npz", "data_covariance.json"}
    )


# Compute the external GRANIITTI I/O installation path
def default_library_path() -> str:
    return os.environ.get("GRANIITTI_IO_PATH") or os.path.join(os.environ["HOME"], "local")


# Build GRANIITTI runtime variables and optionally apply them
def environment(*, cdir: str, libdir: str | None, python_version: str, apply: bool = False) -> dict[str, str]:
    root = libdir or default_library_path()
    old_ld = os.environ.get("LD_LIBRARY_PATH", "")
    old_python = os.environ.get("PYTHONPATH", "")
    variables = {
        "LD_LIBRARY_PATH": (f"{root}/HEPMC3/lib:{root}/HEPMC3/lib64:{root}/LHAPDF/lib:{root}/LHAPDF/lib64:{old_ld}"),
        "PYTHONPATH": (f"{cdir}:{root}/HEPMC3/lib/python{python_version}/site-packages:{old_python}"),
    }
    if apply:
        os.environ.update(variables)
    for path in variables["PYTHONPATH"].split(os.pathsep):
        if path and path not in sys.path:
            sys.path.append(path)
    return variables


# Compute the GRANIITTI version number
def version(cdir: str) -> float:
    return round(float(load_json_file(os.path.join(cdir, "VERSION.json"))["version"]), 3)


# Restore verified generator caches and initialization outputs
def stage_bootstrap(*, cdir, bootstrap):
    from core.tune.cache import stage_content_archive

    stage_content_archive(root=pathlib.Path(cdir).resolve(), bootstrap=bootstrap, kind="graniitti",
                          cache_dirs=("eikonal", "sudakov", "vgrid"), allowed=allowed_bootstrap_path, label="GRANIITTI")

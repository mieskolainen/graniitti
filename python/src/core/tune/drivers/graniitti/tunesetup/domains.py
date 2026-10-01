# Shared GRANIITTI tunesetup optimizer-domain helpers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from functools import cache
from pathlib import Path

import numpy as np
import pyjson5 as json5
from ray import tune

from core.io.serialize import load_json_file
from core.tune.parameters import tools

ANGULAR_COUPLING_FIELDS = {"alpha_ls", "g_ls", "helicity"}
TP_RES_LABELS = {
    (0, 1, 1): "0,0 2,2",
    (0, -1, 1): "1,1 3,3",
    (2, -1, -1): "Gamma0 Gamma2",
    (2, 1, 1): "2,2 4,4",
    (4, 1, 1): "0,2 (2,0)-(2,2) (2,0)+(2,2) 2,4 4,2 4,4 6,4",
}
SOURCE = ContextVar("graniitti_tuning_context")


# Bind explicit source cards and fit settings for one thread or task
@contextmanager
def context(*, cdir, model_path, settings):
    token = SOURCE.set((Path(cdir).resolve(), Path(model_path).resolve(), settings))
    try:
        yield
    finally:
        SOURCE.reset(token)


# Resolve the simulator input directory for this tune construction
def source_root() -> Path:
    return SOURCE.get()[0]


# Locate the selected source model for tuning domains
def model_dir() -> Path:
    return SOURCE.get()[1]


# Read explicit sampling choices supplied by the tuning card
def settings() -> dict:
    return SOURCE.get()[2]


# Compute the coherent MP tuning sectors connected to M=0 in the chosen frame
def mp_spin_sectors(j: int, frame: str) -> list[int]:
    if frame == "CM":
        return [0]
    if frame in {"CS", "HX"}:
        return list(range(0, j + 1, 2))
    raise ValueError(f"Unknown MP spin frame {frame}")


# Load and cache one JSON5 model-data card
@cache
def load_model_json(path) -> dict:
    return load_json_file(path, loader=json5.load)


# Load the shared GENERAL source model card
def load_general_data() -> dict:
    return load_model_json(model_dir() / "GENERAL.json")


# Compute one canonical integer or half-integer angular quantum number
def _angular_number(value: object, *, half_integer: bool, field: str) -> str:
    number = float(value)
    scale = 2 if half_integer else 1
    discrete = int(round(scale * number))
    if not np.isfinite(number) or not np.isclose(scale * number, discrete, rtol=0.0, atol=1.0e-9):
        unit = "integer or half-integer" if half_integer else "integer"
        raise ValueError(f"{field} row label must be {unit}, got {value}")
    if not half_integer or discrete % 2 == 0:
        return str(discrete if not half_integer else discrete // 2)
    return f"{discrete}/2"


# Compute one physical angular coupling row name from its discrete columns
def angular_row_name(field: str, row: list | tuple) -> str:
    field = str(field)
    if field not in ANGULAR_COUPLING_FIELDS:
        raise ValueError(f'Unknown angular coupling field "{field}"')
    if not isinstance(row, (list, tuple)) or len(row) < 2:
        raise ValueError(f"{field} row must contain two angular labels")
    half_integer = field == "helicity"
    first = _angular_number(row[0], half_integer=half_integer, field=field)
    second = _angular_number(row[1], half_integer=half_integer, field=field)
    if len(row) == 5 and field in {"g_ls", "helicity"}:
        m = _angular_number(row[2], half_integer=False, field=field)
        return f"{field}({first},{second},{m})"
    return f"{field}({first},{second})"


# Compute physical Tensor Pomeron resonance coefficient names from JPC
def tp_res_names(param: dict) -> tuple[str, ...]:
    key = tuple(int(param[field]) for field in ("spinX2", "P", "C"))
    if key not in TP_RES_LABELS:
        raise ValueError(f"Unsupported TP resonance quantum numbers {key}")
    return tuple(f"g_tensor({label})" for label in TP_RES_LABELS[key].split())


# Compute covariant Tensor Pomeron continuum coefficient names
def tp_con_names(couplings: list) -> tuple[str, ...]:
    if not isinstance(couplings, list) or len(couplings) not in {1, 2}:
        raise ValueError("TP continuum g_tensor must contain one or two coefficients")
    return ("g_tensor",) if len(couplings) == 1 else ("g_tensor(Gamma0)", "g_tensor(Gamma2)")


# Add signed direction angles with a fixed first row and the selected angular range
def add_direction_params(parameters: dict, keys: list[str], *, geometry: str = "projective") -> None:
    for index, key in enumerate(keys):
        limit = np.pi if geometry == "spherical" and index + 1 == len(keys) else tools.PROJECTIVE_ANGLE_LIMIT
        parameters[key] = tune.uniform(-limit, limit)


# Add an overall norm or Cartesian coupling followed by the signed row direction
def add_residue_params(
    parameters: dict, *, rows: list[str], bounds: tuple[float, float], mode: str, geometry: str = "projective"
) -> None:
    if mode == "mag_phase_cartesian":
        for key in tools.projective_coefficient_keys(rows[0]):
            parameters[key] = tune.uniform(-bounds[1], bounds[1])
    else:
        parameters[tools.projective_norm_key(rows[0])] = tune.uniform(*bounds)
    suffix = tools.SPHERICAL_VECTOR_SUFFIX if geometry == "spherical" else tools.PROJECTIVE_VECTOR_SUFFIX
    add_direction_params(parameters, [f"{row}{suffix}" for row in rows[1:]], geometry=geometry)


# Add one magnitude-phase or Cartesian coupling domain to a parameter space
def add_magnitude_phase_params(
    parameters: dict,
    *,
    base_key: str,
    magnitude_key: str,
    phase_key: str,
    mode: str,
    magnitude_bounds: tuple[float, float],
    cartesian_bounds: tuple[tuple[float, float], tuple[float, float]] | None = None,
    mode_name: str = "mode",
) -> None:
    mode = tools.validate_tune_mode(mode, name=mode_name)
    magnitude_min, magnitude_max = magnitude_bounds
    if not np.all(np.isfinite(magnitude_bounds)) or not 0.0 <= magnitude_min <= magnitude_max:
        raise ValueError("Coupling magnitude bounds must be finite, nonnegative and ordered")
    if mode == "none" or magnitude_max <= 0.0:
        return
    if tools.tune_mode_is_cartesian(mode):
        limits = cartesian_bounds or ((-magnitude_max, magnitude_max),) * 2
        if np.shape(limits) != (2, 2) or not np.all(np.isfinite(limits)):
            raise ValueError("Cartesian coupling bounds require two finite intervals")
        limits = tuple(tuple(sorted(float(value) for value in bounds)) for bounds in limits)
        if all(lower >= upper for lower, upper in limits):
            return
        if any(lower >= upper for lower, upper in limits):
            raise ValueError("Both Cartesian coupling components must have a nonzero tuning range")
        keys = (tools.complex_re_key(base_key), tools.complex_im_key(base_key))
        for key, bounds in zip(keys, limits, strict=True):
            parameters[key] = tune.uniform(*sorted(float(value) for value in bounds))
        return
    if tools.tune_mode_has_magnitude(mode) and magnitude_min < magnitude_max:
        parameters[magnitude_key] = tune.uniform(magnitude_min, magnitude_max)
    if mode.endswith("phase_raw"):
        parameters[tools.raw_phase_key(phase_key)] = tune.uniform(-np.pi, np.pi)
    elif mode.endswith("phase_cayley"):
        parameters[tools.phase_u_key(phase_key)] = tune.uniform(-1.0, 1.0)
        parameters[tools.phase_v_key(phase_key)] = tune.uniform(-1.0, 1.0)

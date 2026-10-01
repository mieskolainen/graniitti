# Validated ICEBO JSON configuration loading
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import copy
import math
import pathlib
import re

from core.io.serialize import load_json_file

ICEBO_SECTIONS = {
    "search": {"device", "direction", "dtype", "duplicate_tolerance"},
    "gp": {
        "kernel",
        "cholesky_attempts",
        "fit_restarts",
        "fit_steps",
        "fit_tolerance_grad",
        "fit_tolerance_change",
        "gradient_clip",
        "initial_jitter",
        "initial_lengthscale",
        "initial_noise",
        "initial_outputscale",
        "jitter_multiplier",
        "learning_rate",
        "lengthscale_prior_log_std",
        "lengthscale_prior_median",
        "noise_prior_log_std",
        "noise_prior_median",
        "outputscale_prior_log_std",
        "outputscale_prior_median",
        "posterior_variance_floor",
    },
    "acquisition": {
        "gradient_steps",
        "lbfgs_tolerance_grad",
        "lbfgs_tolerance_change",
        "mc_cholesky_jitter",
        "mc_samples",
        "optimization_batch_size",
        "raw_samples",
        "restarts",
        "prune_baseline",
        "tau_max",
        "tau_relu",
    },
    "turbo": {
        "enabled",
        "failure_base",
        "failure_dimension_scale",
        "failure_length_multiplier",
        "improvement_relative_tolerance",
        "length_initial",
        "length_max",
        "length_min",
        "restart_min_points",
        "restart_points_per_dimension",
        "success_length_multiplier",
        "success_tolerance",
    },
}


# Require a plain integer without accepting Boolean values
def _integer(value, *, name: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f'ICEBO setting "{name}" must be an integer >= {minimum}')
    return value


# Require one finite real number with an optional strict lower limit
def _real(value, *, name: str, lower: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'ICEBO setting "{name}" must be numeric')
    numeric = float(value)
    if not math.isfinite(numeric) or (lower is not None and numeric <= lower):
        raise ValueError(f'ICEBO setting "{name}" is outside its valid range')
    return numeric


# Validate section names and exact nested keys
def _validate_layout(payload: dict, *, source: pathlib.Path) -> None:
    if payload.get("schema_version") != 1:
        raise ValueError(f"Unsupported ICEBO configuration schema: {source}")
    sections = set(payload) - {"schema_version"}
    if sections != set(ICEBO_SECTIONS):
        raise ValueError(f"ICEBO configuration sections must be exactly {sorted(ICEBO_SECTIONS)}")
    for section, expected in ICEBO_SECTIONS.items():
        values = payload[section]
        if not isinstance(values, dict) or set(values) != expected:
            raise ValueError(f'ICEBO section "{section}" keys must be exactly {sorted(expected)}')


# Validate all nested ICEBO algorithm settings
def validate_icebo_config(payload: object, *, source: pathlib.Path) -> dict:
    if not isinstance(payload, dict):
        raise ValueError(f"ICEBO configuration must be a mapping: {source}")
    _validate_layout(payload, source=source)
    output = copy.deepcopy(payload)
    search = output["search"]
    if search["direction"] not in {"minimize", "maximize"}:
        raise ValueError("ICEBO search.direction must be minimize or maximize")
    if search["dtype"] not in {"auto", "float32", "float64"}:
        raise ValueError("ICEBO search.dtype must be auto, float32 or float64")
    if not isinstance(search["device"], str) or re.fullmatch(r"auto|cpu|cuda(?::[0-9]+)?", search["device"]) is None:
        raise ValueError("ICEBO search.device must be auto, cpu or cuda[:index]")
    _real(search["duplicate_tolerance"], name="search.duplicate_tolerance", lower=0.0)

    gp = output["gp"]
    if gp["kernel"] not in {"rbf", "matern52"}:
        raise ValueError("ICEBO gp.kernel must be rbf or matern52")
    for key in ("fit_steps", "fit_restarts", "cholesky_attempts"):
        _integer(gp[key], name=f"gp.{key}", minimum=1)
    for key in (
        "fit_tolerance_grad",
        "fit_tolerance_change",
        "gradient_clip",
        "initial_jitter",
        "initial_lengthscale",
        "initial_noise",
        "initial_outputscale",
        "learning_rate",
        "jitter_multiplier",
        "lengthscale_prior_log_std",
        "lengthscale_prior_median",
        "noise_prior_log_std",
        "noise_prior_median",
        "outputscale_prior_log_std",
        "outputscale_prior_median",
        "posterior_variance_floor",
    ):
        _real(gp[key], name=f"gp.{key}", lower=0.0)
    if gp["jitter_multiplier"] <= 1.0:
        raise ValueError("ICEBO gp.jitter_multiplier must exceed one")

    acquisition = output["acquisition"]
    for key, minimum in (
        ("gradient_steps", 0),
        ("mc_samples", 8),
        ("optimization_batch_size", 1),
        ("raw_samples", 1),
        ("restarts", 1),
    ):
        _integer(acquisition[key], name=f"acquisition.{key}", minimum=minimum)
    for key in (
        "lbfgs_tolerance_grad",
        "lbfgs_tolerance_change",
        "mc_cholesky_jitter",
        "tau_max",
        "tau_relu",
    ):
        _real(acquisition[key], name=f"acquisition.{key}", lower=0.0)
    if not isinstance(acquisition["prune_baseline"], bool):
        raise ValueError("ICEBO acquisition.prune_baseline must be Boolean")

    turbo = output["turbo"]
    if not isinstance(turbo["enabled"], bool):
        raise ValueError("ICEBO turbo.enabled must be Boolean")
    for key in ("restart_min_points", "restart_points_per_dimension", "success_tolerance"):
        _integer(turbo[key], name=f"turbo.{key}", minimum=1)
    for key in (
        "failure_base",
        "failure_dimension_scale",
        "failure_length_multiplier",
        "improvement_relative_tolerance",
        "length_initial",
        "length_max",
        "length_min",
        "success_length_multiplier",
    ):
        _real(turbo[key], name=f"turbo.{key}", lower=0.0)
    if not turbo["length_min"] < turbo["length_initial"] <= turbo["length_max"]:
        raise ValueError("ICEBO TuRBO lengths require min < initial <= max")
    if turbo["failure_length_multiplier"] >= 1.0:
        raise ValueError("ICEBO TuRBO failures must contract the trust region")
    if turbo["success_length_multiplier"] <= 1.0:
        raise ValueError("ICEBO TuRBO successes must expand the trust region")
    return output


# Load and validate one ICEBO JSON configuration file
def load_icebo_config(path: str | pathlib.Path) -> dict:
    source = pathlib.Path(path).resolve()
    payload = load_json_file(source)
    return validate_icebo_config(payload, source=source)


# Compute a validated copy with explicit dotted-path overrides for tests and studies
def override_icebo_config(settings: dict, overrides: dict[str, object]) -> dict:
    output = copy.deepcopy(settings)
    for path, value in overrides.items():
        parts = path.split(".")
        target = output
        for part in parts[:-1]:
            if not isinstance(target, dict) or part not in target:
                raise ValueError(f'Unknown ICEBO override path "{path}"')
            target = target[part]
        leaf = parts[-1]
        if not isinstance(target, dict) or leaf not in target:
            raise ValueError(f'Unknown ICEBO override path "{path}"')
        target[leaf] = value
    return validate_icebo_config(output, source=pathlib.Path("<ICEBO overrides>"))

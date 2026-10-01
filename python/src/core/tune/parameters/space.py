# Shared typed optimizer parameter-space utilities
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import math
from dataclasses import dataclass
from numbers import Real

import numpy as np


@dataclass(frozen=True)
class ParameterSpec:
    """One finite uniform optimizer parameter"""

    name: str
    kind: str
    lower: float | int
    upper: float | int

    # Build one typed specification from a dictionary or Ray domain
    @classmethod
    def from_domain(cls, name: str, domain) -> ParameterSpec:
        kind = domain.get("type") if isinstance(domain, dict) else None
        if kind is not None and kind not in {"float", "int", "integer", "num", "uniform"}:
            raise ValueError(f'Parameter "{name}" is not a uniform numeric domain')
        lower = domain.get("lower") if isinstance(domain, dict) else getattr(domain, "lower", None)
        upper = domain.get("upper") if isinstance(domain, dict) else getattr(domain, "upper", None)
        if lower is None or upper is None:
            raise ValueError(f'Parameter "{name}" is not a uniform numeric domain')
        is_integer = kind in {"int", "integer"} or domain.__class__.__name__ == "Integer"
        if is_integer:
            lower_i = int(math.ceil(float(lower)))
            upper_i = int(math.floor(float(upper))) if kind else int(upper) - 1
            if upper_i < lower_i:
                raise ValueError(f'Parameter "{name}" has an empty integer domain')
            return cls(str(name), "int", lower_i, upper_i)
        lower_f, upper_f = float(lower), float(upper)
        if not np.isfinite(lower_f) or not np.isfinite(upper_f) or lower_f >= upper_f:
            raise ValueError(f'Parameter "{name}" requires finite lower < upper bounds')
        return cls(str(name), "uniform", lower_f, upper_f)

    # Compute the parameter dictionary used by optimizers
    def as_dict(self) -> dict:
        return {"type": self.kind, "lower": self.lower, "upper": self.upper}

    # Coerce and validate one optimizer configuration value
    def coerce(self, value):
        numeric = float(value)
        if not math.isfinite(numeric):
            raise ValueError(f'Parameter "{self.name}" received non-finite value {value!r}')
        if self.kind != "int":
            if numeric < self.lower or numeric > self.upper:
                raise ValueError(
                    f'Parameter "{self.name}" received out-of-bounds value {value!r}; '
                    f"expected [{self.lower}, {self.upper}]"
                )
            return numeric
        integer = int(round(numeric))
        if not math.isclose(numeric, integer, rel_tol=0.0, abs_tol=1.0e-9):
            raise ValueError(f'Integer parameter "{self.name}" received non-integer value {value!r}')
        if integer < self.lower or integer > self.upper:
            raise ValueError(
                f'Integer parameter "{self.name}" received out-of-bounds value {value!r}; '
                f"expected [{self.lower}, {self.upper}]"
            )
        return integer


# Normalize a mixed dictionary or Ray Tune parameter space
def normalize_param_space(param_space: dict) -> dict:
    return {str(name): ParameterSpec.from_domain(str(name), domain).as_dict() for name, domain in param_space.items()}


# Normalize a parameter space and require continuous parameters
def normalize_continuous_param_space(param_space: dict) -> dict:
    try:
        normalized = normalize_param_space(param_space)
    except ValueError as exc:
        raise ValueError("Gradient optimizers support only continuous uniform search spaces") from exc
    if not normalized:
        raise ValueError("Gradient optimizers require at least one tunable parameter")
    for name, spec in normalized.items():
        if spec["type"] != "uniform":
            raise ValueError(f'Gradient optimizers support only continuous uniform search spaces; key "{name}" is unsupported')
    return {name: {"lower": spec["lower"], "upper": spec["upper"]} for name, spec in sorted(normalized.items())}


# Compute whether a normalized parameter is integer-valued
def is_integer_bound(spec: dict) -> bool:
    return spec.get("type") == "int"


# Coerce a complete optimizer configuration
def typed_config(config: dict, bounds: dict) -> dict:
    return {
        name: ParameterSpec(name, spec["type"], spec["lower"], spec["upper"]).coerce(config[name])
        for name, spec in sorted(bounds.items())
    }


# Validate and type one complete initial optimizer configuration
def validate_initial_config(config: dict, param_space: dict) -> dict:
    if not isinstance(config, dict):
        raise TypeError("Initial fit configuration must be a parameter mapping")
    bounds = normalize_param_space(param_space)
    missing = sorted(set(bounds) - set(config))
    extra = sorted(set(config) - set(bounds))
    if missing or extra:
        raise ValueError(
            "Initial fit configuration keys differ from the parameter space: "
            f"missing={missing}, extra={extra}"
        )
    try:
        return typed_config(config, bounds)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Initial fit configuration is invalid: {exc}") from exc


# Draw one configuration using caller-supplied continuous and inclusive-integer samplers
def sample_config(bounds: dict, *, uniform, integer) -> dict:
    return {
        name: (
            integer(int(spec["lower"]), int(spec["upper"]))
            if is_integer_bound(spec)
            else uniform(spec["lower"], spec["upper"])
        )
        for name, spec in sorted(bounds.items())
    }


# Compute whether two scalar configuration values have the same optimizer meaning
def scalar_values_equal(left, right) -> bool:
    if hasattr(left, "item"):
        left = left.item()
    if hasattr(right, "item"):
        right = right.item()
    if isinstance(left, Real) and isinstance(right, Real):
        return math.isclose(float(left), float(right), rel_tol=1.0e-12, abs_tol=1.0e-12)
    return left == right


# Compute whether a candidate config contains one complete reference point
def config_contains(candidate: dict | None, reference: dict | None) -> bool:
    return bool(
        isinstance(candidate, dict)
        and isinstance(reference, dict)
        and all(key in candidate and scalar_values_equal(candidate[key], value) for key, value in reference.items())
    )


# Map physical coordinates into the unit cube
def to_unit(values: np.ndarray, bounds: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)
    lower = bounds[:, 0]
    width = bounds[:, 1] - lower
    return np.divide(values - lower, width, out=np.zeros_like(values), where=width > 0.0)


# Map unit-cube coordinates into physical coordinates
def from_unit(values: np.ndarray, bounds: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    bounds = np.asarray(bounds, dtype=np.float64)
    return bounds[:, 0] + values * (bounds[:, 1] - bounds[:, 0])

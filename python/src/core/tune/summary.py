# Shared best fit JSON output
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import math
from collections.abc import Mapping
from numbers import Real

BEST_FIT_SCHEMA_VERSION = 1


# Convert one required parameter value into a finite JSON scalar
def _finite_parameter_value(value, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f'Best-fit parameter "{name}" is not numeric')
    output = float(value)
    if not math.isfinite(output):
        raise ValueError(f'Best-fit parameter "{name}" is not finite')
    return output


# Convert one optional uncertainty or objective value into a JSON scalar
def _optional_finite_value(value) -> float | None:
    if value is None or isinstance(value, bool) or not isinstance(value, Real):
        return None
    output = float(value)
    return output if math.isfinite(output) else None


# Build the common best-fit block used by icetune, icescape and iceproxy
def build_best_fit(
    *,
    values: Mapping,
    objective_name: str,
    objective_value,
    source: str,
    uncertainties: Mapping | None = None,
    uncertainty_method: str | None = None,
) -> dict:
    if not isinstance(values, Mapping) or not values:
        raise ValueError("Best-fit values must be a non-empty parameter map")
    if not str(objective_name):
        raise ValueError("Best-fit objective name must be non-empty")
    if not str(source):
        raise ValueError("Best-fit source must be non-empty")
    if uncertainties is not None and not isinstance(uncertainties, Mapping):
        raise ValueError("Best-fit uncertainties must be a parameter map")

    parameters = {}
    for name, value in values.items():
        parameter_name = str(name)
        if not parameter_name:
            raise ValueError("Best-fit parameter name must be non-empty")
        uncertainty = None if uncertainties is None else uncertainties.get(name)
        uncertainty_value = _optional_finite_value(uncertainty)
        if uncertainty_value is not None and uncertainty_value < 0.0:
            raise ValueError(f'Best-fit uncertainty for "{parameter_name}" is negative')
        parameters[parameter_name] = {
            "value": _finite_parameter_value(value, parameter_name),
            "uncertainty": uncertainty_value,
        }

    return {
        "schema_version": BEST_FIT_SCHEMA_VERSION,
        "source": str(source),
        "objective": {"name": str(objective_name), "value": _optional_finite_value(objective_value)},
        "uncertainty_method": (None if uncertainty_method is None else str(uncertainty_method)),
        "parameters": parameters,
    }


# Extract and validate physical parameter values from one common best-fit block
def best_fit_values(best_fit: dict) -> dict[str, float]:
    if not isinstance(best_fit, dict):
        raise ValueError("best_fit is not an object")
    if best_fit.get("schema_version") != BEST_FIT_SCHEMA_VERSION:
        raise ValueError(f"best_fit schema_version must be {BEST_FIT_SCHEMA_VERSION}")
    parameters = best_fit.get("parameters")
    if not isinstance(parameters, dict) or not parameters:
        raise ValueError("best_fit.parameters is not a non-empty parameter map")

    values = {}
    for name, record in parameters.items():
        parameter_name = str(name)
        if not parameter_name:
            raise ValueError("best_fit parameter name must be non-empty")
        if not isinstance(record, dict) or "value" not in record:
            raise ValueError(f'best_fit parameter "{parameter_name}" has no value')
        values[parameter_name] = _finite_parameter_value(record["value"], parameter_name)
    return values

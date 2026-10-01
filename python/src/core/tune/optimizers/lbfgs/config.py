# Bounded L-BFGS settings for distributed simulator evaluations
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import pyjson5

from core.io.serialize import load_json_file
from core.numerics.lbfgsb import validate_line_search


# Validate all controls before starting distributed trials
def load_settings(path):
    settings = load_json_file(path, loader=pyjson5.load)
    integers = {"maxiter", "history_size", "line_search_steps"}
    positive = {"learning_rate", "line_search_decay", "line_search_armijo", "grad_step_rel", "grad_step_abs"}
    tolerances = {"gradient_tolerance", "relative_gradient_tolerance"}
    if not isinstance(settings, dict) or set(settings) != integers | positive | tolerances | {"retry_failed_search"}:
        raise ValueError("L-BFGS settings must specify all optimizer and finite-difference controls")
    for name in integers:
        if type(settings[name]) is not int or settings[name] < 1:
            raise ValueError(f"L-BFGS {name} must be a positive integer")
    for name in positive | tolerances:
        value = settings[name]
        if type(value) not in (float, int) or not math.isfinite(value) or value < 0 or (name in positive and value <= 0):
            raise ValueError(f"Invalid L-BFGS {name}")
    if not isinstance(settings["retry_failed_search"], bool):
        raise ValueError("L-BFGS retry_failed_search must be Boolean")
    validate_line_search(**{name: settings[name] for name in ("learning_rate", "line_search_decay", "line_search_armijo")})
    return settings

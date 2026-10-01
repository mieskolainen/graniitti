# icetune amplitude optimizer and preparation settings
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
import pyjson5

from core.io.serialize import load_json_file
from core.numerics.lbfgsb import validate_line_search


# Validate bounded L-BFGS controls before scheduling any physics evaluations
def load_settings(path):
    settings = load_json_file(path, loader=pyjson5.load)
    integers = {"starts", "maxiter", "history_size", "line_search_steps"}
    nonnegative = {"gradient_tolerance", "relative_gradient_tolerance", "start_relative_range", "covariance_interval"}
    steps = {"learning_rate", "line_search_decay", "line_search_armijo"}
    booleans = {"retry_failed_search", "covariance"}
    if not isinstance(settings, dict) or set(settings) != integers | nonnegative | steps | booleans | {"bank"}:
        raise ValueError("Amplitude optimizer settings must specify all bounded L-BFGS controls")
    validate_line_search(**{name: settings[name] for name in steps})
    for name in integers:
        if type(settings[name]) is not int or settings[name] < 1:
            raise ValueError(f"Amplitude optimizer {name} must be a positive integer")
    for name in nonnegative:
        value = settings[name]
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not np.isfinite(value) or value < 0:
            raise ValueError(f"Amplitude optimizer {name} must be finite and nonnegative")
    for name in booleans:
        if not isinstance(settings[name], bool):
            raise ValueError(f"Amplitude optimizer {name} must be Boolean")
    bank = settings["bank"]
    if not isinstance(bank, dict) or set(bank) != {"nodes", "closure_rtol", "batch_events", "decay_zero", "decay_steps", "mc_stat"}:
        raise ValueError("Amplitude bank settings require nodes, closure_rtol, batch_events, decay_zero, decay_steps and mc_stat")
    if bank["mc_stat"] not in {"source", "reweighted"}:
        raise ValueError("Amplitude MC statistics must be source or reweighted")
    if type(bank["nodes"]) is not int or bank["nodes"] < 2:
        raise ValueError("Amplitude bank nodes must be an integer of at least two")
    if type(bank["batch_events"]) is not int or bank["batch_events"] < 1:
        raise ValueError("Amplitude bank batch_events must be a positive integer")
    if type(bank["decay_steps"]) is not int or bank["decay_steps"] < 0:
        raise ValueError("Amplitude bank decay_steps must be a nonnegative integer")
    for name in ("closure_rtol", "decay_zero"):
        if isinstance(bank[name], bool) or not isinstance(bank[name], (int, float)) or not np.isfinite(bank[name]) or bank[name] <= 0:
            raise ValueError(f"Amplitude bank {name} must be finite and positive")
    return settings

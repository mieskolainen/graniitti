# Function call caching mechanism
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import math
from collections.abc import Mapping, Set
from functools import wraps
from threading import RLock


# Own one event-local value store and its reentrant lock
class EventCache:
    # Initialize an empty event-local cache
    def __init__(self):
        self.values = {}
        self.lock = RLock()


# Convert supported values into an exact event-view cache signature
def _freeze_strict(value):
    value_type = type(value)
    if value_type is dict:
        return ("dict", tuple((_freeze_strict(key), _freeze_strict(item)) for key, item in value.items()))
    if value_type in {tuple, list, set, frozenset}:
        return value_type.__name__, tuple(_freeze_strict(item) for item in value)
    if value is None:
        return ("none",)
    if value_type is bool:
        return "bool", value
    if value_type is int:
        return "int", value
    if value_type is float:
        if not math.isfinite(value):
            raise ValueError("Non-finite floating-point values cannot share an event cache")
        return "float", value.hex()
    if value_type is complex:
        if not math.isfinite(value.real) or not math.isfinite(value.imag):
            raise ValueError("Non-finite complex values cannot share an event cache")
        return "complex", value.real.hex(), value.imag.hex()
    if value_type in {str, bytes}:
        return value_type.__name__, value
    raise ValueError(f"Arguments of type {value_type.__name__} cannot share an event cache")


# Convert nested configuration containers into stable cache-key values
def freeze(value, *, strict: bool = False):
    if strict:
        return _freeze_strict(value)
    if isinstance(value, Mapping):
        return (
            "mapping",
            type(value),
            frozenset((freeze(key, strict=strict), freeze(item, strict=strict)) for key, item in value.items()),
        )
    if isinstance(value, (tuple, list)):
        return (type(value).__name__, type(value), tuple(freeze(item, strict=strict) for item in value))
    if isinstance(value, Set):
        return ("set", type(value), frozenset(freeze(item, strict=strict) for item in value))
    if value is None:
        return ("none",)
    if isinstance(value, bool):
        return "bool", type(value), value
    if isinstance(value, int):
        return "int", type(value), value
    if isinstance(value, float):
        if math.isnan(value):
            raise ValueError("NaN values cannot be used with cache")
        return "float", type(value), float(value).hex()
    if isinstance(value, complex):
        if math.isnan(value.real) or math.isnan(value.imag):
            raise ValueError("Complex NaN values cannot be used with cache")
        return "complex", type(value), float(value.real).hex(), float(value.imag).hex()
    if isinstance(value, (str, bytes)):
        return type(value).__name__, type(value), value
    try:
        hash(value)
    except TypeError as exc:
        raise ValueError(f"Arguments of type {type(value).__name__} cannot be used with cache") from exc
    return "value", type(value), value


# Compute the current value-semantic context of one event wrapper
def _event_context(event):
    try:
        attributes = event.__dict__
    except AttributeError:
        return None

    fields = tuple((name, value) for name, value in attributes.items() if name not in {"evt", "_icecache"})
    try:
        return type(event), id(event.evt), freeze(fields, strict=True)
    except (AttributeError, RecursionError, TypeError, ValueError):
        return None


# Copy one cached result so callers never receive the stored snapshot
def _copy_result(value):
    try:
        return True, copy.deepcopy(value)
    except Exception:
        return False, None


# Cache repeated pure function calls and event projector evaluations
def cache(function):
    cache = EventCache()
    missing = object()

    # Compute or create the cache state owned by one event view
    def event_state(event):
        state = getattr(event, "_icecache", None)
        if state is not None:
            return state
        try:
            state = event.__dict__.setdefault("_icecache", EventCache())
        except AttributeError:
            with cache.lock:
                state = getattr(event, "_icecache", None)
                if state is None:
                    state = EventCache()
                    event._icecache = state
        return state

    # Compute a cached value or evaluate and store the function result
    @wraps(function)
    def wrapper(*args, **kwargs):
        if len(args) > 0 and hasattr(args[0], "evt"):
            state = event_state(args[0])
            context = _event_context(args[0])
            try:
                call_key = (
                    function
                    if len(args) == 1 and not kwargs
                    else (function, freeze(args[1:], strict=True), freeze(kwargs, strict=True))
                )
            except (RecursionError, TypeError, ValueError):
                context = None
            if context is None:
                return function(*args, **kwargs)

            key = context, call_key
            with state.lock:
                result = state.values.get(key, missing)
                if result is not missing:
                    copied, output = _copy_result(result)
                    if copied:
                        return output
                    del state.values[key]

                output = function(*args, **kwargs)
                copied, snapshot = _copy_result(output)
                if copied:
                    state.values[key] = snapshot
                return output

        key = (freeze(args), freeze(kwargs))
        with cache.lock:
            result = cache.values.get(key, missing)
        if result is not missing:
            return result

        result = function(*args, **kwargs)
        with cache.lock:
            cached = cache.values.get(key, missing)
            if cached is missing:
                cache.values[key] = result
                return result
            return cached

    return wrapper

# Regression tests for nested and event-local iceplot cache keys
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import time
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace

from core.io.cache import EventCache
from core.io.cache import cache as event_cache


# Check equivalent nested configuration dictionaries share one cache entry
def test_cache_accepts_nested_configuration():
    calls = []

    # Record evaluations while returning a stable derived value
    @event_cache
    def read_configuration(*, dataset):
        """Return one nested dataset field."""
        calls.append(dataset)
        return dataset["histograms"][0]["scale"]

    first = {
        "histograms": [{"observable": "Rap", "scale": 1.0}],
        "states": ["Upsilon(1S)", "Upsilon(2S)"],
    }
    reordered = {
        "states": ["Upsilon(1S)", "Upsilon(2S)"],
        "histograms": [{"scale": 1.0, "observable": "Rap"}],
    }

    assert read_configuration(dataset=first) == 1.0
    assert read_configuration(dataset=reordered) == 1.0
    assert len(calls) == 1


# Check same-named projectors cannot collide in one event-local cache
def test_separates_same_named_event_projectors():
    # Build one independently decorated event projector
    def build_projector(offset):
        # Evaluate an event field with one projector-specific offset
        @event_cache
        def project(event):
            """Return one offset event value."""
            return event.value + offset

        return project

    event = SimpleNamespace(evt=object(), value=4.0)
    first = build_projector(1.0)
    second = build_projector(2.0)

    assert first(event) == 5.0
    assert second(event) == 6.0


# Check concurrent calls share one completed event-projector evaluation
def test_serializes_event_projector():
    calls = []

    # Delay one projection so concurrent callers overlap
    @event_cache
    def project(event):
        """Return one event value after a short calculation."""
        calls.append(event.value)
        time.sleep(0.01)
        return event.value

    event = SimpleNamespace(evt=object(), value=7.0)
    with ThreadPoolExecutor(max_workers=8) as pool:
        values = list(pool.map(project, [event] * 16))

    assert values == [7.0] * 16
    assert calls == [7.0]


# Check distinct event-wrapper state cannot collide in a shared cache
def test_separates_dynamic_event_view_state():
    calls = []

    # Compute one value controlled by state added by a selection
    @event_cache
    def project(event):
        """Return the selection-local tag."""
        calls.append(event.tag)
        return event.tag

    source = object()
    cache = EventCache()
    first = SimpleNamespace(evt=source, tag="first", _icecache=cache)
    second = SimpleNamespace(evt=source, tag="second", _icecache=cache)

    assert project(first) == "first"
    assert project(second) == "second"
    assert calls == ["first", "second"]


# Check temporary parameter changes cannot contaminate another event view
def test_cache_keys_current_event_view_state():
    calls = []

    # Compute one value controlled by a mutable cut parameter
    @event_cache
    def project(event):
        """Return the current cut parameter."""
        value = event.cut_param["value"]
        calls.append(value)
        return value

    source = object()
    cache = EventCache()
    first = SimpleNamespace(evt=source, cut_param={"value": 1}, _icecache=cache)
    second = SimpleNamespace(evt=source, cut_param={"value": 1}, _icecache=cache)

    first.cut_param["value"] = 2
    assert project(first) == 2
    first.cut_param["value"] = 1
    assert project(second) == 1
    assert calls == [2, 1]


# Check callers cannot modify the private cached result snapshot
def test_cache_copies_mutable_event_results():
    calls = []

    # Compute one mutable projector result
    @event_cache
    def project(event):
        """Return a fresh mutable result."""
        calls.append(event.value)
        return {"values": [event.value]}

    event = SimpleNamespace(evt=object(), value=3.0)
    first = project(event)
    first["values"].append(4.0)

    assert project(event) == {"values": [3.0]}
    assert calls == [3.0]


# Check unsupported event state conservatively disables event caching
def test_no_cache_unsupported_event_state():
    calls = []

    # Count evaluations when the event context cannot be frozen exactly
    @event_cache
    def project(event):
        """Return the number of projector evaluations."""
        calls.append(event.value)
        return len(calls)

    event = SimpleNamespace(evt=object(), value=object())

    assert project(event) == 1
    assert project(event) == 2


# Check results without safe copy support are never stored
def test_no_cache_uncopyable_event_results():
    calls = []

    # Refuse copying to exercise the conservative result path
    class Uncopyable:
        # Reject creation of a cache snapshot
        def __deepcopy__(self, _memo):
            raise TypeError("copy disabled")

    # Compute a result that cannot have a private cache snapshot
    @event_cache
    def project(event):
        """Return one uncopyable result."""
        calls.append(event.value)
        return Uncopyable()

    event = SimpleNamespace(evt=object(), value=5.0)

    assert isinstance(project(event), Uncopyable)
    assert isinstance(project(event), Uncopyable)
    assert calls == [5.0, 5.0]

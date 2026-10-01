# Lazy discovery of simulator drivers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import inspect
from functools import cache
from importlib.metadata import EntryPoint, entry_points


# Resolve driver registrations without importing their simulation dependencies
def driver_entries() -> dict:
    entries = entry_points(group="core.tune.drivers")
    names = [entry.name.upper() for entry in entries]
    if len(set(names)) != len(names):
        raise RuntimeError("Duplicate core driver registrations")
    drivers = dict(zip(names, entries, strict=True))
    for name, module, cls in (
        ("GRANIITTI", "graniitti", "GraniittiDriver"),
        ("PANDORA", "pandora", "PandoraDriver"),
    ):
        drivers.setdefault(name, EntryPoint(
            name=name.lower(), value=f"core.tune.drivers.{module}.driver:{cls}", group="core.tune.drivers"))
    return drivers


# Load only the selected simulator driver
@cache
def get_driver_type(name: str) -> type:
    entries = driver_entries()
    normalized = name.upper()
    if normalized not in entries:
        raise ValueError(f"Unknown simulator driver {name!r}; available drivers: {sorted(entries)}")
    driver = entries[normalized].load()
    if not inspect.isclass(driver) or inspect.isabstract(driver):
        raise TypeError(f"Invalid simulator driver registration: {name}")
    return driver


# Construct one selected simulator driver
def create_driver(name: str):
    return get_driver_type(name)()


# Select the supplied default without loading unrelated drivers
def default_driver_name() -> str:
    return "GRANIITTI"


# Infer a simulator from persisted metadata or its own summary recognizer
def infer_driver_name(summary: dict) -> str:
    if summary.get("simdriver"):
        return get_driver_type(str(summary["simdriver"])).driver_name()
    matches = [name for name in driver_entries() if get_driver_type(name).matches_summary(summary)]
    if len(matches) != 1:
        raise ValueError(f"Summary matches {len(matches)} simulator drivers: {matches}")
    return matches[0]

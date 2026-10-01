# Export exact pytest node IDs and durations for Condor sharding
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import os
import pathlib


# Store timing state for one Pytest process
class TimingPlugin:
    # Initialize per process timing storage
    def __init__(self, output: pathlib.Path):
        self.output = output
        self.durations: dict[str, float] = {}

    # Accumulate setup, call, and teardown duration for one node
    def pytest_runtest_logreport(self, report):
        self.durations[report.nodeid] = self.durations.get(report.nodeid, 0.0) + float(
            report.duration
        )

    # Write timings through an atomic rename after the session
    def pytest_sessionfinish(self, session, exitstatus):
        self.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.output.with_name(f"{self.output.name}.tmp.{os.getpid()}")
        temporary.write_text(
            json.dumps(self.durations, indent=2, sort_keys=True) + "\n", encoding="ascii"
        )
        temporary.replace(self.output)


# Register isolated timing state only when an output path was requested
def pytest_configure(config):
    output = os.environ.get("GRANIITTI_PYTEST_TIMINGS")
    if output:
        config.pluginmanager.register(TimingPlugin(pathlib.Path(output)), "condor-timing")


# Write every collected node ID to the requested manifest
def pytest_collection_finish(session):
    output = os.environ.get("GRANIITTI_PYTEST_COLLECTION")
    if not output:
        return
    path = pathlib.Path(output)
    nodes = [item.nodeid for item in session.items]
    path.write_text(json.dumps(nodes, indent=2) + "\n", encoding="utf-8")

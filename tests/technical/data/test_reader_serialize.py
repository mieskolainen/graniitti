# Check source-defined observables across Ray and HepMC analysis processes
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import multiprocessing
import os
import subprocess
import sys
from pathlib import Path

import pytest
from core.io import readers, steering
from core.tune.drivers.graniitti import runtime
from pyHepMC3 import HepMC3 as h
from pyHepMC3 import std

from tests.technical.support.hepmc import make_pair_event

ROOT = Path(__file__).resolve().parents[3]


# Read real STAR pion events after deserialization and return histograms through a pool
@pytest.mark.parametrize("start_method", ["fork", "spawn"])
def test_observable_process_serialization(tmp_path, start_method):
    if start_method not in multiprocessing.get_all_start_methods():
        pytest.skip(f"{start_method} is unavailable")
    pickler = pytest.importorskip("ray.cloudpickle")
    source = ROOT / "icepack/SOFTCEP/STAR_1792394/pipi/obs.py"
    dataset, path = steering.load_dataset(str(source.with_name("dataset.json")), cdir=ROOT)
    observables, _ = steering.load_observables(str(source), dataset_path=path, cdir=ROOT)
    _, observables = readers.read_hepdata(
        dataset["sets"][0], dataset["datapath"], dataset["type"], observables,
        reader=dataset["reader"], dataset_path=path, cdir=str(ROOT))
    events = tmp_path / "pipi.hepmc3"
    writer = h.WriterAscii(str(events))
    try:
        for index, mass in enumerate((0.8, 1.0, 1.2)):
            event = make_pair_event(mass=mass, pid=211).evt
            event.set_event_number(index)
            data = h.GenEventData()
            event.write_data(data)
            data.weights = std.vector_double([1.0])
            event.read_data(data)
            cross_section = h.GenCrossSection()
            cross_section.set_cross_section(3.0, 0.3, index + 1, index + 1)
            event.set_cross_section(cross_section)
            writer.write_event(event)
            assert not writer.failed()
    finally:
        writer.close()
    param = dict(hepmc3file=str(events), obs=[observables], pid=[[211, -211]],
                 cuts=[str(source.with_name("cuts.py"))], chunk_range=None, xsmode="auto",
                 header_xsection=None, verbose=False, scales=[1.0], k=0, label="STAR",
                 density=False, density_uncertainty="scaled", covariance_mode="full")
    payload = tmp_path / "input.pkl"
    payload.write_bytes(pickler.dumps(param))
    script = '''import multiprocessing, pathlib, pickle, sys
import numpy as np
from core.io import readers

param = pickle.loads(pathlib.Path(sys.argv[1]).read_bytes())
reference = readers.parallel_wrapper(param)
assert reference[0]["M"]["hdata"].entries.sum() == 3
multiprocessing.set_start_method(sys.argv[2], force=True)
with multiprocessing.Pool(2) as pool:
    result = pool.map(readers.parallel_wrapper, [param])[0]
for name, record in result[0].items():
    expected, actual = reference[0][name]["hdata"], record["hdata"]
    for field in ("counts", "errs", "bins", "binscale"):
        np.testing.assert_allclose(getattr(actual, field), getattr(expected, field))
    assert callable(record["obs"]["func"])
'''
    environment = {**os.environ, **runtime.environment(
        cdir=str(ROOT), libdir=None, python_version=f"{sys.version_info.major}.{sys.version_info.minor}")}
    subprocess.run([sys.executable, "-c", script, str(payload), start_method],
                   cwd=tmp_path, env=environment, check=True, timeout=120)

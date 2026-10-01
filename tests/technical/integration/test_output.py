# Generator output paths and automatic folder creation
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import pathlib
import subprocess
import tempfile

import pyjson5
import pytest

from tests.technical.support.hepmc import read_events


# Generate real events with default, relative and absolute output stems
@pytest.mark.parametrize("mode", ["default", "relative", "absolute"])
def test_output_paths(mode):
    """Check CLI and card output paths, parent creation and readable outputs"""
    root = pathlib.Path(__file__).resolve().parents[3]
    (root / "tmp").mkdir(exist_ok=True)
    work = pathlib.Path(tempfile.mkdtemp(prefix="test_output_", dir=root / "tmp"))
    stem = work / "new folder" / "nested" / "events"
    output = str(stem) if mode == "absolute" else "./new folder/nested/events"
    if mode == "default":
        output = work.name
    card = pyjson5.decode((root / "gencard/test.json").read_text())
    card["GENERIC"].update(OUTPUT=output, NEVENTS=2, CORES=1, HIST=1, WEIGHTED=True)
    card["SCATTERING"].update(PROCESS="GP[CON]<F> -> pi+ pi-", LOOPSCREEN=False)
    card["INTEGRATOR"].update(min_samples=1000, max_samples=2000, precision=1.0)
    card["INTEGRATOR"]["VEGAS"].update(ncall=1000, rounds=2)
    inputfile = work / "card.json"
    inputfile.write_text(json.dumps(card))
    command = [str(root / "bin/gr"), "-i", str(inputfile), "-g", "VEGAS", "-f", "hepmc3"]
    if mode == "relative":
        card["GENERIC"]["OUTPUT"] = "unused"
        inputfile.write_text(json.dumps(card))
        command += ["-o", output]
    result = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
    event_stem = root / "output" / output if mode == "default" else stem
    grid_stem = root / "vgrid" / output if mode == "default" else stem
    events = read_events(event_stem.with_suffix(".hepmc3"))
    assert len(events) == 2
    assert all(event.particles() for event in events)
    hist = json.loads(event_stem.with_suffix(".hfast").read_text())
    grid = json.loads(grid_stem.with_suffix(".vgrid").read_text())
    assert hist["h1"]
    assert grid["STAT"]["sigma"] > 0.0

# Regression tests for standalone program input and event output
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import subprocess
import tempfile
from pathlib import Path

import pytest
from pyHepMC3 import HepMC3 as hepmc

ROOT = Path(__file__).resolve().parents[3]


# Keep all generated event samples inside the project temporary directory
@pytest.fixture
def workdir():
    directory = ROOT / "tmp"
    directory.mkdir(exist_ok=True)
    return Path(tempfile.mkdtemp(prefix="program_test_", dir=directory))


# Exercise the installed executable from a normal build without changing its environment
def run_program(name, arguments, workdir):
    binary = ROOT / "bin" / name
    assert binary.is_file(), f"Missing normal build executable: {binary}"
    return subprocess.run(
        [str(binary), *map(str, arguments)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )


# Preserve PID weights when the source path contains directories and spaces
def test_data_converter_nested_input(workdir):
    source = workdir / "input data"
    source.mkdir()
    path = source / "pions.csv"
    path.write_text("0.3,0,0,0,-0.3,0,0,0,0.25\n0.4,0,0,0,-0.4,0,0,0,2.5\n")
    result = run_program("data2hepmc3", [path], workdir)
    assert result.returncode == 0, result.stderr
    reader = hepmc.ReaderAscii(str(workdir / "output" / "pions.csv.hepmc3"))
    events = []
    while not reader.failed():
        event = hepmc.GenEvent()
        reader.read_event(event)
        if not reader.failed():
            events.append(event)
    reader.close()
    assert len(events) == 2
    assert [event.event_number() for event in events] == [0, 1]
    assert [event.weights()[0] for event in events] == pytest.approx([0.25, 2.5])
    for event, momentum in zip(events, (0.3, 0.4), strict=True):
        pions = sorted((p for p in event.particles() if p.status() == 1), key=lambda p: p.pid())
        assert [p.pid() for p in pions] == [-211, 211]
        assert [p.momentum().px() for p in pions] == pytest.approx([-momentum, momentum])
        assert all(p.momentum().e() > abs(p.momentum().px()) for p in pions)


# Fail when the requested output name is already a directory
def test_data_converter_output_failure(workdir):
    path = workdir / "pions.csv"
    path.write_text("0.3,0,0,0,-0.3,0,0,0,1\n")
    (workdir / "output" / "pions.csv.hepmc3").mkdir(parents=True)
    result = run_program("data2hepmc3", [path], workdir)
    assert result.returncode > 0


# Reject an absent event stream without reporting a successful conversion
def test_hepevt_converter_missing_input(workdir):
    result = run_program("hepevt2hepmc3", [workdir / "absent.hepevt"], workdir)
    assert result.returncode > 0


# Convert missing input errors to a normal failure exit instead of aborting
def test_lhe_converter_missing_input(workdir):
    result = run_program("hepmc3tolhe", [workdir / "absent.hepmc3"], workdir)
    assert result.returncode > 0


# Detect an unwritable diffraction output before evaluating the aperture integral
def test_diffraction_output_failure(workdir):
    (workdir / "3D.ascii").mkdir()
    result = run_program("sommerfeld", [], workdir)
    assert result.returncode > 0


# Preserve the input when an existing output is a symbolic link to it
@pytest.mark.parametrize("program", ["data2hepmc3", "hepevt2hepmc3"])
def test_converter_output_alias(program, workdir):
    path = workdir / "pions.csv"
    text = "0.3,0,0,0,-0.3,0,0,0,1\n"
    path.write_text(text)
    output = workdir / "output" if program == "data2hepmc3" else workdir
    output.mkdir(exist_ok=True)
    (output / "pions.csv.hepmc3").symlink_to(path)
    result = run_program(program, [path], workdir)
    assert result.returncode > 0
    assert path.read_text() == text


# Use default transport parameters and reject additional samples instead of ignoring them
def test_transport_defaults_and_sample_count(workdir):
    path = workdir / "pions.csv"
    path.write_text("0.3,0,0,0,-0.3,0,0,0,1\n0.4,0,0,0,-0.4,0,0,0,3\n")
    converted = run_program("data2hepmc3", [path], workdir)
    assert converted.returncode == 0, converted.stderr
    events = workdir / "output" / "pions.csv.hepmc3"
    result = run_program("ot", ["-i", f"{events},{events}"], workdir)
    assert result.returncode == 0, result.stderr
    result = run_program("ot", ["-i", f"{events},{events},{events}", "-a", "1", "-r", "10"], workdir)
    assert result.returncode > 0


# Reject unsafe lattice inputs with a normal error exit instead of a signal
@pytest.mark.parametrize(
    ("program", "arguments"),
    [
        ("pdebench", [64, 128, 16, 0, 1, 0]),
        ("pdebench", [64, 0, 16, 1, 1, 0]),
        ("pdebench", [64, 128, 0, 1, 1, 0]),
        ("pdebench", [64, 127, 16, 1, 1, 1]),
        ("pathmark", [4, 1, 1.5, 0, 0.02, 1]),
        ("pathmark", [4, 0, 1.5, 2, 0.02, 1]),
        ("pathmark", [0, 1, 1.5, 2, 0.02, 1]),
        ("pathmark", [4, 1, 1.5, 4, 0.02, 1]),
        ("pathmark", [4, 1, 1.5, 2, 0.02, -1]),
    ],
)
def test_benchmark_invalid_input(program, arguments, workdir):
    result = run_program(program, arguments, workdir)
    assert result.returncode > 0, result.stderr


# Fit weighted HepMC events through both estimators with scaled response weights
@pytest.mark.skipif(not (ROOT / "bin/fitharmonic").is_file(), reason="Requires a build with ROOT support")
@pytest.mark.parametrize("estimator", ["ALGEBRAIC", "EML"])
def test_harmonic_weighted_measurement(workdir, estimator):
    source = workdir / "pions.csv"
    source.write_text("0.3,0,0,0,-0.3,0,0,0,2\n0.4,0,0,0,-0.4,0,0,0,3\n")
    converted = run_program("data2hepmc3", [source], workdir)
    assert converted.returncode == 0, converted.stderr
    events = str(workdir / "output" / "pions.csv.hepmc3")
    output = workdir / "measurement.json"
    card = {
        "schema": "GRANIITTI_HARMONIC_ANALYSIS_V1",
        "measurement": "CENTRAL", "frame": "CM", "sqrt_s": 13000.0,
        "fiducial": {"central": {"pion_eta": [-5.0, 5.0], "pion_pt": [0.0, 10.0]}},
        "axes": [
            {"coordinate": "M", "bins": 1, "min": 0.0, "max": 2.0},
            {"coordinate": "PT", "bins": 1, "min": 0.0, "max": 1.0},
            {"coordinate": "Y", "bins": 1, "min": -1.0, "max": 1.0},
        ],
        "fit": {
            "estimator": estimator, "lmax": 0, "svd_relative_cut": 0.0,
            "response_jackknife_bins": 0, "eml": {},
        },
        "responses": [{
            "name": "response", "mode": "PAIRED",
            "samples": [{"truth": events, "reco": events, "reference": "ANGULAR_FLAT", "scale": 4e307}],
        }],
        "samples": [{
            "name": "data", "type": "DATA", "response": "response", "paths": [events],
            "normalization": "CROSS_SECTION_PB", "integrated_luminosity_pb": 10.0,
            "luminosity_relative_uncertainty": 0.1,
        }],
        "output": str(output),
    }
    path = workdir / "card.json"
    path.write_text(json.dumps(card))
    fitted = run_program("fitharmonic", ["--card", path], workdir)
    assert fitted.returncode == 0, fitted.stderr
    result = json.loads(output.read_text())["samples"][0]["result"]
    assert result["data_sum_weight"] == pytest.approx(5.0)
    assert result["data_sum_weight2"] == pytest.approx(13.0)
    for level in ("angular_flat", "fiducial", "detector"):
        estimate = result[level]["cells"][0]["estimate"]
        assert estimate["coefficients"] == pytest.approx([0.5], rel=1e-5)
        assert estimate["covariance"][0] == pytest.approx([0.1325], rel=1e-5)
        assert result[level]["global_covariance"][0] == pytest.approx([0.1325], rel=1e-5)

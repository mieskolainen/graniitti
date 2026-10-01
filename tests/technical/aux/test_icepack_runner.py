# Test the icepack command line runner
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path

import pyjson5
import pytest
from core import resource

from tests.technical.support.hepmc import read_events, write_muon_events

ROOT = Path(__file__).resolve().parents[3]
RUNNER = ROOT / "icepack" / "run.sh"


# Select exact dataset cards without also scheduling their nested datasets
@pytest.mark.parametrize(("entries", "expected"), [
    (["GAMMA/continuum/dataset.json"], ["GAMMA/continuum"]),
    (["icepack/PHOTOPROD/H1_1798511/dataset.json"], ["PHOTOPROD/H1_1798511"]),
    (["GAMMA/continuum", "GAMMA/continuum/dataset.json"], [
        "GAMMA/continuum", "GAMMA/continuum/closure/highmass", "GAMMA/continuum/closure/lowmass",
    ]),
])
def test_dataset_selection(entries, expected, tmp_path):
    result = subprocess.run(
        ["bash", str(RUNNER), "--list", *entries],
        cwd=tmp_path, capture_output=True, text=True, check=True,
    )
    assert sorted(result.stdout.splitlines()) == expected


# Continue real HepMC3 analysis after a missing sample and retain the failing exit status
@pytest.mark.parametrize("failure", ["missing_input", "mc_precision"])
def test_plotting_continues_after_failed_dataset(tmp_path, failure):
    from tests.technical.support.hepmc import write_muon_events

    project = tmp_path / "project"
    packs = project / "icepack"
    packs.mkdir(parents=True)
    shutil.copyfile(RUNNER, packs / "run.sh")
    for name in ("python", "tests", "install"):
        (project / name).symlink_to(ROOT / name, target_is_directory=True)
    (project / "output").mkdir()
    write_muon_events(project / "output/good.hepmc3", energies=[0.8, 1.0, 1.2])
    if failure == "mc_precision":
        write_muon_events(project / "output/missing.hepmc3", energies=[0.8, 1.0, 1.2])
    for name in ("missing", "good"):
        entry = packs / name
        entry.mkdir()
        dataset = {
            "active": True, "type": "MC_ONLY",
            "samples": [{"name": "muons", "label": "muons", "parameters": {},
                         "gencard": str(ROOT / "icepack/GAMMA/ATLAS_1377585/mumu/gencard.json")}],
            "plot": {"normalization": "cross_section", "ratio_uncertainty": "combined",
                     "stack": False, "data_style": "hist"},
            "fit": {"normalization": "cross_section"},
            "sets": [{"name": "muons", "data": False, "pid": [13, -13],
                      "cuts": str(ROOT / "icepack/_common/cuts_inclusive.py"),
                      "obs": str(resource("analysis/observables/default.py")), "hist": [{"obs": "M"}]}],
        }
        if name == "missing" and failure == "mc_precision":
            dataset["validation"] = {"min_mc_effective_events": 4}
        (entry / "dataset.json").write_text(json.dumps(dataset))
    environment = os.environ.copy()
    for name in ("DENSITY", "INTEGRATOR", "NEVENTS", "LOOPSCREEN", "WEIGHTED"):
        environment.pop(name, None)
    outputs, reports = project / "results", project / "reports"
    result = subprocess.run(
        ["bash", str(packs / "run.sh"), "missing", "good", "--skip-generate",
         "--output-dir", str(outputs), "--report-dir", str(reports)],
        cwd=project, env=environment, capture_output=True, text=True, timeout=120, check=False,
    )
    assert result.returncode == 1, result.stdout + result.stderr
    assert (reports / "missing.json").exists() == (failure == "mc_precision")
    report = json.loads((reports / "good.json").read_text())
    sample = report["sets"][0]["observables"][0]["samples"][0]
    assert sample["mc_selected_events"] == 3
    assert sample["integral"] == pytest.approx(3.0)
    assert list((outputs / "good/plots").rglob("*.pdf"))


# Prepare a complete isolated pion-pair icepack with the normal runtime and analysis
@pytest.fixture
def pion_pack(tmp_path, monkeypatch):
    project = tmp_path / "project"
    entry = project / "icepack/pions"
    entry.mkdir(parents=True)
    shutil.copyfile(RUNNER, project / "icepack/run.sh")
    for name in ("python", "tests", "install", "modeldata", "VERSION.json"):
        (project / name).symlink_to(ROOT / name, target_is_directory=(ROOT / name).is_dir())
    (project / "bin").mkdir()
    shutil.copyfile(ROOT / "bin/gr", project / "bin/gr")
    shutil.copymode(ROOT / "bin/gr", project / "bin/gr")
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path.parent / "matplotlib"))
    for name in ("DENSITY", "INTEGRATOR", "NEVENTS", "LOOPSCREEN", "WEIGHTED"):
        monkeypatch.delenv(name, raising=False)
    card = pyjson5.loads((ROOT / "gencard/test.json").read_text())
    card["SCATTERING"].update(PROCESS="GP[CON]<F> -> pi+ pi-", RES=[])
    card["GENERIC"].update(INTEGRATOR="VEGAS", WEIGHTED=True, CORES=1, HIST=0)
    card["INTEGRATOR"].update(min_samples=1000, max_samples=2000, precision=1.0)
    card["INTEGRATOR"]["VEGAS"].update(ncall=1000, rounds=2)
    card["GENCUTS"]["<F>"]["M"] = [0.5, 2.0]
    (entry / "gencard.json").write_text(json.dumps(card))
    dataset = {
        "active": True, "type": "MC_ONLY",
        "samples": [{"name": "pions", "label": "pions", "parameters": {}, "gencard": "./gencard.json"}],
        "plot": {"normalization": "cross_section", "ratio_uncertainty": "combined", "stack": False, "data_style": "hist"},
        "fit": {"normalization": "cross_section"},
        "validation": {"nevents": 16, "loopscreen": False},
        "sets": [{"name": "pions", "data": False, "pid": [211, -211],
                  "cuts": str(ROOT / "icepack/_common/cuts_inclusive.py"),
                  "obs": str(resource("analysis/observables/default.py")), "hist": [{"obs": "M"}]}],
    }
    (entry / "dataset.json").write_text(json.dumps(dataset))
    return project, entry, dataset


# Check launcher density controls through the measured histogram integral
@pytest.mark.parametrize(("density", "plan_density"), [(1, False), (0, True), (None, True), (None, False)])
def test_density_follows_dataset_override(pion_pack, density, plan_density):
    project, entry, dataset = pion_pack
    dataset["plot"]["normalization"] = "unit_density" if plan_density else "cross_section"
    if plan_density:
        dataset["plot"]["density_uncertainty"] = "scaled"
    dataset["sets"][0]["pid"] = [13, -13]
    (entry / "dataset.json").write_text(json.dumps(dataset))
    output = project / "output/pions.hepmc3"
    output.parent.mkdir()
    write_muon_events(output, energies=[0.8, 1.0, 1.2])
    result = subprocess.run(
        ["bash", str(project / "icepack/run.sh"), "pions", "--skip-generate",
         *([] if density is None else ["--DENSITY", str(density)]),
         "--output-dir", str(project / "results")],
        cwd=project, capture_output=True, text=True, timeout=240, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads((project / "results/pions/reports/pions.json").read_text())
    sample = report["sets"][0]["observables"][0]["samples"][0]
    normalized = plan_density if density is None else bool(density)
    plot_dir = Path(report["plot_directory"])
    assert plot_dir.parent == project / "results/pions/plots"
    assert list(plot_dir.rglob("*.pdf"))
    assert sample["integral"] == pytest.approx(1.0 if normalized else 3.0)
    assert sample["mc_selected_events"] == 3


# Reject an unknown density value before starting an icepack
def test_density_override_unknown_value():
    environment = os.environ.copy()
    for name in ("DENSITY", "INTEGRATOR", "NEVENTS", "LOOPSCREEN", "WEIGHTED"):
        environment.pop(name, None)
    result = subprocess.run(
        [
            "bash",
            str(RUNNER),
            "icepack/DURHAM/partons",
            "--DENSITY",
            "2",
        ],
        cwd=ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )

    assert result.returncode == 64
    assert "DENSITY must be 0, 1, true or false" in result.stderr


# Check event count and screening overrides in real generated output
@pytest.mark.parametrize(("override", "nevents"), [(None, None), ("0", "24"), ("1", None), ("false", "24")])
def test_generation_follows_dataset_override(pion_pack, override, nevents):
    project, entry, dataset = pion_pack
    result = subprocess.run(
        ["bash", str(project / "icepack/run.sh"), "pions",
         *([] if override is None else ["--LOOPSCREEN", override]),
         *([] if nevents is None else ["--NEVENTS", nevents]),
         "--output-dir", str(project / "results")],
        cwd=project, capture_output=True, text=True, timeout=240, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    commands = list((project / "results/pions/generation").rglob("command.sh"))
    assert len(commands) == len(dataset["samples"])
    assert commands[0].with_name("generator.log").is_file()
    args = shlex.split(commands[0].read_text().splitlines()[-1])
    output = args[args.index("-o") + 1]
    expected_events = int(nevents or dataset["validation"]["nevents"])
    events = read_events(project / "output" / f"{output}.hepmc3")
    assert len(events) == expected_events
    assert all(event.cross_section().xsec() > 0.0 for event in events)
    grid = json.loads((project / "vgrid" / f"{output}.vgrid").read_text())
    assert grid["INTEGRATOR"] == "VEGAS"
    assert grid["LOOPSCREEN"] is (override == "1")
    report = json.loads((project / "results/pions/reports/pions.json").read_text())
    sample = report["sets"][0]["observables"][0]["samples"][0]
    assert sample["mc_selected_events"] == expected_events
    assert sample["integral"] > 0.0


# Repeating an analysis must preserve the previous production plots
def test_plot_productions_are_preserved(pion_pack):
    project, entry, dataset = pion_pack
    dataset["sets"][0]["pid"] = [13, -13]
    (entry / "dataset.json").write_text(json.dumps(dataset))
    output = project / "output/pions.hepmc3"
    output.parent.mkdir()
    write_muon_events(output, energies=[0.8, 1.0, 1.2])
    directories = []
    for _ in range(2):
        subprocess.run(
            ["bash", str(project / "icepack/run.sh"), "pions", "--skip-generate",
             "--output-dir", str(project / "results")],
            cwd=project, capture_output=True, text=True, timeout=240, check=True)
        report = json.loads((project / "results/pions/reports/pions.json").read_text())
        directories.append(Path(report["plot_directory"]))
    assert directories[0] != directories[1]
    assert all(list(directory.rglob("*.pdf")) for directory in directories)

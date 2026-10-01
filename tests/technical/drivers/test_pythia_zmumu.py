# Tests for the PYTHIA dimuon HepMC3 driver
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
import re
import subprocess
from pathlib import Path

import pytest

from tests.technical.support.hepmc import read_events

ROOT = Path(__file__).resolve().parents[3]
DRIVER = ROOT / "bin/pythia_zmumu_hepmc3"
WORKFLOW = ROOT / "tests/external/pythia/cases/hard_pomeron_z/run.sh"
pytestmark = pytest.mark.skipif(not DRIVER.is_file(), reason="Build the Pythia dimuon driver first")

# Write one fast process-only neutral-current Pythia card
def write_process_card(path: Path, lepton_id: int = 13) -> None:
    path.write_text(
        "! Fast process-only card for the native dimuon driver test\n"
        "Beams:idA = 2212\n"
        "Beams:idB = 2212\n"
        "Beams:eCM = 13000.0\n"
        "WeakSingleBoson:ffbar2gmZ = on\n"
        "PhaseSpace:mHatMin = 80.0\n"
        "PhaseSpace:mHatMax = 100.0\n"
        "23:onMode = off\n"
        f"23:onIfAny = {lepton_id}\n"
        "PartonLevel:all = off\n"
        "HadronLevel:all = off\n"
        "Print:quiet = on\n"
    )


# Run the native Pythia driver with the matching XML data tree
def run_driver(card: Path, output: Path, nevents: int, *arguments: str, cwd=ROOT):
    return subprocess.run(
        [str(DRIVER), str(card), str(output), str(nevents), "24680", *arguments],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )


# Parse the compact final driver summary from standard output
def parse_run_summary(output: str) -> dict:
    match = re.search(
        r"PYTHIA_ZMUMU_HEPMC3 written (\d+) events .*? tried (\d+) "
        r"selected (\d+) accepted (\d+) aborted (\d+) sigma_pb "
        r"([-+0-9.eE]+) sigma_err_pb ([-+0-9.eE]+)",
        output,
    )
    assert match is not None
    return {
        "written": int(match.group(1)),
        "tried": int(match.group(2)),
        "selected": int(match.group(3)),
        "accepted": int(match.group(4)),
        "aborted": int(match.group(5)),
        "sigma_pb": float(match.group(6)),
        "sigma_err_pb": float(match.group(7)),
    }


# Parse the process-hook accounting from standard output
def parse_filter_summary(output: str) -> dict:
    match = re.search(
        r"PYTHIA_ZMUMU_HARD_MUON_FILTER seen (\d+) passed (\d+) "
        r"kinematic_vetoed (\d+) malformed (\d+)",
        output,
    )
    assert match is not None
    return {
        "seen": int(match.group(1)),
        "passed": int(match.group(2)),
        "kinematic_vetoed": int(match.group(3)),
        "malformed": int(match.group(4)),
    }


# Check invalid values and duplicate CLI options are rejected
def test_pythia_zmumu_cli_validation():
    invalid_result = subprocess.run(
        [str(DRIVER), "card.cmnd", "output.hepmc3", "1", "--hard-muon-cuts", "maybe"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    assert invalid_result.returncode > 0, invalid_result.stdout + invalid_result.stderr

    duplicate_result = subprocess.run(
        [
            str(DRIVER),
            "card.cmnd",
            "output.hepmc3",
            "1",
            "--hard-muon-cuts",
            "on",
            "--hard-muon-cuts=off",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    assert duplicate_result.returncode > 0, duplicate_result.stdout + duplicate_result.stderr

    wrapper_result = subprocess.run(
        ["bash", str(WORKFLOW), "--hard-muon-cuts", "on", "--hard-muon-cuts=off"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    assert wrapper_result.returncode == 2


# Check filtered kinematics, counters, and post-veto absolute normalization
def test_filter_and_hepmc_xs(tmp_path):
    card = tmp_path / "zmumu.cmnd"
    filtered_path = tmp_path / "filtered.hepmc3"
    unfiltered_path = tmp_path / "unfiltered.hepmc3"
    write_process_card(card)

    filtered_result = run_driver(
        card,
        filtered_path,
        60,
        "--hard-muon-cuts=on",
    )
    unfiltered_result = run_driver(card, unfiltered_path, 60)
    assert filtered_result.returncode == 0, filtered_result.stderr
    assert unfiltered_result.returncode == 0, unfiltered_result.stderr

    filtered_summary = parse_run_summary(filtered_result.stdout)
    unfiltered_summary = parse_run_summary(unfiltered_result.stdout)
    hook_summary = parse_filter_summary(filtered_result.stdout)
    filtered_events = read_events(filtered_path)
    unfiltered_events = read_events(unfiltered_path)

    assert len(filtered_events) == len(unfiltered_events) == 60
    assert hook_summary["passed"] == filtered_summary["accepted"] == 60
    assert hook_summary["seen"] == filtered_summary["selected"]
    assert hook_summary["kinematic_vetoed"] == hook_summary["seen"] - 60
    assert hook_summary["malformed"] == 0

    previous_attempted = 0
    for index, event in enumerate(filtered_events):
        assert list(event.weights()) == pytest.approx([1.0])
        assert int(event.attribute_as_string("hard_process_muon_cuts_enabled")) == 1
        assert float(event.attribute_as_string("hard_process_mu_minus_pt_GeV")) > 20.0
        assert float(event.attribute_as_string("hard_process_mu_plus_pt_GeV")) > 20.0
        assert abs(float(event.attribute_as_string("hard_process_mu_minus_eta"))) < 2.5
        assert abs(float(event.attribute_as_string("hard_process_mu_plus_eta"))) < 2.5
        assert 80.0 <= float(event.attribute_as_string("hard_process_mumu_mass_GeV")) <= 100.0

        cross_section = event.cross_section()
        assert cross_section.get_accepted_events() == index + 1
        attempted = cross_section.get_attempted_events()
        assert attempted >= previous_attempted
        previous_attempted = attempted

    assert all(
        int(event.attribute_as_string("hard_process_muon_cuts_enabled")) == 0
        for event in unfiltered_events
    )

    terminal_cross_section = filtered_events[-1].cross_section()
    assert terminal_cross_section.xsec() == pytest.approx(filtered_summary["sigma_pb"], rel=5e-4)
    assert terminal_cross_section.xsec_err() == pytest.approx(
        filtered_summary["sigma_err_pb"], rel=5e-4
    )
    assert terminal_cross_section.get_accepted_events() == filtered_summary["accepted"]
    assert terminal_cross_section.get_attempted_events() == filtered_summary["tried"]

    cross_section_ratio = filtered_summary["sigma_pb"] / unfiltered_summary["sigma_pb"]
    filter_efficiency = hook_summary["passed"] / hook_summary["seen"]
    ratio_uncertainty = cross_section_ratio * math.sqrt(
        (filtered_summary["sigma_err_pb"] / filtered_summary["sigma_pb"]) ** 2
        + (unfiltered_summary["sigma_err_pb"] / unfiltered_summary["sigma_pb"]) ** 2
    )
    efficiency_uncertainty = math.sqrt(
        filter_efficiency * (1.0 - filter_efficiency) / hook_summary["seen"]
    )
    assert cross_section_ratio < 1.0
    assert abs(cross_section_ratio - filter_efficiency) < 4.0 * math.hypot(
        ratio_uncertainty,
        efficiency_uncertainty,
    )


# Check incompatible non-muon steering fails instead of entering an infinite veto loop
def test_filter_rejects_incompatible_card(tmp_path):
    card = tmp_path / "zee.cmnd"
    output = tmp_path / "zee.hepmc3"
    write_process_card(card, lepton_id=11)

    result = run_driver(card, output, 1, "--hard-muon-cuts", "on")

    assert result.returncode > 0, result.stdout + result.stderr


# Fail on unreadable steering without producing a valid-looking event file
@pytest.mark.parametrize("damage", ["missing", "unknown-setting"])
def test_invalid_steering(tmp_path, damage):
    card, output = tmp_path / "bad.cmnd", tmp_path / "out.hepmc3"
    if damage != "missing":
        write_process_card(card)
        card.write_text(card.read_text() + "UnknownPhysics:setting = on\n")
    result = run_driver(card, output, 1)
    assert result.returncode != 0
    assert not output.exists()


# Protect input aliases and resolve Pythia libraries and XML data from another directory
@pytest.mark.parametrize("alias", ["same", "symlink", "hardlink"])
def test_driver_paths(tmp_path, alias):
    card = tmp_path / "source.cmnd"
    write_process_card(card)
    output = tmp_path / "out.hepmc3"
    if alias == "same":
        output = card
    elif alias == "symlink":
        output.symlink_to(card)
    else:
        output.hardlink_to(card)
    before = card.read_bytes()
    result = run_driver(card, output, 1, cwd=tmp_path)
    assert result.returncode != 0
    assert card.read_bytes() == before
    output = tmp_path / "valid.hepmc3"
    result = run_driver(card, output, 1, cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(read_events(output)) == 1


# Reject seeds that Pythia would otherwise replace with a different value
@pytest.mark.parametrize("seed", ["900000001", "1000000000", "1x"])
def test_driver_seed_range(tmp_path, seed):
    card, output = tmp_path / "source.cmnd", tmp_path / "out.hepmc3"
    write_process_card(card)
    result = subprocess.run([str(DRIVER), str(card), str(output), "1", seed], capture_output=True, timeout=10)
    assert result.returncode != 0
    assert not output.exists()

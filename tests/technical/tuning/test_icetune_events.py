# Check event retention steering and shared paths for uploaded Ray runtimes
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import sys
from pathlib import Path

import pytest
from core import icetune
from core.tune import events
from core.tune.drivers.graniitti.driver import GraniittiDriver
from core.tune.runtime.ray import worker_runtime

from submit import CAMPAIGN_DIR, campaign, runtime


# Pass event retention through submission preflight, initialization and fitting
@pytest.mark.parametrize("phase", ["preflight", "init", "fit"])
def test_save_events_submission(tmp_path, monkeypatch, phase):
    catalog = campaign.load_campaign_catalog(CAMPAIGN_DIR / "campaigns.yml")
    environment = campaign.resolve(catalog, campaign_name="tune-gpom-res-con",
        environment_overrides=campaign.parse_overrides(["fit.save_events=1", f"runtime.shared_output_dir={tmp_path}"]))
    environment.update(RAY_CAMPAIGN_FINGERPRINT="a" * 64, RAY_ADDRESS="lxplus-bootstrap",
        RAY_INIT_COORD_DIR=str(tmp_path / "init"), RAY_INIT_BOOTSTRAP_CACHE=str(tmp_path / "cache"),
        RAY_INIT_RUNTIME_SHA256="b" * 64)
    command = runtime.icetune_command(environment, cdir=tmp_path, phase=phase, tunesetup_path=tmp_path / "tune.json")
    monkeypatch.setattr(sys, "argv", ["icetune", *command[3:]])
    args = icetune.parse_arguments()
    assert args.save_events
    if phase == "fit":
        assert Path(args.events_dir) == tmp_path / "runs/icetune" / environment["RUN_NAME"] / "campaigns" / ("a" * 16) / "events"


# Reject event retention when the optimizer does not generate GRANIITTI trial events
@pytest.mark.parametrize("name", ["tune-gpom-star-cms-ampfit", "tune-pandora-v0"])
def test_save_events_unsupported_driver(name):
    catalog = campaign.load_campaign_catalog(CAMPAIGN_DIR / "campaigns.yml")
    with pytest.raises(ValueError, match="save_events requires"):
        campaign.resolve(catalog, campaign_name=name, environment_overrides={"SAVE_EVENTS": "1"})


# Keep shared event storage absolute when Ray rebases an uploaded simulator checkout
def test_uploaded_runtime_event_storage(tmp_path, monkeypatch):
    head = tmp_path / "head"
    worker = tmp_path / "worker"
    worker.mkdir()
    monkeypatch.chdir(worker)
    destination = str(head / "runs/icetune/training/events")
    param = worker_runtime({"cdir": str(head), "ray_upload_runtime": True, "events_dir": destination}, GraniittiDriver())
    assert param["cdir"] == str(worker)
    assert param["events_dir"] == destination


# Verify real HepMC3 samples and full manifests behind compact event directory names
@pytest.mark.parametrize("damage", ["sample", "manifest"])
def test_compact_event_storage(tmp_path, damage):
    from tests.technical.support.hepmc import read_events, write_muon_events

    source = tmp_path / "source"
    source.mkdir()
    sample = Path(write_muon_events(source / "events.hepmc3"))
    options = dict(source=source, destination=tmp_path / "shared", trial_id="trial", metadata={"run": "test"})
    result = events.publish_samples(**options)
    manifest = Path(result["manifest"])
    saved = manifest.parent / sample.name
    assert manifest.parent.name == result["attempt"][:16]
    assert len(read_events(saved)) == len(read_events(sample)) > 0
    assert saved.read_bytes() == sample.read_bytes()
    assert events.publish_samples(**options) == result
    if damage == "sample":
        saved.write_bytes(b"incomplete transfer")
    else:
        content = json.loads(manifest.read_text())
        content["metadata"]["run"] = "another run"
        manifest.write_text(json.dumps(content))
    with pytest.raises(events.EventTransferError, match="Incomplete event sample|Conflicting event manifest"):
        events.publish_samples(**options)

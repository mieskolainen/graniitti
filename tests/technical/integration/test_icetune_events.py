# Verify retained icetune event samples and shared output transfer with real generation
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import math
from pathlib import Path

import pytest
from core.tune import core, events
from core.tune.drivers.graniitti.driver import GraniittiDriver
from core.tune.drivers.graniitti.tunesetup import card

from tests.technical.integration.test_ampfit import dataset_card
from tests.technical.support.hepmc import read_events

ROOT = Path(__file__).resolve().parents[3]


# Generate an actual trial with the ordinary GRANIITTI driver and VEGAS
@pytest.fixture(scope="module")
def retained_trial(tmp_path_factory):
    directory = tmp_path_factory.mktemp("icetune_events")
    driver = GraniittiDriver()
    dataset_path = Path(dataset_card(directory, "STAR_1792394/KK"))
    dataset = json.loads(dataset_path.read_text())
    dataset["samples"][0]["parameters"]["NUMERICS.json:NUMERICS_VEGAS.automatic_convergence"] = False
    dataset_path.write_text(json.dumps(dataset))
    cards = [card(str(dataset_path), nevents=32, loopscreen=False,
                  xsmode="sample", swap_process="GP[RES+CON]<F>", integrator="VEGAS")]
    driver.init_data(str(directory / "run"), cards, "default", str(ROOT), pickle_dump=False)
    key = "RES|f0_980:GP:phi"
    initial = driver.get_initial_param({key: {}}, {}, str(ROOT))
    config = {key: initial[key] + 0.1}
    param = dict(cdir=str(ROOT), run_name=directory.name, datacards=cards, aux_param_space={},
                 mc_steer=dict(tune_default="TUNE0", tunesetup_name=directory.name), save_events=True,
                 chunksize=10000, processes=1, max_t=180, rngseed=31415, cost="chi2",
                 cost_rho="quadratic", cost_avg="global-mean", optimization={})
    outputs = core.evaluate_trial_outputs(config=config, param=param, simdriver=driver, trial_id=directory.name)
    assert outputs["error"] is None, outputs["error"]
    core.require_finite_trial_cost(outputs=outputs, cost_key="chi2")
    return driver, param, outputs


# Preserve physical events, weights, cross sections, cards and parameter labels through cleanup
def test_saved_trial_events(retained_trial, tmp_path):
    driver, param, outputs = retained_trial
    source = Path(outputs["event_samples"])
    core.cleanup_trial_outputs(simdriver=driver, param=param, tunename=outputs["tunename"])
    metadata = {"config": outputs["config"], "card_config": outputs["card_config"]}
    published = events.publish_samples(source=source, destination=tmp_path, trial_id="trial1", metadata=metadata)
    manifest = json.loads(Path(published["manifest"]).read_text())
    saved = Path(published["manifest"]).parent
    events.verify_samples(saved, manifest)
    assert manifest["metadata"]["config"] == outputs["config"]
    assert (saved / "modeldata/RES/f0_980.json").is_file()
    sample = json.loads((saved / "sample_0000/sample.json").read_text())
    assert sample["rngseed"] > 0 and sample["weighted"] is True
    assert sample["nevents"] == 32
    generated = read_events(saved / "sample_0000/events.hepmc3")
    assert len(generated) == 32
    assert all(event.particles() and event.weights() and math.isfinite(event.weights()[0]) for event in generated)
    assert all(event.cross_section() is not None and math.isfinite(event.cross_section().xsec()) for event in generated)
    assert events.publish_samples(source=source, destination=tmp_path, trial_id="trial1", metadata=metadata) == published
    assert (source / "sample_0000/events.hepmc3").stat().st_ino != (saved / "sample_0000/events.hepmc3").stat().st_ino
    (saved / "sample_0000/events.hepmc3").write_bytes(b"interrupted")
    with pytest.raises(events.EventTransferError, match="Incomplete event sample"):
        events.publish_samples(source=source, destination=tmp_path, trial_id="trial1", metadata=metadata)


# Retry a failed shared destination without publishing a completed manifest
def test_event_transfer_retry(retained_trial, tmp_path):
    _, _, outputs = retained_trial
    source = Path(outputs["event_samples"])
    destination = tmp_path / "shared"
    destination.write_text("unavailable")
    with pytest.raises(events.EventTransferError):
        events.publish_samples(source=source, destination=destination, trial_id="trial2", metadata={})
    assert not list(tmp_path.rglob("manifest.json"))
    destination.rename(tmp_path / "shared._old")
    result = events.publish_samples(source=source, destination=destination, trial_id="trial2", metadata={})
    assert Path(result["manifest"]).is_file()


# Preserve the default removal of ordinary temporary trial events
def test_event_retention_disabled(retained_trial):
    driver, param, outputs = retained_trial
    plain = {**param, "save_events": False}
    result = core.evaluate_trial_outputs(config=outputs["config"], param=plain, simdriver=driver, trial_id="no-events")
    assert result["error"] is None, result["error"]
    assert result["event_samples"] is None
    core.cleanup_trial_outputs(simdriver=driver, param=plain, tunename=result["tunename"])

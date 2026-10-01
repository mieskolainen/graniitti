# Test icetune parameter push workflows
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import cmath
import json
import re
import runpy
import shutil
import sys
from argparse import Namespace
from pathlib import Path

import pytest
from core.io.serialize import load_json_file

ROOT = Path(__file__).resolve().parents[3]

json5 = pytest.importorskip("pyjson5")
from core.tune import push as icetune_push
from core.tune.drivers.graniitti import driver as graniitti
from core.tune.drivers.pandora import driver as pandora
from core.tune.tunesetup import load_tunesetup


# Copy one pristine tune and record every source card byte
def _copy_tune(tmp_path):
    cdir = tmp_path / "graniitti"
    target = cdir / "modeldata" / "TUNE0"
    target.parent.mkdir(parents=True)
    shutil.copytree(ROOT / "modeldata" / "TUNE0", target)
    before = {path.relative_to(target): path.read_bytes() for path in target.rglob("*.json")}
    return cdir, target, before


# Build the shared fitted card config used by push tests
def _card_summary(config, parameters=None):
    return {
        "config": config,
        "card_config": {
            "schema_version": 1,
            "parameters": config if parameters is None else parameters,
            "tables": {},
        },
    }


# Check scalar replacement supports comments, unquoted keys and trailing commas
def test_render_json5_scalars_unrelated_text():
    source = """// card comment
{
  alpha: 1.00, // inline comment
  table: [2, 3,],
}
"""
    desired = json5.loads(source)
    desired["alpha"] = 1.25
    desired["table"][1] = 4

    rendered, changed = icetune_push.render_json5_scalars(source, desired)

    assert changed == [("alpha",), ("table", 1)]
    assert rendered == source.replace("1.00", "1.25").replace("3,", "4,")


# Check an explicitly allowed optional object can be added without reformatting its parent
def test_render_json5_scalars_adds_allowed_object():
    source = """{
  model: {
    phi: 0.0
  }
}
"""
    desired = json5.loads(source)
    desired["model"]["FF_prod"] = {
        "type": "gaussian",
        "norm": "pole",
        "Lambda2": 2.0,
    }

    rendered, changed = icetune_push.render_json5_scalars(
        source,
        desired,
        allowed_object_additions={("model", "FF_prod")},
    )

    assert changed == [("model", "FF_prod")]
    assert rendered == source.replace(
        "    phi: 0.0\n",
        '    phi: 0.0,\n    "FF_prod": {"type": "gaussian", "norm": "pole", "Lambda2": 2.0}\n',
    )


# Check a successful multi-file transaction updates all targets and keeps modes
def test_atomic_write_texts_commits_all_targets(tmp_path):
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    first.write_text("old first\n", encoding="utf-8")
    first.chmod(0o640)

    icetune_push.atomic_write_texts({first: "new first\n", second: "new second\n"})

    assert first.read_text(encoding="utf-8") == "new first\n"
    assert second.read_text(encoding="utf-8") == "new second\n"
    assert first.stat().st_mode & 0o777 == 0o640
    assert {path.name for path in tmp_path.iterdir()} == {"first.json", "second.json"}


# Check a mid-commit failure restores prior bytes, absence and sidecar cleanup
def test_atomic_write_texts_rolls_back_all_targets(tmp_path, monkeypatch):
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    missing = tmp_path / "missing.json"
    first.write_bytes(b"\xffold first\n")
    second.write_bytes(b"old second\n")
    first.chmod(0o640)
    original = {first: first.read_bytes(), second: second.read_bytes()}
    real_replace = icetune_push.os.replace
    calls = 0

    # Fail only the second target replacement so rollback operations can run
    def fail_second(source, target):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected commit failure")
        return real_replace(source, target)

    monkeypatch.setattr(icetune_push.os, "replace", fail_second)
    with pytest.raises(OSError, match="injected commit failure"):
        icetune_push.atomic_write_texts(
            {first: "new first\n", second: "new second\n", missing: "new third\n"}
        )

    assert {path: path.read_bytes() for path in original} == original
    assert first.stat().st_mode & 0o777 == 0o640
    assert not missing.exists()
    assert {path.name for path in tmp_path.iterdir()} == {"first.json", "second.json"}


# Check push input, target and optional driver options parse as one operation
@pytest.mark.parametrize(
    "values",
    [
        ["figs/icetune/run/summary.json", "modeldata/TUNE0"],
        ["figs/icetune/run/summary.json", "modeldata/TUNE0", '{"cards":["GENERAL.json"]}'],
    ],
)
def test_icetune_push_arguments(monkeypatch, values):
    module = runpy.run_module("core.icetune")
    monkeypatch.setattr(sys, "argv", ["icetune", "--push", *values])
    args = module["parse_arguments"]()
    assert args.push == values
    assert args.simdriver is None


# Check the push arguments reject value counts other than two or three
@pytest.mark.parametrize(
    "values",
    [
        ["summary.json"],
        ["summary.json", "modeldata/TUNE0", "{}", "extra"],
    ],
)
def test_push_invalid_arg_count(monkeypatch, values):
    module = runpy.run_module("core.icetune")
    monkeypatch.setattr(sys, "argv", ["icetune", "--push", *values])

    with pytest.raises(SystemExit):
        module["parse_arguments"]()


# Check driver push options require a strict JSON object
@pytest.mark.parametrize("value", ["[]", "null", "not-json"])
def test_push_driver_options_reject_non_objects(value):
    with pytest.raises(ValueError, match="JSON object"):
        icetune_push.parse_driver_options(value)


# Check driver push options preserve nested driver-owned values
def test_push_driver_options_parse_object():
    value = '{"cards":["GENERAL.json","RES*.json"]}'

    assert icetune_push.parse_driver_options(value) == {"cards": ["GENERAL.json", "RES*.json"]}


# Keep data covariance controls on the GRANIITTI driver parser only
def test_icetune_driver_specific_arguments(monkeypatch):
    module = runpy.run_module("core.icetune")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "icetune",
            "--simdriver",
            "GRANIITTI",
            "--mc_correlation_events",
            "123",
            "--mc_correlation_weighting",
            "unweighted",
            "--data_covariance_mode",
            "diagonal",
        ],
    )
    graniitti_args = module["parse_arguments"]()
    assert graniitti_args.mc_correlation_events == 123
    assert graniitti_args.mc_correlation_weighting == "unweighted"
    assert graniitti_args.data_covariance_mode == "diagonal"

    monkeypatch.setattr(
        sys,
        "argv",
        ["icetune", "--simdriver", "GRANIITTI", "--data_covariance_mode", "full"],
    )
    assert module["parse_arguments"]().data_covariance_mode == "full"

    monkeypatch.setattr(
        sys,
        "argv",
        ["icetune", "--simdriver", "GRANIITTI", "--data_covariance_mode", "on"],
    )
    with pytest.raises(SystemExit):
        module["parse_arguments"]()

    monkeypatch.setattr(sys, "argv", ["icetune", "--simdriver", "PANDORA"])
    pandora_args = module["parse_arguments"]()
    assert not hasattr(pandora_args, "data_covariance_mode")
    assert not hasattr(pandora_args, "mc_correlation_events")
    assert not hasattr(pandora_args, "mc_correlation_weighting")


# Check optimizer best-fit outputs map to the shared push configuration
@pytest.mark.parametrize(
    ("payload", "source"),
    [
        (
            {
                "best_fit": {
                    "schema_version": 1,
                    "source": "surrogate",
                    "objective": {"name": "Z", "value": 123.0},
                    "uncertainty_method": "torch_autograd_hessian",
                    "parameters": {"REGGE|omega.GP": {"value": 0.8, "uncertainty": 0.1}},
                },
                "config": {"REGGE|omega.GP": 0.7},
            },
            "best_fit.parameters",
        ),
        ({"config": {"REGGE|omega.GP": 0.8}}, "config"),
        ({"best_config": {"REGGE|omega.GP": 0.8}, "config": {"REGGE|omega.GP": 0.7}}, "best_config"),
        ({"objective_summary": {"config": {"REGGE|omega.GP": 0.8}}}, "objective_summary.config"),
        ({"best_result": {"config": {"REGGE|omega.GP": 0.8}}}, "best_result.config"),
        (
            {"surrogate_best": {"REGGE|omega.GP": {"value": 0.8, "uncertainty": 0.1}}},
            "surrogate_best",
        ),
        (
            {"optimized_best": {"Z": 123.0, "theta": {"REGGE|omega.GP": 0.8}}},
            "optimized_best.theta",
        ),
    ],
)
def test_load_summary_normalizes_optimizer_outputs(tmp_path, payload, source):
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(json.dumps(payload), encoding="utf-8")

    summary = icetune_push.load_summary(summary_path)

    assert summary["config"] == {"REGGE|omega.GP": 0.8}
    assert summary["_push_parameter_source"] == source


# Check malformed or non-finite best-fit values are rejected before card rendering
@pytest.mark.parametrize(
    "payload",
    [
        {"best_fit": {"schema_version": 1, "parameters": {}}},
        {"surrogate_best": {"REGGE|omega.GP": {"uncertainty": 0.1}}},
        {"optimized_best": {"Z": 123.0}},
        {"optimized_best": {"theta": {"REGGE|omega.GP": float("nan")}}},
        {"best_config": {"PWRAP_cal": float("inf")}},
        {"objective_summary": {"config": []}},
        {"best_result": {}},
    ],
)
def test_load_summary_invalid_optimizer_outputs(tmp_path, payload):
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError):
        icetune_push.load_summary(summary_path)


# Read direct parameter baselines while keeping push input validation strict
@pytest.mark.parametrize("config", [{}, {"PWRAP_cal": 1.2}])
@pytest.mark.parametrize("wrapped", [False, True])
def test_load_baseline(tmp_path, config, wrapped):
    source = tmp_path / "baseline.json"
    source.write_text(json.dumps({"config": config} if wrapped else config))
    assert icetune_push.load_summary(source, baseline=True)["config"] == config
    if not wrapped or not config:
        with pytest.raises(ValueError):
            icetune_push.load_summary(source)


# Check confirmation retries ambiguous input and cancels on no or closed input
@pytest.mark.parametrize(("answers", "prompt_count"), [(["maybe", "no"], 2), ([], 1)])
def test_push_confirmation_cancels_safely(monkeypatch, answers, prompt_count):
    responses = iter(answers)
    prompts = []

    # Model explicit answers followed by a closed standard input stream
    def response(prompt):
        prompts.append(prompt)
        try:
            return next(responses)
        except StopIteration:
            raise EOFError from None

    monkeypatch.setattr("builtins.input", response)
    accepted = icetune_push.confirm_change_table([("REGGE|omega.GP", 0.7, 0.8)])
    assert not accepted
    assert len(prompts) == prompt_count


# Check push preview reports changed named couplings without duplicate card rows
def test_push_preview_hides_duplicate_angular_tables():
    base = "RES|f0_980:GP:g_ls"
    old = {
        "parameters": {
            f"{base}(0,0)@MAG": 0.4,
            f"{base}(0,0)@PHASE": 0.0,
        },
        "tables": {base: {"rows": [[0, 0, 0.4, 0.0]]}},
    }
    new = {
        "parameters": {
            f"{base}(0,0)@MAG": 0.2,
            f"{base}(0,0)@PHASE": 0.0,
        },
        "tables": {base: {"rows": [[0, 0, 0.2, 0.0]]}},
    }

    assert icetune_push.card_config_rows(old, new) == [
        (f"{base}(0,0)@MAG", 0.4, 0.2)
    ]


# Check a decoded table remains visible when it has no named card parameters
def test_push_preview_changed_table_only_values():
    base = "RES|f2_1270:MP:polarization.rho_population"
    old = {
        "parameters": {},
        "tables": {base: {"rho_mag": [[0.5, 0.0], [0.0, 0.5]]}},
    }
    new = {
        "parameters": {},
        "tables": {base: {"rho_mag": [[0.4, 0.0], [0.0, 0.6]]}},
    }

    assert icetune_push.card_config_rows(old, new) == [
        (f"{base}:rho_mag[0][0]", 0.5, 0.4),
        (f"{base}:rho_mag[1][1]", 0.5, 0.6),
    ]


# Check the icetune push runner passes an icescape optimum to the selected driver
def test_push_runs_with_icescape_params(tmp_path, monkeypatch):
    module = runpy.run_module("core.icetune")
    summary_path = tmp_path / "parameters.json"
    summary_path.write_text(
        json.dumps({"surrogate_best": {"REGGE|omega.GP": {"value": 0.8, "uncertainty": 0.1}}}),
        encoding="utf-8",
    )
    received = {}

    class RecordingDriver:
        # Record and confirm the normalized summary without modifying a steering card
        def push_parameters(self, *, summary, target_path, cdir, options, confirm):
            received.update(summary)
            assert target_path == "modeldata/TUNE0"
            assert cdir == str(tmp_path)
            assert options == {"cards": ["GENERAL.json"]}
            rows = [("REGGE|omega.GP", 0.7, summary["config"]["REGGE|omega.GP"])]
            assert confirm(rows)
            return rows

    monkeypatch.setitem(
        module["run_push"].__globals__,
        "create_driver",
        lambda _name: RecordingDriver(),
    )
    args = Namespace(
        push=["parameters.json", "modeldata/TUNE0", '{"cards":["GENERAL.json"]}'],
        cdir=str(tmp_path),
        simdriver=None,
    )
    monkeypatch.setattr("builtins.input", lambda _prompt: "yes")

    assert module["run_push"](args)

    assert received["config"] == {"REGGE|omega.GP": 0.8}
    assert received["_push_parameter_source"] == "surrogate_best"
    assert args.simdriver == "GRANIITTI"


# Check the icetune push runner reports cancellation before applying driver changes
def test_push_decline_stops_before_apply(tmp_path, monkeypatch):
    module = runpy.run_module("core.icetune")
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(
        json.dumps({"config": {"REGGE|omega.GP": 0.8}}),
        encoding="utf-8",
    )
    applied = False

    class RecordingDriver:
        # Model the required driver confirmation boundary before applying changes
        def push_parameters(self, *, summary, target_path, cdir, options, confirm):
            nonlocal applied
            assert options == {}
            rows = [("REGGE|omega.GP", 0.7, summary["config"]["REGGE|omega.GP"])]
            if not confirm(rows):
                raise icetune_push.PushCancelled
            applied = True
            return rows

    monkeypatch.setitem(
        module["run_push"].__globals__,
        "create_driver",
        lambda _name: RecordingDriver(),
    )
    monkeypatch.setattr("builtins.input", lambda _prompt: "no")
    args = Namespace(
        push=["summary.json", "modeldata/TUNE0"],
        cdir=str(tmp_path),
        simdriver=None,
    )

    assert not module["run_push"](args)

    assert not applied


# Check GRANIITTI uses its driver decoder but changes only the fitted JSON5 token
def test_graniitti_push_tune_card_format(tmp_path):
    cdir, target, before = _copy_tune(tmp_path)
    summary = _card_summary({"REGGE|omega.GP": 0.8})

    rows = graniitti.GraniittiDriver().push_parameters(
        summary=summary,
        target_path="modeldata/TUNE0",
        cdir=str(cdir),
        options={},
        confirm=lambda _rows: True,
    )

    after = {path.relative_to(target): path.read_bytes() for path in target.rglob("*.json")}
    general_before = before.pop(Path("GENERAL.json")).decode()
    general_after = after.pop(Path("GENERAL.json")).decode()
    omega_before = json5.loads(general_before)["PARAM_REGGE"]["omega"]["GP"]
    omega_token = r'("omega"\s*:\s*\{[^}]*"GP"\s*:\s*)[+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?'
    expected, count = re.subn(omega_token, r'\g<1>0.8', general_before)
    assert count == 1
    assert after == before
    assert general_after == expected
    assert rows == [("REGGE|omega.GP", pytest.approx(omega_before), 0.8)]


# Check GRANIITTI push options restrict updates to selected card families
def test_graniitti_push_updates_only_selected_cards(tmp_path):
    cdir, target, before = _copy_tune(tmp_path)
    resonance = json5.loads((target / "RES" / "f0_980.json").read_text(encoding="utf-8"))
    lambda2 = float(resonance["PARAM_RES"]["MODELS"]["TP"]["FF_prod"]["Lambda2"]) + 0.125
    config = {
        "REGGE|omega.GP": 0.8,
        "RES|f0_980:TP:FF_prod.Lambda2": lambda2,
        "DECAY|9010221:[321,-321]:zeta.GP@PHASE": 0.25,
    }
    parameters = {name.removesuffix("@PHASE"): value for name, value in config.items()}
    summary = _card_summary(config, parameters)

    rows = graniitti.GraniittiDriver().push_parameters(
        summary=summary,
        target_path="modeldata/TUNE0",
        cdir=str(cdir),
        options={"cards": ["RES*.json", "DECAYS.json"]},
        confirm=lambda _rows: True,
    )

    after = {path.relative_to(target): path.read_bytes() for path in target.rglob("*.json")}
    changed = {path for path in before if before[path] != after[path]}
    assert changed == {Path("RES/f0_980.json"), Path("DECAYS.json")}
    assert {row[0] for row in rows} == {
        "RES|f0_980:TP:FF_prod.Lambda2",
        "DECAY|9010221:[321,-321]:zeta.GP",
    }


# Check TP continuum couplings publish only into CON_TP.json
def test_graniitti_push_updates_tensor_con_card(tmp_path):
    cdir, target, before = _copy_tune(tmp_path)
    tensor = json5.loads((target / "CON_TP.json").read_text(encoding="utf-8"))
    key = "CON_TP|995:[211,211]:g_tensor"
    value = 0.85 * float(tensor["995"]["[211,211]"]["g_tensor"][0])
    summary = _card_summary({key: value})

    rows = graniitti.GraniittiDriver().push_parameters(
        summary=summary,
        target_path="modeldata/TUNE0",
        cdir=str(cdir),
        options={"cards": ["CON_TP.json"]},
        confirm=lambda _rows: True,
    )

    after = {path.relative_to(target): path.read_bytes() for path in target.rglob("*.json")}
    changed = {path for path in before if before[path] != after[path]}
    assert changed == {Path("CON_TP.json")}
    assert rows == [(key, pytest.approx(tensor["995"]["[211,211]"]["g_tensor"][0]), value)]


# Check push updates the explicit continuum form factor in the current cards
def test_graniitti_push_updates_con_form_factor(tmp_path):
    cdir, target, before = _copy_tune(tmp_path)
    path = target / "CON_GP.json"
    expected = load_json_file(path, loader=json5.load)["990"]["[111,111]"]["FF_offshell"]
    expected["b"] = pytest.approx(1.25)
    key = "CON_GP|990:[111,111]:FF_offshell.b"
    summary = _card_summary({key: 1.25})

    rows = graniitti.GraniittiDriver().push_parameters(
        summary=summary,
        target_path="modeldata/TUNE0",
        cdir=str(cdir),
        options={"cards": ["CON_GP.json"]},
        confirm=lambda _rows: True,
    )

    card = json5.loads((target / "CON_GP.json").read_text(encoding="utf-8"))
    assert card["990"]["[111,111]"]["FF_offshell"] == expected
    assert rows[0][0] == key
    assert before[Path("CON_GP.json")] != (target / "CON_GP.json").read_bytes()


# Check fitted GP Pomeron and Reggeon couplings survive summary decoding and push
@pytest.mark.parametrize("suffix", ["", "-orear"])
@pytest.mark.parametrize("decoded", [False, True])
def test_graniitti_push_gp_ampfit_couplings(tmp_path, suffix, decoded):
    cdir, target, before = _copy_tune(tmp_path)
    study = load_tunesetup(cdir=str(ROOT), simdriver="graniitti", tune_default="TUNE0",
        name=f"submit/tunecards/graniitti/tune-gpom-star-cms-ampfit{suffix}.json")
    scales = {f"CON_GP|{exchange}:[{pdg},{pdg}]/opposite:helicity(0,0,0)@NORM": scale
              for exchange, scale in ((990, 0.85), (9910, 1.15)) for pdg in (211, 321)}
    assert scales.keys() <= study.param_space.keys()
    parameters = {key: study.param_space[key] for key in scales}
    driver = graniitti.GraniittiDriver()
    initial = driver.get_initial_param(parameters, {}, str(cdir))
    config = {key: initial[key] * scale for key, scale in scales.items()}
    summary = {"config": config}
    if decoded:
        driver.create_steering_card(param_space=config, tunename="FIT", cdir=str(cdir))
        summary["card_config"] = driver._decoded_card_config(config=config, path=str(cdir / "modeldata/FIT"))
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(json.dumps(summary))
    rows = driver.push_parameters(summary=icetune_push.load_summary(summary_path),
        target_path=str(target), cdir=str(cdir), options={"cards": ["CON_GP.json"]}, confirm=lambda _: True)

    expected = json5.loads(before[Path("CON_GP.json")].decode())
    pushed = json5.loads((target / "CON_GP.json").read_text())
    for exchange, scale in (("990", 0.85), ("9910", 1.15)):
        for pair in ("[211,211]", "[321,321]"):
            original = expected[exchange][pair]["opposite"]["helicity"]
            updated = pushed[exchange][pair]["opposite"]["helicity"]
            for old, new in zip(original, updated, strict=True):
                assert cmath.rect(new[-2], new[-1]) == pytest.approx(scale * cmath.rect(old[-2], old[-1]))
            expected[exchange][pair]["opposite"]["helicity"] = updated
    assert pushed == expected
    assert len(rows) == len(scales)
    assert driver.get_initial_param(parameters, {}, str(cdir)) == pytest.approx(config)
    assert all((target / path).read_bytes() == value for path, value in before.items() if path != Path("CON_GP.json"))


# Check GRANIITTI rejects unknown card family selectors before staging
def test_graniitti_push_unknown_card_selector(tmp_path):
    summary = {"config": {"REGGE|omega.GP": 0.8}}

    with pytest.raises(ValueError, match="Unknown GRANIITTI push cards"):
        graniitti.GraniittiDriver().push_parameters(
            summary=summary,
            target_path=str(tmp_path),
            cdir=str(tmp_path),
            options={"cards": ["UNKNOWN.json"]},
            confirm=lambda _rows: True,
        )


# Check declining a GRANIITTI push leaves every target card byte unchanged
def test_graniitti_push_decline_all_target_cards(tmp_path):
    cdir, target, before = _copy_tune(tmp_path)
    summary = _card_summary({"REGGE|omega.GP": 0.8})
    prepared_rows = []

    # Capture the preview while declining the write
    def decline(rows):
        prepared_rows.extend(rows)
        return False

    with pytest.raises(icetune_push.PushCancelled):
        graniitti.GraniittiDriver().push_parameters(
            summary=summary,
            target_path="modeldata/TUNE0",
            cdir=str(cdir),
            options={},
            confirm=decline,
        )

    after = {path.relative_to(target): path.read_bytes() for path in target.rglob("*.json")}
    omega_before = json5.loads(before[Path("GENERAL.json")].decode())["PARAM_REGGE"]["omega"]["GP"]
    assert after == before
    assert prepared_rows == [("REGGE|omega.GP", pytest.approx(omega_before), 0.8)]


# Check Pandora updates catalogued XML and Python values without reformatting either file
def test_pandora_push_xml_and_python_format(tmp_path):
    target = tmp_path / "pandora"
    target.mkdir()
    xml_path = target / "PandoraSettings.xml"
    python_path = target / "run_reco_pandora.py"
    xml_source = """<?xml version="1.0"?>
<!-- keep XML comment -->
<pandora>
  <algorithm type="Example">
    <Gain>1.00</Gain> <!-- keep inline comment -->
    <Enabled>true</Enabled>
  </algorithm>
</pandora>
"""
    python_source = """# keep Python comment
pandora.Parameters = {
    "Scale": ["1.25"],  # keep inline comment
}
"""
    xml_path.write_text(xml_source, encoding="utf-8")
    python_path.write_text(python_source, encoding="utf-8")
    catalog = pandora.build_catalog(
        xml_path,
        optional_xml_numeric_defaults={"Example": {"Threshold": 3.0}},
        wrapper_defaults={"Scale": [1.25]},
    )
    gain = next(entry for entry in catalog["xml_numeric"] if entry["tag"] == "Gain")
    threshold = next(entry for entry in catalog["xml_numeric"] if entry["tag"] == "Threshold")
    enabled = next(entry for entry in catalog["xml_boolean"] if entry["tag"] == "Enabled")
    wrapper = catalog["wrapper_numeric"][0]
    config = {
        gain["key"]: 2.5,
        threshold["key"]: 4.5,
        enabled["key"]: 0.0,
        wrapper["key"]: 1.5,
    }
    summary = {
        "config": config,
        "card_config": pandora._card_config(config, catalog),
    }

    rows = pandora.PandoraDriver().push_parameters(
        summary=summary,
        target_path=str(target),
        cdir=str(tmp_path),
        options={},
        confirm=lambda _rows: True,
    )

    assert xml_path.read_text(encoding="utf-8") == xml_source.replace(
        "<Gain>1.00</Gain>",
        "<Gain>2.5</Gain>",
    ).replace("<Enabled>true</Enabled>", "<Enabled>false</Enabled>").replace(
        "\n  </algorithm>",
        "\n    <Threshold>4.5</Threshold>\n  </algorithm>",
    )
    assert python_path.read_text(encoding="utf-8") == python_source.replace(
        '"Scale": ["1.25"]',
        '"Scale": ["1.5"]',
    )
    assert rows == [
        (gain["key"], 1.0, 2.5),
        (threshold["key"], 3.0, 4.5),
        (enabled["key"], 1.0, 0.0),
        (wrapper["key"], 1.25, 1.5),
    ]


# Keep baseline files unchanged on decline and retain untuned parameters on approval
@pytest.mark.parametrize("existing", [False, True])
def test_pandora_push_baseline(tmp_path, existing):
    target = tmp_path / "baseline" / "parameters.json"
    old = {"PWRAP_cal": 1.0, "PWRAP_other": 2.0} if existing else {}
    before = json.dumps({"config": old})
    if existing:
        target.parent.mkdir()
        target.write_text(before)
    summary = icetune_push.normalize_summary({"config": {"PWRAP_cal": 1.2}}, tmp_path / "summary.json")
    args = dict(summary=summary, target_path=str(target), cdir=str(tmp_path), options={})
    driver = pandora.PandoraDriver()
    with pytest.raises(icetune_push.PushCancelled):
        driver.push_parameters(**args, confirm=lambda _rows: False)
    assert target.read_text() == before if existing else not target.parent.exists()
    rows = driver.push_parameters(**args, confirm=lambda _rows: True)
    assert rows == [("PWRAP_cal", old.get("PWRAP_cal"), 1.2)]
    assert icetune_push.load_summary(target, baseline=True)["config"] == {**old, "PWRAP_cal": 1.2}


# Check selected matrix updates support additions and changed row counts without changing other fields
def test_render_json5_arrays():
    source = '// Keep this comment\n{\n  model: {\n    phase: 0.27,\n    rows: [[1, 2]], // Keep this comment too\n  }\n}\n'
    updates = {('model', 'rows'): [[1, 3], [2, 4]], ('model', 'curve'): [[1, 5, 6]]}
    rendered, changed = icetune_push.render_json5_arrays(source, updates)
    expected = json5.loads(source)
    expected['model'].update(rows=updates['model', 'rows'], curve=updates['model', 'curve'])
    assert json5.loads(rendered) == expected
    assert all(comment in rendered for comment in ('// Keep this comment', '// Keep this comment too'))
    assert changed == list(updates)
    assert icetune_push.render_json5_arrays(rendered, updates) == (rendered, [])
    cleared, _ = icetune_push.render_json5_arrays(rendered, {('model', 'curve'): []})
    assert json5.loads(cleared)['model']['curve'] == []

# Tests for encapsulated icepack steering references
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier

import pytest
from core.io import steering

ROOT = Path(__file__).resolve().parents[3]


# Concurrent dataset imports must share a fully initialized module
def test_python_module_concurrent_import(tmp_path):
    source = tmp_path / "selection.py"
    source.write_text("import time\ntime.sleep(0.1)\ncut_param = {'value': 7}\n")
    barrier = Barrier(8)

    # Start all dataset imports together and access the completed definition immediately
    def load(_):
        barrier.wait(timeout=10)
        module = steering.load_python_module(str(source))
        assert module.cut_param == {"value": 7}
        return module

    with ThreadPoolExecutor(max_workers=8) as pool:
        modules = list(pool.map(load, range(8)))
    assert all(module is modules[0] for module in modules)


# Failed imports must be retried rather than exposing a partially initialized module
def test_python_module_failed_import(tmp_path):
    source = tmp_path / "failed.py"
    source.write_text("raise RuntimeError('incomplete selection')\n")
    for _ in range(2):
        with pytest.raises(RuntimeError, match="incomplete selection"):
            steering.load_python_module(str(source))


# Compute one complete plotting block for dataset schema tests
def valid_plot(normalization="cross_section"):
    plot = {
        "normalization": normalization,
        "ratio_uncertainty": "separate",
        "stack": False,
        "data_style": "hist",
    }
    if normalization == "unit_density":
        plot["density_uncertainty"] = "shape"
    return plot


# Compute one complete fit block for dataset schema tests
def valid_fit(normalization="cross_section"):
    return {"normalization": normalization}


# Compute one complete nominal generator sample for dataset schema tests
def valid_samples():
    return [
        {
            "name": "nominal",
            "label": "GRANIITTI",
            "gencard": "./gencard.json",
            "parameters": {},
        }
    ]


# Compute one minimal valid dataset card for strict schema tests
def valid_dataset_card():
    return {
        "active": True,
        "type": "HEPDATA_TEST",
        "reader": "./reader.py",
        "datapath": "HEPData",
        "samples": valid_samples(),
        "plot": valid_plot(),
        "fit": valid_fit(),
        "sets": [
            {
                "name": "Test",
                "pid": [211, -211],
                "cuts": "./cuts.py",
                "obs": "./obs.py",
                "hist": [{"file": "table.csv", "obs": "M", "scale": 1.0}],
            }
        ],
    }


# Resolve built-in and dataset-local Python definitions deterministically
def test_python_reference_contract(tmp_path):
    project = tmp_path / "project"
    bundle = project / "icepack" / "measurement"
    local = bundle / "defs"
    local.mkdir(parents=True)
    dataset_path = bundle / "dataset.json"
    cut_path = local / "cuts.py"
    dataset_path.touch()
    cut_path.touch()

    assert (
        steering.resolve_python_reference(
            "default",
            package="core.analysis.cuts",
            dataset_path=dataset_path,
            cdir=project,
        )
        == "core.analysis.cuts.default"
    )
    assert steering.resolve_python_reference(
        "./defs/cuts.py",
        package="core.analysis.cuts",
        dataset_path=dataset_path,
        cdir=project,
    ) == str(cut_path)


# Resolve bare data files from datapath and explicit references from the dataset bundle
def test_data_reference_contract(tmp_path):
    project = tmp_path / "project"
    bundle = project / "icepack" / "measurement" / "fiducial"
    shared = project / "HEPData"
    publication = bundle.parent / "publication"
    bundle.mkdir(parents=True)
    shared.mkdir()
    publication.mkdir()
    dataset_path = bundle / "dataset.json"
    shared_table = shared / "table.csv"
    local_table = bundle / "local.csv"
    publication_table = publication / "scalar.csv"
    for path in (dataset_path, shared_table, local_table, publication_table):
        path.touch()

    assert steering.resolve_data_reference(
        "table.csv",
        datapath="HEPData",
        dataset_path=dataset_path,
        cdir=project,
    ) == str(shared_table)
    assert steering.resolve_data_reference(
        "./local.csv",
        datapath="HEPData",
        dataset_path=dataset_path,
        cdir=project,
    ) == str(local_table)
    assert steering.resolve_data_reference(
        "../publication/scalar.csv",
        datapath="HEPData",
        dataset_path=dataset_path,
        cdir=project,
    ) == str(publication_table)


# Reject ambiguous slash-containing Python references without a source suffix
def test_python_reference_rejects_ambiguous_paths(tmp_path):
    with pytest.raises(ValueError, match="must end in '.py'"):
        steering.resolve_python_reference(
            "icepack/measurement/cuts",
            package="core.analysis.cuts",
            dataset_path=None,
            cdir=tmp_path,
        )


# Load observable definitions from a dataset-relative Python file
def test_dataset_relative_observable_loading(tmp_path):
    project = tmp_path / "project"
    bundle = project / "icepack" / "measurement"
    bundle.mkdir(parents=True)
    dataset_path = bundle / "dataset.json"
    obs_path = bundle / "observables.py"
    dataset_path.touch()
    obs_path.write_text(
        "def value(event):\n"
        "    return event\n"
        "obs_value = {'tag': 'value', 'func': value, 'density': False}\n",
        encoding="utf-8",
    )

    observables, resolved = steering.load_observables(
        "./observables.py",
        dataset_path=dataset_path,
        cdir=project,
    )

    assert resolved == str(obs_path)
    assert observables["value"]["func"](3) == 3


# Resolve sample cards, finite scales and generator JSON overrides into one plan
def test_multi_sample_generation_plan_explicit(tmp_path):
    project = tmp_path / "project"
    bundle = project / "icepack" / "measurement"
    bundle.mkdir(parents=True)
    dataset_path = bundle / "dataset.json"
    base_card = bundle / "gencard.json"
    alternate_card = bundle / "alternate.json"
    for path in (dataset_path, base_card, alternate_card):
        path.touch()

    dataset = {
        "samples": [
            {
                "name": "elastic",
                "label": "Elastic",
                "gencard": "./gencard.json",
                "parameters": {"SCATTERING.NSTARS": 0},
            },
            {
                "name": "dissociative",
                "label": "Dissociative",
                "gencard": "./alternate.json",
                "parameters": {
                    "SCATTERING.NSTARS": 1,
                    "VETOCUTS.active": True,
                },
                "scale": 1.25,
            },
        ],
        "sets": [{"sample": "elastic"}, {"sample": "dissociative"}],
        "plot": {
            "normalization": "unit_density",
            "density_uncertainty": "shape",
            "ratio_uncertainty": "separate",
            "stack": True,
            "data_style": "errorbar",
        },
        "fit": valid_fit("unit_density"),
    }

    steering.validate_dataset_samples(dataset, path=str(dataset_path))
    plan = steering.build_generation_plan(
        dataset,
        dataset_path=str(dataset_path),
        cdir=str(project),
        output_prefix="measurement",
    )

    assert plan["stack"]
    assert plan["density"]
    assert [sample["output"] for sample in plan["samples"]] == [
        "measurement_elastic",
        "measurement_dissociative",
    ]
    assert plan["samples"][0]["gencard"] == str(base_card)
    assert plan["samples"][1]["gencard"] == str(alternate_card)
    assert plan["samples"][1]["scale"] == pytest.approx(1.25)
    assert [sample["assignment"] for sample in plan["samples"]] == [
        "elastic",
        "dissociative",
    ]
    assert plan["samples"][1]["overrides"] == [
        "SCATTERING.NSTARS=1",
        "VETOCUTS.active=true",
    ]


# Route several generator samples to one dataset set
def test_generation_plan_supports_sample_overlays(tmp_path):
    project = tmp_path / "project"
    bundle = project / "icepack" / "measurement"
    bundle.mkdir(parents=True)
    dataset_path = bundle / "dataset.json"
    gencard = bundle / "gencard.json"
    dataset_path.touch()
    gencard.touch()
    dataset = {
        "samples": [
            {
                "name": "impulse",
                "label": "Impulse",
                "gencard": "./gencard.json",
                "parameters": {},
            },
            {
                "name": "shadowing",
                "label": "Shadowing",
                "gencard": "./gencard.json",
                "parameters": {},
            },
        ],
        "sets": [{"samples": ["impulse", "shadowing"]}],
        "plot": valid_plot(),
        "fit": valid_fit(),
    }

    steering.validate_dataset_samples(dataset, path=str(dataset_path))
    plan = steering.build_generation_plan(
        dataset,
        dataset_path=str(dataset_path),
        cdir=str(project),
        output_prefix="measurement",
    )

    assert [sample["assignment"] for sample in plan["samples"]] == [
        "impulse",
        "shadowing",
    ]


# Preserve physics notation in a directory-derived generation output prefix
def test_generation_plan_accepts_plus_output_prefix(tmp_path):
    dataset_path = tmp_path / "dataset.json"
    gencard = tmp_path / "gencard.json"
    dataset_path.touch()
    gencard.touch()
    dataset = {
        "samples": valid_samples(),
        "sets": [],
        "plot": valid_plot(),
        "fit": valid_fit(),
    }

    plan = steering.build_generation_plan(
        dataset,
        dataset_path=str(dataset_path),
        cdir=str(tmp_path),
        output_prefix="SOFTCEP__processes__RES+CON_all_pipi",
    )

    assert plan["samples"][0]["output"] == "SOFTCEP__processes__RES+CON_all_pipi"


# Reject ambiguous sample identities and unsupported plotting keys
def test_ambiguous_samples():
    duplicate_samples = {
        "samples": [
            {"name": "same", "label": "First", "gencard": "./a.json", "parameters": {}},
            {"name": "same", "label": "Second", "gencard": "./b.json", "parameters": {}},
        ],
        "plot": valid_plot(),
        "fit": valid_fit(),
    }
    with pytest.raises(ValueError, match="duplicate sample name"):
        steering.validate_dataset_samples(duplicate_samples, path="dataset.json")

    unknown_plot_key = {
        "samples": valid_samples(),
        "plot": {**valid_plot(), "normalize": True},
        "fit": valid_fit(),
    }
    with pytest.raises(KeyError, match="unknown keys"):
        steering.validate_dataset_samples(unknown_plot_key, path="dataset.json")

    unknown_sample_key = {
        "samples": [
            {
                "name": "sample",
                "label": "Sample",
                "gencard": "./gencard.json",
                "parameters": {},
                "extra": True,
            },
        ],
        "plot": valid_plot(),
        "fit": valid_fit(),
    }
    with pytest.raises(KeyError, match="unknown keys"):
        steering.validate_dataset_samples(unknown_sample_key, path="dataset.json")

    missing_density_uncertainty = {
        "samples": valid_samples(),
        "plot": {
            "normalization": "unit_density",
            "ratio_uncertainty": "separate",
            "stack": False,
            "data_style": "hist",
        },
        "fit": valid_fit(),
    }
    with pytest.raises(KeyError, match="density_uncertainty is required"):
        steering.validate_dataset_samples(missing_density_uncertainty, path="dataset.json")

    cross_section_density_uncertainty = {
        "samples": valid_samples(),
        "plot": {**valid_plot(), "density_uncertainty": "scaled"},
        "fit": valid_fit(),
    }
    with pytest.raises(KeyError, match="only valid for unit_density"):
        steering.validate_dataset_samples(
            cross_section_density_uncertainty,
            path="dataset.json",
        )

    invalid_data_style = {
        "samples": valid_samples(),
        "plot": {**valid_plot(), "data_style": "points"},
        "fit": valid_fit(),
    }
    with pytest.raises(ValueError, match="data_style"):
        steering.validate_dataset_samples(invalid_data_style, path="dataset.json")

    invalid_unit = {
        "samples": valid_samples(),
        "plot": {**valid_plot(), "unit": "barn"},
        "fit": valid_fit(),
    }
    with pytest.raises(ValueError, match="plot.unit"):
        steering.validate_dataset_samples(invalid_unit, path="dataset.json")

    invalid_ratio_uncertainty = {
        "samples": valid_samples(),
        "plot": {**valid_plot(), "ratio_uncertainty": "denominator"},
        "fit": valid_fit(),
    }
    with pytest.raises(ValueError, match="ratio_uncertainty"):
        steering.validate_dataset_samples(invalid_ratio_uncertainty, path="dataset.json")

    invalid_ratio_plot = {
        "samples": valid_samples(),
        "plot": {**valid_plot(), "ratio_plot": "yes"},
        "fit": valid_fit(),
    }
    with pytest.raises(TypeError, match="ratio_plot"):
        steering.validate_dataset_samples(invalid_ratio_plot, path="dataset.json")

    invalid_fit = {
        "samples": valid_samples(),
        "plot": valid_plot(),
        "fit": {"normalization": "events"},
    }
    with pytest.raises(ValueError, match="fit.normalization"):
        steering.validate_dataset_samples(invalid_fit, path="dataset.json")

    partially_routed = {
        "samples": [
            {"name": "sd", "label": "SD", "gencard": "./sd.json", "parameters": {}},
            {"name": "dd", "label": "DD", "gencard": "./dd.json", "parameters": {}},
        ],
        "sets": [{"sample": "sd"}, {}],
        "plot": valid_plot(),
        "fit": valid_fit(),
    }
    with pytest.raises(ValueError, match="must all define sample"):
        steering.validate_dataset_samples(partially_routed, path="dataset.json")

    uncovered_sample = {
        "samples": [
            {"name": "sd", "label": "SD", "gencard": "./sd.json", "parameters": {}},
        ],
        "sets": [{"sample": "dd"}],
        "plot": valid_plot(),
        "fit": valid_fit(),
    }
    with pytest.raises(ValueError, match="unknown sample"):
        steering.validate_dataset_samples(uncovered_sample, path="dataset.json")

    empty_set_sample = {
        "samples": [
            {"name": "sd", "label": "SD", "gencard": "./sd.json", "parameters": {}},
        ],
        "sets": [{"sample": ""}],
        "plot": valid_plot(),
        "fit": valid_fit(),
    }
    with pytest.raises(ValueError, match="non-empty string"):
        steering.validate_dataset_samples(empty_set_sample, path="dataset.json")


# Reject uppercase and unknown keys at every strict dataset schema level
@pytest.mark.parametrize(
    ("level", "key"),
    (
        ("top", "ACTIVE"),
        ("set", "NAME"),
        ("hist", "FILE"),
        ("plot", "NORMALIZATION"),
        ("fit", "NORMALIZATION"),
        ("sample", "NAME"),
    ),
)
def test_dataset_schema_rejects_unknown_keys(tmp_path, level, key):
    card = valid_dataset_card()
    if level == "top":
        card[key] = True
    elif level == "set":
        card["sets"][0][key] = True
    elif level == "hist":
        card["sets"][0]["hist"][0][key] = True
    elif level in {"plot", "fit"}:
        card[level][key] = True
    else:
        card["samples"] = [
            {
                "name": "sample",
                "label": "Sample",
                "gencard": "./gencard.json",
                "parameters": {},
                key: True,
            }
        ]
    card_path = tmp_path / f"{level}.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    with pytest.raises(KeyError, match="unknown keys"):
        steering.load_dataset(str(card_path), cdir=tmp_path)


# Reject the removed dataset level generator card and incomplete sample blocks
def test_dataset_schema_requires_sample_gen_cards(tmp_path):
    top_level = valid_dataset_card()
    top_level["gencard"] = "./gencard.json"
    top_level_path = tmp_path / "top_level_gencard.json"
    top_level_path.write_text(json.dumps(top_level), encoding="utf-8")
    with pytest.raises(KeyError, match="unknown keys.*gencard"):
        steering.load_dataset(str(top_level_path), cdir=tmp_path)

    missing_sample_card = valid_dataset_card()
    missing_sample_card["samples"][0].pop("gencard")
    missing_sample_path = tmp_path / "missing_sample_gencard.json"
    missing_sample_path.write_text(json.dumps(missing_sample_card), encoding="utf-8")
    with pytest.raises(KeyError, match="missing keys.*gencard"):
        steering.load_dataset(str(missing_sample_path), cdir=tmp_path)


# Require lowercase scalar histogram bounds only for RAW_SCALAR datasets
def test_raw_scalar_schema_lowercase_hist_bounds(tmp_path):
    card = valid_dataset_card()
    card["type"] = "RAW_SCALAR"
    card.pop("reader")
    card["sets"][0]["hist"][0].update({"xmin": 0.0, "xmax": 1.0, "nbins": 10})
    card_path = tmp_path / "raw.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    loaded, _ = steering.load_dataset(str(card_path), cdir=tmp_path)
    assert loaded["sets"][0]["hist"][0]["nbins"] == 10

    del card["sets"][0]["hist"][0]["xmin"]
    card_path.write_text(json.dumps(card), encoding="utf-8")
    with pytest.raises(KeyError, match="RAW_SCALAR keys"):
        steering.load_dataset(str(card_path), cdir=tmp_path)


# Accept set-specific titles and positive finite MC scales
def test_dataset_schema_accepts_plot_fields(tmp_path):
    card = valid_dataset_card()
    card["sets"][0].update(
        {
            "title": "$2.0 < y < 4.5$",
            "mc_scale": 2.5,
        }
    )
    card_path = tmp_path / "valid_set_fields.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    loaded, _ = steering.load_dataset(str(card_path), cdir=tmp_path)

    assert loaded["sets"][0]["title"] == "$2.0 < y < 4.5$"
    assert loaded["sets"][0]["mc_scale"] == pytest.approx(2.5)


# Accept an MC-only set without data file or scaling fields
def test_dataset_schema_accepts_mc_only_set(tmp_path):
    card = valid_dataset_card()
    card["sets"][0]["data"] = False
    card["sets"][0]["hist"] = [{"obs": "M"}]
    card_path = tmp_path / "mc_only.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    loaded, _ = steering.load_dataset(str(card_path), cdir=tmp_path)

    assert loaded["sets"][0]["data"] is False
    assert loaded["sets"][0]["hist"] == [{"obs": "M"}]


# Accept a complete MC-only card without data steering
def test_dataset_schema_accepts_data_free_card(tmp_path):
    card = valid_dataset_card()
    card.pop("reader")
    card.pop("datapath")
    card["type"] = "MC_ONLY"
    card["sets"][0]["data"] = False
    card["sets"][0]["hist"] = [{"obs": "M"}]
    card_path = tmp_path / "data_free.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    loaded, _ = steering.load_dataset(str(card_path), cdir=tmp_path)

    assert "reader" not in loaded
    assert "datapath" not in loaded


# Accept an omitted integral limit or a unit factor
@pytest.mark.parametrize("factor", [None, 1.0])
def test_dataset_schema_optional_integral_factor(tmp_path, factor):
    card = valid_dataset_card()
    card["validation"] = {"max_integral_factor": factor}
    path = tmp_path / "validation.json"
    path.write_text(json.dumps(card), encoding="utf-8")
    loaded, _ = steering.load_dataset(str(path), cdir=tmp_path)
    assert loaded["validation"]["max_integral_factor"] == factor


# Accept generic report bounds and reject a subunit factor
def test_dataset_schema_validates_report_bounds(tmp_path):
    card = valid_dataset_card()
    card["validation"] = {
        "loopscreen": 1,
        "require_comparison": True,
        "require_differential_comparison": True,
        "require_fiducial_integral_comparison": True,
        "max_integral_factor": 1.75,
        "max_mc_rel_uncertainty": 0.25,
        "min_mc_effective_events": 100.0,
        "min_mc_ess_fraction": 0.01,
    }
    card_path = tmp_path / "validation.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    loaded, _ = steering.load_dataset(str(card_path), cdir=tmp_path)
    assert loaded["validation"]["max_integral_factor"] == pytest.approx(1.75)
    assert loaded["validation"]["max_mc_rel_uncertainty"] == pytest.approx(0.25)
    assert loaded["validation"]["min_mc_effective_events"] == pytest.approx(100.0)
    assert loaded["validation"]["min_mc_ess_fraction"] == pytest.approx(0.01)
    assert loaded["validation"]["require_differential_comparison"] is True
    assert loaded["validation"]["require_fiducial_integral_comparison"] is True

    card["validation"]["max_integral_factor"] = 0.9
    card_path.write_text(json.dumps(card), encoding="utf-8")
    with pytest.raises(ValueError, match="must be at least one"):
        steering.load_dataset(str(card_path), cdir=tmp_path)

    card["validation"]["max_integral_factor"] = 1.75
    card["validation"]["min_mc_ess_fraction"] = 1.01
    card_path.write_text(json.dumps(card), encoding="utf-8")
    with pytest.raises(ValueError, match="must not exceed one"):
        steering.load_dataset(str(card_path), cdir=tmp_path)

    card["validation"]["min_mc_ess_fraction"] = 0.01
    card["plot"]["stack"] = True
    card["samples"].append(
        {
            "name": "second",
            "label": "Second process",
            "gencard": "./gencard.json",
            "parameters": {},
        }
    )
    card_path.write_text(json.dumps(card), encoding="utf-8")
    with pytest.raises(ValueError, match="do not support process stacking"):
        steering.load_dataset(str(card_path), cdir=tmp_path)


# Require detailed comparison requests to enable the common comparison path
def test_detached_comparison(tmp_path):
    card = valid_dataset_card()
    card["validation"] = {"require_differential_comparison": True}
    card_path = tmp_path / "detached_comparison.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    with pytest.raises(ValueError, match="need validation.require_comparison"):
        steering.load_dataset(str(card_path), cdir=tmp_path)


# Accept a named MC reference for data-free, non-stacked closure tests
def test_dataset_schema_accepts_mc_reference(tmp_path):
    card = valid_dataset_card()
    card.pop("reader")
    card.pop("datapath")
    card["type"] = "MC_ONLY"
    card["samples"] = [
        {
            "name": "factorized",
            "label": "Factorized",
            "gencard": "./gencard.json",
            "parameters": {},
        },
        {"name": "central", "label": "Central", "gencard": "./gencard.json", "parameters": {}},
    ]
    card["sets"][0]["data"] = False
    card["sets"][0]["hist"] = [{"obs": "M"}]
    card["validation"] = {"mc_reference": "factorized", "require_comparison": True}
    card_path = tmp_path / "mc_reference.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")
    (tmp_path / "gencard.json").touch()

    loaded, _ = steering.load_dataset(str(card_path), cdir=tmp_path)

    assert loaded["validation"]["mc_reference"] == "factorized"


# Reject an MC reference name which is not one of the generated samples
def test_dataset_schema_unknown_mc_reference(tmp_path):
    card = valid_dataset_card()
    card["validation"] = {"mc_reference": "missing"}
    card_path = tmp_path / "invalid_mc_reference.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    with pytest.raises(ValueError, match="unknown sample"):
        steering.load_dataset(str(card_path), cdir=tmp_path)


# Reject invalid set-specific titles and MC scales
@pytest.mark.parametrize(
    ("field", "value", "error"),
    (
        ("title", "", "non-empty string"),
        ("title", "bad\ntitle", "control characters"),
        ("mc_scale", True, "finite and positive"),
        ("mc_scale", 0.0, "finite and positive"),
        ("mc_scale", -1.0, "finite and positive"),
        ("mc_scale", float("inf"), "finite and positive"),
        ("mc_scale", float("nan"), "finite and positive"),
        ("mc_scale", {}, "finite and positive"),
    ),
)
def test_dataset_schema_invalid_plot_fields(tmp_path, field, value, error):
    card = valid_dataset_card()
    card["sets"][0][field] = value
    card_path = tmp_path / "invalid_set_fields.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    with pytest.raises(ValueError, match=error):
        steering.load_dataset(str(card_path), cdir=tmp_path)


# Reject the removed reference normalization field, including positive divisors
@pytest.mark.parametrize("value", [True, 0.0, -1.0, 2.0, float("inf")])
def test_dataset_schema_rejects_data_divisor(tmp_path, value):
    card = valid_dataset_card()
    card["sets"][0]["hist"][0]["data_divisor"] = value
    card_path = tmp_path / "invalid_data_divisor.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    with pytest.raises(KeyError, match="data_divisor"):
        steering.load_dataset(str(card_path), cdir=tmp_path)


# Reject invalid selections of published table rows
@pytest.mark.parametrize("value", [[], [0, 0], [-1], [True], [1.0]])
def test_dataset_schema_rejects_invalid_rows(tmp_path, value):
    card = valid_dataset_card()
    card["sets"][0]["hist"][0]["rows"] = value
    card_path = tmp_path / "invalid_rows.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    with pytest.raises(ValueError, match="unique nonnegative integers"):
        steering.load_dataset(str(card_path), cdir=tmp_path)


# Reject a non-boolean MC-only set selector
def test_nonboolean_data_flag(tmp_path):
    card = valid_dataset_card()
    card["sets"][0]["data"] = 1
    card_path = tmp_path / "invalid_data_field.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    with pytest.raises(TypeError, match="data must be boolean"):
        steering.load_dataset(str(card_path), cdir=tmp_path)


# Reject set names that map to the same effective plot output directory
def test_duplicate_plot_outputs(tmp_path):
    card = valid_dataset_card()
    duplicate = dict(card["sets"][0])
    duplicate["name"] = "Second"
    duplicate["plotname"] = "Test"
    card["sets"].append(duplicate)
    card_path = tmp_path / "duplicate_plot_output.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    with pytest.raises(ValueError, match="duplicate plot output 'Test'"):
        steering.load_dataset(str(card_path), cdir=tmp_path)


# Accept matching observables in one shared plot group
def test_dataset_schema_accepts_shared_plot_group(tmp_path):
    card = valid_dataset_card()
    card["plot"]["data_linestyle"] = "--"
    card["plot"]["mc_linestyle"] = "-"
    card["sets"][0]["plot_group"] = "comparison"
    duplicate = dict(card["sets"][0])
    duplicate["name"] = "Second"
    card["sets"].append(duplicate)
    card_path = tmp_path / "shared_plot_group.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    dataset, _ = steering.load_dataset(str(card_path), cdir=tmp_path)

    assert [entry["plot_group"] for entry in dataset["sets"]] == [
        "comparison",
        "comparison",
    ]


# Reject mismatched observables in one shared plot group
def test_incompatible_plot_group(tmp_path):
    card = valid_dataset_card()
    card["sets"][0]["plot_group"] = "comparison"
    duplicate = dict(card["sets"][0])
    duplicate["name"] = "Second"
    duplicate["hist"] = [{"file": "table.csv", "obs": "Rap", "scale": 1.0}]
    card["sets"].append(duplicate)
    card_path = tmp_path / "incompatible_plot_group.json"
    card_path.write_text(json.dumps(card), encoding="utf-8")

    with pytest.raises(ValueError, match="same observables"):
        steering.load_dataset(str(card_path), cdir=tmp_path)


# Check plot flags precede sample fields in the shell generation plan
@pytest.mark.parametrize("loopscreen", [None, 0, 1])
@pytest.mark.parametrize("nevents", [None, 25000])
def test_generation_plan_serialize_orders_plot_flags(capsysbinary, loopscreen, nevents):
    plan = {
        "samples": [
            {
                "output": "measurement_elastic",
                "label": "Elastic",
                "gencard": "/project/gencard.json",
                "scale": 1.0,
                "assignment": "elastic",
                "fragmentation": None,
                "overrides": ["SCATTERING.NSTARS=0"],
            }
        ],
        "stack": True,
        "density": True,
        "loopscreen": loopscreen,
        "nevents": nevents,
    }

    steering.write_generation_plan(plan)

    assert capsysbinary.readouterr().out.split(b"\0") == [
        b"1",
        b"1",
        b"1",
        b"" if loopscreen is None else str(loopscreen).encode(),
        b"" if nevents is None else str(nevents).encode(),
        b"measurement_elastic",
        b"Elastic",
        b"/project/gencard.json",
        b"1.0",
        b"elastic",
        b"",
        b"",
        b"",
        b"",
        b"",
        b"1",
        b"SCATTERING.NSTARS=0",
        b"",
    ]




# Reject mistyped observable arguments before event processing
@pytest.mark.parametrize('arguments', [{'target_mass': 1.0}, {'event': None}])
def test_histogram_argument_signature(arguments):
    from core.analysis import obs

    with pytest.raises(TypeError):
        steering.select_histogram_observables(
            {'hist': [{'obs': 'mass', 'args': arguments}]},
            {'mass': {'func': obs.proj_1D_4body_M_A}},
        )

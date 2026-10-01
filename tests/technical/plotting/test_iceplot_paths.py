# Tests for iceplot path and external cut-definition handling
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import subprocess
import sys
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from core.io import steering

ROOT = Path(__file__).resolve().parents[3]

from core.io import readers
from core.plot import plot
from core.stats import hist as hist_tools


# Load the installed iceplot module for direct helper testing
def load_iceplot_module():
    from importlib import import_module

    return import_module('core.iceplot')


# Multiply the sample-wide scale by every dataset-set MC scale
def test_effective_mc_scales_combine_sample_scales():
    iceplot = load_iceplot_module()

    assert iceplot.effective_mc_scales(2.0, [1.0, 2.5]) == pytest.approx([2.0, 5.0])


# Apply publication units independently to mixed-unit observables
def test_obs_unit_overrides_standard_iceplot_scaling():
    iceplot = load_iceplot_module()
    obs = {
        "mass": {"kind": "histogram", "units": {"x": "GeV", "y": "pb"}},
        "angle": {"kind": "histogram", "units": {"x": "rad", "y": "pb"}},
        "photonuclear": {"kind": "point", "units": {"x": "GeV", "y": "mb"}},
        "suppression": {"kind": "point", "units": {"x": "GeV", "y": "1"}},
    }

    configured, scales = iceplot.change_scale(
        all_obs=obs,
        args=SimpleNamespace(unit="nb"),
        units={"angle": "pb"},
    )

    assert configured["mass"]["units"]["y"] == "nb"
    assert configured["angle"]["units"]["y"] == "pb"
    assert configured["photonuclear"]["units"]["y"] == "nb"
    assert configured["suppression"]["units"]["y"] == "1"
    assert scales == pytest.approx({"mass": 1e-3, "angle": 1.0, "photonuclear": 1e-3, "suppression": 1e-3})
    assert iceplot.apply_unit_scale(2.0, scales) == pytest.approx(
        {"mass": 2e-3, "angle": 2.0, "photonuclear": 2e-3, "suppression": 2e-3}
    )


# Preserve published values and asymmetric errors through the actual iceplot unit and density conversion
def test_publication_plots_reader_values_errors():
    iceplot = load_iceplot_module()
    for path in sorted((ROOT / "icepack").rglob("dataset.json")):
        dataset, resolved = steering.load_dataset(str(path), cdir=str(ROOT))
        if not dataset.get("active", True) or not dataset.get("datapath"):
            continue
        density = dataset["plot"]["normalization"] == "unit_density"
        density_uncertainty = dataset["plot"].get("density_uncertainty", "shape")
        for entry in dataset["sets"]:
            if not entry.get("data", True):
                continue
            observables, _ = steering.load_observables(entry["obs"], dataset_path=resolved, cdir=ROOT)
            data, observables = readers.read_hepdata(
                entry, dataset["datapath"], dataset["type"], observables,
                cdir=str(ROOT), reader=dataset["reader"], dataset_path=resolved,
            )
            originals = deepcopy(data)
            observables, scales = iceplot.change_scale(
                observables, SimpleNamespace(unit=dataset["plot"].get("unit", "pb")),
                units={hist["obs"]: hist["unit"] for hist in entry["hist"] if "unit" in hist},
            )
            plots = plot.histhepdata(data, observables, scale=scales, density=density,
                                       density_uncertainty=density_uncertainty)
            for obs, raw in originals.items():
                for sub in raw["y"] if isinstance(raw["y"], dict) else (None,):
                    key = obs if sub is None else f"{obs}_{hist_tools.bins2txt(raw['bins'][sub][0])}_{hist_tools.bins2txt(raw['bins'][sub][1])}"
                    plotted = plots[key]
                    values = np.asarray(raw["y"] if sub is None else raw["y"][sub])
                    np.testing.assert_array_equal(data[obs]["y"] if sub is None else data[obs]["y"][sub], values)
                    valid = plotted["hdata"].valid
                    divisor = np.sum(values[valid] * plotted["hdata"].binwidth[valid]) if density else 1.0
                    units = {"b": 1.0, "mb": 1e3, "ub": 1e6, "nb": 1e9, "pb": 1e12, "fb": 1e15}
                    factor = 1.0 if density else raw["scale"] * units.get(observables[obs]["units"]["y"], 1.0)
                    np.testing.assert_allclose(plotted["hdata"].counts_scaled[valid], factor * values[valid] / divisor,
                                               err_msg=f"{path}: {key}")
                    if density:
                        assert plotted["hdata"].integral() == pytest.approx(1.0)
                    if (not density or density_uncertainty == "scaled") and "y_err_up" in raw:
                        errors = np.asarray([raw[f"y_err_{side}"] if sub is None else raw[f"y_err_{side}"][sub]
                                             for side in ("down", "up")])
                        np.testing.assert_allclose(plotted["total_errs"][:, valid], factor * errors[:, valid] / divisor)


# Preserve the CDF asymmetric statistical errors when converting pb to fb
def test_asym_error_unit_conversion():
    iceplot = load_iceplot_module()
    path = ROOT / "icepack/GAMMA/integrated/CDF_2007/dataset.json"
    dataset, resolved = steering.load_dataset(str(path), cdir=str(ROOT))
    entry = dataset["sets"][0]
    observables, _ = steering.load_observables(entry["obs"], dataset_path=resolved, cdir=ROOT)
    data, observables = readers.read_hepdata(entry, dataset["datapath"], dataset["type"], observables,
                                         cdir=str(ROOT), reader=dataset["reader"], dataset_path=resolved)
    observables, scale = iceplot.change_scale(observables, SimpleNamespace(unit="fb"))
    plotted = plot.histhepdata(data, observables, scale=scale)["cross_section"]
    assert plotted["hdata"].integral() == pytest.approx(1600.0)
    np.testing.assert_allclose(plotted["stat_errs"], [[300.0], [500.0]])
    np.testing.assert_allclose(plotted["total_errs"], [[np.hypot(300.0, 300.0)], [np.hypot(500.0, 300.0)]])


# Keep mixed-unit histogram summaries explicit per observable
def test_mixed_units_hist_summary(tmp_path):
    iceplot = load_iceplot_module()
    records = {
        "mass": {
            "hdata": make_histogram([1.0], 0.1),
            "obs": {"units": {"y": "nb"}},
        },
        "angle": {
            "hdata": make_histogram([2.0], 0.1),
            "obs": {"units": {"y": "pb"}},
        },
    }

    iceplot.print_mc_data_comparison_tables(
        source=None,
        mc_sets=None,
        data_sets=[records],
        names=["publication"],
        unit="nb",
        output_dir=tmp_path,
        table_prefix="mixed",
    )

    payload = json.loads((tmp_path / "mixed_set_000.json").read_text(encoding="utf-8"))
    assert payload["metadata"]["unit"] == "mixed"
    assert {row["observable"]: row["unit"] for row in payload["rows"]} == {
        "mass": "nb",
        "angle": "pb",
    }


# Verify explicit HepMC3 files and analysis directories resolve from project paths
def test_resolve_absolute_and_project_inputs(tmp_path):
    iceplot = load_iceplot_module()
    project = tmp_path / "project"
    icepack = project / "icepack"
    external = tmp_path / "external"
    project.mkdir(parents=True)
    icepack.mkdir()
    external.mkdir()

    project_analysis = icepack / "analysis"
    project_analysis.mkdir()
    project_card = project_analysis / "dataset.json"
    absolute_mc = external / "events.hepmc3"
    absolute_analysis = external / "analysis"
    absolute_analysis.mkdir()
    absolute_card = absolute_analysis / "dataset.json"
    for path in [project_card, absolute_mc, absolute_card]:
        path.touch()

    assert iceplot.resolve_analysis_input("icepack/analysis", str(project)) == (
        str(project_analysis),
        str(project_card),
    )
    assert iceplot.resolve_hepmc3_source(str(absolute_mc), str(project)) == [str(absolute_mc)]
    assert iceplot.resolve_analysis_input(str(absolute_analysis), str(project)) == (
        str(absolute_analysis),
        str(absolute_card),
    )


# Verify explicit project paths use cdir
def test_resolve_explicit_relative_inputs(tmp_path):
    iceplot = load_iceplot_module()
    project = tmp_path / "project"
    samples = project / "samples"
    analysis = project / "cards" / "measurement"
    samples.mkdir(parents=True)
    analysis.mkdir(parents=True)
    mc_path = samples / "events.hepmc3"
    card_path = analysis / "dataset.json"
    mc_path.touch()
    card_path.touch()

    assert iceplot.resolve_hepmc3_source("samples/events.hepmc3", str(project)) == [str(mc_path)]
    assert iceplot.resolve_analysis_input("cards/measurement", str(project)) == (
        str(analysis),
        str(card_path),
    )


# Require the analysis directory to contain its dataset steering card
def test_missing_dataset_input(tmp_path):
    iceplot = load_iceplot_module()
    analysis = tmp_path / "analysis"
    analysis.mkdir()

    with pytest.raises(FileNotFoundError, match="dataset.json"):
        iceplot.resolve_analysis_input(str(analysis), str(tmp_path))


# Verify a quoted glob expands into one ordered input group
def test_resolve_hepmc3_glob_source(tmp_path):
    iceplot = load_iceplot_module()
    output = tmp_path / "output"
    output.mkdir()
    first = output / "myprocess_001.hepmc3"
    second = output / "myprocess_002.hepmc3"
    ignored = output / "other.hepmc3"
    for path in (second, ignored, first):
        path.touch()

    source = iceplot.resolve_hepmc3_source("output/myprocess_*.hepmc3", str(tmp_path))

    assert source == [str(first), str(second)]
    assert iceplot.hepmc3_source_tag("myprocess_*") == "myprocess"


# Verify an unmatched HepMC3 glob reports the resolved pattern
def test_resolve_hepmc3_glob_requires_matches(tmp_path):
    iceplot = load_iceplot_module()

    with pytest.raises(FileNotFoundError, match="matched no files"):
        iceplot.resolve_hepmc3_source("missing_*", str(tmp_path))


# Verify cut paths broadcast cleanly and external files load as modules
def test_external_cut_definition_loading(tmp_path):
    iceplot = load_iceplot_module()
    cut_file = tmp_path / "fiducial.py"
    cut_file.write_text(
        "cut_param = {'THRESHOLD': 0.2}\ndef cut_func(event):\n    return event == 'accepted'\n"
    )

    references = iceplot.resolve_cli_cuts(
        cuts=[str(cut_file)], hepmc3_count=2, base_dir=str(tmp_path)
    )
    assert references == [str(cut_file), str(cut_file)]

    module = readers.load_cut_module(references[0])
    assert module.cut_param == {"THRESHOLD": 0.2}
    assert module.cut_func("accepted")
    assert not module.cut_func("rejected")


# Verify invalid cut multiplicities fail instead of creating nested references
def test_cut_definition_count_must_match_inputs(tmp_path):
    iceplot = load_iceplot_module()
    with pytest.raises(ValueError, match="one definition per HepMC3 input"):
        iceplot.resolve_cli_cuts(cuts=["first", "second"], hepmc3_count=3, base_dir=str(tmp_path))


# Verify default tags contain basenames rather than absolute directory components
def test_input_tag_uses_only_the_filename():
    iceplot = load_iceplot_module()
    assert iceplot.input_tag("/data/run/events.hepmc3") == "events"
    assert iceplot.input_tag("/data/cards/measurement.json") == "measurement"


# Verify dataset tags retain their unique icepack bundle identity
def test_dataset_tag_uses_bundle_path(tmp_path):
    iceplot = load_iceplot_module()
    card = tmp_path / "icepack" / "SOFTCEP" / "STAR_1792394" / "pipi" / "dataset.json"
    card.parent.mkdir(parents=True)
    card.touch()

    assert iceplot.dataset_tag(str(card), str(tmp_path)) == "SOFTCEP__STAR_1792394__pipi"


# Verify the CLI stores resolved absolute inputs and creates a safe output tag
def test_parse_args_accepts_absolute_inputs(monkeypatch, tmp_path):
    iceplot = load_iceplot_module()
    mc_path = tmp_path / "events.hepmc3"
    cut_path = tmp_path / "fiducial.py"
    mc_path.touch()
    cut_path.write_text("cut_param = {}\ndef cut_func(event):\n    return True\n")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "iceplot",
            "--hepmc3",
            str(mc_path),
            "--cuts",
            str(cut_path),
            "--cores",
            "1",
        ],
    )
    args = iceplot.parse_args()

    assert args.hepmc3files == [str(mc_path)]
    assert args.cut_references == [str(cut_path)]
    assert args.output == "MC__events"


# Verify the CLI retains one MC source when a glob resolves to several files
def test_parse_args_groups_glob_matches(monkeypatch, tmp_path):
    iceplot = load_iceplot_module()
    output = tmp_path / "output"
    output.mkdir()
    files = [output / "grid_001.hepmc3", output / "grid_002.hepmc3"]
    for path in files:
        path.touch()

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "iceplot",
            "--cdir",
            str(tmp_path),
            "--hepmc3",
            "output/grid_*.hepmc3",
            "--cores",
            "1",
        ],
    )
    args = iceplot.parse_args()

    assert args.hepmc3_sources == [[str(path) for path in files]]
    assert args.hepmc3files == [str(output / "grid_*.hepmc3")]
    assert args.output == "MC__grid"


# Keep analysis output names independent of the number of MC samples
def test_parse_args_compact_analysis_output(monkeypatch, tmp_path):
    iceplot = load_iceplot_module()
    inputs = [tmp_path / f"sample_{index}.hepmc3" for index in range(3)]
    for path in inputs:
        path.touch()

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "iceplot",
            "--cdir",
            str(ROOT),
            "--analysis",
            "icepack/DURHAM/partons",
            "--hepmc3",
            *(str(path) for path in inputs),
            "--cores",
            "1",
        ],
    )

    args = iceplot.parse_args()

    assert args.output == "DURHAM__partons"


# Verify process stacking supports final density normalization
def test_parse_args_validates_process_stacking(monkeypatch, tmp_path):
    iceplot = load_iceplot_module()
    first = tmp_path / "first.hepmc3"
    second = tmp_path / "second.hepmc3"
    first.touch()
    second.touch()
    base_args = [
        "iceplot",
        "--hepmc3",
        str(first),
        str(second),
        "--stack",
        "--cores",
        "1",
    ]

    monkeypatch.setattr(sys, "argv", base_args)
    args = iceplot.parse_args()
    assert args.stack
    assert args.output.endswith("__[stack]")

    monkeypatch.setattr(sys, "argv", [*base_args, "--density"])
    args = iceplot.parse_args()
    assert args.stack
    assert args.density
    assert args.output.endswith("__[density]__[stack]")

    monkeypatch.setattr(sys, "argv", ["iceplot", "--hepmc3", str(first), "--stack"])
    with pytest.raises(ValueError, match="at least two"):
        iceplot.parse_args()


# Verify dataset plot policies become the iceplot defaults
def test_parse_args_uses_dataset_plot_policy(monkeypatch, tmp_path):
    iceplot = load_iceplot_module()
    inputs = [tmp_path / f"sample_{index}.hepmc3" for index in range(3)]
    for path in inputs:
        path.touch()

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "iceplot",
            "--cdir",
            str(ROOT),
            "--analysis",
            "icepack/GAMMA/ATLAS_1377585/mumu",
            "--hepmc3",
            *(str(path) for path in inputs),
            "--cores",
            "1",
        ],
    )
    args = iceplot.parse_args()

    assert args.density
    assert args.density_uncertainty == "shape"
    assert args.ratio_uncertainty == "combined"
    assert args.stack
    assert args.unit == "pb"
    assert args.dataset_card["fit"] == {"normalization": "unit_density"}


# Verify a dataset cross-section unit is used unless the CLI overrides it
@pytest.mark.parametrize(("option", "expected"), [(None, "mb"), ("nb", "nb")])
def test_parse_args_dataset_xs_unit(monkeypatch, option, expected):
    iceplot = load_iceplot_module()
    argv = [
        "iceplot",
        "--cdir",
        str(ROOT),
        "--analysis",
        "icepack/UPC/PHOTOPROD/ALICE_1840600/jpsi_coherent",
        "--cores",
        "1",
    ]
    if option is not None:
        argv.extend(("--unit", option))
    monkeypatch.setattr(sys, "argv", argv)

    args = iceplot.parse_args()

    assert args.unit == expected


# Verify both density CLI states override the dataset plot policy
@pytest.mark.parametrize(
    ("analysis", "option", "expected"),
    [
        ("icepack/DURHAM/partons", "--density", True),
        ("icepack/DURHAM/chic/figure8a", "--no-density", False),
    ],
)
def test_parse_args_density_overrides_dataset_policy(
    monkeypatch, tmp_path, analysis, option, expected
):
    iceplot = load_iceplot_module()
    source = tmp_path / "sample.hepmc3"
    source.touch()

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "iceplot",
            "--cdir",
            str(ROOT),
            "--analysis",
            analysis,
            "--hepmc3",
            str(source),
            option,
            "--cores",
            "1",
        ],
    )

    args = iceplot.parse_args()

    assert args.density is expected
    assert ("__[density]" in args.output) is expected


# Keep a dataset stack policy valid when plotting HEPData without MC inputs
def test_cli_analysis_stack(monkeypatch):
    iceplot = load_iceplot_module()
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "iceplot",
            "--cdir",
            str(ROOT),
            "--analysis",
            "icepack/GAMMA/ATLAS_1377585/mumu",
            "--cores",
            "1",
        ],
    )

    args = iceplot.parse_args()

    assert args.hepmc3 == []
    assert args.stack


# Verify explicit plotting policies override the dataset card defaults
def test_cli_plot_overrides(monkeypatch, tmp_path):
    iceplot = load_iceplot_module()
    inputs = [tmp_path / f"sample_{index}.hepmc3" for index in range(3)]
    for path in inputs:
        path.touch()

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "iceplot",
            "--cdir",
            str(ROOT),
            "--analysis",
            "icepack/GAMMA/ATLAS_1377585/mumu",
            "--hepmc3",
            *(str(path) for path in inputs),
            "--density-uncertainty",
            "scaled",
            "--ratio-uncertainty",
            "separate",
            "--cores",
            "1",
        ],
    )
    args = iceplot.parse_args()

    assert args.dataset_card["plot"]["density_uncertainty"] == "shape"
    assert args.dataset_card["plot"]["ratio_uncertainty"] == "combined"
    assert args.density_uncertainty == "scaled"
    assert args.ratio_uncertainty == "separate"


# Build one complete synthetic HepMC3 weight summary
def make_weight_summary(weights, xsection=10.0, overflow_count=0):
    reference = weights[0]
    is_constant = all(weight == reference for weight in weights)
    first_nonconstant = next(
        (index for index, weight in enumerate(weights) if weight != reference),
        None,
    )
    return {
        "nevents": len(weights),
        "maxevents_reached": False,
        "reference_weight": reference,
        "first_nonconstant_event": first_nonconstant,
        "first_nonconstant_weight": (
            None if first_nonconstant is None else weights[first_nonconstant]
        ),
        "min_weight": min(weights),
        "max_weight": max(weights),
        "mean_weight": sum(weights) / len(weights),
        "max_abs_deviation": max(abs(weight - reference) for weight in weights),
        "finite_weight_count": len(weights),
        "missing_weight_count": 0,
        "nonfinite_weight_count": 0,
        "maximum_weight_overflow_count": overflow_count,
        "last_xsection_pb": xsection,
        "last_xsection_pb_err": 0.1,
        "is_constant": is_constant,
        "is_nonconstant": not is_constant,
        "rel_tol": 1e-12,
        "abs_tol": 0.0,
    }


# Verify split unweighted files remain one header-normalized source
def test_unweighted_split_summary():
    iceplot = load_iceplot_module()
    combined = iceplot.combine_weight_summaries(
        [make_weight_summary([1.0, 1.0]), make_weight_summary([1.0, 1.0, 1.0])]
    )

    assert combined["nevents"] == 5
    assert combined["is_constant"]
    assert iceplot.resolve_xsmode("auto", combined)[0] == "header"


# Verify split weighted files are detected even when each split is internally constant
def test_weighted_split_summary():
    iceplot = load_iceplot_module()
    combined = iceplot.combine_weight_summaries(
        [make_weight_summary([1.0, 1.0]), make_weight_summary([2.0, 2.0])]
    )

    assert combined["nevents"] == 4
    assert combined["first_nonconstant_event"] == 2
    assert combined["is_nonconstant"]
    assert iceplot.resolve_xsmode("auto", combined)[0] == "sample"


# Verify split overflow samples retain header normalization and overflow counts
def test_overflow_split_summary():
    iceplot = load_iceplot_module()
    combined = iceplot.combine_weight_summaries(
        [
            make_weight_summary([1.0, 1.0]),
            make_weight_summary([2.5, 1.0], overflow_count=1),
        ]
    )

    assert combined["nevents"] == 4
    assert combined["maximum_weight_overflow_count"] == 1
    assert combined["is_nonconstant"]
    assert iceplot.resolve_xsmode("auto", combined)[0] == "header"


# Build one synthetic worker diagnostic record
def make_chunk_diagnostic(
    filename,
    *,
    nevents,
    attempted,
    wsum,
    wsum2,
    xsection,
    xsection_error,
    after,
):
    return {
        "filename": filename,
        "before": {
            "nevents": nevents,
            "attempted_events": attempted,
            "wsum": wsum,
            "wsum2": wsum2,
            "xsection_pb": xsection,
            "xsection_pb_err": xsection_error,
        },
        "after": [after],
    }


# Verify unweighted chunk summaries preserve one header cross section
def test_unweighted_source_summary():
    iceplot = load_iceplot_module()
    chunks = [
        make_chunk_diagnostic(
            "events.hepmc3",
            nevents=2,
            attempted=2,
            wsum=2.0,
            wsum2=2.0,
            xsection=10.0,
            xsection_error=0.5,
            after={"nevents": 1, "wsum": 1.0, "wsum2": 1.0},
        ),
        make_chunk_diagnostic(
            "events.hepmc3",
            nevents=3,
            attempted=3,
            wsum=3.0,
            wsum2=3.0,
            xsection=10.0,
            xsection_error=0.5,
            after={"nevents": 2, "wsum": 2.0, "wsum2": 2.0},
        ),
    ]

    summary = iceplot.aggregate_source_diagnostics(chunks, xsmode="header")

    assert summary["before"]["xsection_pb"] == pytest.approx(10.0)
    assert summary["before"]["xsection_pb_err"] == pytest.approx(0.5)
    assert summary["before"]["ess_fraction"] == pytest.approx(1.0)
    assert summary["after"][0]["event_acceptance"] == pytest.approx(0.6)
    assert summary["after"][0]["acceptance"] == pytest.approx(0.6)
    assert summary["after"][0]["xsection_pb"] == pytest.approx(6.0)


# Verify weighted chunk summaries combine sums and attempted events globally
def test_weighted_source_summary():
    iceplot = load_iceplot_module()
    chunks = [
        make_chunk_diagnostic(
            "grid_1.hepmc3",
            nevents=2,
            attempted=2,
            wsum=4.0,
            wsum2=10.0,
            xsection=2.0e12,
            xsection_error=0.0,
            after={"nevents": 1, "wsum": 3.0, "wsum2": 9.0},
        ),
        make_chunk_diagnostic(
            "grid_2.hepmc3",
            nevents=2,
            attempted=2,
            wsum=6.0,
            wsum2=20.0,
            xsection=3.0e12,
            xsection_error=0.0,
            after={"nevents": 1, "wsum": 2.0, "wsum2": 4.0},
        ),
    ]

    summary = iceplot.aggregate_source_diagnostics(chunks, xsmode="sample")

    assert summary["before"]["xsection_pb"] == pytest.approx(2.5e12)
    assert summary["before"]["xsection_pb_err"] == pytest.approx(np.sqrt(30.0) / 4.0 * 1e12)
    assert summary["before"]["ess_fraction"] == pytest.approx(10.0**2 / (4.0 * 30.0))
    assert summary["after"][0]["wsum2"] == pytest.approx(13.0)
    assert summary["after"][0]["event_acceptance"] == pytest.approx(0.5)
    assert summary["after"][0]["acceptance"] == pytest.approx(0.5)
    assert summary["after"][0]["xsection_pb_err"] == pytest.approx(np.sqrt(13.0) / 4.0 * 1e12)


# Verify the shared table utility writes the structured report outputs
def test_print_and_save_table(tmp_path):
    iceplot = load_iceplot_module()
    json_path, text_path = iceplot.print_and_save_table(
        title="iceplot: Test table",
        headers=["Name", "value"],
        rows=[("alpha", "1.250E+00")],
        json_rows=[{"name": "alpha", "value": 1.25}],
        output_dir=tmp_path,
        table_name="table_test",
        metadata={"kind": "test"},
    )

    assert json_path == tmp_path / "table_test.json"
    assert text_path == tmp_path / "table_test.txt"
    assert json.loads(json_path.read_text(encoding="utf-8"))["rows"] == [
        {"name": "alpha", "value": 1.25}
    ]
    assert text_path.is_file()


# Construct real histogram data with fixed differential bins and uncertainties
def make_histogram(counts, error, bins=None):
    bins = np.arange(len(counts) + 1, dtype=float) if bins is None else np.asarray(bins, dtype=float)
    return hist_tools.hobj(counts=np.asarray(counts, dtype=float), errs=np.full(len(counts), error),
                        bins=bins, cbins=0.5 * (bins[1:] + bins[:-1]))


# Select MC-only observable definitions directly from the set histogram list
def test_select_mc_only_observables():
    iceplot = load_iceplot_module()
    observable = {"tag": "M", "func": lambda event: event, "bins": np.array([0.0, 1.0])}

    selected = iceplot.select_mc_only_observables(
        dataset_set={"hist": [{"obs": "M"}]},
        all_obs={"M": observable},
    )

    assert list(selected) == ["M"]
    assert selected["M"] is not observable


# Verify unit-density tables use dimensionless integral terminology
def test_density_table_norm(tmp_path):
    iceplot = load_iceplot_module()
    histogram = make_histogram([0.4, 0.6], 0.1)

    iceplot.print_mc_data_comparison_tables(
        source="events.hepmc3",
        mc_sets=[{"M": {"hdata": histogram}}],
        data_sets=None,
        names=[""],
        unit="pb",
        normalization="unit_density",
        output_dir=tmp_path,
        table_prefix="density",
    )

    payload = json.loads((tmp_path / "density_set_000.json").read_text(encoding="utf-8"))
    assert payload["metadata"]["normalization"] == "unit_density"
    assert payload["metadata"]["unit"] == "1"
    assert "mc_integral" in payload["rows"][0]
    assert "mc_xs" not in payload["rows"][0]


# Verify a zero data cross section produces an explicit undefined ratio
def test_xs_ratio_rejects_zero_denominator():
    iceplot = load_iceplot_module()

    assert iceplot.cross_section_ratio(1.0, 0.1, 0.0, 0.2) == (None, None)


# Verify density reports retain the normalization-induced MC anticorrelations
def test_shape_density_report_integral_error_zero():
    iceplot = load_iceplot_module()
    histogram = hist_tools.hobj(
        counts=np.asarray([3.0, 7.0]),
        errs=np.sqrt(np.asarray([3.0, 7.0])),
        bins=np.asarray([0.0, 1.0, 3.0]),
        cbins=np.asarray([0.5, 2.0]),
        density=True,
        density_uncertainty="shape",
    )

    assert iceplot.histogram_integral(histogram) == pytest.approx(1.0)
    assert iceplot.histogram_integral_error(histogram) == pytest.approx(0.0, abs=1.0e-15)


# Verify stacked reports compare the normalized total rather than raw components
def test_unit_density_stack_report_normalized_total(tmp_path):
    iceplot = load_iceplot_module()
    bins = np.asarray([0.0, 1.0, 2.0])
    centers = np.asarray([0.5, 1.5])

    # Build one synthetic MC histogram record
    def record(label, counts):
        return {
            "hdata": hist_tools.hobj(
                counts=np.asarray(counts, dtype=float),
                errs=np.sqrt(np.asarray(counts, dtype=float)),
                bins=bins,
                cbins=centers,
            ),
            "hfunc": "hist",
            "label": label,
            "color": "blue",
            "style": dict(iceplot.plots.hist_style_step),
            "obs": {
                "xlabel": "$M$",
                "ylabel": "$d\\sigma/dM$",
                "units": {"x": "1", "y": "pb"},
            },
        }

    mclist = [
        [{"M": record("process A", [2.0, 1.0])}],
        [{"M": record("process B", [1.0, 1.0])}],
    ]
    data = [
        {
            "M": {
                "hdata": hist_tools.hobj(
                    counts=np.asarray([0.6, 0.4]),
                    errs=np.zeros(2),
                    bins=bins,
                    cbins=centers,
                ),
                "uncertainties": [],
            }
        }
    ]
    args = SimpleNamespace(
        stack=True,
        density=True,
        density_uncertainty="shape",
        ratio_uncertainty="combined",
        mc_hepdata_sample=None,
        hepmc3files=["a.hepmc3", "b.hepmc3"],
        hepmc3_sources=[["a.hepmc3"], ["b.hepmc3"]],
        output="validation",
        output_dir=str(tmp_path),
        dataset_file="icepack/example/dataset.json",
    )

    report = iceplot.build_validation_report(
        args=args,
        mclist=mclist,
        data=data,
        names=["set/"],
        dataset_card=None,
    )
    samples = report["sets"][0]["observables"][0]["samples"]

    assert report["schema_version"] == 1
    assert report["stack"] is True
    assert len(samples) == 1
    assert samples[0]["label"] == "MC total"
    assert len(samples[0]["components"]) == 2
    assert samples[0]["integral"] == pytest.approx(1.0)
    assert samples[0]["data_integral"] == pytest.approx(1.0)
    assert samples[0]["integral_error"] == pytest.approx(0.0, abs=1.0e-15)


# Verify a HEPData-only report retains set metadata without indexing an MC sample
@pytest.mark.parametrize("stack", [False, True])
def test_data_only_report(tmp_path, stack):
    iceplot = load_iceplot_module()
    data = [
        {
            "M": {
                "hdata": make_histogram([1.0, 2.0], 0.1),
                "uncertainties": [],
            }
        }
    ]
    args = SimpleNamespace(
        stack=stack,
        density=False,
        density_uncertainty="scaled",
        ratio_uncertainty="separate",
        mc_hepdata_sample=None,
        hepmc3files=[],
        hepmc3_sources=[],
        output="hepdata_only",
        output_dir=str(tmp_path),
        dataset_file="icepack/example/dataset.json",
    )
    dataset_card = {
        "sets": [
            {
                "name": "published",
                "title": "Fiducial selection",
                "region": "fiducial",
                "mc_scale": 2.5,
            }
        ]
    }

    report = iceplot.build_validation_report(
        args=args,
        mclist=[],
        data=data,
        names=["published/"],
        dataset_card=dataset_card,
    )

    assert report["schema_version"] == 1
    assert report["stack"] is stack
    assert report["sets"] == [
        {
            "name": "published",
            "plotname": "published",
            "title": "Fiducial selection",
            "region": "fiducial",
            "mc": True,
            "mc_scale": 2.5,
            "observables": [
                {
                    "observable": "M",
                    "kind": "histogram",
                    "samples": [],
                }
            ],
        }
    ]


# Report an MC-only set without creating data comparison fields
def test_validation_report_supports_mc_without_data(tmp_path):
    iceplot = load_iceplot_module()
    mclist = [[{"M": {"hdata": make_histogram([1.0, 2.0], 0.1), "label": "MC"}}]]
    args = SimpleNamespace(
        stack=False,
        density=False,
        density_uncertainty="scaled",
        ratio_uncertainty="separate",
        mc_hepdata_sample=None,
        hepmc3files=["events.hepmc3"],
        hepmc3_sources=[["events.hepmc3"]],
        output="mc_only",
        output_dir=str(tmp_path),
        dataset_file="icepack/example/dataset.json",
    )

    report = iceplot.build_validation_report(
        args=args,
        mclist=mclist,
        data=[None],
        names=["diagnostic/"],
        dataset_card={"sets": [{"name": "diagnostic", "data": False}]},
        source_diagnostics=[{"after": [{"nevents": 80, "ess_fraction": 0.25}]}],
    )
    sample = report["sets"][0]["observables"][0]["samples"][0]

    assert sample["integral"] == pytest.approx(3.0)
    assert sample["mc_selected_events"] == 80
    assert sample["mc_ess_fraction"] == pytest.approx(0.25)
    assert sample["mc_effective_events"] == pytest.approx(20.0)
    assert "data_integral" not in sample
    assert sample["comparison_status"] == "prediction"

    # An explicit prediction remains valid beside required measured comparisons
    report["validation"]["require_comparison"] = True
    with pytest.raises(ValueError, match="no MC-to-reference comparison"):
        iceplot.check_physics(report)
    measured = dict(sample, comparison_status="compared")
    report["sets"][0]["observables"][0]["samples"].append(measured)
    iceplot.check_physics(report)
    sample["empty"] = True
    with pytest.raises(ValueError, match="empty MC histogram"):
        iceplot.check_physics(report)


# Compare every non-reference MC output against a named MC closure reference
def test_validation_report_compares_mc_reference(tmp_path):
    iceplot = load_iceplot_module()
    mclist = [
        [{"M": {"hdata": make_histogram([1.0, 2.0], 0.1), "label": "Factorized"}}],
        [{"M": {"hdata": make_histogram([1.1, 1.9], 0.1), "label": "Central"}}],
    ]
    args = SimpleNamespace(
        stack=False,
        density=False,
        density_uncertainty="scaled",
        ratio_uncertainty="separate",
        mc_hepdata_sample=None,
        hepmc3files=["factorized.hepmc3", "central.hepmc3"],
        hepmc3_sources=[["factorized.hepmc3"], ["central.hepmc3"]],
        output="mc_closure",
        output_dir=str(tmp_path),
        dataset_file="icepack/example/dataset.json",
    )
    dataset_card = {
        "samples": [
            {"name": "factorized", "label": "Factorized"},
            {"name": "central", "label": "Central"},
        ],
        "sets": [{"name": "closure", "data": False}],
        "validation": {"mc_reference": "factorized"},
    }

    report = iceplot.build_validation_report(
        args=args,
        mclist=mclist,
        data=[None],
        names=["closure/"],
        dataset_card=dataset_card,
    )
    reference, comparison = report["sets"][0]["observables"][0]["samples"]

    assert reference["comparison_status"] == "reference"
    assert comparison["comparison_status"] == "compared"
    assert comparison["comparison_kind"] == "mc_reference"
    assert comparison["comparison_reference"] == "Factorized"
    assert comparison["shape_l1"] == pytest.approx(0.2 / 3.0)
    assert comparison["chi2"] == pytest.approx(2 * 0.1**2 / (2 * 0.1**2))


# Verify numerical literature metrics use real histogram integrals and shapes
def test_validation_report_comparison_metrics():
    iceplot = load_iceplot_module()
    mc = make_histogram([2.0, 3.0, 5.0], 0.2)
    data = make_histogram([1.0, 4.0, 5.0], 0.3)

    metrics = iceplot.comparison_metrics(mc, data)

    assert metrics["data_integral"] == 10.0
    assert metrics["chi2_ndf"] == pytest.approx(2.0 / (0.2**2 + 0.3**2) / 3)
    assert metrics["ndf"] == 3
    assert np.isclose(metrics["shape_l1"], 0.2)
    assert np.isfinite(metrics["integral_pull"])
    assert metrics["mc_max_rel_uncertainty"] == pytest.approx(0.1)
    assert metrics["mc_zero_bin_count"] == 0
    assert metrics["mc_valid_bin_count"] == 3


# Report zero MC bins separately from the finite populated-bin precision
def test_report_mc_zero_bins():
    iceplot = load_iceplot_module()
    histogram = make_histogram([4.0, 0.0, -2.0], 0.2)

    summary = iceplot.histogram_summary(histogram, allow_empty=True)

    assert summary["mc_max_rel_uncertainty"] == pytest.approx(0.1)
    assert summary["mc_zero_bin_count"] == 1
    assert summary["mc_valid_bin_count"] == 3


# Reject equal length MC and data histograms with different physical bins
def test_report_bin_mismatch():
    iceplot = load_iceplot_module()
    mc = make_histogram([2.0, 3.0], 0.2, bins=[0.0, 1.0, 2.0])
    data = make_histogram([2.0, 3.0], 0.3, bins=[0.0, 0.5, 2.0])

    with pytest.raises(ValueError, match="different bin edges"):
        iceplot.comparison_metrics(mc, data)


# Verify an unavailable chi-square remains a nonfatal report field
def test_validation_report_nonfinite_chi2_nonfatal():
    iceplot = load_iceplot_module()
    mc = make_histogram([2.0, 3.0], 0.0)
    data = make_histogram([1.0, 4.0], 0.0)

    metrics = iceplot.comparison_metrics(mc, data)

    assert metrics["chi2"] is None
    assert metrics["chi2_ndf"] is None
    assert metrics["chi2_status"] == "unavailable"
    assert metrics["shape_l1"] == pytest.approx(0.4)


# Verify report integrals and chi-square retain collective data correlations
def test_validation_report_source_cov():
    iceplot = load_iceplot_module()
    mc = make_histogram([2.0, 2.0], 0.0)
    data = make_histogram([1.0, 1.0], 1.0)
    collective = [
        {
            "name": "normalization",
            "category": "systematic",
            "correlation": "collective",
            "effect": "multiplicative",
            "up": np.ones(2),
            "down": np.ones(2),
            "shift": np.ones(2),
        }
    ]

    metrics = iceplot.comparison_metrics(
        mc,
        data,
        data_uncertainties=collective,
    )

    assert metrics["data_integral_error"] == pytest.approx(2.0)
    assert metrics["chi2"] == pytest.approx(1.0)
    assert metrics["ndf"] == 1
    assert metrics["chi2_ndf"] == pytest.approx(1.0)


# Verify paired report metrics exclude data intervals omitted from publication
def test_validation_report_common_valid_intervals():
    iceplot = load_iceplot_module()
    mc = make_histogram([1.0, 100.0], 1.0)
    data = make_histogram([1.0, 0.0], 1.0)
    data.valid = np.asarray([True, False])

    metrics = iceplot.comparison_metrics(mc, data)

    assert metrics["integral"] == pytest.approx(1.0)
    assert metrics["data_integral"] == pytest.approx(1.0)
    assert metrics["integral_pull"] == pytest.approx(0.0)
    assert metrics["shape_l1"] == pytest.approx(0.0)


# Verify one routed MC input selects only dataset sets with its sample name
def test_dataset_indices_mc_select_matching_samples():
    iceplot = load_iceplot_module()
    datasets = [
        {"sample": "sd"},
        {"sample": "dd"},
        {"sample": "dd"},
    ]

    assert list(iceplot.dataset_indices_for_mc(None, datasets)) == [0, 1, 2]
    assert iceplot.dataset_indices_for_mc("sd", datasets) == [0]
    assert iceplot.dataset_indices_for_mc("dd", datasets) == [1, 2]
    with pytest.raises(ValueError, match="no dataset sets assigned"):
        iceplot.dataset_indices_for_mc("elastic", datasets)


# Exclude explicit data-only observables from event projection and ratio panels
def test_data_only_dataset_skip_mc_routing_ratios():
    iceplot = load_iceplot_module()
    datasets = [
        {"name": "event", "sample": "mc"},
        {"name": "point", "mc": False},
    ]

    assert iceplot.dataset_indices_for_mc(None, datasets) == [0]
    assert iceplot.dataset_indices_for_mc("mc", datasets) == [0]
    assert iceplot.mc_indices_for_dataset(["mc"], datasets[1], 1) == []
    card = {"sets": datasets, "plot": {}}
    assert not iceplot.plot_group_ratio_enabled(card, [1])


# Keep data-only measurements outside an independent MC process stack
def test_process_stack_data_only():
    iceplot = load_iceplot_module()
    args = SimpleNamespace(mc_hepdata_sample=None, density=False, density_uncertainty="scaled")
    card = {"sets": [{"name": "data", "mc": False}]}
    assert iceplot.build_summed_mc_process_sets(args, [[{}], [{}]], ["data"], card) == [None]


# Analyze real events while retaining published measurements without MC projections
def test_data_only_measurements_survive_mc_analysis(tmp_path):
    from tests.technical.support.hepmc import write_muon_events

    events = tmp_path / "muons.hepmc3"
    write_muon_events(events)
    output = tmp_path / "plots"
    source = ROOT / "icepack/UPC/PHOTOPROD/ALICE_1840600/jpsi_coherent"
    card, _ = steering.load_dataset(str(source / "dataset.json"), cdir=ROOT)
    card.pop("validation")
    card["sets"][0]["mc"] = False
    for entry in card["sets"]:
        entry.pop("samples")
        for key in ("cuts", "obs"):
            entry[key] = str(source / entry[key])
    card["sets"].append(dict(name="muons", data=False, pid=[-13, 13],
                             cuts=str(ROOT / "icepack/_common/cuts_inclusive.py"),
                             obs="default", hist=[dict(obs="M")]))
    (tmp_path / "dataset.json").write_text(json.dumps(card))
    subprocess.run(
        [sys.executable, "-m", "core.iceplot", "--hepmc3", str(events),
         "--analysis", str(tmp_path),
         "--output-dir", str(output), "--output", "mixed", "--cores", "1"],
        check=True, capture_output=True, text=True,
    )
    assert list((output / "mixed").glob("**/hplot__M.pdf"))
    assert list((output / "mixed").glob("**/hplot__pt2.pdf"))


# Route several MC inputs to one dataset set
def test_overlay_sample_selection():
    iceplot = load_iceplot_module()
    assignments = ["impulse", "glauber", "shadowing", "other"]
    dataset = {"samples": ["impulse", "glauber", "shadowing"]}
    datasets = [dataset, {"sample": "other"}]

    iceplot.validate_mc_hepdata_samples(assignments, datasets, len(assignments))
    assert list(iceplot.mc_indices_for_dataset(assignments, dataset, len(assignments))) == [
        0,
        1,
        2,
    ]
    assert iceplot.dataset_indices_for_mc("shadowing", datasets) == [0]


# Group only sets that explicitly share a plot group
def test_dataset_plot_groups_order():
    iceplot = load_iceplot_module()
    datasets = [
        {"name": "first", "plot_group": "comparison"},
        {"name": "standalone", "plotname": "single"},
        {"name": "second", "plot_group": "comparison"},
    ]

    assert iceplot.dataset_plot_groups(datasets) == [
        ("comparison/", [0, 2]),
        ("single/", [1]),
    ]


# Keep direct CLI plots below the requested output directory
def test_plot_group_empty_name():
    iceplot = load_iceplot_module()
    assert iceplot.dataset_plot_groups([{"name": ""}]) == [("", [0])]


# Assign palette colors to MC and black to data
def test_append_colored_plot_channel_black_data():
    iceplot = load_iceplot_module()
    first = {"M": {"color": None}}
    second = {"M": {"color": None}}
    reference = {"M": {"color": None}}

    grouped_mc = []
    grouped_data = []
    next_index = iceplot.append_colored_plot_channel(
        grouped_mc=grouped_mc,
        grouped_data=grouped_data,
        mc_sets=[first, second],
        data_set=reference,
        color_index=0,
        data_color="black",
    )

    assert next_index == 2
    assert [record["M"]["color"] for record in grouped_mc + grouped_data] == [
        iceplot.plots.colors(0),
        iceplot.plots.colors(1),
        "black",
    ]
    assert len({record["M"]["color"] for record in grouped_mc + grouped_data}) == 3

    grouped_mc = []
    grouped_data = []
    next_index = iceplot.append_colored_plot_channel(
        grouped_mc=grouped_mc,
        grouped_data=grouped_data,
        mc_sets=[first],
        data_set=reference,
        color_index=2,
    )

    assert next_index == 4
    assert [record["M"]["color"] for record in grouped_mc + grouped_data] == [
        iceplot.plots.colors(2),
        iceplot.plots.colors(3),
    ]

    grouped_mc = []
    grouped_data = []
    next_index = iceplot.append_colored_plot_channel(
        grouped_mc, grouped_data, [first], reference, 2, data_color="black", match_colors=True
    )
    assert next_index == 3
    assert grouped_mc[0]["M"]["color"] == grouped_data[0]["M"]["color"] == iceplot.plots.colors(2)
    assert first["M"]["color"] is None
    assert reference["M"]["color"] is None


# Keep grouped ratios opt in while allowing an explicit override for every group
def test_plot_group_ratio_enabled_dataset_policy():
    iceplot = load_iceplot_module()

    assert iceplot.plot_group_ratio_enabled(None, [0]) is True
    assert iceplot.plot_group_ratio_enabled(None, [0, 1]) is False
    assert iceplot.plot_group_ratio_enabled({"plot": {"ratio_plot": True}}, [0, 1]) is True
    assert iceplot.plot_group_ratio_enabled({"plot": {"ratio_plot": False}}, [0]) is False


# Apply line styles without replacing other histogram drawing settings
def test_hist_linestyle_record_style():
    iceplot = load_iceplot_module()
    histogram_sets = [{"M": {"style": {"lw": 2, "ls": "-"}}}]

    iceplot.set_histogram_linestyle(histogram_sets, "--")

    assert histogram_sets[0]["M"]["style"] == {"lw": 2, "ls": "--"}


# Verify source diagnostics retain only the routed cut sets
def test_route_source_diagnostics_selects_matching():
    iceplot = load_iceplot_module()
    summary = {
        "before": {"nevents": 20},
        "after": [
            {"nevents": 10},
            {"nevents": 4},
            {"nevents": 3},
        ],
    }

    routed = iceplot.route_source_diagnostics(summary, [1, 2])

    assert routed["before"] == {"nevents": 20}
    assert routed["after"] == [{"nevents": 4}, {"nevents": 3}]
    assert len(summary["after"]) == 3


# Verify empty MC histograms retain absolute metrics without fake normalization
def test_validation_report_empty_mc_metrics():
    iceplot = load_iceplot_module()
    mc = make_histogram([0.0, 0.0, 0.0], 0.0)
    data = make_histogram([1.0, 4.0, 5.0], 0.3)

    metrics = iceplot.comparison_metrics(mc, data)

    assert metrics["comparison_status"] == "empty_mc"
    assert metrics["chi2_ndf"] == pytest.approx((1.0**2 + 4.0**2 + 5.0**2) / 0.3**2 / 3)
    assert metrics["shape_l1"] is None
    assert metrics["integral_pull"] < 0.0
    assert metrics["mc_max_rel_uncertainty"] is None
    assert metrics["mc_zero_bin_count"] == 3


# Verify signed bin cancellation does not mark an occupied histogram as empty
def test_report_signed_cancellation():
    iceplot = load_iceplot_module()
    mc = make_histogram([-1.0, 1.0], 0.2)
    data = make_histogram([1.0, 1.0], 0.3)

    summary = iceplot.histogram_summary(mc)
    metrics = iceplot.comparison_metrics(mc, data)
    histogram = hist_tools.hobj(
        counts=np.array([-1.0, 1.0]),
        errs=np.array([0.2, 0.2]),
        bins=np.array([0.0, 1.0, 2.0]),
        cbins=np.array([0.5, 1.5]),
    )

    assert summary["integral"] == pytest.approx(0.0)
    assert summary["empty"] is False
    assert summary["mc_max_rel_uncertainty"] == pytest.approx(0.2)
    assert summary["mc_zero_bin_count"] == 0
    assert histogram.is_empty is False
    assert metrics["comparison_status"] == "compared"
    assert metrics["shape_l1"] is None


# Verify negative cross sections cannot enter a validation report
def test_validation_report_rejects_negative_hist():
    iceplot = load_iceplot_module()
    histogram = make_histogram([-1.0, 0.0], 0.1)

    with pytest.raises(ValueError, match="invalid histogram integral"):
        iceplot.histogram_summary(histogram, allow_empty=True)


# Verify iceplot persists identical JSON and tab-separated report content
def test_write_validation_report(tmp_path):
    iceplot = load_iceplot_module()
    report = {
        "schema_version": 1,
        "output": "validation",
        "plot_directory": str(tmp_path / "plots"),
        "hepdata": "icepack/example/dataset.json",
        "normalization": "cross_section",
        "density_uncertainty": None,
        "ratio_uncertainty": "combined",
        "stack": False,
        "sets": [
            {
                "name": "published",
                "title": None,
                "region": None,
                "mc_scale": 1.0,
                "observables": [
                    {
                        "observable": "M",
                        "kind": "histogram",
                        "samples": [
                            {
                                "label": "GRANIITTI",
                                "comparison_status": "compared",
                                "empty": False,
                                "integral": 1.2,
                                "integral_error": 0.1,
                                "data_integral": 1.0,
                                "data_integral_error": 0.1,
                                "integral_pull": 1.414,
                                "chi2": 2.0,
                                "ndf": 4,
                                "chi2_ndf": 0.5,
                                "shape_l1": 0.1,
                            }
                        ],
                    }
                ],
            }
        ],
    }
    path = tmp_path / "tables" / "validation.json"

    iceplot.write_validation_report(report, path)

    assert path.is_file()
    text = path.with_suffix(".txt").read_text(encoding="utf-8")
    assert "chi2_ndf" in text
    assert "mc_max_rel_uncertainty" in text
    assert "mc_effective_events" in text
    assert "published\tM\tGRANIITTI\tcompared\tfalse" in text

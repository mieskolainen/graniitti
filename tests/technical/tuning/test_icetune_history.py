# Test the icetune history command
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import math
import runpy
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]

from core.tune import history as icetune_history


# Compute one valid compact history payload
def _history() -> dict:
    return {
        "history_schema_version": 1,
        "optimization": {"cost": "loss"},
        "parameter_space": [
            {"name": "raw_scale", "lower": 0.0, "upper": 1.0},
            {"name": "raw_phase", "lower": -3.2, "upper": 3.2},
        ],
        "trials": [
            {
                "card_config": {"schema_version": 1, "parameters": {"scale": 4.0, "phase": 1.0}, "tables": {}},
                "config": {"raw_scale": 0.75, "raw_phase": 1.0},
                "trial_id": "trial-10",
                "metrics": {"loss": 0.25},
                "search_payload": {"kind": "hebo"},
                "started_at_unix": 9.0,
                "completed_at_datetime": "2026-08-23T10:00:10Z",
                "completed_at_unix": 10.0,
                "node_id": "node-b",
            },
            {
                "card_config": {"schema_version": 1, "parameters": {"scale": 2.0, "phase": -1.0}, "tables": {}},
                "config": {"raw_scale": 0.25, "raw_phase": -1.0},
                "trial_id": "trial-2",
                "metrics": {"loss": 1.5},
                "search_payload": {"kind": "random"},
                "started_at_unix": 1.0,
                "completed_at_datetime": "2026-08-23T10:00:02Z",
                "completed_at_unix": 2.0,
                "node_id": "node-a",
            },
        ],
    }


# Write one history payload for path and schema tests
def _write_history(path: Path, payload: dict | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_history() if payload is None else payload), encoding="utf-8")
    return path


# Write one resolved campaign environment snapshot
def _write_snapshot(path: Path, **environment) -> None:
    lines = ["environment:", *(f"  {key}: {value}" for key, value in environment.items())]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# Write one minimal Ray Tune trial result stream
def _write_ray_result(
    run_dir: Path,
    trial_id: str,
    *,
    cost: float,
    timestamp: float,
) -> Path:
    path = run_dir / f"CFunc_{trial_id}" / "result.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    base = {
        "config": {"x": 0.4},
        "cost_definition": "symmetric_log_ratio_chi2_v2",
        "hostname": "ray-node",
        "trial_id": trial_id,
    }
    rows = [
        base | {"ratio2": cost + offset, "time_total_s": runtime, "timestamp": timestamp - offset}
        for offset, runtime in ((1.0, 0.5), (0.0, 2.0))
    ]
    path.write_text("\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8")
    return path


# Resolve direct files, run directories and bare run names
def test_history_inputs(tmp_path):
    history_path = _write_history(
        tmp_path / "runs" / "icetune" / "GP_DENSITY_STAR_CMS" / "history.json"
    )
    run_dir = history_path.parent

    assert icetune_history.resolve_history_path(history_path, tmp_path) == history_path
    assert icetune_history.resolve_history_path(run_dir, tmp_path) == history_path
    assert (
        icetune_history.resolve_history_path(Path("GP_DENSITY_STAR_CMS"), tmp_path)
        == history_path
    )


# Resolve a bare run name to its most recently updated nested campaign history
def test_nested_campaign_history(tmp_path):
    older = _write_history(tmp_path / "runs" / "icetune" / "GP520" / "campaigns" / ("a" * 64) / "history.json")
    latest = _write_history(tmp_path / "runs" / "icetune" / "GP520" / "campaigns" / ("b" * 64) / "history.json")
    older.touch()
    latest.touch()

    assert icetune_history.resolve_history_path(Path("GP520"), tmp_path) == latest


# Resolve a bare Ray run name to its Tune experiment directory
def test_resolve_history_path_accepts_ray_run(tmp_path):
    run_dir = tmp_path / "runs" / "icetune" / "RAY_GP"
    _write_ray_result(run_dir, "a1b2c3d4", cost=4.0, timestamp=20.0)

    assert icetune_history.resolve_history_path(Path("RAY_GP"), tmp_path) == run_dir


# Load the final finite result from each Ray trial stream
def test_load_history_reads_ray_results(tmp_path):
    run_dir = tmp_path / "runs" / "icetune" / "RAY_GP"
    _write_ray_result(run_dir, "trial-b", cost=2.0, timestamp=20.0)
    _write_ray_result(run_dir, "trial-a", cost=3.0, timestamp=10.0)
    history, cost_key, points = icetune_history.load_history(run_dir)

    expected = {"backend": "ray", "cost": "ratio2"}
    assert history["optimization"] == expected
    assert cost_key == "ratio2"
    assert [point.trial_id for point in points] == ["trial-a", "trial-b"]
    assert [point.cost for point in points] == pytest.approx([3.0, 2.0])
    assert all(point.node_id == "ray-node" for point in points)
    assert icetune_history.random_regime_end(history, points) is None


# Infer the Wasserstein metric from the definition saved by the real trial writer
def test_load_ray_wasserstein_history(tmp_path):
    from core.tune.core import optimizer_cost_definition

    path = _write_ray_result(tmp_path, "trial-w", cost=2.5, timestamp=10.0)
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    for row in rows:
        row["wasserstein"] = row.pop("ratio2")
        row["cost_definition"] = optimizer_cost_definition("wasserstein")
    path.write_text("\n".join(json.dumps(row) for row in rows) + "\n")
    _, cost_key, points = icetune_history.load_history(tmp_path)
    assert cost_key == "wasserstein"
    assert [point.cost for point in points] == pytest.approx([2.5])


# Prefer immutable Ray identity metadata over launcher snapshots
def test_load_ray_history_campaign_identity(tmp_path):
    run_dir = tmp_path / "RAY_GP"
    result_path = _write_ray_result(run_dir, "trial-a", cost=3.0, timestamp=10.0)
    _write_snapshot(
        run_dir / "campaign.resolved.yml", ALGORITHM="basic", COST="loss", RAND_TRIALS=99
    )
    identity = {
        "algorithm": "hebo",
        "optimization": {"cost": "ratio2", "rand_trials": 7},
        "physics": {"plot_brand": "CUSTOM", "simdriver": "PANDORA"},
    }
    manifest = {
        "fingerprint": icetune_history.json_fingerprint(identity),
        "identity": identity,
        "schema_version": 1,
    }
    (run_dir / "icetune_campaign.json").write_text(json.dumps(manifest), encoding="utf-8")

    history, cost_key, _ = icetune_history.load_history(result_path)

    assert cost_key == "ratio2"
    expected = {"backend": "ray", "cost": "ratio2", "optimizer": "hebo", "rand_trials": 7}
    assert history["optimization"] == expected
    assert history["simdriver"] == "PANDORA"
    assert history["plot_brand"] == "CUSTOM"


# Load finite costs and order trial identifiers by their numeric parts
def test_load_history_natural_trial_order(tmp_path):
    history_path = _write_history(tmp_path / "history.json")

    _, cost_key, points = icetune_history.load_history(history_path)

    assert cost_key == "loss"
    assert [point.trial_id for point in points] == ["trial-2", "trial-10"]
    assert [point.cost for point in points] == pytest.approx([1.5, 0.25])


# Keep equal cost rows in natural trial order when sorting by cost
def test_cost_sort_is_stable(tmp_path, capsys):
    payload = _history()
    payload["trials"][0]["metrics"]["loss"] = 1.5
    history_path = _write_history(tmp_path / "history.json", payload)

    icetune_history.main(
        [
            str(history_path),
            "--sort",
            "--cdir",
            str(tmp_path),
            "--run-name",
            "test-run",
        ]
    )

    output = capsys.readouterr().out
    assert output.index("trial-2") < output.index("trial-10")


# Sort descending costs with the minimum trial last
def test_sorted_costs_puts_minimum_last():
    points = [
        icetune_history.CostPoint("trial-10", 0.25, "last", 10.0, "node-b"),
        icetune_history.CostPoint("trial-2", 1.5, "first", 2.0, "node-a"),
        icetune_history.CostPoint("trial-3", 0.75, "middle", 6.0, "node-a"),
    ]

    ordered = icetune_history.sorted_costs(points)

    assert [point.cost for point in ordered] == pytest.approx([1.5, 0.75, 0.25])
    assert ordered[-1].trial_id == "trial-10"


# Keep decoded physical values separate from raw optimizer coordinates
def test_param_evolution_distinct_bases(tmp_path):
    history_path = _write_history(tmp_path / "history.json")
    history, _, points = icetune_history.load_history(history_path)

    _, physical_names, physical = icetune_history.parameter_evolution(history, points, physical=True)
    ordered, optimizer_names, optimizer = icetune_history.parameter_evolution(history, points, physical=False)

    assert physical_names == ["phase", "scale"]
    assert physical["scale"] == pytest.approx([2.0, 4.0])
    assert optimizer_names == ["raw_scale", "raw_phase"]
    assert optimizer["raw_scale"] == pytest.approx([0.25, 0.75])
    assert icetune_history.trial_stages(history, ordered) == [False, True]


# Select fitted physical couplings without table indices or fixed unselected rows
@pytest.mark.parametrize(
    "base, labels, rows",
    [
        ("RES|f0_test:GP:g_ls", "0,2", [[0, 2, 0.8, 0.0], [2, 2, 0.5, 0.0]]),
        ("RES|f0_test:GP:helicity", "0,0", [[0, 0, 0.8, 0.0], [1, 1, 0.5, 0.0]]),
        ("CON_GP|990:[211,211]/opposite:helicity", "0,0,0",
         [[0, 0, 0, 0.8, 0.0], [1, 1, 0, 0.5, 0.0]]),
    ],
)
def test_param_evolution_selects_named_couplings(base, labels, rows):
    from core.tune.push import flatten_card_config

    parameters = {f"{base}({labels})@MAG": 0.8, f"{base}({labels})@PHASE": 0.0}
    card = {"parameters": parameters, "tables": {base: {"rows": rows}}}
    point = icetune_history.CostPoint("trial-0", 1.0, "initial", 0.0, "node", card_config=card)
    _, names, values = icetune_history.parameter_evolution({}, [point], physical=True)

    assert names == sorted(parameters)
    assert values == {name: [value] for name, value in parameters.items()}
    assert flatten_card_config(card)[f"{base}:rows[1][0]"] == rows[1][0]


# Reject malformed schema fields and costs before printing trial data
@pytest.mark.parametrize(
    "mutate, message",
    [
        (lambda data: data.update(history_schema_version=2), "history_schema_version"),
        (lambda data: data.update(history_schema_version=1.0), "history_schema_version"),
        (lambda data: data["optimization"].pop("cost"), "optimization.cost"),
        (lambda data: data["trials"][0]["metrics"].update(loss=True), "numeric"),
        (lambda data: data["trials"][0]["metrics"].update(loss=float("inf")), "nonfinite"),
        (lambda data: data["trials"][0]["metrics"].update(loss=10**400), "nonfinite"),
        (
            lambda data: data["trials"][1].update(trial_id="trial-10"),
            "duplicate trial identifier",
        ),
        (lambda data: data["trials"][0].update(node_id=""), "node_id"),
    ],
)
def test_load_history_invalid_history(tmp_path, mutate, message):
    payload = _history()
    mutate(payload)
    history_path = _write_history(tmp_path / "history.json", payload)

    with pytest.raises(ValueError, match=message):
        icetune_history.load_history(history_path)


# Derive elapsed time from the completion text when its Unix companion is absent
def test_load_history_printed_completion_time(tmp_path):
    payload = _history()
    for trial in payload["trials"]:
        trial.pop("completed_at_unix")
    history_path = _write_history(tmp_path / "history.json", payload)

    _, _, points = icetune_history.load_history(history_path)
    ordered, elapsed, _, _, _, _ = icetune_history.cost_evolution(points)

    assert [point.trial_id for point in ordered] == ["trial-2", "trial-10"]
    assert elapsed == pytest.approx([0.0, 8.0])


# Compute cumulative time from irregular chronological completion times
def test_cost_evolution_uses_completion_time():
    points = [
        icetune_history.CostPoint("trial-3", 3.0, "later", 130.0, "node-a"),
        icetune_history.CostPoint("trial-1", 5.0, "first", 100.0, "node-a"),
        icetune_history.CostPoint("trial-2", 7.0, "middle", 106.0, "node-b"),
    ]

    ordered, elapsed, costs, average, sigma, best = icetune_history.cost_evolution(points)

    assert [point.trial_id for point in ordered] == ["trial-1", "trial-2", "trial-3"]
    assert elapsed == pytest.approx([0.0, 6.0, 30.0])
    assert costs == pytest.approx([5.0, 7.0, 3.0])
    assert average == pytest.approx([5.0, 6.0, 5.0])
    assert sigma == pytest.approx([0.0, 1.0, math.sqrt(8.0 / 3.0)])
    assert best == pytest.approx([5.0, 5.0, 3.0])


# Preserve ampfit starts through loading and classify initial evaluations as optimizer work
@pytest.mark.parametrize("starts", [2, 128])
def test_ampfit_start_metadata(tmp_path, starts):
    payload = _history()
    payload["optimization"].update(optimizer="ampfit", rand_trials=500, ampfit={"starts": starts})
    payload["trials"][0]["search_payload"] = {"kind": "ampfit", "start": 1}
    payload["trials"][1]["search_payload"] = {"kind": "initial", "start": 0}
    history, _, points = icetune_history.load_history(_write_history(tmp_path / "history.json", payload))
    assert [point.start for point in points] == [0, 1]
    assert icetune_history.trial_stages(history, points) == [True, True]
    assert icetune_history.random_regime_end(history, points) is None
    groups = icetune_history.start_trajectories(points, history)
    assert [(label, indices) for label, _, indices in groups] == [("Start 0", [0]), ("Start 1", [1])]
    assert icetune_history.start_trajectories(points[1:], history)[0][1] == groups[1][1]


# Locate the last random completion on the same time axis as plotted costs
def test_random_regime_end_random_completion_time():
    points = [
        icetune_history.CostPoint("trial-1", 5.0, "first", 100.0, "node-a", "initial", None, 90.0),
        icetune_history.CostPoint("trial-2", 4.0, "random", 115.0, "node-a", "random", None, 101.0),
        icetune_history.CostPoint("trial-3", 3.0, "adaptive", 140.0, "node-a", "hebo", None, 125.0),
    ]

    assert icetune_history.random_regime_end({}, points) == pytest.approx(15.0)


# Include random trials that finish after the first adaptive result
def test_random_regime_end_includes_pending_warmup():
    point = icetune_history.CostPoint
    points = [
        point("trial-000000", 5.0, "initial", 100.0, "node-a", "initial", None, 90.0),
        point("trial-000500", 3.0, "adaptive", 130.0, "node-b", "hebo", None, 110.0),
        point("trial-000001", 4.0, "random", 180.0, "node-c", "random", None, 95.0),
    ]
    history = {"optimization": {"optimizer": "hebo", "rand_trials": 500}}

    assert icetune_history.trial_stages(history, points) == [False, True, False]
    assert icetune_history.random_regime_end(history, points) == pytest.approx(80.0)


# Keep ICEBO Sobol proposals inside the random trial regime
def test_random_regime_end_recognizes_icebo_warmup():
    point = icetune_history.CostPoint
    points = [
        point("trial-1", 5.0, "warmup", 100.0, "node-a", "icebo", "sobol_warmup", 90.0),
        point("trial-2", 3.0, "adaptive", 130.0, "node-a", "icebo", "qlognei", 120.0),
    ]

    assert icetune_history.random_regime_end({}, points) == pytest.approx(0.0)


# Accept a valid snapshot before its first completed trial
def test_empty_history_is_printable(tmp_path, capsys):
    payload = _history()
    payload["trials"] = []
    history_path = _write_history(tmp_path / "history.json", payload)

    icetune_history.main([str(history_path), "--cdir", str(tmp_path), "--run-name", "empty-run"])

    assert capsys.readouterr().out.startswith("trial_id")
    assert (tmp_path / "figs" / "icetune" / "empty-run" / "cost_evolution.png").is_file()
    assert (tmp_path / "figs" / "icetune" / "empty-run" / "cost_evolution_logy.png").is_file()
    assert (tmp_path / "figs" / "icetune" / "empty-run" / "cost_sorted.png").is_file()
    assert (tmp_path / "figs" / "icetune" / "empty-run" / "cost_sorted_logy.png").is_file()
    assert (tmp_path / "figs" / "icetune" / "empty-run" / "parameter_evolution_physical.pdf").is_file()
    assert (tmp_path / "figs" / "icetune" / "empty-run" / "parameter_evolution_optimizer.pdf").is_file()


# Infer the run name and write the plot beside existing run figures
def test_main_writes_cost_plot_bare_run_name(tmp_path, capsys):
    _write_history(
        tmp_path / "runs" / "icetune" / "GP_DENSITY_STAR_CMS" / "history.json"
    )

    icetune_history.main(["GP_DENSITY_STAR_CMS", "--cdir", str(tmp_path)])

    output = capsys.readouterr().out
    plot_path = tmp_path / "figs" / "icetune" / "GP_DENSITY_STAR_CMS" / "cost_evolution.png"
    log_path = plot_path.with_name("cost_evolution_logy.png")
    sorted_path = plot_path.with_name("cost_sorted.png")
    sorted_log_path = plot_path.with_name("cost_sorted_logy.png")
    assert plot_path.is_file()
    assert log_path.is_file()
    assert sorted_path.is_file()
    assert sorted_log_path.is_file()
    assert str(plot_path) in output
    assert str(sorted_path) in output
    assert str(sorted_log_path) in output


# Mirror nested campaign histories below the canonical campaign figure directory
def test_main_writes_cost_plots_nested_campaign(tmp_path, capsys):
    campaign_id = "f" * 64
    _write_history(tmp_path / "runs" / "icetune" / "GP520" / "campaigns" / campaign_id / "history.json")

    icetune_history.main(["GP520", "--cdir", str(tmp_path)])

    output = capsys.readouterr().out
    plot_dir = tmp_path / "figs" / "icetune" / "GP520" / "campaigns" / campaign_id
    for name in (
        "cost_evolution.png",
        "cost_evolution_logy.png",
        "cost_sorted.png",
        "cost_sorted_logy.png",
        "parameter_evolution_physical.pdf",
        "parameter_evolution_optimizer.pdf",
    ):
        path = plot_dir / name
        assert path.is_file()
        assert str(path) in output


# Print and plot a Ray run through the same history command
def test_main_writes_cost_plot_for_ray_run(tmp_path, capsys):
    run_dir = tmp_path / "runs" / "icetune" / "RAY_GP"
    _write_ray_result(run_dir, "ray-1", cost=3.0, timestamp=10.0)
    _write_ray_result(run_dir, "ray-2", cost=2.0, timestamp=20.0)

    icetune_history.main(["RAY_GP", "--cdir", str(tmp_path)])

    output = capsys.readouterr().out
    plot_path = tmp_path / "figs" / "icetune" / "RAY_GP" / "cost_evolution.png"
    log_path = plot_path.with_name("cost_evolution_logy.png")
    sorted_path = plot_path.with_name("cost_sorted.png")
    sorted_log_path = plot_path.with_name("cost_sorted_logy.png")
    assert "ray-1" in output
    assert "ray-2" in output
    assert plot_path.is_file()
    assert log_path.is_file()
    assert sorted_path.is_file()
    assert sorted_log_path.is_file()


# Dispatch history before normal tuning argument parsing and driver creation
def test_main_dispatches_history_command(monkeypatch):
    module = runpy.run_module("core.icetune")
    received = []

    monkeypatch.setattr(sys, "argv", ["icetune", "history", "run-name", "--sort"])
    monkeypatch.setattr(icetune_history, "main", lambda argv: received.extend(argv))
    monkeypatch.setitem(
        module,
        "create_driver",
        lambda _name: pytest.fail("history command created a simulator driver"),
    )

    module["main"]()

    assert received == ["run-name", "--sort"]

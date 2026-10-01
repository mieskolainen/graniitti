# Distributed amplitude L-BFGS replay using the shared histogram objective and Torch autograd
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch
from core import resource
from core.numerics import lbfgsb
from core.stats import hist
from core.stats import objective as icecost
from core.tune.optimizers.ampfit.config import load_settings
from core.tune.search import SearchState


# Preserve physical MC errors while differentiating a source variance objective
@pytest.mark.parametrize("density", [False, True])
def test_amplitude_source_errors(density):
    from core.tune.drivers.graniitti.ampfit.amplitude import source_mc

    bins = np.arange(3)
    centers = hist.edge2centerbins(bins)
    reference = hist.hobj(np.array([2.0, 1.0]), np.array([0.3, 0.2]), bins, centers,
                            density=density, density_uncertainty="shape")
    data = hist.hobj(np.array([2.5, 1.0]), np.array([0.2, 0.1]), bins, centers,
                       density=density, density_uncertainty="shape")

    # Compute weighted bins with a parameter dependent MC variance
    def cost(value):
        weights = torch.stack((value, value, value.new_tensor(1.0)))
        counts, errors, _, _ = hist.hist(np.array([0.2, 0.4, 1.2]), bins=bins, weights=weights)
        mc = hist.hobj(counts, errors, bins, centers, density=density, density_uncertainty="shape")
        results = {"mc": [[{"M": {"hdata": mc, "source_statistics": reference}}], []],
                   "data": [[{"M": {"hdata": data}}], []]}
        converted = source_mc(results)
        assert converted["mc"][-1] == []
        fixed = converted["mc"][0][0]["M"]["hdata"]
        torch.testing.assert_close(fixed.counts_scaled, mc.counts_scaled)
        np.testing.assert_allclose(fixed.errs_scaled, reference.errs_scaled)
        assert results["mc"][0][0]["M"]["hdata"] is mc
        assert not np.allclose(mc.errs_scaled.detach().numpy(), reference.errs_scaled)
        expected = ((mc.counts_scaled - value.new_tensor(data.counts_scaled))**2 /
                    value.new_tensor(reference.errs_scaled**2 + data.errs_scaled**2)).sum()
        actual = icecost.chi2_cost(fixed, data)[0]
        torch.testing.assert_close(actual, expected)
        return actual

    value = torch.tensor(3.0, dtype=torch.float64, requires_grad=True)
    assert torch.autograd.gradcheck(cost, (value,))
    assert torch.autograd.gradgradcheck(cost, (value,))


# Preserve raw weighted MC statistics under cross section and density normalization
@pytest.mark.parametrize("density", [False, True])
def test_amplitude_histogram_statistics(density):
    from core.tune.drivers.graniitti.ampfit.ampcheck import histogram_statistics

    counts, errors, bins, centers = hist.hist(
        np.array([0.1, 0.2, 1.2]), bins=np.arange(5), weights=np.array([1.0, 3.0, 2.0]))
    mc = hist.hobj(counts, errors, bins, centers, density=density, density_uncertainty="shape")
    data = hist.hobj(np.array([4.0, 2.0, 1.0, 10.0]), np.ones(4), bins, centers,
                        valid=np.array([True, True, True, False]), density=density)
    result, = histogram_statistics([{"M": {"hdata": mc}}], [{"M": {"hdata": data, "fitw": 1.0}}])
    assert result["valid_bins"] == 3
    assert result["positive_data_bins"] == 3
    assert result["positive_data_bins_with_mc"] == 2
    assert result["empty_mc_bins"] == [2]
    assert result["nonfinite_mc_bins"] == []
    assert result["mc_error_larger_than_data_bins"] == [0, 1]
    assert result["median_mc_to_data_error"] > 1.0
    assert result["min_effective_events_occupied_bin"] == pytest.approx(1.0)
    assert result["median_effective_events_occupied_bin"] == pytest.approx((16.0 / 10.0 + 1.0) / 2.0)


# Resolve continuous amplitude coordinates with unequal physical scales
@pytest.fixture
def search_inputs():
    settings = load_settings(resource("tune/settings/ampfit.json"))
    settings.update(starts=2, maxiter=100)
    return dict(
        args=SimpleNamespace(
            algorithm="ampfit",
            cost="chi2",
            rngseed=7,
            ampfit_settings=settings,
        ),
        bounds={
            "g": {"type": "float", "lower": 0.1, "upper": 3.0},
            "phi": {"type": "float", "lower": -0.7, "upper": 0.7},
        },
        initial_points={"g": 0.8, "phi": 0.15},
        async_proposals=True,
    )


# Preserve the source start while seeding the remaining ampfit trajectories
def test_amplitude_start_seed(search_inputs):
    first = SearchState(**search_inputs).ask_many(0, 2)
    assert first == SearchState(**search_inputs).ask_many(0, 2)
    changed = copy.deepcopy(search_inputs)
    changed["args"].rngseed += 1
    other = SearchState(**changed).ask_many(0, 2)
    assert first[0] == other[0]
    assert first[1][0] != other[1][0]


# Bound relative starts for positive, negative, zero and nearly bounded parameters
@pytest.mark.parametrize("relative", [0.0, 0.1, 0.25])
def test_amplitude_start_range(search_inputs, relative):
    search_inputs["args"].ampfit_settings.update(starts=32, start_relative_range=relative)
    search_inputs["bounds"] = {
        name: {"type": "float", "lower": lower, "upper": upper}
        for name, lower, upper in [("g", 0.1, 3.0), ("negative", -3.0, -0.1), ("zero", -1.0, 1.0),
                                    ("low", 0.99, 2.0), ("high", -2.0, -0.99)]}
    initial = search_inputs["initial_points"] = dict(g=0.8, negative=-0.8, zero=0.0, low=1.0, high=-1.0)
    proposals = SearchState(**search_inputs).ask_many(0, 32)
    assert len(proposals) == (32 if relative > 0.0 else 1)
    assert proposals[0][0] == pytest.approx(initial)
    for name, value in initial.items():
        values = np.array([config[name] for config, _ in proposals])
        bounds = search_inputs["bounds"][name]
        assert np.all(values >= bounds["lower"]) and np.all(values <= bounds["upper"])
        assert np.all(np.abs(values - value) <= relative * abs(value) + 1e-14)
        if relative > 0.0 and abs(value) > 0.0:
            assert np.ptp(values) > 0.0
    np.testing.assert_allclose([config["zero"] for config, _ in proposals], 0.0, atol=1e-14)


# Preserve discrete real row phases while perturbing the active production coordinates
@pytest.mark.parametrize("initial_phase", [0.0, -np.pi])
def test_amplitude_start_phases(search_inputs, initial_phase):
    from core.tune.drivers.graniitti.driver import GraniittiDriver
    from core.tune.parameters import tools

    base = "RES|f0_980:GP:g_ls"
    phase_key = tools.raw_phase_key("RES|f0_980:GP:phi")
    initial = {tools.projective_norm_key(base + "(0,0)"): 1.2,
               base + "(2,2)@PROJECTIVE": -0.4, base + "(4,4)@PROJECTIVE": 0.7,
               phase_key: initial_phase}
    search_inputs["args"].ampfit_settings.update(starts=8, start_relative_range=0.1)
    search_inputs["initial_points"] = initial
    search_inputs["bounds"] = {
        name: {"type": "float", "lower": lower, "upper": upper}
        for name in initial
        for lower, upper in [(0.1, 3.0) if tools.is_projective_norm_key(name) else
                             (-np.pi, np.pi) if name == phase_key else
                             (-tools.PROJECTIVE_ANGLE_LIMIT, tools.PROJECTIVE_ANGLE_LIMIT)]}
    driver = GraniittiDriver()
    path = str(Path(__file__).resolve().parents[3] / "modeldata/TUNE0")
    proposals = SearchState(**search_inputs).ask_many(0, 8)
    for config, _ in proposals:
        decoded, _, groups = driver._prepare_card_params(config, path=path)
        entry = next(iter(groups[3].values()))
        phases = np.array([decoded[key] for key in entry["phase_keys"]])
        assert np.all(np.isclose(phases, 0.0) | np.isclose(phases, -np.pi))
        magnitudes = np.array([decoded[key] for key in entry["mag_keys"]])
        np.testing.assert_allclose(entry["amplitudes"], magnitudes * np.exp(1j * (phases + config[phase_key])),
                                   rtol=1e-12, atol=1e-12)
    if initial_phase < 0.0:
        assert np.ptp([config[phase_key] for config, _ in proposals]) > 0.0


# Reject missing reference values through the optimizer and the real CLI
def test_amplitude_start_requires_initial(search_inputs, monkeypatch, capsys):
    from core.icetune import parse_arguments

    search_inputs["initial_points"] = None
    with pytest.raises(ValueError, match="require initial parameter values"):
        SearchState(**search_inputs)
    monkeypatch.setattr(sys, "argv", ["icetune", "--simdriver", "graniitti", "--algorithm", "ampfit",
                                      "--no_initial_point"])
    with pytest.raises(SystemExit) as error:
        parse_arguments()
    assert error.value.code == 2
    assert "require initial parameter values" in capsys.readouterr().err


# Compute coherent intensities in four bins through the shared histogram cost
def objective(value):
    bank = value.new_tensor([0.3, 0.7, 1.0, 1.3])
    phase = value.new_tensor([0.0, 0.8, 1.7, 2.4])
    weights = (bank * value[0] * torch.exp(1j * value[1]) + torch.exp(1j * phase)).abs().square()
    target = (bank * 1.4 * torch.exp(1j * value.new_tensor(-0.2)) + torch.exp(1j * phase)).abs().square()
    counts, _, bins, centers = hist.hist(np.arange(4) + 0.5, bins=np.arange(5), weights=weights)
    mc = hist.hobj(counts, torch.zeros_like(counts), bins, centers)
    data = hist.hobj(target.detach().numpy(), np.full(4, 0.2), bins, centers)
    result = icecost.evaluate_cost_bundle(
        results={"mc": [[{"M": {"hdata": mc}}]], "data": [[{"M": {"hdata": data, "fitw": 1.0}}]]},
        selected_cost="chi2",
        cost_rho="quadratic",
        cost_avg="global-mean",
        covariance_payload=None,
        rngseed=0,
        wasserstein_cache=None,
    )
    return result["metrics"]["chi2"]


# Complete actual autograd worker evaluations and rebuild the search from JSON history
@pytest.mark.parametrize("controls", [{}, {
    "learning_rate": 0.8, "line_search_decay": 0.4, "line_search_armijo": 1e-3,
}])
def test_amplitude_search_replay(search_inputs, controls):
    search_inputs["args"].ampfit_settings.update(controls)
    records = []
    names = sorted(search_inputs["bounds"])
    for _iteration in range(180):
        search = SearchState(**search_inputs)
        search.observe(json.loads(json.dumps(records)))
        proposals = search.ask_many(len(records), 2)
        replay = SearchState(**search_inputs)
        replay.observe(copy.deepcopy(records))
        assert proposals == replay.ask_many(len(records), 2)
        if not proposals:
            assert search.amplitude.finished
            break
        search.set_pending([config for config, _ in proposals])
        assert search.ask_many(len(records), 2) == []
        for config, payload in proposals:
            values = torch.tensor([config[name] for name in names], dtype=torch.float64)
            cost, gradient = lbfgsb.value_gradient(objective, values)
            records.append(
                dict(
                    trial_id=f"trial-{len(records):06d}",
                    config=config,
                    search_payload=payload,
                    metrics={"chi2": float(cost)},
                    gradient=dict(zip(names, gradient.tolist(), strict=True)),
                )
            )
    else:
        pytest.fail("Amplitude optimizer did not finish")
    best = min(records, key=lambda record: record["metrics"]["chi2"])
    assert best["metrics"]["chi2"] < 1e-12
    np.testing.assert_allclose([best["config"][name] for name in names], [1.4, -0.2], atol=1e-6)
    assert all(item["success"] for item in search.amplitude.diagnostics)


# Reject invalid optimizer steering before scheduling worker evaluations
@pytest.mark.parametrize(("name", "value"), [
    ("learning_rate", 0.0), ("learning_rate", float("nan")),
    ("line_search_decay", 0.0), ("line_search_decay", 1.0),
    ("line_search_armijo", 1.0), ("line_search_armijo", True),
    ("start_relative_range", -0.1), ("start_relative_range", True),
    ("start_relative_range", float("nan")), ("start_relative_range", float("inf")),
    ("start_relative_range", "0.1"), ("covariance", 1), ("covariance", "false"),
    ("covariance_interval", -1.0), ("covariance_interval", True),
    ("covariance_interval", float("nan")), ("covariance_interval", float("inf")),
])
def test_amplitude_steering(tmp_path, name, value):
    settings = load_settings(resource("tune/settings/ampfit.json"))
    settings[name] = value
    card = tmp_path / "ampfit.json"
    card.write_text(json.dumps(settings))
    with pytest.raises(ValueError, match=name):
        load_settings(card)


# Refuse scalar-only history instead of silently substituting approximate gradients
def test_amplitude_search_requires_gradient(search_inputs):
    search = SearchState(**search_inputs)
    config, payload = search.ask_many(0, 1)[0]
    search.observe([dict(trial_id="trial-000000", config=config, search_payload=payload, metrics={"chi2": 1.0})])
    with pytest.raises(ValueError, match="autograd gradient"):
        search.ask_many(1, 1)


# Drive actual asynchronous Ray search proposals to convergence and resume a saved search
def test_amplitude_ray_search_resume(search_inputs, tmp_path):
    import time

    from core.tune.backends.ray import PercentileSearch
    from ray.tune.search import Searcher

    args = search_inputs["args"]
    args.max_concurrent_trials = 2
    args.num_trials = 300
    args.rand_trials = 0
    args.proposal_batch_size = 2
    args.proposal_batch_fraction = 1.0
    args.surrogate_fit_percentile = 0.8
    arguments = {key: value for key, value in search_inputs.items() if key != "async_proposals"}
    search = PercentileSearch(**arguments)
    names = sorted(arguments["bounds"])
    checkpoint = str(tmp_path / "amplitude_search.pkl")
    deadline = time.monotonic() + 60
    resumed = False
    while time.monotonic() < deadline:
        trial = f"worker-{search.next_index:06d}"
        config = search.suggest(trial)
        if config == Searcher.FINISHED:
            break
        if config is None:
            time.sleep(0.01)
            continue
        values = torch.tensor([config[name] for name in names], dtype=torch.float64)
        cost, gradient = lbfgsb.value_gradient(objective, values)
        result = {"chi2": float(cost), "gradient": dict(zip(names, gradient.tolist(), strict=True))}
        search.on_trial_complete(trial, result=result)
        if len(search.completed_records) >= 4 and not resumed:
            search._collect_fit()
            search.save(checkpoint)
            if search.fit_executor is not None:
                search.fit_executor.shutdown(wait=True)
            search = PercentileSearch(**arguments)
            search.restore(checkpoint)
            resumed = True
    else:
        pytest.fail("Ray amplitude search did not finish")
    if search.fit_executor is not None:
        search.fit_executor.shutdown(wait=True)
    assert resumed
    assert search.optimizer_finished
    assert search.next_index < args.num_trials
    assert min(record["metrics"]["chi2"] for record in search.completed_records) < 1e-12


# Keep Cartesian and signed projective production derivatives at vanishing components
@pytest.mark.parametrize("geometry", ["cartesian", "projective"])
def test_driver_amp_coords_autograd(geometry):

    from core.tune.drivers.graniitti.driver import GraniittiDriver
    from core.tune.parameters import tools

    driver = GraniittiDriver()
    path = str(Path(__file__).resolve().parents[3] / "modeldata" / "TUNE0")
    values = torch.tensor([0.0, 0.0] if geometry == "cartesian" else [1.2, 0.0, 0.0, 0.3],
                          dtype=torch.float64, requires_grad=True)

    # Evaluate the same prepared coordinate groups used when writing physical steering cards
    def evaluate(theta):
        if geometry == "cartesian":
            base = "RES|f0_980:MP:g"
            config = {tools.complex_re_key(base): theta[0], tools.complex_im_key(base): theta[1]}
            _, _, groups = driver._prepare_card_params(config, path=path)
            amplitudes = [group.amplitude for group in groups[2]]
        else:
            base = "RES|f0_980:GP:g_ls"
            config = {tools.projective_norm_key(base + "(0,0)"): theta[0],
                      base + "(2,2)@PROJECTIVE": theta[1],
                      base + "(4,4)@PROJECTIVE": theta[2],
                      tools.raw_phase_key("RES|f0_980:GP:phi"): theta[3]}
            _, _, groups = driver._prepare_card_params(config, path=path)
            amplitudes = next(iter(groups[3].values()))["amplitudes"]
        return torch.view_as_real(torch.stack(amplitudes))

    assert torch.autograd.gradcheck(evaluate, (values,), eps=1e-5)
    assert torch.autograd.gradgradcheck(evaluate, (values,), eps=1e-5)


# Serialize vanishing LS norms with canonical zero phases while retaining amplitude derivatives
@pytest.mark.parametrize("geometry", ["PROJECTIVE", "SPHERICAL"])
def test_zero_ls_norm_coordinates(geometry):

    from core.tune.drivers.graniitti.driver import GraniittiDriver
    from core.tune.parameters import tools

    driver = GraniittiDriver()
    path = str(Path(__file__).resolve().parents[3] / "modeldata" / "TUNE0")
    base = "RES|f0_980:GP:g_ls"
    norm = torch.tensor(0.0, dtype=torch.float64, requires_grad=True)
    config = {tools.projective_norm_key(base + "(0,0)"): norm,
              base + f"(2,2)@{geometry}": -0.4,
              base + f"(4,4)@{geometry}": -0.7}
    decoded, _, groups = driver._prepare_card_params(config, path=path)
    entry = next(iter(groups[3].values()))
    assert all(float(decoded[key]) == 0.0 for key in entry["mag_keys"] + entry["phase_keys"])
    amplitudes = torch.stack(entry["amplitudes"])
    derivative = torch.autograd.grad(amplitudes.real.sum(), norm)[0]
    assert float(derivative) == pytest.approx(sum(float(x) for x in entry["decoded_vector"]))


# Check complex interference curvature and propagate the full covariance into couplings
def test_amplitude_covariance(tmp_path):
    from core.stats.transform import ParameterTransform
    from core.tune.optimizers.ampfit.report import save_covariance

    values = torch.tensor([1.4, -0.2], dtype=torch.float64)
    names = ["magnitude", "phase"]

    # Decode one complex coupling into its real and imaginary components
    def decode(parameters):
        amplitude = parameters["magnitude"] * torch.exp(1j * parameters["phase"])
        return {"real": amplitude.real, "imaginary": amplitude.imag}

    transform = ParameterTransform(names, decode, values)
    result = save_covariance(objective=objective, values=values, names=names,
                             bounds=np.array([[0.1, 3.0], [-np.pi, np.pi]]), transform=transform,
                             output_root=tmp_path, metadata={"trial_id": "coherent"})
    bank, phases = np.array([0.3, 0.7, 1.0, 1.3]), np.array([0.0, 0.8, 1.7, 2.4])
    magnitude, phase = values.numpy()
    jacobian = np.column_stack((2 * bank**2 * magnitude + 2 * bank * np.cos(phase - phases),
                                -2 * bank * magnitude * np.sin(phase - phases)))
    expected = np.linalg.inv(jacobian.T @ jacobian / 0.2**2)
    optimizer = json.loads((tmp_path / "optimizer/parameters.json").read_text())
    physical = json.loads((tmp_path / "physical/parameters.json").read_text())
    np.testing.assert_allclose(optimizer["covariance"], expected, rtol=1e-12)
    jacobian = np.array([[np.cos(phase), -magnitude * np.sin(phase)],
                         [np.sin(phase), magnitude * np.cos(phase)]])
    np.testing.assert_allclose(physical["covariance"], jacobian @ expected @ jacobian.T, rtol=1e-12)
    np.testing.assert_allclose(np.square(physical["errors"]), np.diag(physical["covariance"]), rtol=1e-12)
    assert result["status"] == "local"
    assert optimizer["trial_id"] == physical["trial_id"] == "coherent"


# Keep bound and unresolved curvature uncertainties unavailable instead of regularized errors
@pytest.mark.parametrize("curvature", [0.0, -1.0, 1.0])
def test_amplitude_covariance_limits(tmp_path, curvature):
    from core.stats.transform import ParameterTransform
    from core.tune.optimizers.ampfit.report import save_covariance

    values = torch.zeros(2, dtype=torch.float64)
    names = ["bound", "free"]

    # Evaluate a quadratic with controlled free curvature and one active bound
    def quadratic(point):
        return point[0].square() + curvature * point[1].square()

    result = save_covariance(objective=quadratic, values=values, names=names,
                             bounds=np.array([[0.0, 1.0], [-1.0, 1.0]]),
                             transform=ParameterTransform(names, dict, values),
                             output_root=tmp_path, metadata={"trial_id": "curvature"})
    payload = json.loads((tmp_path / "optimizer/parameters.json").read_text())
    assert payload["errors"][0] is None
    if curvature > 0.0:
        assert payload["errors"][1] == pytest.approx(1.0)
        assert result["status"] == "local"
    else:
        assert payload["errors"][1] is None
        assert result["status"] == "unavailable"


# Preserve complete covariance transfers when a better fit is published before or after them
@pytest.mark.parametrize("late", [False, True])
def test_covariance_publication(tmp_path, late):
    from core.io.serialize import write_json_file
    from core.tune.backends.ray import BestState, TrialOutputCallback, pack_ray_figure_tree

    state = BestState(experiment_dir=str(tmp_path / "run"), cdir=str(tmp_path), run_name="fit", cost="chi2")
    report = tmp_path / "report"
    write_json_file(report / "covariance.json", {"trial_id": "first", "status": "local"})
    for basis in ("optimizer", "physical"):
        write_json_file(report / basis / "parameters.json",
                        {"trial_id": "first", "covariance": [[1.0, 0.2], [0.2, 2.0]]})
    expected = {path.relative_to(report): path.read_bytes() for path in report.rglob("*") if path.is_file()}
    payload = pack_ray_figure_tree(report)
    figures = tmp_path / "histograms"
    for index, trial_id in enumerate(("first", "better", "best")):
        cost = float(3 - index)
        summary = {"trial_id": trial_id, "metrics": {"chi2": cost}}
        write_json_file(figures / "summary.json", summary)
        state.publish_render({"trial_id": trial_id, "cost": cost, "render_initial": False, "render_best": True},
                             {"figures": pack_ray_figure_tree(figures), "pickle": None})
        if index == int(late):
            state.publish_covariance({"trial_id": "first"}, payload)
        assert json.loads((state.figure_dir / "summary.json").read_text()) == summary
    target = state.figure_dir / "covariance"
    for name, data in expected.items():
        assert (target / name).read_bytes() == data
    restored = TrialOutputCallback(experiment_dir=str(tmp_path / "run"), param=state.param, global_state=state)
    assert restored.rendered_trial_id == "best"
    assert restored.covariance_trial_id == "first"
    with pytest.raises(RuntimeError, match="wrong trial identity"):
        state.publish_covariance({"trial_id": "other"}, payload)
    for name, data in expected.items():
        assert (target / name).read_bytes() == data

    # Replace the complete report when curvature fails without retaining previous uncertainties
    write_json_file(tmp_path / "failed" / "covariance.json", {"trial_id": "best", "status": "unavailable"})
    state.publish_covariance({"trial_id": "best"}, pack_ray_figure_tree(tmp_path / "failed"))
    assert json.loads((target / "covariance.json").read_text())["status"] == "unavailable"
    assert not (target / "optimizer").exists()
    assert not (target / "physical").exists()

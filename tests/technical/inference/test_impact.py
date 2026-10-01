# Tests for surrogate parameter impact diagnostics
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path

import numpy as np
import pytest
import torch
from core.inference import impact
from core.plot import plot
from core.tune import likelihood as icetune_likelihood
from sklearn.linear_model import Ridge


# Build one internally consistent synthetic histogram likelihood
def _likelihood(records):
    valid_records = [record for record in records if record["valid"]]
    total = float(sum(record["chi2"] for record in valid_records))
    return {
        "schema_version": 1,
        "valid": True,
        "kind": "gaussian_chi2",
        "covariance_mode": "diagonal",
        "objective": {"name": "chi2", "value": total, "ndf": 4.0 * len(valid_records)},
        "observables": records,
        "mc_uncertainty": {
            "source": "single_trial",
            "sigma_logL": None,
            "replicate_count": 1,
        },
    }


# Build one compact observable record
def _record(dataset, subset, observable, chi2, valid=True):
    return {
        "dataset": dataset,
        "subset": subset,
        "observable": observable,
        "chi2": chi2,
        "ndf": 4.0,
        "fitw": 1.0,
        "valid": valid,
    }


# Check history loading aligns reordered and missing observable records
def test_history_loader_aligned_obs_chi2(tmp_path):
    history_path = tmp_path / "history.json"
    history = {
        "history_schema_version": 1,
        "simdriver": "PANDORA",
        "plot_brand": "CUSTOM",
        "likelihood_default": "logL",
        "parameter_space": [
            {
                "name": "x",
                "lower": 0.0,
                "upper": 1.0,
                "prior": {"type": "uniform", "lower": 0.0, "upper": 1.0},
            },
            {
                "name": "y",
                "lower": -1.0,
                "upper": 1.0,
                "prior": {"type": "uniform", "lower": -1.0, "upper": 1.0},
            },
        ],
        "trials": [
            {"theta": [0.0], "likelihood": {"covariance_mode": "full"}},
            {
                "theta": [0.1, -0.5],
                "theta_hash": "a",
                "likelihood": _likelihood(
                    [
                        _record(1, 0, "M", 7.0),
                        _record(0, 0, "Rap", 3.0),
                    ]
                ),
            },
            {
                "theta": [0.4, 0.0],
                "theta_hash": "b",
                "likelihood": _likelihood(
                    [
                        _record(0, 0, "Rap", 2.0),
                        _record(1, 0, "M", 6.0),
                    ]
                ),
            },
            {
                "theta": [0.8, 0.5],
                "theta_hash": "c",
                "likelihood": _likelihood(
                    [
                        _record(0, 0, "Rap", 5.0),
                        _record(1, 0, "M", 9.0, valid=False),
                    ]
                ),
            },
        ],
    }
    history_path.write_text(json.dumps(history), encoding="utf-8")

    loaded = icetune_likelihood.load_likelihood_history(history_path)

    assert loaded["simdriver"] == "PANDORA"
    assert loaded["plot_brand"] == "CUSTOM"
    assert loaded["covariance_mode"] == "diagonal"
    assert [
        (item["dataset"], item["subset"], item["observable"]) for item in loaded["observables"]
    ] == [(0, 0, "Rap"), (1, 0, "M")]
    assert np.allclose(
        loaded["observable_chi2"][:2],
        [[3.0, 7.0], [2.0, 6.0]],
    )
    assert loaded["observable_chi2"][2, 0] == pytest.approx(5.0)
    assert np.isnan(loaded["observable_chi2"][2, 1])


# Check smooth impact fits retain signed correlated observable responses
@pytest.mark.parametrize("offset", [10.0, -10.0])
def test_obs_surrogate_param_impacts_directional(offset):
    rng = np.random.default_rng(42)
    X = rng.uniform(0.0, 1.0, size=(240, 2))
    Y = np.column_stack(
        (
            offset + 4.0 * X[:, 0] - X[:, 1],
            3.0 - 0.5 * X[:, 0] + 2.0 * X[:, 1],
        )
    )
    model = impact.fit_observable_surrogate(
        X_model=X,
        observable_chi2=Y,
        n_features=24,
        scale=0.5,
        alpha=1.0e-4,
        validation_fraction=0.2,
        rngseed=11,
        min_trials=30,
    )
    shifts = impact.fixed_parameter_shifts(
        best_fit=np.array([0.5, 0.5]),
        errors=np.array([0.1, 0.2]),
        bounds=np.array([[0.0, 1.0], [0.0, 1.0]]),
    )
    impacts = impact.compute_parameter_impacts(
        model=model,
        best_fit=np.array([0.5, 0.5]),
        bounds=np.array([[0.0, 1.0], [0.0, 1.0]]),
        feature_transform=lambda values: values,
        shifts=shifts,
    )

    assert impacts["mode"] == "fixed"
    assert np.all(model["validation_r2"] > 0.99)
    assert model["target_transform"] == "standard"
    assert (impacts["baseline"][0] > 0.0) == (offset > 0.0)
    nearest = np.argmin(np.linalg.norm(X - np.array([0.5, 0.5]), axis=1))
    assert impacts["baseline_reference"]["kind"] == "nearest_realized_trial"
    assert np.allclose(impacts["baseline"], Y[nearest])
    assert impacts["delta_minus"][0, 0] == pytest.approx(-0.4, abs=0.03)
    assert impacts["delta_plus"][0, 0] == pytest.approx(0.4, abs=0.03)
    assert impacts["delta_minus"][0, 1] == pytest.approx(0.2, abs=0.03)
    assert impacts["delta_plus"][0, 1] == pytest.approx(-0.2, abs=0.03)
    assert impact.ranked_parameter_indices(impacts, 0, 2).tolist()[0] == 0


# Check constrained impacts follow a correlated global likelihood valley
def test_profile_correlated_shifts():
    best_fit = np.array([0.5, 0.5])
    errors = np.array([0.2, 0.2])
    bounds = np.array([[0.0, 1.0], [0.0, 1.0]])

    # Evaluate a correlated quadratic global likelihood
    def global_surrogate(values):
        values = np.asarray(values, dtype=np.float64)
        x = values[:, 0]
        y = values[:, 1]
        return 100.0 * (x + y - 1.0) ** 2 + (x - y) ** 2

    gradient_calls = []

    # Compute the exact physical global-likelihood gradient
    def value_gradient(values):
        gradient_calls.append(1)
        x, y = np.asarray(values, dtype=np.float64)
        common = 200.0 * (x + y - 1.0)
        return (
            float(global_surrogate(np.asarray([[x, y]]))[0]),
            np.asarray([common + 2.0 * (x - y), common - 2.0 * (x - y)]),
        )

    fixed = impact.fixed_parameter_shifts(
        best_fit=best_fit,
        errors=errors,
        bounds=bounds,
    )
    profiled = impact.profiled_parameter_shifts(
        func_hat=global_surrogate,
        best_fit=best_fit,
        errors=errors,
        bounds=bounds,
        covariance=np.eye(2) * errors[0] ** 2,
        maxiter=50,
        value_gradient_func=value_gradient,
    )

    assert profiled["mode"] == "profiled"
    assert profiled["plus_points"][0, 0] == pytest.approx(0.7)
    assert profiled["plus_points"][0, 1] == pytest.approx(0.3039604, abs=1.0e-5)
    assert profiled["minus_points"][0, 1] == pytest.approx(0.6960396, abs=1.0e-5)
    fixed_plus_z = global_surrogate(fixed["plus_points"][[0]])[0]
    profile_plus = profiled["profile"]["plus"][0]
    assert profile_plus["success"]
    assert gradient_calls
    assert profile_plus["value"] < fixed_plus_z
    assert profile_plus["delta_z"] == pytest.approx(profile_plus["value"])


# Check impact layout uses rendered label widths and signed shift colors
def test_plot_compacts_left_margin_colors_shifts(tmp_path, monkeypatch):
    shifts = impact.fixed_parameter_shifts(
        best_fit=np.array([0.5, 0.5]),
        errors=np.array([0.1, 0.1]),
        bounds=np.array([[0.0, 1.0], [0.0, 1.0]]),
    )

    # Compute a non-zero directional response for both parameters
    def evaluate(points):
        return 4.0 + 2.0 * np.asarray(points)[:, 0:1] - np.asarray(points)[:, 1:2]

    impacts = impact.compute_direct_parameter_impacts(
        evaluate_observable_chi2=evaluate,
        best_fit=np.array([0.5, 0.5]),
        bounds=np.array([[0.0, 1.0], [0.0, 1.0]]),
        shifts=shifts,
    )
    original_close = impact.plt.close
    monkeypatch.setattr(impact.plt, "close", lambda _fig: None)
    impact.plot_observable_impact(
        output_dir=tmp_path,
        observable={
            "dataset": 0,
            "subset": 0,
            "observable": "dPhi_pp",
            "ndf": 4,
        },
        observable_index=0,
        param_names=[
            "RES|f2_1950:GP:g_ls[0]@RE",
            "CON_GP|990:[211,211]/opposite:helicity[0]@RE",
        ],
        impacts=impacts,
        model={
            "validation_r2": np.array([0.75]),
            "validation_mae": np.array([0.1]),
            "trial_count": np.array([100]),
        },
        top_n=2,
        plot_brand="PANDORA",
    )
    figure = impact.plt.gcf()
    figure.canvas.draw()
    labels = figure.axes[0].get_yticklabels()
    minimum_x = (
        min(label.get_window_extent(figure.canvas.get_renderer()).x0 for label in labels)
        / figure.bbox.width
    )
    assert 0.005 <= minimum_x <= 0.03
    minus_bars, plus_bars = figure.axes[1].containers
    assert minus_bars.patches[0].get_facecolor()[:3] == pytest.approx(plot.colors(1))
    assert plus_bars.patches[0].get_facecolor()[:3] == pytest.approx(plot.colors(0))
    original_close(figure)


# Check rendering publishes CMS-style plots and machine-readable impacts
def test_impact_outputs(tmp_path):
    rng = np.random.default_rng(7)
    X = rng.uniform(0.0, 1.0, size=(160, 2))
    Y = np.column_stack(
        (
            5.0 + 2.0 * X[:, 0] + X[:, 1],
            8.0 - X[:, 0] + 3.0 * X[:, 1],
        )
    )
    observables = [
        {"dataset": 0, "subset": 0, "observable": "M", "ndf": 8.0, "fitw": 1.0},
        {"dataset": 0, "subset": 0, "observable": "Rap", "ndf": 6.0, "fitw": 1.0},
    ]

    result = impact.render_parameter_impacts(
        X_model=X,
        observable_chi2=Y,
        observables=observables,
        param_names=["coupling", "form_factor"],
        best_fit=np.array([0.5, 0.5]),
        errors=np.array([0.1, 0.1]),
        bounds=np.array([[0.0, 1.0], [0.0, 1.0]]),
        feature_transform=lambda values: values,
        output_dir=tmp_path,
        n_features=16,
        scale=0.5,
        alpha=1.0e-4,
        validation_fraction=0.2,
        rngseed=9,
        min_trials=30,
        top_n=2,
    )

    assert result["published_observable_count"] == 2
    assert result["mode"] == "fixed"
    assert result["surrogate"]["target_transform"] == "standard"
    assert result["surrogate"]["impact_response"] == "surrogate_differences"
    assert result["profiled_shift_count"] == 0
    assert Path(result["manifest"]).is_file()
    for output in result["outputs"]:
        assert Path(output["png"]).is_file()
        assert Path(output["pdf"]).is_file()
        payload = json.loads(Path(output["json"]).read_text(encoding="utf-8"))
        assert payload["impact_schema_version"] == 1
        assert payload["mode"] == "fixed"
        assert len(payload["parameters"]) == 2
        assert any(
            abs(value) > 1.0e-8
            for parameter in payload["parameters"]
            for value in (
                parameter["delta_chi2_minus"],
                parameter["delta_chi2_plus"],
            )
            if value is not None
        )

    # Evaluate a smooth global objective for profiled output serialization
    def global_surrogate(values):
        values = np.asarray(values, dtype=np.float64)
        return np.sum((values - 0.5) ** 2, axis=1)

    profiled_shifts = impact.profiled_parameter_shifts(
        func_hat=global_surrogate,
        best_fit=np.array([0.5, 0.5]),
        errors=np.array([0.1, 0.1]),
        bounds=np.array([[0.0, 1.0], [0.0, 1.0]]),
        covariance=np.eye(2) * 0.01,
        maxiter=20,
    )
    profiled_result = impact.render_parameter_impacts(
        X_model=X,
        observable_chi2=Y[:, :1],
        observables=observables[:1],
        param_names=["coupling", "form_factor"],
        best_fit=np.array([0.5, 0.5]),
        errors=np.array([0.1, 0.1]),
        bounds=np.array([[0.0, 1.0], [0.0, 1.0]]),
        feature_transform=lambda values: values,
        output_dir=tmp_path / "profiled",
        shifts=profiled_shifts,
        n_features=16,
        scale=0.5,
        alpha=1.0e-4,
        validation_fraction=0.2,
        rngseed=9,
        min_trials=30,
        top_n=2,
    )
    profiled_payload = json.loads(
        Path(profiled_result["outputs"][0]["json"]).read_text(encoding="utf-8")
    )
    assert profiled_result["mode"] == "profiled"
    assert profiled_result["profiled_shift_count"] == 4
    assert profiled_payload["parameters"][0]["profile_minus"] is not None


# Check direct replica impacts do not fit a redundant observable surrogate
def test_direct_impact_response(tmp_path):
    shifts = impact.fixed_parameter_shifts(
        best_fit=np.array([0.5, 0.5]),
        errors=np.array([0.1, 0.2]),
        bounds=np.array([[0.0, 1.0], [0.0, 1.0]]),
    )

    # Evaluate two exact observable chi2 responses
    def evaluate(points):
        return np.column_stack(
            (
                10.0 + 4.0 * points[:, 0] - points[:, 1],
                3.0 - 0.5 * points[:, 0] + 2.0 * points[:, 1],
            )
        )

    result = impact.render_direct_parameter_impacts(
        evaluate_observable_chi2=evaluate,
        observables=[
            {"dataset": 0, "subset": 0, "observable": "M", "ndf": 3},
            {"dataset": 0, "subset": 0, "observable": "Rap", "ndf": 3},
        ],
        param_names=["coupling", "form_factor"],
        best_fit=np.array([0.5, 0.5]),
        bounds=np.array([[0.0, 1.0], [0.0, 1.0]]),
        shifts=shifts,
        output_dir=tmp_path,
        validation_r2=np.array([0.99, 0.98]),
        validation_mae=np.array([0.1, 0.2]),
        trial_count=100,
        top_n=2,
    )

    assert result["surrogate"]["kind"] == "direct_histogram_replica"
    assert result["published_observable_count"] == 2
    payload = json.loads(Path(result["outputs"][0]["json"]).read_text(encoding="utf-8"))
    assert payload["baseline_chi2"] == pytest.approx(11.5)
    assert payload["baseline_reference"]["kind"] == "histogram_surrogate_optimum"


# Match the regularized solution for tall, wide and rank-deficient feature matrices
@pytest.mark.parametrize("shape", [(17, 4), (5, 13)])
@pytest.mark.parametrize("alpha", [0.0, 0.3])
def test_torch_ridge_matches_svd_reference(shape, alpha):
    rng = np.random.default_rng(103)
    features = rng.normal(size=shape)
    features[:, -1] = features[:, 0]
    target = rng.normal(size=(shape[0], 2)) + np.array([3.0, -2.0])
    phi = np.column_stack((np.ones(shape[0]), features))
    fitted = Ridge(alpha=alpha, solver="svd").fit(features, target)
    expected = np.vstack((fitted.intercept_, fitted.coef_.T))
    actual = impact.fit_ridge_tensor(torch.tensor(phi), torch.tensor(target), alpha).numpy()
    np.testing.assert_allclose(actual, expected, atol=1.0e-12)
    if alpha > 0.0:
        normal_residual = features.T @ (phi @ actual - target) + alpha * actual[1:]
        np.testing.assert_allclose(normal_residual, 0.0, atol=1.0e-12)


# Retain ridge shrinkage when a constant design contains no estimable slopes
def test_torch_ridge_constant_features():
    phi = torch.ones((5, 4), dtype=torch.float64)
    target = torch.arange(5, dtype=torch.float64)
    for alpha in (0.0, 0.2):
        coefficients = impact.fit_ridge_tensor(phi, target, alpha)
        torch.testing.assert_close(coefficients, torch.tensor([[2.0], [0.0], [0.0], [0.0]], dtype=torch.float64))

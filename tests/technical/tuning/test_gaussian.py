# Gaussian amplitude fit objectives, derivatives and covariance normalization
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from types import SimpleNamespace

import numpy as np
import pytest
import torch
from core import resource
from core.numerics import lbfgsb
from core.stats import cov, hist, objective
from core.tune import likelihood
from core.tune.drivers.graniitti.driver import GraniittiDriver
from core.tune.optimizers.ampfit.config import load_settings
from core.tune.search import SearchState


# Compare dense and structured Gaussian likelihoods including physical covariance units
@pytest.mark.parametrize("scale", [1e-80, 1.0, 1e80])
def test_gaussian_covariance(scale):
    block = np.array([[0.8, 0.2, 0.0], [0.2, 0.6, 0.0], [0.0, 0.0, 1.2]])
    vectors = np.array([[0.3], [0.1], [-0.2]])
    covariance = block + vectors @ vectors.T
    mc_error = np.array([0.1, 0.2, 0.3])
    residual = np.array([0.5, -0.3, 0.1])
    total = covariance + np.diag(mc_error**2)
    expected = residual @ np.linalg.solve(total, residual) + np.linalg.slogdet(total)[1] + 6 * np.log(scale)
    params = dict(mc_prediction=residual * scale, data_values=np.zeros(3),
                  mc_stat_uncertainty=mc_error * scale, fit_weights=np.ones(3), gaussian=True)
    structure = cov.prepare_covariance_structure(cov.covariance_decomposition(block * scale**2, vectors * scale))
    for actual in (objective.correlated_chi2_arrays(**params, data_total_covariance=covariance * scale**2),
                   objective.structured_chi2_arrays(**params, covariance_structure=structure)):
        assert actual is not None and actual["valid"]
        assert actual["gaussian"] == pytest.approx(expected)


# Differentiate a normalized spectrum with an exact zero mode and repeated positive eigenvalues
@pytest.mark.parametrize("structured", [False, True])
def test_gaussian_singular_derivatives(structured):
    projector = np.eye(3) - np.ones((3, 3)) / 3
    residual = np.array([0.2, -0.1, -0.1])

    # Vary the covariance on a fixed normalization subspace
    def cost(value):
        covariance = value[0].square() * torch.tensor(projector)
        params = dict(mc_prediction=value[1] * torch.tensor(residual), data_values=np.zeros(3),
                      mc_stat_uncertainty=np.zeros(3), fit_weights=np.ones(3), gaussian=True)
        if structured:
            structure = dict(size=3, diagonal=np.diag(projector).copy(), collective_vectors=np.empty((3, 0)),
                             components=[dict(indices=np.arange(3), covariance=projector)])
            actual = objective.structured_covariance_cost(residual=params["mc_prediction"], structure=structure,
                fit_weights=np.ones(3), covariance_transform=value[0].expand(3), diagonal_uncertainty=np.zeros(3),
                gaussian=True, return_contributions=True)
            result = actual["cost"] + actual["logdet"]
        else:
            actual = objective.correlated_chi2_arrays(**params, data_total_covariance=covariance)
            result = actual["gaussian"]
        assert actual["rank"] == 2
        expected = value[1].square() * float(residual @ residual) / value[0].square() + 4 * value[0].log()
        torch.testing.assert_close(result, expected)
        return result

    value = torch.tensor([0.8, 1.2], dtype=torch.float64, requires_grad=True)
    assert torch.autograd.gradcheck(cost, (value,))
    assert torch.autograd.gradgradcheck(cost, (value,))


# Reject a residual outside the support of a normalized Gaussian
def test_gaussian_outside_support():
    projector = np.eye(3) - np.ones((3, 3)) / 3
    result = objective.correlated_chi2_arrays(mc_prediction=np.ones(3), data_values=np.zeros(3),
        mc_stat_uncertainty=np.zeros(3), fit_weights=np.ones(3), data_total_covariance=projector, gaussian=True)
    assert not result["valid"]


# Evaluate complex interference and model dependent MC variance through the actual driver
@pytest.mark.parametrize("source", [False, True])
@pytest.mark.parametrize("full", [False, True])
def test_gaussian_amplitude(source, full, tmp_path):
    driver = GraniittiDriver()
    bins = np.arange(4)
    centers = hist.edge2centerbins(bins)
    data = hist.hobj(np.array([1.4, 0.7, 1.1]), np.array([0.2, 0.1, 0.2]), bins, centers)
    reference = hist.hobj(np.ones(3), np.full(3, 0.15), bins, centers)
    driver.data = [[{"M": {"hdata": data, "fitw": 1.0}}]]
    driver.data_covariance_payload = dict(layout=cov.build_layout(driver.data),
                                         data_total_covariance=np.diag(data.errs_scaled**2))
    param = dict(cost="gaussian", cost_rho="quadratic", cost_avg="sum", rngseed=0,
                 data_covariance_mode="full" if full else "diagonal",
                 mc_steer={"ampfit": {"mc_stat": "source" if source else "reweighted"}})
    driver.prepare_trial_runtime(param)

    # Construct coherent amplitudes with a physical relative phase
    def cost(value, save=False):
        amplitude = value[0] * torch.exp(1j * value[1]) + value.new_tensor([0.4, 0.6, 0.8])
        prediction = amplitude.abs().square()
        errors = prediction * 0.2
        mc = hist.hobj(prediction, errors, bins, centers)
        result = driver.trial_costs(dict(mc=[[{"M": {"hdata": mc, "source_statistics": reference}}]],
                                        data=driver.data), param)
        variance = value.new_tensor(data.errs_scaled**2) + (0.15**2 if source else errors.square())
        expected = ((prediction - value.new_tensor(data.counts_scaled))**2 / variance + variance.log()).sum()
        torch.testing.assert_close(result["metrics"]["gaussian"], expected)
        assert result["likelihood"]["objective"]["name"] == "gaussian"
        assert result["likelihood"]["components"]["observable_diagnostics"] == "chi2_only"
        assert result["likelihood"]["objective"]["reduced"] is None
        assert result["likelihood"]["mc_uncertainty"]["sigma_objective"] is None
        if save:
            path = tmp_path / "history.json"
            path.write_text(json.dumps(dict(parameter_space=[dict(name="g", lower=0.0, upper=3.0),
                dict(name="phi", lower=-3.0, upper=3.0)], trials=[dict(theta=value.detach().tolist(),
                likelihood=result["likelihood"])])))
            history = likelihood.load_likelihood_history(path)
            assert history["two_nll"][0] == pytest.approx(float(expected.detach()))
            assert history["logL"][0] == pytest.approx(-0.5 * float(expected.detach()))
            assert np.isnan(history["chi2"][0])
            from core.icescape import load_history_surface

            surface = load_history_surface(SimpleNamespace(history_path=str(path), run_name="gaussian",
                                                          cdir=str(tmp_path), max_trials=None))
            assert surface["Z"][0] == pytest.approx(float(expected.detach()))
            assert surface["logL"][0] == pytest.approx(-0.5 * float(expected.detach()))
        return result["metrics"]["gaussian"]

    value = torch.tensor([0.5, 0.3], dtype=torch.float64, requires_grad=True)
    cost(value, save=True)
    assert torch.autograd.gradcheck(cost, (value,))
    assert torch.autograd.gradgradcheck(cost, (value,))


# Fit a variance parameter whose chi square alone would prefer the largest allowed value
def test_gaussian_ampfit_minimum():
    settings = load_settings(resource("tune/settings/ampfit.json"))
    settings.update(starts=1, maxiter=100)
    search = SearchState(args=SimpleNamespace(algorithm="ampfit", cost="gaussian", rngseed=7,
        ampfit_settings=settings), bounds={"sigma": dict(type="uniform", lower=0.2, upper=4.0)},
        initial_points={"sigma": 2.0}, async_proposals=True)
    records = []

    # Evaluate the actual diagonal histogram likelihood
    def cost(value):
        mc = hist.hobj(value.new_tensor([1.0]), value, np.arange(2), np.array([0.5]))
        data = hist.hobj(np.zeros(1), np.zeros(1), np.arange(2), np.array([0.5]))
        return objective.gaussian_cost(mc, data)[0]

    for _ in range(200):
        search.observe(records)
        proposals = search.ask_many(len(records), 1)
        if not proposals:
            break
        search.set_pending([config for config, _ in proposals])
        for config, payload in proposals:
            value, gradient = lbfgsb.value_gradient(cost, torch.tensor([config["sigma"]], dtype=torch.float64))
            records.append(dict(trial_id=str(len(records)), config=config, search_payload=payload,
                metrics={"gaussian": float(value)}, gradient={"sigma": float(gradient[0])}))
    assert search.amplitude.finished
    best = min(records, key=lambda row: row["metrics"]["gaussian"])
    assert best["config"]["sigma"] == pytest.approx(1.0, abs=1e-5)


# Reject nonunit fit weights before amplitude evaluation
def test_gaussian_fit_weights():
    driver = GraniittiDriver()
    driver.data = [[{"M": {"fitw": 2.0}}]]
    with pytest.raises(ValueError, match="unit positive fit weights"):
        driver.prepare_trial_runtime({"cost": "gaussian"})



# Validate Gaussian steering through the actual driver parser checks
@pytest.mark.parametrize("change", [{}, {"algorithm": "hebo"}, {"cost_avg": "global-mean"}, {"cost_rho": "huber"}])
def test_gaussian_cli(change):
    import argparse

    args = SimpleNamespace(algorithm="ampfit", precompute=True, cost="gaussian", cost_rho="quadratic",
                           cost_avg="sum", data_covariance_mode="full", mc_correlation_events=None)
    args.__dict__.update(change)
    if change:
        with pytest.raises(SystemExit):
            GraniittiDriver.validate_cli_arguments(argparse.ArgumentParser(), args)
    else:
        GraniittiDriver.validate_cli_arguments(argparse.ArgumentParser(), args)


# Do not reconstruct a different likelihood from incomplete Gaussian trial histograms
def test_gaussian_histogram_input(tmp_path):
    import pickle

    from core.inference import mcmc

    directory = tmp_path / "runs/icetune/gaussian/results"
    directory.mkdir(parents=True)
    histogram = hist.hobj(np.ones(1), np.ones(1), np.arange(2), np.array([0.5]))
    records = [[{"M": {"hdata": histogram, "fitw": 1.0}}]]
    with (directory / "TUNE_icetune_trial-000000.pkl").open("wb") as stream:
        pickle.dump({"replica_schema_version": 1, "param": {"cost": "gaussian"},
                     "results": {"mc": records, "data": records}}, stream)
    with pytest.raises(ValueError, match="icescape history input"):
        mcmc.collect_simu(run_name="gaussian", cdir=str(tmp_path))


# Keep finite Gaussian values when physical variances would underflow if squared
@pytest.mark.parametrize("scale", [1e-200, 1.0, 1e200])
def test_gaussian_diagonal_scale(scale):
    bins = np.arange(2)
    centers = np.array([0.5])
    mc = hist.hobj(np.array([2.0]) * scale, np.array([0.5]) * scale, bins, centers)
    data = hist.hobj(np.array([1.0]) * scale, np.array([0.5]) * scale, bins, centers)
    value, ndf = objective.gaussian_cost(mc, data)
    assert ndf == 1
    assert value == pytest.approx(2.0 + np.log(0.5) + 2.0 * np.log(scale))



# Reject an invalid Gaussian histogram instead of fitting only the remaining observables
@pytest.mark.parametrize("bad", [np.nan, np.inf])
def test_gaussian_invalid_histogram(bad):
    bins = np.arange(2)
    centers = np.array([0.5])
    good = hist.hobj(np.ones(1), np.ones(1), bins, centers)
    invalid = hist.hobj(np.array([bad]), np.ones(1), bins, centers)
    results = dict(mc=[[{"good": {"hdata": good}, "bad": {"hdata": invalid}}]],
                   data=[[{name: {"hdata": good, "fitw": 1.0} for name in ("good", "bad")}]])
    result = objective.evaluate_cost_bundle(results=results, selected_cost="gaussian", cost_rho="quadratic",
        cost_avg="sum", covariance_payload=None, rngseed=0, wasserstein_cache=None)
    assert np.isposinf(result["metrics"]["gaussian"])
    assert not result["likelihood"]["valid"]


# Honor covariance in the shared per-histogram Gaussian API
def test_gaussian_histogram_covariance():
    bins = np.arange(3)
    centers = hist.edge2centerbins(bins)
    data = hist.hobj(np.array([1.0, 2.0]), np.ones(2), bins, centers)
    mc = hist.hobj(np.array([1.2, 1.7]), np.full(2, 0.1), bins, centers)
    covariance = np.array([[1.0, 0.5], [0.5, 1.0]])
    results = dict(mc=[[{"M": {"hdata": mc}}]], data=[[{"M": {"hdata": data, "fitw": 1.0}}]])
    payload = dict(layout=cov.build_layout(results["data"]), data_total_covariance=covariance)
    costs, _, _, _ = objective.compute_costs(results, "gaussian", "quadratic", covariance_payload=payload)
    residual = mc.counts_scaled - data.counts_scaled
    total = covariance + np.diag(mc.errs_scaled**2)
    assert costs[0][0]["M"] == pytest.approx(residual @ np.linalg.solve(total, residual) + np.linalg.slogdet(total)[1])


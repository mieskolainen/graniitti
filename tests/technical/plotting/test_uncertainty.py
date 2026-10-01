# Weighted covariance, MC support and differentiable normalization tests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
import pytest
import torch
from core.plot import plot
from core.stats import cov, hist, objective, uncertainty


# Fill the real MC histogram path with known event identities and weights
def sample(weights, values=(0.2, 0.7, 1.4), ids=(0, 1, 2), density=False, full=True):
    data = dict(data={"x": np.asarray(values)}, weights=weights, event_ids=np.asarray(ids),
                source="sample", xsection_pb=weights.sum())
    obs = {"x": dict(bins=np.array([0., 1., 2.]), xlabel="x", units={"x": "", "y": ""}, xlim=[0., 2.])}
    return plot.histmc(data, obs, density=density, density_uncertainty="shape",
                       covariance_mode="full" if full else "diagonal")["x"]["hdata"]


# Match both covariance modes to the explicit event normalization Jacobian
@pytest.mark.parametrize("density", [False, True])
def test_weighted_covariance_and_autograd(density):
    weights = torch.tensor([1., 2., 3.], dtype=torch.float64, requires_grad=True)
    histogram = sample(weights, density=density)
    jacobian = torch.autograd.functional.jacobian(lambda w: sample(w, density=density).counts_scaled, weights)
    expected = jacobian @ torch.diag(weights**2) @ jacobian.T
    torch.testing.assert_close(histogram.covariance_scaled, expected)
    torch.testing.assert_close(histogram.errs_scaled**2, expected.diag())
    assert torch.autograd.gradcheck(lambda w: sample(w, density=density).errs_scaled, (weights,))
    numpy_histogram = sample(weights.detach().numpy(), density=density)
    np.testing.assert_allclose(numpy_histogram.covariance_scaled, expected.detach().numpy(), atol=1e-15)


# Shared event identities generate cross-observable covariance including normalization
@pytest.mark.parametrize("density", [False, True])
def test_shared_event_covariance(density):
    weights = np.array([1., 2., 3.])
    first = sample(weights, density=density)
    second = sample(weights, values=(1.2, 0.3, 1.6), density=density)
    joint = uncertainty.joint_covariance([first, second], [np.arange(2), np.arange(2)])
    contributions = np.array([[1., 0., 0., 1.], [2., 0., 2., 0.], [0., 3., 0., 3.]])
    if density:
        for start in (0, 2):
            block = contributions[:, start:start+2]
            contributions[:, start:start+2] = (block - weights[:, None] * block.sum(0) / weights.sum()) / weights.sum()
    np.testing.assert_allclose(joint, contributions.T @ contributions, atol=1e-15)


# Empty chunks retain exposure and occupancy without changing merged uncertainties
def test_chunk_merging_and_support():
    first = sample(np.array([1., 2.]), values=(.2, .7), ids=(0, 1))
    second = sample(np.array([3.]), values=(1.4,), ids=(2,))
    # Identical chunk scales represent equal generation exposure
    merged = first.fuse_independent_chunk(second)
    whole = sample(np.array([1., 2., 3.]))
    np.testing.assert_array_equal(merged.entries, whole.entries)
    np.testing.assert_allclose(merged.counts_scaled, whole.counts_scaled / 2)
    np.testing.assert_allclose(merged.covariance_scaled, whole.covariance_scaled / 4)
    assert uncertainty.mc_statistics(merged)["mc_empty_bin_count"] == 0
    assert uncertainty.mc_statistics(first)["mc_empty_bins"] == [1]
    value, ndf = objective.chi2_cost(first, whole)
    assert np.isfinite(value) and ndf == 2


# Weighted cancellations retain uncertainty without being counted as empty MC
def test_cancellation_is_not_empty():
    histogram = sample(np.array([1., -1., 2.]))
    assert histogram.entries[0] == 2
    assert histogram.errs_scaled[0]**2 == pytest.approx(2.)
    assert uncertainty.mc_statistics(histogram)["mc_empty_bin_count"] == 0


# Full covariance affects both the fit objective and its MC fluctuation estimate
def test_full_fit_covariance():
    mc = np.array([2., 4.])
    data = np.array([1., 1.])
    mc_cov = np.array([[2., 1.], [1., 3.]])
    data_cov = np.eye(2)
    result = objective.correlated_chi2_arrays(mc_prediction=mc, data_values=data,
        mc_stat_uncertainty=np.sqrt(mc_cov.diagonal()), fit_weights=np.ones(2),
        data_total_covariance=data_cov, mc_covariance=mc_cov)
    gradient = np.linalg.solve(data_cov + mc_cov, mc - data)
    assert result["chi2"] == pytest.approx((mc-data) @ gradient)
    assert result["sigma_chi2"] == pytest.approx(2*np.sqrt(gradient @ mc_cov @ gradient))


# Cross sections and histogram integrals use the same weighted Poisson convention
def test_integral_error():
    weights = np.array([1., 2., 3.])
    value, error = uncertainty.sample_cross_section(weights, attempted=10, scale=2)
    histogram = sample(weights)
    assert value == pytest.approx(histogram.integral() * .2)
    assert error == pytest.approx(histogram.integral_error() * .2)



# Measured empty MC bins contribute to the fit while acceptance gaps remain excluded
@pytest.mark.parametrize("full", [False, True])
def test_empty_bin_fit_diagnostics(full):
    prediction = sample(np.array([2.]), values=(.2,), ids=(0,))
    data = hist.hobj(counts=np.array([3., 4.]), errs=np.array([1., 2.]), bins=prediction.bins,
                     cbins=prediction.cbins)
    results = {"mc": [[{"x": {"hdata": prediction}}]],
               "data": [[{"x": {"hdata": data, "fitw": 1.}}]]}
    payload = {"layout": cov.build_layout(results["data"]), "data_total_covariance": data.covariance_scaled} if full else None
    bundle = objective.evaluate_cost_bundle(results=results, selected_cost="chi2", cost_rho="quadratic",
        cost_avg="sum", covariance_payload=payload, rngseed=0, wasserstein_cache=None)
    assert bundle["metrics"]["chi2"] == pytest.approx(1/5 + 16/4)
    record = bundle["likelihood"]["mc_statistics"][0]
    assert record["observable"] == "x"
    assert record["mc_empty_bin_count"] == 1
    assert record["mc_empty_bins"] == [1]
    from core import iceplot
    plotted = iceplot.comparison_metrics(prediction, data)
    assert plotted["mc_empty_bin_count"] == record["mc_empty_bin_count"]
    data.valid[1] = False
    assert uncertainty.comparison_statistics(results)[0]["mc_empty_bin_count"] == 0


# Full normalized covariance remains finite for very small generator weights
@pytest.mark.parametrize("magnitude", [1e-200, 1., 1e150])
def test_shared_density_weight_scale(magnitude):
    weights = np.array([1., 2., 3.])
    first = sample(weights*magnitude, density=True)
    second = sample(weights*magnitude, values=(1.2, .3, 1.6), density=True)
    expected = uncertainty.joint_covariance([sample(weights, density=True),
        sample(weights, values=(1.2, .3, 1.6), density=True)], [np.arange(2), np.arange(2)])
    actual = uncertainty.joint_covariance([first, second], [np.arange(2), np.arange(2)])
    np.testing.assert_allclose(actual, expected, atol=1e-15)


# Frozen source statistics preserve complete covariance without freezing predictions
@pytest.mark.parametrize("density", [False, True])
def test_frozen_source_covariance(density):
    from core.tune.drivers.graniitti.ampfit.amplitude import source_mc

    weights = np.array([1., 2., 3.])
    reference = [sample(weights, density=density), sample(weights, values=(1.2, .3, 1.6), density=density)]
    predictions = [sample(weights*weights, density=density), sample(weights*weights, values=(1.2, .3, 1.6), density=density)]
    results = {"mc": [[{str(i): {"hdata": h, "source_statistics": r}
                         for i, (h, r) in enumerate(zip(predictions, reference, strict=True))}]]}
    frozen = [item["hdata"] for item in source_mc(results)["mc"][0][0].values()]
    selections = [np.arange(2), np.arange(2)]
    np.testing.assert_allclose(uncertainty.joint_covariance(frozen, selections),
                               uncertainty.joint_covariance(reference, selections))
    np.testing.assert_allclose(frozen[0].counts_scaled, predictions[0].counts_scaled)
    assert not np.allclose(frozen[0].counts_scaled, reference[0].counts_scaled)
    assert frozen[0].integral_error() == pytest.approx(reference[0].integral_error())


# Resolving a common error preserves both marginals and creates the expected covariance
def test_split_uncertainty():
    from core.stats.uncertainty import uncertainty_source
    marginal = uncertainty_source('syst', [5., 10.], down=[4., 8.], category='systematic',
                                  correlation='uncorrelated', provenance='Published systematic error')
    common = uncertainty_source('lumi', [3., 6.], category='systematic', correlation='collective')
    sources = uncertainty.split_uncertainty([marginal], [common], 'syst')
    result = uncertainty.finalize_uncertainties({'y': np.array([10., 20.])}, sources)
    np.testing.assert_allclose(result['y_err_up'], [5., 10.])
    np.testing.assert_allclose(result['y_err_down'], [4., 8.])
    assert sum(uncertainty.source_covariance(s)[0, 1] for s in sources) == pytest.approx(18.)
    np.testing.assert_array_equal(marginal['up'], [5., 10.])
    with pytest.raises(ValueError, match='exceed'):
        uncertainty.split_uncertainty([marginal], [common, common], 'syst')


# Signed offsets retain their envelope while the Gaussian covariance follows the Jacobian
@pytest.mark.parametrize('transform', [np.diag([2., 3.]), np.array([[.5, .5]]), np.array([[.3, -.7], [-.3, .7]])])
def test_offset_propagation(transform):
    offsets = np.array([[2., -1.], [-3., 4.]])
    source = uncertainty.offset_source('shift', *offsets.T, 'Published signed variations')
    result = uncertainty.linear_source(source, transform)
    varied = transform @ offsets
    np.testing.assert_allclose(result['covariance'], transform @ source['covariance'] @ transform.T)
    np.testing.assert_allclose(result['up'], np.maximum(varied.max(axis=1), 0))
    np.testing.assert_allclose(result['down'], np.maximum(-varied.min(axis=1), 0))
    np.testing.assert_array_equal(source['offsets'], offsets)


# Preserve the Gaussian covariance and asymmetric margins independently during rebinning
@pytest.mark.parametrize("transform", [np.array([[0.25, 0.75]]), np.array([[0.5, 0.5], [0.0, 1.0]])])
def test_asymmetric_linear_covariance(transform):
    from core.stats.uncertainty import uncertainty_source
    source = uncertainty_source("syst", [2.0, 5.0], down=[1.0, 2.0], category="systematic", correlation="uncorrelated")
    result = uncertainty.linear_source(source, transform, full=True)
    expected = transform @ uncertainty.source_covariance(source) @ transform.T
    np.testing.assert_allclose(uncertainty.source_covariance(result), expected)
    for side in ("up", "down"):
        np.testing.assert_allclose(result[side], np.sqrt(transform**2 @ source[side]**2))
    scale = np.diag(np.arange(2.0, 2.0 + len(transform)))
    converted = uncertainty.linear_source(result, scale)
    np.testing.assert_allclose(uncertainty.source_covariance(converted), scale @ expected @ scale)
    for side in ("up", "down"):
        np.testing.assert_allclose(converted[side], scale @ result[side])


# Integrate differential densities and bin cross sections with their covariance after fusion
@pytest.mark.parametrize("differential", [False, True])
@pytest.mark.parametrize("density", [False, True])
def test_bin_cross_sections(differential, density):
    from core import iceplot

    weights = np.array([1., 2., 3.])
    obs = {"x": dict(bins=np.array([2., 5., 9.]), xlabel="x", units={"x": "", "y": ""},
                     xlim=[2., 9.], differential=differential)}
    events = dict(data={"x": np.array([3., 4., 6.])}, weights=weights, xsection_pb=weights.sum())
    mc = plot.histmc(events, obs)["x"]["hdata"]
    data = uncertainty.finalize_uncertainties(dict(x=mc.cbins, bins=mc.bins, binwidth=mc.binwidth,
        xlim=mc.bins[[0, -1]], y=mc.counts_scaled, scale=1e-12, fitw=1.),
        [uncertainty.uncertainty_source("stat", mc.errs_scaled, category="statistical", correlation="uncorrelated")])
    mc = plot.histmc(events, obs, density=density, density_uncertainty="scaled")["x"]["hdata"]
    measured = plot.histhepdata({"x": data}, obs, density=density, density_uncertainty="scaled")["x"]
    total = 1.0 if density else weights.sum()
    error = np.linalg.norm(weights) / (weights.sum() if density else 1.0)
    np.testing.assert_allclose(mc.counts_scaled, measured["hdata"].counts_scaled)
    for histogram in (mc, measured["hdata"]):
        assert histogram.integral() == pytest.approx(total)
        assert histogram.integral_error() == pytest.approx(error)
        assert iceplot.histogram_integral(histogram) == pytest.approx(total)
        assert iceplot.histogram_integral_error(histogram) == pytest.approx(error)
    fused = mc.fuse_independent_chunk(mc)
    assert fused.integral() == pytest.approx(total)
    assert fused.integral_error() == pytest.approx(error / np.sqrt(2))
    if density:
        return
    assert mc.sum_independent_process(mc).integral() == pytest.approx(2 * weights.sum())
    assert hist.hobj.sum_independent_processes([mc, mc]).integral() == pytest.approx(2 * weights.sum())


# Preserve the null density integral with highly unequal event weights
@pytest.mark.parametrize("small", [1e-3, 1e-8, 1e-12])
@pytest.mark.parametrize("magnitude", [1e-150, 1.0, 1e100])
def test_density_dominant_weight(small, magnitude):
    weights = np.array([1.0, small, small / 10]) * magnitude
    histogram = sample(weights, values=(0.2, 1.2, 1.4), density=True)
    covariance = histogram.covariance_scaled
    response = (np.eye(2)[:, [0, 1, 1]] - histogram.counts_scaled[:, None]) * weights / weights.sum()
    np.testing.assert_allclose(covariance, response @ response.T, rtol=1e-12, atol=64 * np.finfo(float).eps * small)
    np.testing.assert_allclose(histogram.errs_scaled**2, covariance.diagonal())
    assert histogram.integral_error() == pytest.approx(0.0, abs=1e-14)


# Cluster weight reduction sums repeated entries within bins before squaring
@pytest.mark.parametrize("order", [[0, 1, 2, 3, 4], [4, 2, 0, 3, 1]])
def test_cluster_weights(order):
    bins = np.array([0, 0, 1, 0, -1])[order]
    weights = np.array([1., 2., 4., 5., 100.])[order]
    events = np.array([7, 7, 7, 9, 9])[order]
    np.testing.assert_allclose(uncertainty.cluster_weight_squares(bins, weights, events, 3), [34, 16, 0])
    np.testing.assert_array_equal(uncertainty.cluster_weight_squares([], [], [], 3), np.zeros(3))

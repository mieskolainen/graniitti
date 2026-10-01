# Derivatives of weighted histograms, normalization and fit residuals
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
import pytest
import torch
from core.numerics import array
from core.plot import plot
from core.stats import hist as histogram
from core.stats import objective


# Construct weighted events with overflow, bin edges and two physical shape parameters
def prediction(theta, *, density=False, uncertainty="shape", category=False):
    x = np.array([-1.4, -1.0, -0.7, -0.1, 0.0, 0.2, 0.6, 1.5, 1.7])
    basis = np.column_stack((x, x**2))
    xp = array.namespace(theta)
    weights = xp.exp(array.asarray(basis, like=theta) @ theta)
    bins = np.array([-1.0, 0.0, 0.6, 1.5])
    obs = dict(
        bins=bins,
        xlim=[-1.0, 1.5],
        valid=np.array([True, False, True]),
        mc_scale=np.array([0.7, 1.4, 1.1]),
        units={"x": "GeV", "y": "pb"},
        xlabel="Mass",
    )
    values = x
    if category:
        obs.update(
            bins={0: [[0.0, 1.0], [0.0, 1.0], bins]},
            xlim={0: [-1.0, 1.5]},
            valid={0: obs["valid"]},
            symmetrized_fill=False,
            category_labels=["p1", "p2"],
            category_units=["GeV", "GeV"],
        )
        values = np.column_stack((np.linspace(0.1, 1.1, len(x)), np.full(len(x), 0.2), x))
    scale = xp.exp(0.2 * theta[0])
    mcdata = dict(data={"M": values}, weights=weights, xsection_pb=scale * weights.sum())
    return next(
        iter(plot.histmc(mcdata, {"M": obs}, density=density, density_uncertainty=uncertainty, scale=1.3).values())
    )["hdata"]


# Check shared NumPy and Torch equations and automatic first and second derivatives
@pytest.mark.parametrize("density,uncertainty", [(False, "scaled"), (True, "scaled"), (True, "shape")])
@pytest.mark.parametrize("category", [False, True])
@pytest.mark.parametrize("cost", ["chi2", "ratio2"])
def test_histogram_autograd(density, uncertainty, category, cost):
    theta = torch.tensor([0.13, -0.09], dtype=torch.float64, requires_grad=True)
    options = dict(density=density, uncertainty=uncertainty, category=category)
    hist = prediction(theta.detach().numpy(), **options)
    data = histogram.hobj(
        counts=np.array([1.2, 0.0, -0.3]),
        errs=np.array([0.2, 0.0, 0.1]),
        bins=hist.bins,
        cbins=hist.cbins,
        valid=hist.valid,
    )

    # Evaluate the real histogram constructor and cost without explicit derivative formulas
    def evaluate(value):
        result = prediction(value, **options)
        return result.counts_scaled, result.errs_scaled, getattr(objective, cost + "_cost")(result, data)[0]

    actual = evaluate(theta)
    expected = evaluate(theta.detach().numpy())
    for first, second in zip(actual, expected, strict=True):
        np.testing.assert_allclose(first.detach().numpy(), second, rtol=2e-13, atol=1e-14)
    assert torch.autograd.gradcheck(evaluate, (theta,), eps=1e-5)
    assert torch.autograd.gradgradcheck(evaluate, (theta,), eps=1e-5)


# Preserve autograd through independent process sums and event chunk fusion
@pytest.mark.parametrize(
    "operation,density",
    [("sum_independent_processes", False), ("fuse_independent_chunks", False), ("fuse_independent_chunks", True)],
)
def test_combined_histogram_autograd(operation, density):
    theta = torch.tensor([0.13, -0.09], dtype=torch.float64, requires_grad=True)

    # Combine independent samples whose weights have different shape dependences
    def evaluate(value):
        samples = [prediction(value, density=density), prediction(2 * value, density=density)]
        hist = getattr(histogram.hobj, operation)(samples)
        return hist.counts_scaled, hist.errs_scaled

    actual, expected = evaluate(theta), evaluate(theta.detach().numpy())
    for first, second in zip(actual, expected, strict=True):
        np.testing.assert_allclose(first.detach().numpy(), second, rtol=2e-13, atol=1e-14)
    assert torch.autograd.gradcheck(evaluate, (theta,), eps=1e-5)
    assert torch.autograd.gradgradcheck(evaluate, (theta,), eps=1e-5)


# Empty bins and vanishing amplitudes must retain finite automatic derivatives
@pytest.mark.parametrize("values", [[0.0, 0.0, 0.0], [1.0, 0.0, 2.0]])
def test_empty_bin_autograd(values):
    theta = torch.tensor(values, dtype=torch.float64, requires_grad=True)

    # Histogram physical squared amplitudes including a bin with no generated events
    def evaluate(value):
        counts, errors, _, _ = histogram.hist(
            np.array([0.1, 0.2, 0.8]), bins=np.array([0.0, 0.5, 1.0, 1.5]), weights=value**2
        )
        return counts + errors

    gradient = torch.autograd.grad(evaluate(theta).sum(), theta)[0]
    assert torch.all(torch.isfinite(gradient))
    assert torch.autograd.gradcheck(evaluate, (theta,), eps=1e-5)


# Differentiate joint histogram costs through the shared structured covariance and dataset weights
@pytest.mark.parametrize("covariance", [False, True])
@pytest.mark.parametrize("cost", ["chi2", "ratio2"])
@pytest.mark.parametrize("average", ["global-mean", "dataset-mean"])
def test_joint_cost_autograd(covariance, cost, average):
    from core.stats import cov
    from core.stats import objective as icecost

    theta = torch.tensor([0.13, -0.09], dtype=torch.float64, requires_grad=True)
    reference = prediction(theta.detach().numpy())
    data = [
        [
            {
                "M": {
                    "hdata": histogram.hobj(
                        counts=np.array([1.2, 0.0, 0.3]),
                        errs=np.array([0.2, 0.0, 0.1]),
                        bins=reference.bins,
                        cbins=reference.cbins,
                        valid=reference.valid,
                    ),
                    "fitw": 0.7 + index,
                }
            }
        ]
        for index in range(2)
    ]
    payload = None
    if covariance:
        block = np.diag([0.03, 0.007, 0.02, 0.006])
        block[0, 1] = block[1, 0] = 0.002
        vectors = np.array([[0.02], [0.01], [0.03], [0.02]])
        payload = dict(
            layout=cov.build_layout(data),
            data_total_covariance=block + vectors @ vectors.T,
            **cov.covariance_decomposition(block, vectors),
        )

    # Compute the same joint objective used by the generator driver
    def evaluate(value):
        mc = [[{"M": {"hdata": prediction((index + 1) * value)}}] for index in range(2)]
        result = icecost.evaluate_cost_bundle(
            results=dict(mc=mc, data=data),
            selected_cost=cost,
            cost_rho="quadratic",
            cost_avg=average,
            covariance_payload=payload,
            rngseed=0,
            wasserstein_cache=None,
        )
        return result["metrics"][cost]

    np.testing.assert_allclose(evaluate(theta).detach().numpy(), evaluate(theta.detach().numpy()), rtol=1e-12)
    assert torch.autograd.gradcheck(evaluate, (theta,), eps=1e-5)
    assert torch.autograd.gradgradcheck(evaluate, (theta,), eps=1e-5)

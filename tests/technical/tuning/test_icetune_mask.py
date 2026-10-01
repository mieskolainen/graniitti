# Tests for icetune histogram cost validity and reproducibility
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import importlib
import json
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
from core.stats import objective
from core.tune import core as icetune

pytest.importorskip("pyjson5")
from core.plot import loss, plot
from core.stats import cov, hist, uncertainty
from core.stats import objective as icecost
from core.tune import likelihood as icetune_likelihood
from core.tune.drivers.graniitti import driver as graniitti_driver

plot = importlib.reload(plot)
icetune = importlib.reload(icetune)
graniitti_driver = importlib.reload(graniitti_driver)


# Check the changed ratio2 semantics have a stable persisted definition
def test_ratio2_cost_definition_is_versioned():
    assert icetune.optimizer_cost_definition("ratio2") == "symmetric_asinh_ratio_chi2_v3"
    assert icetune.optimizer_cost_definition("chi2") == "chi2"


# Check numerical thread limits match one trial's assigned CPU count
def test_numeric_thread_limits():
    result = subprocess.run([sys.executable, '-c', """
import os
import numpy as np
from threadpoolctl import threadpool_info
from core.tune.runtime.process import configure_numerical_threads
import torch
torch.set_num_threads(1)
assert configure_numerical_threads(3) == 3
assert torch.get_num_threads() == 3
np.eye(2) @ np.eye(2)
pools = [pool for pool in threadpool_info() if pool['user_api'] == 'blas']
assert pools and all(pool['num_threads'] == 3 for pool in pools)
for name in ('BLIS', 'MKL', 'NUMEXPR', 'OMP', 'OPENBLAS', 'VECLIB_MAXIMUM'):
    assert os.environ[name + '_NUM_THREADS' if name != 'VECLIB_MAXIMUM' else name + '_THREADS'] == '3'
"""], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr


# Build a compact histogram fixture
def _hobj(counts, errs, valid=None):
    counts = np.asarray(counts, dtype=float)
    errs = np.asarray(errs, dtype=float)
    bins = np.arange(len(counts) + 1, dtype=float)
    cbins = 0.5 * (bins[:-1] + bins[1:])
    return hist.hobj(counts=counts, errs=errs, bins=bins, cbins=cbins, valid=valid)


def _hist_entry(hdata, label):
    return {
        "hdata": hdata,
        "obs": {
            "xlabel": "x",
            "ylabel": "y",
            "units": {"x": "", "y": ""},
            "ylim": None,
            "xlim": (0.0, float(len(hdata.counts))),
            "density": False,
            "figsize": (4, 3),
        },
        "label": label,
        "color": None,
        "hfunc": "hist",
        "style": {},
    }


# Check explicit observable validity for finite, empty, zero-MC and malformed inputs
def test_compute_costs_explicit_validity_mask():
    results = {
        "mc": [
            [
                {
                    "valid_obs": {"hdata": _hobj([10.0, 20.0], [1.0, 1.0])},
                    "zero_mc_obs": {"hdata": _hobj([0.0, 0.0], [1.0, 1.0])},
                    "zero_ndf_obs": {"hdata": _hobj([0.0, 0.0], [0.0, 0.0])},
                    "nan_cost_obs": {"hdata": _hobj([10.0, 20.0], [np.nan, 1.0])},
                }
            ]
        ],
        "data": [
            [
                {
                    "valid_obs": {"hdata": _hobj([12.0, 18.0], [1.0, 1.0]), "fitw": 1.0},
                    "zero_mc_obs": {"hdata": _hobj([5.0, 5.0], [1.0, 1.0]), "fitw": 1.0},
                    "zero_ndf_obs": {"hdata": _hobj([0.0, 0.0], [0.0, 0.0]), "fitw": 1.0},
                    "nan_cost_obs": {"hdata": _hobj([12.0, 18.0], [1.0, 1.0]), "fitw": 1.0},
                }
            ]
        ],
    }

    cost_all, ndf_all, weight_all, valid_all = icecost.compute_costs(
        results=results,
        cost_func="chi2",
        cost_rho="quadratic",
    )

    assert valid_all[0][0]["valid_obs"] == 1
    assert valid_all[0][0]["zero_mc_obs"] == 1
    assert valid_all[0][0]["zero_ndf_obs"] == 0
    assert valid_all[0][0]["nan_cost_obs"] == 0

    assert ndf_all[0][0]["valid_obs"] > 0
    assert ndf_all[0][0]["zero_mc_obs"] == 2
    assert cost_all[0][0]["zero_mc_obs"] == pytest.approx(25.0)
    assert ndf_all[0][0]["zero_ndf_obs"] == 0
    assert np.isnan(cost_all[0][0]["nan_cost_obs"])
    assert weight_all[0][0]["valid_obs"] == pytest.approx(1.0)


# Check covariance-aware histogram diagnostics use each marginal covariance block
def test_compute_costs_marginal_cov_for_chi2():
    results = {
        "mc": [[{"obs": {"hdata": _hobj([11.0, 18.0], [0.5, 1.0])}}]],
        "data": [[{"obs": {"hdata": _hobj([10.0, 20.0], [1.0, 2.0]), "fitw": 1.5}}]],
    }
    covariance = np.array([[1.0, 0.4], [0.4, 4.0]], dtype=float)
    payload = {
        "layout": [
            {
                "dataset": 0,
                "subset": 0,
                "observable": "obs",
                "bins": [0, 1],
                "start": 0,
                "stop": 2,
            }
        ],
        "data_total_covariance": covariance,
    }

    cost, ndf, weight, valid = icecost.compute_costs(
        results=results,
        cost_func="chi2",
        cost_rho="quadratic",
        covariance_payload=payload,
    )
    residual = np.array([1., -2.])
    expected = residual @ np.linalg.solve(covariance + np.diag([.25, 1.]), residual)
    assert cost[0][0]['obs'] == pytest.approx(expected)
    assert ndf[0][0]['obs'] == 2
    assert weight[0][0]["obs"] == pytest.approx(1.5)
    assert valid[0][0]["obs"] == 1


def test_combine_costs_mask_inf_when_all_invalid():
    cost = [[{"valid_obs": 10.0, "masked_obs": 1000.0}]]
    ndf = [[{"valid_obs": 5.0, "masked_obs": 10.0}]]
    weight = [[{"valid_obs": 1.0, "masked_obs": 1.0}]]
    valid = [[{"valid_obs": 1, "masked_obs": 0}]]

    value = icetune_likelihood.combine_costs(
        cost=cost,
        ndf=ndf,
        weight=weight,
        valid=valid,
        cost_avg="global-mean-unweighted",
    )

    assert value == pytest.approx(2.0)

    all_invalid = [[{"valid_obs": 0, "masked_obs": 0}]]
    value = icetune_likelihood.combine_costs(
        cost=cost,
        ndf=ndf,
        weight=weight,
        valid=all_invalid,
        cost_avg="global-mean-unweighted",
    )

    assert np.isinf(value)


# Check dataset normalization preserves relative weights and ignores histogram multiplicity
@pytest.mark.parametrize("copies", [1, 20])
def test_dataset_mean_balances_costs(copies):
    cost = [[{"a": 8.0, "b": 2.0}], [{"c": 900.0, "masked": 1e6}] * copies]
    ndf = [[{"a": 2.0, "b": 1.0}], [{"c": 100.0, "masked": 1000.0}] * copies]
    weight = [[{"a": 3.0, "b": 1.0}], [{"c": 7.0, "masked": 0.0}] * copies]
    valid = [[{"a": 1, "b": 1}], [{"c": 1, "masked": 1}] * copies]

    value = icetune_likelihood.combine_costs(cost, ndf, weight, valid, "dataset-mean")
    assert value == pytest.approx(((3 * 8 + 2) / (3 * 2 + 1) + 9) / 2)
    cost[0][0]["a"] = np.inf
    assert np.isinf(icetune_likelihood.combine_costs(cost, ndf, weight, valid, "dataset-mean"))
    valid = [[{name: 0 for name in subset} for subset in dataset] for dataset in valid]
    assert np.isinf(icetune_likelihood.combine_costs(cost, ndf, weight, valid, "dataset-mean"))


# Check balanced driver costs against independent dataset costs with real histogram covariance
@pytest.mark.parametrize("selected_cost", ["chi2", "ratio2"])
@pytest.mark.parametrize("full", [False, True])
def test_cost_bundle_balances_datasets(selected_cost, full):
    results = {
        "mc": [[{"a": {"hdata": _hobj([12.0], [0.4])}}],
               [{"b": {"hdata": _hobj([7.0, 9.0, 11.0], [0.3, 0.3, 0.3])}}]],
        "data": [[{"a": {"hdata": _hobj([10.0], [1.0]), "fitw": 2.0}}],
                 [{"b": {"hdata": _hobj([8.0, 8.0, 8.0], [0.8, 0.8, 0.8]), "fitw": 4.0}}]],
    }
    payload = {
        "layout": cov.build_layout(results["data"]),
        "data_total_covariance": np.diag([1.0, 0.8**2, 0.8**2, 0.8**2]),
    }
    cost, ndf, _, _ = icecost.compute_costs(results=results, cost_func=selected_cost, cost_rho="quadratic")
    expected = (cost[0][0]["a"] / ndf[0][0]["a"] + cost[1][0]["b"] / ndf[1][0]["b"]) / 2
    bundle = icecost.evaluate_cost_bundle(
        results=results, selected_cost=selected_cost, cost_rho="quadratic", cost_avg="dataset-mean",
        covariance_payload=payload if full else None, rngseed=1234, wasserstein_cache=None,
    )
    assert bundle["metrics"][selected_cost] == pytest.approx(expected)
    chi2, _, _, _ = icecost.compute_costs(results=results, cost_func="chi2", cost_rho="quadratic")
    assert bundle["likelihood"]["objective"]["value"] == pytest.approx(
        2 * chi2[0][0]["a"] + 4 * chi2[1][0]["b"]
    )


def test_likelihood_totals_and_obs():
    cost = [[{"a": 10.0, "b": 100.0, "bad": np.nan}]]
    ndf = [[{"a": 5.0, "b": 10.0, "bad": 2.0}]]
    weight = [[{"a": 2.0, "b": 3.0, "bad": 1.0}]]
    valid = [[{"a": 1, "b": 0, "bad": 0}]]

    likelihood = icetune_likelihood.build_likelihood_payload(
        cost=cost,
        ndf=ndf,
        weight=weight,
        valid=valid,
        rho="quadratic",
    )

    assert likelihood["valid"] is True
    assert likelihood["kind"] == "gaussian_chi2"
    assert likelihood["covariance_mode"] == "diagonal"
    assert likelihood["residual_model"] == "quadratic"
    assert likelihood["objective"]["value"] == pytest.approx(20.0)
    assert likelihood["objective"]["ndf"] == pytest.approx(5.0)
    assert likelihood["objective"]["reduced"] == pytest.approx(4.0)
    assert len(likelihood["observables"]) == 3
    assert likelihood["observables"][2]["objective_value"] is None
    assert "nll_total" not in likelihood
    assert "logL" not in likelihood
    assert likelihood["mc_uncertainty"]["replicate_count"] == 1


def test_metric_no_implied_logl():
    likelihood = icetune_likelihood.build_metric_likelihood_payload(
        metrics={"loss": 3.0}, cost_key="loss"
    )

    assert likelihood["valid"] is False
    assert likelihood["kind"] == "optimizer_cost"
    assert likelihood["objective"]["name"] == "loss"
    assert likelihood["objective"]["value"] == pytest.approx(3.0)
    assert "nll_total" not in likelihood
    assert "logL" not in likelihood


def test_param_space_theta_and_hash_stable():
    bounds = {
        "z": {"lower": -1.0, "upper": 1.0},
        "a": {"lower": 0.0, "upper": 2.0},
    }
    parameter_space = icetune.build_parameter_space(bounds)
    theta, fingerprint = icetune.build_theta_record(
        {"z": 0.5, "a": 1.25},
        parameter_space,
    )

    assert [item["name"] for item in parameter_space] == ["a", "z"]
    assert theta == pytest.approx([1.25, 0.5])
    assert fingerprint == icetune.theta_hash([1.25, 0.5])
    assert fingerprint != icetune.theta_hash([0.5, 1.25])


def test_ratio2_cost_finite_sparse_mc_bins():
    h_mc = _hobj([0.0], [0.0])
    h_data = _hobj([10.0], [np.sqrt(10.0)])

    cost, ndf = objective.ratio2_cost(h_mc=h_mc, h_data=h_data, rho="quadratic")
    nearby, _ = objective.ratio2_cost(
        h_mc=_hobj([1.0e-12], [0.0]),
        h_data=h_data,
        rho="quadratic",
    )

    assert ndf == 1
    assert np.isfinite(cost)
    assert nearby == pytest.approx(cost, rel=1.0e-12)


def test_ratio2_common_scale_invariance():
    h_mc_a = _hobj([10.0], [1.0])
    h_data_a = _hobj([12.0], [1.2])
    h_mc_b = _hobj([1000.0], [100.0])
    h_data_b = _hobj([1200.0], [120.0])

    cost_a, ndf_a = objective.ratio2_cost(h_mc=h_mc_a, h_data=h_data_a, rho="quadratic")
    cost_b, ndf_b = objective.ratio2_cost(h_mc=h_mc_b, h_data=h_data_b, rho="quadratic")

    assert ndf_a == 1
    assert ndf_b == 1
    assert cost_a == pytest.approx(cost_b)


# Check that exchanging MC and data changes only the asinh-ratio pull sign
def test_ratio2_cost_sym_under_mc_data_exchange():
    h_mc = _hobj([0.0, 4.0, 15.0], [0.0, 0.8, 1.5])
    h_data = _hobj([5.0, 2.0, 3.0], [1.0, 0.4, 0.3])

    forward, forward_ndf = objective.ratio2_cost(
        h_mc=h_mc,
        h_data=h_data,
        rho="quadratic",
    )
    reverse, reverse_ndf = objective.ratio2_cost(
        h_mc=h_data,
        h_data=h_mc,
        rho="quadratic",
    )

    assert forward_ndf == 3
    assert reverse_ndf == forward_ndf
    assert reverse == pytest.approx(forward)


# Check that the multiplicative pull approaches ordinary chi-square near agreement
def test_ratio2_cost_additive_chi2_local_limit():
    h_mc = _hobj([100.001], [2.0])
    h_data = _hobj([100.0], [3.0])

    ratio2, ratio_ndf = objective.ratio2_cost(h_mc=h_mc, h_data=h_data)
    chi2, chi_ndf = objective.chi2_cost(h_mc=h_mc, h_data=h_data)

    assert ratio_ndf == chi_ndf == 1
    assert ratio2 == pytest.approx(chi2, rel=2e-5)


# Check that a certain nonzero disagreement rejects the trial
def test_ratio2_zero_variance_conflict():
    h_mc = _hobj([1.0], [0.0])
    h_data = _hobj([2.0], [0.0])

    cost, ndf = objective.ratio2_cost(h_mc=h_mc, h_data=h_data)

    assert ndf == 1
    assert np.isposinf(cost)


# Check that chi2 keeps zero MC against positive data in the objective
def test_chi2_zero_mc_bins():
    results = {
        "mc": [
            [
                {
                    "obs": {"hdata": _hobj([0.0, 0.0], [0.0, 0.0])},
                }
            ]
        ],
        "data": [
            [
                {
                    "obs": {"hdata": _hobj([5.0, 0.0], [1.0, 0.0]), "fitw": 1.0},
                }
            ]
        ],
    }

    cost_all, ndf_all, weight_all, valid_all = icecost.compute_costs(
        results=results,
        cost_func="chi2",
        cost_rho="quadratic",
    )

    assert valid_all[0][0]["obs"] == 1
    assert ndf_all[0][0]["obs"] == 1
    assert cost_all[0][0]["obs"] == pytest.approx(25.0)
    assert weight_all[0][0]["obs"] == pytest.approx(1.0)


# Check that additive and ratio costs retain either direction of zero mismatch
@pytest.mark.parametrize("cost_func", ["chi2", "ratio2"])
@pytest.mark.parametrize(
    ("mc_counts", "mc_errs", "data_counts", "data_errs"),
    [
        ([0.0, 0.0], [0.0, 0.0], [5.0, 0.0], [1.0, 0.0]),
        ([5.0, 0.0], [1.0, 0.0], [0.0, 0.0], [0.0, 0.0]),
    ],
)
def test_compute_costs_sided_zero_bins(
    cost_func,
    mc_counts,
    mc_errs,
    data_counts,
    data_errs,
):
    results = {
        "mc": [
            [
                {
                    "obs": {"hdata": _hobj(mc_counts, mc_errs)},
                }
            ]
        ],
        "data": [
            [
                {
                    "obs": {"hdata": _hobj(data_counts, data_errs), "fitw": 1.0},
                }
            ]
        ],
    }

    cost_all, ndf_all, weight_all, valid_all = icecost.compute_costs(
        results=results,
        cost_func=cost_func,
        cost_rho="quadratic",
    )

    assert valid_all[0][0]["obs"] == 1
    assert ndf_all[0][0]["obs"] == 1
    assert np.isfinite(cost_all[0][0]["obs"])
    assert cost_all[0][0]["obs"] > 0.0
    assert weight_all[0][0]["obs"] == pytest.approx(1.0)


# Check that Wasserstein retains zero mismatches and reports informative bins
@pytest.mark.parametrize(
    ("mc_counts", "data_counts"),
    [
        ([0.0, 0.0], [5.0, 0.0]),
        ([5.0, 0.0], [0.0, 0.0]),
    ],
)
def test_wasserstein_one_sided_zeros(mc_counts, data_counts):
    results = {
        "mc": [
            [
                {
                    "obs": {"hdata": _hobj(mc_counts, [0.0, 0.0])},
                }
            ]
        ],
        "data": [
            [
                {
                    "obs": {"hdata": _hobj(data_counts, [0.0, 0.0]), "fitw": 1.0},
                }
            ]
        ],
    }

    cost_all, ndf_all, weight_all, valid_all = icecost.compute_costs(
        results=results,
        cost_func="wasserstein",
        cost_rho="quadratic",
    )

    assert valid_all[0][0]["obs"] == 1
    assert ndf_all[0][0]["obs"] == 1
    # Integrate the linear CDF rise over the first bin and its plateau over the second
    assert cost_all[0][0]["obs"] == pytest.approx(7.5)
    assert weight_all[0][0]["obs"] == pytest.approx(1.0)


# Check Wasserstein toys use reproducible common random numbers
def test_wasserstein_cost_reproducible_seed():
    h_mc = _hobj([8.0, 12.0], [0.8, 1.2])
    h_data = _hobj([10.0, 10.0], [1.0, 1.0])

    first, _ = objective.wasserstein_cost(h_mc, h_data, N=64, rngseed=1234)
    repeated, _ = objective.wasserstein_cost(h_mc, h_data, N=64, rngseed=1234)
    changed, _ = objective.wasserstein_cost(h_mc, h_data, N=64, rngseed=4321)

    assert first == repeated
    assert first != changed


# Check Wasserstein rejects signed measures instead of omitting them
def test_wasserstein_cost_rejects_negative_bin():
    cost, ndf = objective.wasserstein_cost(
        _hobj([-1.0], [0.2]),
        _hobj([1.0], [0.2]),
    )

    assert ndf == 1
    assert np.isposinf(cost)


# Check Wasserstein toys retain the full within-histogram data covariance
def test_wasserstein_cost_data_total_cov():
    h_mc = _hobj([8.0, 12.0], [0.8, 1.2])
    h_data = _hobj([10.0, 10.0], [1.0, 1.0])
    diagonal = np.eye(2)
    correlated = np.array([[1.0, 0.9], [0.9, 1.0]])

    diagonal_cost, _ = objective.wasserstein_cost(
        h_mc,
        h_data,
        N=128,
        rngseed=1234,
        data_total_covariance=diagonal,
    )
    correlated_cost, _ = objective.wasserstein_cost(
        h_mc,
        h_data,
        N=128,
        rngseed=1234,
        data_total_covariance=correlated,
    )

    assert diagonal_cost != correlated_cost


# Check cached Wasserstein replicas avoid repeated covariance factorization
def test_wasserstein_replica_cache(monkeypatch):
    h_mc = _hobj([8.0, 12.0], [0.8, 1.2])
    h_data = _hobj([10.0, 10.0], [1.0, 1.0])
    covariance = np.array([[1.0, 0.9], [0.9, 1.0]])
    state = objective.build_wasserstein_replica_state(
        data_mean=h_data.counts_scaled,
        data_total_covariance=covariance,
        sample_count=64,
        rngseed=1234,
    )
    expected, _ = objective.wasserstein_cost(
        h_mc,
        h_data,
        N=64,
        replica_state=state,
    )

    def reject_factorization(_):
        raise AssertionError("cached Wasserstein evaluation refactorized covariance")

    monkeypatch.setattr(uncertainty, "covariance_square_root", reject_factorization)
    repeated, _ = objective.wasserstein_cost(
        h_mc,
        h_data,
        N=64,
        replica_state=state,
    )

    assert repeated == expected
    assert not state["mc_standard"].flags.writeable
    assert not state["data_replicas"].flags.writeable


# Check each correlated trial evaluates only its required covariance objectives
@pytest.mark.parametrize(
    ("selected_cost", "expected_cost_calls", "expected_eigh_calls"),
    [
        ("chi2", [], 1),
        ("ratio2", ["ratio2"], 2),
        ("wasserstein", ["wasserstein"], 1),
    ],
)
def test_cost_reuses_factorization(
    monkeypatch,
    selected_cost,
    expected_cost_calls,
    expected_eigh_calls,
):
    results = {
        "mc": [[{"obs": {"hdata": _hobj([11.0, 18.0], [0.5, 1.0])}}]],
        "data": [[{"obs": {"hdata": _hobj([10.0, 20.0], [1.0, 2.0]), "fitw": 1.0}}]],
    }
    payload = {
        "layout": cov.build_layout(results["data"]),
        "data_total_covariance": np.array([[1.0, 0.4], [0.4, 4.0]]),
    }
    cache = icecost.build_wasserstein_replica_cache(
        data=results["data"],
        covariance_payload=payload,
        rngseed=1234,
    )
    cost_calls = []
    eigh_calls = 0
    original_compute_costs = icecost.compute_costs
    original_eigh = np.linalg.eigh

    def tracked_compute_costs(*args, **kwargs):
        cost_calls.append(kwargs["cost_func"])
        return original_compute_costs(*args, **kwargs)

    def tracked_eigh(*args, **kwargs):
        nonlocal eigh_calls
        eigh_calls += 1
        return original_eigh(*args, **kwargs)

    monkeypatch.setattr(icecost, "compute_costs", tracked_compute_costs)
    monkeypatch.setattr(np.linalg, "eigh", tracked_eigh)
    bundle = icecost.evaluate_cost_bundle(
        results=results,
        selected_cost=selected_cost,
        cost_rho="quadratic",
        cost_avg="global-mean",
        covariance_payload=payload,
        rngseed=1234,
        wasserstein_cache=cache,
    )

    assert cost_calls == expected_cost_calls
    assert eigh_calls == expected_eigh_calls
    assert set(bundle["metrics"]) == {"chi2", selected_cost}
    if selected_cost == "ratio2":
        joint_ratio2 = objective.correlated_ratio2(results=results, payload=payload)
        assert bundle["metrics"]["ratio2"] == pytest.approx(
            joint_ratio2["cost"] / joint_ratio2["ndf"]
        )


# Check full and diagonal ratio2 use the same relative fit weight convention
@pytest.mark.parametrize(
    "cost_avg",
    ["global-mean", "local-mean", "global-mean-unweighted", "local-mean-unweighted", "dataset-mean"],
)
def test_ratio2_diagonal_fit_weights(cost_avg):
    results = {
        "mc": [
            [
                {
                    "a": {"hdata": _hobj([12.0, 18.0], [0.4, 0.6])},
                    "b": {"hdata": _hobj([7.0], [0.3])},
                }
            ]
        ],
        "data": [
            [
                {
                    "a": {"hdata": _hobj([10.0, 20.0], [1.0, 1.5]), "fitw": 2.0},
                    "b": {"hdata": _hobj([8.0], [0.8]), "fitw": 4.0},
                }
            ]
        ],
    }
    payload = {
        "layout": cov.build_layout(results["data"]),
        "data_total_covariance": np.diag([1.0, 1.5**2, 0.8**2]),
    }
    diagonal = icecost.evaluate_cost_bundle(
        results=results,
        selected_cost="ratio2",
        cost_rho="quadratic",
        cost_avg=cost_avg,
        covariance_payload=None,
        rngseed=1234,
        wasserstein_cache=None,
    )
    full = icecost.evaluate_cost_bundle(
        results=results,
        selected_cost="ratio2",
        cost_rho="quadratic",
        cost_avg=cost_avg,
        covariance_payload=payload,
        rngseed=1234,
        wasserstein_cache=None,
    )

    assert full["metrics"]["ratio2"] == pytest.approx(diagonal["metrics"]["ratio2"])


# Check ratio2 fit weights are relative while Gaussian likelihood weights are absolute
@pytest.mark.parametrize("cost_rho", ["quadratic", "huber"])
def test_full_cost_fit_weight_scaling(cost_rho):
    results = {
        "mc": [
            [
                {
                    "a": {"hdata": _hobj([12.0], [0.4])},
                    "b": {"hdata": _hobj([7.0], [0.3])},
                }
            ]
        ],
        "data": [
            [
                {
                    "a": {"hdata": _hobj([10.0], [1.0]), "fitw": 1.0},
                    "b": {"hdata": _hobj([8.0], [0.8]), "fitw": 2.0},
                }
            ]
        ],
    }
    payload = {
        "layout": cov.build_layout(results["data"]),
        "data_total_covariance": np.array([[1.0, 0.25], [0.25, 0.8**2]]),
    }

    # Evaluate one bundle using the current mutable fit weights
    def evaluate():
        return icecost.evaluate_cost_bundle(
            results=results,
            selected_cost="ratio2",
            cost_rho=cost_rho,
            cost_avg="global-mean",
            covariance_payload=payload,
            rngseed=1234,
            wasserstein_cache=None,
        )

    unit = evaluate()
    results["data"][0][0]["a"]["fitw"] *= 3.0
    results["data"][0][0]["b"]["fitw"] *= 3.0
    scaled = evaluate()

    assert scaled["metrics"]["ratio2"] == pytest.approx(unit["metrics"]["ratio2"])
    assert scaled["metrics"]["chi2"] == pytest.approx(3.0 * unit["metrics"]["chi2"])
    assert scaled["likelihood"]["objective"]["value"] == pytest.approx(
        3.0 * unit["likelihood"]["objective"]["value"]
    )


# Check zero fit weight removes a signed ratio bin without invalidating the joint cost
@pytest.mark.parametrize("cost_avg", ["global-mean", "global-mean-unweighted"])
def test_full_ratio2_ignores_zero_weight_signed_obs(cost_avg):
    results = {
        "mc": [
            [
                {
                    "active": {"hdata": _hobj([12.0], [0.4])},
                    "off": {"hdata": _hobj([-2.0], [0.3])},
                }
            ]
        ],
        "data": [
            [
                {
                    "active": {"hdata": _hobj([10.0], [1.0]), "fitw": 1.0},
                    "off": {"hdata": _hobj([-1.0], [0.8]), "fitw": 0.0},
                }
            ]
        ],
    }
    payload = {
        "layout": cov.build_layout(results["data"]),
        "data_total_covariance": np.diag([1.0, 0.8**2]),
    }
    bundle = icecost.evaluate_cost_bundle(
        results=results,
        selected_cost="ratio2",
        cost_rho="quadratic",
        cost_avg=cost_avg,
        covariance_payload=payload,
        rngseed=1234,
        wasserstein_cache=None,
    )

    assert np.isfinite(bundle["metrics"]["ratio2"])
    assert bundle["likelihood"]["objective"]["ndf"] == pytest.approx(1.0)


# Check signed measurements remain in the asinh ratio cost
def test_full_ratio2_signed_obs_joint_cost():
    results = {
        "mc": [
            [
                {
                    "active": {"hdata": _hobj([12.0], [0.4])},
                    "signed": {"hdata": _hobj([-2.0], [0.3])},
                }
            ]
        ],
        "data": [
            [
                {
                    "active": {"hdata": _hobj([10.0], [1.0]), "fitw": 1.0},
                    "signed": {"hdata": _hobj([-1.0], [0.8]), "fitw": 1.0},
                }
            ]
        ],
    }
    payload = {
        "layout": cov.build_layout(results["data"]),
        "data_total_covariance": np.diag([1.0, 0.8**2]),
    }
    joint = objective.correlated_ratio2(results=results, payload=payload)
    bundle = icecost.evaluate_cost_bundle(
        results=results,
        selected_cost="ratio2",
        cost_rho="quadratic",
        cost_avg="global-mean",
        covariance_payload=payload,
        rngseed=1234,
        wasserstein_cache=None,
    )
    independent, ndf = objective.ratio2_cost(
        h_mc=results["mc"][0][0]["active"]["hdata"],
        h_data=results["data"][0][0]["active"]["hdata"],
    )
    signed, signed_ndf = objective.ratio2_cost(
        h_mc=results["mc"][0][0]["signed"]["hdata"],
        h_data=results["data"][0][0]["signed"]["hdata"],
    )

    assert joint["valid"]
    assert joint["cost"] == pytest.approx(independent + signed)
    assert joint["ndf"] == ndf + signed_ndf
    assert [item["valid"] for item in joint["observables"]] == [True, True]
    assert bundle["metrics"]["ratio2"] == pytest.approx(joint["cost"] / joint["ndf"])
    assert bundle["likelihood"]["objective"]["ndf"] == pytest.approx(2.0)


# Check empty joint records retain their data fit weight without becoming active
def test_joint_cost_arrays_handle_empty_obs_weight():
    results = {
        "data": [
            [
                {
                    "empty": {"hdata": _hobj([], []), "fitw": 0.2},
                }
            ]
        ],
    }
    joint = {
        "observables": [
            {
                "dataset": 0,
                "subset": 0,
                "observable": "empty",
                "ndf": None,
                "valid": False,
                "objective_value": None,
                "weight": None,
            }
        ]
    }

    cost, ndf, weight, valid = icecost.joint_observable_cost_arrays(results, joint)

    assert np.isnan(cost[0][0]["empty"])
    assert ndf[0][0]["empty"] == 0.0
    assert weight[0][0]["empty"] == pytest.approx(0.2)
    assert valid[0][0]["empty"] == 0


# Check that every histogram cost ignores only explicitly masked or empty bins
@pytest.mark.parametrize(
    ("cost_func", "expected_ndf"),
    [("chi2", 1), ("ratio2", 1), ("wasserstein", 1)],
)
def test_compute_costs_share_explicit_bin_validity(cost_func, expected_ndf):
    results = {
        "mc": [
            [
                {
                    "obs": {
                        "hdata": _hobj([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], valid=[True, True, False])
                    },
                }
            ]
        ],
        "data": [
            [
                {
                    "obs": {
                        "hdata": _hobj([5.0, 0.0, 7.0], [1.0, 0.0, 1.0], valid=[True, True, True]),
                        "fitw": 1.0,
                    },
                }
            ]
        ],
    }

    cost_all, ndf_all, _, valid_all = icecost.compute_costs(
        results=results,
        cost_func=cost_func,
        cost_rho="quadratic",
    )

    assert valid_all[0][0]["obs"] == 1
    assert ndf_all[0][0]["obs"] == expected_ndf
    assert np.isfinite(cost_all[0][0]["obs"])
    assert cost_all[0][0]["obs"] > 0.0


# Check that malformed active bins invalidate every histogram cost
@pytest.mark.parametrize("cost_func", ["chi2", "ratio2", "wasserstein"])
def test_compute_costs_reject_nonfinite_active_bins(cost_func):
    results = {
        "mc": [
            [
                {
                    "obs": {"hdata": _hobj([5.0, 1.0], [np.nan, 1.0])},
                }
            ]
        ],
        "data": [
            [
                {
                    "obs": {"hdata": _hobj([4.0, 1.0], [1.0, 1.0]), "fitw": 1.0},
                }
            ]
        ],
    }

    cost_all, ndf_all, _, valid_all = icecost.compute_costs(
        results=results,
        cost_func=cost_func,
        cost_rho="quadratic",
    )

    assert valid_all[0][0]["obs"] == 0
    assert ndf_all[0][0]["obs"] == 0
    assert np.isnan(cost_all[0][0]["obs"])


# Check certain disagreement rejects the model
def test_chi2_zero_error_disagreement_infinite():
    cost, ndf = objective.chi2_cost(
        h_mc=_hobj([0.0], [0.0]),
        h_data=_hobj([5.0], [0.0]),
        rho="quadratic",
    )

    assert ndf == 1
    assert np.isposinf(cost)


# Check a certain disagreement dominates usable chi-square bins
def test_chi2_zero_variance_disagreement():
    cost, ndf = objective.chi2_cost(
        h_mc=_hobj([0.0, 1.0], [0.0, 1.0]),
        h_data=_hobj([5.0, 2.0], [0.0, 1.0]),
        rho="quadratic",
    )

    assert ndf == 2
    assert np.isposinf(cost)


# Check a certain disagreement rejects an icetune trial
def test_cost_zero_variance_conflict():
    results = {
        "mc": [
            [
                {
                    "obs": {"hdata": _hobj([0.0, 1.0], [0.0, 1.0])},
                }
            ]
        ],
        "data": [
            [
                {
                    "obs": {"hdata": _hobj([5.0, 2.0], [0.0, 1.0]), "fitw": 1.0},
                }
            ]
        ],
    }

    cost_all, ndf_all, weight_all, valid_all = icecost.compute_costs(
        results=results,
        cost_func="chi2",
        cost_rho="quadratic",
    )
    combined = icetune_likelihood.combine_costs(
        cost=cost_all,
        ndf=ndf_all,
        weight=weight_all,
        valid=valid_all,
        cost_avg="global-mean",
    )

    assert valid_all[0][0]["obs"] == 1
    assert ndf_all[0][0]["obs"] == 2
    assert np.isposinf(cost_all[0][0]["obs"])
    assert np.isposinf(combined)


def test_loss_validity_null_chi2red(tmp_path):
    cdir = tmp_path / "graniitti"
    cdir.mkdir()

    results = {
        "datasets": [{"type": "test", "sets": [{"name": "ZERO_NDF"}]}],
        "mc": [
            [
                {
                    "obs": _hist_entry(_hobj([0.0, 0.0], [0.0, 0.0]), "MC"),
                }
            ]
        ],
        "data": [
            [
                {
                    "obs": {**_hist_entry(_hobj([0.0, 0.0], [0.0, 0.0]), "Data"), "fitw": 1.0},
                }
            ]
        ],
    }

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        loss.visualize_losses(
            results=results,
            run_name="validity_mask_test",
            cdir=str(cdir),
            tunename="TUNE_TEST",
            valid_arr=[[{"obs": 0}]],
            summary={"cost": "chi2", "metrics": {"chi2": 1.0}, "tunename": "TUNE_TEST"},
        )

    with open(
        cdir / "figs" / "icetune" / "validity_mask_test" / "summary.json", encoding="utf-8"
    ) as f:
        summary = json.load(f)

    label = "test__ZERO_NDF"
    assert summary["valid"][label]["obs"] == 0
    assert summary["nbin"][label]["obs"] == pytest.approx(0.0)
    assert summary["chi2red"][label]["obs"] is None
    assert isinstance(summary["created_at_unix"], float)
    assert isinstance(summary["created_at_datetime"], str)


# Plot real correlated likelihoods without replacing marginal contributions or fractional ranks
def test_joint_likelihood_plot(tmp_path):
    data, mc = {}, {}
    for name, count in [('first', 2.), ('second', 3.)]:
        data[name] = {**_hist_entry(_hobj([count], [1.]), 'Data'), 'fitw': 1.}
        mc[name] = _hist_entry(_hobj([count + 1], [0.]), 'MC')
    results = {'datasets': [{'type': 'test', 'sets': [{'name': 'FIDUCIAL', 'title': 'Fiducial muon cuts'}]}],
               'mc': [[mc]], 'data': [[data]]}
    covariance = {'layout': cov.build_layout(results['data']), 'data_total_covariance': np.ones((2, 2))}
    bundle = icecost.evaluate_cost_bundle(
        results=results, selected_cost='chi2', cost_rho='quadratic', cost_avg='global-mean',
        covariance_payload=covariance, rngseed=1, wasserstein_cache=None)
    payload = graniitti_driver.GraniittiDriver().render_trial_figures_to_dir(
        outputs={'results': results, **bundle, 'tunename': 'TUNE_FULL'},
        param={'run_name': 'joint_plot_test', 'cdir': str(tmp_path)},
        summary_payload={'cost': 'chi2', 'metrics': bundle['metrics']},
        output_dir=str(tmp_path / 'figures'), summary_file=str(tmp_path / 'summary.json'))
    # The supported residual [1, 1] has total chi-square 1 and rank 1, shared equally between bins
    for name in data:
        assert payload['chi2']['test__FIDUCIAL'][name] == pytest.approx(.5)
        assert payload['nbin']['test__FIDUCIAL'][name] == pytest.approx(.5)
        assert payload['chi2red']['test__FIDUCIAL'][name] == pytest.approx(1.)
        assert payload['valid']['test__FIDUCIAL'][name] == 1
    assert list((tmp_path / 'figures').rglob('*.pdf'))


# Check that a publish transaction with only a summary is rejected
def test_missing_rendered_figures(tmp_path):
    cdir = tmp_path / "graniitti"
    cdir.mkdir()
    run_name = "publish_summary_only_test"

    # Write only the summary file without any rendered outputs
    def render_summary_only(output_dir, summary_file):
        Path(summary_file).write_text(
            json.dumps({"cost": "chi2", "metrics": {"chi2": 3.0}, "trial_id": "trial-empty"}),
            encoding="utf-8",
        )
        return {"cost": "chi2", "metrics": {"chi2": 3.0}, "trial_id": "trial-empty"}

    with pytest.raises(RuntimeError, match="produced no figure outputs"):
        icetune._publish_with_lock(
            cdir=str(cdir),
            run_name=run_name,
            render_fn=render_summary_only,
        )

    assert not Path(icetune.icetune_figure_dir(cdir=str(cdir), run_name=run_name)).exists()
    assert not Path(
        icetune.icetune_figure_publish_tmp_dir(cdir=str(cdir), run_name=run_name)
    ).exists()


# Reject an impossible observable even when another observable fits exactly
@pytest.mark.parametrize("fitw", [0.0, 1.0])
def test_cost_bundle_rejects_zero_error_disagreement(fitw):
    results = {
        "mc": [[{"good": {"hdata": _hobj([1], [1])}, "bad": {"hdata": _hobj([1], [0])}}]],
        "data": [[{"good": {"hdata": _hobj([1], [1]), "fitw": 1.0},
                  "bad": {"hdata": _hobj([2], [0]), "fitw": fitw}}]],
    }
    bundle = icecost.evaluate_cost_bundle(
        results=results, selected_cost="chi2", cost_rho="quadratic", cost_avg="global-mean",
        covariance_payload=None, rngseed=0, wasserstein_cache=None,
    )
    assert bundle["likelihood"]["valid"] is (fitw <= 0.0)
    if fitw > 0.0:
        assert np.isinf(bundle["metrics"]["chi2"])
        assert bundle["likelihood"]["objective"]["value"] is None
    else:
        assert bundle["metrics"]["chi2"] == pytest.approx(0.0)


# Integrate physical densities independently of subdivision and unequal bin widths
@pytest.mark.parametrize("bins", [[0, 1, 2], [0, 0.5, 1, 1.5, 2], [0, 0.2, 1, 1.3, 2]])
def test_wasserstein_integrates_bin_density(bins):
    bins = np.asarray(bins, dtype=float)
    centers = 0.5 * (bins[:-1] + bins[1:])
    first = (centers < 1.0).astype(float)
    second = (centers > 1.0).astype(float)
    mc = hist.hobj(counts=first, errs=np.zeros_like(first), bins=bins, cbins=centers)
    data = hist.hobj(counts=second, errs=np.zeros_like(second), bins=bins, cbins=centers)
    value, _ = objective.wasserstein_cost(mc, data, N=1)
    assert value == pytest.approx(1.0)


# Integrate a CDF difference that changes sign inside a bin
@pytest.mark.parametrize("bins", [[0, 1, 3, 4], [0, 1, 2, 3, 4]])
def test_wasserstein_integrates_cdf_zero_crossing(bins):
    bins = np.asarray(bins, dtype=float)
    centers = 0.5 * (bins[:-1] + bins[1:])
    first = ((centers < 1.0) | (centers > 3.0)).astype(float)
    second = ((centers > 1.0) & (centers < 3.0)).astype(float)
    mc = hist.hobj(counts=first, errs=np.zeros_like(first), bins=bins, cbins=centers)
    data = hist.hobj(counts=second, errs=np.zeros_like(second), bins=bins, cbins=centers)
    value, _ = objective.wasserstein_cost(mc, data, N=1)
    assert value == pytest.approx(2.0)

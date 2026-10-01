# Histogram objectives for plotting, generator fits and differentiable amplitude fits
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import torch

import numpy as np
from core.numerics import array
from core.stats import cov, hist, uncertainty
from core.tune import likelihood as icetune_likelihood


# Select histogram keys after verifying that MC and data observables match
def get_histogram_keys(mc_subset: dict, data_subset: dict, context: str = "") -> list:
    mc_keys = set(mc_subset.keys())
    data_keys = set(data_subset.keys())

    if mc_keys != data_keys:
        missing_in_mc = sorted(data_keys - mc_keys)
        missing_in_data = sorted(mc_keys - data_keys)
        raise Exception(
            __name__
            + f".get_histogram_keys: MC/Data histogram key mismatch{context} | missing_in_mc={missing_in_mc} | missing_in_data={missing_in_data}"
        )

    if len(mc_keys) == 0:
        raise Exception(__name__ + f".get_histogram_keys: Empty histogram key set{context}")

    return sorted(mc_keys)


# Compute one deterministic per-observable Wasserstein seed
def wasserstein_seed(*, rngseed: int, dataset: int, subset: int, observable_index: int):
    return np.random.SeedSequence([int(rngseed) % (2**32), int(dataset), int(subset), int(observable_index)])


# Precompute fixed data replicas and MC standard-normal draws for one tune
def build_wasserstein_replica_cache(
    *, data: list, covariance_payload: dict | None, rngseed: int, sample_count: int = 100
) -> dict:
    cache = {}
    for dataset, dataset_data in enumerate(data):
        for subset, subset_data in enumerate(dataset_data):
            for observable_index, observable in enumerate(sorted(subset_data)):
                h_data = subset_data[observable]["hdata"]
                counts = np.asarray(h_data.counts_scaled, dtype=float)
                errors = np.asarray(h_data.errs_scaled, dtype=float)
                active = (
                    np.asarray(h_data.valid, dtype=bool) & np.isfinite(counts) & np.isfinite(errors) & (errors >= 0.0)
                )
                mean = np.where(active, counts, 0.0)
                if covariance_payload is None:
                    covariance = np.diag(np.where(active, errors**2, 0.0))
                else:
                    covariance = cov.observable_covariance_block(
                        payload=covariance_payload,
                        dataset=dataset,
                        subset=subset,
                        observable=observable,
                        bin_count=len(counts),
                    )
                    covariance[~active, :] = 0.0
                    covariance[:, ~active] = 0.0
                cache[(dataset, subset, observable)] = build_wasserstein_replica_state(
                    data_mean=mean,
                    data_total_covariance=covariance,
                    sample_count=sample_count,
                    rngseed=wasserstein_seed(
                        rngseed=rngseed, dataset=dataset, subset=subset, observable_index=observable_index
                    ),
                )
    return cache


# Evaluate one requested histogram cost over all observables
def compute_costs(
    results: dict,
    cost_func: str,
    cost_rho: str,
    covariance_payload: dict | None = None,
    rngseed: int = 0,
    wasserstein_cache: dict | None = None,
) -> tuple[list, list, list, list]:
    mc = results["mc"]
    data = results["data"]

    output = {name: [] for name in ("cost", "ndf", "weight", "valid")}

    for i in range(len(mc)):
        dataset_output = {name: [] for name in output}

        for j in range(len(mc[i])):
            subset_output = {name: {} for name in output}

            hist_keys = get_histogram_keys(
                mc_subset=mc[i][j], data_subset=data[i][j], context=f" for dataset={i}, subset={j}"
            )

            for observable_index, obs_key in enumerate(hist_keys):
                if cost_func in {"chi2", "gaussian"}:
                    h_mc = mc[i][j][obs_key]["hdata"]
                    h_data = data[i][j][obs_key]["hdata"]
                    if covariance_payload is None:
                        c, n = (gaussian_cost(h_mc, h_data) if cost_func == "gaussian" else
                                chi2_cost(h_mc=h_mc, h_data=h_data, rho=cost_rho))
                    else:
                        bins, covariance = cov.observable_covariance_selection(
                            payload=covariance_payload,
                            dataset=i,
                            subset=j,
                            observable=obs_key,
                            bin_count=len(h_mc.counts_scaled),
                        )
                        mc_valid = np.asarray(h_mc.valid, dtype=bool)[bins]
                        bins = bins[mc_valid]
                        covariance = covariance[np.ix_(mc_valid, mc_valid)]
                        result = correlated_chi2_arrays(
                            mc_prediction=h_mc.counts_scaled[bins],
                            data_values=np.asarray(h_data.counts_scaled, dtype=float)[bins],
                            mc_stat_uncertainty=h_mc.errs_scaled[bins],
                            fit_weights=np.ones(len(bins), dtype=float),
                            data_total_covariance=covariance,
                            mc_covariance=h_mc.covariance_scaled[bins][:, bins],
                            gaussian=cost_func == "gaussian",
                        )
                        c = result[cost_func] if result["valid"] else np.inf
                        n = result["ndf"]
                elif cost_func == "ratio2":
                    c, n = ratio2_cost(
                        h_mc=mc[i][j][obs_key]["hdata"], h_data=data[i][j][obs_key]["hdata"], rho=cost_rho
                    )
                elif cost_func == "wasserstein":
                    data_total_covariance = None
                    replica_state = None
                    if wasserstein_cache is not None:
                        replica_state = wasserstein_cache[(i, j, obs_key)]
                    elif covariance_payload is not None:
                        bin_count = len(mc[i][j][obs_key]["hdata"].counts_scaled)
                        data_total_covariance = cov.observable_covariance_block(
                            payload=covariance_payload, dataset=i, subset=j, observable=obs_key, bin_count=bin_count
                        )
                    c, n = wasserstein_cost(
                        h_mc=mc[i][j][obs_key]["hdata"],
                        h_data=data[i][j][obs_key]["hdata"],
                        rngseed=wasserstein_seed(
                            rngseed=rngseed, dataset=i, subset=j, observable_index=observable_index
                        ),
                        data_total_covariance=data_total_covariance,
                        replica_state=replica_state,
                    )
                else:
                    raise Exception(__name__ + f".compute_costs: Unknown cost_func = {cost_func}")

                fitw = data[i][j][obs_key]["fitw"]
                # Keep positive-infinite costs active so impossible disagreements reject the trial
                value = array.to_numpy(c)
                cost_is_valid = np.isfinite(value) or np.isposinf(value)
                is_valid = int(cost_is_valid and np.isfinite(n) and np.isfinite(fitw) and (n > 0) and (fitw > 0))

                for name, value in zip(output, (c, n, fitw, is_valid), strict=True):
                    subset_output[name][obs_key] = value

            for name in output:
                dataset_output[name].append(subset_output[name])
        for name in output:
            output[name].append(dataset_output[name])

    return tuple(output.values())


# Convert additive joint contributions into the nested observable cost structure
def joint_observable_cost_arrays(results: dict, joint: dict) -> tuple[list, list, list, list]:
    records = {
        (int(item["dataset"]), int(item["subset"]), str(item["observable"])): item
        for item in joint.get("observables", [])
    }
    output = {name: [] for name in ("cost", "ndf", "weight", "valid")}
    for dataset, dataset_data in enumerate(results["data"]):
        dataset_output = {name: [] for name in output}
        for subset, subset_data in enumerate(dataset_data):
            subset_output = {name: {} for name in output}
            for observable in sorted(subset_data):
                record = records.get((dataset, subset, observable))
                data_fit_weight = float(subset_data[observable]["fitw"])
                record_fit_weight = None if record is None else record.get("weight")
                fit_weight = float(record_fit_weight) if record_fit_weight is not None else data_fit_weight
                valid = bool(record is not None and record.get("valid", False) and fit_weight > 0.0)
                objective = None if record is None else record.get("objective_value")
                ndf = None if record is None else record.get("ndf")
                cost = objective / fit_weight if valid and objective is not None else np.nan
                values = (cost, 0.0 if ndf is None else float(ndf), fit_weight, int(valid))
                for name, value in zip(output, values, strict=True):
                    subset_output[name][observable] = value
            for name in output:
                dataset_output[name].append(subset_output[name])
        for name in output:
            output[name].append(dataset_output[name])
    return tuple(output.values())


# Combine one nested observable cost tuple using the configured convention
def combined_observable_cost(arrays: tuple[list, list, list, list], cost_avg: str) -> float:
    cost, ndf, weight, valid = arrays
    return icetune_likelihood.combine_costs(cost=cost, ndf=ndf, weight=weight, valid=valid, cost_avg=cost_avg)


# Combine one joint cost and retain robust global mode support
def combined_joint_cost(*, results: dict, joint: dict, layout: list[dict], cost_avg: str) -> float:
    if joint.get("observables"):
        return combined_observable_cost(joint_observable_cost_arrays(results, joint), cost_avg)
    ndf = float(joint.get("ndf", 0.0))
    cost = joint.get("cost")
    if not joint.get("valid", False) or cost is None or ndf <= 0.0:
        return np.inf
    if cost_avg == "dataset-mean":
        raise ValueError("Dataset normalization requires additive observable costs, use quadratic residuals")
    if joint.get("relative_fitw", False):
        return cost / ndf
    if cost_avg.endswith("-unweighted"):
        return cost / ndf
    fitw = []
    for record in layout:
        item = results["data"][record["dataset"]][record["subset"]][record["observable"]]
        weight = float(item["fitw"])
        if len(record["bins"]) > 0 and np.isfinite(weight) and weight > 0.0:
            fitw.append(weight)
    if not fitw:
        return np.inf
    return len(fitw) / float(np.sum(fitw)) * cost / ndf


# Evaluate only the selected optimizer cost and canonical Gaussian likelihood
def evaluate_cost_bundle(
    *,
    results: dict,
    selected_cost: str,
    cost_rho: str,
    cost_avg: str,
    covariance_payload: dict | None,
    rngseed: int,
    wasserstein_cache: dict | None,
) -> dict:
    arrays = {name: {} for name in ("cost", "ndf", "weight", "valid")}
    metrics = {}

    if covariance_payload is not None:
        joint = correlated_chi2(results=results, payload=covariance_payload, gaussian=selected_cost == "gaussian")
        chi2_arrays = joint_observable_cost_arrays(results, joint)
        likelihood = icetune_likelihood.build_joint_likelihood_payload(joint=joint, covariance_mode="full")
        metrics["chi2"] = joint["chi2"] if joint["valid"] else np.inf
    else:
        chi2_arrays = compute_costs(results=results, cost_func="chi2", cost_rho="quadratic")
        likelihood = icetune_likelihood.build_likelihood_payload(
            cost=chi2_arrays[0], ndf=chi2_arrays[1], weight=chi2_arrays[2], valid=chi2_arrays[3], rho="quadratic"
        )
        metrics["chi2"] = combined_observable_cost(chi2_arrays, "sum")
    for name, values in zip(arrays, chi2_arrays, strict=True):
        arrays[name]["chi2"] = values

    if selected_cost == "chi2" and cost_avg == "dataset-mean" and likelihood["valid"]:
        metrics["chi2"] = combined_observable_cost(chi2_arrays, cost_avg)

    if selected_cost == "gaussian":
        if covariance_payload is None:
            selected_arrays = compute_costs(results, "gaussian", "quadratic")
            metrics["gaussian"] = combined_observable_cost(selected_arrays, "sum")
            for name, values in zip(arrays, selected_arrays, strict=True):
                arrays[name]["gaussian"] = values
        else:
            metrics["gaussian"] = joint["gaussian"] if joint["valid"] else np.inf
        value = float(array.to_numpy(metrics["gaussian"]))
        chi2 = likelihood["objective"]["value"]
        likelihood["valid"] = bool(likelihood["valid"] and np.isfinite(value))
        likelihood["kind"] = "gaussian"
        likelihood["objective"].update(name="gaussian", value=value if likelihood["valid"] else None, reduced=None)
        likelihood["components"] = {**likelihood.get("components", {}), "chi2": chi2,
            "logdet": value - chi2 if likelihood["valid"] else None, "observable_diagnostics": "chi2_only"}
        likelihood["mc_uncertainty"]["sigma_objective"] = None
    elif selected_cost != "chi2":
        selected_arrays = compute_costs(
            results=results,
            cost_func=selected_cost,
            cost_rho=cost_rho,
            covariance_payload=covariance_payload,
            rngseed=rngseed,
            wasserstein_cache=wasserstein_cache,
        )
        for name, values in zip(arrays, selected_arrays, strict=True):
            arrays[name][selected_cost] = values
        metrics[selected_cost] = combined_observable_cost(selected_arrays, cost_avg)
        if selected_cost == "ratio2" and covariance_payload is not None:
            ratio_fitw = None
            if cost_avg.endswith("-unweighted"):
                data_fitw = cov.comparison_vectors(results, covariance_payload["layout"])[-1]
                ratio_fitw = np.where(np.isfinite(data_fitw) & (data_fitw > 0.0), 1.0, 0.0)
            ratio = correlated_ratio2(
                results=results, payload=covariance_payload, rho=cost_rho, fitw=ratio_fitw, relative=True
            )
            metrics["ratio2"] = combined_joint_cost(
                results=results, joint=ratio, layout=covariance_payload["layout"], cost_avg=cost_avg
            )

    likelihood["mc_statistics"] = uncertainty.comparison_statistics(results)
    output = {
        "metrics": metrics,
        "cost_arr": arrays["cost"],
        "ndf_arr": arrays["ndf"],
        "weight_arr": arrays["weight"],
        "valid_arr": arrays["valid"],
        "likelihood": likelihood,
    }
    return output


# Evaluate additive diagonal or generalized histogram chi-square terms
def histogram_chi2_contributions_torch(
    *,
    counts: torch.Tensor,
    errors: torch.Tensor,
    data_counts: torch.Tensor,
    data_errors: torch.Tensor,
    mask: torch.Tensor,
    fit_weights: torch.Tensor,
    data_total_covariance: torch.Tensor | None = None,
    covariance_indices: torch.Tensor | None = None,
) -> torch.Tensor:
    import torch

    active = mask & torch.isfinite(fit_weights) & (fit_weights > 0.0)
    if data_total_covariance is None:
        residual = _additive_histogram_residual(
            torch.where(active, counts, 0.0),
            torch.where(active, errors, 0.0),
            torch.where(active, data_counts, 0.0),
            torch.where(active, data_errors, 0.0),
        )
        residual = torch.where(torch.isnan(residual), 0.0, residual)
        return torch.where(active, fit_weights, 0.0) * residual.square()
    if covariance_indices is None:
        raise ValueError("Correlated Torch chi-square requires covariance indices")
    selected = torch.nonzero(active[covariance_indices], as_tuple=True)[0]
    covariance_indices = covariance_indices[selected]
    if not covariance_indices.numel():
        return torch.where(active, counts + errors, 0.0) * 0.0
    covariance = data_total_covariance[selected][:, selected]
    contributions = []
    for prediction, error in zip(counts, errors, strict=True):
        result = correlated_chi2_arrays(
            mc_prediction=prediction[covariance_indices],
            data_values=data_counts[covariance_indices],
            mc_stat_uncertainty=error[covariance_indices],
            fit_weights=fit_weights[covariance_indices],
            data_total_covariance=covariance,
        )
        terms = result.get("contributions") if result["valid"] else None
        if terms is None:
            terms = prediction.new_full((len(covariance_indices),), torch.inf)
        contributions.append(torch.zeros_like(prediction).index_copy(0, covariance_indices, terms))
    return torch.stack(contributions)


# Evaluate diagonal or correlated histogram chi-square with Torch autograd
def histogram_chi2_torch(
    *,
    counts: torch.Tensor,
    errors: torch.Tensor,
    data_counts: torch.Tensor,
    data_errors: torch.Tensor,
    mask: torch.Tensor,
    fit_weights: torch.Tensor,
    data_total_covariance: torch.Tensor | None = None,
    covariance_indices: torch.Tensor | None = None,
) -> torch.Tensor:
    contribution = histogram_chi2_contributions_torch(
        counts=counts,
        errors=errors,
        data_counts=data_counts,
        data_errors=data_errors,
        mask=mask,
        fit_weights=fit_weights,
        data_total_covariance=data_total_covariance,
        covariance_indices=covariance_indices,
    )
    return contribution.sum(dim=1)


# Compute compatible histogram arrays and one shared statistical validity mask
def _cost_bin_inputs(h_mc: hist.hobj, h_data: hist.hobj):
    h_mc.check_compatible_bins(other=h_data, operation="histogram cost")

    counts_mc = array.asarray(h_mc.counts_scaled, dtype=float)
    xp = array.namespace(counts_mc)
    errs_mc = array.asarray(h_mc.errs_scaled, like=counts_mc, dtype=float)
    counts_data = array.asarray(h_data.counts_scaled, like=counts_mc, dtype=float)
    errs_data = array.asarray(h_data.errs_scaled, like=counts_mc, dtype=float)
    arrays = (counts_mc, errs_mc, counts_data, errs_data)

    if any(values.ndim != 1 for values in arrays):
        raise ValueError(__name__ + ".cost: histogram counts and errors must be one-dimensional")
    if any(values.shape != counts_mc.shape for values in arrays[1:]):
        raise ValueError(__name__ + ".cost: histogram counts and errors must have identical shapes")

    bins = np.asarray(h_mc.bins, dtype=float)
    if bins.ndim != 1 or len(bins) != len(counts_mc) + 1:
        raise ValueError(__name__ + ".cost: histogram bin edges must have length nbin + 1")
    if not np.all(np.isfinite(bins)) or np.any(np.diff(bins) <= 0.0):
        raise ValueError(__name__ + ".cost: histogram bin edges must be finite and strictly increasing")

    valid_mc = array.asarray(h_mc.valid, like=counts_mc, dtype=bool)
    valid_data = array.asarray(h_data.valid, like=counts_mc, dtype=bool)
    if valid_mc.shape != counts_mc.shape or valid_data.shape != counts_mc.shape:
        raise ValueError(__name__ + ".cost: histogram validity masks must match the bin counts")

    active = valid_mc & valid_data
    if (
        any(not xp.all(xp.isfinite(values[active])) for values in arrays)
        or xp.any(errs_mc[active] < 0.0)
        or xp.any(errs_data[active] < 0.0)
    ):
        return None

    return counts_mc, errs_mc, counts_data, errs_data, active



# Compute standardized residuals and mark bins without variance as unavailable
def _standardized_residual(delta, sigma):
    xp = array.namespace(delta)
    uncertain = sigma > 0.0
    ratio = delta / xp.where(uncertain, sigma, 1.0)
    certain = xp.where(xp.abs(delta) <= np.finfo(float).tiny, xp.full_like(delta, xp.nan), xp.full_like(delta, xp.inf))
    return xp.where(uncertain, ratio, certain)



# Build deterministic Wasserstein replicas from one fixed data covariance
def build_wasserstein_replica_state(
    *,
    data_mean: np.ndarray,
    data_total_covariance: np.ndarray,
    sample_count: int = 100,
    rngseed: int | np.random.SeedSequence = 0,
) -> dict:
    data_mean = np.asarray(data_mean, dtype=float)
    data_total_covariance = np.asarray(data_total_covariance, dtype=float)
    if data_mean.ndim != 1:
        raise ValueError(__name__ + ".build_wasserstein_replica_state: mean must be a vector")
    if data_total_covariance.shape != (len(data_mean), len(data_mean)):
        raise ValueError(__name__ + ".build_wasserstein_replica_state: covariance has incompatible dimensions")
    if not isinstance(sample_count, (int, np.integer)) or sample_count < 1:
        raise ValueError(__name__ + ".build_wasserstein_replica_state: sample_count must be positive")

    rng = np.random.default_rng(rngseed)
    mc_standard = rng.standard_normal((sample_count, len(data_mean)))
    data_standard = rng.standard_normal((sample_count, len(data_mean)))
    factor = uncertainty.covariance_square_root(data_total_covariance)
    data_replicas = np.maximum(data_mean + data_standard @ factor.T, 0.0)
    mc_standard.setflags(write=False)
    data_replicas.setflags(write=False)
    return {
        "sample_count": int(sample_count),
        "bin_count": int(len(data_mean)),
        "mc_standard": mc_standard,
        "data_replicas": data_replicas,
    }



# Validate one reusable Wasserstein common-random-number state
def _validated_wasserstein_replica_state(state: dict, *, sample_count: int, bin_count: int) -> dict:
    if not isinstance(state, dict):
        raise TypeError(__name__ + ".wasserstein_cost: replica_state must be a dictionary")
    if int(state.get("sample_count", -1)) != int(sample_count):
        raise ValueError(__name__ + ".wasserstein_cost: replica sample count does not match N")
    if int(state.get("bin_count", -1)) != int(bin_count):
        raise ValueError(__name__ + ".wasserstein_cost: replica bin count is incompatible")
    mc_standard = np.asarray(state.get("mc_standard"), dtype=float)
    data_replicas = np.asarray(state.get("data_replicas"), dtype=float)
    expected = (int(sample_count), int(bin_count))
    if mc_standard.shape != expected or data_replicas.shape != expected:
        raise ValueError(__name__ + ".wasserstein_cost: malformed replica state arrays")
    return {"mc_standard": mc_standard, "data_replicas": data_replicas}



# Compute the unbalanced Wasserstein-1 histogram cost
def wasserstein_cost(
    h_mc: hist.hobj,
    h_data: hist.hobj,
    N=100,
    *,
    rngseed: int | np.random.SeedSequence = 0,
    data_total_covariance: np.ndarray | None = None,
    replica_state: dict | None = None,
):
    inputs = _cost_bin_inputs(h_mc=h_mc, h_data=h_data)
    if inputs is None:
        return np.nan, 0

    counts_mc, errs_mc, counts_data, errs_data, active = inputs
    occupied = (
        (np.abs(counts_mc) > np.finfo(float).tiny)
        | (errs_mc > 0.0)
        | (np.abs(counts_data) > np.finfo(float).tiny)
        | (errs_data > 0.0)
    )
    nbin = int(np.sum(active & occupied))
    if nbin == 0:
        return 0.0, 0
    if not isinstance(N, (int, np.integer)) or N < 1:
        raise ValueError(__name__ + ".wasserstein_cost: N must be a positive integer")
    if np.any(counts_mc[active] < 0.0) or np.any(counts_data[active] < 0.0):
        return np.inf, nbin

    counts_mc = np.where(active, counts_mc, 0.0)
    errs_mc = np.where(active, errs_mc, 0.0)
    counts_data = np.where(active, counts_data, 0.0)
    errs_data = np.where(active, errs_data, 0.0)
    dx = hist.bins2binwidth(h_mc.bins)
    if replica_state is None:
        if data_total_covariance is None:
            data_total_covariance = np.diag(errs_data**2)
        else:
            data_total_covariance = np.asarray(data_total_covariance, dtype=float)
            if data_total_covariance.shape != (len(counts_data), len(counts_data)):
                raise ValueError(__name__ + ".wasserstein_cost: data covariance has incompatible dimensions")
            data_total_covariance = data_total_covariance.copy()
            data_total_covariance[~active, :] = 0.0
            data_total_covariance[:, ~active] = 0.0
        replica_state = build_wasserstein_replica_state(
            data_mean=counts_data, data_total_covariance=data_total_covariance, sample_count=N, rngseed=rngseed
        )
    state = _validated_wasserstein_replica_state(replica_state, sample_count=N, bin_count=len(counts_data))
    c_mc = np.maximum(counts_mc + state["mc_standard"] * errs_mc, 0.0)
    c_data = np.where(active[None, :], state["data_replicas"], 0.0)

    # Integrate differential bin heights to cumulative cross sections
    right = np.cumsum((c_mc - c_data) * dx, axis=1)
    left = np.concatenate((np.zeros_like(right[:, :1]), right[:, :-1]), axis=1)

    # Integrate the absolute linear CDF difference within each bin
    total = np.abs(left) + np.abs(right)
    crossing = (left < 0.0) & (right > 0.0) | (left > 0.0) & (right < 0.0)
    area = 0.5 * total
    np.divide(left**2 + right**2, 2.0 * total, out=area, where=crossing)
    W1_mean = np.mean(np.sum(dx * area, axis=1))

    return W1_mean, nbin



# Compute a robust or quadratic residual loss
def rho_wrapper(r: np.ndarray, rho: str = "quadratic"):
    """Residual rho function"""

    xp = array.namespace(r)

    # Quadratic
    if rho == "quadratic":
        return r**2
    # Absolute
    if rho == "absolute":
        return xp.abs(r)
    # Cauchy type
    if rho == "cauchy":
        return xp.log(1 + r**2)
    # Huber type
    if rho == "huber":
        abs_r = xp.abs(r)
        return xp.where(abs_r <= 1.0, 0.5 * r**2, abs_r - 0.5)
    raise Exception(__name__ + ".rho_wrapper: Unknown residual rho function mode")



# Reduce one histogram residual definition over informative bins
def _histogram_residual_cost(h_mc: hist.hobj, h_data: hist.hobj, rho: str, residual_fn):
    inputs = _cost_bin_inputs(h_mc=h_mc, h_data=h_data)
    if inputs is None:
        return np.nan, 0
    counts_mc, errs_mc, counts_data, errs_data, informative = inputs
    xp = array.namespace(counts_mc)
    nbin = int(informative.sum())
    if nbin == 0:
        return 0.0, 0
    selected = tuple(values[informative] for values in inputs[:4])
    residual = residual_fn(*selected)
    if residual is None:
        return np.nan, 0
    available = ~xp.isnan(residual)
    nbin = int(available.sum())
    if nbin == 0:
        return 0.0, 0
    if xp.any(xp.isinf(residual[available])):
        return np.inf, nbin
    return rho_wrapper(residual[available], rho=rho).sum(), nbin



# Compute additive residuals with MC and data uncertainties in quadrature
def _additive_histogram_residual(counts_mc, errs_mc, counts_data, errs_data):
    return _standardized_residual(delta=counts_mc - counts_data, sigma=array.hypot(errs_mc, errs_data))



# Compute symmetric asinh coordinates and analytic Jacobians
def asinh_ratio_coordinates(counts_mc, errs_mc, counts_data, errs_data):
    # [REFERENCE: Lupton, Gunn and Szalay, AJ 118 (1999) 1406, arXiv:astro-ph/9903081]
    xp = array.namespace(counts_mc)
    soft = array.hypot(errs_mc, errs_data)
    uncertain = soft > 0.0
    scale = xp.where(uncertain, soft, 1.0)
    residual = xp.arcsinh(counts_mc / scale) - xp.arcsinh(counts_data / scale)
    residual = xp.where(uncertain, residual, counts_mc - counts_data)
    mc_scale = array.hypot(counts_mc, soft)
    data_scale = array.hypot(counts_data, soft)
    mc_jacobian = xp.where(uncertain, 1.0 / xp.where(uncertain, mc_scale, 1.0), 0.0)
    data_jacobian = xp.where(uncertain, -1.0 / xp.where(uncertain, data_scale, 1.0), 0.0)
    return residual, mc_jacobian, data_jacobian



# Compute symmetric multiplicative residuals with propagated asinh uncertainty
def _multiplicative_histogram_residual(counts_mc, errs_mc, counts_data, errs_data):
    residual, mc_jacobian, data_jacobian = asinh_ratio_coordinates(counts_mc, errs_mc, counts_data, errs_data)
    sigma = array.hypot(mc_jacobian * errs_mc, data_jacobian * errs_data)
    return _standardized_residual(delta=residual, sigma=sigma)



# Compute the additive uncertainty weighted histogram cost
def chi2_cost(h_mc: hist.hobj, h_data: hist.hobj, rho: str = "quadratic"):
    return _histogram_residual_cost(h_mc, h_data, rho, _additive_histogram_residual)



# Compute the Gaussian 2NLL without its constant normalization on independent bins
def gaussian_cost(h_mc: hist.hobj, h_data: hist.hobj):
    inputs = _cost_bin_inputs(h_mc, h_data)
    if inputs is None:
        return np.inf, int(np.count_nonzero(h_data.valid))
    mc, mc_error, data, data_error, active = inputs
    xp = array.namespace(mc)
    residual = mc[active] - data[active]
    sigma = array.hypot(mc_error[active], data_error[active])
    positive = sigma > 0.0
    if xp.any((~positive) & (xp.abs(residual) > 0.0)):
        return np.inf, int(active.sum())
    value = ((residual[positive] / sigma[positive])**2 + 2.0 * xp.log(sigma[positive])).sum()
    return value, int(positive.sum())


# Compute a differentiable determinant on the retained covariance subspace
def covariance_logdet(covariance, basis):
    xp = array.namespace(covariance)
    reduced = covariance if basis.shape[1] == len(covariance) else basis.T @ covariance @ basis
    return xp.linalg.slogdet(reduced)[1]


# Compute the multiplicative uncertainty weighted histogram cost
def ratio2_cost(h_mc: hist.hobj, h_data: hist.hobj, rho: str = "quadratic"):
    return _histogram_residual_cost(h_mc, h_data, rho, _multiplicative_histogram_residual)



# Evaluate one residual vector in covariance eigenmodes
def covariance_cost(
    *,
    residual: np.ndarray,
    covariance: np.ndarray,
    fit_weights: np.ndarray,
    rho: str = "quadratic",
    mc_stat_uncertainty: np.ndarray | None = None,
    mc_covariance=None,
    return_contributions: bool = False,
    gaussian: bool = False,
) -> dict:
    residual = array.asarray(residual, dtype=float)
    covariance = array.asarray(covariance, like=residual)
    fit_weights = array.asarray(fit_weights, like=residual)
    xp = array.namespace(residual)
    original = residual
    if covariance.shape != (len(residual), len(residual)) or fit_weights.shape != residual.shape:
        raise ValueError("Residual covariance or fit weights have incompatible dimensions")
    if mc_stat_uncertainty is not None:
        mc_stat_uncertainty = array.asarray(mc_stat_uncertainty, like=residual)
        if mc_stat_uncertainty.shape != residual.shape:
            raise ValueError("Residual MC errors have incompatible dimensions")

    diagonal = xp.diag(covariance)
    selected = xp.isfinite(residual) & xp.isfinite(diagonal) & xp.isfinite(fit_weights) & (fit_weights > 0.0)
    certain = selected & (diagonal <= 0.0)
    if xp.any(certain & (xp.abs(residual) > np.finfo(float).tiny)):
        return {"valid": False, "cost": None, "ndf": int(certain.sum()), "rank": 0}
    active = selected & (diagonal > 0.0)
    residual = residual[active] * xp.sqrt(fit_weights[active])
    covariance = covariance[active][:, active]
    mc_error = None if mc_stat_uncertainty is None else mc_stat_uncertainty[active]
    if not len(residual) or not xp.all(xp.isfinite(covariance)):
        return {"valid": False, "cost": None, "ndf": 0, "rank": 0}

    # Classify retained modes without differentiating eigenvectors at repeated eigenvalues
    support = covariance if xp is np else covariance.detach()
    eigenvalues, eigenvectors = xp.linalg.eigh(support)
    maximum = float(xp.abs(eigenvalues).max())
    threshold = max(maximum * 1.0e-12, 1.0e-30)
    kept = eigenvalues > threshold
    rank = int(kept.sum())
    if float(eigenvalues[0]) < -threshold or rank == 0:
        return {"valid": False, "cost": None, "ndf": 0, "rank": 0}
    discarded = eigenvectors[:, ~kept].T @ residual
    scale = max(float(array.to_numpy(xp.linalg.norm(residual))), np.finfo(float).tiny)
    if float(array.to_numpy(xp.linalg.norm(discarded))) > 1.0e-12 * scale:
        return {"valid": False, "cost": None, "ndf": rank, "rank": rank}

    logdet = None
    if gaussian:
        basis = eigenvectors[:, kept]
        logdet = covariance_logdet(covariance, basis)
        reduced = covariance if rank == len(residual) else basis.T @ covariance @ basis
        rhs = residual if rank == len(residual) else basis.T @ residual
        solved = xp.linalg.solve(reduced, rhs)
        inverse_residual = solved if rank == len(residual) else basis @ solved
    elif rho == "quadratic":
        # Solve full rank covariance directly and retain the pseudoinverse for singular modes
        inverse_residual = (xp.linalg.solve(covariance, residual) if rank == len(residual) else
                            xp.linalg.pinv(covariance, hermitian=True, rcond=threshold / maximum) @ residual)
    else:
        if xp is not np:
            raise ValueError("Autograd with full covariance requires quadratic residuals")
        standardized = (eigenvectors[:, kept].T @ residual) / np.sqrt(eigenvalues[kept])
        mode_cost = rho_wrapper(standardized, rho=rho)
        nonzero = standardized != 0.0
        allocation = np.where(nonzero, mode_cost / np.where(nonzero, standardized, 1.0), 0.0)
        inverse_residual = eigenvectors[:, kept] @ (allocation / np.sqrt(eigenvalues[kept]))
    local_cost = residual * inverse_residual
    cost = local_cost.sum()
    contributions = ndf_contributions = sigma_cost = None
    if return_contributions:
        contributions = xp.full_like(original, xp.nan)
        contributions[active] = local_cost
        ndf_contributions = xp.zeros_like(original)
        ndf_contributions[active] = (eigenvectors[:, kept] ** 2).sum(1)
    if rho == "quadratic" and mc_error is not None:
        gradient = xp.sqrt(fit_weights[active]) * inverse_residual
        sigma_cost = 2.0 * (xp.linalg.norm(gradient * mc_error) if mc_covariance is None else
                            array.sqrt(gradient @ mc_covariance[active][:, active] @ gradient))
    return {"valid": bool(xp.isfinite(cost)), "cost": cost, "sigma_cost": sigma_cost, "logdet": logdet,
            "ndf": rank, "rank": rank, "bins": len(residual), "contributions": contributions,
            "ndf_contributions": ndf_contributions}



# Evaluate one quadratic residual using exact covariance blocks and collective sources
def structured_covariance_cost(
    *,
    residual: np.ndarray,
    structure: dict,
    fit_weights: np.ndarray,
    covariance_transform: np.ndarray,
    diagonal_uncertainty: np.ndarray,
    mc_stat_uncertainty: np.ndarray | None = None,
    return_contributions: bool = False,
    gaussian: bool = False,
) -> dict | None:
    residual = array.asarray(residual, dtype=float)
    xp = array.namespace(residual)
    fit_weights = array.asarray(fit_weights, like=residual)
    transform = array.asarray(covariance_transform, like=residual)
    diagonal_uncertainty = array.asarray(diagonal_uncertainty, like=residual)
    size = int(structure["size"])
    if any(value.shape != (size,) for value in (residual, fit_weights, transform, diagonal_uncertainty)):
        raise ValueError("Structured covariance vectors have incompatible dimensions")
    if mc_stat_uncertainty is not None and mc_stat_uncertainty.shape != (size,):
        raise ValueError("Structured covariance MC errors have incompatible dimensions")

    vectors = array.asarray(structure["collective_vectors"], like=residual)
    covariance_diagonal = array.asarray(structure["diagonal"], like=residual)
    total_diagonal = transform**2 * (covariance_diagonal + xp.sum(vectors**2, axis=1)) + diagonal_uncertainty**2
    selected = (
        xp.isfinite(residual)
        & xp.isfinite(total_diagonal)
        & xp.isfinite(fit_weights)
        & (fit_weights > 0.0)
        & xp.isfinite(transform)
        & xp.isfinite(diagonal_uncertainty)
    )
    certain = selected & (total_diagonal <= 0.0)
    if xp.any(certain & ~(xp.abs(residual) <= np.finfo(float).tiny)):
        return {"valid": False, "cost": None, "ndf": int(xp.sum(certain)), "rank": 0}
    active = selected & (total_diagonal > 0.0)
    active_indices = np.flatnonzero(array.to_numpy(active))
    if len(active_indices) == 0:
        return {"valid": False, "cost": None, "ndf": 0, "rank": 0}

    active_position = np.full(size, -1, dtype=np.int64)
    active_position[active_indices] = np.arange(len(active_indices), dtype=np.int64)
    weighted_residual = residual[active] * xp.sqrt(fit_weights[active])
    active_vectors = vectors[active] * transform[active, None]
    blocks = []
    row_bound = xp.zeros_like(weighted_residual)
    for component in structure["components"]:
        indices = np.asarray(component["indices"], dtype=np.int64)
        retained = array.to_numpy(active)[indices]
        if not np.any(retained):
            continue
        selected = indices[retained]
        positions = active_position[selected]
        block = array.asarray(component["covariance"][np.ix_(retained, retained)], like=residual)
        local_transform = transform[selected]
        block = local_transform[:, None] * block * local_transform[None, :]
        block = block + xp.diag(diagonal_uncertainty[selected] ** 2)
        blocks.append((positions, block))
        row_bound[positions] = xp.sum(xp.abs(block), axis=1)
    if active_vectors.shape[1] > 0:
        row_bound += xp.abs(active_vectors) @ xp.sum(xp.abs(active_vectors), axis=0)
    eigenvalue_bound = float(row_bound.max())

    right_hand_side = xp.column_stack((weighted_residual, active_vectors))
    inverse_rhs = xp.zeros_like(right_hand_side)
    projected_residual = xp.zeros_like(weighted_residual)
    projector_diagonal = xp.zeros_like(weighted_residual)
    rank = 0
    logdet = 0.0 if gaussian else None
    for positions, block in blocks:
        if len(positions) == 1:
            eigenvalue = block[0, 0]
            threshold = max(abs(eigenvalue) * 1.0e-12, 1.0e-30)
            if not xp.isfinite(eigenvalue):
                return None
            if eigenvalue <= threshold:
                if float(xp.sum(active_vectors[positions] ** 2)) > threshold:
                    return None
                continue
            if gaussian:
                logdet = logdet + xp.log(eigenvalue)
            inverse_rhs[positions] = right_hand_side[positions] / eigenvalue
            projected_residual[positions] = weighted_residual[positions]
            projector_diagonal[positions] = 1.0
            rank += 1
            continue
        eigenvalues, eigenvectors = xp.linalg.eigh(block if xp is np else block.detach())
        threshold = max(float(xp.abs(eigenvalues).max()) * 1.0e-12, 1.0e-30)
        if not xp.all(xp.isfinite(eigenvalues)):
            return None
        if float(eigenvalues[0]) < -threshold:
            return None
        kept = eigenvalues > threshold
        if xp.any(~kept) and active_vectors.shape[1] > 0:
            discarded_sources = eigenvectors[:, ~kept].T @ active_vectors[positions]
            if float(xp.linalg.norm(discarded_sources)) ** 2 > threshold:
                return None
        if not xp.any(kept):
            continue
        if gaussian:
            basis = eigenvectors[:, kept]
            logdet = logdet + covariance_logdet(block, basis)
            reduced = basis.T @ block @ basis
            inverse_rhs[positions] = basis @ xp.linalg.solve(reduced, basis.T @ right_hand_side[positions])
        else:
            precision = xp.linalg.pinv(block, hermitian=True, rcond=threshold / float(xp.abs(eigenvalues).max()))
            inverse_rhs[positions] = precision @ right_hand_side[positions]
        residual_projection = eigenvectors[:, kept].T @ weighted_residual[positions]
        projected_residual[positions] = eigenvectors[:, kept] @ residual_projection
        projector_diagonal[positions] = xp.sum(eigenvectors[:, kept] ** 2, axis=1)
        rank += int(xp.sum(kept))

    if rank == 0:
        return {"valid": False, "cost": None, "ndf": 0, "rank": 0}
    residual_scale = max(float(xp.linalg.norm(weighted_residual)), np.finfo(float).tiny)
    if float(xp.linalg.norm(weighted_residual - projected_residual)) > 1.0e-12 * residual_scale:
        return {"valid": False, "cost": None, "ndf": rank, "rank": rank}

    inverse_residual = inverse_rhs[:, 0]
    inverse_vectors = inverse_rhs[:, 1:]
    if active_vectors.shape[1] > 0:
        correction = array.asarray(np.eye(active_vectors.shape[1]), like=residual) + active_vectors.T @ inverse_vectors
        correction = 0.5 * (correction + correction.T)
        correction_eigenvalues, _ = xp.linalg.eigh(correction if xp is np else correction.detach())
        correction_threshold = max(float(xp.abs(correction_eigenvalues).max()) * 1.0e-12, 1.0e-30)
        if not xp.all(xp.isfinite(correction_eigenvalues)) or float(correction_eigenvalues[0]) <= correction_threshold:
            return None
        if gaussian:
            logdet = logdet + xp.linalg.slogdet(correction)[1]
        source_projection = active_vectors.T @ inverse_residual
        source_solution = xp.linalg.solve(correction, source_projection)
        inverse_residual = inverse_residual - inverse_vectors @ source_solution

    reconstructed = xp.zeros_like(weighted_residual)
    for positions, block in blocks:
        reconstructed[positions] = block @ inverse_residual[positions]
    if active_vectors.shape[1] > 0:
        reconstructed += active_vectors @ (active_vectors.T @ inverse_residual)
    residual_error = float(xp.linalg.norm(reconstructed - projected_residual, ord=np.inf))
    residual_scale = max(
        float(xp.linalg.norm(projected_residual, ord=np.inf)),
        eigenvalue_bound * float(xp.linalg.norm(inverse_residual, ord=np.inf)),
        1.0e-30,
    )
    if not np.isfinite(residual_error) or residual_error > 1.0e-9 * residual_scale:
        return None

    cost = weighted_residual @ inverse_residual
    contributions = None
    ndf_contributions = None
    if return_contributions:
        contributions = xp.full_like(residual, xp.nan)
        contributions[active] = weighted_residual * inverse_residual
        ndf_contributions = xp.zeros_like(residual)
        ndf_contributions[active] = projector_diagonal
    sigma_cost = None
    if mc_stat_uncertainty is not None:
        active_mc_error = array.asarray(mc_stat_uncertainty, like=residual)[active]
        gradient_scale = xp.sqrt(fit_weights[active]) * inverse_residual
        sigma_cost = 2.0 * xp.linalg.norm(gradient_scale * active_mc_error)
    return {
        "valid": bool(xp.isfinite(cost)),
        "cost": cost,
        "logdet": logdet,
        "sigma_cost": sigma_cost,
        "ndf": rank,
        "rank": rank,
        "bins": int(len(active_indices)),
        "contributions": contributions,
        "ndf_contributions": ndf_contributions,
    }



# Compute the numerical scale from finite bins with positive fit weight
def _chi2_array_scale(mc, data, errors, weights, covariance, *, diagonal=False) -> tuple[np.ndarray, float]:
    weights = array.asarray(weights, like=mc)
    xp = array.namespace(mc)
    expected_shape = (len(mc),) if diagonal else (len(mc), len(mc))
    if weights.shape != mc.shape or covariance.shape != expected_shape:
        raise ValueError("Correlated chi-square weights or covariance have incompatible dimensions")
    active = xp.isfinite(weights) & (weights > 0.0)
    if not xp.any(active):
        return active, 0.0
    selected_covariance = covariance if xp.all(active) else (
        covariance[active] if covariance.ndim == 1 else covariance[active][:, active])
    if (
        any(not xp.all(xp.isfinite(value[active])) for value in (mc, data, errors))
        or not xp.all(xp.isfinite(selected_covariance))
    ):
        return active, np.nan
    maxima = (xp.abs(mc[active]).max(), xp.abs(data[active]).max(), xp.abs(errors[active]).max(),
              xp.sqrt(xp.abs(selected_covariance).max()))
    scale = max(float(array.to_numpy(value)) for value in maxima)
    return active, scale



# Restore physical covariance units in the Gaussian objective
def gaussian_result(result, scale):
    if result.get("valid", False):
        logdet = result["logdet"] + 2.0 * result["rank"] * np.log(scale)
        result["gaussian"] = result["chi2"] + logdet
        result["valid"] = bool(array.namespace(logdet).isfinite(result["gaussian"]))
    return result


# Evaluate covariance-aware chi-square arrays with a cached structured covariance
def structured_chi2_arrays(
    *,
    mc_prediction: np.ndarray,
    data_values: np.ndarray,
    mc_stat_uncertainty: np.ndarray,
    fit_weights: np.ndarray,
    data_total_covariance: np.ndarray | None = None,
    covariance_structure: dict,
    gaussian: bool = False,
) -> dict | None:
    mc_prediction = array.asarray(mc_prediction, dtype=float)
    xp = array.namespace(mc_prediction)
    data_values = array.asarray(data_values, like=mc_prediction)
    mc_stat_uncertainty = array.asarray(mc_stat_uncertainty, like=mc_prediction)
    if data_total_covariance is None:
        data_total_covariance = (np.asarray(covariance_structure["diagonal"]) +
                                 np.sum(np.asarray(covariance_structure["collective_vectors"]) ** 2, axis=1))
    data_total_covariance = array.asarray(data_total_covariance, like=mc_prediction)
    if mc_prediction.shape != data_values.shape or mc_prediction.shape != mc_stat_uncertainty.shape:
        raise ValueError("Structured chi-square arrays have incompatible dimensions")
    active, scale = _chi2_array_scale(
        mc_prediction, data_values, mc_stat_uncertainty, fit_weights, data_total_covariance,
        diagonal=data_total_covariance.ndim == 1,
    )
    if not np.isfinite(scale) or scale <= 0.0:
        return None
    result = structured_covariance_cost(
        residual=xp.where(active, mc_prediction, 0.0) / scale - xp.where(active, data_values, 0.0) / scale,
        structure=covariance_structure,
        fit_weights=fit_weights,
        covariance_transform=xp.where(active, xp.ones_like(mc_prediction) / scale, 0.0),
        diagonal_uncertainty=xp.where(active, mc_stat_uncertainty, 0.0) / scale,
        mc_stat_uncertainty=xp.where(active, mc_stat_uncertainty, 0.0) / scale,
        return_contributions=True,
        gaussian=gaussian,
    )
    if result is None:
        return None
    result["chi2"] = result.pop("cost")
    result["sigma_chi2"] = result.pop("sigma_cost", None)
    return gaussian_result(result, scale) if gaussian else result



# Evaluate covariance-aware chi-square directly from flat arrays
def correlated_chi2_arrays(
    *,
    mc_prediction: np.ndarray,
    data_values: np.ndarray,
    mc_stat_uncertainty: np.ndarray,
    fit_weights: np.ndarray,
    data_total_covariance: np.ndarray,
    mc_covariance=None,
    gaussian: bool = False,
) -> dict:
    mc_prediction = array.asarray(mc_prediction, dtype=float)
    xp = array.namespace(mc_prediction)
    data_values = array.asarray(data_values, like=mc_prediction)
    mc_stat_uncertainty = array.asarray(mc_stat_uncertainty, like=mc_prediction)
    data_total_covariance = array.asarray(data_total_covariance, like=mc_prediction)
    if mc_prediction.shape != data_values.shape or mc_prediction.shape != mc_stat_uncertainty.shape:
        raise ValueError("Correlated chi-square arrays have incompatible dimensions")
    active, scale = _chi2_array_scale(
        mc_prediction, data_values, mc_stat_uncertainty, fit_weights, data_total_covariance
    )
    if not np.isfinite(scale) or scale <= 0.0:
        return {"valid": False, "chi2": None, "sigma_chi2": None, "ndf": 0, "rank": 0}
    scaled_mc = xp.where(active, mc_prediction, 0.0) / scale
    scaled_data = xp.where(active, data_values, 0.0) / scale
    scaled_mc_error = xp.where(active, mc_stat_uncertainty, 0.0) / scale
    residual_covariance = xp.where(active[:, None] & active[None, :], data_total_covariance, 0.0) / scale / scale
    scaled_mc_covariance = (xp.diag(scaled_mc_error**2) if mc_covariance is None else
                            array.asarray(mc_covariance, like=mc_prediction) / scale / scale)
    if scaled_mc_covariance.shape != residual_covariance.shape:
        raise ValueError("MC covariance has incompatible dimensions")
    residual_covariance = residual_covariance + scaled_mc_covariance
    result = covariance_cost(
        residual=scaled_mc - scaled_data,
        covariance=residual_covariance,
        fit_weights=fit_weights,
        rho="quadratic",
        mc_stat_uncertainty=scaled_mc_error,
        mc_covariance=scaled_mc_covariance,
        return_contributions=True,
        gaussian=gaussian,
    )
    result["chi2"] = result.pop("cost")
    result["sigma_chi2"] = result.pop("sigma_cost", None)
    return gaussian_result(result, scale) if gaussian else result



# Evaluate covariance-aware chi-square for a batch of MC predictions
def correlated_chi2_batch(
    *,
    mc_prediction: np.ndarray,
    data_values: np.ndarray,
    mc_stat_uncertainty: np.ndarray,
    fit_weights: np.ndarray,
    data_total_covariance: np.ndarray,
    covariance_structure: dict | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    mc_prediction = np.atleast_2d(np.asarray(mc_prediction, dtype=float))
    mc_stat_uncertainty = np.atleast_2d(np.asarray(mc_stat_uncertainty, dtype=float))
    if mc_prediction.shape != mc_stat_uncertainty.shape:
        raise ValueError("Batched MC values and errors have incompatible dimensions")
    values = np.full(len(mc_prediction), np.inf, dtype=float)
    errors = np.full(len(mc_prediction), np.inf, dtype=float)
    for index in range(len(mc_prediction)):
        result = None
        if covariance_structure is not None:
            result = structured_chi2_arrays(
                mc_prediction=mc_prediction[index],
                data_values=data_values,
                mc_stat_uncertainty=mc_stat_uncertainty[index],
                fit_weights=fit_weights,
                data_total_covariance=data_total_covariance,
                covariance_structure=covariance_structure,
            )
        if result is None:
            result = correlated_chi2_arrays(
                mc_prediction=mc_prediction[index],
                data_values=data_values,
                mc_stat_uncertainty=mc_stat_uncertainty[index],
                fit_weights=fit_weights,
                data_total_covariance=data_total_covariance,
            )
        if result["valid"]:
            values[index] = float(result["chi2"])
            sigma = result.get("sigma_chi2")
            errors[index] = 0.0 if sigma is None else float(sigma)
    return values, errors



# Partition one joint precision-weighted chi-square into observable records
def correlated_observable_records(*, result: dict, layout: list[dict], fit_weights: np.ndarray) -> list[dict]:
    contributions = result.get("contributions")
    ndf_contributions = result.get("ndf_contributions")
    if contributions is None or ndf_contributions is None:
        return []
    contributions = array.asarray(contributions, dtype=float)
    xp = array.namespace(contributions)
    ndf_contributions = array.to_numpy(ndf_contributions)
    fit_weights = array.to_numpy(fit_weights)
    if contributions.shape != fit_weights.shape or ndf_contributions.shape != fit_weights.shape:
        raise ValueError("Joint observable contribution arrays have incompatible dimensions")

    records = []
    for record in layout:
        indices = cov.record_indices(record)
        active = np.isfinite(array.to_numpy(contributions[indices]))
        valid = bool(result.get("valid", False) and np.any(active))
        objective = xp.sum(contributions[indices][active]) if valid else None
        ndf = float(np.sum(ndf_contributions[indices][active])) if valid else None
        weights = fit_weights[indices]
        weight = float(weights[0]) if len(weights) else None
        records.append(
            {
                "dataset": int(record["dataset"]),
                "subset": int(record["subset"]),
                "observable": str(record["observable"]),
                "ndf": ndf,
                "valid": valid,
                "objective_value": objective,
                "weight": weight,
            }
        )
    return records



# Normalize positive fit weights to unit observable mean
def relative_fitw(fitw: np.ndarray, layout: list[dict]) -> tuple[np.ndarray, float]:
    weights = np.asarray(fitw, dtype=float).copy()
    observable_weights = []
    for record in layout:
        local = weights[cov.record_indices(record)]
        positive = local[np.isfinite(local) & (local > 0.0)]
        if len(positive) > 0:
            observable_weights.append(float(positive[0]))
    if not observable_weights:
        return weights, 1.0
    scale = len(observable_weights) / float(np.sum(observable_weights))
    weights[np.isfinite(weights) & (weights > 0.0)] *= scale
    return weights, scale



# Evaluate the asinh ratio cost with the full data covariance
def correlated_ratio2(
    results: dict, payload: dict, *, rho: str = "quadratic", fitw: np.ndarray | None = None, relative: bool = False
) -> dict:
    mc, data, mc_error, data_error, data_fitw = cov.comparison_vectors(results, payload["layout"])
    xp = array.namespace(mc)
    data, mc_error, data_error = (array.asarray(value, like=mc) for value in (data, mc_error, data_error))
    fit_weights = data_fitw if fitw is None else np.asarray(fitw, dtype=float)
    if fit_weights.shape != mc.shape:
        raise ValueError("Ratio fit weights have incompatible dimensions")
    valid = (
        array.asarray(cov.comparison_mc_valid(results, payload["layout"]), like=mc, dtype=bool)
        & xp.isfinite(mc)
        & xp.isfinite(data)
        & xp.isfinite(mc_error)
        & xp.isfinite(data_error)
        & (mc_error >= 0.0)
        & (data_error >= 0.0)
    )
    fit_weights = np.where(array.to_numpy(valid), fit_weights, 0.0)
    fitw_scale = 1.0
    if relative:
        fit_weights, fitw_scale = relative_fitw(fit_weights, payload["layout"])
    active = array.asarray(np.isfinite(fit_weights) & (fit_weights > 0.0), like=mc, dtype=bool)

    coordinate_mc = xp.where(active, mc, 0.0)
    coordinate_data = xp.where(active, data, 0.0)
    coordinate_mc_error = xp.where(active, mc_error, 0.0)
    coordinate_data_error = xp.where(active, data_error, 0.0)
    residual, mc_jacobian, data_jacobian = asinh_ratio_coordinates(
        coordinate_mc, coordinate_mc_error, coordinate_data, coordinate_data_error
    )
    data_total_covariance = array.asarray(payload["data_total_covariance"], like=mc)
    residual_covariance = data_jacobian[:, None] * data_total_covariance * data_jacobian[None, :]
    mc_covariance = cov.comparison_covariance(results, payload["layout"])
    mc_covariance = xp.diag(mc_error**2) if mc_covariance is None else array.asarray(mc_covariance, like=mc)
    residual_covariance = residual_covariance + mc_jacobian[:, None] * mc_covariance * mc_jacobian[None, :]
    result = covariance_cost(
        residual=residual, covariance=residual_covariance, fit_weights=fit_weights, rho=rho, return_contributions=True
    )
    result["relative_fitw"] = bool(relative)
    result["fitw_scale"] = float(fitw_scale)
    result["observables"] = correlated_observable_records(
        result=result, layout=payload["layout"], fit_weights=fit_weights
    )
    return result



# Evaluate the joint generalized chi-square with the fixed data covariance
def correlated_chi2(results: dict, payload: dict, *, gaussian: bool = False) -> dict:
    mc, data, mc_error, _, fitw = cov.comparison_vectors(results, payload["layout"])
    fitw = np.where(cov.comparison_mc_valid(results, payload["layout"]), fitw, 0.0)
    mc_covariance = cov.comparison_covariance(results, payload["layout"])
    structure = cov.prepare_covariance_structure(payload)
    result = None
    if mc_covariance is None and structure is not None:
        result = structured_chi2_arrays(mc_prediction=mc, data_values=data, mc_stat_uncertainty=mc_error,
                                       fit_weights=fitw, covariance_structure=structure, gaussian=gaussian)
    if result is None:
        result = correlated_chi2_arrays(mc_prediction=mc, data_values=data, mc_stat_uncertainty=mc_error,
                                       fit_weights=fitw, data_total_covariance=payload["data_total_covariance"],
                                       mc_covariance=mc_covariance, gaussian=gaussian)
    result["observables"] = correlated_observable_records(result=result, layout=payload["layout"], fit_weights=fitw)
    return result



# Evaluate the shared independent-bin Gaussian chi2 surface
def global_histogram_chi2(
    Y_data: np.ndarray,
    E_data: np.ndarray,
    Y_hat: np.ndarray,
    E_hat: np.ndarray,
    mask: np.ndarray | None = None,
    fit_weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    Y_data = np.asarray(Y_data, dtype=np.float64)
    E_data = np.asarray(E_data, dtype=np.float64)
    Y_hat = np.asarray(Y_hat, dtype=np.float64)
    E_hat = np.asarray(E_hat, dtype=np.float64)
    active = ((Y_data != 0.0) | (Y_hat != 0.0)) if mask is None else np.asarray(mask, dtype=bool)
    weights = (
        np.ones_like(Y_data, dtype=np.float64) if fit_weights is None else np.asarray(fit_weights, dtype=np.float64)
    )
    if weights.shape != Y_data.shape:
        raise ValueError("Histogram fit weights have incompatible dimensions")

    active = active & np.isfinite(weights) & (weights > 0.0)
    weights = np.where(active, weights, 0.0)
    residual = np.where(active, Y_hat, 0.0) - np.where(active, Y_data, 0.0)
    E_hat = np.where(active, E_hat, 0.0)
    variance = np.where(active, E_data, 0.0)**2 + E_hat**2
    contribution = np.zeros_like(residual, dtype=np.float64)
    uncertain = variance > 0.0
    np.divide(residual**2, variance, out=contribution, where=uncertain)
    contribution[(~uncertain) & (residual != 0.0)] = np.inf
    contribution = np.where(active, weights * contribution, 0.0)
    chi2 = np.sum(contribution, axis=-1)

    slope = np.zeros_like(residual, dtype=np.float64)
    np.divide(residual, variance, out=slope, where=uncertain)
    sigma_term = (weights * slope) ** 2 * E_hat**2
    chi2_error = 2.0 * np.sqrt(np.sum(np.where(active, sigma_term, 0.0), axis=-1))
    return chi2, chi2_error

# Likelihood construction and history loading for icetune
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import copy
import math
import os

import numpy as np

from core.io.serialize import finite_or_none
from core.numerics import array
from core.tune.io import read_json
from core.tune.parameters.space import to_unit


# Compute one canonical optimizer-objective block
def _objective_block(*, name: str, value, ndf, valid: bool) -> dict:
    value = finite_or_none(value)
    ndf = finite_or_none(ndf)
    return {
        "name": str(name),
        "value": value if valid else None,
        "ndf": ndf if valid else None,
        "reduced": finite_or_none(value / ndf) if valid and ndf is not None and ndf > 0.0 else None,
    }


# Compute the shared schema envelope for every icetune likelihood kind
def _likelihood_payload(
    *,
    kind: str,
    covariance_mode: str | None,
    objective: dict,
    valid: bool,
    observables: list[dict] | None = None,
    sigma_objective=None,
    **fields,
) -> dict:
    return {
        "schema_version": 1,
        "valid": bool(valid),
        "kind": str(kind),
        "covariance_mode": covariance_mode,
        "objective": objective,
        "observables": observables or [],
        "mc_uncertainty": {
            "source": "single_trial",
            "sigma_objective": finite_or_none(sigma_objective),
            "replicate_count": 1,
        },
        **fields,
    }


def ll2array(x, dtype=np.float64):
    """
    Convert a nested structure (list of datasets -> list of subsets -> dict of observables)
    into a flattened numpy array.
    """

    return np.asarray([array.to_numpy(value) for dataset in x for subset in dataset for value in subset.values()], dtype=dtype)


def iter_observable_costs(cost, ndf, weight, valid):
    """Yield compact observable-level cost records from nested icetune arrays"""

    if cost is None or ndf is None or weight is None or valid is None:
        return
    for i, dataset in enumerate(cost):
        for j, subset in enumerate(dataset):
            for obs_key, chi2 in subset.items():
                yield {
                    "dataset": int(i),
                    "subset": int(j),
                    "observable": str(obs_key),
                    "chi2": finite_or_none(chi2),
                    "ndf": finite_or_none(ndf[i][j].get(obs_key)),
                    "fitw": finite_or_none(weight[i][j].get(obs_key)),
                    "valid": bool(valid[i][j].get(obs_key, 0) == 1),
                }


# Convert nested histogram costs into canonical observable likelihood records
def build_observable_records(*, cost, ndf, weight, valid) -> list[dict]:
    records = []
    for record in iter_observable_costs(cost, ndf, weight, valid):
        item = dict(record)
        chi2 = item.pop("chi2")
        fitw = item.pop("fitw")
        item["objective_value"] = finite_or_none(chi2 * fitw) if chi2 is not None and fitw is not None else None
        item["weight"] = fitw
        records.append(item)
    return records


def build_likelihood_payload(*, cost, ndf, weight, valid, rho: str = "quadratic") -> dict:
    """Build a canonical independent-bin Gaussian objective payload"""

    per_observable = list(iter_observable_costs(cost, ndf, weight, valid))
    valid_records = [
        r
        for r in per_observable
        if r["valid"]
        and r["ndf"] is not None
        and r["fitw"] is not None
        and r["ndf"] > 0.0
        and r["fitw"] > 0.0
    ]
    objective_total = sum(r["fitw"] * r["chi2"] if r["chi2"] is not None else math.inf for r in valid_records)
    ndf_total = sum(r["ndf"] for r in valid_records)
    is_valid = len(valid_records) > 0 and math.isfinite(objective_total) and math.isfinite(ndf_total)
    observables = build_observable_records(cost=cost, ndf=ndf, weight=weight, valid=valid)

    return _likelihood_payload(
        kind="gaussian_chi2",
        covariance_mode="diagonal",
        objective=_objective_block(name="chi2", value=objective_total, ndf=ndf_total, valid=is_valid),
        valid=is_valid,
        observables=observables,
        residual_model=str(rho),
    )


# Build one canonical covariance-aware Gaussian objective payload
def build_joint_likelihood_payload(*, joint: dict, covariance_mode: str) -> dict:
    value = finite_or_none(joint.get("chi2"))
    sigma_value = finite_or_none(joint.get("sigma_chi2"))
    ndf = finite_or_none(joint.get("ndf"))
    valid = bool(joint.get("valid", False) and value is not None and ndf is not None)
    observables = [{**record, "objective_value": finite_or_none(array.to_numpy(record.get("objective_value")))}
                   for record in joint.get("observables", [])]
    components = {"bins": int(joint.get("bins", 0)), "covariance_rank": int(joint.get("rank", 0))}
    if observables:
        components["observable_diagnostics"] = "joint_precision_partition_additive"
    return _likelihood_payload(
        kind="gaussian_chi2",
        covariance_mode=str(covariance_mode),
        objective=_objective_block(name="chi2", value=value, ndf=ndf, valid=valid),
        valid=valid,
        observables=observables,
        sigma_objective=sigma_value,
        components=components,
    )


def build_metric_likelihood_payload(*, metrics: dict, cost_key: str, error=None) -> dict:
    """Build a scalar cost payload without assigning likelihood semantics"""

    value = finite_or_none(metrics.get(cost_key))
    payload = _likelihood_payload(
        kind="optimizer_cost",
        covariance_mode=None,
        objective=_objective_block(name=cost_key, value=value, ndf=None, valid=True),
        valid=False,
    )
    if error is not None:
        payload["error"] = str(error)
    return payload


# Convert one likelihood payload to the four canonical scalar conventions
def _likelihood_scalars(likelihood: dict) -> tuple:
    if likelihood.get("schema_version") != 1:
        raise ValueError("Likelihood history requires likelihood schema version 1")
    objective = likelihood.get("objective", {}) or {}
    objective_name = objective.get("name")
    if objective_name in {"chi2", "gaussian", "nll"}:
        objective_value = finite_or_none(objective.get("value"))
        value_chi2 = objective_value if objective_name == "chi2" else None
        value_nll = (
            objective_value / 2.0 if objective_name in {"chi2", "gaussian"} and objective_value is not None else objective_value
        )
        value_logL = -value_nll if value_nll is not None else None
        value_two_nll = 2.0 * value_nll if value_nll is not None else None
    else:
        value_logL = value_nll = value_two_nll = value_chi2 = None
    return value_logL, value_nll, value_two_nll, value_chi2, objective_name


# Normalize trial Monte Carlo uncertainty to the canonical 2NLL convention
def _normalized_mc_uncertainty(likelihood: dict, objective_name: str | None):
    uncertainty = copy.deepcopy(likelihood.get("mc_uncertainty"))
    if not isinstance(uncertainty, dict):
        return uncertainty
    sigma_two_nll = finite_or_none(uncertainty.get("sigma_two_nll"))
    sigma_objective = finite_or_none(uncertainty.get("sigma_objective"))
    sigma_logl = finite_or_none(uncertainty.get("sigma_logL"))
    if sigma_two_nll is None and sigma_objective is not None:
        scale = {"chi2": 1.0, "gaussian": 1.0, "nll": 2.0}.get(objective_name)
        sigma_two_nll = sigma_objective * scale if scale is not None else None
    if sigma_two_nll is None and sigma_logl is not None:
        sigma_two_nll = 2.0 * sigma_logl
    uncertainty["sigma_two_nll"] = sigma_two_nll
    return uncertainty


# Extract one sparse observable row while accumulating shared column metadata
def _observable_history_row(likelihood: dict, metadata: dict) -> dict:
    row = {}
    records = likelihood.get("observables", likelihood.get("per_observable", []))
    for record in records or []:
        if not isinstance(record, dict):
            continue
        try:
            key = (int(record["dataset"]), int(record["subset"]), str(record["observable"]))
        except (KeyError, TypeError, ValueError):
            continue
        metadata.setdefault(
            key,
            {
                "dataset": key[0],
                "subset": key[1],
                "observable": key[2],
                "ndf": finite_or_none(record.get("ndf")),
                "fitw": finite_or_none(record.get("weight", record.get("fitw"))),
            },
        )
        value = finite_or_none(record.get("objective_value", record.get("chi2")))
        if bool(record.get("valid", False)) and value is not None:
            row[key] = value
    return row


# Densify sparse observable histories in stable metadata order
def _observable_history_matrix(rows: list[dict], keys: list[tuple]) -> np.ndarray:
    matrix = np.full((len(rows), len(keys)), np.nan, dtype=np.float64)
    index = {key: column for column, key in enumerate(keys)}
    for row_number, row in enumerate(rows):
        for key, value in row.items():
            matrix[row_number, index[key]] = value
    return matrix


# Load scalar and observable-level likelihood arrays from one history snapshot
def load_likelihood_history(path: str | os.PathLike[str]) -> dict:
    path = os.fspath(path)
    if os.path.isdir(path):
        path = os.path.join(path, "history.json")
    history = read_json(path, default={}) or {}
    campaign = read_json(os.path.join(os.path.dirname(path), "icetune_campaign.json"), default={}) or {}
    physics = campaign.get("identity", {}).get("physics", {})
    parameter_space = history.get("parameter_space", [])
    param_names = [item["name"] for item in parameter_space]
    lower = np.array([float(item["lower"]) for item in parameter_space], dtype=np.float64)
    upper = np.array([float(item["upper"]) for item in parameter_space], dtype=np.float64)

    scalar_names = ("logL", "nll", "two_nll", "chi2")
    scalars = {name: [] for name in scalar_names}
    X, valid = [], []
    hashes = []
    mc_uncertainty = []
    covariance_modes = []
    observable_rows = []
    observable_metadata = {}
    for trial in history.get("trials", []):
        likelihood = trial.get("likelihood", {}) or {}
        theta = trial.get("theta")
        if theta is None or len(theta) != len(param_names):
            continue
        covariance_modes.append(likelihood.get("covariance_mode"))
        values = _likelihood_scalars(likelihood)
        value_logL, value_nll, value_two_nll, value_chi2, objective_name = values
        is_valid = bool(likelihood.get("valid", False)) and value_logL is not None
        X.append([float(x) for x in theta])
        for name, value in zip(scalar_names, values[:4], strict=True):
            scalars[name].append(value if value is not None else np.nan)
        valid.append(is_valid)
        hashes.append(trial.get("theta_hash"))
        mc_uncertainty.append(_normalized_mc_uncertainty(likelihood, objective_name))
        observable_rows.append(_observable_history_row(likelihood, observable_metadata))

    X = np.array(X, dtype=np.float64)
    if X.size == 0:
        X = np.empty((0, len(param_names)), dtype=np.float64)
    observable_keys = sorted(observable_metadata)
    observable_chi2 = _observable_history_matrix(observable_rows, observable_keys)
    X_unit = to_unit(X, np.column_stack((lower, upper)))
    distinct_covariance_modes = {str(mode) for mode in covariance_modes if mode is not None}
    if len(distinct_covariance_modes) > 1:
        raise ValueError("Likelihood history mixes incompatible covariance modes")
    covariance_mode = next(iter(distinct_covariance_modes)) if distinct_covariance_modes else None
    return {
        "X": X,
        "X_unit": X_unit,
        **{name: np.array(scalars[name], dtype=np.float64) for name in scalar_names},
        "valid": np.array(valid, dtype=bool),
        "theta_hash": hashes,
        "mc_uncertainty": mc_uncertainty,
        "covariance_mode": covariance_mode,
        "observable_chi2": observable_chi2,
        "observables": [copy.deepcopy(observable_metadata[key]) for key in observable_keys],
        "param_names": param_names,
        "parameter_space": copy.deepcopy(parameter_space),
        "optimization": copy.deepcopy(history.get("optimization", {})),
        "parameter_topology": copy.deepcopy(history.get("parameter_topology", {})),
        "simdriver": copy.deepcopy(history.get("simdriver")),
        "mc_steer": copy.deepcopy(history.get("mc_steer") or physics.get("mc_steer", {})),
        "plot_brand": copy.deepcopy(history.get("plot_brand")),
    }


# Compute observable coefficients shared by scalar costs and their parameter derivatives
def cost_weights(cost: list, ndf: list, weight: list, valid: list, cost_avg: str) -> np.ndarray:
    flat_cost, flat_ndf, flat_weight = (ll2array(value) for value in (cost, ndf, weight))
    flat_valid = ll2array(valid, dtype=np.int64)
    if len({len(value) for value in (flat_cost, flat_ndf, flat_weight, flat_valid)}) != 1:
        raise ValueError("Observable cost input lengths differ")
    active = ((flat_valid == 1) & (np.isfinite(flat_cost) | np.isposinf(flat_cost))
              & np.isfinite(flat_ndf) & np.isfinite(flat_weight) & (flat_ndf > 0.0) & (flat_weight > 0.0))
    output = np.zeros_like(flat_cost)
    ndf, weights = flat_ndf[active], flat_weight[active]
    if not len(ndf):
        return output
    if cost_avg == "dataset-mean":
        dataset = np.repeat(np.arange(len(cost)), [sum(len(subset) for subset in item) for item in cost])[active]
        bins = np.bincount(dataset, weights=weights * ndf)
        output[active] = weights / bins[dataset] / np.count_nonzero(bins > 0.0)
    elif cost_avg in {"local-mean", "local-mean-unweighted", "global-mean", "global-mean-unweighted"}:
        factors = np.ones_like(weights) if cost_avg.endswith("-unweighted") else weights * len(weights) / np.sum(weights)
        output[active] = factors / (len(ndf) * ndf if cost_avg.startswith("local") else np.sum(ndf))
    elif cost_avg == "sum":
        output[active] = weights
    else:
        raise ValueError(f"Unknown observable cost normalization {cost_avg!r}")
    return output


# Combine observable costs with bin, histogram, or dataset normalization
def combine_costs(cost: list, ndf: list, weight: list, valid: list, cost_avg: str) -> float:
    factors = cost_weights(cost, ndf, weight, valid, cost_avg)
    active = factors > 0.0
    if not np.any(active):
        return np.inf
    values = array.stack([value for dataset in cost for subset in dataset for value in subset.values()])
    return array.asarray(factors[active], like=values) @ values[active]

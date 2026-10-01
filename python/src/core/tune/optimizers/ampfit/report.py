# Local autograd covariance and shared icescape parameter plots for ampfit
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import time
from pathlib import Path

import numpy as np
import torch

from core.inference import fit
from core.io.serialize import json_safe, write_json_file
from core.numerics import lbfgsb
from core.stats.transform import covariance_summary


# Save local curvature, both parameter bases and the same figures used by icescape
def save_covariance(*, objective, values, names, bounds, transform, output_root, metadata, plot_brand="GRANIITTI"):
    output_root = Path(output_root)
    point = values.detach().cpu().numpy()
    value, gradient = lbfgsb.value_gradient(objective, values)
    gradient = gradient.detach().cpu().numpy()
    diagnostics = fit.parameter_bound_diagnostics(point, bounds, gradient=gradient)
    free = torch.tensor([not item["at_bound"] for item in diagnostics], dtype=torch.bool)
    # Differentiate one direction at a time without duplicating event tensors for every parameter
    started = time.monotonic()
    hessian, covariance, curvature = lbfgsb.hessian_covariance(objective, values, free_mask=free, vectorize=False)
    curvature["elapsed_seconds"] = time.monotonic() - started
    valid = curvature["free_positive_definite"] and curvature["regularized_eigenvalues"] == 0
    if not valid:
        covariance[:] = torch.nan
    covariance = covariance.cpu().numpy()
    errors, correlation = covariance_summary(covariance)
    optimizer = dict(basis="optimizer", names=names, values=point, covariance=covariance,
                     errors=errors, correlation=correlation, bound_diagnostics=diagnostics)
    fit_label = f"ampfit {metadata['trial_id']}"
    fit.plot_parameter_uncertainties(param_names=names, best_fit=point, errors=errors, bounds=bounds,
                                     output_root=output_root / "optimizer", realized_best=point,
                                     plot_brand=plot_brand, fit_label=fit_label)
    if np.any(np.isfinite(correlation)):
        fit.plot_correlation_matrix(correlation=correlation, param_names=names,
                                    output_root=output_root / "optimizer", plot_brand=plot_brand)
    physical = fit.save_physical_parameters(transform=transform, best_fit=point, covariance=covariance,
                                            output_root=output_root, plot_brand=plot_brand, fit_label=fit_label)
    metadata = {**metadata, "updated_at_unix": time.time(), "objective": {"name": "gaussian" if metadata.get("selected_cost") == "gaussian" else "chi2", "value": float(value)},
                "uncertainty_method": "torch_autograd_hessian", "curvature": curvature,
                "status": "local" if valid and curvature["free_dimension"] else "unavailable",
                "gradient": gradient, "hessian": hessian.cpu().numpy(), "hessian_basis": "optimizer", "hessian_names": names}
    for result in (optimizer, physical):
        directory = output_root / result["basis"]
        result["relative_error_percent"] = 100.0 * fit.parameter_relative_uncertainties(result["values"], result["errors"])
        write_json_file(directory / "parameters.json", json_safe({**result, **metadata}), indent=4)
        rows = [[name, f"{value:.4g}", f"{error:.4g}" if np.isfinite(error) else "N/A",
                 f"{relative:.4g}" if np.isfinite(relative) else "N/A"]
                for name, value, error, relative in zip(result["names"], result["values"], result["errors"],
                                                       result["relative_error_percent"], strict=True)]
        title = f"ampfit {result['basis']} parameters at {metadata['trial_id']} ({metadata['status']})"
        table = fit.format_table(["Parameter", "Value", "Uncertainty", "Relative uncertainty [%]"], rows)
        text = title + "\n" + table + "\n"
        (directory / "parameters.txt").write_text(text, encoding="utf-8")
        print(text, flush=True)
    write_json_file(output_root / "covariance.json", json_safe(metadata), indent=4)
    return metadata


# Evaluate prepared amplitude banks with the recorded statistical objective
def amplitude_objective(driver, param, names):
    settings = {**param, "cost": "gaussian" if param["cost"] == "gaussian" else "chi2", "cost_avg": "sum"}
    amplitude = dict(directory=driver.amplitude_directory(cdir=param["cdir"], run_name=param["run_name"]),
                     prepare=False, controls=param["mc_steer"]["ampfit"], workers=param.get("ampfit_workers", {}))

    # Preserve the Torch graph through the existing histogram and covariance calculations
    def objective(values):
        parameters = dict(zip(names, values, strict=True)) | param["aux_param_space"]
        mc = driver.compute(tunename=param["mc_steer"]["tune_default"], datacards=param["datacards"],
                            mc_steer=param["mc_steer"], cdir=param["cdir"], processes=1,
                            max_t=param["max_t"], rngseed=param["rngseed"],
                            amplitude={**amplitude, "parameters": parameters})
        results = dict(mc=mc, data=driver.data, obs=driver.obs, datasets=driver.datasets)
        return driver.trial_costs(results, settings)["metrics"][settings["cost"]]

    return objective


# Compute one accepted trial covariance independently of its histogram figures
def trial_covariance(*, driver, record, param, output_root):
    names = sorted(record["config"])
    values = torch.tensor([record["config"][name] for name in names], dtype=torch.float64)
    limits = {item["name"]: (item["lower"], item["upper"]) for item in param["parameter_space"]}
    bounds = np.array([limits[name] for name in names], dtype=float)
    transform = driver.parameter_transform(names, values.tolist(), cdir=param["cdir"], metadata=param)
    return save_covariance(objective=amplitude_objective(driver, param, names), values=values, names=names,
                           bounds=bounds, transform=transform, output_root=output_root,
                           metadata={"trial_id": record["trial_id"], "theta_hash": record.get("theta_hash"),
                                     "metrics": record["metrics"], "selected_cost": param["cost"]},
                           plot_brand=param.get("plot_brand") or "GRANIITTI")

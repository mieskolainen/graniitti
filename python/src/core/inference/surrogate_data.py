# Shared full-histogram surface loading for surrogate workflows
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import numpy as np

from core.inference import fit, mcmc
from core.stats import cov, objective
from core.tune.parameters import tools
from core.tune.parameters import topology as parameter_topology


# Load aligned full-trial histograms and reconstruct their global likelihood
def load_histogram_surface(args, *, active_nonzero_bins: bool) -> dict:
    trials = mcmc.collect_simu(
        run_name=args.run_name, cdir=args.cdir, max_trials=args.max_trials, use_cached=args.cached
    )
    X = np.asarray(trials["X"], dtype=np.float64)
    mc = np.asarray(trials["Y"], dtype=np.float64)
    mc_error = np.asarray(trials["E"], dtype=np.float64)
    data = np.asarray(trials["Y_data"], dtype=np.float64)
    data_error = np.asarray(trials["E_data"], dtype=np.float64)
    manifest = trials["histogram_manifest"]
    valid = np.concatenate([np.asarray(item["valid"], dtype=bool) for item in manifest])
    active = valid & ((data != 0.0) | np.any(mc != 0.0, axis=0)) if active_nonzero_bins else valid
    fit_weights = cov.manifest_fit_weights(manifest)
    covariance_payload = covariance = covariance_indices = None
    if cov.uses_full_data_covariance(trials):
        covariance_payload, covariance_indices, covariance = cov.load_manifest_covariance(
            cdir=args.cdir, run_name=args.run_name, manifest=manifest, fit_weights=fit_weights
        )
        Z, Z_error = objective.correlated_chi2_batch(
            data_values=data[covariance_indices],
            mc_prediction=mc[:, covariance_indices],
            mc_stat_uncertainty=mc_error[:, covariance_indices],
            fit_weights=fit_weights[covariance_indices],
            data_total_covariance=covariance,
            covariance_structure=covariance_payload.get("_selected_covariance_structure"),
        )
        covariance_mode = "full"
    else:
        Z, Z_error = objective.global_histogram_chi2(
            Y_data=data, E_data=data_error, Y_hat=mc, E_hat=mc_error, mask=active, fit_weights=fit_weights
        )
        covariance_mode = "diagonal"

    param_names = [str(name) for name in trials["config_keys"]]
    bounds = fit.parameter_bounds(X)
    topology = tools.build_parameter_topology(param_names)
    return {
        "X": X,
        "X_model": parameter_topology.model_features(X, bounds, param_names, topology),
        "Y_data": data,
        "E_data": data_error,
        "Y_hat": mc,
        "E_hat": mc_error,
        "Z": np.asarray(Z, dtype=np.float64),
        "Z_err": np.asarray(Z_error, dtype=np.float64),
        "param_names": param_names,
        "bounds": bounds,
        "parameter_topology": topology,
        "model_feature_names": parameter_topology.model_feature_names(param_names, topology),
        "histogram_manifest": manifest,
        "likelihood_mask": active,
        "fit_weights": fit_weights,
        "covariance_mode": covariance_mode,
        "covariance_payload": covariance_payload,
        "data_total_covariance": covariance,
        "covariance_indices": covariance_indices,
        "simdriver": trials.get("simdriver"),
        "mc_steer": trials.get("mc_steer", {}),
        "plot_brand": trials.get("plot_brand"),
    }

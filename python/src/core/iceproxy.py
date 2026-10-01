# Fit a histogram surrogate from full icetune trial outputs
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import json
import pathlib
import sys

import numpy as np
import torch
from termcolor import cprint
from tqdm import tqdm

from core import __AUTHOR__, __RELEASE__, __version__
from core.inference import fit as fitting
from core.inference import gp as gp_model
from core.inference import impact, mcmc, neural, surrogate_cli, surrogate_data
from core.io.files import ensure_dir
from core.io.serialize import json_safe
from core.stats import hist, objective
from core.stats import objective as icecost
from core.tune import summary as fit_summary
from core.tune.parameters import topology as parameter_topology


# Parse command line arguments using icescape-compatible names
def parse_args():
    return surrogate_cli.parse_args(
        description=f"%(prog)s {__version__} {__RELEASE__} [{__AUTHOR__}]",
        replica=True,
    )


# Compute the persisted backend name for one CLI model choice
def replica_backend_name(model: str) -> str:
    return {
        "rff": "pca_rff_ridge",
        "gp": "gp_pca",
        "neural": "neural",
    }[model]


# Re-export shared validation summaries
print_gp_validation_summary = fitting.print_gp_validation_summary
print_neural_validation_summary = fitting.print_neural_validation_summary


# Re-export shared topology, runtime, split and likelihood machinery
model_features = parameter_topology.model_features
model_features_torch = parameter_topology.model_features_torch
configure_torch_runtime = fitting.configure_torch_runtime
training_validation_indices = fitting.training_validation_indices
global_2NLL = objective.global_histogram_chi2


# Convert one optional array into a typed Torch tensor
def optional_tensor(value, *, dtype: torch.dtype, device=None):
    return None if value is None else torch.as_tensor(value, dtype=dtype, device=device)


# Load full histogram tensors from core.icetune pickle outputs
def load_replica_surface(args) -> dict:
    cprint("Using full histogram replica input path", "yellow")
    template = mcmc.collect_replica_template(
        run_name=args.run_name,
        cdir=args.cdir,
    )
    surface = surrogate_data.load_histogram_surface(args, active_nonzero_bins=True)
    surface.update(
        {
            "plot_results": template["results"],
            "plot_template_source": template["source"],
            "simdriver": surface.get("simdriver") or template.get("simdriver"),
            "plot_brand": surface.get("plot_brand") or template.get("plot_brand"),
        }
    )
    return surface


# Apply positive log transforms and training-only standardization
def prepare_targets(
    Y_hat: np.ndarray,
    E_hat: np.ndarray,
    Y_data: np.ndarray,
    E_data: np.ndarray,
    fit_weights: np.ndarray,
    data_total_covariance: np.ndarray | None,
    covariance_indices: np.ndarray | None,
    training_indices: np.ndarray | None = None,
) -> dict:
    Y_hat = np.asarray(Y_hat, dtype=np.float64)
    E_hat = np.asarray(E_hat, dtype=np.float64)
    floor = 1.0e-12
    transformed = np.hstack(
        (
            np.log(np.maximum(Y_hat, floor)),
            np.log(np.maximum(E_hat, floor)),
        )
    )
    training = (
        np.arange(len(transformed), dtype=np.int64)
        if training_indices is None
        else np.asarray(training_indices, dtype=np.int64)
    )
    scaled, target_mu, target_std = fitting.standardize(transformed, training)

    return {
        "Y_train_target": scaled,
        "E_train_target": np.ones_like(transformed),
        "target_mu": target_mu,
        "target_std": target_std,
        "n_bins": int(Y_hat.shape[1]),
        "target_floor": floor,
        "data_counts": np.asarray(Y_data, dtype=np.float64),
        "data_errs": np.asarray(E_data, dtype=np.float64),
        "fit_weights": np.asarray(fit_weights, dtype=np.float64),
        "data_total_covariance": (
            None
            if data_total_covariance is None
            else np.asarray(data_total_covariance, dtype=np.float64)
        ),
        "covariance_indices": (
            None if covariance_indices is None else np.asarray(covariance_indices, dtype=np.int64)
        ),
    }


# Build the covariance-aware likelihood tensors used at runtime and in checkpoints
def likelihood_tensor_state(
    preprocess: dict, likelihood_mask: np.ndarray, dtype: torch.dtype, device=None
) -> dict:
    return {
        "data_counts": torch.as_tensor(preprocess["data_counts"], dtype=dtype, device=device),
        "data_errs": torch.as_tensor(preprocess["data_errs"], dtype=dtype, device=device),
        "fit_weights": torch.as_tensor(preprocess["fit_weights"], dtype=dtype, device=device),
        "likelihood_mask": torch.as_tensor(likelihood_mask, dtype=torch.bool, device=device),
        "data_total_covariance": optional_tensor(
            preprocess.get("data_total_covariance"), dtype=dtype, device=device
        ),
        "covariance_indices": optional_tensor(
            preprocess.get("covariance_indices"), dtype=torch.long, device=device
        ),
    }


# Move replica preprocessing arrays onto the selected Torch device
def torch_preprocess(
    preprocess: dict, likelihood_mask: np.ndarray, device: str | torch.device, dtype: torch.dtype
) -> dict:
    return {
        "target_mu": torch.as_tensor(preprocess["target_mu"], dtype=dtype, device=device),
        "target_std": torch.as_tensor(preprocess["target_std"], dtype=dtype, device=device),
        "n_bins": int(preprocess["n_bins"]),
        **likelihood_tensor_state(preprocess, likelihood_mask, dtype, device),
    }


# Decode standardized log targets into positive counts and MC errors
def decode_targets_torch(
    prediction_scaled: torch.Tensor, preprocess_torch: dict
) -> tuple[torch.Tensor, torch.Tensor]:
    transformed = prediction_scaled * preprocess_torch["target_std"] + preprocess_torch["target_mu"]
    n_bins = int(preprocess_torch["n_bins"])
    counts = torch.exp(transformed[:, :n_bins])
    errors = torch.exp(transformed[:, n_bins:])
    return counts, errors


# Propagate standardized target variance through the positive log inverse
def decode_target_variance_torch(
    prediction_scaled: torch.Tensor, variance_scaled: torch.Tensor, preprocess_torch: dict
) -> tuple[torch.Tensor, torch.Tensor]:
    counts, errors = decode_targets_torch(prediction_scaled, preprocess_torch)
    transformed_variance = variance_scaled * preprocess_torch["target_std"].square()
    n_bins = int(preprocess_torch["n_bins"])
    count_variance = counts.square() * transformed_variance[:, :n_bins]
    error_variance = errors.square() * transformed_variance[:, n_bins:]
    return count_variance, error_variance


# Bind predicted histograms to the shared covariance-aware likelihood tensors
def histogram_likelihood_tensors(
    counts: torch.Tensor, errors: torch.Tensor, preprocess_torch: dict
) -> dict:
    return {
        "counts": counts,
        "errors": errors,
        "data_counts": preprocess_torch["data_counts"],
        "data_errors": preprocess_torch["data_errs"],
        "mask": preprocess_torch["likelihood_mask"],
        "fit_weights": preprocess_torch["fit_weights"],
        "data_total_covariance": preprocess_torch["data_total_covariance"],
        "covariance_indices": preprocess_torch["covariance_indices"],
    }


# Evaluate the differentiable tune chi2 from predicted histograms
def replica_Z_torch(
    counts: torch.Tensor, errors: torch.Tensor, preprocess_torch: dict
) -> torch.Tensor:
    return icecost.histogram_chi2_torch(
        **histogram_likelihood_tensors(counts, errors, preprocess_torch)
    )


# Evaluate one chi2 column per original observable
def observable_chi2_torch(
    counts: torch.Tensor, errors: torch.Tensor, preprocess_torch: dict, manifest: list[dict]
) -> torch.Tensor:
    contribution = icecost.histogram_chi2_contributions_torch(
        **histogram_likelihood_tensors(counts, errors, preprocess_torch)
    )
    return torch.stack(
        [torch.sum(contribution[:, item["start"] : item["stop"]], dim=1) for item in manifest],
        dim=1,
    )


# Estimate held-out standardized residual uncertainty per target component
def validation_residual_scale(
    target_scaled: np.ndarray, prediction_scaled: np.ndarray
) -> np.ndarray:
    residual = np.asarray(prediction_scaled) - np.asarray(target_scaled)
    scale = np.sqrt(np.mean(residual**2, axis=0))
    return np.where(np.isfinite(scale) & (scale > 1.0e-6), scale, 1.0e-6)


# Collect the held-out residual statistics shared by every replica backend
def validation_statistics(target_scaled: np.ndarray, prediction_scaled: np.ndarray) -> dict:
    residual = np.asarray(prediction_scaled) - np.asarray(target_scaled)
    return {
        "validation_mae_scaled": float(np.mean(np.abs(residual))),
        "validation_rmse_scaled": float(np.sqrt(np.mean(residual**2))),
        "residual_scale": validation_residual_scale(target_scaled, prediction_scaled),
    }


# Fit one joint PCA basis over all histogram bins and observables
def fit_joint_pca(
    Y_scaled: np.ndarray, variance: float, max_modes: int, device: str = "cpu"
) -> tuple[np.ndarray, dict]:
    Y_scaled = np.asarray(Y_scaled, dtype=np.float64)
    pca_mean = np.mean(Y_scaled, axis=0)
    centered = Y_scaled - pca_mean
    if str(device).startswith("cuda"):
        centered_tensor = torch.as_tensor(
            centered,
            dtype=torch.float32,
            device=device,
        )
        _, singular_tensor, components_tensor = torch.linalg.svd(
            centered_tensor,
            full_matrices=False,
        )
        singular_values = singular_tensor.cpu().numpy().astype(np.float64)
        components_t = components_tensor.cpu().numpy().astype(np.float64)
    else:
        _, singular_values, components_t = np.linalg.svd(
            centered,
            full_matrices=False,
        )

    power = singular_values**2
    total_power = np.sum(power)
    if total_power > 0.0:
        explained_ratio = power / np.sum(power)
        cumulative = np.cumsum(explained_ratio)
        n_modes = int(np.searchsorted(cumulative, variance, side="left") + 1)
    else:
        explained_ratio = np.zeros_like(power)
        if len(explained_ratio) > 0:
            explained_ratio[0] = 1.0
        n_modes = 1
    n_modes = max(1, min(n_modes, int(max_modes), components_t.shape[0]))

    components = components_t[:n_modes]
    scores = centered @ components.T
    pca = {
        "mean": pca_mean,
        "components": components,
        "singular_values": singular_values[:n_modes],
        "explained_variance_ratio": explained_ratio[:n_modes],
        "retained_variance": float(np.sum(explained_ratio[:n_modes])),
        "n_modes": int(n_modes),
    }
    return scores, pca


# Reconstruct standardized log-histogram vectors from PCA scores
def reconstruct_joint_pca(scores: np.ndarray, pca: dict) -> np.ndarray:
    scores = np.asarray(scores, dtype=np.float64)
    return scores @ np.asarray(pca["components"], dtype=np.float64) + np.asarray(
        pca["mean"], dtype=np.float64
    )


# Project standardized log-histogram vectors onto an existing joint PCA basis
def project_joint_pca(Y_scaled: np.ndarray, pca: dict) -> np.ndarray:
    centered = np.asarray(Y_scaled, dtype=np.float64) - np.asarray(pca["mean"], dtype=np.float64)
    return centered @ np.asarray(pca["components"], dtype=np.float64).T


# Build RFF maps through the shared inference implementation
build_rff_map = impact.build_rff_map
# Evaluate RFF features through the shared inference implementation
rff_feature_matrix = impact.rff_feature_matrix
# Fit ridge coefficients through the shared inference implementation
fit_ridge_coefficients = impact.fit_ridge_coefficients


# Build the affine and random Fourier feature tensor used by every RFF backend
def rff_feature_tensor(
    X_model: torch.Tensor,
    omega: torch.Tensor,
    phase: torch.Tensor,
) -> torch.Tensor:
    nonlinear = np.sqrt(2.0 / omega.shape[1]) * torch.cos(X_model @ omega + phase)
    intercept = torch.ones(
        (len(X_model), 1),
        dtype=X_model.dtype,
        device=X_model.device,
    )
    return torch.cat((intercept, X_model, nonlinear), dim=1)


# Fit one RFF ridge coefficient matrix on the selected arithmetic device
def fit_rff_coefficients(
    X_model: np.ndarray, scores: np.ndarray, rff: dict, alpha: float, device: str
) -> np.ndarray:
    if not str(device).startswith("cuda"):
        return fit_ridge_coefficients(
            rff_feature_matrix(X_model, rff),
            scores,
            alpha,
        )
    X_tensor = torch.as_tensor(X_model, dtype=torch.float64, device=device)
    omega = torch.as_tensor(rff["omega"], dtype=torch.float64, device=device)
    phase = torch.as_tensor(rff["phase"], dtype=torch.float64, device=device)
    Phi = rff_feature_tensor(X_tensor, omega, phase)
    target = torch.as_tensor(scores, dtype=torch.float64, device=device)
    coefficients = impact.fit_ridge_tensor(Phi, target, alpha)
    return coefficients.detach().cpu().numpy()


# Fit one coefficient matrix for every member of an existing RFF ensemble
def fit_rff_ensemble(
    feature_maps: list[dict],
    X_model: np.ndarray,
    scores: np.ndarray,
    alpha: float,
    device: str,
) -> list[dict]:
    return [
        {
            "rff": rff,
            "coef": fit_rff_coefficients(X_model, scores, rff, alpha, device).astype(np.float64),
        }
        for rff in feature_maps
    ]


# Train a validated PCA-latent surrogate and refit it with all rows
def train_pca_rff_ridge(
    X_model: np.ndarray,
    Y_scaled: np.ndarray,
    args,
    training_indices: np.ndarray,
    validation_indices: np.ndarray,
) -> tuple[list[dict], dict, dict]:
    scores_train, pca = fit_joint_pca(
        Y_scaled=Y_scaled[training_indices],
        variance=float(args.pca_variance),
        max_modes=int(args.pca_max_modes),
        device=args.device,
    )
    cprint(
        f"Joint PCA: modes = {pca['n_modes']}, retained variance = {pca['retained_variance']:0.5f}",
        "yellow",
    )
    cprint(
        f"PCA RFF ridge: features = {args.rff_features}, ensemble = {args.rff_ensemble}, "
        f"scale = {args.rff_scale}, alpha = {args.ridge_alpha}",
        "yellow",
    )

    rng = np.random.default_rng(int(args.rngseed) + 104729)
    feature_maps = []
    for _ in range(int(args.rff_ensemble)):
        feature_maps.append(
            build_rff_map(
                input_dim=X_model.shape[1],
                n_features=int(args.rff_features),
                scale=float(args.rff_scale),
                rng=rng,
            )
        )

    validation_models = fit_rff_ensemble(
        feature_maps,
        X_model[training_indices],
        scores_train,
        float(args.ridge_alpha),
        args.device,
    )
    validation_scores = predict_scores_pca_rff_ridge(
        validation_models,
        X_model[validation_indices],
    )
    validation_prediction = reconstruct_joint_pca(validation_scores, pca)
    scores_full = project_joint_pca(Y_scaled, pca)
    models = fit_rff_ensemble(
        feature_maps,
        X_model,
        scores_full,
        float(args.ridge_alpha),
        args.device,
    )
    stats = {
        "validation_count": int(len(validation_indices)),
        **validation_statistics(
            Y_scaled[validation_indices],
            validation_prediction,
        ),
        "_validation_prediction_scaled": validation_prediction,
    }
    return models, pca, stats


# Predict PCA latent scores from a PCA RFF ridge ensemble
def predict_scores_pca_rff_ridge(models: list[dict], X_model: np.ndarray) -> np.ndarray:
    scores = None
    for model in models:
        Phi = rff_feature_matrix(X_model, model["rff"])
        pred = Phi @ np.asarray(model["coef"], dtype=np.float64)
        scores = pred if scores is None else scores + pred
    return scores / float(len(models))


# Move one PCA decoder and its calibrated residual scale onto Torch
def pca_torch_state(pca: dict, residual_scale: np.ndarray, device: str, dtype: torch.dtype) -> dict:
    return {
        "mean": torch.as_tensor(pca["mean"], dtype=dtype, device=device),
        "components": torch.as_tensor(pca["components"], dtype=dtype, device=device),
        "residual_scale": torch.as_tensor(residual_scale, dtype=dtype, device=device),
    }


# Move one fitted PCA-RFF ensemble onto Torch for inference and autograd
def rff_torch_state(
    models: list[dict], pca: dict, residual_scale: np.ndarray, device: str, dtype: torch.dtype
) -> dict:
    return {
        "models": [
            {
                "omega": torch.as_tensor(model["rff"]["omega"], dtype=dtype, device=device),
                "phase": torch.as_tensor(model["rff"]["phase"], dtype=dtype, device=device),
                "coef": torch.as_tensor(model["coef"], dtype=dtype, device=device),
            }
            for model in models
        ],
        "pca": pca_torch_state(pca, residual_scale, device, dtype),
    }


# Predict standardized targets and calibrated variance with Torch RFF
def predict_scaled_rff_torch(
    state: dict, X_model: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor]:
    predictions = []
    for model in state["models"]:
        features = rff_feature_tensor(X_model, model["omega"], model["phase"])
        scores = features @ model["coef"]
        predictions.append(scores @ state["pca"]["components"] + state["pca"]["mean"])
    stacked = torch.stack(predictions, dim=0)
    mean = torch.mean(stacked, dim=0)
    ensemble_variance = torch.var(stacked, dim=0, unbiased=False)
    variance = ensemble_variance + state["pca"]["residual_scale"].square()
    return mean, variance


# Train a scalable independent-output GP and refit its posterior with all rows
def train_gp_pca(
    X_model: np.ndarray,
    Y_scaled: np.ndarray,
    args,
    training_indices: np.ndarray | None = None,
    validation_indices: np.ndarray | None = None,
) -> tuple[gp_model.GaussianProcess, dict, dict]:
    if training_indices is None or validation_indices is None:
        training_indices, validation_indices = training_validation_indices(
            n_rows=len(X_model),
            validation_fraction=args.validation_fraction,
            rngseed=args.rngseed,
            required_training=[int(np.argmin(np.sum(Y_scaled**2, axis=1)))],
        )
    X_train = X_model[training_indices]
    X_val = X_model[validation_indices]
    Y_train_scaled = Y_scaled[training_indices]
    Y_val_scaled = Y_scaled[validation_indices]
    scores_train, pca = fit_joint_pca(
        Y_scaled=Y_train_scaled,
        variance=float(args.pca_variance),
        max_modes=int(args.pca_max_modes),
        device=args.device,
    )
    scores_val = project_joint_pca(Y_val_scaled, pca)
    cprint(
        f"Joint PCA: modes = {pca['n_modes']}, retained variance = {pca['retained_variance']:0.5f}",
        "yellow",
    )
    fitting.print_table(
        "GP-PCA training split:",
        ["Array", "Shape"],
        [
            ["X_train", X_train.shape],
            ["X_val", X_val.shape],
            ["scores_train", scores_train.shape],
            ["scores_val", scores_val.shape],
        ],
    )

    gp = gp_model.GaussianProcess(
        X_train=torch.as_tensor(X_train),
        Y_train=torch.as_tensor(scores_train),
        scale=float(args.gp_scale),
        noise=float(args.gp_noise),
        device=args.device,
        dtype=getattr(args, "gp_dtype", "auto"),
        kernel=getattr(args, "gp_kernel", "Matern"),
        independent_outputs=True,
    )
    gp_stats = {
        "train_loss": [],
        "validation_loss": [],
        "initial_validation_loss": None,
        "best_validation_loss": None,
        "best_step": None,
        "stopped_step": None,
    }
    if args.gp_optimize and args.gp_opt_steps > 0:
        gp_stats = gp.optimize_hyperparameters(
            num_steps=int(args.gp_opt_steps),
            lr=float(args.gp_lr),
            X_val=torch.as_tensor(X_val),
            Y_val=torch.as_tensor(scores_val),
            patience=args.gp_patience,
            min_delta=args.gp_min_delta,
            validation_interval=getattr(args, "gp_validation_interval", 5),
            log_interval=getattr(args, "gp_log_interval", 20),
        )
    with torch.inference_mode():
        validation_scores = (
            gp.predict_mean(torch.as_tensor(X_val, dtype=gp.dtype, device=gp.device)).cpu().numpy()
        )
    validation_prediction = reconstruct_joint_pca(validation_scores, pca)
    gp_stats.update(validation_statistics(Y_val_scaled, validation_prediction))
    gp_stats["_validation_prediction_scaled"] = validation_prediction
    scores_full = project_joint_pca(Y_scaled, pca)
    gp.set_training_data(
        X_train=torch.as_tensor(X_model),
        Y_train=torch.as_tensor(scores_full),
    )
    gp.freeze_for_inference()
    print_gp_validation_summary(gp_stats, title="GP-PCA validation summary:")
    return gp, pca, gp_stats


# Predict standardized targets and calibrated variance with the GP-PCA model
def predict_scaled_gp_torch(
    gp, pca_torch: dict, X_model: torch.Tensor, include_variance: bool = True
) -> tuple[torch.Tensor, torch.Tensor]:
    if include_variance:
        score_mean, score_variance = gp.predict_marginals(X_model)
    else:
        score_mean = gp.predict_mean(X_model)
        score_variance = torch.zeros_like(score_mean)
    mean = score_mean @ pca_torch["components"] + pca_torch["mean"]
    variance = score_variance @ pca_torch["components"].square()
    if include_variance:
        variance = variance + pca_torch["residual_scale"].square()
    return mean, variance


# Construct the histogram-replica neural model parameters
def neural_model_parameters(in_dim: int, out_dim: int, args) -> dict:
    hidden_dim = int(args.nn_hidden_dim)
    if hidden_dim <= 0:
        hidden_dim = min(max(4 * in_dim, 64), 512)
    return neural.lzmlp_parameters(
        in_dim=in_dim,
        out_dim=out_dim,
        hidden_dim=hidden_dim,
        hidden_layers=3,
    )


# Predict standardized neural targets with held-out residual variance
def predict_scaled_neural_torch(
    model, residual_scale: torch.Tensor, X_model: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor]:
    prediction = model(X_model)
    variance = residual_scale.square().expand_as(prediction)
    return prediction, variance


# Find a good starting point with batched on-device random scanning
def random_search(
    surrogate_Z_torch, bounds: np.ndarray, args, dtype: torch.dtype
) -> tuple[np.ndarray, float]:
    total = int(args.random_samples)
    batch_size = int(args.random_batch_size)
    D = bounds.shape[0]
    best_x = None
    best_z = float("inf")
    device = torch.device(args.device)
    limits = torch.as_tensor(bounds, dtype=dtype, device=device)
    generator = torch.Generator(device=device)
    generator.manual_seed(int(args.rngseed))

    cprint(f"Replica random search: {total} surrogate samples, batch size {batch_size}", "yellow")
    with tqdm(total=total, desc="Replica random search") as pbar:
        for start in range(0, total, batch_size):
            n = min(batch_size, total - start)
            unit = torch.rand(
                (n, D),
                dtype=dtype,
                device=device,
                generator=generator,
            )
            X_trial = limits[:, 0] + unit * (limits[:, 1] - limits[:, 0])
            with torch.inference_mode():
                Z = surrogate_Z_torch(X_trial)
            finite = torch.isfinite(Z)
            if bool(torch.any(finite)):
                candidate = torch.where(
                    finite,
                    Z,
                    torch.full_like(Z, torch.inf),
                )
                local_index = int(torch.argmin(candidate))
                local_z = float(candidate[local_index].detach().cpu())
                if local_z < best_z:
                    best_z = local_z
                    best_x = X_trial[local_index].detach().cpu().numpy().copy()
                    pbar.set_postfix(best_Z=f"{best_z:.3g}")
            pbar.update(n)

    if best_x is None:
        raise RuntimeError("Replica random search did not find a finite surrogate point")
    return best_x, best_z


# Propagate model and fitted-parameter uncertainty to predictive histograms
def predictive_histogram_summary(
    predict_histograms_torch,
    point: np.ndarray,
    covariance: np.ndarray,
    device: str,
    dtype: torch.dtype,
) -> dict:
    point_tensor = (
        torch.as_tensor(
            point,
            dtype=dtype,
            device=device,
        )
        .clone()
        .requires_grad_(True)
    )

    # Compute only the central count prediction for Jacobian propagation
    def count_prediction(values: torch.Tensor) -> torch.Tensor:
        counts, _, _, _ = predict_histograms_torch(
            values[None, :],
            include_model_variance=False,
        )
        return counts[0]

    with torch.enable_grad():
        jacobian = torch.autograd.functional.jacobian(
            count_prediction,
            point_tensor,
            vectorize=True,
        )
    with torch.inference_mode():
        counts, mc_errors, model_count_variance, model_error_variance = predict_histograms_torch(
            point_tensor.detach()[None, :],
            include_model_variance=True,
        )
    covariance_tensor = torch.as_tensor(
        covariance,
        dtype=dtype,
        device=device,
    )
    covariance_tensor = torch.where(
        torch.isfinite(covariance_tensor),
        covariance_tensor,
        torch.zeros_like(covariance_tensor),
    )
    parameter_variance = torch.einsum(
        "bi,ij,bj->b",
        jacobian,
        covariance_tensor,
        jacobian,
    )
    parameter_variance = torch.clamp(parameter_variance, min=0.0)
    model_count_variance = torch.clamp(model_count_variance[0], min=0.0)
    total_variance = mc_errors[0].square() + model_count_variance + parameter_variance
    return {
        "counts": counts[0].detach().cpu().numpy(),
        "mc_errors": mc_errors[0].detach().cpu().numpy(),
        "model_count_sigma": torch.sqrt(model_count_variance).cpu().numpy(),
        "model_error_sigma": torch.sqrt(torch.clamp(model_error_variance[0], min=0.0))
        .cpu()
        .numpy(),
        "parameter_sigma": torch.sqrt(parameter_variance).detach().cpu().numpy(),
        "total_sigma": torch.sqrt(torch.clamp(total_variance, min=0.0)).detach().cpu().numpy(),
        "count_jacobian": jacobian.detach().cpu().numpy(),
    }


# Replace one full trial template with surrogate-predictive histogram arrays
def fill_predictive_results(template_results: dict, manifest: list[dict], prediction: dict) -> dict:
    results = copy.deepcopy(template_results)
    for item in manifest:
        start = int(item["start"])
        stop = int(item["stop"])
        record = results["mc"][item["dataset"]][item["subset"]][item["observable"]]
        hdata = record["hdata"]
        counts_scaled = np.asarray(prediction["counts"][start:stop], dtype=np.float64)
        errors_scaled = np.asarray(prediction["total_sigma"][start:stop], dtype=np.float64)
        # Predictions already include the histogram normalization and its uncertainty
        record["hdata"] = hist.hobj(
            counts=counts_scaled,
            errs=errors_scaled,
            bins=hdata.bins,
            cbins=hdata.cbins,
            valid=np.asarray(item["valid"], dtype=np.bool_),
        )
        record["label"] = "Replica optimum"
    return results


# Render predictive histograms with the same driver used by the tune
def render_predictive_histograms(
    *,
    surface: dict,
    prediction: dict,
    point: np.ndarray,
    backend_name: str,
    output_root: pathlib.Path,
    args,
) -> dict:
    output_dir = pathlib.Path(output_root) / "predictive_histograms"
    ensure_dir(output_dir)
    arrays_path = output_dir / "predictive_histograms.npz"
    np.savez_compressed(
        arrays_path,
        optimal_parameters=np.asarray(point, dtype=np.float64),
        **{key: np.asarray(value) for key, value in prediction.items()},
    )
    outputs = {
        "arrays": str(arrays_path),
        "figure_root": str(output_dir),
        "template_source": str(surface["plot_template_source"]),
    }
    template = surface["plot_results"]
    if template.get("datasets") is None:
        outputs["status"] = "skipped: plotting template has no datasets"
        return outputs

    results = fill_predictive_results(
        template_results=template,
        manifest=surface["histogram_manifest"],
        prediction=prediction,
    )
    _, _, _, valid_arr = icecost.compute_costs(
        results=results,
        cost_func="chi2",
        cost_rho="quadratic",
    )
    summary_file = output_dir / "summary.json"
    from core.plot.loss import visualize_losses

    visualize_losses(
        results=results,
        run_name=args.run_name,
        cdir=args.cdir,
        tunename=f"iceproxy {backend_name}",
        summary={
            "run_name": args.run_name,
            "surrogate_model": backend_name,
            "optimal_parameters": {
                name: float(point[index]) for index, name in enumerate(surface["param_names"])
            },
            "uncertainty": (
                "quadrature of predicted MC statistics, held-out surrogate "
                "residual, and fitted-parameter covariance"
            ),
        },
        valid_arr=valid_arr,
        output_dir=str(output_dir),
        summary_file=str(summary_file),
    )
    outputs["summary"] = str(summary_file)
    outputs["status"] = "rendered"
    return outputs


# Measure per-observable fidelity from a held-out standardized prediction
def observable_validation_metrics_from_scaled(
    *,
    prediction_scaled: np.ndarray,
    Y_validation: np.ndarray,
    E_validation: np.ndarray,
    preprocess_torch: dict,
    manifest: list[dict],
) -> tuple[np.ndarray, np.ndarray]:
    prediction = torch.as_tensor(
        prediction_scaled,
        dtype=preprocess_torch["data_counts"].dtype,
        device=preprocess_torch["data_counts"].device,
    )
    with torch.inference_mode():
        counts, errors = decode_targets_torch(prediction, preprocess_torch)
        predicted = observable_chi2_torch(
            counts,
            errors,
            preprocess_torch,
            manifest,
        )
        realized = observable_chi2_torch(
            torch.as_tensor(
                Y_validation,
                dtype=prediction.dtype,
                device=prediction.device,
            ),
            torch.as_tensor(
                E_validation,
                dtype=prediction.dtype,
                device=prediction.device,
            ),
            preprocess_torch,
            manifest,
        )
    predicted_np = predicted.cpu().numpy()
    realized_np = realized.cpu().numpy()
    residual = predicted_np - realized_np
    denominator = np.sum(
        (realized_np - np.mean(realized_np, axis=0, keepdims=True)) ** 2,
        axis=0,
    )
    r2 = 1.0 - np.divide(
        np.sum(residual**2, axis=0),
        denominator,
        out=np.full(realized_np.shape[1], np.nan, dtype=np.float64),
        where=denominator > 0.0,
    )
    return r2, np.mean(np.abs(residual), axis=0)


# Build one physical-coordinate histogram predictor from a backend closure
def make_histogram_predictor(*, predict_scaled_torch, surface: dict, preprocess_torch: dict):
    # Predict positive bin values and their calibrated model variances
    def predict_histograms_torch(X_physical: torch.Tensor, include_model_variance: bool = True):
        X_model = model_features_torch(
            X_physical,
            surface["bounds"],
            surface["param_names"],
            surface["parameter_topology"],
        )
        prediction_scaled, variance_scaled = predict_scaled_torch(
            X_model,
            include_model_variance,
        )
        counts, errors = decode_targets_torch(
            prediction_scaled,
            preprocess_torch,
        )
        if include_model_variance:
            count_variance, error_variance = decode_target_variance_torch(
                prediction_scaled,
                variance_scaled,
                preprocess_torch,
            )
        else:
            count_variance = torch.zeros_like(counts)
            error_variance = torch.zeros_like(errors)
        return counts, errors, count_variance, error_variance

    return predict_histograms_torch


# Serialize the likelihood tensors shared by every replica backend
def replica_likelihood_checkpoint(surface: dict, preprocess: dict) -> dict:
    """Return covariance-aware likelihood state for one saved replica model"""
    return {
        "covariance_mode": str(surface["covariance_mode"]),
        **likelihood_tensor_state(preprocess, surface["likelihood_mask"], torch.float32),
    }


# Serialize metadata and target preprocessing shared by every replica backend
def replica_checkpoint(surface: dict, preprocess: dict, args, model: str) -> dict:
    return {
        "replica_schema_version": 1,
        "surrogate_model": model,
        "training_args": vars(args),
        "param_names": surface["param_names"],
        "parameter_topology": surface["parameter_topology"],
        "model_feature_names": surface["model_feature_names"],
        "bounds": torch.tensor(surface["bounds"], dtype=torch.float32),
        "target_mu": torch.tensor(preprocess["target_mu"], dtype=torch.float32),
        "target_std": torch.tensor(preprocess["target_std"], dtype=torch.float32),
        "n_bins": int(preprocess["n_bins"]),
        "target_floor": float(preprocess["target_floor"]),
        **replica_likelihood_checkpoint(surface, preprocess),
        "histogram_manifest": surface["histogram_manifest"],
    }


# Serialize the PCA state shared by the RFF and GP replica backends
def pca_checkpoint(pca: dict) -> dict:
    return {
        "pca_mean": torch.tensor(pca["mean"], dtype=torch.float32),
        "pca_components": torch.tensor(pca["components"], dtype=torch.float32),
        "pca_singular_values": torch.tensor(pca["singular_values"], dtype=torch.float32),
        "pca_explained_variance_ratio": torch.tensor(
            pca["explained_variance_ratio"], dtype=torch.float32
        ),
        "pca_retained_variance": float(pca["retained_variance"]),
    }


# Persist one replica checkpoint below its backend output directory
def save_replica_checkpoint(checkpoint: dict, filename: str, args) -> pathlib.Path:
    save_dir = pathlib.Path(args.output_root) / replica_backend_name(args.model)
    ensure_dir(save_dir)
    save_path = save_dir / filename
    torch.save(checkpoint, save_path)
    return save_path


# Persist a neural replica surrogate and all preprocessing needed for inference
def save_neural_checkpoint(
    model, model_param: dict, surface: dict, preprocess: dict, stats: dict, args
) -> pathlib.Path:
    checkpoint = {
        **replica_checkpoint(surface, preprocess, args, "neural"),
        "model_state_dict": model.state_dict(),
        "model_param": model_param,
        "train_stats": json_safe(stats),
    }
    return save_replica_checkpoint(checkpoint, "replica_neural_model.pt", args)


# Persist a PCA RFF ridge replica surrogate and all preprocessing needed for inference
def save_pca_rff_ridge_checkpoint(
    models: list[dict], pca: dict, stats: dict, surface: dict, preprocess: dict, args
) -> pathlib.Path:
    checkpoint = {
        **replica_checkpoint(surface, preprocess, args, "pca_rff_ridge"),
        "models": [
            {
                "omega": torch.tensor(model["rff"]["omega"], dtype=torch.float32),
                "phase": torch.tensor(model["rff"]["phase"], dtype=torch.float32),
                "coef": torch.tensor(model["coef"], dtype=torch.float32),
                "scale": float(model["rff"]["scale"]),
                "n_features": int(model["rff"]["n_features"]),
            }
            for model in models
        ],
        "ridge_alpha": float(args.ridge_alpha),
        **pca_checkpoint(pca),
        "train_stats": json_safe(stats),
    }
    return save_replica_checkpoint(checkpoint, "replica_pca_rff_ridge_model.pt", args)


# Persist a GP-PCA replica surrogate and all preprocessing needed for inference
def save_gp_pca_checkpoint(
    gp, pca: dict, gp_stats: dict, surface: dict, preprocess: dict, args
) -> pathlib.Path:
    checkpoint = {
        **replica_checkpoint(surface, preprocess, args, "gp_pca"),
        "gp_state_dict": gp.state_dict(),
        "gp_train_X": gp.X_train.detach().cpu(),
        "gp_train_Y": gp.Y_train.detach().cpu(),
        "gp_scale": float(args.gp_scale),
        "gp_noise": float(args.gp_noise),
        "gp_kernel": str(args.gp_kernel),
        "gp_dtype": str(gp.dtype),
        "gp_independent_outputs": True,
        "gp_train_stats": json_safe(gp_stats),
        **pca_checkpoint(pca),
    }
    return save_replica_checkpoint(checkpoint, "replica_gp_pca_model.pt", args)


# Train and package the PCA-RFF replica backend
def build_rff_replica_backend(
    surface: dict,
    preprocess: dict,
    training_indices: np.ndarray,
    validation_indices: np.ndarray,
    args,
) -> dict:
    models, pca, stats = train_pca_rff_ridge(
        X_model=surface["X_model"],
        Y_scaled=preprocess["Y_train_target"],
        args=args,
        training_indices=training_indices,
        validation_indices=validation_indices,
    )
    validation_prediction = stats.pop("_validation_prediction_scaled")
    dtype = torch.float64
    state = rff_torch_state(
        models=models,
        pca=pca,
        residual_scale=stats["residual_scale"],
        device=args.device,
        dtype=dtype,
    )

    # Evaluate the RFF ensemble mean and optional calibrated variance
    def predict_scaled(X_model: torch.Tensor, include_variance: bool):
        mean, variance = predict_scaled_rff_torch(state, X_model)
        return mean, variance if include_variance else torch.zeros_like(mean)

    checkpoint = save_pca_rff_ridge_checkpoint(
        models=models,
        pca=pca,
        stats=stats,
        surface=surface,
        preprocess=preprocess,
        args=args,
    )
    return {
        "predict_scaled_torch": predict_scaled,
        "dtype": dtype,
        "checkpoint": checkpoint,
        "validation_prediction_scaled": validation_prediction,
        "summary": {
            "pca": {
                "n_modes": int(pca["n_modes"]),
                "retained_variance": float(pca["retained_variance"]),
            },
            "rff_ridge": {
                "features": int(args.rff_features),
                "ensemble": int(args.rff_ensemble),
                "scale": float(args.rff_scale),
                "ridge_alpha": float(args.ridge_alpha),
                **stats,
            },
        },
    }


# Train and package the independent-output GP-PCA replica backend
def build_gp_replica_backend(
    surface: dict,
    preprocess: dict,
    training_indices: np.ndarray,
    validation_indices: np.ndarray,
    args,
) -> dict:
    gp, pca, stats = train_gp_pca(
        X_model=surface["X_model"],
        Y_scaled=preprocess["Y_train_target"],
        args=args,
        training_indices=training_indices,
        validation_indices=validation_indices,
    )
    validation_prediction = stats.pop("_validation_prediction_scaled")
    pca_torch = pca_torch_state(pca, stats["residual_scale"], args.device, gp.dtype)

    # Evaluate the GP-PCA posterior mean and optional calibrated variance
    def predict_scaled(X_model: torch.Tensor, include_variance: bool):
        return predict_scaled_gp_torch(
            gp,
            pca_torch,
            X_model,
            include_variance=include_variance,
        )

    checkpoint = save_gp_pca_checkpoint(
        gp=gp,
        pca=pca,
        gp_stats=stats,
        surface=surface,
        preprocess=preprocess,
        args=args,
    )
    return {
        "predict_scaled_torch": predict_scaled,
        "dtype": gp.dtype,
        "checkpoint": checkpoint,
        "validation_prediction_scaled": validation_prediction,
        "summary": {
            "pca": {
                "n_modes": int(pca["n_modes"]),
                "retained_variance": float(pca["retained_variance"]),
            },
            "gp_validation": stats,
        },
    }


# Train, validate and full-refit the neural replica backend
def build_neural_replica_backend(
    surface: dict,
    preprocess: dict,
    training_indices: np.ndarray,
    validation_indices: np.ndarray,
    args,
) -> dict:
    X_train = surface["X_model"][training_indices]
    X_val = surface["X_model"][validation_indices]
    y_train = preprocess["Y_train_target"][training_indices]
    y_val = preprocess["Y_train_target"][validation_indices]
    err_train = preprocess["E_train_target"][training_indices]
    err_val = preprocess["E_train_target"][validation_indices]
    fitting.print_table(
        "Neural training split:",
        ["Array", "Shape"],
        [
            ["X_train", X_train.shape],
            ["X_val", X_val.shape],
            ["y_train", y_train.shape],
            ["y_val", y_val.shape],
        ],
    )
    model_param = neural_model_parameters(
        in_dim=X_train.shape[1],
        out_dim=y_train.shape[1],
        args=args,
    )
    selected_model, full_model, stats = neural.train_validated_lzmlp(
        model_param=model_param,
        training=(X_train, y_train, err_train),
        validation=(X_val, y_val, err_val),
        full=(
            surface["X_model"],
            preprocess["Y_train_target"],
            preprocess["E_train_target"],
        ),
        options=surrogate_cli.neural_training_options(args),
        use_errors=args.mc_errors,
    )
    print_neural_validation_summary(stats)
    selected_model.eval()
    with torch.inference_mode():
        validation_prediction = (
            selected_model(torch.as_tensor(X_val, dtype=torch.float64, device=args.device))
            .cpu()
            .numpy()
        )
    stats.update(validation_statistics(y_val, validation_prediction))

    full_epochs = int(stats.get("best_epoch") or 0)
    residual_tensor = torch.as_tensor(
        stats["residual_scale"],
        dtype=torch.float64,
        device=args.device,
    )

    # Evaluate the full-refit network and optional held-out residual variance
    def predict_scaled(X_model: torch.Tensor, include_variance: bool):
        mean, variance = predict_scaled_neural_torch(
            full_model,
            residual_tensor,
            X_model.to(dtype=torch.float64),
        )
        return mean, variance if include_variance else torch.zeros_like(mean)

    checkpoint = save_neural_checkpoint(
        model=full_model,
        model_param=model_param,
        surface=surface,
        preprocess=preprocess,
        stats=stats,
        args=args,
    )
    return {
        "predict_scaled_torch": predict_scaled,
        "dtype": torch.float64,
        "checkpoint": checkpoint,
        "validation_prediction_scaled": validation_prediction,
        "summary": {
            "neural_validation": stats,
            "full_refit_epochs": full_epochs,
        },
    }


# Dispatch one fully validated and refitted replica backend
def build_replica_backend(
    surface: dict,
    preprocess: dict,
    training_indices: np.ndarray,
    validation_indices: np.ndarray,
    args,
) -> dict:
    builders = {
        "rff": build_rff_replica_backend,
        "gp": build_gp_replica_backend,
        "neural": build_neural_replica_backend,
    }
    if args.model not in builders:
        raise ValueError(f"Unknown model = {args.model}")
    return builders[args.model](surface, preprocess, training_indices, validation_indices, args)


# Write a compact JSON summary for the run
def save_summary(
    surface: dict,
    best_random: tuple[np.ndarray, float],
    best_opt: tuple[np.ndarray, float],
    parameter_uncertainties: np.ndarray,
    uncertainty_method: str | None,
    args,
    extra: dict | None = None,
) -> pathlib.Path:
    output_dir = pathlib.Path(args.output_root)
    ensure_dir(output_dir)
    summary_path = output_dir / "summary.json"

    best_trial_index = int(np.argmin(surface["Z"]))
    optimized_values = {
        name: float(best_opt[0][i]) for i, name in enumerate(surface["param_names"])
    }
    optimized_uncertainties = {
        name: parameter_uncertainties[i] for i, name in enumerate(surface["param_names"])
    }
    summary = {
        "replica_schema_version": 1,
        "output_run_id": pathlib.Path(args.output_root).name,
        "run_manifest": str(pathlib.Path(args.output_root) / "manifest.json"),
        "best_fit": fit_summary.build_best_fit(
            values=optimized_values,
            uncertainties=optimized_uncertainties,
            objective_name="Z",
            objective_value=best_opt[1],
            source="surrogate",
            uncertainty_method=uncertainty_method,
        ),
        "run_name": args.run_name,
        "simdriver": surface.get("simdriver"),
        "plot_brand": args.plot_brand,
        "surrogate_model": replica_backend_name(args.model),
        "objective_convention": {
            "symbol": "Z",
            "definition": (
                "fit-weighted generalized chi2"
                if surface["covariance_mode"] == "full"
                else "fit-weighted sum of informative-bin chi2 terms"
            ),
            "covariance_mode": surface["covariance_mode"],
        },
        "quality_mode": args.quality_mode,
        "quality_percentile": args.quality_percentile,
        "max_delta_Z": args.max_delta_Z,
        "max_Z": args.max_Z,
        "n_trials": int(len(surface["Z"])),
        "n_bins": int(surface["Y_hat"].shape[1]),
        "param_names": surface["param_names"],
        "parameter_topology": surface["parameter_topology"],
        "model_feature_names": surface["model_feature_names"],
        "realized_best": {
            "Z": float(surface["Z"][best_trial_index]),
            "theta": {
                name: float(surface["X"][best_trial_index, i])
                for i, name in enumerate(surface["param_names"])
            },
        },
        "random_best": {
            "Z": float(best_random[1]),
            "theta": {
                name: float(best_random[0][i]) for i, name in enumerate(surface["param_names"])
            },
        },
        "optimized_best": {
            "Z": float(best_opt[1]),
            "theta": optimized_values,
        },
    }
    if extra is not None:
        summary.update(json_safe(extra))

    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(json_safe(summary), handle, indent=4)
    return summary_path


# Build NumPy and scalar Torch views of the replica chi2 objective
def build_replica_objective(
    predict_histograms_torch, preprocess_torch: dict, device: str, dtype: torch.dtype
):
    # Evaluate a batch of physical coordinates on the active Torch device
    def surrogate_Z_torch(X_physical: torch.Tensor) -> torch.Tensor:
        counts, errors, _, _ = predict_histograms_torch(
            X_physical,
            include_model_variance=False,
        )
        return replica_Z_torch(counts, errors, preprocess_torch)

    # Evaluate a NumPy batch for shared plotting and posterior machinery
    def surrogate_Z(X_physical: np.ndarray) -> np.ndarray:
        values = np.asarray(X_physical, dtype=np.float64)
        values = np.atleast_2d(values)
        with torch.inference_mode():
            prediction = surrogate_Z_torch(torch.as_tensor(values, dtype=dtype, device=device))
        return prediction.detach().cpu().numpy()

    # Evaluate one scalar physical-coordinate objective for Torch autograd
    def scalar_objective(values: torch.Tensor) -> torch.Tensor:
        return surrogate_Z_torch(values[None, :])[0]

    return surrogate_Z_torch, surrogate_Z, scalar_objective


# Render direct fixed or profiled histogram impacts from the replica backend
def render_replica_impacts(
    *,
    m,
    best_fit: np.ndarray,
    covariance: np.ndarray,
    errors: np.ndarray,
    surface: dict,
    predict_histograms_torch,
    preprocess_torch: dict,
    surrogate_Z,
    value_gradient,
    torch_config: dict | None,
    validation_r2: np.ndarray,
    validation_mae: np.ndarray,
    output_root: pathlib.Path,
    args,
) -> dict | None:
    if m is None or not np.any(np.isfinite(errors)):
        cprint("Skipping impacts because fitted Hessian errors are unavailable", "yellow")
        return None
    if args.impact_mode == "fixed":
        shifts = impact.fixed_parameter_shifts(
            best_fit=best_fit,
            errors=errors,
            bounds=surface["bounds"],
        )
    else:
        shifts = impact.profiled_parameter_shifts(
            func_hat=surrogate_Z,
            best_fit=best_fit,
            errors=errors,
            bounds=surface["bounds"],
            covariance=covariance,
            maxiter=int(args.impact_profile_maxiter),
            value_gradient_func=value_gradient,
            torch_config=torch_config,
            workers=int(args.profile_workers),
        )

    # Evaluate direct per-observable chi2 values for all impact points
    def evaluate_observables(points: np.ndarray) -> np.ndarray:
        with torch.inference_mode():
            counts, mc_errors, _, _ = predict_histograms_torch(
                torch.as_tensor(
                    points,
                    dtype=preprocess_torch["data_counts"].dtype,
                    device=preprocess_torch["data_counts"].device,
                ),
                include_model_variance=False,
            )
            values = observable_chi2_torch(
                counts,
                mc_errors,
                preprocess_torch,
                surface["histogram_manifest"],
            )
        return values.cpu().numpy()

    return impact.render_direct_parameter_impacts(
        evaluate_observable_chi2=evaluate_observables,
        observables=surface["histogram_manifest"],
        param_names=surface["param_names"],
        best_fit=best_fit,
        bounds=surface["bounds"],
        shifts=shifts,
        output_dir=pathlib.Path(output_root) / "optimizer" / "impacts" / args.impact_mode,
        physical={"transform": surrogate_cli.parameter_transform(args, surface), "best_fit": best_fit,
                  "covariance": covariance, "bounds": surface["bounds"], "mode": args.impact_mode,
                  "output_dir": pathlib.Path(output_root) / "physical" / "impacts" / args.impact_mode,
                  "maxiter": args.impact_profile_maxiter, "value_gradient_func": value_gradient},
        validation_r2=validation_r2,
        validation_mae=validation_mae,
        trial_count=len(surface["Z"]),
        top_n=int(args.impact_top),
        plot_brand=args.plot_brand,
    )


# Load, clean and quality-filter the full replica trial surface
def prepare_replica_surface(args) -> dict:
    surface = load_replica_surface(args)
    surface, _ = surrogate_cli.quality_filter_surface(
        surface,
        finite_keys=("X", "X_model", "Y_hat", "E_hat", "Z"),
        row_keys=("X", "X_model", "Y_hat", "E_hat", "Z", "Z_err"),
        args=args,
    )
    return surface


# Run device scanning and the optional torch L-BFGS refinement
def optimize_replica_surface(
    *,
    surface: dict,
    surrogate_Z_torch,
    scalar_objective,
    surrogate_Z,
    value_gradient,
    hessian,
    dtype: torch.dtype,
    output_root: pathlib.Path,
    args,
) -> dict:
    best_random = random_search(
        surrogate_Z_torch=surrogate_Z_torch,
        bounds=surface["bounds"],
        args=args,
        dtype=dtype,
    )
    best_trial_index = int(np.argmin(surface["Z"]))
    realized_best = surface["X"][best_trial_index]
    realized_Z = float(surface["Z"][best_trial_index])
    torch_config = surrogate_cli.torch_minimizer_config(
        args=args,
        objective_torch=scalar_objective,
        objective_batch_torch=surrogate_Z_torch,
        device=args.device,
        dtype=dtype,
    )
    if args.optimize:
        cprint(f"Starting replica {args.optimizer_backend} refinement", "yellow")
        m, surrogate_best = surrogate_cli.optimize_surrogate(
            args=args,
            X=surface["X"],
            Z=surface["Z"],
            param_names=surface["param_names"],
            func_hat=surrogate_Z,
            bounds=surface["bounds"],
            output_root=output_root,
            realized_best=realized_best,
            realized_Z=realized_Z,
            x0=best_random[0],
            objective_torch=scalar_objective,
            objective_batch_torch=surrogate_Z_torch,
            device=args.device,
            dtype=dtype,
            transform=surrogate_cli.parameter_transform(args, surface),
        )
        best_opt = (surrogate_best, float(m.fmin.fval))
    else:
        cprint("Skipping replica surrogate refinement; using random-search best", "yellow")
        m = None
        best_opt = (
            np.asarray(best_random[0], dtype=np.float64),
            fitting.evaluate_scalar(surrogate_Z, best_random[0]),
        )
    fitting.print_table(
        "Replica optimization summary:",
        ["Candidate", "Z"],
        [
            ["random search", best_random[1]],
            ["surrogate optimum", best_opt[1]],
        ],
        align=["left", "right"],
    )
    _, parameter_uncertainties = fitting.fitted_covariance(
        m,
        surface["param_names"],
    )
    return {
        "m": m,
        "best_random": best_random,
        "best_opt": best_opt,
        "parameter_uncertainties": parameter_uncertainties,
        "uncertainty_method": (
            None
            if m is None
            else surrogate_cli.optimizer_uncertainty_method(args)
        ),
        "realized_best": realized_best,
        "realized_Z": realized_Z,
        "torch_config": torch_config,
    }


# Produce optional profiles, impacts, posterior and predictive histograms
def render_replica_postfit(
    *,
    fit: dict,
    surface: dict,
    predict_histograms_torch,
    preprocess_torch: dict,
    surrogate_Z,
    scalar_objective,
    value_gradient,
    validation_r2: np.ndarray,
    validation_mae: np.ndarray,
    dtype: torch.dtype,
    backend_name: str,
    output_root: pathlib.Path,
    args,
) -> dict:
    surrogate_cli.render_profiles(
        args=args,
        X=surface["X"],
        Z=surface["Z"],
        param_names=surface["param_names"],
        func_hat=surrogate_Z,
        output_root=output_root,
        realized_best=fit["realized_best"],
        surrogate_best=fit["best_opt"][0],
        realized_Z=fit["realized_Z"],
        surrogate_Z=float(fit["best_opt"][1]),
        bounds=surface["bounds"],
        value_gradient_func=value_gradient,
        torch_config=fit["torch_config"],
        transform=surrogate_cli.parameter_transform(args, surface),
    )
    covariance, errors = fitting.fitted_covariance(
        fit["m"],
        surface["param_names"],
    )
    impacts = render_replica_impacts(
        m=fit["m"],
        best_fit=fit["best_opt"][0],
        covariance=covariance,
        errors=errors,
        surface=surface,
        predict_histograms_torch=predict_histograms_torch,
        preprocess_torch=preprocess_torch,
        surrogate_Z=surrogate_Z,
        value_gradient=value_gradient,
        torch_config=fit["torch_config"],
        validation_r2=validation_r2,
        validation_mae=validation_mae,
        output_root=output_root,
        args=args,
    )
    surrogate_cli.run_posterior(
        args=args,
        objective_torch=scalar_objective,
        device=args.device,
        dtype=dtype,
        bounds=surface["bounds"],
        param_names=surface["param_names"],
        output_root=output_root,
        start=fit["best_opt"][0],
        transform=surrogate_cli.parameter_transform(args, surface),
    )
    predictive = None
    if args.predictive_plots:
        prediction = predictive_histogram_summary(
            predict_histograms_torch=predict_histograms_torch,
            point=fit["best_opt"][0],
            covariance=covariance,
            device=args.device,
            dtype=dtype,
        )
        predictive = render_predictive_histograms(
            surface=surface,
            prediction=prediction,
            point=fit["best_opt"][0],
            backend_name=backend_name,
            output_root=output_root,
            args=args,
        )
        cprint(
            f"Predictive histogram output: {predictive['figure_root']}",
            "yellow",
        )
    return {
        "impacts": impacts,
        "predictive_histograms": predictive,
    }


# Main executable path
def main():
    args = parse_args()
    start_time = surrogate_cli.start_run(
        args,
        tool="iceproxy",
        version=__version__,
        command=sys.argv,
    )

    surface = prepare_replica_surface(args)
    surrogate_cli.record_surface(
        args,
        surface,
        inputs={
            "mode": "full",
            "run_directory": str(pathlib.Path(args.cdir) / "runs" / "icetune" / args.run_name),
            "plot_template_source": surface["plot_template_source"],
            "simdriver": surface.get("simdriver"),
            "covariance_mode": surface["covariance_mode"],
        },
        histogram_bin_count=int(surface["Y_hat"].shape[1]),
    )
    best_trial_index = int(np.argmin(surface["Z"]))
    training_indices, validation_indices = training_validation_indices(
        n_rows=len(surface["Z"]),
        validation_fraction=args.validation_fraction,
        rngseed=args.rngseed,
        required_training=[best_trial_index],
    )
    preprocess = prepare_targets(
        Y_hat=surface["Y_hat"],
        E_hat=surface["E_hat"],
        Y_data=surface["Y_data"],
        E_data=surface["E_data"],
        fit_weights=surface["fit_weights"],
        data_total_covariance=surface["data_total_covariance"],
        covariance_indices=surface["covariance_indices"],
        training_indices=training_indices,
    )
    backend = build_replica_backend(
        surface=surface,
        preprocess=preprocess,
        training_indices=training_indices,
        validation_indices=validation_indices,
        args=args,
    )
    cprint(
        f"Saved replica surrogate checkpoint to: {backend['checkpoint']}",
        "yellow",
    )
    dtype = backend["dtype"]
    preprocess_torch = torch_preprocess(
        preprocess=preprocess,
        likelihood_mask=surface["likelihood_mask"],
        device=args.device,
        dtype=dtype,
    )
    predict_histograms_torch = make_histogram_predictor(
        predict_scaled_torch=backend["predict_scaled_torch"],
        surface=surface,
        preprocess_torch=preprocess_torch,
    )
    surrogate_Z_torch, surrogate_Z, scalar_objective = build_replica_objective(
        predict_histograms_torch=predict_histograms_torch,
        preprocess_torch=preprocess_torch,
        device=args.device,
        dtype=dtype,
    )
    value_gradient, hessian = fitting.build_torch_derivative_callbacks(
        objective=scalar_objective,
        device=args.device,
        dtype=dtype,
    )
    validation_r2, validation_mae = observable_validation_metrics_from_scaled(
        prediction_scaled=backend["validation_prediction_scaled"],
        Y_validation=surface["Y_hat"][validation_indices],
        E_validation=surface["E_hat"][validation_indices],
        preprocess_torch=preprocess_torch,
        manifest=surface["histogram_manifest"],
    )
    backend_name = replica_backend_name(args.model)
    output_root = pathlib.Path(args.output_root) / backend_name
    fit = optimize_replica_surface(
        surface=surface,
        surrogate_Z_torch=surrogate_Z_torch,
        scalar_objective=scalar_objective,
        surrogate_Z=surrogate_Z,
        value_gradient=value_gradient,
        hessian=hessian,
        dtype=dtype,
        output_root=output_root,
        args=args,
    )
    postfit = render_replica_postfit(
        fit=fit,
        surface=surface,
        predict_histograms_torch=predict_histograms_torch,
        preprocess_torch=preprocess_torch,
        surrogate_Z=surrogate_Z,
        scalar_objective=scalar_objective,
        value_gradient=value_gradient,
        validation_r2=validation_r2,
        validation_mae=validation_mae,
        dtype=dtype,
        backend_name=backend_name,
        output_root=output_root,
        args=args,
    )

    summary_extra = {
        **backend["summary"],
        "training_split": {
            "training_count": int(len(training_indices)),
            "validation_count": int(len(validation_indices)),
            "validation_indices": validation_indices.tolist(),
        },
        "observable_validation": {
            "r2": validation_r2,
            "mae_chi2": validation_mae,
        },
        "derivatives": {
            "gradient": "torch autograd",
            "hessian": "torch autograd",
            "device": str(args.device),
            "dtype": str(dtype),
        },
        **postfit,
    }
    summary_path = save_summary(
        surface=surface,
        best_random=fit["best_random"],
        best_opt=fit["best_opt"],
        parameter_uncertainties=fit["parameter_uncertainties"],
        uncertainty_method=fit["uncertainty_method"],
        args=args,
        extra=summary_extra,
    )
    cprint(f"Saved replica summary to: {summary_path}", "yellow")
    surrogate_cli.complete_run(
        args,
        tool="iceproxy",
        started=start_time,
        model=backend_name,
    )


if __name__ == "__main__":
    main()

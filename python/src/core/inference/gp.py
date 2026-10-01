# Exact scalar and multi-output Gaussian processes with input uncertainty
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import math
import time
from dataclasses import dataclass

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as functional
import torch.optim as optim
from torch import Tensor
from tqdm import tqdm

from core.numerics import lbfgsb
from core.numerics.linalg import solve_triangular

# Share one local ISO timestamp across surrogate trainers
from core.tune.runtime.process import training_timestamp

__all__ = ["ExactGaussianProcess", "GaussianPosterior", "exact_gp_defaults", "robust_cholesky"]


# Resolve the GP arithmetic precision with float64 as the stable default
def resolve_gp_dtype(device: str, dtype: str | torch.dtype | None) -> torch.dtype:
    if isinstance(dtype, torch.dtype):
        return dtype
    dtype_name = "auto" if dtype is None else str(dtype).lower()
    if dtype_name == "auto":
        return torch.float64
    if dtype_name == "float32":
        return torch.float32
    if dtype_name == "float64":
        return torch.float64
    raise ValueError("GP dtype must be auto, float32 or float64")


# Compute pairwise squared distances without cancellation and with double-backward support
def pairwise_squared_distance(X1: Tensor, X2: Tensor) -> Tensor:
    return (X1[:, None, :] - X2[None, :, :]).square().sum(dim=-1)


# Compute Matern five halves correlation with analytic curvature at coincidence
def matern52_correlation(squared_distance: Tensor) -> Tensor:
    scaled_squared = 5.0 * squared_distance
    scaled = torch.sqrt(torch.clamp(scaled_squared, min=1.0e-8))
    correlation = (1.0 + scaled + scaled.square() / 3.0) * torch.exp(-scaled)
    expansion = 1.0 - scaled_squared / 6.0 + scaled_squared.square() / 24.0
    return torch.where(scaled_squared < 1.0e-8, expansion, correlation)


# 1. Multi-output RBF kernel with input uncertainty


class MultiOutputKernel(nn.Module):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        lengthscale: float = 1.0,
        variance: float = 1.0,
        kind: str = "RBF",
        dtype: torch.dtype = torch.float64,
    ) -> None:
        """
        Multi-output RBF kernel with support for input uncertainties.
        The output covariance is represented via a Cholesky factor.
        """
        super().__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        self.kind = str(kind).upper()
        if self.kind not in {"RBF", "MATERN"}:
            raise ValueError("GP kernel kind must be RBF or Matern")

        # Optimizable positive kernel parameters
        self.lengthscale = nn.Parameter(torch.tensor(lengthscale, dtype=dtype))
        self.variance = nn.Parameter(torch.tensor(variance, dtype=dtype))

        # Lower-triangular matrix for output covariance (Cholesky factor)
        self.L = nn.Parameter(torch.eye(output_dim, dtype=dtype))

    # Evaluate the selected stationary kernel from squared distances
    def from_squared_distance(self, squared_distance: Tensor) -> Tensor:
        if self.kind == "RBF":
            return self.variance * torch.exp(-0.5 * squared_distance / self.lengthscale.square())
        return self.variance * matern52_correlation(squared_distance / self.lengthscale.square())

    def base_kernel(self, X1: Tensor, X2: Tensor, X1_err: Tensor | None = None, X2_err: Tensor | None = None) -> Tensor:
        """
        Computes the selected stationary kernel with optional input uncertainty.
        """
        if X1_err is None and X2_err is None:
            squared_distance = pairwise_squared_distance(X1, X2)
            return self.from_squared_distance(squared_distance)
        if self.kind != "RBF":
            raise NotImplementedError("Input-uncertainty integration is implemented only for the RBF kernel")

        X1_err = torch.zeros_like(X1) if X1_err is None else X1_err
        X2_err = torch.zeros_like(X2) if X2_err is None else X2_err

        # Scale inputs and uncertainties
        X1_scaled = X1 / self.lengthscale
        X2_scaled = X2 / self.lengthscale
        X1_err_scaled = X1_err / self.lengthscale
        X2_err_scaled = X2_err / self.lengthscale

        # Compute pairwise differences and total variance from uncertainties
        # diff: shape (N, M, D_in)
        diff = X1_scaled[:, None, :] - X2_scaled[None, :, :]
        total_var = (X1_err_scaled[:, None, :] ** 2) + (X2_err_scaled[None, :, :] ** 2)

        # RBF exponent with uncertainty adjustment
        exponent = -0.5 * torch.sum(diff**2 / (1 + total_var), dim=-1)
        # Jacobian adjustment from input uncertainty
        det_term = -0.5 * torch.sum(torch.log(1 + total_var), dim=-1)

        return self.variance * torch.exp(exponent + det_term)

    # Expand one input kernel into the multi-output block structure
    def from_base_kernel(self, K_input: Tensor) -> Tensor:
        output_cov = self.L @ self.L.T
        return (K_input[:, None, :, None] * output_cov[None, :, None, :]).permute(0, 2, 1, 3)

    def forward(self, X1: Tensor, X2: Tensor, X1_err: Tensor | None = None, X2_err: Tensor | None = None) -> Tensor:
        """
        Compute the multi-output kernel matrix.
        Returns a 4D tensor of shape (N, M, output_dim, output_dim).
        """
        K_input = self.base_kernel(X1, X2, X1_err, X2_err)
        return self.from_base_kernel(K_input)


# 2. Utility to reshape 4D kernel to 2D block matrix


def blockify_kernel(K_4d: Tensor) -> Tensor:
    N, M, D_out, _ = K_4d.shape
    return K_4d.permute(0, 2, 1, 3).reshape(N * D_out, M * D_out)


# 3. Robust Cholesky with adaptive jitter


def robust_multioutput_cholesky(K: Tensor, max_attempts: int = 5, init_jitter: float = 1e-6) -> Tensor:
    identity = torch.eye(K.shape[0], dtype=K.dtype, device=K.device)
    jitter = float(init_jitter)
    for _ in range(max_attempts):
        L, info = torch.linalg.cholesky_ex(K + jitter * identity)
        if int(torch.max(info).detach().cpu()) == 0:
            return L
        jitter *= 10.0

    # Shift the spectrum and return a genuine triangular Cholesky factor
    min_eigenvalue = float(torch.min(torch.linalg.eigvalsh(K)).detach().cpu())
    spectral_jitter = max(jitter, float(init_jitter) - min_eigenvalue)
    return torch.linalg.cholesky(K + spectral_jitter * identity)


# 4. GaussianProcess class


class GaussianProcess(nn.Module):
    def __init__(
        self,
        X_train: Tensor,
        Y_train: Tensor,
        X_err: Tensor | None = None,
        Y_err: Tensor | None = None,
        scale: float = 0.25,
        noise: float = 1e-6,
        num_inducing: int | None = None,
        device: str | None = None,
        dtype: str | torch.dtype | None = None,
        kernel: str = "RBF",
        independent_outputs: bool = False,
    ) -> None:
        """
        Multi-output Gaussian Process Regression with a stationary kernel.
        Independent-output mode shares one N by N kernel across target columns,
        while correlated mode uses the original output-block covariance.
        """
        super().__init__()
        self.device = device if device is not None else ("cuda" if torch.cuda.is_available() else "cpu")
        if str(self.device).startswith("cuda") and not torch.cuda.is_available():
            raise RuntimeError("CUDA was requested for the GP but no CUDA device is available")
        self.dtype = resolve_gp_dtype(self.device, dtype)
        self.independent_outputs = bool(independent_outputs)
        if self.independent_outputs and Y_err is not None:
            raise NotImplementedError("Independent-output GP currently requires Y_err=None")
        self.X_train = self._model_tensor(X_train, copy=True)
        self.Y_train = self._model_tensor(Y_train, copy=True)
        self.X_err = self._model_tensor(X_err, copy=True)
        self.Y_err = self._model_tensor(Y_err, copy=True)

        self.D_in = self.X_train.shape[1]
        self.N, self.D_out = self.Y_train.shape
        self.num_inducing = num_inducing

        # Initialize the multi-output kernel
        self.kernel_module = MultiOutputKernel(
            input_dim=self.D_in,
            output_dim=1 if self.independent_outputs else self.D_out,
            lengthscale=scale,
            variance=1.0,
            kind=kernel,
            dtype=self.dtype,
        ).to(self.device)

        # Noise parameter added to the diagonal
        self.noise = nn.Parameter(torch.tensor(noise, dtype=self.dtype, device=self.device))
        self._refresh_training_cache()

        device_name = str(self.device)
        if str(self.device).startswith("cuda"):
            device_name = f"{self.device} ({torch.cuda.get_device_name(self.device)})"
        print(f"{training_timestamp()} GaussianProcess: device = {device_name}, dtype = {self.dtype}")
        print(
            f"{training_timestamp()} GaussianProcess: "
            f"N = {self.N}, D_in = {self.D_in}, D_out = {self.D_out}, "
            f"kernel = {self.kernel_module.kind}, "
            f"independent_outputs = {self.independent_outputs}"
        )

        self.update_covariance()

    # Move one optional tensor to model arithmetic while optionally detaching a private copy
    def _model_tensor(self, value: Tensor | None, *, copy: bool = False) -> Tensor | None:
        if value is None:
            return None
        tensor = value.detach().clone() if copy else value
        return tensor.to(device=self.device, dtype=self.dtype)

    # Refresh flattened targets and cached training distances
    def _refresh_training_cache(self) -> None:
        self.N, self.D_out = self.Y_train.shape
        self.Y_train_flat = self.Y_train if self.independent_outputs else self.Y_train.reshape(self.N * self.D_out, 1)
        self.Y_err_flat = self.Y_err.square().reshape(self.N * self.D_out) if self.Y_err is not None else None
        self._train_squared_distance = (
            pairwise_squared_distance(self.X_train, self.X_train) if self.X_err is None else None
        )

    # Replace training data while preserving selected hyperparameters
    def set_training_data(
        self, X_train: Tensor, Y_train: Tensor, X_err: Tensor | None = None, Y_err: Tensor | None = None
    ) -> None:
        X_new = self._model_tensor(X_train, copy=True)
        Y_new = self._model_tensor(Y_train, copy=True)
        if X_new.ndim != 2 or Y_new.ndim != 2 or len(X_new) != len(Y_new):
            raise ValueError("GP training inputs and targets must be aligned matrices")
        if X_new.shape[1] != self.D_in or Y_new.shape[1] != self.D_out:
            raise ValueError("Replacement GP training dimensions do not match the model")
        if self.independent_outputs and Y_err is not None:
            raise NotImplementedError("Independent-output GP currently requires Y_err=None")
        self.X_train = X_new
        self.Y_train = Y_new
        self.X_err = self._model_tensor(X_err, copy=True)
        self.Y_err = self._model_tensor(Y_err, copy=True)
        self._refresh_training_cache()
        self.update_covariance()

    # Build the training kernel from cached distances when possible
    def _training_kernel(self) -> Tensor:
        if self.independent_outputs:
            if self._train_squared_distance is None:
                return self.kernel_module.base_kernel(self.X_train, self.X_train, self.X_err, self.X_err)
            return self.kernel_module.from_squared_distance(self._train_squared_distance)
        if self._train_squared_distance is None:
            return blockify_kernel(self.kernel_module(self.X_train, self.X_train, self.X_err, self.X_err))
        K_input = self.kernel_module.from_squared_distance(self._train_squared_distance)
        return blockify_kernel(self.kernel_module.from_base_kernel(K_input))

    # Build the noisy symmetric covariance for the current training state
    def _training_covariance(self) -> Tensor:
        K = self._training_kernel()
        diag_noise = torch.zeros(K.shape[0], dtype=K.dtype, device=K.device)
        if self.Y_err_flat is not None:
            diag_noise += self.Y_err_flat
        diag_noise += self.noise.expand(K.shape[0])
        K = K + torch.diag(diag_noise)
        return 0.5 * (K + K.T)

    def update_covariance(self) -> None:
        """Compute the full training covariance matrix and perform a robust Cholesky decomposition."""
        if self.num_inducing is not None:
            raise NotImplementedError("Sparse GP is not implemented for multi-output block matrices.")
        self.K = self._training_covariance()
        self.L = robust_multioutput_cholesky(self.K)
        self.alpha = torch.cholesky_solve(self.Y_train_flat, self.L)

    # Freeze hyperparameters and detach posterior state for input autograd
    def freeze_for_inference(self) -> None:
        for parameter in self.parameters():
            parameter.requires_grad_(False)
        self.K = self.K.detach()
        self.L = self.L.detach()
        self.alpha = self.alpha.detach()

    def compute_cross_cov(self, X_test: Tensor, X_test_err: Tensor | None = None) -> Tensor:
        """Compute the cross-covariance between test inputs and training inputs."""
        return self._covariance(X_test, self.X_train, X_test_err, self.X_err)

    def compute_test_cov(self, X_test: Tensor, X_test_err: Tensor | None = None) -> Tensor:
        """Compute the covariance matrix for test inputs."""
        return self._covariance(X_test, X_test, X_test_err, X_test_err)

    # Build one independent or output-block covariance from the shared input kernel
    def _covariance(self, X1: Tensor, X2: Tensor, X1_err: Tensor | None, X2_err: Tensor | None) -> Tensor:
        base = self.kernel_module.base_kernel(X1, X2, X1_err, X2_err)
        return base if self.independent_outputs else blockify_kernel(self.kernel_module.from_base_kernel(base))

    # Compute the input-kernel variance at each test point
    def _input_variance(self, X_test: Tensor, X_test_err: Tensor | None) -> Tensor:
        if X_test_err is None:
            return self.kernel_module.variance.expand(len(X_test))
        return torch.diag(self.kernel_module.base_kernel(X_test, X_test, X_test_err, X_test_err))

    # Predict only the latent mean without constructing test covariance
    def predict_mean(self, X_test: Tensor, X_test_err: Tensor | None = None) -> Tensor:
        X_test = self._model_tensor(X_test)
        X_test_err = self._model_tensor(X_test_err)
        K_s = self.compute_cross_cov(X_test, X_test_err)
        if self.independent_outputs:
            return K_s @ self.alpha
        return (K_s @ self.alpha).reshape(X_test.shape[0], self.D_out)

    # Predict marginal means and variances without a dense test covariance
    def predict_marginals(
        self, X_test: Tensor, X_test_err: Tensor | None = None, Y_test_err: Tensor | None = None
    ) -> tuple[Tensor, Tensor]:
        X_test = self._model_tensor(X_test)
        X_test_err = self._model_tensor(X_test_err)
        Y_test_err = self._model_tensor(Y_test_err)

        K_s = self.compute_cross_cov(X_test, X_test_err)
        mean_flat = K_s @ self.alpha
        solved = torch.cholesky_solve(K_s.T, self.L)
        reduction = torch.sum(K_s * solved.T, dim=1)

        if self.independent_outputs:
            variance = torch.clamp(self._input_variance(X_test, X_test_err) + self.noise - reduction, min=1e-12)[
                :, None
            ].expand(-1, self.D_out)
            if Y_test_err is not None:
                variance = variance + Y_test_err.square()
            return mean_flat, variance

        output_variance = torch.diag(self.kernel_module.L @ self.kernel_module.L.T)
        input_variance = self._input_variance(X_test, X_test_err)
        prior_variance = (input_variance[:, None] * output_variance[None, :]).reshape(-1)
        prior_variance += self.noise
        if Y_test_err is not None:
            prior_variance += Y_test_err.square().reshape(-1)
        variance_flat = torch.clamp(prior_variance - reduction, min=1e-12)
        return (mean_flat.reshape(len(X_test), self.D_out), variance_flat.reshape(len(X_test), self.D_out))

    def forward(
        self, X_test: Tensor, X_test_err: Tensor | None = None, Y_test_err: Tensor | None = None
    ) -> (Tensor, Tensor):
        """
        Compute the predictive mean and variance for test inputs.
        Returns:
            mu_s:  Predictive mean of shape (N_test, D_out)
            var_s: Predictive variance of shape (N_test, D_out)
        """
        try:
            return self.predict_marginals(X_test, X_test_err=X_test_err, Y_test_err=Y_test_err)
        except RuntimeError:
            # Preserve a robust fallback for unusual linear algebra backends
            X_test = self._model_tensor(X_test)
            X_test_err = self._model_tensor(X_test_err)
            Y_test_err = self._model_tensor(Y_test_err)
            K_s = self.compute_cross_cov(X_test, X_test_err)
            K_ss = self.compute_test_cov(X_test, X_test_err)
            K_ss += self.noise * torch.eye(K_ss.shape[0], dtype=self.dtype, device=self.device)
            if Y_test_err is not None:
                K_ss += torch.diag(Y_test_err.square().reshape(-1))
            covariance = K_ss - K_s @ torch.linalg.solve(self.K, K_s.T)
            if self.independent_outputs:
                mean = K_s @ self.alpha
                variance = torch.clamp(torch.diag(covariance), min=1e-12)
                return mean, variance[:, None].expand(-1, self.D_out)
            mean = (K_s @ self.alpha).reshape(len(X_test), self.D_out)
            variance = torch.clamp(torch.diag(covariance), min=1e-12)
            return mean, variance.reshape(len(X_test), self.D_out)

    def predict(
        self, X_test: Tensor, X_test_err: Tensor | None = None, Y_test_err: Tensor | None = None
    ) -> (Tensor, Tensor):
        return self.forward(X_test, X_test_err, Y_test_err)

    # Compute a detached snapshot of the current trainable state
    def _state_snapshot(self) -> dict:
        return {key: value.detach().clone() for key, value in self.state_dict().items()}

    # Compute validation negative log predictive density for holdout monitoring
    def validation_loss(self, X_val: Tensor, Y_val: Tensor, Y_val_err: Tensor | None = None) -> Tensor:
        X_val = self._model_tensor(X_val)
        Y_val = self._model_tensor(Y_val)
        Y_val_err = self._model_tensor(Y_val_err)

        mean, var = self.forward(X_val, Y_test_err=Y_val_err)
        total_var = torch.clamp(var, min=1e-12)
        residual = Y_val - mean
        log_norm = torch.log(torch.tensor(2.0 * np.pi, dtype=total_var.dtype, device=self.device) * total_var)
        return torch.mean(0.5 * (residual**2 / total_var + log_norm))

    # Optimize GP hyperparameters with holdout-based model selection
    def optimize_hyperparameters(
        self,
        num_steps: int = 100,
        lr: float = 0.01,
        optimize_params: list | None = None,
        X_val: Tensor | None = None,
        Y_val: Tensor | None = None,
        Y_val_err: Tensor | None = None,
        patience: int | None = None,
        min_delta: float = 0.0,
        validation_interval: int = 1,
        log_interval: int = 20,
    ) -> dict:
        if validation_interval < 1:
            raise ValueError("GP validation interval must be at least 1")
        if log_interval < 1:
            raise ValueError("GP log interval must be at least 1")
        if optimize_params is None:
            # By default optimize the lengthscale, variance, and noise
            optimize_params = [self.kernel_module.lengthscale, self.kernel_module.variance, self.noise]
            if not self.independent_outputs and self.D_out > 1:
                optimize_params.append(self.kernel_module.L)
        optimizer = optim.Adam(optimize_params, lr=lr)
        has_validation = X_val is not None and Y_val is not None and X_val.shape[0] > 0
        best_state = None
        best_val_loss = float("inf")
        best_step = None
        optimization_start = time.perf_counter()
        stats = {
            "train_loss": [],
            "validation_loss": [],
            "validation_steps": [],
            "initial_validation_loss": None,
            "best_validation_loss": None,
            "best_step": None,
            "stopped_step": None,
            "elapsed_seconds": None,
            "device": str(self.device),
            "dtype": str(self.dtype),
        }

        if has_validation:
            with torch.no_grad():
                initial_val_loss = float(self.validation_loss(X_val, Y_val, Y_val_err).detach().cpu())
            stats["initial_validation_loss"] = initial_val_loss
            if np.isfinite(initial_val_loss):
                best_state = self._state_snapshot()
                best_val_loss = initial_val_loss
                best_step = -1

        for step in range(num_steps):
            optimizer.zero_grad()
            K = self._training_covariance()
            L = robust_multioutput_cholesky(K, max_attempts=3, init_jitter=1e-6)
            alpha = torch.cholesky_solve(self.Y_train_flat, L)

            log_det_factor = self.D_out if self.independent_outputs else 1
            log_likelihood = (
                -0.5 * torch.sum(self.Y_train_flat * alpha)
                - log_det_factor * torch.sum(torch.log(torch.diag(L)))
                - 0.5
                * (self.N * self.D_out)
                * torch.log(torch.tensor(2.0 * np.pi, dtype=self.dtype, device=self.device))
            )
            loss = -log_likelihood
            loss.backward()
            optimizer.step()

            val_loss = None
            with torch.no_grad():
                self.kernel_module.lengthscale.clamp_(1e-6, 1e6)
                self.kernel_module.variance.clamp_(1e-6, 1e6)
                self.noise.clamp_(1e-9, 1e6)
                should_validate = has_validation and (step % validation_interval == 0 or step == num_steps - 1)
                if should_validate:
                    self.update_covariance()
                    val_loss_tensor = self.validation_loss(X_val, Y_val, Y_val_err)
                    val_loss = float(val_loss_tensor.detach().cpu())
                    stats["validation_loss"].append(val_loss)
                    stats["validation_steps"].append(step)

                    if np.isfinite(val_loss) and val_loss < best_val_loss - min_delta:
                        best_state = self._state_snapshot()
                        best_val_loss = val_loss
                        best_step = step

            stats["train_loss"].append(float(loss.detach().cpu()))

            if step % log_interval == 0 or step == num_steps - 1:
                val_msg = f", val_nlpd = {val_loss:.4f}" if val_loss is not None else ""
                elapsed = time.perf_counter() - optimization_start
                print(
                    f"{training_timestamp()} Step {step}: "
                    f"elapsed = {elapsed:.1f} s, loss = {loss.item():.4f}"
                    f"{val_msg}, "
                    f"lengthscale = {self.kernel_module.lengthscale.item():.4f}, "
                    f"variance = {self.kernel_module.variance.item():.4f}, "
                    f"noise = {self.noise.item():.6f}"
                )

            waited_steps = step - best_step if best_step is not None else 0
            if (
                has_validation
                and val_loss is not None
                and patience is not None
                and patience > 0
                and waited_steps >= patience
            ):
                stats["stopped_step"] = step
                best_step_label = "initial" if best_step == -1 else str(best_step)
                print(
                    f"{training_timestamp()} Early stopping GP hyperparameter "
                    f"optimization at step {step}; best validation NLPD was "
                    f"{best_val_loss:.4f} at step {best_step_label}"
                )
                break

        if has_validation and best_state is not None:
            self.load_state_dict(best_state)
            stats["best_validation_loss"] = best_val_loss
            stats["best_step"] = best_step

        # Update the covariance with optimized or validation-selected hyperparameters
        self.update_covariance()
        stats["elapsed_seconds"] = time.perf_counter() - optimization_start
        print(
            f"{training_timestamp()} Optimized lengthscale: "
            f"{self.kernel_module.lengthscale.item():.4f}, "
            f"variance: {self.kernel_module.variance.item():.4f}, "
            f"noise: {self.noise.item():.6f}"
        )
        if has_validation and best_state is not None:
            best_step_label = "initial" if best_step == -1 else str(best_step)
            print(f"{training_timestamp()} Selected GP validation NLPD: {best_val_loss:.4f} at step {best_step_label}")
        return stats


@dataclass(frozen=True)
class GaussianPosterior:
    """Marginal latent and observation moments"""

    mean: torch.Tensor
    epistemic_variance: torch.Tensor
    aleatoric_variance: torch.Tensor

    # Compute the full predictive variance
    @property
    def predictive_variance(self) -> torch.Tensor:
        return self.epistemic_variance + self.aleatoric_variance


class EuclideanGeometry(nn.Module):
    """ARD Euclidean geometry for unconstrained model coordinates"""

    # Initialize a Euclidean input geometry
    def __init__(self, dimension: int) -> None:
        super().__init__()
        self.dimension = int(dimension)
        self.effective_dim = int(dimension)

    # Compute pairwise ARD squared distances
    def squared_distance(
        self, first: torch.Tensor, second: torch.Tensor, lengthscale: torch.Tensor, *, diag: bool = False
    ) -> torch.Tensor:
        scale = lengthscale.reshape(-1)
        if scale.numel() != self.dimension:
            raise ValueError("GP length scales do not match the input dimension")
        if diag:
            return ((first - second) / scale).square().sum(dim=-1)
        difference = first.unsqueeze(-2) - second.unsqueeze(-3)
        return (difference / scale).square().sum(dim=-1)


# Compute the default exact GP settings without hidden numerical choices
def exact_gp_defaults() -> dict:
    return {
        "kernel": "matern52",
        "fit_steps": 100,
        "fit_restarts": 3,
        "learning_rate": 1.0,
        "fit_tolerance_grad": 1.0e-7,
        "fit_tolerance_change": 1.0e-9,
        "gradient_clip": 20.0,
        "initial_lengthscale": 0.25,
        "initial_outputscale": 1.0,
        "initial_noise": 1.0e-4,
        "lengthscale_prior_median": 0.25,
        "lengthscale_prior_log_std": 1.0,
        "outputscale_prior_median": 1.0,
        "outputscale_prior_log_std": 1.0,
        "noise_prior_median": 1.0e-3,
        "noise_prior_log_std": 2.0,
        "cholesky_attempts": 7,
        "initial_jitter": 1.0e-8,
        "jitter_multiplier": 10.0,
        "posterior_variance_floor": 1.0e-12,
    }


class ExactGaussianProcess(nn.Module):
    """Exact scalar radial GP with known and learned noise"""

    # Initialize one exact scalar GP
    def __init__(
        self,
        X_train: torch.Tensor,
        Y_train: torch.Tensor,
        *,
        Y_err: torch.Tensor | None = None,
        geometry: nn.Module | None = None,
        settings: dict | None = None,
        device: str | torch.device = "cpu",
        dtype: torch.dtype = torch.float64,
        seed: int = 0,
    ) -> None:
        super().__init__()
        self.device = torch.device(device)
        self.dtype = dtype
        self.seed = int(seed)
        self.settings = exact_gp_defaults() if settings is None else copy.deepcopy(settings)
        if self.settings["kernel"] not in {"rbf", "matern52"}:
            raise ValueError("Exact GP kernel must be rbf or matern52")
        X = self._tensor(X_train)
        if X.ndim != 2:
            raise ValueError("Exact GP training inputs must be a matrix")
        self.geometry = copy.deepcopy(EuclideanGeometry(X.shape[1]) if geometry is None else geometry).to(
            device=self.device, dtype=self.dtype
        )
        if self.geometry.dimension != X.shape[1]:
            raise ValueError("Exact GP geometry does not match the training inputs")
        initial_lengthscale = torch.full(
            (self.geometry.effective_dim,),
            float(self.settings["initial_lengthscale"]),
            device=self.device,
            dtype=self.dtype,
        )
        self.raw_lengthscale = nn.Parameter(_inverse_softplus(initial_lengthscale))
        self.raw_outputscale = nn.Parameter(_inverse_softplus(self._scalar(self.settings["initial_outputscale"])))
        self.raw_noise = nn.Parameter(_inverse_softplus(self._scalar(self.settings["initial_noise"])))
        self.constant_mean = nn.Parameter(torch.zeros((), device=self.device, dtype=self.dtype))
        self.X = torch.empty((0, X.shape[1]), device=self.device, dtype=self.dtype)
        self.y = torch.empty((0,), device=self.device, dtype=self.dtype)
        self.known_variance = torch.empty((0,), device=self.device, dtype=self.dtype)
        self.cholesky: torch.Tensor | None = None
        self.alpha: torch.Tensor | None = None
        self.jitter = 0.0
        self.set_training_data(X, Y_train, Y_err=Y_err)

    # Convert values to model arithmetic
    def _tensor(self, values: torch.Tensor) -> torch.Tensor:
        return torch.as_tensor(values, device=self.device, dtype=self.dtype)

    # Create one scalar in model arithmetic
    def _scalar(self, value: float) -> torch.Tensor:
        return torch.tensor(float(value), device=self.device, dtype=self.dtype)

    # Compute positive ARD length scales
    @property
    def lengthscale(self) -> torch.Tensor:
        return functional.softplus(self.raw_lengthscale)

    # Compute the positive latent process variance
    @property
    def outputscale(self) -> torch.Tensor:
        return functional.softplus(self.raw_outputscale)

    # Compute the learned homoscedastic observation variance
    @property
    def noise(self) -> torch.Tensor:
        return functional.softplus(self.raw_noise)

    # Evaluate the configured radial covariance on the input geometry
    def kernel(self, first: torch.Tensor, second: torch.Tensor) -> torch.Tensor:
        squared = self.geometry.squared_distance(first, second, self.lengthscale)
        correlation = torch.exp(-0.5 * squared) if self.settings["kernel"] == "rbf" else matern52_correlation(squared)
        return self.outputscale * correlation

    # Replace training observations while preserving hyperparameters
    def set_training_data(
        self, X_train: torch.Tensor, Y_train: torch.Tensor, *, Y_err: torch.Tensor | None = None
    ) -> None:
        X = self._tensor(X_train).detach().clone()
        y = self._tensor(Y_train).reshape(-1).detach().clone()
        if X.ndim != 2 or y.shape != (len(X),):
            raise ValueError("Exact GP training inputs and targets are not aligned")
        if X.shape[1] != self.geometry.dimension:
            raise ValueError("Replacement GP inputs do not match the model dimension")
        errors = torch.zeros_like(y) if Y_err is None else self._tensor(Y_err).reshape(-1).detach().clone()
        if errors.shape != y.shape or not bool(torch.all(torch.isfinite(errors))):
            raise ValueError("Exact GP observation errors must be finite and aligned")
        if bool(torch.any(errors < 0.0)):
            raise ValueError("Exact GP observation errors must be nonnegative")
        if not bool(torch.all(torch.isfinite(X))) or not bool(torch.all(torch.isfinite(y))):
            raise ValueError("Exact GP training data must be finite")
        self.X = X
        self.y = y
        self.known_variance = errors.square()
        self._refresh()

    # Build the noisy training covariance
    def _training_covariance(self) -> torch.Tensor:
        covariance = self.kernel(self.X, self.X)
        diagonal = self.known_variance + self.noise
        return 0.5 * (covariance + covariance.T) + torch.diag(diagonal)

    # Refresh the exact posterior factorization
    def _refresh(self) -> None:
        covariance = self._training_covariance()
        self.cholesky, self.jitter = robust_cholesky(
            covariance,
            attempts=self.settings["cholesky_attempts"],
            initial_jitter=self.settings["initial_jitter"],
            jitter_multiplier=self.settings["jitter_multiplier"],
        )
        centered = (self.y - self.constant_mean).unsqueeze(-1)
        self.alpha = torch.cholesky_solve(centered, self.cholesky).squeeze(-1)

    # Compute one LogNormal log density including its normalization
    @staticmethod
    def _lognormal_log_prob(values: torch.Tensor, median: float, log_std: float) -> torch.Tensor:
        sigma = float(log_std)
        log_values = torch.log(values)
        standardized = (log_values - math.log(float(median))) / sigma
        return -0.5 * standardized.square() - log_values - math.log(sigma * math.sqrt(2.0 * math.pi))

    # Compute the negative log posterior of all GP hyperparameters
    def negative_log_posterior(self) -> torch.Tensor:
        covariance = self._training_covariance()
        cholesky, _ = robust_cholesky(
            covariance,
            attempts=self.settings["cholesky_attempts"],
            initial_jitter=self.settings["initial_jitter"],
            jitter_multiplier=self.settings["jitter_multiplier"],
        )
        centered = (self.y - self.constant_mean).unsqueeze(-1)
        alpha = torch.cholesky_solve(centered, cholesky)
        loss = 0.5 * torch.sum(centered * alpha)
        loss = loss + torch.log(torch.diag(cholesky)).sum()
        loss = loss + 0.5 * len(self.X) * math.log(2.0 * math.pi)
        loss = (
            loss
            - self._lognormal_log_prob(
                self.lengthscale, self.settings["lengthscale_prior_median"], self.settings["lengthscale_prior_log_std"]
            ).sum()
        )
        loss = loss - self._lognormal_log_prob(
            self.outputscale, self.settings["outputscale_prior_median"], self.settings["outputscale_prior_log_std"]
        )
        loss = loss - self._lognormal_log_prob(
            self.noise, self.settings["noise_prior_median"], self.settings["noise_prior_log_std"]
        )
        return loss

    # Initialize one MAP optimization restart from the configured priors
    def _initialize_restart(self, restart: int) -> None:
        if restart == 0:
            return
        generator = torch.Generator(device=self.device)
        generator.manual_seed(self.seed + 104729 * restart)

        # Draw one positive parameter from a configured LogNormal prior
        def draw(median_key: str, std_key: str, shape: tuple[int, ...]) -> torch.Tensor:
            normal = torch.randn(shape, generator=generator, device=self.device, dtype=self.dtype)
            return torch.exp(math.log(float(self.settings[median_key])) + float(self.settings[std_key]) * normal)

        self.raw_lengthscale.data.copy_(
            _inverse_softplus(
                draw("lengthscale_prior_median", "lengthscale_prior_log_std", (self.geometry.effective_dim,))
            )
        )
        self.raw_outputscale.data.copy_(
            _inverse_softplus(draw("outputscale_prior_median", "outputscale_prior_log_std", ()))
        )
        self.raw_noise.data.copy_(_inverse_softplus(draw("noise_prior_median", "noise_prior_log_std", ())))
        self.constant_mean.data.copy_(torch.mean(self.y))

    # Fit MAP hyperparameters with multistart initialization or warm continuation
    def fit(self, *, steps: int | None = None, restarts: int | None = None, warm_start: bool = False) -> dict:
        fit_steps = self.settings["fit_steps"] if steps is None else int(steps)
        fit_restarts = self.settings["fit_restarts"] if restarts is None else int(restarts)
        initial_state = self._parameter_snapshot()
        best_state = None
        best_loss = math.inf
        successful = 0
        converged = 0
        completed_steps = 0
        evaluations = gradients = 0
        parameters = tuple(self.parameters())
        started = time.perf_counter()

        # Evaluate the MAP objective and defer derivatives until line search acceptance
        def evaluate(values):
            nonlocal evaluations
            torch.nn.utils.vector_to_parameters(values.detach(), parameters)
            evaluations += 1
            try:
                with torch.enable_grad():
                    loss = self.negative_log_posterior()
            except torch.linalg.LinAlgError:
                return values.new_tensor(torch.inf), torch.zeros_like(values)

            # Differentiate the accepted covariance without repeating its factorization
            def derivative():
                nonlocal best_loss, best_state, gradients
                gradients += 1
                try:
                    gradient = torch.cat([part.reshape(-1) for part in torch.autograd.grad(loss, parameters)])
                except torch.linalg.LinAlgError:
                    return torch.full_like(values, torch.nan)
                numeric = float(loss.detach())
                if numeric < best_loss and bool(torch.isfinite(gradient).all()):
                    best_loss, best_state = numeric, self._parameter_snapshot()
                return gradient.detach()

            return loss.detach(), derivative

        for restart in range(max(1, fit_restarts)):
            self._restore_parameters(initial_state)
            self._initialize_restart(restart)
            start = torch.nn.utils.parameters_to_vector(parameters).detach().clone()
            _, loss, _, diagnostics = lbfgsb.minimize(
                None, start, evaluator=evaluate, bounds=None, maxiter=max(1, fit_steps),
                learning_rate=self.settings["learning_rate"],
                gradient_tolerance=self.settings["fit_tolerance_grad"],
                function_tolerance=self.settings["fit_tolerance_change"],
            )
            completed_steps += diagnostics["iterations"]
            successful += int(bool(torch.isfinite(loss)))
            converged += int(diagnostics["success"])
            # Continue a converged fit and retain random restarts as recovery
            if warm_start and diagnostics["success"]:
                break
        if best_state is None:
            raise RuntimeError("Exact GP MAP optimization did not produce a finite state")
        self._restore_parameters(best_state)
        self._refresh()
        return {
            "best_negative_log_posterior": best_loss,
            "steps": completed_steps,
            "restarts": restart + 1,
            "successful_restarts": successful,
            "converged_restarts": converged,
            "evaluations": evaluations,
            "gradient_evaluations": gradients,
            "elapsed_seconds": time.perf_counter() - started,
        }

    # Compute mean negative log predictive density on a holdout sample
    def validation_loss(
        self, X_val: torch.Tensor, Y_val: torch.Tensor, Y_val_err: torch.Tensor | None = None
    ) -> torch.Tensor:
        target = self._tensor(Y_val).reshape(-1)
        mean, variance = self.predict_marginals(X_val)
        total = variance.reshape(-1)
        if Y_val_err is not None:
            total = total + self._tensor(Y_val_err).reshape(-1).square()
        total = torch.clamp(total, min=self.settings["posterior_variance_floor"])
        residual = target - mean.reshape(-1)
        return 0.5 * torch.mean(residual.square() / total + torch.log(2.0 * math.pi * total))

    # Rebuild the posterior and evaluate one optional holdout loss
    def _holdout_loss(
        self, X_val: torch.Tensor | None, Y_val: torch.Tensor | None, Y_val_err: torch.Tensor | None
    ) -> float | None:
        if X_val is None or Y_val is None or len(X_val) == 0:
            return None
        with torch.no_grad():
            self._refresh()
            value = self.validation_loss(X_val, Y_val, Y_val_err)
        return float(value.detach().cpu())

    # Fit one multistart MAP realization with optional validation selection
    def _fit_validation_restart(
        self,
        *,
        initial_state: dict[str, torch.Tensor],
        restart: int,
        total_restarts: int,
        num_steps: int,
        learning_rate: float,
        X_val: torch.Tensor | None,
        Y_val: torch.Tensor | None,
        Y_val_err: torch.Tensor | None,
        patience: int | None,
        min_delta: float,
        validation_interval: int,
        log_interval: int,
    ) -> dict:
        self._restore_parameters(initial_state)
        self._initialize_restart(restart)
        optimizer = torch.optim.Adam(self.parameters(), lr=learning_rate)
        started = time.perf_counter()
        has_validation = X_val is not None and Y_val is not None and len(X_val) > 0
        initial_validation = self._holdout_loss(X_val, Y_val, Y_val_err)
        valid_initial = initial_validation is not None and math.isfinite(initial_validation)
        best_state = self._parameter_snapshot() if valid_initial else None
        best_loss = initial_validation if valid_initial else math.inf
        best_step = -1 if valid_initial else None
        train_loss = []
        validation_loss = []
        validation_steps = []
        stopped_step = None

        label = f"Exact GP restart {restart + 1}/{total_restarts}"
        with tqdm(range(num_steps), desc=label, unit="step", leave=False) as progress:
            for step in progress:
                optimizer.zero_grad(set_to_none=True)
                loss = self.negative_log_posterior()
                if not bool(torch.isfinite(loss)):
                    break
                numeric_loss = float(loss.detach().cpu())
                if not has_validation and numeric_loss < best_loss:
                    best_state = self._parameter_snapshot()
                    best_loss = numeric_loss
                    best_step = step - 1
                loss.backward()
                if any(p.grad is not None and not torch.isfinite(p.grad).all() for p in self.parameters()):
                    break
                torch.nn.utils.clip_grad_norm_(self.parameters(), self.settings["gradient_clip"])
                optimizer.step()
                train_loss.append(numeric_loss)

                validation = None
                should_validate = has_validation and (step % validation_interval == 0 or step == num_steps - 1)
                if should_validate:
                    validation = self._holdout_loss(X_val, Y_val, Y_val_err)
                    validation_loss.append(validation)
                    validation_steps.append(step)
                    if math.isfinite(validation) and validation < best_loss - min_delta:
                        best_state = self._parameter_snapshot()
                        best_loss = validation
                        best_step = step

                if step % log_interval == 0 or step == num_steps - 1:
                    fields = {"loss": f"{numeric_loss:.4g}"}
                    if validation is not None:
                        fields["val_nlpd"] = f"{validation:.4g}"
                    progress.set_postfix(fields)

                waited_steps = step - best_step if best_step is not None else 0
                if validation is not None and patience is not None and patience > 0 and waited_steps >= patience:
                    stopped_step = step
                    break

        with torch.no_grad():
            self._refresh()
            final_loss = float(self.negative_log_posterior().detach().cpu())
        if not has_validation and math.isfinite(final_loss) and final_loss < best_loss:
            best_state = self._parameter_snapshot()
            best_loss = final_loss
            best_step = len(train_loss) - 1
        successful = best_state is not None and math.isfinite(best_loss)
        elapsed = time.perf_counter() - started
        status = "early stop" if stopped_step is not None else "complete"
        tqdm.write(
            f"{label}: {status}, steps = {len(train_loss)}, elapsed = {elapsed:.1f} s, selected loss = {best_loss:.6g}"
        )
        return {
            "best_state": best_state,
            "best_loss": best_loss,
            "best_step": best_step,
            "stopped_step": stopped_step,
            "train_loss": train_loss,
            "validation_loss": validation_loss,
            "validation_steps": validation_steps,
            "initial_validation_loss": initial_validation,
            "steps": len(train_loss),
            "elapsed_seconds": elapsed,
            "successful": successful,
        }

    # Fit the shared GP and report optional holdout predictive density
    def optimize_hyperparameters(
        self,
        num_steps: int = 100,
        lr: float = 0.03,
        optimize_params: list | None = None,
        X_val: torch.Tensor | None = None,
        Y_val: torch.Tensor | None = None,
        Y_val_err: torch.Tensor | None = None,
        patience: int | None = None,
        min_delta: float = 0.0,
        validation_interval: int = 1,
        log_interval: int = 20,
        restarts: int | None = None,
    ) -> dict:
        del optimize_params
        if num_steps < 1:
            raise ValueError("GP optimization steps must be positive")
        if lr <= 0.0:
            raise ValueError("GP learning rate must be positive")
        if restarts is not None and restarts < 1:
            raise ValueError("GP optimization restarts must be positive")
        if min_delta < 0.0:
            raise ValueError("GP minimum validation improvement must be non-negative")
        if validation_interval < 1:
            raise ValueError("GP validation interval must be at least 1")
        if log_interval < 1:
            raise ValueError("GP log interval must be at least 1")

        fit_restarts = int(self.settings["fit_restarts"] if restarts is None else restarts)
        initial_state = self._parameter_snapshot()
        started = time.perf_counter()
        best_state = None
        best_loss = math.inf
        best_restart = None
        best_result = None
        train_loss = []
        validation_loss = []
        validation_steps = []
        restart_stats = []
        has_validation = X_val is not None and Y_val is not None and len(X_val) > 0

        with tqdm(total=fit_restarts, desc="Exact GP restarts", unit="restart") as progress:
            for restart in range(fit_restarts):
                result = self._fit_validation_restart(
                    initial_state=initial_state,
                    restart=restart,
                    total_restarts=fit_restarts,
                    num_steps=int(num_steps),
                    learning_rate=float(lr),
                    X_val=X_val,
                    Y_val=Y_val,
                    Y_val_err=Y_val_err,
                    patience=patience,
                    min_delta=float(min_delta),
                    validation_interval=int(validation_interval),
                    log_interval=int(log_interval),
                )
                train_loss.extend(result["train_loss"])
                step_offset = restart * int(num_steps)
                validation_loss.extend(result["validation_loss"])
                validation_steps.extend(step_offset + step for step in result["validation_steps"])
                restart_stats.append(
                    {
                        "restart": restart,
                        "steps": result["steps"],
                        "stopped_step": result["stopped_step"],
                        "elapsed_seconds": result["elapsed_seconds"],
                        "selected_loss": result["best_loss"],
                        "successful": result["successful"],
                    }
                )
                if result["successful"] and result["best_loss"] < best_loss:
                    best_state = result["best_state"]
                    best_loss = result["best_loss"]
                    best_restart = restart
                    best_result = result
                progress.update(1)
                progress.set_postfix(loss=f"{result['best_loss']:.4g}")

        if best_state is None or best_result is None:
            raise RuntimeError("Exact GP MAP optimization did not produce a finite state")
        self._restore_parameters(best_state)
        self._refresh()
        with torch.no_grad():
            selected_nlp = float(self.negative_log_posterior().detach().cpu())
        return {
            "best_negative_log_posterior": selected_nlp,
            "steps": len(train_loss),
            "restarts": fit_restarts,
            "successful_restarts": sum(result["successful"] for result in restart_stats),
            "restart_stats": restart_stats,
            "train_loss": train_loss,
            "validation_loss": validation_loss,
            "validation_steps": validation_steps,
            "initial_validation_loss": best_result["initial_validation_loss"],
            "best_validation_loss": best_loss if has_validation else None,
            "best_restart": best_restart,
            "best_step": best_result["best_step"] if has_validation else None,
            "stopped_restart": best_restart if best_result["stopped_step"] is not None else None,
            "stopped_step": best_result["stopped_step"],
            "elapsed_seconds": time.perf_counter() - started,
            "device": str(self.device),
            "dtype": str(self.dtype),
        }

    # Compute a detached trainable parameter snapshot
    def _parameter_snapshot(self) -> dict[str, torch.Tensor]:
        return {name: parameter.detach().clone() for name, parameter in self.named_parameters()}

    # Restore only trainable parameters without replacing data buffers
    def _restore_parameters(self, snapshot: dict[str, torch.Tensor]) -> None:
        parameters = dict(self.named_parameters())
        with torch.no_grad():
            for name, value in snapshot.items():
                parameters[name].copy_(value)

    # Compute latent posterior marginal moments
    def posterior(self, values: torch.Tensor) -> GaussianPosterior:
        points = self._tensor(values)
        cross = self.kernel(points, self.X)
        mean = self.constant_mean + cross @ self.alpha
        solved = solve_triangular(self.cholesky, cross.transpose(-1, -2))
        prior = self.outputscale.expand(mean.shape)
        epistemic = torch.clamp(prior - solved.square().sum(dim=-2), min=self.settings["posterior_variance_floor"])
        return GaussianPosterior(mean, epistemic, self.noise.expand(mean.shape))

    # Compute a dense latent posterior covariance
    def posterior_covariance(self, values: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        points = self._tensor(values)
        cross = self.kernel(points, self.X)
        mean = self.constant_mean + cross @ self.alpha
        solved = solve_triangular(self.cholesky, cross.transpose(-1, -2))
        covariance = self.kernel(points, points) - solved.transpose(-1, -2) @ solved
        return mean, 0.5 * (covariance + covariance.transpose(-1, -2))

    # Compute posterior cross covariance between two point sets
    def posterior_cross_covariance(self, first: torch.Tensor, second: torch.Tensor) -> torch.Tensor:
        left = self._tensor(first)
        right = self._tensor(second)
        cross_left = self.kernel(left, self.X)
        cross_right = self.kernel(self.X, right)
        solved = torch.cholesky_solve(cross_right, self.cholesky)
        return self.kernel(left, right) - cross_left @ solved

    # Compute the posterior mean as one column
    def predict_mean(self, values: torch.Tensor) -> torch.Tensor:
        points = self._tensor(values)
        return (self.constant_mean + self.kernel(points, self.X) @ self.alpha).reshape(-1, 1)

    # Compute predictive marginal mean and variance as columns
    def predict_marginals(self, values: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        posterior = self.posterior(values)
        return posterior.mean.reshape(-1, 1), posterior.predictive_variance.reshape(-1, 1)

    # Evaluate predictive marginal mean and variance
    def forward(self, values: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        return self.predict_marginals(values)

    # Freeze hyperparameters while retaining input autograd
    def freeze_for_inference(self) -> None:
        for parameter in self.parameters():
            parameter.requires_grad_(False)
        self.X = self.X.detach()
        self.y = self.y.detach()
        self.known_variance = self.known_variance.detach()
        self.cholesky = self.cholesky.detach()
        self.alpha = self.alpha.detach()

    # Compute compact fitted model diagnostics
    def diagnostics(self) -> dict:
        return {
            "surrogate": f"exact_{self.settings['kernel']}_map_gp",
            "training_points": len(self.X),
            "lengthscale_min": float(torch.min(self.lengthscale).detach().cpu()),
            "lengthscale_median": float(torch.median(self.lengthscale).detach().cpu()),
            "lengthscale_max": float(torch.max(self.lengthscale).detach().cpu()),
            "outputscale": float(self.outputscale.detach().cpu()),
            "noise_variance": float(self.noise.detach().cpu()),
            "known_heteroscedastic": bool(torch.any(self.known_variance > 0.0)),
            "kernel_dimension": self.geometry.effective_dim,
            "cholesky_jitter": self.jitter,
            "settings": copy.deepcopy(self.settings),
        }


# Compute a stable inverse softplus for positive parameter initialization
def _inverse_softplus(value: torch.Tensor) -> torch.Tensor:
    return value + torch.log(-torch.expm1(-value))


# Compute a robust Cholesky factor with the smallest successful jitter
def robust_cholesky(
    matrix: torch.Tensor, *, attempts: int, initial_jitter: float, jitter_multiplier: float
) -> tuple[torch.Tensor, float]:
    symmetric = 0.5 * (matrix + matrix.transpose(-1, -2))
    identity = torch.eye(matrix.shape[-1], device=matrix.device, dtype=matrix.dtype)
    scale = torch.clamp(torch.mean(torch.diagonal(symmetric, dim1=-2, dim2=-1).abs()), min=1.0)
    jitter = float(initial_jitter)
    for _ in range(max(1, int(attempts))):
        try:
            return torch.linalg.cholesky(symmetric + jitter * scale * identity), jitter
        except torch.linalg.LinAlgError:
            pass
        jitter *= float(jitter_multiplier)
    minimum = torch.min(torch.linalg.eigvalsh(symmetric))
    shift = torch.clamp(-minimum + float(initial_jitter) * scale, min=jitter * scale)
    return torch.linalg.cholesky(symmetric + shift * identity), jitter

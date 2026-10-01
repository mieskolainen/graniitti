# Gaussian process input and target uncertainty closure
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import numpy as np
import torch
from core.inference.gp import GaussianProcess
from torch import Tensor

OUTPUT_DIR = pathlib.Path(__file__).resolve().parents[4] / "figs"


# Compute the smooth reference function used by the uncertainty study
def true_function(x: np.ndarray) -> np.ndarray:
    return np.sin(x)


# Render one fitted GP against the reference function and training sample
def plot_gp_results(
    X_train: Tensor,
    Y_train: Tensor,
    X_train_err: Tensor | None,
    Y_train_err: Tensor | None,
    gp_model: GaussianProcess,
    title: str,
    filename: pathlib.Path,
) -> None:
    X_test = np.linspace(0, 10, 200).reshape(-1, 1)
    mu_pred, var_pred = gp_model.predict(torch.tensor(X_test, dtype=torch.float64))
    mu_pred = mu_pred.detach().cpu().numpy()
    std_pred = np.sqrt(var_pred.detach().cpu().numpy())
    y_true = true_function(X_test.flatten())

    plt.figure(figsize=(8, 5))
    plt.plot(X_test, y_true, "k--", label="True function")
    if mu_pred.ndim == 1 or mu_pred.shape[1] == 1:
        mu_pred = mu_pred.flatten()
        std_pred = std_pred.flatten()
        plt.plot(X_test, mu_pred, "b", label="GP mean")
        plt.fill_between(
            X_test.flatten(), mu_pred - std_pred, mu_pred + std_pred, alpha=0.2, label="68% CL"
        )
    else:
        plt.plot(X_test, mu_pred[:, 0], "b", label="GP mean (dim0)")
        plt.fill_between(
            X_test.flatten(),
            mu_pred[:, 0] - std_pred[:, 0],
            mu_pred[:, 0] + std_pred[:, 0],
            alpha=0.2,
            label="68% CL (dim0)",
        )

    x_train = X_train.cpu().numpy().flatten()
    y_train = Y_train.cpu().numpy()
    if y_train.ndim > 1 and y_train.shape[1] == 1:
        y_train = y_train.flatten()
    xerr = X_train_err.cpu().numpy().flatten() if X_train_err is not None else None
    yerr = Y_train_err.cpu().numpy().flatten() if Y_train_err is not None else None

    plt.errorbar(x_train, y_train, xerr=xerr, yerr=yerr, fmt="ro", label="Training data")
    plt.title(title)
    plt.xlabel("x")
    plt.xlim(X_test.min(), X_test.max())
    plt.ylabel("y")
    plt.legend()
    plt.tight_layout()
    filename.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(filename)
    plt.close()
    print(f"Saved figure to {filename}")


# Run one GP uncertainty configuration
def run_gp_uncertainty(
    kernel: str, *, noise_std_x: float | None, noise_std_y: float | None
) -> None:
    np.random.seed(0)
    torch.manual_seed(0)
    sample_count = 30
    x_true = np.sort(np.random.uniform(0, 10, sample_count)).reshape(-1, 1)
    x_observed = (
        x_true
        if noise_std_x is None
        else x_true + np.random.normal(0, noise_std_x, size=x_true.shape)
    )
    y_observed = true_function(x_true.flatten())
    if noise_std_y is not None:
        y_observed += np.random.normal(0, noise_std_y, size=sample_count)
    x_train = torch.tensor(x_observed, dtype=torch.float64)
    y_train = torch.tensor(y_observed.reshape(-1, 1), dtype=torch.float64)
    x_error = (
        None
        if noise_std_x is None
        else torch.full((sample_count, 1), noise_std_x, dtype=torch.float64)
    )
    y_error = (
        None
        if noise_std_y is None
        else torch.full((sample_count, 1), noise_std_y, dtype=torch.float64)
    )
    mode = ("x" if noise_std_x is not None else "") + ("y" if noise_std_y is not None else "")
    descriptions = {
        "": "no measurement uncertainty",
        "x": "measurement uncertainty on X",
        "y": "measurement uncertainty on Y",
        "xy": "measurement uncertainty on X and Y",
    }
    gp = GaussianProcess(
        x_train,
        y_train,
        X_err=x_error,
        Y_err=y_error,
        noise=1e-6,
        device="cpu",
        kernel=kernel,
    )
    gp.optimize_hyperparameters(num_steps=100, lr=0.01)
    filename = OUTPUT_DIR / f"gp_{mode + '_' if mode else 'no_'}uncertainty_{kernel}.png"
    plot_gp_results(
        x_train,
        y_train,
        x_error,
        y_error,
        gp,
        f"GP [{kernel}] ({descriptions[mode]})",
        filename,
    )


# Run all standalone GP uncertainty configurations
def main() -> None:
    for kernel in ("RBF",):
        for noise_std_x, noise_std_y in (
            (None, None),
            (0.15, None),
            (None, 0.2),
            (0.15, 0.2),
        ):
            run_gp_uncertainty(
                kernel,
                noise_std_x=noise_std_x,
                noise_std_y=noise_std_y,
            )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# Plot Sommerfeld Rayleigh scalar diffraction benchmark output
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
"""
Plot Sommerfeld-Rayleigh scalar diffraction benchmark output

mikael.mieskolainen@cern.ch, 2026
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parents[3]
DEFAULT_INPUT = REPO_DIR / "3D.ascii"
DEFAULT_OUTPUT = SCRIPT_DIR / "figs" / "sommerfeld.pdf"
Y_LIMITS = (-1.0, 1.0)
Z_LIMITS = (-1.0, 5.0)
SCREEN_SEGMENTS = ((-1.0, -0.8), (-0.6, -0.2), (0.2, 0.6), (0.8, 1.0))


# Parse command line arguments for the diffraction analysis
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=DEFAULT_INPUT, help="sommerfeld 3D.ascii output"
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output PDF path")
    parser.add_argument("--x-index", type=int, default=0, help="Stored x-grid index")
    parser.add_argument("--time-index", type=int, default=0, help="Stored time-grid index")
    parser.add_argument("--show", action="store_true", help="Show the figure interactively")
    return parser.parse_args()


# Load the flattened complex tensor written by the sommerfeld executable
def load_tensor(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Sommerfeld array not found: {path}")

    data = np.loadtxt(path, delimiter=",", comments="#")
    data = np.atleast_2d(data)
    if data.shape[1] != 6:
        raise ValueError(f"Expected six columns in {path}, got {data.shape[1]}")
    if not np.all(np.isfinite(data)):
        raise ValueError(f"Non-finite values found in {path}")

    raw_indices = data[:, :4]
    indices = np.rint(raw_indices).astype(np.intp)
    if not np.allclose(raw_indices, indices):
        raise ValueError(f"Non-integer tensor indices found in {path}")
    if np.any(indices < 0):
        raise ValueError(f"Negative tensor indices found in {path}")

    shape = tuple(int(np.max(indices[:, axis])) + 1 for axis in range(4))
    flat_indices = np.ravel_multi_index(tuple(indices.T), shape)
    if np.unique(flat_indices).size != data.shape[0]:
        raise ValueError(f"Duplicate tensor indices found in {path}")
    if data.shape[0] != int(np.prod(shape)):
        raise ValueError(f"Incomplete tensor grid in {path}: shape {shape}, rows {data.shape[0]}")

    tensor = np.empty(shape, dtype=np.complex128)
    tensor[tuple(indices.T)] = data[:, 4] + 1j * data[:, 5]
    return tensor


# Compute a physical coordinate axis for one stored grid dimension
def physical_axis(limits: tuple[float, float], size: int) -> np.ndarray:
    return np.linspace(limits[0], limits[1], size)


# Extract one real field slice in z by y display order
def real_field_slice(tensor: np.ndarray, x_index: int, time_index: int) -> np.ndarray:
    if not 0 <= x_index < tensor.shape[0]:
        raise IndexError(f"x-index {x_index} outside [0, {tensor.shape[0] - 1}]")
    if not 0 <= time_index < tensor.shape[3]:
        raise IndexError(f"time-index {time_index} outside [0, {tensor.shape[3] - 1}]")
    return tensor[x_index, :, :, time_index].real.T


# Keep a square plot box across matplotlib versions
def square_axes(ax) -> None:
    try:
        ax.set_box_aspect(1)
    except AttributeError:
        ax.set_aspect("auto")


# Draw the opaque parts of the aperture plane
def draw_aperture_screen(ax) -> None:
    for lower, upper in SCREEN_SEGMENTS:
        ax.plot((lower, upper), (0.0, 0.0), "k-", linewidth=2.0)


# Save a figure into the requested output path
def save_figure(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved {path}")


# Plot one real Sommerfeld-Rayleigh field slice
def plot_diffraction(
    tensor: np.ndarray,
    x_index: int,
    time_index: int,
    output: Path,
    show: bool,
) -> None:
    import matplotlib.pyplot as plt

    field = real_field_slice(tensor, x_index, time_index)
    y_values = physical_axis(Y_LIMITS, tensor.shape[1])
    z_values = physical_axis(Z_LIMITS, tensor.shape[2])
    bound = float(np.max(np.abs(tensor)))
    if not bound > 0.0:
        raise ValueError("Sommerfeld tensor has no nonzero field values")

    fig, ax = plt.subplots()
    ax.imshow(
        field,
        extent=(y_values[0], y_values[-1], z_values[0], z_values[-1]),
        origin="lower",
        aspect="auto",
        cmap="hot",
        vmin=-bound,
        vmax=bound,
        interpolation="nearest",
    )
    draw_aperture_screen(ax)
    ax.set_xlabel(r"$y$")
    ax.set_ylabel(r"$z$")
    ax.set_xticks(np.arange(Y_LIMITS[0], Y_LIMITS[1] + 0.1, 0.2))
    ax.set_yticks(np.arange(Z_LIMITS[0], Z_LIMITS[1] + 0.1, 0.5))
    square_axes(ax)
    save_figure(fig, output)

    if show:
        plt.show()
    else:
        plt.close(fig)


# Run the Sommerfeld-Rayleigh diffraction analysis
def main() -> int:
    args = parse_args()
    tensor = load_tensor(args.input)
    print(f"Loaded {args.input}: tensor shape {tensor.shape}")
    plot_diffraction(tensor, args.x_index, args.time_index, args.output, args.show)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

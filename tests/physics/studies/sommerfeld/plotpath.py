#!/usr/bin/env python3
# Plot the Feynman path-integral benchmark distribution
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = SCRIPT_DIR / "pathintegral_reference.csv"
DEFAULT_OUTPUT = SCRIPT_DIR / "figs" / "pathintegral.pdf"
BIN_LOW_COLUMN = 0
BIN_HIGH_COLUMN = 1
VALUE_COLUMN = 2


# Parse command line arguments for the path-integral analysis
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="CSV histogram or captured pathmark output",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output PDF path")
    parser.add_argument("--show", action="store_true", help="Show the figure interactively")
    return parser.parse_args()


# Parse numeric comma-separated histogram rows while ignoring pathmark log text
def load_histogram(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Path-integral histogram not found: {path}")

    rows: list[list[float]] = []
    width: int | None = None
    for line_number, line in enumerate(path.read_text(encoding="ascii").splitlines(), start=1):
        text = line.strip()
        if not text or text.startswith("#") or "," not in text:
            continue
        try:
            row = [float(value.strip()) for value in text.split(",")]
        except ValueError:
            continue
        if len(row) < 3:
            continue
        if width is None:
            width = len(row)
        if len(row) != width:
            raise ValueError(f"Inconsistent column count at {path}:{line_number}")
        rows.append(row)

    if not rows:
        raise ValueError(f"No numeric histogram rows found in {path}")
    data = np.asarray(rows, dtype=float)
    validate_histogram(data, path)
    return data


# Validate path-integral histogram edges and positive-definite values
def validate_histogram(data: np.ndarray, path: Path) -> None:
    if not np.all(np.isfinite(data)):
        raise ValueError(f"Non-finite histogram values found in {path}")
    if np.any(data[:, BIN_HIGH_COLUMN] <= data[:, BIN_LOW_COLUMN]):
        raise ValueError(f"Non-positive histogram bin width found in {path}")
    if np.any(data[:, VALUE_COLUMN] < 0.0):
        raise ValueError(f"Negative positive-definite histogram value found in {path}")
    if not float(np.max(data[:, VALUE_COLUMN])) > 0.0:
        raise ValueError(f"Histogram has no nonzero values: {path}")


# Compute bin centers and a unit-normalized positive-definite distribution
def normalized_distribution(data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    centers = 0.5 * (data[:, BIN_LOW_COLUMN] + data[:, BIN_HIGH_COLUMN])
    values = data[:, VALUE_COLUMN]
    return centers, values / np.max(values)


# Save a figure into the requested output path
def save_figure(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved {path}")


# Plot the normalized path-integral detector distribution
def plot_path_integral(data: np.ndarray, output: Path, show: bool) -> None:
    import matplotlib.pyplot as plt

    centers, values = normalized_distribution(data)
    fig, ax = plt.subplots()
    ax.plot(centers, values, "k-", linewidth=1.5)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$|A|^2 / \max |A|^2$")
    ax.set_xlim(float(data[0, BIN_LOW_COLUMN]), float(data[-1, BIN_HIGH_COLUMN]))
    ax.set_ylim(bottom=0.0)
    save_figure(fig, output)

    if show:
        plt.show()
    else:
        plt.close(fig)


# Run the Feynman path-integral analysis
def main() -> int:
    args = parse_args()
    data = load_histogram(args.input)
    print(f"Loaded {args.input}: {data.shape[0]} bins")
    plot_path_integral(data, args.output, args.show)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

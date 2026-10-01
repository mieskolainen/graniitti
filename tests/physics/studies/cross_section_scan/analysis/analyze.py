#!/usr/bin/env python3
# Analyze cross sections as a function of energy with screening on and off
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
"""
Cross sections as a function of energy with screening on and off

mikael.mieskolainen@cern.ch, 2026
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_ENERGY_COUNT = 10
BARN_TO_MICROBARN = 1e6


@dataclass(frozen=True)
class ScanColumns:
    sqrts: int = 0
    xs_el: int = 4
    xs_sd: int = 5
    xs_dd: int = 6
    xs_el_screened: int = 7
    xs_sd_screened: int = 8
    xs_dd_screened: int = 9


# Parse command line arguments for the analysis paths
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--false-scan", type=Path, default=SCRIPT_DIR / "scan_screening_false.csv")
    parser.add_argument("--true-scan", type=Path, default=SCRIPT_DIR / "scan_screening_true.csv")
    parser.add_argument("--fig-dir", type=Path, default=SCRIPT_DIR / "figs")
    parser.add_argument("--show", action="store_true", help="Show figures interactively")
    parser.add_argument("--no-ratio-only", action="store_true", help="Do not write xsecratios.pdf")
    return parser.parse_args()


# Print the energy list used by the scan steering script
def print_energy_values() -> None:
    sqrts_values = np.logspace(1.3, 6.3, DEFAULT_ENERGY_COUNT)
    values = ",".join(f"{value:0.6E}" for value in sqrts_values)
    print(f"E={values}\n")


# Load one tab separated xscan table
def load_scan(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Input scan file not found: {path}")
    data = np.loadtxt(path, delimiter="\t", skiprows=1)
    data = np.atleast_2d(data)
    if data.shape[1] < 7:
        raise ValueError(f"Expected at least 7 columns in {path}, got {data.shape[1]}")
    return data


# Load false and true screening scans and combine the process columns
def load_combined_scan(false_path: Path, true_path: Path) -> np.ndarray:
    false_scan = load_scan(false_path)
    true_scan = load_scan(true_path)
    if false_scan.shape[0] != true_scan.shape[0]:
        raise ValueError("False and true screening scans have different row counts")
    if false_scan.shape[1] != true_scan.shape[1]:
        raise ValueError("False and true screening scans have different column counts")
    if not np.allclose(false_scan[:, 0], true_scan[:, 0]):
        raise ValueError("False and true screening scans have different energy values")

    combined = np.column_stack((false_scan[:, :7], true_scan[:, 4:7]))
    combined[:, 1:] *= BARN_TO_MICROBARN
    return combined


# Save a figure with tight bounding box into the requested path
def save_figure(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")


# Plot cross sections and their screening ratios in one stacked figure
def plot_cross_sections(
    scan: np.ndarray, columns: ScanColumns, fig_dir: Path, ratio_only: bool
) -> None:
    import matplotlib.pyplot as plt

    fig, (ax_top, ax_ratio) = plt.subplots(
        2,
        1,
        figsize=(6.0, 6.2),
        sharex=True,
        gridspec_kw={"height_ratios": [3.8, 1.0], "hspace": 0.05},
    )

    plot_upper_panel(ax_top, scan, columns)
    plot_ratio_panel(ax_ratio, scan, columns)
    save_figure(fig, fig_dir / "xsec.pdf")

    if ratio_only:
        fig_ratio, ax = plt.subplots(figsize=(5.8, 3.6))
        plot_ratio_panel(ax, scan, columns)
        save_figure(fig_ratio, fig_dir / "xsecratios.pdf")


# Plot the unscreened and screened exclusive process cross sections
def plot_upper_panel(ax, scan: np.ndarray, columns: ScanColumns) -> None:
    x_values = scan[:, columns.sqrts]
    series = [
        (columns.xs_el, "k-", r"$\pi^+\pi^-_{\rm EL}$"),
        (columns.xs_sd, "k--", r"$\pi^+\pi^-_{\rm SD}$"),
        (columns.xs_dd, "k:", r"$\pi^+\pi^-_{\rm DD}$"),
        (columns.xs_el_screened, "r-", r"$\pi^+\pi^-_{\rm EL}$ ($S^2$)"),
        (columns.xs_sd_screened, "r--", r"$\pi^+\pi^-_{\rm SD}$ ($S^2$)"),
        (columns.xs_dd_screened, "r:", r"$\pi^+\pi^-_{\rm DD}$ ($S^2$)"),
    ]
    for column, style, label in series:
        ax.plot(x_values, scan[:, column], style, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(r"$\sigma$ ($\mu$b)")
    ax.set_xlim(10.0, float(np.max(x_values)))
    ax.set_ylim(1e-1, 1e2)
    ax.legend(loc="upper left", frameon=False)
    ax.tick_params(axis="x", which="both", labelbottom=False)


# Plot the screened to unscreened cross section ratios
def plot_ratio_panel(ax, scan: np.ndarray, columns: ScanColumns) -> None:
    x_values = scan[:, columns.sqrts]
    pairs = [
        (columns.xs_el, columns.xs_el_screened, "k-", r"EL"),
        (columns.xs_sd, columns.xs_sd_screened, "k--", r"SD"),
        (columns.xs_dd, columns.xs_dd_screened, "k:", r"DD"),
    ]
    for unscreened, screened, style, label in pairs:
        ax.plot(x_values, scan[:, screened] / scan[:, unscreened], style, label=label)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\sqrt{s}$ (GeV)")
    ax.set_ylabel(r"$\langle S^2 \rangle$")
    ax.set_xlim(10.0, float(np.max(x_values)))
    ax.set_ylim(0.0, 0.5)
    ax.set_yticks(np.arange(0.0, 0.5, 0.1))


# Print the average cross section ratios used as diagnostics
def print_ratio_summary(scan: np.ndarray, columns: ScanColumns) -> None:
    print("Before screening:")
    sd_to_el = float(np.mean(scan[:, columns.xs_sd] / scan[:, columns.xs_el]))
    dd_to_sd = float(np.mean(scan[:, columns.xs_dd] / scan[:, columns.xs_sd]))
    print(f"SD_to_EL_S  = {sd_to_el:0.6g}")
    print(f"DD_to_SD_S  = {dd_to_sd:0.6g}")

    print("After screening:")
    sd_to_el_screened = float(
        np.mean(scan[:, columns.xs_sd_screened] / scan[:, columns.xs_el_screened])
    )
    dd_to_sd_screened = float(
        np.mean(scan[:, columns.xs_dd_screened] / scan[:, columns.xs_sd_screened])
    )
    print(f"SD_to_EL_S2 = {sd_to_el_screened:0.6g}")
    print(f"DD_to_SD_S2 = {dd_to_sd_screened:0.6g}")


# Execute the scan cross section analysis workflow
def main() -> None:
    args = parse_args()
    print_energy_values()

    columns = ScanColumns()
    scan = load_combined_scan(args.false_scan, args.true_scan)
    plot_cross_sections(scan, columns, args.fig_dir, not args.no_ratio_only)
    print_ratio_summary(scan, columns)

    if args.show:
        import matplotlib.pyplot as plt

        plt.show()


if __name__ == "__main__":
    main()

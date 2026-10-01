#!/usr/bin/env python3
# Plot eikonal densities in impact parameter space
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
"""
Read and plot eikonal densities in impact-parameter space

mikael.mieskolainen@cern.ch, 2026
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from core.tune.drivers.graniitti.eikonal import EikonalMatrixFile, read_matrix_output

from tests.physics.studies.elastic.analysis.steering import MODELS, scan_outputs

X_IND = 0
RE_IND = 1
IM_IND = 2
AMP_RE_IND = 3
AMP_IM_IND = 4

DEFAULT_SQRTS_PP = [200.0, 7000.0, 13000.0, 60000.0]
DEFAULT_SQRTS_PPBAR = [546.0, 1960.0]
GEV2FM = 0.1973
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parents[4]
IMPACT_SPIRAL_XLIM = (0.0, 0.045)
IMPACT_SPIRAL_YLIM = (0.0, 0.70)
IMPACT_SPIRAL_STEP = 20
IMPACT_SPIRAL_LABEL_B = (0.05, 0.2, 0.5, 1.0, 1.5, 2.0)


# Compute default sqrt(s) values for the selected beam mode
def default_sqrts_for_beam(mode: str) -> list[float]:
    if mode == "pp":
        return list(DEFAULT_SQRTS_PP)
    if mode == "ppbar":
        return list(DEFAULT_SQRTS_PPBAR)
    return sorted(DEFAULT_SQRTS_PP + DEFAULT_SQRTS_PPBAR)


# Parse command line arguments for the impact-space eikonal analysis
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=tuple(MODELS), default="single")
    parser.add_argument("--scan-dir", type=Path, help="Directory containing the steered pp and ppbar scans")
    parser.add_argument("--fig-dir", type=Path, default=REPO_DIR / "figs" / "elastic")
    parser.add_argument("--sqrts", type=float, nargs="+")
    parser.add_argument("--beam", choices=("auto", "pp", "ppbar"), default="auto")
    parser.add_argument("--channel", default="0,0", help="Channel pair i,j for multichannel files")
    parser.add_argument("--show", action="store_true", help="Show figures interactively")
    args = parser.parse_args()
    if args.scan_dir is None:
        args.scan_dir = REPO_DIR / "tmp" / "elastic" / args.model
    if args.sqrts is None:
        args.sqrts = default_sqrts_for_beam(args.beam)
    return args


# Parse one selected physical final state channel pair
def parse_channel(text: str) -> tuple[int, int]:
    try:
        first, second = text.split(",", 1)
        return int(first), int(second)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("channel must have form i,j") from exc


# Compute the signed PDG beams expected for one energy
def expected_beams(sqrts: float, mode: str) -> set[tuple[int, int]]:
    if mode == "pp":
        return {(2212, 2212)}
    if mode == "ppbar":
        return {(2212, -2212), (-2212, 2212)}
    if any(np.isclose(sqrts, value) for value in DEFAULT_SQRTS_PPBAR):
        return {(2212, -2212), (-2212, 2212)}
    return {(2212, 2212)}


# Compute the canonical signed PDG beam pair for labels
def canonical_beam(sqrts: float, mode: str) -> tuple[int, int]:
    if (2212, -2212) in expected_beams(sqrts, mode):
        return (2212, -2212)
    return (2212, 2212)


# Compute a compact label for the signed beam pair
def beam_label(beam: tuple[int, int]) -> str:
    if beam == (2212, 2212):
        return r"$pp$"
    if beam in {(2212, -2212), (-2212, 2212)}:
        return r"$p\bar{p}$"
    return rf"${beam[0]} {beam[1]}$"


# Compute the explicitly selected eikonal model after checking availability
def select_nchannels(
    files: list[EikonalMatrixFile],
    requested: int,
    sqrts_values: list[float],
    beam_mode: str,
    channel: tuple[int, int],
) -> int:
    if not files:
        raise FileNotFoundError("No matrix eikonal outputs found")
    if min(channel) < 0 or max(channel) >= requested:
        raise ValueError(f"Channel {channel} is outside the N={requested} Good Walker basis")
    matches = [
        item
        for item in files
        if item.nchannels == requested
        and (item.beam1, item.beam2) in expected_beams(item.sqrts, beam_mode)
        and any(np.isclose(item.sqrts, sqrts) for sqrts in sqrts_values)
    ]
    if matches:
        return requested
    raise FileNotFoundError(
        f"No matrix eikonal outputs match N{requested} and the requested energies"
    )


# Load one selected impact-space eikonal channel over the requested energies
def load_channel_data(
    files: list[EikonalMatrixFile],
    sqrts_values: list[float],
    beam_mode: str,
    nchannels: int,
    channel: tuple[int, int],
) -> list[tuple[float, np.ndarray]]:
    datasets: list[tuple[float, np.ndarray]] = []
    for sqrts in sqrts_values:
        matches = [
            item
            for item in files
            if item.nchannels == nchannels
            and (item.beam1, item.beam2) in expected_beams(sqrts, beam_mode)
            and np.isclose(item.sqrts, sqrts)
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected exactly one steered output for N{nchannels} sqrt(s)={sqrts:g}, "
                f"found {len(matches)}. Run the elastic study with the requested model first"
            )
        chosen = matches[0]
        output = read_matrix_output(chosen)
        opacity = output.eigen_opacity(channel)
        amplitude = output.physical_impact_amplitude(channel)
        data = np.column_stack(
            (output.b_node, opacity.real, opacity.imag, amplitude.real, amplitude.imag)
        )
        datasets.append((sqrts, data))
    if not datasets:
        raise FileNotFoundError(f"No matrix eikonal data loaded for N{nchannels} channel {channel}")
    return datasets


# Compute a compact legend label for one energy point
def legend_label(
    sqrts: float,
    beam: tuple[int, int],
) -> str:
    return rf"{beam_label(beam)} $\sqrt{{s}} = {sqrts:g}$ GeV"


# Compute a filename suffix for the selected channel
def suffix(nchannels: int, channel: tuple[int, int]) -> str:
    if nchannels == 1 and channel == (0, 0):
        return ""
    return f"_ch{channel[0]}{channel[1]}"


# Compute the model-specific output directory
def model_fig_dir(fig_dir: Path, nchannels: int) -> Path:
    return fig_dir / f"N{nchannels}"


# Keep square plot boxes across matplotlib versions
def square_axes(ax) -> None:
    try:
        ax.set_box_aspect(1)
    except AttributeError:
        ax.set_aspect("equal", adjustable="box")


# Compute finite points inside a rectangular plot window
def plot_window_mask(
    x: np.ndarray,
    y: np.ndarray,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> np.ndarray:
    return (
        np.isfinite(x)
        & np.isfinite(y)
        & (x >= xlim[0])
        & (x <= xlim[1])
        & (y >= ylim[0])
        & (y <= ylim[1])
    )


# Save a figure with tight bounding box into the requested path
def save_figure(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")


# Plot impact-space figures for one selected eikonal channel
def plot_channel(
    datasets: list[tuple[float, np.ndarray]],
    fig_dir: Path,
    beam_mode: str,
    nchannels: int,
    channel: tuple[int, int],
) -> None:
    import matplotlib.pyplot as plt

    labels = [legend_label(sqrts, canonical_beam(sqrts, beam_mode)) for sqrts, _ in datasets]
    out_suffix = suffix(nchannels, channel)
    step = 2

    fig, ax = plt.subplots()
    for _, data in datasets:
        ax.plot(GEV2FM * data[::step, X_IND], data[::step, IM_IND], "-")
    ax.set_prop_cycle(None)
    for _, data in datasets:
        ax.plot(GEV2FM * data[::step, X_IND], data[::step, RE_IND], "--")
    ax.set_title(r"Re[$\Omega$] (dashed), Im[$\Omega$] (solid)")
    ax.legend(labels, frameon=False)
    ax.set_xlabel(r"$b$ (fm)")
    ax.set_ylabel(r"$\Omega(s,b)$")
    ax.set_xlim(1e-2, 2.5)
    ax.set_ylim(bottom=0.0)
    square_axes(ax)
    save_figure(fig, fig_dir / f"b_space_omega{out_suffix}.pdf")

    fig, ax = plt.subplots()
    for _, data in datasets:
        ax.plot(GEV2FM * data[::step, X_IND], data[::step, RE_IND] / data[::step, IM_IND], "-")
    ax.legend(labels, frameon=False, loc="lower right")
    ax.set_xlabel(r"$b$ (fm)")
    ax.set_ylabel(r"Re $[\Omega(s,b)]$ / Im $[\Omega(s,b)]$")
    ax.axis([0.0, 4.0, 0.12, 0.24])
    square_axes(ax)
    save_figure(fig, fig_dir / f"b_space_omega_ratio{out_suffix}.pdf")

    fig, ax = plt.subplots()
    for _, data in datasets:
        amplitude = data[::step, AMP_RE_IND] + 1j * data[::step, AMP_IM_IND]
        ax.plot(GEV2FM * data[::step, X_IND], amplitude.imag, "-")
    ax.set_prop_cycle(None)
    for _, data in datasets:
        amplitude = data[:, AMP_RE_IND] + 1j * data[:, AMP_IM_IND]
        ax.plot(GEV2FM * data[:, X_IND], amplitude.real, "--")
    ax.set_title(r"Re[A] (dashed), Im[A] (solid)")
    ax.legend(labels, frameon=False)
    ax.set_xlabel(r"$b$ (fm)")
    ax.set_ylabel(r"$A_{el}(s,b)$")
    ax.set_xticks(np.linspace(0.0, 4.0, 11))
    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(bottom=0.0)
    square_axes(ax)
    save_figure(fig, fig_dir / f"b_space_amplitude{out_suffix}.pdf")

    fig, ax = plt.subplots()
    for index, (_, data) in enumerate(datasets):
        amplitude = data[:, AMP_RE_IND] + 1j * data[:, AMP_IM_IND]
        sample = np.arange(0, len(data), IMPACT_SPIRAL_STEP)
        mask = plot_window_mask(
            amplitude.real[sample],
            amplitude.imag[sample],
            IMPACT_SPIRAL_XLIM,
            IMPACT_SPIRAL_YLIM,
        )
        ax.plot(
            amplitude.real[sample][mask],
            amplitude.imag[sample][mask],
            ".",
            markersize=1.5,
            rasterized=True,
        )
        if index == len(datasets) - 1:
            for label_index, target_b in enumerate(IMPACT_SPIRAL_LABEL_B):
                k = int(np.argmin(np.abs(GEV2FM * data[:, X_IND] - target_b)))
                if plot_window_mask(
                    amplitude.real[k : k + 1],
                    amplitude.imag[k : k + 1],
                    IMPACT_SPIRAL_XLIM,
                    IMPACT_SPIRAL_YLIM,
                )[0]:
                    ax.annotate(
                        f"{GEV2FM * data[k, X_IND]:.2g}",
                        xy=(amplitude.real[k], amplitude.imag[k]),
                        xytext=(6, 5 if label_index % 2 == 0 else -8),
                        textcoords="offset points",
                        fontsize=6,
                        bbox={
                            "boxstyle": "round,pad=0.12",
                            "fc": "white",
                            "ec": "none",
                            "alpha": 0.75,
                        },
                        arrowprops={"arrowstyle": "-", "lw": 0.35, "color": "0.35"},
                        clip_on=True,
                        zorder=5,
                    )
    ax.set_xlabel(r"Re [$A_{el}(s,b)$]")
    ax.set_ylabel(r"Im [$A_{el}(s,b)$]")
    ax.axis([*IMPACT_SPIRAL_XLIM, *IMPACT_SPIRAL_YLIM])
    square_axes(ax)
    save_figure(fig, fig_dir / f"b_space_amplitude_spiral{out_suffix}.pdf")

    fig, ax = plt.subplots()
    for _, data in datasets:
        amplitude = data[:, AMP_RE_IND] + 1j * data[:, AMP_IM_IND]
        ax.plot(GEV2FM * data[:, X_IND], amplitude.real / amplitude.imag, "-")
    ax.legend(labels, frameon=False, loc="upper left")
    ax.set_xlabel(r"$b$ (fm)")
    ax.set_ylabel(r"Re [$A_{el}(s,b)$] / Im [$A_{el}(s,b)$]")
    ax.axis([0.0, GEV2FM * max(data[-1, X_IND] for _, data in datasets), 1e-2, 1.0])
    square_axes(ax)
    save_figure(fig, fig_dir / f"b_amplitude_re_im{out_suffix}.pdf")


# Execute the impact-space eikonal analysis workflow
def main() -> None:
    args = parse_args()
    matrix_files = scan_outputs(args.scan_dir, args.model, args.beam)
    channel = parse_channel(args.channel)

    nchannels = select_nchannels(matrix_files, MODELS[args.model], args.sqrts, args.beam, channel)
    datasets = load_channel_data(matrix_files, args.sqrts, args.beam, nchannels, channel)
    plot_channel(datasets, model_fig_dir(args.fig_dir, nchannels), args.beam, nchannels, channel)

    if args.show:
        import matplotlib.pyplot as plt

        plt.show()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# Compare single, double and triple channel eikonal arrays
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
"""
Compare single, double and triple channel eikonal arrays

mikael.mieskolainen@cern.ch, 2026
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from core.tune.drivers.graniitti.eikonal import read_matrix_output

from tests.physics.studies.elastic.analysis.steering import scan_outputs

X_IND = 0
RE_IND = 1
IM_IND = 2
AMP_RE_IND = 3
AMP_IM_IND = 4

GEV2FM = 0.1973
GEV2MB = 0.389
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parents[4]
MODEL_TO_NCHANNELS = {"single": 1, "double": 2, "triple": 3}
NCHANNELS_TO_MODEL = {value: key for key, value in MODEL_TO_NCHANNELS.items()}
MODEL_TO_FOLDER = {"single": "N1", "double": "N2", "triple": "N3"}
MODEL_STYLES = {
    "single": {"color": "black"},
    "double": {"color": "tab:red"},
    "triple": {"color": "tab:blue"},
}
CHANNEL_COLORS = [
    "black",
    "tab:red",
    "tab:blue",
    "tab:green",
    "tab:purple",
    "tab:orange",
    "tab:brown",
    "tab:pink",
    "tab:gray",
]
TSPACE_XLIM = (1e-2, 4.0)


@dataclass(frozen=True)
class ChannelData:
    model: str
    nchannels: int
    channel: tuple[int, int]
    sqrts: float
    impact: np.ndarray
    momentum: np.ndarray


# Parse command line arguments for the multichannel comparison
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-dir", type=Path, default=REPO_DIR / "tmp" / "multichannel")
    parser.add_argument("--fig-dir", type=Path, default=REPO_DIR / "figs" / "multichannel")
    parser.add_argument("--sqrts", type=float, nargs="+", default=[13000.0])
    parser.add_argument("--beam1", type=int, default=2212)
    parser.add_argument("--beam2", type=int, default=2212)
    parser.add_argument(
        "--channels",
        nargs="+",
        default=["auto"],
        help="Good Walker pair labels i,j or auto for upper triangular pairs",
    )
    parser.add_argument(
        "--ordered-channels",
        action="store_true",
        help="Use every ordered i,j pair instead of i<=j pairs",
    )
    parser.add_argument(
        "--models",
        choices=tuple(MODEL_TO_NCHANNELS),
        nargs="+",
        default=list(MODEL_TO_NCHANNELS),
    )
    parser.add_argument("--show", action="store_true", help="Show figures interactively")
    return parser.parse_args()


# Parse an eigenchannel pair given as i,j
def parse_channel(text: str) -> tuple[int, int]:
    try:
        first, second = text.split(",", 1)
        return int(first), int(second)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("channel must have form i,j") from exc


# Compute physically unique channel pairs for one N-channel model
def physical_channel_pairs(nchannels: int, ordered: bool) -> list[tuple[int, int]]:
    if ordered:
        return [(i, j) for i in range(nchannels) for j in range(nchannels)]
    return [(i, j) for i in range(nchannels) for j in range(i, nchannels)]


# Compute requested channel pairs for one model
def requested_channels(texts: list[str], nchannels: int, ordered: bool) -> list[tuple[int, int]]:
    if len(texts) == 1 and texts[0] == "auto":
        return physical_channel_pairs(nchannels, ordered)
    channels = [parse_channel(text) for text in texts]
    bad = [
        channel
        for channel in channels
        if channel[0] < 0 or channel[1] < 0 or channel[0] >= nchannels or channel[1] >= nchannels
    ]
    if bad:
        model = NCHANNELS_TO_MODEL.get(nchannels, f"N{nchannels}")
        print(f"Skipping channels outside {model}: {bad}")
    if ordered:
        return [
            channel
            for channel in channels
            if channel[0] >= 0
            and channel[1] >= 0
            and channel[0] < nchannels
            and channel[1] < nchannels
        ]
    return [
        channel
        for channel in channels
        if channel[0] >= 0
        and channel[1] >= 0
        and channel[0] < nchannels
        and channel[1] < nchannels
        and channel[0] <= channel[1]
    ]


# Load impact and momentum arrays for every requested model/channel pair
def load_channel_data(
    scan_dir: Path,
    models: list[str],
    channel_texts: list[str],
    ordered_channels: bool,
    beam: tuple[int, int],
    sqrts: float,
) -> list[ChannelData]:
    loaded_outputs = {}
    datasets: list[ChannelData] = []
    for model in models:
        nchannels = MODEL_TO_NCHANNELS[model]
        if beam not in ((2212, 2212), (2212, -2212)):
            raise ValueError('Recorded scans require pp or ppbar beams')
        files = scan_outputs(scan_dir / model, model, 'pp' if beam[1] == 2212 else 'ppbar')
        matches = [item for item in files if np.isclose(item.sqrts, sqrts)]
        if len(matches) != 1:
            raise ValueError(f'Expected one recorded {model} output at {sqrts:g} GeV')
        selected = matches[0]
        if selected.path not in loaded_outputs:
            loaded_outputs[selected.path] = read_matrix_output(selected)
        output = loaded_outputs[selected.path]
        for channel in requested_channels(channel_texts, nchannels, ordered_channels):
            opacity = output.eigen_opacity(channel)
            impact_amplitude = output.physical_impact_amplitude(channel)
            momentum_amplitude = output.physical_momentum_amplitude(channel)
            datasets.append(
                ChannelData(
                    model=model,
                    nchannels=nchannels,
                    channel=channel,
                    sqrts=sqrts,
                    impact=np.column_stack(
                        (
                            output.b_node,
                            opacity.real,
                            opacity.imag,
                            impact_amplitude.real,
                            impact_amplitude.imag,
                        )
                    ),
                    momentum=np.column_stack(
                        (output.q2_node, momentum_amplitude.real, momentum_amplitude.imag)
                    ),
                )
            )
    if not datasets:
        raise FileNotFoundError("No multichannel eikonal arrays matched the requested selection")
    return datasets


# Save a figure with tight bounding box into the requested path
def save_figure(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    print(f"Saved {path}")


# Keep square plot boxes across matplotlib versions
def square_axes(ax) -> None:
    try:
        ax.set_box_aspect(1)
    except AttributeError:
        ax.set_aspect("equal", adjustable="box")


# Compute a compact legend label for one model
def model_label(dataset: ChannelData) -> str:
    return f"{dataset.model} (N={dataset.nchannels})"


# Compute a compact legend label for one channel pair
def channel_label(dataset: ChannelData) -> str:
    return f"({dataset.channel[0]},{dataset.channel[1]})"


# Compute a filename-safe label for one channel pair
def channel_slug(channel: tuple[int, int]) -> str:
    return f"{channel[0]}{channel[1]}"


# Compute plotting style for one model overlay
def model_style(dataset: ChannelData, _size: int, imaginary: bool) -> dict[str, object]:
    _ = imaginary
    style = MODEL_STYLES[dataset.model]
    return {
        "color": style["color"],
        "linewidth": 1.2,
    }


# Compute plotting style for one channel pair overlay
def channel_style(dataset: ChannelData, _size: int, imaginary: bool) -> dict[str, object]:
    _ = imaginary
    index = dataset.channel[0] * dataset.nchannels + dataset.channel[1]
    color = CHANNEL_COLORS[index % len(CHANNEL_COLORS)]
    return {
        "color": color,
        "linewidth": 1.2,
    }


# Group loaded arrays by channel pair
def group_by_channel(datasets: list[ChannelData]) -> dict[tuple[int, int], list[ChannelData]]:
    groups: dict[tuple[int, int], list[ChannelData]] = {}
    for dataset in datasets:
        groups.setdefault(dataset.channel, []).append(dataset)
    return {channel: groups[channel] for channel in sorted(groups)}


# Group loaded arrays by model name
def group_by_model(datasets: list[ChannelData]) -> dict[str, list[ChannelData]]:
    groups: dict[str, list[ChannelData]] = {}
    for dataset in datasets:
        groups.setdefault(dataset.model, []).append(dataset)
    return {model: groups[model] for model in MODEL_TO_NCHANNELS if model in groups}


# Compute the elastic dsigma/dt from the stored momentum-space amplitude
def dsigma_dt(momentum: np.ndarray, sqrts: float) -> np.ndarray:
    amplitude = momentum[:, RE_IND] + 1j * momentum[:, IM_IND]
    return np.abs(amplitude) ** 2 / (16.0 * np.pi * sqrts**4) * GEV2MB


# Compute y-axis limits for momentum-space amplitude plots
def momentum_amplitude_ylim(datasets: list[ChannelData]) -> tuple[float, float]:
    values: list[np.ndarray] = []
    for dataset in datasets:
        x_values = dataset.momentum[:, X_IND]
        amplitude = (
            dataset.momentum[:, RE_IND] + 1j * dataset.momentum[:, IM_IND]
        ) / dataset.sqrts**2
        mask = (
            np.isfinite(amplitude.real)
            & np.isfinite(amplitude.imag)
            & (x_values >= TSPACE_XLIM[0])
            & (x_values <= TSPACE_XLIM[1])
        )
        values.extend([amplitude.real[mask], amplitude.imag[mask]])

    non_empty = [value for value in values if value.size > 0]
    if not non_empty:
        return -20.0, 20.0
    finite = np.concatenate(non_empty)
    ymin = min(float(np.min(finite)), 0.0)
    ymax = max(float(np.max(finite)), 0.0)
    span = max(ymax - ymin, 1.0)
    return ymin - 0.08 * span, ymax + 0.08 * span


# Plot impact-space eikonal densities for the requested dataset group
def plot_impact_omega(
    datasets: list[ChannelData],
    fig_dir: Path,
    stem: str,
    style_func,
    label_func,
    title: str | None,
) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    for dataset in datasets:
        b_values = GEV2FM * dataset.impact[:, X_IND]
        ax.plot(
            b_values,
            dataset.impact[:, IM_IND],
            "-",
            label=f"{label_func(dataset)} Im",
            **style_func(dataset, len(b_values), imaginary=True),
        )
        ax.plot(
            b_values,
            dataset.impact[:, RE_IND],
            "--",
            label=f"{label_func(dataset)} Re",
            **style_func(dataset, len(b_values), imaginary=False),
        )
    if title is not None:
        ax.set_title(title)
    ax.set_xlabel(r"$b$ (fm)")
    ax.set_ylabel(r"$\Omega(s,b)$")
    ax.set_xlim(1e-2, 2.5)
    ax.set_ylim(bottom=0.0)
    ax.legend(frameon=False)
    square_axes(ax)
    save_figure(fig, fig_dir / f"{stem}_b_space_omega.pdf")
    plt.close(fig)


# Plot the eikonalized impact-space elastic amplitude
def plot_impact_amplitude(
    datasets: list[ChannelData],
    fig_dir: Path,
    stem: str,
    style_func,
    label_func,
    title: str | None,
) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    for dataset in datasets:
        b_values = GEV2FM * dataset.impact[:, X_IND]
        amplitude = dataset.impact[:, AMP_RE_IND] + 1j * dataset.impact[:, AMP_IM_IND]
        ax.plot(
            b_values,
            amplitude.imag,
            "-",
            label=f"{label_func(dataset)} Im",
            **style_func(dataset, len(b_values), imaginary=True),
        )
        ax.plot(
            b_values,
            amplitude.real,
            "--",
            label=f"{label_func(dataset)} Re",
            **style_func(dataset, len(b_values), imaginary=False),
        )
    if title is not None:
        ax.set_title(title)
    ax.set_xlabel(r"$b$ (fm)")
    ax.set_ylabel(r"$A_{el}(s,b)$")
    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(bottom=0.0)
    ax.legend(frameon=False)
    square_axes(ax)
    save_figure(fig, fig_dir / f"{stem}_b_space_amplitude.pdf")
    plt.close(fig)


# Plot momentum-space elastic differential cross sections
def plot_dsigma_dt(
    datasets: list[ChannelData],
    fig_dir: Path,
    stem: str,
    style_func,
    label_func,
    title: str | None,
) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    for dataset in datasets:
        ax.plot(
            dataset.momentum[:, X_IND],
            dsigma_dt(dataset.momentum, dataset.sqrts),
            "-",
            label=label_func(dataset),
            **style_func(dataset, len(dataset.momentum), imaginary=True),
        )
    ax.set_xscale("log")
    if title is not None:
        ax.set_title(title)
    ax.set_yscale("log")
    ax.set_xlabel(r"$-t$ (GeV$^2$)")
    ax.set_ylabel(r"$d\sigma/dt$ (mb/GeV$^2$)")
    ax.set_xlim(1e-2, 4.0)
    ax.set_ylim(1e-10, 1e6)
    ax.legend(frameon=False)
    square_axes(ax)
    save_figure(fig, fig_dir / f"{stem}_t_space_dsigmadt.pdf")
    plt.close(fig)


# Plot momentum-space elastic amplitudes normalized by s
def plot_momentum_amplitude(
    datasets: list[ChannelData],
    fig_dir: Path,
    stem: str,
    style_func,
    label_func,
    title: str | None,
) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    for dataset in datasets:
        amplitude = (
            dataset.momentum[:, RE_IND] + 1j * dataset.momentum[:, IM_IND]
        ) / dataset.sqrts**2
        ax.plot(
            dataset.momentum[:, X_IND],
            amplitude.imag,
            "-",
            label=f"{label_func(dataset)} Im",
            **style_func(dataset, len(dataset.momentum), imaginary=True),
        )
        ax.plot(
            dataset.momentum[:, X_IND],
            amplitude.real,
            "--",
            label=f"{label_func(dataset)} Re",
            **style_func(dataset, len(dataset.momentum), imaginary=False),
        )
    ax.axhline(0.0, color="0.82", linewidth=0.8, zorder=1)
    ax.set_xscale("log")
    if title is not None:
        ax.set_title(title)
    ax.set_xlabel(r"$-t$ (GeV$^2$)")
    ax.set_ylabel(r"$A_{el}(s,t) / s$")
    ax.set_xlim(*TSPACE_XLIM)
    ax.set_ylim(*momentum_amplitude_ylim(datasets))
    ax.legend(frameon=False)
    square_axes(ax)
    save_figure(fig, fig_dir / f"{stem}_t_space_amplitude.pdf")
    plt.close(fig)


# Plot the fundamental multichannel comparison figures for one dataset group
def plot_group(
    datasets: list[ChannelData],
    fig_dir: Path,
    stem: str,
    style_func,
    label_func,
    title: str | None,
) -> None:
    plot_impact_omega(datasets, fig_dir, stem, style_func, label_func, title)
    plot_impact_amplitude(datasets, fig_dir, stem, style_func, label_func, title)
    plot_dsigma_dt(datasets, fig_dir, stem, style_func, label_func, title)
    plot_momentum_amplitude(datasets, fig_dir, stem, style_func, label_func, title)


# Configure matplotlib for non-interactive PDF output when needed
def configure_matplotlib(show: bool) -> None:
    if show:
        return
    import matplotlib

    matplotlib.use("Agg")


# Plot model and channel comparisons at one energy
def plot_comparisons(datasets: list[ChannelData], fig_dir: Path) -> None:
    for channel, channel_datasets in group_by_channel(datasets).items():
        stem = f"model_compare_ch{channel_slug(channel)}"
        title = f"channel ({channel[0]},{channel[1]})"
        plot_group(
            channel_datasets, fig_dir / "compare", stem, model_style, model_label, title
        )

    for model, model_datasets in group_by_model(datasets).items():
        plot_group(
            model_datasets,
            fig_dir / MODEL_TO_FOLDER[model],
            "channel_compare",
            channel_style,
            channel_label,
            None,
        )


# Plot energy overlays for each model and channel pair
def plot_energies(datasets: list[ChannelData], fig_dir: Path) -> None:
    energies = list(dict.fromkeys(f"{dataset.sqrts:g}" for dataset in datasets))
    colors = {energy: f"C{index % 10}" for index, energy in enumerate(energies)}

    # Compute one color per energy with solid imaginary and dashed real amplitudes
    def style(dataset: ChannelData, _size: int, imaginary: bool) -> dict[str, object]:
        return {"color": colors[f"{dataset.sqrts:g}"], "linewidth": 1.2}

    # Compute the collision energy legend label
    def label(dataset: ChannelData) -> str:
        return rf"$\sqrt{{s}} = {dataset.sqrts:g}$ GeV"

    for model, model_datasets in group_by_model(datasets).items():
        for channel, channel_datasets in group_by_channel(model_datasets).items():
            plot_group(
                channel_datasets,
                fig_dir / MODEL_TO_FOLDER[model],
                f"energy_compare_ch{channel_slug(channel)}",
                style,
                label,
                f"{model}, channel ({channel[0]},{channel[1]})",
            )


# Execute the multichannel comparison workflow
def main() -> None:
    args = parse_args()
    configure_matplotlib(args.show)

    datasets = []
    for sqrts in args.sqrts:
        selected = load_channel_data(
            scan_dir=args.scan_dir,
            models=args.models,
            channel_texts=args.channels,
            ordered_channels=args.ordered_channels,
            beam=(args.beam1, args.beam2),
            sqrts=sqrts,
        )
        fig_dir = args.fig_dir / f"{sqrts:g}" if len(args.sqrts) > 1 else args.fig_dir
        plot_comparisons(selected, fig_dir)
        datasets.extend(selected)
    if len(args.sqrts) > 1:
        plot_energies(datasets, args.fig_dir)

    if args.show:
        import matplotlib.pyplot as plt

        plt.show()


if __name__ == "__main__":
    main()

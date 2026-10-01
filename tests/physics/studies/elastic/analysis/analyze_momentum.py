#!/usr/bin/env python3
# Plot eikonal amplitudes and elastic differential cross sections
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
"""
Plot eikonal amplitudes in momentum space & dsigma/dt
Plot elastic HEPData overlaid

mikael.mieskolainen@cern.ch, 2026
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from core.iceplot import comparison_metrics, write_validation_report
from core.io import readers, steering
from core.plot.plot import histhepdata
from core.stats.hist import hobj
from core.stats.validation import check_measurements
from core.tune.drivers.graniitti.eikonal import (
    EikonalMatrixFile,
    discover_matrix_outputs,
    read_matrix_output,
    select_matrix_output,
)
from scipy.integrate import trapezoid

from icepack._common.hepdata import read_table
from tests.physics.studies.elastic.analysis.steering import MODELS, scan_outputs

X_IND = 0
RE_IND = 1
IM_IND = 2

DEFAULT_SQRTS_PP = [200.0, 7000.0, 13000.0, 60000.0]
DEFAULT_SQRTS_PPBAR = [546.0, 1960.0]
DEFAULT_FACTOR_BY_SQRTS = {
    200.0: 4000.0,
    546.0: 400.0,
    1960.0: 20.0,
    7000.0: 1.0,
    13000.0: 1.0 / 20.0,
    60000.0: 1.0 / 1000.0,
}
DSIGMA_LEGEND_FONTSIZE = 6
PNG_DPI = 300
GEV2FM = 0.1973
GEV2MB = 0.389
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parents[4]
DEFAULT_ANALYSES = [
    REPO_DIR / "icepack/ELASTIC/STAR_1791591",
    REPO_DIR / "icepack/ELASTIC/CDF_359411",
    REPO_DIR / "icepack/ELASTIC/D0_1117021",
    REPO_DIR / "icepack/ELASTIC/TOTEM_1220862",
    REPO_DIR / "icepack/ELASTIC/TOTEM_922651",
    REPO_DIR / "icepack/ELASTIC/TOTEM_1710340",
]
MOMENTUM_SPIRAL_XLIM = (-5.0, 5.0)
MOMENTUM_SPIRAL_YLIM = (-10.0, 100.0)
MOMENTUM_AMPLITUDE_XLIM = (1e-2, 4.0)
QUASI_ENTROPY_T_SCALE = 1.0


@dataclass(frozen=True)
class HEPDataOverlay:
    path: Path
    energy: float
    label: str
    dataset: Path
    set_index: int


# Compute default sqrt(s) values for the selected beam mode
def default_sqrts_for_beam(mode: str) -> list[float]:
    if mode == "pp":
        return list(DEFAULT_SQRTS_PP)
    if mode == "ppbar":
        return list(DEFAULT_SQRTS_PPBAR)
    return sorted(DEFAULT_SQRTS_PP + DEFAULT_SQRTS_PPBAR)


# Compute default plotting scale factors for the requested energies
def default_factors_for_sqrts(sqrts_values: list[float]) -> list[float]:
    factors: list[float] = []
    for sqrts in sqrts_values:
        factor = 1.0
        for reference_sqrts, reference_factor in DEFAULT_FACTOR_BY_SQRTS.items():
            if np.isclose(sqrts, reference_sqrts):
                factor = reference_factor
                break
        factors.append(factor)
    return factors


# Parse command line arguments for the momentum-space eikonal analysis
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=tuple(MODELS), default="single")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--scan-dir", type=Path, help="Directory containing the steered pp and ppbar scans")
    source.add_argument("--eikonal-dir", type=Path, help="Read the newest matrix output for each model, beam and energy")
    parser.add_argument("--analysis", type=Path, nargs="+", default=DEFAULT_ANALYSES)
    parser.add_argument("--fig-dir", type=Path, default=REPO_DIR / "figs" / "elastic")
    parser.add_argument("--sqrts", type=float, nargs="+")
    parser.add_argument("--factor", type=float, nargs="+")
    parser.add_argument("--beam", choices=("auto", "pp", "ppbar"), default="auto")
    parser.add_argument("--channel", default="0,0", help="Channel pair i,j for multichannel files")
    parser.add_argument("--show", action="store_true", help="Show figures interactively")
    args = parser.parse_args()
    if args.scan_dir is None and args.eikonal_dir is None:
        args.scan_dir = REPO_DIR / "tmp" / "elastic" / args.model
    if args.sqrts is None:
        args.sqrts = default_sqrts_for_beam(args.beam)
    if args.factor is None:
        args.factor = default_factors_for_sqrts(args.sqrts)
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
        raise ValueError(f"Channel {channel} is outside the N={requested} physical basis")
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


# Load one selected momentum-space eikonal channel over the requested energies
def load_channel_data(
    files: list[EikonalMatrixFile],
    sqrts_values: list[float],
    factors: list[float],
    beam_mode: str,
    nchannels: int,
    channel: tuple[int, int],
) -> list[tuple[float, float, np.ndarray]]:
    datasets: list[tuple[float, float, np.ndarray]] = []
    for sqrts, factor in zip(sqrts_values, factors, strict=False):
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
        amplitude = output.physical_momentum_amplitude(channel)
        forward = output.physical_forward_amplitude(channel)
        q2 = np.concatenate(([0.0], output.q2_node))
        complete_amplitude = np.concatenate(([forward], amplitude))
        data = np.column_stack((q2, complete_amplitude.real, complete_amplitude.imag))
        datasets.append((sqrts, factor, data))
    if not datasets:
        raise FileNotFoundError(f"No matrix eikonal data loaded for N{nchannels} channel {channel}")
    return datasets


# Compute a compact legend label for one energy point
def legend_label(
    sqrts: float,
    beam: tuple[int, int],
) -> str:
    value, unit = energy_value_unit(sqrts)
    return rf"{beam_label(beam)} $\sqrt{{s}} = {value}$ {unit}"


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


# Save a figure with tight bounds and emit a matching PNG beside each PDF
def save_figure(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    if path.suffix.lower() == ".pdf":
        fig.savefig(path.with_suffix(".png"), dpi=PNG_DPI, bbox_inches="tight")


# Compute the scale factor only for a matching theory energy
def factor_for_energy(
    sqrts_values: list[float], factors: list[float], target: float
) -> float | None:
    if not sqrts_values:
        return None
    matches = np.flatnonzero(np.isclose(sqrts_values, target))
    return factors[int(matches[0])] if len(matches) else None


# Compute elastic differential cross section from momentum-space amplitude
def dxs_from_amplitude(sqrts: float, data: np.ndarray) -> np.ndarray:
    amplitude = data[:, RE_IND] + 1j * data[:, IM_IND]
    return np.abs(amplitude) ** 2 / (16.0 * np.pi * sqrts**4) * GEV2MB


# Compute quasi-entropy of the normalized differential t-spectrum
def quasi_entropy(dxs: np.ndarray, t_values: np.ndarray) -> float:
    sigma_el = trapezoid(dxs, t_values)
    if not np.isfinite(sigma_el) or sigma_el <= 0.0:
        return float("nan")
    probability_density = dxs / sigma_el
    density = np.zeros_like(probability_density)
    mask = np.isfinite(probability_density) & (probability_density > 0.0)
    density[mask] = -probability_density[mask] * np.log(
        probability_density[mask] * QUASI_ENTROPY_T_SCALE
    )
    return float(trapezoid(density, t_values))


# Compute y-axis limits for momentum-space amplitude plots
def momentum_amplitude_ylim(datasets: list[tuple[float, float, np.ndarray]]) -> tuple[float, float]:
    values: list[np.ndarray] = []
    for sqrts, _, data in datasets:
        amplitude = (data[:, RE_IND] + 1j * data[:, IM_IND]) / sqrts**2
        mask = (
            np.isfinite(amplitude.real)
            & np.isfinite(amplitude.imag)
            & (data[:, X_IND] >= MOMENTUM_AMPLITUDE_XLIM[0])
            & (data[:, X_IND] <= MOMENTUM_AMPLITUDE_XLIM[1])
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


# Compute the unique collision energy from the original HEPData qualifiers
def hepdata_energy(table: dict) -> float | None:
    values = [entry["value"] for name, entries in table["qualifiers"].items()
              if "SQRT(S)" in name.upper() for entry in entries]
    if not values:
        values = table.get("keywords", {}).get("cmenergies", [])
    energies = []
    for value in values:
        numbers = re.findall(r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?", str(value))
        if len(numbers) != 1:
            raise ValueError(f"Elastic table needs one collision energy: {value}")
        energies.append(float(numbers[0]) * (1000.0 if "TEV" in str(value).upper() else 1.0))
    if energies and not np.allclose(energies, energies[0]):
        raise ValueError("Elastic table contains multiple collision energies")
    return energies[0] if energies else None


# Compute the elastic beam label from the original HEPData reaction metadata
def hepdata_process(table: dict) -> str:
    reactions = list(table.get("keywords", {}).get("reactions", []))
    reactions.extend(entry["value"] for entry in table["qualifiers"].get("RE", []))
    for reaction in reactions:
        text = reaction.upper()
        if "PBAR P" in text:
            return r"$p\bar{p}$"
        if "P P --> P P" in text:
            return r"$pp$"
    return ""


# Compute compact experiment name from a HEPData record folder
def hepdata_experiment(path: Path) -> str:
    names = {
        "HEPData-ins1117021-v1-json": "D0",
        "HEPData-ins1220862-v1-json": "TOTEM low-$|t|$",
        "HEPData-ins1710340-v1-json": "TOTEM",
        "HEPData-ins1791591-v1-json": "STAR",
        "HEPData-ins201990-v1-json": "UA4",
        "HEPData-ins214689-v1-json": "ISR",
        "HEPData-ins359411-v1-json": "CDF",
        "HEPData-ins84176-v1-json": "ISR",
        "HEPData-ins922651-v1-json": "TOTEM high-$|t|$",
    }
    return names.get(path.parent.name, path.parent.name.removesuffix("-v1-json"))


# Compute a numeric legend value with a decimal marker for integer values
def energy_number(value: float) -> str:
    return f"{value:.3g}"


# Compute compact center-of-mass energy value and plain text unit
def energy_value_unit(sqrts: float) -> tuple[str, str]:
    if sqrts < 1000.0:
        return energy_number(sqrts), "GeV"
    return energy_number(sqrts / 1000.0), "TeV"


# Read and scale one elastic measurement through its configured icepack reader
def read_measurement(overlay):
    dataset, resolved = steering.load_dataset(str(overlay.dataset), cdir=REPO_DIR)
    entry = dataset["sets"][overlay.set_index]
    obs, _ = steering.load_observables(entry["obs"], dataset_path=resolved, cdir=REPO_DIR)
    data, obs = readers.read_hepdata(entry, dataset["datapath"], dataset["type"], obs,
                                   cdir=str(REPO_DIR), reader=dataset["reader"], dataset_path=resolved)
    record = histhepdata(data, obs, MC_XS_SCALE=1e3)["mandelstam_t"]
    return dataset, data["mandelstam_t"]["x"], record


# Compute all icepack selected HEPData overlays used in the d sigma / dt figure
def analysis_overlays(analysis_dirs: list[Path]) -> list[HEPDataOverlay]:
    overlays: list[HEPDataOverlay] = []
    for analysis_dir in analysis_dirs:
        dataset_file = analysis_dir.expanduser()
        if not dataset_file.is_absolute():
            dataset_file = REPO_DIR / dataset_file
        dataset_file = dataset_file.resolve() / "dataset.json"
        dataset, dataset_path = steering.load_dataset(str(dataset_file), cdir=str(REPO_DIR))
        for set_index, dataset_set in enumerate(dataset["sets"]):
            if not dataset_set.get("data", True):
                continue
            for histogram in dataset_set["hist"]:
                if histogram["obs"] != "mandelstam_t":
                    continue
                path = Path(
                    steering.resolve_data_reference(
                        histogram["file"],
                        datapath=dataset["datapath"],
                        dataset_path=dataset_path,
                        cdir=str(REPO_DIR),
                    )
                )
                table = read_table(path)
                if not table["values"]:
                    continue
                energy = hepdata_energy(table)
                if energy is None:
                    continue
                value, unit = energy_value_unit(energy)
                label = (
                    f"{hepdata_experiment(path)} {hepdata_process(table)} "
                    rf"$\sqrt{{s}} = {value}$ {unit}"
                )
                overlays.append(
                    HEPDataOverlay(
                        path=path,
                        energy=energy,
                        label=label,
                        dataset=dataset_file,
                        set_index=set_index,
                    )
                )
    return sorted(overlays, key=lambda item: (item.energy, item.path.parent.name, item.path.name))


# Overlay experimental elastic data and return legend handles with labels
def overlay_experiment(
    ax, analysis_dirs: list[Path], sqrts_values: list[float], factors: list[float], beam_mode: str
) -> tuple[list[object], list[str]]:
    from matplotlib import colormaps

    handles: list[object] = []
    labels: list[str] = []
    beams = {beam_label(canonical_beam(sqrts, beam_mode)) for sqrts in sqrts_values}
    overlays = analysis_overlays(analysis_dirs)
    colors = colormaps["plasma"](np.linspace(0.85, 0.05, len(overlays)))
    for index, overlay in enumerate(overlays):
        factor = factor_for_energy(sqrts_values, factors, overlay.energy)
        if factor is None:
            continue
        table = read_table(overlay.path)
        if hepdata_process(table) not in beams:
            continue
        _, x, record = read_measurement(overlay)
        reference = record["hdata"]
        valid = np.asarray(reference.valid, dtype=bool)
        x = np.asarray(x)[valid]
        y = factor * reference.counts_scaled[valid]
        errors = np.asarray(record.get("total_errs", reference.errs_scaled))
        yerr = factor * errors[..., valid]
        plotted_handle = ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt=".",
            markersize=2.4,
            color=colors[index],
        )
        handles.append(plotted_handle)
        labels.append(overlay.label)
    return handles, labels


# Integrate the native prediction over complete measured bins without extrapolation
def bin_average(x, values, bins):
    if np.any(np.diff(x) <= 0) or bins[0] < x[0] or bins[-1] > x[-1]:
        raise ValueError("Measured bins must lie inside the ordered amplitude grid")
    averages = []
    for lower, upper in zip(bins[:-1], bins[1:], strict=True):
        nodes = np.r_[lower, x[(x > lower) & (x < upper)], upper]
        averages.append(trapezoid(np.interp(nodes, x, values), nodes) / (upper - lower))
    return np.asarray(averages)


# Test the physical elastic amplitude against the same data and covariance as its icepack
def measurement_report(datasets, analysis_dirs, beam_mode):
    report = {"normalization": "cross_section", "sets": [], "excluded": [], "validation": {}}
    for overlay in analysis_overlays(analysis_dirs):
        matches = [(energy, values) for energy, _, values in datasets if np.isclose(energy, overlay.energy)
                   and hepdata_process(read_table(overlay.path)) == beam_label(canonical_beam(energy, beam_mode))]
        if not matches:
            report["excluded"].append({"measurement": overlay.label, "reason": "no matching beam and energy"})
            continue
        dataset, _, record = read_measurement(overlay)
        reference = record["hdata"]
        energy, values = matches[0]
        prediction = bin_average(values[:, X_IND], dxs_from_amplitude(energy, values), reference.bins)
        model = hobj(counts=prediction, errs=np.zeros_like(prediction), bins=reference.bins, valid=reference.valid)
        sample = {"label": "elastic amplitude", **comparison_metrics(model, reference, record["uncertainties"])}
        report["sets"].append({"name": overlay.label, "source": str(overlay.path),
                              "observables": [{"observable": "mandelstam_t", "kind": "histogram", "samples": [sample]}]})
        report["validation"] = {"measurement": dataset["validation"]["measurement"]}
    return report


# Compute B(t) on the native momentum grid, omitting the inserted near duplicate forward point
def local_slope(data, dxs):
    q2 = data[1:, X_IND]
    return 0.5 * (q2[:-1] + q2[1:]), -np.diff(np.log(dxs[1:])) / np.diff(q2)


# Plot momentum-space figures for one selected eikonal channel
def plot_channel(
    datasets: list[tuple[float, float, np.ndarray]],
    fig_dir: Path,
    analysis_dirs: list[Path],
    beam_mode: str,
    nchannels: int,
    channel: tuple[int, int],
) -> None:
    import matplotlib.pyplot as plt

    sqrts_values = [sqrts for sqrts, _, _ in datasets]
    factors = [factor for _, factor, _ in datasets]
    labels = [legend_label(sqrts, canonical_beam(sqrts, beam_mode)) for sqrts, _, _ in datasets]
    out_suffix = suffix(nchannels, channel)
    dxs_values = [dxs_from_amplitude(sqrts, data) for sqrts, _, data in datasets]

    fig, ax = plt.subplots()
    theory_handles = []
    tval = 3.52
    for (_sqrts, factor, data), dxs in zip(datasets, dxs_values, strict=False):
        ind = int(np.argmin(np.abs(data[:, X_IND] - tval)))
        if factor >= 1.0:
            txt = rf"$\times \, {factor:.0f}$"
        elif factor > 0.001:
            txt = rf"$\times \, {factor:.2f}$"
        else:
            txt = rf"$\times \, {factor:.3f}$"
        ax.text(tval, factor * dxs[ind] * 0.92, txt, fontsize=7)

    for (sqrts, factor, data), dxs in zip(datasets, dxs_values, strict=False):
        sigma_tot = data[0, IM_IND] / sqrts**2 * GEV2MB
        sigma_el = trapezoid(dxs, data[:, X_IND])
        q_entropy = quasi_entropy(dxs, data[:, X_IND])
        print(
            f"{sqrts:.1f}, tot: {sigma_tot:.1f} mb, el: {sigma_el:.1f} mb, "
            f"quasi-entropy: {q_entropy:.3f}"
        )
        (handle,) = ax.plot(data[:, X_IND], factor * dxs, "-", linewidth=1.1)
        theory_handles.append(handle)

    data_handles, data_labels = (overlay_experiment(ax, analysis_dirs, sqrts_values, factors, beam_mode)
                                 if channel == (0, 0) else ([], []))
    t_three_gluon = np.linspace(1.5, 8.0, 1000)

    scale = 5e-6  # Arbitrary scale
    (three_gluon_handle,) = ax.plot(t_three_gluon, scale * t_three_gluon ** (-8), "k-.")
    theory_legend = ax.legend(
        theory_handles + [three_gluon_handle],
        labels + [r"$t^{-8}$"],
        frameon=False,
        fontsize=DSIGMA_LEGEND_FONTSIZE,
        loc="upper right",
    )
    ax.add_artist(theory_legend)
    if data_handles:
        ax.legend(
            data_handles,
            data_labels,
            frameon=False,
            fontsize=DSIGMA_LEGEND_FONTSIZE,
            loc="upper center",
            bbox_to_anchor=(0.46, 0.98),
            borderaxespad=0.0,
            ncol=2 if len(data_handles) > 24 else 1,
            columnspacing=0.8,
            handletextpad=0.3,
        )
    ax.set_yscale("log")
    ax.set_xlabel(r"$-t$ (GeV$^2$)")
    ax.set_ylabel(r"$d\sigma/dt$ (mb/GeV$^2$)")
    ax.set_xticks(np.linspace(0.0, 8.0, 17))
    ax.axis([0.0, 3.5, 1e-10, 1e6])
    square_axes(ax)
    save_figure(fig, fig_dir / f"t_space_dsigmadt{out_suffix}.pdf")

    partial: list[np.ndarray] = []
    for (_, _, data), dxs in zip(datasets, dxs_values, strict=False):
        partial.append(local_slope(data, dxs)[1])

    for plot_index in (1, 2):
        fig, ax = plt.subplots()
        for (_, _, data), values in zip(datasets, partial, strict=False):
            ax.plot(0.5 * (data[1:-1, X_IND] + data[2:, X_IND]), values)
        if plot_index == 1:
            ax.plot(np.linspace(1e-3, 10.0, 10), np.zeros(10), "k-")
            ax.axis([0.0, 4.0, -20.0, 40.0])
            ax.set_xticks(np.linspace(0.0, 4.0, 11))
        else:
            ax.axis([0.0, 0.25, 5.0, 40.0])
        ax.legend(labels, frameon=False)
        ax.set_xlabel(r"$-t$ (GeV$^2$)")
        ax.set_ylabel(r"$B(t) \equiv \frac{d}{dt}\ln(d\sigma/dt)$ (GeV$^{-2}$)")
        square_axes(ax)
        save_figure(fig, fig_dir / f"t_space_bslope_{plot_index}{out_suffix}.pdf")

    fig, ax = plt.subplots()
    legs: list[str] = []
    first_t = 0.5 * (datasets[0][2][1:-1, X_IND] + datasets[0][2][2:, X_IND])
    for target_t in [0.0, 0.02, 0.05]:
        bin_index = int(np.argmin(np.abs(first_t - target_t)))
        b_values = np.asarray([values[bin_index] for values in partial])
        ax.plot(sqrts_values, b_values, "s-")
        legs.append(rf"$|t| = {first_t[bin_index]:.2f}$ GeV$^2$")
    ax.legend(legs, frameon=False, loc="lower right")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\sqrt{s}$ (GeV)")
    ax.set_ylabel(r"$B(t)$ (GeV$^{-2}$)")
    square_axes(ax)
    save_figure(fig, fig_dir / f"t_space_bt0{out_suffix}.pdf")

    fig, ax = plt.subplots()
    step = 2
    for sqrts, _, data in datasets:
        ax.plot(data[::step, X_IND], data[::step, IM_IND] / sqrts**2, "-")
    ax.set_prop_cycle(None)
    for sqrts, _, data in datasets:
        ax.plot(data[::step, X_IND], data[::step, RE_IND] / sqrts**2, "--")
    ax.set_title(r"Re[A] (dashed), Im[A] (solid)")
    ax.axhline(0.0, color="0.82", linewidth=0.8, zorder=1)
    ax.legend(labels, frameon=False)
    ax.set_xscale("log")
    ax.set_xlabel(r"$-t$ (GeV$^2$)")
    ax.set_ylabel(r"$A_{el}(s,t) / s$")
    ax.set_xlim(*MOMENTUM_AMPLITUDE_XLIM)
    ax.set_ylim(*momentum_amplitude_ylim(datasets))
    square_axes(ax)
    save_figure(fig, fig_dir / f"t_space_amplitude{out_suffix}.pdf")

    fig, ax = plt.subplots()
    for index, (sqrts, _, data) in enumerate(datasets):
        amplitude = (data[:, RE_IND] + 1j * data[:, IM_IND]) / sqrts**2
        mask = plot_window_mask(
            amplitude.real, amplitude.imag, MOMENTUM_SPIRAL_XLIM, MOMENTUM_SPIRAL_YLIM
        )
        ax.plot(amplitude.real[mask], amplitude.imag[mask], "-")
        if index == len(datasets) - 1:
            max_index = max(1, len(data) // 8)
            for k in range(0, max_index, 45):
                if mask[k]:
                    ax.text(
                        amplitude.real[k], amplitude.imag[k], f"{-data[k, X_IND]:.2f}", clip_on=True
                    )
    ax.set_xlabel(r"Re [$A_{el}(s,t) / s$]")
    ax.set_ylabel(r"Im [$A_{el}(s,t) / s$]")
    ax.set_xticks(np.linspace(MOMENTUM_SPIRAL_XLIM[0], MOMENTUM_SPIRAL_XLIM[1], 11))
    ax.axis([*MOMENTUM_SPIRAL_XLIM, *MOMENTUM_SPIRAL_YLIM])
    square_axes(ax)
    save_figure(fig, fig_dir / f"t_space_amplitude_spiral{out_suffix}.pdf")

    fig, ax = plt.subplots()
    for sqrts, _, data in datasets:
        amplitude = (data[:, RE_IND] + 1j * data[:, IM_IND]) / sqrts**2
        ax.plot(data[:, X_IND], amplitude.real / amplitude.imag, "-")
    ax.plot(np.linspace(1e-3, 10.0, 10), np.zeros(10), "k-")
    ax.legend(labels, frameon=False)
    ax.set_xlabel(r"$-t$ (GeV$^2$)")
    ax.set_ylabel(r"Re [$A_{el}(s,t)$] / Im [$A_{el}(s,t)$]")
    ax.axis([1e-2, 4.0, -15.0, 15.0])
    square_axes(ax)
    save_figure(fig, fig_dir / f"t_space_amplitude_re_im_ratio{out_suffix}.pdf")

    fig, ax = plt.subplots()
    for sqrts, _, data in datasets:
        amplitude = (data[:, RE_IND] + 1j * data[:, IM_IND]) / sqrts**2
        ax.plot(data[:, X_IND], amplitude.real / amplitude.imag, "-")
    ax.legend(labels, frameon=False)
    ax.set_xlabel(r"$-t$ (GeV$^2$)")
    ax.set_ylabel(r"Re [$A_{el}(s,t)$] / Im [$A_{el}(s,t)$]")
    ax.axis([0.0, 0.1, 0.0, 0.14])
    square_axes(ax)
    save_figure(fig, fig_dir / f"t_space_amplitude_re_im_ratio_zoom{out_suffix}.pdf")

    fig, ax = plt.subplots()
    legs = []
    first_t = 0.5 * (datasets[0][2][1:-1, X_IND] + datasets[0][2][2:, X_IND])
    for target_t in [0.0, 0.02, 0.05, 0.10]:
        bin_index = int(np.argmin(np.abs(first_t - target_t)))
        rho = np.asarray(
            [data[bin_index, RE_IND] / data[bin_index, IM_IND] for _, _, data in datasets]
        )
        ax.plot(sqrts_values, rho, "s-")
        legs.append(rf"$|t| = {first_t[bin_index]:.2f}$ GeV$^2$")
    ax.legend(legs, frameon=False)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\sqrt{s}$ (GeV)")
    ax.set_ylabel(r"$\rho(t) \equiv \mathrm{Re}[A_{el}(s,t)] / \mathrm{Im}[A_{el}(s,t)]$")
    energy_min = min(sqrts_values)
    energy_max = max(sqrts_values)
    if np.isclose(energy_min, energy_max):
        energy_min /= np.sqrt(2.0)
        energy_max *= np.sqrt(2.0)
    ax.axis([energy_min, energy_max, 0.0, 0.14])
    square_axes(ax)
    save_figure(fig, fig_dir / f"t_space_rho{out_suffix}.pdf")


# Execute the momentum-space eikonal analysis workflow
def main() -> None:
    args = parse_args()
    if len(args.factor) != len(args.sqrts):
        raise ValueError("--factor must have the same number of values as --sqrts")

    if args.eikonal_dir is None:
        matrix_files = scan_outputs(args.scan_dir, args.model, args.beam)
    else:
        available = discover_matrix_outputs(args.eikonal_dir)
        matrix_files = []
        for sqrts in args.sqrts:
            for beam in expected_beams(sqrts, args.beam):
                selected = select_matrix_output(available, MODELS[args.model], beam, sqrts)
                if selected is not None:
                    matrix_files.append(selected)
    channel = parse_channel(args.channel)

    nchannels = select_nchannels(matrix_files, MODELS[args.model], args.sqrts, args.beam, channel)
    datasets = load_channel_data(
        matrix_files, args.sqrts, args.factor, args.beam, nchannels, channel
    )
    plot_channel(
        datasets,
        model_fig_dir(args.fig_dir, nchannels),
        args.analysis,
        args.beam,
        nchannels,
        channel,
    )

    if channel == (0, 0):
        report = measurement_report(datasets, args.analysis, args.beam)
        report.update(model=args.model, channel=list(channel), input_files=[str(item.path) for item in matrix_files])
        write_validation_report(report, model_fig_dir(args.fig_dir, nchannels) / "measurements.json")
        check_measurements(report)

    if args.show:
        import matplotlib.pyplot as plt

        plt.show()


if __name__ == "__main__":
    main()

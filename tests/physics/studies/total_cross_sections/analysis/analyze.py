# Plot and fit total and diffractive cross sections as a function of energy
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from core.iceplot import comparison_metrics, write_validation_report
from core.io.serialize import load_json_file
from core.stats.hist import hobj
from core.stats.validation import check_measurements


@dataclass(frozen=True)
class ColumnIndex:
    sqrts: int
    total: int
    inelastic: int
    elastic: int
    sd: int
    dd: int
    sd_error: int = 10
    dd_error: int = 11


@dataclass(frozen=True)
class CouplingFit:
    values: np.ndarray
    chi2: np.ndarray
    best: float
    low_1sigma: float
    up_1sigma: float
    low_2sigma: float
    up_2sigma: float
    best_index: int
    low_index: int
    up_index: int


@dataclass(frozen=True)
class DiffractionBands:
    sd: np.ndarray
    sd_low: np.ndarray
    sd_up: np.ndarray
    dd: np.ndarray
    dd_low: np.ndarray
    dd_up: np.ndarray


@dataclass(frozen=True)
class Uncertainty:
    label: str
    minus: float
    plus: float


@dataclass(frozen=True)
class CrossSectionPoint:
    sqrts_gev: float
    value_mb: float
    uncertainties: tuple[Uncertainty, ...]


@dataclass(frozen=True)
class InclusiveCrossSectionPoint:
    channel: str
    sqrts_gev: float
    value_mb: float
    uncertainties: tuple[Uncertainty, ...]


@dataclass(frozen=True)
class InclusiveCrossSectionSeries:
    experiment: str
    collision: str
    label: str
    marker: str
    source: str
    points: tuple[InclusiveCrossSectionPoint, ...]


@dataclass(frozen=True)
class DiffractionSeries:
    experiment: str
    channel: str
    label: str
    marker: str
    selection: str
    source: str
    fit_observable: str | None
    points: tuple[CrossSectionPoint, ...]
    collision: str = "pp"


SD_XI_005 = "sd_xi_005"
DD_XI_005 = "dd_xi_005"


# Combine independent published uncertainty components in quadrature
def combined_uncertainty(
    point: CrossSectionPoint | InclusiveCrossSectionPoint,
) -> tuple[float, float]:
    minus = np.sqrt(sum(component.minus**2 for component in point.uncertainties))
    plus = np.sqrt(sum(component.plus**2 for component in point.uncertainties))
    return float(minus), float(plus)


# Parse command line arguments for the analysis paths and constants
def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Plot and fit total and diffractive cross sections"
    )
    parser.add_argument(
        "--scan", type=Path, default=script_dir / "scan.csv", help="Input scan.csv file"
    )
    parser.add_argument(
        "--fig-dir", type=Path, default=script_dir / "figs", help="Output PDF directory"
    )
    parser.add_argument(
        "--g3p-orig",
        type=float,
        default=None,
        help="Triple Pomeron coupling used in the simulation",
    )
    parser.add_argument(
        "--quiet-fit", action="store_true", help="Do not print every chi2 scan point"
    )
    parser.add_argument("--beam", choices=("pp", "ppbar"), help="Beam mode for a scan without metadata")
    return parser.parse_args()


# Print the energy list used by the scan steering scripts
def print_energy_values() -> None:
    sqrts_values = np.logspace(1.3, 6.3, 100)
    values = ",".join(f"{value:0.6E}" for value in sqrts_values)
    print(f"E={values}\n")


# Load the scan table and convert cross sections from barn to mb
def load_mc(scan_path: Path) -> np.ndarray:
    if not scan_path.exists():
        raise FileNotFoundError(f"Input scan file not found: {scan_path}")
    mc = np.loadtxt(scan_path, delimiter="\t", skiprows=1)
    mc = np.atleast_2d(mc)
    if mc.shape[1] != 12:
        raise ValueError(f"Expected 12 columns including four integration errors in {scan_path}, got {mc.shape[1]}")
    if not np.isfinite(mc).all() or np.any(mc[:, 0] <= 0) or np.any(np.diff(mc[:, 0]) <= 0) or np.any(mc[:, 8:] < 0):
        raise ValueError("Scan energies must increase and all cross sections and errors must be finite with nonnegative errors")
    mc[:, 1:] *= 1e3
    return mc


# Compute the generated scan column mapping
def scan_columns() -> ColumnIndex:
    return ColumnIndex(sqrts=0, total=1, inelastic=2, elastic=3, sd=6, dd=7)


# Load cited measurements from the reference JSON without duplicating numerical data
def measurement_series(channel):
    path = Path(__file__).resolve().parents[5] / "tests/references/total_cross_sections.json"
    groups = {}
    for row in load_json_file(path).values():
        metadata = dict(row["series"])
        if metadata["channel"].lower() != channel:
            continue
        if channel == "inclusive":
            metadata.pop("channel")
        measurement = row["measurement"]
        errors = tuple(Uncertainty(name, *pair) for name, pair in measurement.items()
                       if name != "value" and pair is not None)
        args = dict(sqrts_gev=row["sqrts"] * 1000, value_mb=measurement["value"], uncertainties=errors)
        point = InclusiveCrossSectionPoint(channel=row["channel"], **args) if channel == "inclusive" else CrossSectionPoint(**args)
        groups.setdefault(metadata["experiment"], (metadata, []))[1].append(point)
    series = InclusiveCrossSectionSeries if channel == "inclusive" else DiffractionSeries
    return tuple(series(**metadata, points=tuple(points)) for metadata, points in groups.values())


# Read the inclusive cross sections and their quoted error components
def total_measurements() -> tuple[InclusiveCrossSectionSeries, ...]:
    return measurement_series("inclusive")


# Read single diffraction cross sections with their original fiducial definitions
def sd_measurements() -> tuple[DiffractionSeries, ...]:
    return measurement_series("sd")


# Read double diffraction cross sections with their original fiducial definitions
def dd_measurements() -> tuple[DiffractionSeries, ...]:
    return measurement_series("dd")


# Compute only measurement series with an exact scan-observable match
def fitted_measurement_series(
    measurements: tuple[DiffractionSeries, ...],
) -> tuple[DiffractionSeries, ...]:
    return tuple(series for series in measurements if series.fit_observable is not None)


# Evaluate the scan prediction matched to one measurement series
def matched_prediction(
    series: DiffractionSeries,
    point: CrossSectionPoint,
    p3value: float,
    mc: np.ndarray,
    columns: ColumnIndex,
    g3p_orig: float,
) -> tuple[float, float]:
    if series.fit_observable == SD_XI_005:
        column, error_column, scale = columns.sd, columns.sd_error, p3value / g3p_orig
    elif series.fit_observable == DD_XI_005:
        column, error_column, scale = columns.dd, columns.dd_error, (p3value / g3p_orig)**2
    else:
        raise ValueError(f"No matched scan observable for {series.experiment} {series.channel}")
    value, error = scan_prediction(mc, point.sqrts_gev, column, error_column)
    return value * scale, error * abs(scale)


# Interpolate a positive cross section and propagate independent scan integration errors
def scan_prediction(mc, energy, column, error_column=None):
    energies = mc[:, 0]
    if np.any(np.diff(energies) <= 0) or not energies[0] <= energy <= energies[-1]:
        raise ValueError("Measurement energy must lie inside the ordered scan")
    if len(energies) == 1:
        indices, weights = [0], np.array([1.0])
    else:
        lower = int(np.clip(np.searchsorted(energies, energy) - 1, 0, len(energies) - 2))
        indices = [lower, lower + 1]
        fraction = np.log(energy / energies[lower]) / np.log(energies[lower + 1] / energies[lower])
        weights = np.array([1.0 - fraction, fraction])
    values = mc[indices, column]
    if np.any(values <= 0) or not np.isfinite(values).all():
        raise ValueError("Cross section interpolation requires positive finite values")
    prediction = float(np.exp(weights @ np.log(values)))
    error = 0.0 if error_column is None else prediction * np.linalg.norm(weights * mc[indices, error_column] / values)
    return prediction, float(error)


# Test the original model rates with matched beams and cuts, before fitting g3P
def measurement_report(mc, columns, inclusive, diffraction):
    criteria = load_json_file(Path(__file__).resolve().parents[5] / "icepack/_common/SETTINGS.json")
    report = {"validation": criteria, "normalization": "cross_section", "sets": [], "excluded": []}
    for series in (*inclusive, *diffraction):
        for point in series.points:
            channel = series.channel.lower() if isinstance(series, DiffractionSeries) else point.channel
            context = f"{series.experiment}/{channel}/{point.sqrts_gev:g}"
            if isinstance(series, DiffractionSeries) and series.fit_observable is None:
                reason = "different fiducial cuts: " + series.selection
            elif not mc[0, columns.sqrts] <= point.sqrts_gev <= mc[-1, columns.sqrts]:
                reason = "measurement energy outside the scan"
            else:
                reason = None
            if reason:
                report["excluded"].append({"measurement": context, "reason": reason, "source": series.source})
                continue
            error_column = getattr(columns, channel + "_error", None)
            value, error = scan_prediction(mc, point.sqrts_gev, getattr(columns, channel), error_column)
            down, up = combined_uncertainty(point)
            data_error = 0.5 * (down + up)
            mc_hist, data_hist = [hobj(counts=np.array([y]), errs=np.array([dy]), bins=np.array([0.0, 1.0]))
                                  for y, dy in ((value, error), (point.value_mb, data_error))]
            sample = {"label": "configured model", **comparison_metrics(mc_hist, data_hist)}
            report["sets"].append({"name": context, "source": series.source,
                                  "observables": [{"observable": "cross_section", "kind": "histogram", "samples": [sample]}]})
    if not report["sets"]:
        report["validation"] = {}
    return report


# Select the asymmetric measurement error on the prediction side
def residual_uncertainty(point: CrossSectionPoint, prediction: float) -> float:
    minus, plus = combined_uncertainty(point)
    return plus if prediction >= point.value_mb else minus


# Calculate the chi2 value for one trial Triple Pomeron coupling
def chi2_for_coupling(
    p3value: float,
    mc: np.ndarray,
    columns: ColumnIndex,
    measurements: tuple[DiffractionSeries, ...],
    g3p_orig: float,
) -> float:
    chi2 = 0.0
    for series in fitted_measurement_series(measurements):
        for point in series.points:
            prediction, mc_error = matched_prediction(series, point, p3value, mc, columns, g3p_orig)
            error = np.hypot(residual_uncertainty(point, prediction), mc_error)
            chi2 += ((prediction - point.value_mb) / error) ** 2
    return float(chi2)


# Find the first grid index away from the minimum above chi2_min + 1
def chi2_interval_index(
    values: np.ndarray, start: int, stop: int, step: int, threshold: float
) -> int:
    boundary = 0 if step < 0 else len(values) - 1
    for index in range(start, stop, step):
        if values[index] > threshold:
            return index
    return boundary


# Fit the Triple Pomeron coupling to matched measurements
def fit_g3p(
    mc: np.ndarray,
    columns: ColumnIndex,
    measurements: tuple[DiffractionSeries, ...],
    g3p_orig: float,
    quiet: bool,
) -> CouplingFit:
    p3values = np.linspace(0.05, 0.25, 1000)
    chi2_values = np.zeros_like(p3values)
    for index, value in enumerate(p3values):
        chi2 = chi2_for_coupling(value, mc, columns, measurements, g3p_orig)
        chi2_values[index] = chi2
        if not quiet:
            print(f"chi2 = {chi2:0.5f}, g3P = {value:0.5f} ")

    best_index = int(np.argmin(chi2_values))
    chi2_min = float(chi2_values[best_index])
    low_index = chi2_interval_index(chi2_values, best_index, -1, -1, chi2_min + 1.0)
    up_index = chi2_interval_index(chi2_values, best_index, len(chi2_values), 1, chi2_min + 1.0)
    low_2index = chi2_interval_index(chi2_values, best_index, -1, -1, chi2_min + 4.0)
    up_2index = chi2_interval_index(chi2_values, best_index, len(chi2_values), 1, chi2_min + 4.0)

    best = float(p3values[best_index])
    low_1sigma = float(p3values[low_index])
    up_1sigma = float(p3values[up_index])
    return CouplingFit(
        values=p3values,
        chi2=chi2_values,
        best=best,
        low_1sigma=low_1sigma,
        up_1sigma=up_1sigma,
        low_2sigma=float(p3values[low_2index]),
        up_2sigma=float(p3values[up_2index]),
        best_index=best_index,
        low_index=low_index,
        up_index=up_index,
    )


# Build SD and DD cross section bands from the fitted coupling interval
def diffraction_bands(
    mc: np.ndarray,
    columns: ColumnIndex,
    fit: CouplingFit,
    g3p_orig: float,
) -> DiffractionBands:
    sd = (mc[:, columns.sd] / g3p_orig) * fit.best
    sd_up = (mc[:, columns.sd] / g3p_orig) * fit.up_2sigma
    sd_low = (mc[:, columns.sd] / g3p_orig) * fit.low_2sigma
    dd = (mc[:, columns.dd] / g3p_orig**2) * fit.best**2
    dd_up = (mc[:, columns.dd] / g3p_orig**2) * fit.up_2sigma**2
    dd_low = (mc[:, columns.dd] / g3p_orig**2) * fit.low_2sigma**2
    return DiffractionBands(sd=sd, sd_low=sd_low, sd_up=sd_up, dd=dd, dd_low=dd_low, dd_up=dd_up)


# Save one figure as a tightly cropped PDF
def save_pdf(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.08)
    print(f"Saved {path}")


# Create a new Matplotlib figure and axes pair
def new_figure():
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    return fig, ax


# Draw a filled uncertainty band between two curves
def plot_fill(ax, x: np.ndarray, y_up: np.ndarray, y_low: np.ndarray, color: str = "0.7") -> None:
    lower = np.minimum(y_low, y_up)
    upper = np.maximum(y_low, y_up)
    ax.fill_between(x, lower, upper, facecolor=color, edgecolor="0.5", alpha=0.35, linewidth=0.4)


# Apply common axis and legend formatting used by the analysis plots
def finish_axes(ax, xlabel: str, ylabel: str, legend_loc: str | None = None) -> None:
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    if legend_loc is not None:
        ax.legend(loc=legend_loc, frameon=False)


# Plot one sourced measurement series with asymmetric published errors
def plot_diffraction_series(ax, series: DiffractionSeries) -> None:
    energy = np.asarray([point.sqrts_gev for point in series.points])
    values = np.asarray([point.value_mb for point in series.points])
    errors = np.asarray([combined_uncertainty(point) for point in series.points]).T
    ax.errorbar(
        energy,
        values,
        yerr=errors,
        fmt=series.marker,
        capsize=2.0,
        markersize=4.0,
        label=series.label,
    )


# Plot one sourced inclusive series with asymmetric published errors
def plot_inclusive_series(ax, series: InclusiveCrossSectionSeries):
    energy = np.asarray([point.sqrts_gev for point in series.points])
    values = np.asarray([point.value_mb for point in series.points])
    errors = np.asarray([combined_uncertainty(point) for point in series.points]).T
    return ax.errorbar(
        energy,
        values,
        yerr=errors,
        fmt=series.marker,
        capsize=2.0,
        markersize=4.0,
        label=series.label,
    )


# Plot total, inelastic, and elastic cross section measurements
def plot_total(
    mc: np.ndarray,
    columns: ColumnIndex,
    measurements: tuple[InclusiveCrossSectionSeries, ...],
    fig_dir: Path,
) -> None:
    import matplotlib.pyplot as plt

    fig, ax = new_figure()
    for column, style in [
        (columns.total, "k-"),
        (columns.inelastic, "r-"),
        (columns.elastic, "k--"),
    ]:
        ax.plot(mc[:, columns.sqrts], mc[:, column], style)

    handles = [plot_inclusive_series(ax, series) for series in measurements]
    ax.legend(
        handles,
        [series.label for series in measurements],
        loc="upper left",
        frameon=False,
    )
    ax.set_xscale("log")
    ax.set_xlim(10.0, float(np.max(mc[:, columns.sqrts])))
    ax.set_ylim(0.0, 110.0)
    finish_axes(ax, r"$\sqrt{s}$ (GeV)", r"$\sigma$ (mb)")
    save_pdf(fig, fig_dir / "total.pdf")
    plt.close(fig)


# Plot the Triple Pomeron coupling chi2 scan
def plot_chi2(fit: CouplingFit, fig_dir: Path) -> None:
    import matplotlib.pyplot as plt

    fig, ax = new_figure()
    y_values = np.linspace(0.5, float(np.max(fit.chi2)), 10)
    ax.plot(np.ones_like(y_values) * fit.low_1sigma, y_values, "r--", label=r"$\chi^2_{min} + 1$")
    ax.plot(np.ones_like(y_values) * fit.up_1sigma, y_values, "r--")
    ax.plot(fit.values, fit.chi2, "k-")
    ax.set_yscale("log")
    ax.set_xlim(float(np.min(fit.values)), float(np.max(fit.values)))
    ax.set_ylim(float(0.9 * np.min(fit.chi2)), float(np.max(fit.chi2)))
    finish_axes(ax, r"$g_{3P}$", r"$\chi^2$", "best")
    save_pdf(fig, fig_dir / "g3p_chi2fit.pdf")
    plt.close(fig)


# Plot diffraction ratios and the Pumplin bound
def plot_ratio(
    mc: np.ndarray,
    columns: ColumnIndex,
    bands: DiffractionBands,
    fig_dir: Path,
) -> None:
    import matplotlib.pyplot as plt

    fig, ax = new_figure()
    energy = mc[:, columns.sqrts]
    elastic = mc[:, columns.elastic]
    total = mc[:, columns.total]
    inelastic = mc[:, columns.inelastic]

    plot_fill(
        ax,
        energy,
        (elastic + bands.sd_up + bands.dd_up) / total,
        (elastic + bands.sd_low + bands.dd_low) / total,
    )
    ax.plot(
        energy,
        (elastic + bands.sd + bands.dd) / total,
        "k-",
        linewidth=1.1,
        label=r"$\sigma_{EL+SD+DD}/\sigma_{TOT}$",
    )
    plot_fill(ax, energy, bands.sd_up / inelastic, bands.sd_low / inelastic)
    ax.plot(energy, bands.sd / inelastic, "r--", label=r"$\sigma_{SD}/\sigma_{INEL}$")
    plot_fill(ax, energy, bands.dd_up / inelastic, bands.dd_low / inelastic)
    ax.plot(energy, bands.dd / inelastic, "k-.", label=r"$\sigma_{DD}/\sigma_{INEL}$")
    ax.plot(energy, elastic / inelastic, "k:", label=r"$\sigma_{EL}/\sigma_{INEL}$")
    ax.plot(energy, np.ones_like(energy) * 0.5, "r-", linewidth=1.1, label="MP bound")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(10.0, float(np.max(energy)))
    ax.set_ylim(1e-2, 1.0)
    finish_axes(ax, r"$\sqrt{s}$ (GeV)", "", "lower left")
    save_pdf(fig, fig_dir / "ratio.pdf")
    plt.close(fig)


# Plot single diffractive cross sections and measurements
def plot_sd(
    mc: np.ndarray,
    columns: ColumnIndex,
    bands: DiffractionBands,
    measurements: tuple[DiffractionSeries, ...],
    fig_dir: Path,
) -> None:
    import matplotlib.pyplot as plt

    fig, ax = new_figure()
    energy = mc[:, columns.sqrts]
    plot_fill(ax, energy, bands.sd_up, bands.sd_low)
    ax.plot(energy, bands.sd, "k-", label=r"GRANIITTI, $\xi_X < 0.05$")
    for series in measurements:
        plot_diffraction_series(ax, series)
    ax.set_xscale("log")
    ax.set_xlim(10.0, float(np.max(energy)))
    ax.set_ylim(0.0, 22.0)
    finish_axes(ax, r"$\sqrt{s}$ (GeV)", r"$\sigma_{SD}$ (mb)", "upper left")
    save_pdf(fig, fig_dir / "sd.pdf")
    plt.close(fig)


# Plot double diffractive cross sections and measurements
def plot_dd(
    mc: np.ndarray,
    columns: ColumnIndex,
    bands: DiffractionBands,
    measurements: tuple[DiffractionSeries, ...],
    fig_dir: Path,
) -> None:
    import matplotlib.pyplot as plt

    fig, ax = new_figure()
    energy = mc[:, columns.sqrts]
    plot_fill(ax, energy, bands.dd_up, bands.dd_low)
    ax.plot(energy, bands.dd, "k-", label=r"GRANIITTI, $\xi_{X,Y} < 0.05$")
    for series in measurements:
        plot_diffraction_series(ax, series)
    ax.set_xscale("log")
    ax.set_xlim(10.0, float(np.max(energy)))
    ax.set_ylim(0.0, 14.0)
    finish_axes(ax, r"$\sqrt{s}$ (GeV)", r"$\sigma_{DD}$ (mb)", "upper left")
    save_pdf(fig, fig_dir / "dd.pdf")
    plt.close(fig)


# Configure Matplotlib for non-interactive PDF output
def configure_matplotlib() -> None:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    plt.rcParams.update(
        {
            "font.size": 11,
            "legend.fontsize": 9,
            "axes.labelsize": 12,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


# Run the full cross section analysis workflow
def main() -> int:
    args = parse_args()
    configure_matplotlib()
    print_energy_values()

    columns = scan_columns()
    mc = load_mc(args.scan)
    metadata_path = args.scan.with_suffix('.json')
    metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
    coupling = args.g3p_orig if args.g3p_orig is not None else metadata.get('g3p')
    beam = args.beam or metadata.get('beam')
    if coupling is None or not np.isfinite(coupling) or coupling <= 0 or beam is None:
        raise ValueError('Scan metadata or explicit --g3p-orig and --beam are required')
    xsd = tuple(series for series in sd_measurements() if series.collision == beam)
    xdd = tuple(series for series in dd_measurements() if series.collision == beam)
    measurements = tuple(replace(series, points=tuple(point for point in series.points
        if mc[0, columns.sqrts] <= point.sqrts_gev <= mc[-1, columns.sqrts]))
        for series in fitted_measurement_series(xsd + xdd))
    measurements = tuple(series for series in measurements if series.points)
    if measurements:
        fit = fit_g3p(mc, columns, measurements, coupling, args.quiet_fit)
        print(f"g3P = {fit.best:.4f} ({fit.low_1sigma:.4f} ... {fit.up_1sigma:.4f})")
        bands = diffraction_bands(mc, columns, fit, coupling)
        plot_chi2(fit, args.fig_dir)
    else:
        print('No matching measurement energies, plotting the generated rates without a fit')
        sd, dd = mc[:, columns.sd], mc[:, columns.dd]
        bands = DiffractionBands(sd, sd, sd, dd, dd, dd)
    inclusive = tuple(series for series in total_measurements()
                      if series.collision.replace('pbar p', 'ppbar') == beam)
    plot_total(mc, columns, inclusive, args.fig_dir)
    plot_ratio(mc, columns, bands, args.fig_dir)
    plot_sd(mc, columns, bands, xsd, args.fig_dir)
    plot_dd(mc, columns, bands, xdd, args.fig_dir)
    report = measurement_report(mc, columns, inclusive, xsd + xdd)
    report.update(scan=str(args.scan.resolve()), model=metadata.get("model"), beam=beam, g3p=coupling)
    write_validation_report(report, args.fig_dir / "measurements.json")
    check_measurements(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

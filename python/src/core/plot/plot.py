# Histogram and ROC plotting functions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import math
from types import MappingProxyType

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import LogLocator
from sklearn.metrics import auc, roc_curve

from core.numerics import array
from core.stats import hist, uncertainty
from core.stats.uncertainty import source_covariance

DEFAULT_LEGEND_FONTSIZE = 5.5
MIN_LEGEND_FONTSIZE = 3.0
DEFAULT_FIGSIZE = (4, 3.75)


# Colors
imperial_dark_blue = (0, 0.24, 0.45)
imperial_light_blue = (0, 0.43, 0.69)
imperial_dark_red = (0.75, 0.10, 0.0)
imperial_green = (0.0, 0.54, 0.23)
imperial_dark_brown = (0.38, 0.05, 0.0)
imperial_light_brown = (1.0, 0.6, 0.2)
imperial_gray = (0.6, 0.6, 0.6)
imperial_light_gray = (0.3, 0.3, 0.3)


def colors(i, power=0.34):
    c = [
        imperial_dark_red,
        imperial_dark_blue,
        imperial_green,
        imperial_light_blue,
        imperial_dark_brown,
        imperial_light_brown,
        imperial_gray,
        imperial_light_gray,
    ]

    if i < len(c):
        return c[i]
    else:
        alpha = 1.0 / power
        color = c[i % len(c)]
        return tuple(np.clip(alpha * np.array(color), 0.0, 1.0))


""" Global marker styles

zorder : approximate plotting order
lw     : linewidth
ls     : linestyle
"""
errorbar_style = MappingProxyType({"zorder": 3, "ls": " ", "lw": 1, "marker": "o", "markersize": 2.5})
plot_style = MappingProxyType({"zorder": 2, "ls": "-", "lw": 1})
hist_style_step = MappingProxyType({"zorder": 0, "ls": "-", "lw": 1, "histtype": "step"})
hist_style_fill = MappingProxyType({"zorder": 0, "ls": "-", "lw": 1, "histtype": "stepfilled"})
hist_style_bar = MappingProxyType({"zorder": 0, "ls": "-", "lw": 1, "histtype": "bar"})
data_total_style = MappingProxyType({"alpha": 0.12, "hatch": "////", "lw": 0})
data_stat_style = MappingProxyType({"alpha": 0.14, "lw": 0})


def ratioerr(A, B, sigma_A, sigma_B, sigma_AB=0, EPS=1e-15):
    """Ratio f(A,B) = A/B error, by Taylor expansion of f."""
    A_safe = np.array(A, copy=True, dtype=float)
    B_safe = np.array(B, copy=True, dtype=float)
    A_safe[np.abs(A_safe) < EPS] = EPS
    B_safe[np.abs(B_safe) < EPS] = EPS
    return np.abs(A_safe / B_safe) * np.sqrt(
        (sigma_A / A_safe) ** 2 + (sigma_B / B_safe) ** 2 - 2 * sigma_AB / (A_safe * B_safe)
    )


# Compute ratio-panel errors according to one explicit display convention
def ratio_plot_errors(A, B, sigma_A, sigma_B, mode: str, reference: bool, EPS=1e-15):
    if mode not in {"combined", "separate", "numerator", "none"}:
        raise ValueError(f"Unknown ratio uncertainty mode '{mode}'")
    if mode == "none":
        return np.zeros_like(np.asarray(A, dtype=float))

    denominator = np.asarray(B, dtype=float).copy()
    denominator[np.abs(denominator) < EPS] = EPS
    if reference:
        return np.abs(np.asarray(sigma_B, dtype=float) / denominator)
    if mode == "combined":
        return ratioerr(A=A, B=B, sigma_A=sigma_A, sigma_B=sigma_B)
    return np.abs(np.asarray(sigma_A, dtype=float) / denominator)


DIMENSIONLESS_DIFFERENTIAL_UNITS = {"", "1", "unit", "rad"}
XS_UNITS = {"b": 1e-12, "mb": 1e-9, "ub": 1e-6, "nb": 1e-3, "pb": 1.0, "fb": 1e3}


# Identify cross section histograms and dimensionful cross section points
def cross_section_observable(observable):
    return observable.get("kind", "histogram") == "histogram" or observable["units"]["y"] in XS_UNITS


# Canonicalize unit strings before display and denominator logic
def normalize_unit_text(unit):
    return str(unit).strip()


# Compute true for units that should be displayed as dimensionless
def is_dimensionless_unit(unit):
    return (
        normalize_unit_text(unit).replace("$", "").replace("{", "").replace("}", "") in DIMENSIONLESS_DIFFERENTIAL_UNITS
    )


# Attach a unit suffix only for dimensional axes
def label_with_unit(label, unit):
    unit = normalize_unit_text(unit)
    if is_dimensionless_unit(unit):
        return label
    return f"{label} [{unit}]"


def parse_unit_power(unit):
    """Return base-unit powers from a simple unit string"""
    unit = normalize_unit_text(unit)
    if is_dimensionless_unit(unit):
        return {}

    compact = unit.replace("$", "").replace("{", "").replace("}", "").replace("**", "^").strip()
    if is_dimensionless_unit(compact):
        return {}

    if "^" in compact:
        base, exponent = compact.split("^", 1)
        return {base: int(exponent)}
    return {compact: 1}


def format_unit_powers(powers):
    """Format base-unit powers as a compact tex-compatible product"""
    pieces = []
    for base, exponent in powers.items():
        if exponent == 0:
            continue
        if exponent == 1:
            pieces.append(base)
        else:
            pieces.append(f"{base}$^{exponent}$")
    return r" ".join(pieces)


def differential_denominator_unit(units, explicit_dimensionless=False):
    """Return product unit for a differential denominator"""
    powers = {}
    for unit in units:
        for base, exponent in parse_unit_power(unit).items():
            powers[base] = powers.get(base, 0) + exponent
    formatted = format_unit_powers(powers)
    if formatted == "" and explicit_dimensionless:
        return "1"
    return formatted


# Build the y-axis unit without displaying dimensionless factors
def ylabel_unit(units):
    numerator = normalize_unit_text(units.get("y", ""))
    if is_dimensionless_unit(numerator):
        numerator = ""
    denominator = units.get("yden", None)

    if denominator is None:
        denominator = differential_denominator_unit([units.get("x", "")])
    elif is_dimensionless_unit(denominator):
        denominator = ""

    denominator = normalize_unit_text(denominator)
    if denominator == "":
        return numerator
    if numerator == "":
        return f"1 / {denominator}"
    return f"{numerator} / {denominator}"


# Compute uniformly spaced bin edges including the upper edge
def stepspace(start: float, stop: float, step: float):
    return np.arange(start, stop + step, step)


# Draw the unity guide behind ratio uncertainties and central values
def plot_horizontal_line(ax, color=(0.5, 0.5, 0.5), linewidth=0.9):
    xlim = ax.get_xlim()
    ax.plot(np.linspace(xlim[0], xlim[1], 2), np.array([1, 1]), color=color, linewidth=linewidth, zorder=-2)


def tick_calc(lim, step, N: int = 6):
    """Tick spacing calculator."""
    return [np.round(lim[0] + i * step, N) for i in range(1 + math.floor((lim[1] - lim[0]) / step))]


def set_axis_ticks(ax, ticks, dim="x"):
    """Set ticks of the axis."""
    if dim == "x":
        ax.set_xticks(ticks)
        ax.set_xticklabels(list(map(str, ticks)))
    elif dim == "y":
        ax.set_yticks(ticks)
        ax.set_yticklabels(list(map(str, ticks)))


def tick_creator(
    ax,
    xtick_step=None,
    ytick_step=None,
    ylim_ratio=(0.5, 1.5),
    ratio_plot=True,
    minorticks_on=True,
    ytick_ratio_step=0.25,
    labelsize=9,
    labelsize_ratio=8,
    **kwargs,
):
    """Axis tick constructor."""

    # Get limits
    xlim = ax[0].get_xlim()
    ylim = ax[0].get_ylim()

    # X-axis
    if xtick_step is not None:
        ticks = tick_calc(lim=xlim, step=xtick_step)
        set_axis_ticks(ax[-1], ticks, "x")

    # Y-axis
    if ytick_step is not None:
        ticks = tick_calc(lim=ylim, step=ytick_step)
        set_axis_ticks(ax[0], ticks, "y")

    # Y-ratio-axis
    if ratio_plot:
        ax[0].tick_params(labelbottom=False)
        ax[1].tick_params(axis="y", labelsize=labelsize_ratio)

        if ylim_ratio is not None:
            ticks = tick_calc(lim=ylim_ratio, step=ytick_ratio_step)
            ticks = ticks[0:-1]  # Remove the last (collapses with upper plot)
            set_axis_ticks(ax[1], ticks, "y")
            ax[1].set_ylim(ylim_ratio)

    # Tick settings
    for a in ax:
        if minorticks_on:
            a.minorticks_on()
        a.tick_params(top=True, bottom=True, right=True, left=True, which="both", direction="in", labelsize=labelsize)

    return ax


# Shrink an axis label if the rendered text exceeds the axis width
def shrink_axis_label_to_width(fig, ax, label, fontsize, min_fontsize=5.5):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    label_width = label.get_window_extent(renderer=renderer).width
    axis_width = ax.get_window_extent(renderer=renderer).width

    if label_width <= axis_width or label_width <= 0:
        return fontsize

    return max(min_fontsize, fontsize * axis_width / label_width)


def create_axes(
    xlabel="$x$",
    ylabel=r"Counts",
    ylabel_ratio="Ratio",
    xlim=(0, 1),
    ylim=None,
    ratio_plot=True,
    figsize=DEFAULT_FIGSIZE,
    fontsize=9,
    units=None,
    xlabel_suffix="",
    **kwargs,
):
    """Axes creator."""
    units = {"x": "", "y": ""} if units is None else units

    # Create subplots
    N = 2 if ratio_plot else 1
    gridspec_kw = {"height_ratios": (3.333, 1) if ratio_plot else (1,), "hspace": 0.0}
    fig, ax = plt.subplots(N, figsize=figsize, gridspec_kw=gridspec_kw)
    ax = [ax] if (N == 1) else ax

    # Axes limits
    xscale = kwargs.get("xscale", "linear")
    if xscale not in {"linear", "log"}:
        raise ValueError(f"Unknown xscale = {xscale}")
    for a in ax:
        a.set_xscale(xscale)
        a.set_xlim(*xlim)

    if ylim is not None:
        ax[0].set_ylim(*ylim)

    # Axes labels
    yden = units.get("yden", differential_denominator_unit([units.get("x", "")], explicit_dimensionless=True))
    if kwargs.get("density", False):
        ylabel = f"$1/N$  {ylabel}"
        if is_dimensionless_unit(yden):
            yden = ""
        if yden != "":
            ylabel = f"{ylabel} / [{yden}]"
    else:
        yunit = ylabel_unit(units)
        if yunit:
            ylabel = f"{ylabel}  [{yunit}]"
    xlabel = label_with_unit(xlabel, units["x"])
    if xlabel_suffix:
        xlabel = f"{xlabel} | {xlabel_suffix}"

    ax[0].set_ylabel(ylabel, fontsize=fontsize)
    ax[-1].set_xlabel(xlabel, fontsize=fontsize)
    ax[-1].xaxis.label.set_fontsize(
        shrink_axis_label_to_width(fig=fig, ax=ax[-1], label=ax[-1].xaxis.label, fontsize=fontsize)
    )

    # Ratio plot
    if ratio_plot:
        ax[1].set_ylabel(ylabel_ratio, fontsize=fontsize)

    # Setup ticks
    ax = tick_creator(ax=ax, ratio_plot=ratio_plot, **kwargs)

    return fig, ax


# Compute ordered legend handles and labels from one axes object
def ordered_handles_labels(ax=None, order=None, unique=False):
    def unique_everseen(seq, key=None):
        seen = set()
        seen_add = seen.add
        return [x for x, k in zip(seq, key, strict=False) if not (k in seen or seen_add(k))]

    if ax is None:
        ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()

    if len(handles) == 0:
        return [], []

    # Sort both labels and handles by labels
    labels, handles = zip(*sorted(zip(labels, handles, strict=False), key=lambda t: t[0]), strict=False)

    # Sort according to a given list, which may be incomplete
    if order is not None:
        keys = dict(zip(order, range(len(order)), strict=False))
        labels, handles = zip(
            *sorted(zip(labels, handles, strict=False), key=lambda t, keys=keys: keys.get(t[0], np.inf)), strict=False
        )

    # Keep only the first of each handle
    if unique:
        labels, handles = zip(*unique_everseen(zip(labels, handles, strict=False), key=labels), strict=False)
    return list(handles), list(labels)


# Draw an ordered legend on one axes object
def ordered_legend(ax=None, order=None, frameon=False, unique=False, **kwargs):
    handles, labels = ordered_handles_labels(ax=ax, order=order, unique=unique)
    if len(handles) == 0:
        return (handles, labels)
    ax.legend(handles, labels, frameon=frameon, **kwargs)

    return (handles, labels)


def generate_colormap():
    """Default colormap."""
    # Take colors
    color = plt.cm.Set1(np.linspace(0, 1, 10))

    # Add black
    black = np.ones((1, 4))
    black[:, 0:3] = 0.0
    color = np.concatenate((black, color))

    return color


# Compute fill properties without line or marker only arguments
def _filled_error_style(style, *, zorder, defaults):
    output = {**defaults, **style, "zorder": zorder}
    for key in ("histtype", "ls", "linestyle", "marker", "markersize"):
        output.pop(key, None)
    output["lw"] = 0
    return output


# Draw full-width boxes for every finite run of histogram uncertainties
def _draw_binned_error(ax, bins, y, err, color, style):
    y = np.asarray(y, dtype=float)
    err = np.asarray(err, dtype=float)
    if err.shape not in {y.shape, (2, len(y))} or len(bins) != len(y) + 1:
        raise ValueError("plot: binned uncertainty arrays must match the histogram bins")
    down, up = (err, err) if err.ndim == 1 else err

    output = []
    finite = np.isfinite(y) & np.isfinite(down) & np.isfinite(up)
    for start, stop in _finite_histogram_runs(np.where(finite, y, np.nan)):
        lower = y[start:stop] - np.abs(down[start:stop])
        upper = y[start:stop] + np.abs(up[start:stop])
        output.append(
            ax.fill_between(
                bins[start : stop + 1],
                np.append(lower, lower[-1]),
                np.append(upper, upper[-1]),
                step="post",
                color=color,
                **style,
            )
        )
    return output


# Draw a step histogram statistical uncertainty band
def hist_filled_error(ax, bins, cbins, y, err, color, **kwargs):
    del cbins
    style = _filled_error_style(kwargs, zorder=float(kwargs.get("zorder", 0.0)) - 1.0, defaults=data_stat_style)
    return _draw_binned_error(ax, bins, y, err, color, style)


# Test one rendered artist path against an inside legend box
def _path_intersects_legend(path, transform, legend_bbox, *, filled):
    try:
        display_path = transform.transform_path(path)
        return display_path.intersects_bbox(legend_bbox, filled=filled)
    except (AttributeError, TypeError, ValueError):
        return False


# Compute whether an inside legend covers rendered physics content
def _legend_overlaps_data(axis, legend, renderer):
    legend_bbox = legend.get_window_extent(renderer)
    for artist in (*axis.lines, *axis.patches):
        if not artist.get_visible():
            continue
        if _path_intersects_legend(artist.get_path(), artist.get_transform(), legend_bbox, filled=False):
            return True
    for collection in axis.collections:
        if not collection.get_visible():
            continue
        if any(
            _path_intersects_legend(path, collection.get_transform(), legend_bbox, filled=True)
            for path in collection.get_paths()
        ):
            return True
    return False


# Draw and shrink one inside legend until it clears plotted content
def _draw_fitted_inside_legend(fig, axis, handles, labels, properties):
    properties = properties.copy()
    fontsize = FontProperties(size=properties.get("fontsize", DEFAULT_LEGEND_FONTSIZE)).get_size_in_points()
    while True:
        properties["fontsize"] = fontsize
        legend = axis.legend(handles, labels, **properties)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        if not _legend_overlaps_data(axis, legend, renderer):
            return legend
        next_fontsize = max(MIN_LEGEND_FONTSIZE, 0.9 * fontsize)
        if math.isclose(next_fontsize, fontsize):
            return legend
        legend.remove()
        fontsize = next_fontsize


# Draw a legend either inside the plot axes or outside the figure frame
def draw_superplot_legend(fig, ax, legend_labels, legend_position="inside", legend_properties=None, legend_bbox=None):
    if legend_labels == []:
        return

    legend_properties = {} if legend_properties is None else legend_properties.copy()
    legend_properties.setdefault("fontsize", DEFAULT_LEGEND_FONTSIZE)

    legend_position = str(legend_position).lower()
    if legend_position == "none":
        return

    if legend_position == "inside":
        if legend_bbox is not None:
            legend_properties["bbox_to_anchor"] = tuple(legend_bbox)
        handles, labels = ordered_handles_labels(ax=ax[0], order=legend_labels)
        if handles:
            frameon = legend_properties.pop("frameon", False)
            _draw_fitted_inside_legend(fig, ax[0], handles, labels, {**legend_properties, "frameon": frameon})
        return

    handles, labels = ordered_handles_labels(ax=ax[0], order=legend_labels)
    if len(handles) == 0:
        return

    if legend_position == "outside-right":
        fig.subplots_adjust(right=0.68)
        loc = "center left"
        bbox = tuple(legend_bbox) if legend_bbox is not None else (0.70, 0.58)
    elif legend_position == "outside-top":
        fig.subplots_adjust(top=0.78)
        loc = "lower center"
        bbox = tuple(legend_bbox) if legend_bbox is not None else (0.50, 0.98)
    else:
        raise ValueError(f"Unknown legend_position = {legend_position}")

    fig.legend(
        handles,
        labels,
        loc=loc,
        bbox_to_anchor=bbox,
        frameon=legend_properties.pop("frameon", False),
        **legend_properties,
    )


# Replace excluded publication bins with gaps only in plotting arrays
def _mask_invalid_plot_bins(histogram, values, errors):
    values = np.asarray(values, dtype=float).copy()
    errors = np.asarray(errors, dtype=float).copy()
    valid = np.asarray(histogram.valid, dtype=np.bool_)
    if values.shape != errors.shape or valid.shape != values.shape:
        raise ValueError("plot: plotted values, errors and validity mask must be aligned")
    values[~valid] = np.nan
    errors[~valid] = np.nan
    return values, errors


# Apply the shared numerator and reference interval mask to ratio plotting arrays
def _mask_invalid_ratio_bins(numerator, reference, values, errors):
    values = np.asarray(values, dtype=float).copy()
    errors = np.asarray(errors, dtype=float).copy()
    numerator_valid = np.asarray(numerator.valid, dtype=np.bool_)
    reference_valid = np.asarray(reference.valid, dtype=np.bool_)
    if values.shape != errors.shape or numerator_valid.shape != values.shape or reference_valid.shape != values.shape:
        raise ValueError("plot: ratio values, errors and validity masks must be aligned")
    valid = numerator_valid & reference_valid
    values[~valid] = np.nan
    errors[~valid] = np.nan
    return values, errors


# Compute half-open runs of finite histogram bins
def _finite_histogram_runs(values):
    finite = np.isfinite(np.asarray(values, dtype=float))
    changes = np.diff(np.pad(finite.astype(np.int8), (1, 1)))
    starts = np.flatnonzero(changes == 1)
    stops = np.flatnonzero(changes == -1)
    return list(zip(starts, stops, strict=True))


# Draw step lines as independent finite runs so one gap cannot hide later bins
def _draw_step_histogram(ax, hdata, values, color, label, style):
    step_style = dict(style)
    step_style.pop("histtype", None)
    for index, (start, stop) in enumerate(_finite_histogram_runs(values)):
        ax.stairs(
            values[start:stop],
            hdata.bins[start : stop + 1],
            baseline=None,
            color=color,
            label=label if index == 0 else None,
            **step_style,
        )


# Draw one histogram record with its selected Matplotlib primitive
def _draw_histogram_record(ax, record, values, errors, color, label=None, data_uncertainty=True):
    hdata = record["hdata"]
    style = record["style"]
    hfunc = record["hfunc"]
    values, errors = _mask_invalid_plot_bins(hdata, values, errors)
    if data_uncertainty and "total_errs" in record:
        errors = np.where(hdata.valid, record["total_errs"], np.nan)
    stat_errors = record.get("stat_errs") if data_uncertainty else None
    if stat_errors is not None:
        stat_errors = np.where(hdata.valid, stat_errors, np.nan)
        total_style = _filled_error_style(
            style, zorder=float(style.get("zorder", 0.0)) - 1.0, defaults=data_total_style
        )
        total_style.update(data_total_style)
        _draw_binned_error(ax, hdata.bins, values, errors, color, total_style)
    if hfunc == "hist":
        if style.get("histtype", "bar") == "step":
            _draw_step_histogram(ax, hdata, values, color, label, style)
        else:
            ax.hist(hdata.cbins, bins=hdata.bins, weights=values, color=color, label=label, **style)
        if not record.get("_stack_component", False):
            if stat_errors is None:
                hist_filled_error(ax=ax, bins=hdata.bins, cbins=hdata.cbins, y=values, err=errors, color=color, **style)
            else:
                stat_style = _filled_error_style(
                    style, zorder=float(style.get("zorder", 0.0)) - 2.0, defaults=data_stat_style
                )
                stat_style.update(data_stat_style)
                _draw_binned_error(ax, hdata.bins, values, stat_errors, color, stat_style)
    elif hfunc == "errorbar":
        ax.errorbar(
            hdata.cbins, values, yerr=errors if stat_errors is None else stat_errors, color=color, label=label, **style
        )
    elif hfunc == "plot":
        ax.plot(hdata.cbins, values, color=color, label=label, **style)
        fill_style = {**style, "lw": 0, "zorder": float(style.get("zorder", 0.0)) - 1.0}
        down, up = (errors, errors) if np.ndim(errors) == 1 else errors
        ax.fill_between(hdata.cbins, values - down, values + up, alpha=0.2, color=color, **fill_style)


# Draw denominator uncertainty around unity as a ratio reference band
def _draw_ratio_reference_band(ax, hdata, errors, color):
    style = _filled_error_style({}, zorder=-1, defaults=data_total_style)
    return _draw_binned_error(
        ax, hdata.bins, np.ones_like(np.asarray(errors, dtype=float)), errors, color=color, style=style
    )


# Compute record copies with a plotting order above every reference record
def foreground_histogram_records(records: list[dict], references: list[dict]) -> list[dict]:
    reference_zorders = [float(record.get("style", {}).get("zorder", 0.0)) for record in references]
    minimum_zorder = max(reference_zorders, default=0.0) + 1.0

    output = []
    for record in records:
        visual = record.copy()
        style = record.get("style", {})
        zorder = max(float(style.get("zorder", 0.0)), minimum_zorder)
        visual["style"] = {**style, "zorder": zorder}
        output.append(visual)
    return output


# Compute validated automatic y-axis limits for one superimposed histogram plot
def _superplot_y_limits(
    *, bottom_count: float, ceiling_count: float, yscale: str, log_y_scale: tuple, requested
) -> tuple[float, float]:
    have_positive_bottom = np.isfinite(bottom_count) and bottom_count > 0.0
    have_positive_ceiling = np.isfinite(ceiling_count) and ceiling_count > 0.0
    if yscale == "log":
        lower_scale, upper_scale = map(float, log_y_scale)
        if lower_scale <= 0.0 or upper_scale <= 0.0:
            raise ValueError("Logarithmic y-limit scale factors must be positive")
        lower = bottom_count * lower_scale if have_positive_bottom else 0.1
        upper = ceiling_count * upper_scale if have_positive_ceiling else 1.0
        if requested is not None:
            requested_lower, requested_upper = map(float, requested)
            if np.isfinite(requested_lower) and requested_lower > 0.0:
                lower = max(lower, requested_lower)
            upper = requested_upper
    elif yscale == "linear":
        if requested is None:
            lower = 0.0
            upper = ceiling_count * 1.5 if have_positive_ceiling else 1.0
        else:
            lower, upper = map(float, requested)
    else:
        raise ValueError(f"Unknown yscale = {yscale}")

    if not np.isfinite(lower) or not np.isfinite(upper) or upper <= lower:
        raise ValueError(f"Invalid {yscale} y-axis limits [{lower}, {upper}]")
    return lower, upper


# Draw superimposed histograms and their optional ratio panel
def superplot(
    data: list,
    observable: dict = None,
    ratio_plot: bool = True,
    ratio_data: list | None = None,
    uncertainty_record: dict | None = None,
    legend_order: list[str] | None = None,
    yscale: str = "linear",
    ratio_uncertainty: str = "combined",
    legend_counts: bool = False,
    color: bool = None,
    legend_properties: dict = None,
    legend_position: str = "inside",
    legend_bbox: list = None,
    bottom_PRC: float = 0.0,
    log_y_scale: tuple = (0.1, 10),
    EPS: float = 1e-20,
    verbose: bool = False,
):
    if observable is None:
        observable = data[0]["obs"]

    if verbose:
        print(observable)

    fig, ax = create_axes(**observable, ratio_plot=ratio_plot)

    if color is None:
        color = generate_colormap()

    legend_labels = []

    # y-axis limit
    bottom_count = np.inf
    ceiling_count = 0.0

    # Plot histograms
    for i in range(len(data)):
        if data[i]["hdata"].is_empty:
            print(__name__ + f".superplot: Skipping empty histogram for entry {i}")
            continue

        c = data[i]["color"]
        if c is None:
            c = color[i]

        counts, errs = _mask_invalid_plot_bins(
            data[i]["hdata"], data[i]["hdata"].counts_scaled, data[i]["hdata"].errs_scaled
        )
        # ** For visualization autolimits **
        # Include positive bin contents and uncertainty edges in the logarithmic range
        down, up = np.where(data[i]["hdata"].valid, data[i]["total_errs"], np.nan) if "total_errs" in data[i] else (errs, errs)
        lower_values = np.concatenate((counts, counts - np.abs(down)))
        positive_counts = lower_values[np.isfinite(lower_values) & (lower_values > EPS)]
        if len(positive_counts):
            bottom_count = min(bottom_count, float(np.percentile(positive_counts, bottom_PRC)))
        upper_edges = counts + np.abs(up)
        finite_upper_edges = upper_edges[np.isfinite(upper_edges)]
        if len(finite_upper_edges):
            ceiling_count = max(ceiling_count, float(np.max(finite_upper_edges)))
        label = data[i]["label"]
        if legend_counts:
            label += f" $N={np.sum(data[i]['hdata'].counts):.1f}$"

        legend_labels.append(label)
        _draw_histogram_record(ax[0], data[i], counts, errs, c, label)

    if uncertainty_record is not None:
        total = uncertainty_record["hdata"]
        _draw_histogram_record(
            ax=ax[0], record=uncertainty_record, values=total.counts_scaled, errors=total.errs_scaled, color="black"
        )

    # Plot ratiohistograms
    if ratio_plot:
        plot_horizontal_line(ax[1])
        ratio_records = data if ratio_data is None else ratio_data

        for i in range(len(ratio_records)):
            if ratio_records[i]["hdata"].is_empty:
                print(__name__ + f".superplot: Skipping empty histogram for entry {i} (ratioplot)")
                continue

            c = ratio_records[i]["color"]
            if c is None:
                c = color[i]

            A = ratio_records[i]["hdata"].counts_scaled
            B = ratio_records[0]["hdata"].counts_scaled
            sigma_A = ratio_records[i]["hdata"].errs_scaled
            sigma_B = ratio_records[0]["hdata"].errs_scaled
            ratio_errs = ratio_plot_errors(
                A=A, B=B, sigma_A=sigma_A, sigma_B=sigma_B, mode=ratio_uncertainty, reference=i == 0
            )
            ratio = np.ones_like(B, dtype=float) if i == 0 else A / (B + 1e-30)
            ratio, ratio_errs = _mask_invalid_ratio_bins(
                ratio_records[i]["hdata"], ratio_records[0]["hdata"], ratio, ratio_errs
            )

            if i == 0:
                if ratio_uncertainty != "none":
                    _draw_ratio_reference_band(ax=ax[1], hdata=ratio_records[i]["hdata"], errors=ratio_errs, color=c)
                continue
            _draw_histogram_record(ax[1], ratio_records[i], ratio, ratio_errs, c, data_uncertainty=False)
    # Upper figure

    # Log y-scale
    ax[0].set_yscale(yscale)
    if yscale == "log":
        ax[0].yaxis.set_major_locator(LogLocator(base=10, numticks=100))
        ax[0].yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))

    # y-limits
    ax[0].set_ylim(
        _superplot_y_limits(
            bottom_count=bottom_count,
            ceiling_count=ceiling_count,
            yscale=yscale,
            log_y_scale=log_y_scale,
            requested=observable.get("ylim"),
        )
    )

    # Legend
    draw_superplot_legend(
        fig=fig,
        ax=ax,
        legend_labels=legend_labels if legend_order is None else legend_order,
        legend_position=legend_position,
        legend_properties=legend_properties,
        legend_bbox=legend_bbox,
    )
    return fig, ax


# Apply the configured unit density label or build the generic count label
def change2density_label(obs):
    density_ylabel = obs.get("density_ylabel")
    if density_ylabel is None:
        xlabel = obs["xlabel"].replace("$", "")
        density_ylabel = "$\\frac{1}{N} \\; " + f"dN/d{xlabel}$"
    obs["ylabel"] = density_ylabel
    obs["units"]["y"] = "1"

    return obs


def bins2label(bins):
    """Format bin interval for human-readable axis labels"""
    return f"[{bins[0]:0.2f}, {bins[1]:0.2f}]"


def category_interval_edges(bins, obs_name):
    """Validate and return the two edges of one category interval"""
    edges = np.asarray(bins, dtype=float)

    if len(edges) != 2:
        raise ValueError(f'plot: observable "{obs_name}" category intervals require exactly two bin edges')

    if not np.all(np.isfinite(edges)):
        raise ValueError(f'plot: observable "{obs_name}" category interval has non-finite edges')

    if edges[1] <= edges[0]:
        raise ValueError(f'plot: observable "{obs_name}" category interval edges are not increasing')

    return edges[0], edges[1]


# Compute the displayed bin area for differential normalization
def category_plane_area(x0_bins, x1_bins, symmetrized, obs_name):
    # Symmetrized filling sums both proton orders per displayed rectangle
    # A diagonal bin is filled once and retains its full rectangular area
    x0_min, x0_max = category_interval_edges(bins=x0_bins, obs_name=obs_name)
    x1_min, x1_max = category_interval_edges(bins=x1_bins, obs_name=obs_name)

    x0_width = x0_max - x0_min
    x1_width = x1_max - x1_min
    return x0_width * x1_width


def validate_category_axis_metadata(obs_config, obs_name):
    """Validate and return 2D category axis label metadata"""
    if "category_labels" not in obs_config or "category_units" not in obs_config:
        raise KeyError(
            f'plot: observable "{obs_name}" with dict bins requires "category_labels" and "category_units"'
        )

    labels = obs_config["category_labels"]
    units = obs_config["category_units"]

    if len(labels) != 2 or len(units) != 2:
        raise ValueError(
            f'plot: observable "{obs_name}" requires exactly two category_labels and two category_units'
        )

    return tuple(labels), tuple(units)


def category_differential_units(obs_config, obs_name):
    """Return observable units with the full category differential denominator"""
    _, category_units = validate_category_axis_metadata(obs_config=obs_config, obs_name=obs_name)
    units = copy.deepcopy(obs_config["units"])
    units["yden"] = differential_denominator_unit([units.get("x", "")] + list(category_units))
    return units


def category_xlabel_suffix(obs_config, obs_name, category_bins):
    """Format category axis labels from two sliced hyperbin dimensions"""
    labels, units = validate_category_axis_metadata(obs_config=obs_config, obs_name=obs_name)

    if len(category_bins) != 2:
        raise ValueError(f'plot: observable "{obs_name}" requires exactly two category bin arrays')

    return ", ".join(
        label_with_unit(f"{label} $\\in$ {bins2label(bins)}", unit)
        for label, unit, bins in zip(labels, units, category_bins, strict=False)
    )


# Validate and return the score metadata of one ROC observable
def roc_score_metadata(observable):
    roc = observable.get("roc", {})
    labels = list(roc.get("score_labels", []))
    directions = [str(value).lower() for value in roc.get("directions", [])]

    if len(labels) == 0:
        raise ValueError("plot: ROC observable requires at least one score label")
    if len(directions) != len(labels):
        raise ValueError("plot: ROC directions must match the number of score labels")
    if any(value not in ["higher", "lower"] for value in directions):
        raise ValueError("plot: ROC score directions must be 'higher' or 'lower'")

    default_styles = ["-", "--", "-.", ":"]
    linestyles = list(
        roc.get("linestyles", [default_styles[index % len(default_styles)] for index in range(len(labels))])
    )
    if len(linestyles) != len(labels):
        raise ValueError("plot: ROC linestyles must match the number of score labels")

    score_colors = list(roc.get("score_colors", [None] * len(labels)))
    if len(score_colors) != len(labels):
        raise ValueError("plot: ROC score colors must match the number of score labels")

    return labels, directions, linestyles, score_colors


# Convert one chunk of ROC scores and weights to validated arrays
def prepare_roc_data(scores, weights, score_count):
    scores = np.asarray(scores, dtype=float)
    weights = np.asarray(weights, dtype=float).reshape(-1)

    if scores.size == 0:
        scores = np.empty((0, score_count), dtype=float)
    elif scores.ndim == 1 and score_count == 1:
        scores = scores.reshape(-1, 1)
    elif scores.ndim == 1 and len(weights) == 1 and scores.shape[0] == score_count:
        scores = scores.reshape(1, score_count)

    if scores.ndim != 2 or scores.shape[1] != score_count:
        raise ValueError(f"plot: ROC projector must return {score_count} score value(s) per event")
    if scores.shape[0] != len(weights):
        raise ValueError("plot: ROC score and event-weight counts do not match")
    if not np.all(np.isfinite(scores)) or not np.all(np.isfinite(weights)):
        raise ValueError("plot: ROC scores and weights must be finite")
    if np.any(weights < 0.0):
        raise ValueError("plot: ROC curves require nonnegative event weights")

    return {"scores": scores, "weights": weights}


# Concatenate two independently processed ROC chunks
def fuse_roc_data(first, second):
    first_scores = np.asarray(first["scores"], dtype=float)
    second_scores = np.asarray(second["scores"], dtype=float)
    if first_scores.ndim != 2 or second_scores.ndim != 2:
        raise ValueError("plot: ROC chunk scores must be two-dimensional")
    if first_scores.shape[1] != second_scores.shape[1]:
        raise ValueError("plot: ROC chunks contain different score dimensions")

    return {
        "scores": np.concatenate((first_scores, second_scores), axis=0),
        "weights": np.concatenate((first["weights"], second["weights"]), axis=0),
    }


# Compute an exact weighted empirical ROC curve with tied scores grouped
def weighted_roc_curve(signal_scores, signal_weights, background_scores, background_weights, higher_is_signal=True):
    signal = prepare_roc_data(signal_scores, signal_weights, score_count=1)
    background = prepare_roc_data(background_scores, background_weights, score_count=1)
    if np.sum(signal["weights"]) <= 0.0 or np.sum(background["weights"]) <= 0.0:
        raise ValueError("plot: ROC signal and background require positive total weight")

    scores = np.concatenate((signal["scores"][:, 0], background["scores"][:, 0]))
    labels = np.concatenate(
        (np.ones(len(signal["weights"]), dtype=int), np.zeros(len(background["weights"]), dtype=int))
    )
    weights = np.concatenate((signal["weights"], background["weights"]))
    fpr, tpr, _ = roc_curve(
        labels, scores if higher_is_signal else -scores, sample_weight=weights, drop_intermediate=False
    )
    return fpr, tpr


# Integrate one ROC curve using the standard linear-FPR trapezoidal area
def weighted_roc_auc(fpr, tpr):
    return float(auc(fpr, tpr))


# Format one ROC legend entry using optional observable steering
def format_roc_curve_label(score_label, comparison_label, auc, roc_config):
    template = roc_config.get("label_template", None)
    if template is not None:
        try:
            return template.format(score=score_label, comparison=comparison_label, auc=auc)
        except (KeyError, ValueError) as exc:
            raise ValueError(f'plot: invalid ROC label_template "{template}"') from exc

    label = f"{score_label}: {comparison_label}"
    if roc_config.get("show_auc", True):
        label += f" (AUC = {auc:.3f})"
    return label


# Draw all configured score and sample comparisons on one ROC axes
def rocplot(
    data,
    observable=None,
    title="",
    title_loc="left",
    legend_properties=None,
    legend_position="inside",
    legend_bbox=None,
):
    if observable is None:
        observable = data[0]["obs"]
    labels, directions, linestyles, score_colors = roc_score_metadata(observable)
    comparisons = observable.get("roc", {}).get("comparisons", [])
    if len(comparisons) == 0:
        raise ValueError("plot: ROC observable requires at least one sample comparison")

    fig, axis = plt.subplots(1, figsize=observable.get("figsize", DEFAULT_FIGSIZE))
    xlim = observable.get("xlim", (1.0e-3, 1.0))
    ylim = observable.get("ylim", (0.0, 1.0))
    xscale = observable.get("xscale", "log")
    yscale = observable.get("yscale", "linear")
    curve_labels = []

    for comparison in comparisons:
        signal_index = int(comparison["signal"])
        background_index = int(comparison["background"])
        color_index = int(comparison.get("color_sample", signal_index))
        if not (0 <= signal_index < len(data) and 0 <= background_index < len(data) and 0 <= color_index < len(data)):
            raise IndexError("plot: ROC comparison sample index is out of range")

        signal = data[signal_index]
        background = data[background_index]
        comparison_label = comparison.get("label", f"{signal['label']} vs {background['label']}")
        comparison_color = comparison.get("color", data[color_index].get("color", None))
        comparison_linestyle = comparison.get("linestyle", None)

        for score_index, score_label in enumerate(labels):
            fpr, tpr = weighted_roc_curve(
                signal_scores=signal["rocdata"]["scores"][:, score_index],
                signal_weights=signal["rocdata"]["weights"],
                background_scores=background["rocdata"]["scores"][:, score_index],
                background_weights=background["rocdata"]["weights"],
                higher_is_signal=directions[score_index] == "higher",
            )
            curve_label = format_roc_curve_label(
                score_label=score_label,
                comparison_label=comparison_label,
                auc=weighted_roc_auc(fpr, tpr),
                roc_config=observable.get("roc", {}),
            )
            curve_labels.append(curve_label)

            plot_fpr = np.maximum(fpr, xlim[0]) if xscale == "log" else fpr
            color = score_colors[score_index]
            if color is None:
                color = comparison_color
            linestyle = comparison_linestyle
            if linestyle is None:
                linestyle = linestyles[score_index]
            axis.step(
                plot_fpr,
                tpr,
                where="post",
                color=color,
                linestyle=linestyle,
                linewidth=observable.get("roc", {}).get("linewidth", 1.4),
                label=curve_label,
            )

    if observable.get("roc", {}).get("show_diagonal", True):
        diagonal_x = np.geomspace(xlim[0], xlim[1], 200) if xscale == "log" else np.linspace(xlim[0], xlim[1], 200)
        axis.plot(diagonal_x, diagonal_x, color="0.55", linestyle=":", linewidth=0.8)

    axis.set_xscale(xscale)
    axis.set_yscale(yscale)
    axis.set_xlim(xlim)
    axis.set_ylim(ylim)
    axis.set_xlabel(observable["xlabel"])
    axis.set_ylabel(observable["ylabel"])
    axis.set_title(title, loc=title_loc, fontsize=6)
    axis.minorticks_on()
    axis.tick_params(top=True, bottom=True, right=True, left=True, which="both", direction="in", labelsize=9)
    draw_superplot_legend(
        fig=fig,
        ax=[axis],
        legend_labels=curve_labels,
        legend_position=legend_position,
        legend_properties=legend_properties,
        legend_bbox=legend_bbox,
    )
    return fig, [axis]


# Histogram ordinary observables or retain raw weighted ROC scores
def histmc(
    mcdata: dict,
    obs: dict,
    density: bool = False,
    density_uncertainty: str = "scaled",
    covariance_mode: str = "diagonal",
    scale: float | dict[str, float] | None = None,
    color: tuple = (0, 0, 1),
    label: str = "none",
    style: dict | None = None,
    verbose: bool = False,
):
    obj = {}
    style = dict(hist_style_step if style is None else style)

    for OBS in obs:
        observable = obs[OBS]
        x = mcdata["data"][OBS]
        weights = mcdata["weights"]
        kind = observable.get("kind", "histogram")

        if kind == "roc":
            labels, _, _, _ = roc_score_metadata(observable)
            obj[OBS] = {
                "rocdata": prepare_roc_data(x, weights, score_count=len(labels)),
                "color": color,
                "label": label,
                "obs": copy.deepcopy(observable),
            }
            continue
        if kind != "histogram":
            raise ValueError(f'plot: unknown observable kind "{kind}" for "{OBS}"')

        bins = observable["bins"]
        xlim = observable["xlim"]
        valid = observable.get("valid", None)
        category_histogram = type(bins) is dict
        sum_weights = weights.sum()
        if category_histogram:
            symmetrized_fill = obs[OBS]["symmetrized_fill"]
            members = hist.assign_2d_hyperbins(x=x[:, 0:2], hyperbins=bins, symmetrized=symmetrized_fill)
        subcategories = bins if category_histogram else (None,)

        # Histogram standard 1D input or each 3D hyperbin category through one path
        for SUB in subcategories:
            if category_histogram:
                x0_bins = bins[SUB][0]
                x1_bins = bins[SUB][1]
                bins_ = bins[SUB][2]
                valid_ = valid[SUB] if isinstance(valid, dict) and SUB in valid else None
                IND = members[SUB]
                values_ = x[IND, 2]
                weights_ = weights[IND]
                category_area = category_plane_area(
                    x0_bins=x0_bins, x1_bins=x1_bins, symmetrized=symmetrized_fill, obs_name=OBS
                )
                binwidth = category_area * hist.bins2binwidth(bins_)
                OBS_NAME = f"{OBS}_{hist.bins2txt(x0_bins)}_{hist.bins2txt(x1_bins)}"
            else:
                bins_ = bins
                valid_ = valid
                values_ = x
                weights_ = weights
                binwidth = hist.bins2binwidth(bins_)
                OBS_NAME = OBS

            # Compute the histogram and differential cross section over its configured range
            # Division by the full weight sum preserves the histogram overflow normalization
            counts, errs, _, cbins = hist.hist(x=values_, bins=bins_, weights=weights_)
            if not observable.get("differential", True):
                binwidth = np.ones_like(binwidth)
            binscale = hist.compute_binscale(
                binwidth=binwidth, totalweight=sum_weights, xsection_pb=mcdata["xsection_pb"]
            )
            # Preserve normalization exposure even when no event survives this selection
            if "sample_weights" in mcdata:
                binscale = hist.compute_binscale(binwidth, np.sum(mcdata["sample_weights"]), mcdata["sample_xsection_pb"])
            # Apply the optional external cross-section scale
            if scale is not None:
                binscale = binscale * (scale.get(OBS, 1.0) if isinstance(scale, dict) else scale)
            # Apply the observable MC normalization
            binscale = binscale * array.asarray(observable.get("mc_scale", 1.0), like=binscale)

            indices = hist.bin_indices(values_, bins_)
            selected = indices >= 0
            entries = np.bincount(indices[selected], minlength=len(bins_) - 1)
            events = {}
            if covariance_mode == "full":
                if "source" not in mcdata or "event_ids" not in mcdata:
                    raise ValueError("Full MC covariance requires source and event identities")
                ids = np.asarray(mcdata["event_ids"])
                ids = ids[IND] if category_histogram else ids
                events[mcdata["source"]] = (ids[selected], indices[selected], weights_[selected])
            elif covariance_mode != "diagonal":
                raise ValueError("MC covariance mode must be diagonal or full")

            obj[OBS_NAME] = {
                "hdata": hist.hobj(
                    counts,
                    errs,
                    bins_,
                    cbins,
                    binscale,
                    valid=valid_,
                    density=density,
                    differential=observable.get("differential", True),
                    density_uncertainty=density_uncertainty,
                    entries=entries,
                    mc_events=events,
                ),
                "hfunc": "hist",
                "color": color,
                "label": label,
                "style": style,
                "obs": copy.deepcopy(obs[OBS]),
            }
            if category_histogram:
                obj[OBS_NAME]["obs"]["units"] = category_differential_units(obs_config=obs[OBS], obs_name=OBS)
                obj[OBS_NAME]["obs"]["bins"] = copy.deepcopy(bins_)
                obj[OBS_NAME]["obs"]["xlim"] = copy.deepcopy(xlim[SUB])
                obj[OBS_NAME]["obs"]["xlabel_suffix"] = category_xlabel_suffix(
                    obs_config=obs[OBS], obs_name=OBS, category_bins=(x0_bins, x1_bins)
                )

            if density:
                obj[OBS_NAME]["obs"] = change2density_label(obj[OBS_NAME]["obs"])

            if verbose:
                print(f"histmc: integral = {obj[OBS_NAME]['hdata'].integral():0.2E} ({OBS_NAME})")

    return obj


# Histogram HEPData with the selected iceplot drawing primitive
def histhepdata(
    hepdata: dict,
    obs: dict,
    scale: float | dict[str, float] | None = None,
    density: bool = False,
    density_uncertainty: str = "shape",
    MC_XS_SCALE: float = 1e12,
    label: str = "Data",
    hfunc: str = "hist",
    style: dict | None = None,
    verbose: bool = False,
):

    # Over all observables
    obj = {}
    if hfunc == "hist":
        default_style = hist_style_step
    elif hfunc == "errorbar":
        default_style = errorbar_style
    else:
        raise ValueError(f'plot: unknown HEPData drawing primitive "{hfunc}"')
    style = dict(default_style if style is None else style)

    for OBS in obs:
        dataset = hepdata[OBS]
        bins = dataset["bins"]
        category_histogram = type(bins) is dict
        subcategories = bins if category_histogram else (None,)

        # Treat a simple 1D histogram as one category without hyperbin metadata
        for SUB in subcategories:
            if category_histogram:
                x0_bins = copy.deepcopy(bins[SUB][0])
                x1_bins = copy.deepcopy(bins[SUB][1])
                bins_ = copy.deepcopy(bins[SUB][2])
                y_source = dataset["y"][SUB]
                yerr_source = dataset["y_err"][SUB]
                binwidth_ = copy.deepcopy(dataset["binwidth"][SUB][2])
                cbins_ = copy.deepcopy(dataset["x"][SUB][2])
                xlim_ = copy.deepcopy(dataset["xlim"][SUB])
                valid = dataset.get("valid")
                valid_ = valid[SUB] if isinstance(valid, dict) and SUB in valid else None
                uncertainty_sources = dataset.get("uncertainties", {}).get(SUB, [])
                OBS_NAME = f"{OBS}_{hist.bins2txt(x0_bins)}_{hist.bins2txt(x1_bins)}"
            else:
                bins_ = bins
                y_source = dataset["y"]
                yerr_source = dataset["y_err"]
                binwidth_ = dataset["binwidth"]
                cbins_ = dataset["x"]
                xlim_ = dataset["xlim"]
                valid_ = dataset.get("valid")
                uncertainty_sources = dataset.get("uncertainties", [])
                OBS_NAME = OBS

            point_data = not cross_section_observable(obs[OBS])
            binscale_ = dataset["scale"] if point_data else dataset["scale"] * MC_XS_SCALE

            if "uncertainties" not in dataset:
                uncertainty_sources = [
                    {
                        "name": "total",
                        "category": "combined",
                        "correlation": "uncorrelated",
                        "effect": "additive",
                        "up": np.abs(np.asarray(yerr_source, dtype=float)),
                        "down": np.abs(np.asarray(yerr_source, dtype=float)),
                    }
                ]

            # Convert cross sections or counts per bin before unit density normalization
            if density and not obs[OBS].get("differential", True):
                transform = np.diag(1.0 / np.asarray(binwidth_, dtype=float))
                y_source, yerr_source = transform @ y_source, transform @ yerr_source
                uncertainty_sources = [uncertainty.linear_source(source, transform) for source in uncertainty_sources]
            y_, yerr_ = copy.deepcopy(y_source), copy.deepcopy(yerr_source)

            # Additional scale factor
            if scale is not None and not point_data:
                binscale_ *= scale.get(OBS, 1.0) if isinstance(scale, dict) else scale

            plotted_sources = uncertainty.transform_sources(
                uncertainty_sources,
                values=np.asarray(y_source, dtype=float),
                binwidth=np.asarray(binwidth_, dtype=float),
                binscale=binscale_,
                density=density,
                density_uncertainty=density_uncertainty,
                valid=valid_,
            )

            # Density integral 1 over the histogram bins
            if density:
                y_, _ = hist.normalize_bin_heights(values=y_, errors=yerr_, binwidth=binwidth_, valid=valid_)
                covariance = sum(
                    (source_covariance(source) for source in plotted_sources),
                    start=np.zeros((len(y_), len(y_)), dtype=float),
                )
                yerr_ = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
                binscale_ = 1.0

            obj[OBS_NAME] = {
                "hdata": hist.hobj(y_, yerr_, bins_, cbins_, binscale_, valid=valid_,
                                   differential=density or obs[OBS].get("differential", True)),
                "hfunc": hfunc,
                "color": (0, 0, 0),
                "label": label,
                "style": style,
                "obs": copy.deepcopy(obs[OBS]),
                "fitw": dataset["fitw"],
                "uncertainties": plotted_sources,
                "stat_errs": np.sqrt(
                    np.clip(
                        np.diag(
                            sum(
                                (
                                    source_covariance(source)
                                    for source in plotted_sources
                                    if source["category"] == "statistical"
                                ),
                                start=np.zeros((len(y_), len(y_)), dtype=float),
                            )
                        ),
                        0.0,
                        None,
                    )
                ),
            }
            # Preserve quoted asymmetric margins in plots independently of the symmetric fit covariance
            if not density or density_uncertainty == "scaled":
                error_scale = binscale_
                if density:
                    error_scale = 1.0 / np.sum(np.where(valid_ if valid_ is not None else True, y_source, 0.0) * binwidth_)
                for key, prefix in (("total_errs", "y_err"), ("stat_errs", "y_err_stat")):
                    if f"{prefix}_down" in dataset and f"{prefix}_up" in dataset:
                        errors = [dataset[f"{prefix}_{side}"] for side in ("down", "up")]
                        if category_histogram:
                            errors = [item[SUB] for item in errors]
                        obj[OBS_NAME][key] = np.asarray(errors, dtype=float) * error_scale
            if category_histogram:
                obj[OBS_NAME]["obs"]["units"] = category_differential_units(obs_config=obs[OBS], obs_name=OBS)

            if density:
                obj[OBS_NAME]["obs"] = change2density_label(obj[OBS_NAME]["obs"])

            # Set the common histogram display parameters
            obj[OBS_NAME]["obs"]["bins"] = copy.deepcopy(bins_)
            obj[OBS_NAME]["obs"]["xlim"] = xlim_
            if category_histogram:
                obj[OBS_NAME]["obs"]["xlabel_suffix"] = category_xlabel_suffix(
                    obs_config=obs[OBS], obs_name=OBS, category_bins=(x0_bins, x1_bins)
                )

            if verbose:
                print(
                    f"histhepdata: integral = {obj[OBS_NAME]['hdata'].integral():0.2E}, bins(min,max) = [{np.min(bins_):0.2f}, {np.max(bins_):0.2f}] [{OBS_NAME}]"
                )

    return obj


# Normalize cumulative cross-section histograms with one shared final total
def normalize_cumulative_process_stack(histograms: list[hist.hobj], density_uncertainty: str = "scaled") -> None:
    if density_uncertainty not in {"scaled", "shape"}:
        raise ValueError("plot: density uncertainty must be either 'scaled' or 'shape'")

    binwidth = np.asarray(histograms[-1].measure, dtype=float)
    denominator_counts = np.asarray(histograms[-1].counts_scaled, dtype=float) * binwidth
    denominator_errs = np.asarray(histograms[-1].errs_scaled, dtype=float) * binwidth
    common_valid = np.asarray(histograms[-1].valid, dtype=np.bool_)
    normalization = float(np.sum(np.where(common_valid, denominator_counts, 0.0)))
    if not np.isfinite(normalization) or normalization <= 0.0:
        raise ValueError("plot: density process stacking requires a finite positive total integral")
    for histogram in histograms:
        histogram.counts = np.asarray(histogram.counts_scaled, dtype=float) * binwidth
        histogram.errs = np.asarray(histogram.errs_scaled, dtype=float) * binwidth
        histogram.binscale = 1.0
        histogram.valid = np.asarray(histogram.valid, dtype=np.bool_) & common_valid
        histogram.density = True
        histogram.density_uncertainty = density_uncertainty
        histogram.density_denominator_counts = denominator_counts.copy()
        histogram.density_denominator_errs = denominator_errs.copy()


# Build filled cumulative plotting records and their summed process total
def build_filled_process_stack_records(
    records: list[dict], density: bool = False, density_uncertainty: str = "scaled"
) -> tuple[list[dict], dict]:
    if not records:
        raise ValueError("plot: process stacking requires at least one histogram record")

    process_histograms = []
    for index, record in enumerate(records):
        if "rocdata" in record or "hdata" not in record:
            raise ValueError("plot: process stacking only supports ordinary histograms")
        if record.get("hfunc") != "hist":
            raise ValueError("plot: process stacking requires MC histogram records")
        values = np.asarray(record["hdata"].counts_scaled, dtype=float)
        if np.any(values < 0.0):
            raise ValueError(f"plot: process stacking requires nonnegative bins, process {index} has negatives")
        process_histograms.append(record["hdata"])

    cumulative_histograms = hist.hobj.stack_independent_processes(process_histograms)
    if density:
        normalize_cumulative_process_stack(cumulative_histograms, density_uncertainty=density_uncertainty)

    cumulative_records = []
    for record, cumulative in zip(records, cumulative_histograms, strict=True):
        visual = record.copy()
        visual["hdata"] = cumulative
        visual["style"] = {**hist_style_fill, "alpha": 0.75}
        visual["_stack_component"] = True
        if density:
            visual["obs"] = change2density_label(copy.deepcopy(record["obs"]))
        cumulative_records.append(visual)

    total_record = records[-1].copy()
    total_record["hdata"] = cumulative_histograms[-1]
    total_record["label"] = "MC total"
    total_record["color"] = "black"
    total_record["style"] = dict(hist_style_step)
    total_record.pop("_stack_component", None)
    if density:
        total_record["obs"] = change2density_label(copy.deepcopy(total_record["obs"]))
    return list(reversed(cumulative_records)), total_record


# Fuse histogram and ROC outputs from independent chunks of the same source
def fuse_worker_chunk_outputs(chunk_outputs: list) -> list:
    if not chunk_outputs:
        raise ValueError("plot: independent chunk fusion requires at least one worker output")

    reference = chunk_outputs[0]
    for chunk_index, chunk in enumerate(chunk_outputs[1:], start=1):
        if len(chunk) != len(reference):
            raise ValueError(f"plot: chunk {chunk_index} has {len(chunk)} sets, expected {len(reference)}")
        for target, source in zip(reference, chunk, strict=True):
            if set(target) != set(source):
                raise ValueError(f"plot: chunk {chunk_index} has incompatible observables")

    combined = copy.deepcopy(reference)
    for set_index, target in enumerate(combined):
        for name, target_record in target.items():
            records = [chunk[set_index][name] for chunk in chunk_outputs]
            record_types = [
                "roc" if "rocdata" in record else "histogram" if "hdata" in record else "malformed"
                for record in records
            ]
            if len(set(record_types)) != 1 or record_types[0] == "malformed":
                raise ValueError(f"plot: chunks disagree on observable {name!r} type")
            if record_types[0] == "roc":
                rocdata = copy.deepcopy(records[0]["rocdata"])
                for record in records[1:]:
                    rocdata = fuse_roc_data(rocdata, record["rocdata"])
                target_record["rocdata"] = rocdata
            else:
                target_record["hdata"] = hist.hobj.fuse_independent_chunks([record["hdata"] for record in records])
    return combined

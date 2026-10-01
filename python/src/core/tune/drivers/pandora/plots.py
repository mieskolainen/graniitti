# Particle flow diagnostic plots for the Pandora driver
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import pathlib
import re
from typing import Any

import numpy as np

from core.io.files import ensure_dir
from core.tune.drivers.pandora import diagnostics


# Compute pyplot configured for noninteractive validation outputs
def _pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


# Save and close one plot in the canonical Pandora raster format
def _save_figure(fig: Any, output_path: pathlib.Path, style, **kwargs) -> None:
    ensure_dir(output_path.parent)
    fig.savefig(output_path, **style["savefig"], **kwargs)
    _pyplot().close(fig)


# Write one histogram with central interval and target markers
def _plot_histogram(*, detail, output_path, title, xlabel, style):
    cfg = style["histogram"]
    colors = style["colors"]

    stats = detail["stats"]
    fig, ax = _pyplot().subplots(**cfg["figure"])
    ax.stairs(detail["counts"], detail["bins"], **cfg["fill"], color=colors[1])
    annotation = f"N = {detail['n']}"
    mean, sigma = stats["mean"], stats["sigma"]
    if mean is not None:
        ax.axvspan(stats["low"], stats["high"], **cfg["interval"], label="smallest 90% interval")
        ax.axvline(mean, color=colors[0], label="Mean90")
        mode = detail["mode"]
        annotation += f"\nMean90 = {mean:{cfg['format']}}\n" + rf"$\sigma_{{90}}$ = {sigma:{cfg['format']}}" + f"\nMode = {mode:{cfg['format']}}"
    if detail.get("target") is not None:
        ax.axvline(detail["target"], color=colors[2], label="target")
    if "calibration" in detail:
        calibration = detail["calibration"]
        tolerance = calibration["tolerance"]
        ax.axvline(calibration["mean"], **style["reference"], label="Full mean")
        ax.axvspan(1 - tolerance, 1 + tolerance, color=colors[2], alpha=cfg["calibration_alpha"], label="calibration interval")
        status = "passed" if calibration["passed"] else "outside tolerance"
        annotation += f"\n\nFull mean = {calibration['mean']:{cfg['format']}}\nCalibration:\n{status}"
    ax.text(*cfg["annotation_position"], annotation, transform=ax.transAxes, **cfg["annotation"])
    ax.set(title=title, xlabel=xlabel, ylabel="Weighted events")
    ax.legend(**style["legend"])
    ax.grid(**style["grid"])
    fig.tight_layout()
    _save_figure(fig, output_path, style)


# Compute class indices in the physics display order for confusion plots
def _pf_confusion_display_indices(class_names: tuple[str, ...]) -> np.ndarray:
    lookup = {name: index for index, name in enumerate(class_names)}
    ordered = [lookup[name] for name in diagnostics.PF_CONFUSION_PLOT_CLASS_NAMES if name in lookup]
    ordered.extend(index for index, name in enumerate(class_names) if name not in diagnostics.PF_CONFUSION_CLASS_NAMES)
    return np.asarray(ordered, dtype=int)


# Reorder one rectangular confusion-matrix payload into the plot display convention
def _pf_confusion_display_payload(row_names, column_names, matrices):
    rows, cols = map(_pf_confusion_display_indices, (row_names, column_names))
    reordered = {key: np.asarray(value)[np.ix_(rows, cols)] if key.startswith(("row_", "column_")) else value
                 for key, value in matrices.items() if key not in {"row_names", "column_names"}}
    reordered.update(truth_denominator=matrices["truth_denominator"][rows],
                     reco_denominator=matrices["reco_denominator"][cols])
    return tuple(row_names[i] for i in rows), tuple(column_names[i] for i in cols), reordered


# Compute the filename key and physics label from the same energy interval
def _energy_bin(low, high):
    if low is None or high is None:
        return "inclusive", "inclusive"
    key = f"{low:g}_{high:g}".replace(".", "p").replace("-", "m")
    return re.sub(r"[^A-Za-z0-9_]+", "_", key).strip("_"), rf"$E = {low:g} - {high:g}$ GeV"


# Write one split-cell PF class confusion matrix for a selected energy range
def _plot_pf_class_confusion_matrix(*, matrices, low, high, output_path, style):
    from matplotlib.patches import Rectangle
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    cfg = style["class_confusion"]

    plt = _pyplot()
    rows, cols, matrices = _pf_confusion_display_payload(matrices["row_names"], matrices["column_names"], matrices)
    fig, ax = plt.subplots(**cfg["figure"])
    styles = (("row", cfg["cmaps"][0], -0.5, r"Efficiency ($N_{ij}/N_{\mathrm{gen},i}$)", cfg["colorbar_pads"][0]),
        ("column", cfg["cmaps"][1], 0.0, r"Purity ($N_{ij}/N_{\mathrm{reco},j}$)", cfg["colorbar_pads"][1]))
    for prefix, color, offset, _, _ in styles:
        cmap = plt.get_cmap(color)
        for (row, col), value in np.ndenumerate(matrices[prefix + "_rate"]):
            available = (matrices["truth_denominator"][row] > 0.0 if prefix == "row"
                         else matrices["reco_denominator"][col] > 0.0)
            ax.add_patch(Rectangle((col - 0.5, row + offset), 1.0, 0.5,
                                  facecolor=cmap(value) if available else cfg["empty_color"], **cfg["cell"]))
            down, up = value - matrices[prefix + "_lower"][row, col], matrices[prefix + "_upper"][row, col] - value
            fmt = cfg["rate_format"]
            text = rf"${value:{fmt}}^{{+{up:{fmt}}}}_{{-{down:{fmt}}}}$" if available else r"$\mathrm{n/a}$"
            ax.text(col, row + offset + cfg["text_offsets"][prefix != "row"], text,
                    **cfg["text"], color=cfg["text_colors"][int(value > cfg["text_threshold"])])
    for row, name in enumerate(rows):
        if name in cols:
            ax.add_patch(Rectangle((cols.index(name) - 0.5, row - 0.5), 1, 1,
                                  **cfg["diagonal"]))
    for axis, names, side in (("x", cols, "reco"), ("y", rows, "truth")):
        particle = "reco" if side == "reco" else "gen"
        labels = ["no PFO\n(lost)" if name == "lost" else "no Gen\n(fake)" if name == "fake" else name.replace("_", " ") +
                  f"\n$N_{{\\mathrm{{{particle}}}}} = {matrices[side + '_denominator'][i]:{cfg['count_format']}}$"
                  for i, name in enumerate(names)]
        getattr(ax, f"set_{axis}ticks")(np.arange(len(names)), labels, **({"rotation": cfg["tick_rotation"], "ha": "right"} if axis == "x" else {}))
    ax.set(xlim=(-0.5, len(cols) - 0.5), ylim=(len(rows) - 0.5, -0.5), title=_energy_bin(low, high)[1],
           xlabel=r"Reco (PFO) class, purity binned in $E_{\mathrm{reco}}$",
           ylabel=r"Gen class, efficiency binned in $E_{\mathrm{gen}}$")
    divider = make_axes_locatable(ax)
    for _, cmap, _, label, pad in styles:
        scale = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
        fig.colorbar(scale, cax=divider.append_axes("right", size=cfg["colorbar_size"], pad=pad)).set_label(label)
    _save_figure(fig, output_path, style, bbox_inches="tight")


# Write one PF object efficiency and purity profile versus one object variable
def _plot_pf_object_performance(*, data, bins, particle_class, output_path, axis, style):
    cfg = style["performance"]
    colors = style["colors"]

    _, _, xlabel, _, log_x = diagnostics.PF_OBJECT_PERFORMANCE_AXES[axis]
    fig, (counts, ax) = _pyplot().subplots(2, 1, **cfg["figure"])
    for side, label, color in (("truth", "Efficiency", colors[0]), ("reco", "Purity", colors[1])):
        result = data["rates"][side]
        selected = result["nonzero"]
        ax.errorbar(result["centers"][selected], result["rate"][selected], yerr=result["yerr"][:, selected],
                    **cfg["errorbar"], color=color, label=label)
        counts.stairs(result["denominator"], bins, color=color, label="Gen" if side == "truth" else "Reco (PFO)")
    ax.set(xlabel=xlabel, ylim=(0.0, 1.0), yticks=style["rate_ticks"], xlim=(bins[0], bins[-1]))
    counts.set(title="PF " + particle_class.replace("_", " "), ylabel="Counts")
    for panel in (ax, counts):
        panel.legend(**style["legend"])
        panel.grid(**style["grid"])
        if log_x:
            panel.set_xscale("log")
    _save_figure(fig, output_path, style)


# Write energy bias and resolution versus generated energy for one truth class
def plot_resolution(*, profiles, bins, particle_class, output_path, style):
    cfg = style["resolution"]
    colors = style["colors"]

    fig, axes = _pyplot().subplots(2, 2, **cfg["figure"])
    panels = (("median", "Median bias", r"$Q_{50}(\delta_E)$"),
              ("q68", "Central 68% half width", r"$(Q_{84}-Q_{16})/2$"),
              ("mean90", "Mean bias in the 90% interval", r"$\langle\delta_E\rangle_{90}$"),
              ("sigma90", r"$\sigma_{90}$", r"$\sigma_{90}(\delta_E)$"))
    for result, label, color in zip(profiles, ("All assigned matches", "Correct particle type"),
                                    (colors[0], colors[1]), strict=True):
        for ax, (key, _, _) in zip(axes.flat, panels, strict=True):
            valid = np.isfinite(result[key])
            ax.errorbar(result["centers"][valid], result[key][valid], yerr=result["errors"][key][valid],
                        **cfg["errorbar"], color=color, label=label)
    for ax, (key, title, ylabel) in zip(axes.flat, panels, strict=True):
        has_matches = any(np.any(np.isfinite(line.get_ydata())) for line in ax.lines)
        ax.set(title=title, xlim=(bins[0], bins[-1]), xscale="log", ylabel=ylabel)
        ax.tick_params(axis="x", labelbottom=True)
        if key in ("median", "mean90"):
            ax.axhline(0, **style["reference"])
        else:
            ax.set_ylim(bottom=0)
        if not has_matches:
            ax.text(0.5, 0.5, "No assigned matches", transform=ax.transAxes, ha="center", va="center")
    axes[0, 0].legend(**style["legend"])
    fig.suptitle(particle_class.replace("_", " ")
                 + r": $\delta_E=(E_{\mathrm{reco}}-E_{\mathrm{gen}})/E_{\mathrm{gen}}$")
    fig.supxlabel(r"Generated particle energy $E_{\mathrm{gen}}$ [GeV]")
    _save_figure(fig, output_path, style)


# Write differential matching outcomes and angular distributions for one truth class
def _plot_pf_matching(*, data, bins, particle_class, axis, max_angle, output_path, style):
    from matplotlib.ticker import MaxNLocator

    cfg = style["matching"]
    colors = style["colors"]

    _, _, xlabel, _, log_x = diagnostics.PF_OBJECT_PERFORMANCE_AXES[axis]
    fig, axes = _pyplot().subplots(2, 2, **cfg["figure"])
    ax = axes[0, 0]
    for result, label, color in zip(data["outcomes"],
        ("Correct class", "Wrong class", "No candidate", "Candidate assigned elsewhere"),
        (colors[2], colors[0], colors[1], colors[3]), strict=True):
        selected = result["nonzero"]
        ax.errorbar(result["centers"][selected], result["rate"][selected], yerr=result["yerr"][:, selected],
                    **cfg["errorbar"], label=label, color=color)
    ax.set(title="Matching outcome", ylabel="Fraction of Gen particles", ylim=(0, 1), yticks=style["rate_ticks"])
    ax.legend(**cfg["legend"])
    ax.grid(**style["grid"])
    for ax, key, title, ylabel in (
        (axes[0, 1], "candidates", "Candidates within 3D angle cut", "Reco candidates / Gen particle"),
        (axes[1, 0], "angle", "Assigned 3D angle", r"$\theta_{ij}$ [rad]"),
        (axes[1, 1], "dr", r"Assigned $\Delta R$", r"$\Delta R_{ij}=\sqrt{(\Delta\eta)^2+(\Delta\phi)^2}$")):
        histogram = data["matching"][key]
        if histogram is not None:
            counts, xbins, ybins = histogram
            image = ax.pcolormesh(xbins, ybins, counts.T, cmap=cfg["cmap"])
            fig.colorbar(image, ax=ax, label="Counts")
        else:
            ax.text(0.5, 0.5, "No particles" if key == "candidates" else "No assigned matches",
                    transform=ax.transAxes, ha="center", va="center")
        if key == "candidates":
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        elif key == "angle":
            ax.axhline(max_angle, color=colors[3], linestyle=cfg["cut_linestyle"], label="3D angle cut")
            ax.set_ylim(0, max_angle * cfg["angle_margin"])
            ax.legend(**cfg["legend"])
        ax.set(title=title, ylabel=ylabel)
    for ax in axes.flat:
        ax.set(xlabel="Gen " + xlabel, xlim=(bins[0], bins[-1]))
        if log_x:
            ax.set_xscale("log")
    fig.suptitle("3D angle matching: Gen " + particle_class.replace("_", " "))
    _save_figure(fig, output_path, style)


# Write composition diagnostics from one set of weighted event profiles
def _plot_pf_composition(composition, root, *, low, high, style):
    from matplotlib.ticker import MaxNLocator

    cfg = style["composition"]
    colors = style["colors"]

    key, label = _energy_bin(low, high)
    suffix, title = ("", "") if low is None else ("_" + key, "\n" + label)
    names = composition["class_names"]
    styles = (("gen", cfg["formats"][0], colors[0], "Gen"), ("pfo", cfg["formats"][1], colors[1], "Reco (PFO)"))
    fig, axes = _pyplot().subplots(2, 3, **cfg["figure"])
    fig2d, axes2d = _pyplot().subplots(2, 3, **cfg["figure"])
    for ax, ax2d, name in zip(axes.flat, axes2d.flat, (*names, "total"), strict=True):
        histogram = composition["multiplicity"][name]
        bins = histogram["bins"]
        for counts, (_, _, color, label) in zip(histogram["counts"], styles, strict=True):
            ax.stairs(counts, bins, color=color, label=label)
        ax.set(title=name.replace("_", " "), xlabel="multiplicity / event", ylabel="Weighted events",
               xlim=(bins[0], bins[-1]))
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid(**cfg["grid"])
        image = ax2d.pcolormesh(bins, bins, histogram["joint"].T, cmap=cfg["cmap"])
        fig2d.colorbar(image, ax=ax2d, label="events")
        ax2d.plot([bins[0], bins[-1]], [bins[0], bins[-1]], **cfg["diagonal"])
        ax2d.set(title=name.replace("_", " "), xlabel="Gen multiplicity / event",
                 ylabel="Reco (PFO) multiplicity / event", xlim=(bins[0], bins[-1]),
                 ylim=(bins[0], bins[-1]), aspect="equal")
        ax2d.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax2d.yaxis.set_major_locator(MaxNLocator(integer=True))
    axes.flat[0].legend(**style["legend"])
    for figure, name in ((fig, "particle_multiplicity"), (fig2d, "particle_multiplicity_2d")):
        figure.suptitle("PF Particle Multiplicity" + title)
        _save_figure(figure, root / f"{name}{suffix}.png", style)
    fig, ax = _pyplot().subplots(**cfg["fraction_figure"])
    for side, marker, color, label in styles:
        means, widths = composition["fractions"][side]
        ax.errorbar(np.arange(len(names)), means, yerr=widths, fmt=marker, **cfg["errorbar"], color=color, label=label)
    ax.set_xticks(np.arange(len(names)), [name.replace("_", " ") for name in names], rotation=cfg["tick_rotation"], ha="right")
    ax.set(ylim=(0, None), ylabel="momentum fraction / event", title="PF Particle Momentum Fraction" + title)
    ax.yaxis.set_major_locator(MaxNLocator(**cfg["fraction_ticks"]))
    ax.grid(axis="y", **style["grid"])
    ax.legend(**style["legend"])
    _save_figure(fig, root / f"particle_momentum_fraction{suffix}.png", style)


# Write calibration and the raw terms entering the logarithmic PF risk
def _plot_pf_loss(*, detail, output_path, style):
    cfg = style["loss"]
    colors = style["colors"]

    terms = detail["loss_components"]
    fig, ax = _pyplot().subplots(**cfg["figure"])
    ax.bar(list(terms), list(terms.values()), color=colors[1])
    state = "passed" if detail["calibration"]["passed"] else "outside tolerance"
    ax.set(ylabel="Scaled loss before combination", title="PF loss = calibration + log(1 + PF risk)\n"
           f"Cost = {detail['cost']:{cfg['format']}}, calibration {state}")
    ax.grid(axis="y", **style["grid"])
    fig.tight_layout()
    _save_figure(fig, output_path, style)


# Write signed loss contributions by particle type or particle energy
def _plot_loss_contributions(data, output_path, *, energy, style):
    cfg = style["contributions"]
    colors = style["colors"]

    names = [name.replace("_", " ") for name in data["classes"]]
    fig, axes = _pyplot().subplots(2, 3, **cfg["figure"], sharex=energy, sharey=not energy)
    for ax, (component, values) in zip(axes.flat, data["energy" if energy else "particle"].items(), strict=True):
        if energy:
            for i, name in enumerate(names):
                ax.stairs(values[i], data["bins"], baseline=None, color=colors[i], label=name,
                          **cfg["line"], linestyle=cfg["linestyles"][i % len(cfg["linestyles"])])
            ax.set(xscale="log", xlim=(data["bins"][0], data["bins"][-1]))
            ax.tick_params(axis="x", labelbottom=True)
            ax.axhline(0, **cfg["reference"])
        else:
            bars = ax.barh(np.arange(len(names)), values, **cfg["bar"], color=[colors[i] for i in range(len(names))])
            ax.bar_label(bars, **cfg["bar_label"])
            ax.set(yticks=np.arange(len(names)), yticklabels=names, ylim=(len(names) - 0.5, -0.5))
            ax.margins(x=cfg["bar_margin"])
            ax.axvline(0, **cfg["reference"])
            ax.tick_params(axis="y", length=0)
        ax.set_title(component.capitalize(), **cfg["title"])
        ax.set_axisbelow(True)
        ax.grid(axis="y" if energy else "x", **cfg["grid"])
        ax.tick_params(**cfg["ticks"])
    if energy:
        handles, labels = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, ncol=len(names), **cfg["legend"])
        fig.supylabel("Signed loss contribution per energy bin", **cfg["label"])
        fig.supxlabel("Gen particle energy [GeV] (unmatched PFOs use reconstructed energy)", **cfg["energy_label"])
    else:
        fig.supxlabel("Signed loss contribution", **cfg["label"])
    fig.tight_layout(rect=cfg["energy_rect" if energy else "particle_rect"], **cfg["layout"])
    _save_figure(fig, output_path, style)


# Write one labeled energy confusion heatmap with GEN rows and PFO columns
def _plot_pf_confusion_heatmap(*, matrix, row_names, column_names, output_path, title,
                               xlabel="Reco (PFO) class", ylabel="Gen class", symmetric_zero=False,
                               colorbar_label="", style):
    from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm

    cfg = style["heatmap"]
    rows, cols = map(_pf_confusion_display_indices, (row_names, column_names))
    values = np.asarray(matrix, dtype=float)[np.ix_(rows, cols)].T
    row_names, column_names = tuple(column_names[i] for i in cols), tuple(row_names[i] for i in rows)
    # Compute the colormap and normalization, retaining white for zero energy and a symmetric residual scale
    finite = values[np.isfinite(values)]
    if symmetric_zero:
        limit = max(float(np.max(np.abs(finite))) if finite.size else 0.0, cfg["empty_scale"])
        # Keep both residual signs light enough for black cell labels
        cmap = LinearSegmentedColormap.from_list("pf_confusion_purple_white_green", cfg["residual_colors"])
        norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    else:
        limit = max(float(np.max(np.clip(finite, 0, None))) if finite.size else 0.0, cfg["empty_scale"])
        cmap = LinearSegmentedColormap.from_list("pf_confusion_white_to_red", cfg["energy_colors"])
        norm = Normalize(vmin=0, vmax=limit, clip=True)
    fig, ax = _pyplot().subplots(**cfg["figure"])
    if values.ndim == 2 and values.size:
        image = ax.imshow(values, cmap=cmap, aspect="auto", norm=norm)
        for (row, col), value in np.ndenumerate(values):
            ax.text(col, row, f"{value:{cfg['format']}}", **cfg["text"])
        # Label absent truth and PFO entries explicitly on their respective axes
        xlabels = ["no PFO\n(lost)" if name == "lost" else name.replace("_", " ") for name in column_names]
        ylabels = ["no Gen\n(fake)" if name == "fake" else name.replace("_", " ") for name in row_names]
        ax.set_xticks(np.arange(values.shape[1]), xlabels, rotation=cfg["tick_rotation"], ha="right")
        ax.set_yticks(np.arange(values.shape[0]), ylabels)
        fig.colorbar(image, ax=ax, **cfg["colorbar"]).set_label(colorbar_label)
    else:
        ax.text(0.5, 0.5, "No confusion matrix", transform=ax.transAxes, ha="center", va="center")
    ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    _save_figure(fig, output_path, style)


# Write the particle flow objective and its weighted physics diagnostics
def write_plots(data, fig_dir):
    from core.plot.plot import colors

    style = data["settings"]["plots"]
    style = dict(style, colors=[colors(i) for i in style["imperial_colors"]])
    detail, root = data["detail"], fig_dir / "pflow"
    paths = {}

    # Record each diagnostic under its existing filename and summary key
    def output(name, subdir="", key=None):
        path = root / subdir / (name + ".png")
        ensure_dir(path.parent)
        paths[key or "pflow_" + name] = str(path)
        return path

    for prefix, name, title, xlabel, target in (
        ("", "pf_total_energy_relative", "PF Total Energy Relative Response", r"$E_{\mathrm{reco}}/E_{\mathrm{gen}}$", 1.0),
        ("visible_total_energy_", "pf_total_energy", "PFO Reconstructed Energy", r"$E_{\mathrm{reco}}$ [GeV]",
         detail["truth_energy_stats"]["mean"]),
        ("truth_energy_", "gen_total_energy", "Generated Visible Energy", r"$E_{\mathrm{gen}}$ [GeV]", None)):
        profile = dict(data["histograms"][prefix], stats=detail[prefix + "stats"], target=target)
        if not prefix:
            profile["calibration"] = detail["calibration"]
        _plot_histogram(style=style, detail=profile, output_path=output(name, key=name), title=title, xlabel=xlabel)
    paths["pflow"] = paths["pf_total_energy_relative"]
    _plot_pf_loss(style=style, detail=detail, output_path=output("loss_components", "components"))
    if "loss_contributions" in data:
        for energy, name in ((False, "loss_by_particle"), (True, "loss_vs_energy")):
            _plot_loss_contributions(data["loss_contributions"], output(name, "components"), energy=energy, style=style)
    confusion = detail["confusion"]
    for field, suffix, title, label in (
        ("matrix", "", "Energy Confusion", diagnostics.PF_ENERGY_CONFUSION_COLORBAR_LABEL),
        ("residual", "_residual", "Energy Confusion Residual", diagnostics.PF_CONFUSION_RESIDUAL_COLORBAR_LABEL)):
        _plot_pf_confusion_heatmap(style=style, matrix=confusion[field], row_names=confusion["row_names"], column_names=confusion["column_names"],
            output_path=output("energy_confusion" + suffix, "confusion"), title=title,
            symmetric_zero=field == "residual", colorbar_label=label)
    for particle, prepared in data["particles"].items():
        plot_resolution(style=style, profiles=prepared["resolution"], bins=data["bins"]["energy"], particle_class=particle,
            output_path=output(f"{particle}_resolution_vs_energy", "resolution"))
        for axis, values in prepared["axes"].items():
            _plot_pf_object_performance(style=style, data=values, bins=data["bins"][axis], particle_class=particle, axis=axis,
                output_path=output(f"{particle}_performance_vs_{axis}", "object_performance"))
            _plot_pf_matching(style=style, data=values, bins=data["bins"][axis], particle_class=particle, axis=axis,
                max_angle=detail["loss"]["matching_max_angle_rad"],
                output_path=output(f"{particle}_matching_vs_{axis}", "matching"))
    for interval in data["slices"]:
        low, high = interval["low"], interval["high"]
        suffix = "" if low is None else "_" + _energy_bin(low, high)[0]
        for name in ("particle_multiplicity", "particle_multiplicity_2d", "particle_momentum_fraction"):
            output(name + suffix, "composition")
        _plot_pf_composition(interval["composition"], root / "composition", low=low, high=high, style=style)
        _plot_pf_class_confusion_matrix(style=style, matrices=interval["matrices"], low=low, high=high,
            output_path=output("class_confusion_" + _energy_bin(low, high)[0], "class_confusion"))
    return paths

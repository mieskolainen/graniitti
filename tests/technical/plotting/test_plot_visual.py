# Basic operation tests for the public plotting helpers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
import pytest
from core.plot.plot import (
    build_filled_process_stack_records,
    create_axes,
    data_total_style,
    errorbar_style,
    hist_style_step,
    ordered_legend,
    plot_horizontal_line,
    plot_style,
    superplot,
)
from core.stats.hist import hist, hist_obj, hobj
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.container import ErrorbarContainer


# Render a dashed step histogram with its uncertainty fill
def test_dashed_hist_error_renders(tmp_path):
    bins = np.array([0.0, 1.0, 2.0])
    record = {
        "hdata": hobj(
            counts=np.array([2.0, 1.0]),
            errs=np.array([0.2, 0.1]),
            bins=bins,
            cbins=np.array([0.5, 1.5]),
        ),
        "hfunc": "hist",
        "color": "black",
        "label": "Reference",
        "style": {**hist_style_step, "ls": "--"},
        "obs": {
            "xlim": (0.0, 2.0),
            "ylim": None,
            "xlabel": "$x$",
            "ylabel": "$d\\sigma/dx$",
            "units": {"x": "GeV", "y": "pb"},
            "ylim_ratio": (0.0, 2.0),
        },
    }

    fig, _ = superplot([record], ratio_plot=False)
    output = tmp_path / "dashed.pdf"
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)

    assert output.stat().st_size > 0


# Check an excluded interval does not erase later histogram step lines
def test_gapped_step_hist_draws_finite_run():
    bins = np.arange(5.0)
    record = {
        "hdata": hobj(
            counts=np.array([0.6, 0.5, 0.25, 0.18]),
            errs=np.array([0.05, 0.04, 0.03, 0.02]),
            bins=bins,
            cbins=np.array([0.5, 1.5, 2.5, 3.5]),
            valid=np.array([True, True, False, True]),
        ),
        "hfunc": "hist",
        "color": "black",
        "label": "Reference",
        "style": dict(hist_style_step),
    }
    fig, axis = plt.subplots()
    try:
        from core.plot import plot

        plot._draw_histogram_record(
            axis,
            record,
            record["hdata"].counts,
            record["hdata"].errs,
            record["color"],
            record["label"],
        )
        step_paths = [patch.get_path().vertices for patch in axis.patches]

        assert len(step_paths) == 2
        assert all(np.all(np.isfinite(vertices)) for vertices in step_paths)
        assert (np.min(step_paths[0][:, 0]), np.max(step_paths[0][:, 0])) == (0.0, 2.0)
        assert (np.min(step_paths[1][:, 0]), np.max(step_paths[1][:, 0])) == (3.0, 4.0)
    finally:
        plt.close(fig)


# Render deterministic plotting examples and require nonempty PDF outputs
def test_plot_visual_basic(tmp_path):
    """Visual unit tests"""

    output_dir = tmp_path / "plot"
    output_dir.mkdir(parents=True)

    # Synthetic input data
    rng = np.random.default_rng(12345)
    r1 = rng.standard_normal(25000) * 0.8
    r2 = rng.standard_normal(25000) * 1
    r3 = rng.standard_normal(25000) * 1.2
    r4 = rng.standard_normal(25000) * 1.5

    # ------------------------------------------------------------------------
    # Mathematical definitions

    # Momentum squared
    def pt2(x):
        return np.power(x, 2)

    # ------------------------------------------------------------------------
    # Observables containers

    obs_pt2 = {
        # Axis limits
        "xlim": (0, 1.5),
        "ylim": None,
        "xlabel": r"$p_t^2$",
        "ylabel": r"Counts",
        "units": {"x": r"GeV$^2$", "y": r"counts"},
        "label": r"Transverse momentum squared",
        # Ratio
        "ylim_ratio": (0.7, 1.3),
        # Histogramming
        "bins": np.linspace(0, 1.5, 60),
        "density": False,
        # Function to calculate
        "func": pt2,
    }

    # ------------------------------------------------------------------------
    # ** Example **

    fig1, ax1 = create_axes(**obs_pt2, ratio_plot=False)
    counts, errs, bins, cbins = hist(
        obs_pt2["func"](r1), bins=obs_pt2["bins"], density=obs_pt2["density"]
    )
    ax1[0].errorbar(
        x=cbins, y=counts, yerr=errs, color=(0, 0, 0), label="Data $\\alpha$", **errorbar_style
    )
    ax1[0].legend(frameon=False)
    fig1.savefig(output_dir / "testplot_1.pdf", bbox_inches="tight")

    # ------------------------------------------------------------------------
    # ** Example **

    fig2, ax2 = create_axes(**obs_pt2, ratio_plot=False)
    counts, errs, bins, cbins = hist(
        obs_pt2["func"](r1), bins=obs_pt2["bins"], density=obs_pt2["density"]
    )
    ax2[0].hist(
        x=cbins,
        bins=bins,
        weights=counts,
        color=(0.5, 0.2, 0.1),
        label="Data $\\alpha$",
        **hist_style_step,
    )
    ax2[0].legend(frameon=False)
    fig2.savefig(output_dir / "testplot_2.pdf", bbox_inches="tight")

    # ------------------------------------------------------------------------
    # ** Example **

    fig3, ax3 = create_axes(**obs_pt2, ratio_plot=True)

    counts1, errs, bins, cbins = hist(
        obs_pt2["func"](r1), bins=obs_pt2["bins"], density=obs_pt2["density"]
    )
    ax3[0].hist(
        x=cbins, bins=bins, weights=counts1, color=(0, 0, 0), label="Data 1", **hist_style_step
    )

    counts2, errs, bins, cbins = hist(
        obs_pt2["func"](r2), bins=obs_pt2["bins"], density=obs_pt2["density"]
    )
    ax3[0].hist(
        x=cbins,
        bins=bins,
        weights=counts2,
        color=(1, 0, 0),
        alpha=0.5,
        label="Data 2",
        **hist_style_step,
    )

    ordered_legend(ax=ax3[0], order=["Data 1", "Data 2"])

    # Ratio
    plot_horizontal_line(ax3[1])
    ax3[1].hist(
        x=cbins,
        bins=bins,
        weights=counts2 / (counts1 + 1e-30),
        color=(1, 0, 0),
        alpha=0.5,
        label="Data $\\beta$",
        **hist_style_step,
    )

    fig3.savefig(output_dir / "testplot_3.pdf", bbox_inches="tight")

    # ------------------------------------------------------------------------
    # ** Example **

    data_template = {
        "data": None,
        "weights": None,
        "label": "Data",
        "hfunc": "errorbar",
        "style": errorbar_style,
        "obs": obs_pt2,
        "hdata": None,
        "color": None,
    }

    # Data source <-> Observable collections
    data1 = data_template.copy()  # Deep copies
    data2 = data_template.copy()
    data3 = data_template.copy()
    data4 = data_template.copy()

    data1.update(
        {
            "data": r1,
            "label": "Data $\\alpha$",
            "hfunc": "errorbar",
            "style": errorbar_style,
        }
    )
    data2.update(
        {
            "data": r2,
            "label": "Data $\\beta$",
            "hfunc": "hist",
            "style": hist_style_step,
        }
    )
    data3.update(
        {
            "data": r3,
            "label": "Data $\\gamma$",
            "hfunc": "hist",
            "style": hist_style_step,
        }
    )
    data4.update(
        {
            "data": r4,
            "label": "Data $\\delta$",
            "hfunc": "plot",
            "style": plot_style,
        }
    )

    data = [data1, data2, data3, data4]

    # Calculate histograms
    for i in range(len(data)):
        data[i]["hdata"] = hist_obj(
            data[i]["obs"]["func"](data[i]["data"]), bins=data[i]["obs"]["bins"]
        )

    # Plot it
    fig4, ax4 = superplot(data, ratio_plot=True, yscale="log")
    fig5, ax5 = superplot(data, ratio_plot=True, yscale="linear", ratio_uncertainty="none")

    fig4.savefig(output_dir / "testplot_4.pdf", bbox_inches="tight")
    fig5.savefig(output_dir / "testplot_5.pdf", bbox_inches="tight")

    stacked, total = build_filled_process_stack_records([data2, data3])
    fig6, _ = superplot(
        stacked,
        ratio_plot=False,
        uncertainty_record=total,
        legend_order=[data2["label"], data3["label"]],
    )
    fig6.savefig(output_dir / "testplot_6.pdf", bbox_inches="tight")

    for index in range(1, 7):
        plot_path = output_dir / f"testplot_{index}.pdf"
        assert plot_path.is_file()
        assert plot_path.stat().st_size > 0
    plt.close("all")


# Check automatic and configured linear histogram limits
def test_linear_ylim_override():
    observable = {
        "xlim": (0.0, 2.0),
        "ylim": None,
        "xlabel": "$x$",
        "ylabel": "Counts",
        "units": {"x": "", "y": ""},
        "ylim_ratio": (0.0, 2.0),
    }
    record = {
        "hdata": hist_obj(np.array([0.5, 1.5]), bins=np.array([0.0, 1.0, 2.0])),
        "hfunc": "hist",
        "color": "black",
        "label": "Test",
        "style": hist_style_step,
        "obs": observable,
    }
    record["hdata"].errs = np.array([1.0, 0.5])

    fig, ax = superplot([record], observable=observable, ratio_plot=False, yscale="linear")
    try:
        assert ax[0].get_ylim()[0] == pytest.approx(0.0)
        assert ax[0].get_ylim()[1] == pytest.approx(3.0)
    finally:
        plt.close(fig)

    observable["ylim"] = (-1.0, 4.0)
    fig, ax = superplot([record], observable=observable, ratio_plot=False, yscale="linear")
    try:
        assert ax[0].get_ylim() == pytest.approx((-1.0, 4.0))
    finally:
        plt.close(fig)


# Check the default ratio errors draw only a hatched reference uncertainty band
@pytest.mark.parametrize("ylim_ratio", [None, (0.5, 1.5)])
def test_ratio_reference_band(ylim_ratio):
    observable = {
        "xlim": (0.0, 2.0),
        "ylim": None,
        "xlabel": "$x$",
        "ylabel": "Counts",
        "units": {"x": "", "y": ""},
        "ylim_ratio": ylim_ratio,
    }
    bins = np.array([0.0, 1.0, 2.0])
    centers = np.array([0.5, 1.5])
    reference = {
        "hdata": hobj(
            counts=np.array([10.0, 20.0]),
            errs=np.array([1.0, 4.0]),
            bins=bins,
            cbins=centers,
        ),
        "hfunc": "errorbar",
        "color": "black",
        "label": "Reference",
        "style": errorbar_style,
        "obs": observable,
    }
    numerator = {
        "hdata": hobj(
            counts=np.array([20.0, 10.0]),
            errs=np.array([2.0, 1.0]),
            bins=bins,
            cbins=centers,
        ),
        "hfunc": "errorbar",
        "color": "red",
        "label": "Numerator",
        "style": errorbar_style,
        "obs": observable,
    }

    fig, ax = superplot(
        [reference, numerator],
        observable=observable,
        ratio_plot=True,
    )
    try:
        if ylim_ratio is None:
            assert ax[1].get_autoscaley_on()
            assert ax[1].get_ylim()[0] < 0.5 and ax[1].get_ylim()[1] > 2.0
        else:
            assert ax[1].get_ylim() == pytest.approx(ylim_ratio)
        bands = [item for item in ax[1].collections if isinstance(item, PolyCollection)]
        errorbars = [item for item in ax[1].containers if isinstance(item, ErrorbarContainer)]

        assert len(bands) == 1
        assert len(errorbars) == 1
        assert bands[0].get_hatch() == data_total_style["hatch"]
        vertices = np.concatenate([path.vertices for path in bands[0].get_paths()])
        assert np.min(vertices[:, 0]) == pytest.approx(0.0)
        assert np.max(vertices[:, 0]) == pytest.approx(2.0)
        assert np.min(vertices[:, 1]) == pytest.approx(0.8)
        assert np.max(vertices[:, 1]) == pytest.approx(1.2)
    finally:
        plt.close(fig)


# Check point data use statistical bars and full-width hatched total errors
def test_point_stat_total_errors():
    observable = {
        "xlim": (0.0, 2.0),
        "ylim": None,
        "xlabel": "$x$",
        "ylabel": "Counts",
        "units": {"x": "", "y": ""},
        "ylim_ratio": (0.5, 1.5),
    }
    bins = np.array([0.0, 1.0, 2.0])
    record = {
        "hdata": hobj(
            counts=np.array([10.0, 20.0]),
            errs=np.array([3.0, 4.0]),
            bins=bins,
            cbins=np.array([0.5, 1.5]),
        ),
        "hfunc": "errorbar",
        "color": "black",
        "label": "Data",
        "style": dict(errorbar_style),
        "stat_errs": np.array([1.0, 2.0]),
        "obs": observable,
    }

    fig, ax = superplot([record], observable=observable, ratio_plot=False)
    try:
        bands = [item for item in ax[0].collections if isinstance(item, PolyCollection)]
        errorbars = [item for item in ax[0].containers if isinstance(item, ErrorbarContainer)]

        assert len(bands) == 1
        assert len(errorbars) == 1
        assert bands[0].get_hatch() == data_total_style["hatch"]
        vertices = np.concatenate([path.vertices for path in bands[0].get_paths()])
        assert np.min(vertices[:, 0]) == pytest.approx(0.0)
        assert np.max(vertices[:, 0]) == pytest.approx(2.0)
        assert np.min(vertices[:, 1]) == pytest.approx(7.0)
        assert np.max(vertices[:, 1]) == pytest.approx(24.0)
        segments = errorbars[0].lines[2][0].get_segments()
        np.testing.assert_allclose(segments[0][:, 1], [9.0, 11.0])
        np.testing.assert_allclose(segments[1][:, 1], [18.0, 22.0])
    finally:
        plt.close(fig)


# Check step data use uniform statistical shading and hatched total errors
def test_step_stat_total_errors():
    observable = {
        "xlim": (0.0, 2.0),
        "ylim": None,
        "xlabel": "$x$",
        "ylabel": "Counts",
        "units": {"x": "", "y": ""},
        "ylim_ratio": (0.5, 1.5),
    }
    record = {
        "hdata": hobj(
            counts=np.array([10.0, 20.0]),
            errs=np.array([3.0, 4.0]),
            bins=np.array([0.0, 1.0, 2.0]),
            cbins=np.array([0.5, 1.5]),
        ),
        "hfunc": "hist",
        "color": "black",
        "label": "Data",
        "style": dict(hist_style_step),
        "stat_errs": np.array([1.0, 2.0]),
        "obs": observable,
    }

    fig, ax = superplot([record], observable=observable, ratio_plot=False)
    try:
        bands = [item for item in ax[0].collections if isinstance(item, PolyCollection)]

        assert len(ax[0].patches) == 1
        assert len(bands) == 2
        total = next(item for item in bands if item.get_hatch())
        statistical = next(item for item in bands if not item.get_hatch())
        assert total.get_hatch() == data_total_style["hatch"]
        total_vertices = np.concatenate([path.vertices for path in total.get_paths()])
        stat_vertices = np.concatenate([path.vertices for path in statistical.get_paths()])
        assert np.min(total_vertices[:, 1]) == pytest.approx(7.0)
        assert np.max(total_vertices[:, 1]) == pytest.approx(24.0)
        assert np.min(stat_vertices[:, 1]) == pytest.approx(9.0)
        assert np.max(stat_vertices[:, 1]) == pytest.approx(22.0)
    finally:
        plt.close(fig)


# Check ratio step curves stay above their uncertainty bands and unity guide
def test_superplot_ratio_hist_lines_foreground():
    observable = {
        "xlim": (0.0, 3.0),
        "ylim": None,
        "xlabel": "$x$",
        "ylabel": "Counts",
        "units": {"x": "", "y": ""},
        "ylim_ratio": (0.0, 2.0),
    }
    bins = np.arange(4.0)
    centers = np.array([0.5, 1.5, 2.5])

    # Build one step histogram ratio input
    def histogram_record(counts, errors, color, label):
        return {
            "hdata": hobj(
                counts=np.asarray(counts),
                errs=np.asarray(errors),
                bins=bins,
                cbins=centers,
            ),
            "hfunc": "hist",
            "color": color,
            "label": label,
            "style": dict(hist_style_step),
            "obs": observable,
        }

    records = [
        histogram_record([2.0, 4.0, 2.0], [0.2, 0.4, 0.2], "red", "Reference"),
        histogram_record([1.0, 4.0, 3.0], [0.1, 0.4, 0.3], "blue", "Numerator"),
    ]
    fig, ax = superplot(records, observable=observable, ratio_plot=True)
    try:
        ratio_axis = ax[1]
        step_depths = [item.get_zorder() for item in ratio_axis.patches]
        band_depths = [item.get_zorder() for item in ratio_axis.collections]
        guide_depths = [item.get_zorder() for item in ratio_axis.lines]

        assert len(step_depths) == 1
        assert len(band_depths) == 2
        assert len(guide_depths) == 1
        assert max(guide_depths) < min(band_depths)
        assert max(band_depths) < min(step_depths)
    finally:
        plt.close(fig)

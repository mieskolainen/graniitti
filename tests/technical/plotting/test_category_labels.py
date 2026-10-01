# Unit tests for categorical iceplot axis labels
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numpy as np
import pytest
from core.plot import plot
from core.stats import cov
from matplotlib import pyplot as plt

EXPECTED_SUFFIX = "p1T $\\in$ [0.00, 0.80] [GeV], p2T $\\in$ [0.00, 0.80] [GeV]"


# Compute a minimal 3D category observable for label formatting tests
def category_obs(x0_bins=None, x1_bins=None, symmetrized_fill=False):
    if x0_bins is None:
        x0_bins = np.array([0.0, 0.8])
    if x1_bins is None:
        x1_bins = np.array([0.0, 0.8])

    bins = {
        "cat0": {
            0: x0_bins,
            1: x1_bins,
            2: np.array([0.0, 1.0, 2.0]),
        }
    }

    return {
        "toy": {
            "tag": "toy",
            "xlim": {"cat0": (0.0, 2.0)},
            "ylim": None,
            "xlabel": "max(t, u)",
            "ylabel": "d sigma",
            "units": {"x": "GeV", "y": "pb"},
            "label": "Toy",
            "ylim_ratio": (0.0, 2.0),
            "ytick_ratio_step": 0.5,
            "bins": bins,
            "density": False,
            "symmetrized_fill": symmetrized_fill,
            "category_labels": ("p1T", "p2T"),
            "category_units": ("GeV", "GeV"),
            "func": None,
        }
    }


# Compute a minimal MC histogram input with one entry inside the category
def category_mcdata():
    return {
        "data": {"toy": np.array([[0.4, 0.4, 0.5]])},
        "weights": np.array([1.0]),
        "xsection_pb": 1.0,
    }


# Compute a minimal HEPData-like histogram input inside one category
def category_hepdata():
    return {
        "toy": {
            "y": {"cat0": np.array([1.0, 2.0])},
            "y_err": {"cat0": np.array([0.1, 0.2])},
            "bins": {
                "cat0": {
                    0: np.array([0.0, 0.8]),
                    1: np.array([0.0, 0.8]),
                    2: np.array([0.0, 1.0, 2.0]),
                }
            },
            "binwidth": {"cat0": {2: np.array([1.0, 1.0])}},
            "x": {"cat0": {2: np.array([0.5, 1.5])}},
            "xlim": {"cat0": (0.0, 2.0)},
            "scale": 1.0,
            "fitw": np.array([1.0, 1.0]),
        }
    }


def test_create_axes_reduces_long_xlabel_fontsize():
    """
    Check that long category labels are automatically shrunk
    """
    long_suffix = ", ".join([EXPECTED_SUFFIX] * 3)
    fig, ax = plot.create_axes(
        xlabel="max(t, u)",
        ylabel="d sigma",
        units={"x": "GeV", "y": "pb"},
        xlabel_suffix=long_suffix,
        ratio_plot=False,
        fontsize=9,
    )

    try:
        assert ax[-1].xaxis.label.get_fontsize() < 9
    finally:
        plt.close(fig)


# Shrink an inside legend when it covers a rendered curve
def test_draw_superplot_legend_shrinks_data_overlap():
    fig, axis = plt.subplots()
    axis.plot([0.0, 1.0], [0.95, 0.95], label="Overlapped sample")
    axis.set_ylim(0.0, 1.0)

    try:
        plot.draw_superplot_legend(
            fig=fig,
            ax=[axis],
            legend_labels=["Overlapped sample"],
            legend_properties={"fontsize": 8.0, "loc": "upper right"},
        )
        legend = axis.get_legend()
        assert legend is not None
        assert legend.get_texts()[0].get_fontsize() < 8.0
    finally:
        plt.close(fig)


# Check HEPData preserves automatic and explicitly configured y limits
@pytest.mark.parametrize("ylim", [None, (-10.0, 10.0)])
def test_histhepdata_config_ylim(ylim):
    obs = category_obs()
    obs["toy"]["ylim"] = ylim
    hepdata = category_hepdata()
    hepdata["toy"]["y"]["cat0"] = np.array([-1.0, 2.0])
    hepdata["toy"]["y_err"]["cat0"] = np.array([0.5, 1.0])

    out = plot.histhepdata(hepdata=hepdata, obs=obs)
    obs_name = "toy_[0.00-0.80]_[0.00-0.80]"

    assert out[obs_name]["obs"]["ylim"] == ylim

    if ylim is None:
        reference = out[obs_name]
        mc = copy.deepcopy(reference)
        mc["hdata"].counts *= 1000.0
        mc["hdata"].errs *= 1000.0
        for records in ([reference, mc], [mc, reference]):
            for yscale in ("linear", "log"):
                fig, axes = plot.superplot(records, ratio_plot=False, yscale=yscale)
                try:
                    lower, upper = axes[0].get_ylim()
                    assert upper > np.max(mc["hdata"].counts_scaled + mc["hdata"].errs_scaled)
                    assert lower < reference["hdata"].counts_scaled[1] - reference["hdata"].errs_scaled[1]
                finally:
                    plt.close(fig)


def test_missing_category_labels():
    """
    Check that dict-binned observables must define category label metadata
    """
    obs = category_obs()
    del obs["toy"]["category_labels"]

    with pytest.raises(KeyError, match="category_labels"):
        plot.histmc(mcdata=category_mcdata(), obs=obs)


def test_histmc_malformed_category_label_metadata():
    """
    Check that category label metadata has exactly two dimensions
    """
    obs = copy.deepcopy(category_obs())
    obs["toy"]["category_units"] = ("GeV",)

    with pytest.raises(ValueError, match="exactly two"):
        plot.histmc(mcdata=category_mcdata(), obs=obs)


def test_symmetrized_category_area():
    """
    Check that differential normalization uses displayed category widths
    """
    assert plot.category_plane_area(
        x0_bins=np.array([0.0, 1.0]),
        x1_bins=np.array([0.0, 1.0]),
        symmetrized=True,
        obs_name="toy",
    ) == pytest.approx(1.0)

    assert plot.category_plane_area(
        x0_bins=np.array([0.0, 1.0]),
        x1_bins=np.array([2.0, 5.0]),
        symmetrized=True,
        obs_name="toy",
    ) == pytest.approx(3.0)

    assert plot.category_plane_area(
        x0_bins=np.array([0.0, 2.0]),
        x1_bins=np.array([1.0, 3.0]),
        symmetrized=True,
        obs_name="toy",
    ) == pytest.approx(4.0)

    assert plot.category_plane_area(
        x0_bins=np.array([0.0, 2.0]),
        x1_bins=np.array([1.0, 3.0]),
        symmetrized=False,
        obs_name="toy",
    ) == pytest.approx(4.0)


# Check proton exchange, diagonal counting, weighted errors and covariance assignments
@pytest.mark.parametrize(
    "x1_bins, points, area",
    [
        ([2.0, 5.0], [[0.5, 2.5, 0.5], [2.5, 0.5, 0.5], [1.0, 5.0, 1.5]], 3.0),
        ([0.0, 1.0], [[0.25, 0.75, 0.5], [0.75, 0.25, 0.5], [1.0, 1.0, 1.5]], 1.0),
    ],
)
def test_symmetrized_binwidth(x1_bins, points, area):
    obs = category_obs(
        x0_bins=np.array([0.0, 1.0]),
        x1_bins=np.array(x1_bins),
        symmetrized_fill=True,
    )
    mcdata = {
        "data": {"toy": np.array(points)},
        "weights": np.array([1.0, 2.0, 3.0]),
        "event_ids": np.arange(3),
        "xsection_pb": 1.0,
    }

    out = plot.histmc(mcdata=mcdata, obs=obs)
    name = next(iter(out))
    hdata = out[name]["hdata"]

    assert hdata.counts == pytest.approx([3.0, 3.0])
    assert hdata.errs == pytest.approx([np.sqrt(5.0), 3.0])
    assert hdata.binscale == pytest.approx(np.full(2, 1.0 / (6.0 * area)))
    assert hdata.integral() * area == pytest.approx(mcdata["xsection_pb"])

    event_ids, bin_ids = cov.subset_assignments(mcdata, obs)[name]
    np.testing.assert_array_equal(event_ids, mcdata["event_ids"])
    np.testing.assert_array_equal(bin_ids, [0, 0, 1])

    mcdata["data"]["toy"] = mcdata["data"]["toy"][:, [1, 0, 2]]
    swapped = plot.histmc(mcdata=mcdata, obs=obs)[name]["hdata"]
    assert swapped.counts == pytest.approx(hdata.counts)
    assert swapped.errs == pytest.approx(hdata.errs)
    np.testing.assert_array_equal(cov.subset_assignments(mcdata, obs)[name][1], bin_ids)

# Unit tests for plot histogram fusion
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np
import pytest
from core.plot import plot
from core.stats import hist
from core.stats.uncertainty import finalize_uncertainties, source_covariance, uncertainty_source


# Compute a one-bin histogram with configurable counts and scale
def one_bin_hobj(counts, errs, binscale, bins=None):
    if bins is None:
        bins = np.array([0.0, 1.0])
    cbins = np.array([(bins[0] + bins[1]) / 2.0])

    return hist.hobj(
        counts=np.array([counts], dtype=float),
        errs=np.array([errs], dtype=float),
        bins=bins,
        cbins=cbins,
        binscale=np.array([binscale], dtype=float),
    )


# Check that weighted histogram errors retain very small representable weights
def test_hist_error_avoids_weight_square_underflow():
    counts, errors, _, _ = hist.hist(
        x=np.array([0.25, 0.75]),
        bins=np.array([0.0, 1.0]),
        weights=np.array([1.0e-200, 1.0e-200]),
    )

    assert counts[0] == pytest.approx(2.0e-200, abs=0.0)
    assert errors[0] == pytest.approx(np.sqrt(2.0) * 1.0e-200, abs=0.0)
    assert errors[0] > 0.0


# Check that weighted histogram errors avoid intermediate square overflow
def test_hist_error_avoids_weight_square_overflow():
    _, errors, _, _ = hist.hist(
        x=np.array([0.25, 0.75]),
        bins=np.array([0.0, 1.0]),
        weights=np.array([1.0e200, 1.0e200]),
    )

    assert np.isfinite(errors[0])
    assert errors[0] == pytest.approx(np.sqrt(2.0) * 1.0e200)


# Check that histogram rebinning uses the same stable error norm
def test_rebin_error_avoids_square_underflow():
    _, _, errors = hist.rebin_histogram(
        bin_edges=np.array([0.0, 1.0, 2.0]),
        bin_contents=np.array([1.0e-200, 1.0e-200]),
        bin_errors=np.array([1.0e-200, 1.0e-200]),
        rebin_factor=2,
        differential=False,
    )

    assert errors[0] == pytest.approx(np.sqrt(2.0) * 1.0e-200, abs=0.0)


def test_zero_count_harmonic_scale():
    """
    Check that empty selected chunks still contribute normalization exposure
    """
    empty = one_bin_hobj(counts=0.0, errs=0.0, binscale=0.01)
    filled = one_bin_hobj(counts=10.0, errs=3.0, binscale=0.01)

    left = empty.fuse_independent_chunk(filled)
    right = filled.fuse_independent_chunk(empty)

    np.testing.assert_allclose(left.counts, [10.0])
    np.testing.assert_allclose(right.counts, [10.0])
    np.testing.assert_allclose(left.errs, [3.0])
    np.testing.assert_allclose(right.errs, [3.0])
    np.testing.assert_allclose(left.binscale, [0.005])
    np.testing.assert_allclose(right.binscale, [0.005])


def test_zero_count_fusion_order():
    """
    Check that in-place fusion keeps empty chunk exposure independent of order
    """
    left = hist.hobj.fuse_independent_chunks(
        [
            one_bin_hobj(counts=0.0, errs=0.0, binscale=0.01),
            one_bin_hobj(counts=10.0, errs=3.0, binscale=0.01),
        ]
    )
    right = hist.hobj.fuse_independent_chunks(
        [
            one_bin_hobj(counts=10.0, errs=3.0, binscale=0.01),
            one_bin_hobj(counts=0.0, errs=0.0, binscale=0.01),
        ]
    )

    np.testing.assert_allclose(left.counts, right.counts)
    np.testing.assert_allclose(left.errs, right.errs)
    np.testing.assert_allclose(left.binscale, right.binscale)
    np.testing.assert_allclose(left.binscale, [0.005])


def test_zero_scale_hist_fusion_identity():
    """
    Check that histograms without finite positive scale are true identities
    """
    empty = one_bin_hobj(counts=0.0, errs=0.0, binscale=0.0)
    filled = one_bin_hobj(counts=10.0, errs=3.0, binscale=0.01)

    left = empty.fuse_independent_chunk(filled)
    right = filled.fuse_independent_chunk(empty)

    np.testing.assert_allclose(left.counts, [10.0])
    np.testing.assert_allclose(right.counts, [10.0])
    np.testing.assert_allclose(left.binscale, [0.01])
    np.testing.assert_allclose(right.binscale, [0.01])


def test_hist_fusion_rejects_different_bin_edges():
    """
    Check that incompatible histograms raise a controlled error
    """
    left = one_bin_hobj(counts=1.0, errs=1.0, binscale=0.01)
    right = one_bin_hobj(
        counts=1.0,
        errs=1.0,
        binscale=0.01,
        bins=np.array([0.0, 2.0]),
    )

    with pytest.raises(ValueError, match="different histogram bins"):
        left.fuse_independent_chunk(right)


# Build one density histogram through the shared iceplot and icetune MC path
def density_histogram(values, weights, density_uncertainty="scaled"):
    observable = {
        "x": {
            "func": None,
            "bins": np.array([0.0, 1.0, 2.0]),
            "xlim": (0.0, 2.0),
            "xlabel": "$x$",
            "units": {"x": "1", "y": "pb"},
        }
    }
    return plot.histmc(
        mcdata={
            "data": {"x": np.asarray(values, dtype=float)},
            "weights": np.asarray(weights, dtype=float),
            "xsection_pb": 1.0,
        },
        obs=observable,
        density=True,
        density_uncertainty=density_uncertainty,
    )


# Check unequal density chunks reproduce one serial histogram exactly
@pytest.mark.parametrize("density_uncertainty", ["scaled", "shape"])
def test_density_unequal_chunks(density_uncertainty):
    first = density_histogram(
        values=[0.2, 0.4],
        weights=[1.0, 1.0],
        density_uncertainty=density_uncertainty,
    )
    second = density_histogram(
        values=[1.1, 1.2, 1.3, 1.4],
        weights=[1.0, 2.0, 3.0, 4.0],
        density_uncertainty=density_uncertainty,
    )
    serial = density_histogram(
        values=[0.2, 0.4, 1.1, 1.2, 1.3, 1.4],
        weights=[1.0, 1.0, 1.0, 2.0, 3.0, 4.0],
        density_uncertainty=density_uncertainty,
    )

    fused = plot.fuse_worker_chunk_outputs([[first], [second]])[0]["x"]["hdata"]
    expected = serial["x"]["hdata"]

    assert fused.density
    np.testing.assert_allclose(fused.counts, expected.counts)
    np.testing.assert_allclose(fused.errs, expected.errs)
    np.testing.assert_allclose(fused.counts_scaled, expected.counts_scaled)
    np.testing.assert_allclose(fused.errs_scaled, expected.errs_scaled)
    np.testing.assert_allclose(fused.covariance_scaled, expected.covariance_scaled)
    assert fused.integral() == pytest.approx(1.0)


# Check one-bin shape normalization removes its only counting mode
def test_bin_shape_density_zero_variance():
    histogram = hist.hobj(
        counts=np.asarray([4.0]),
        errs=np.asarray([2.0]),
        bins=np.asarray([0.0, 2.0]),
        cbins=np.asarray([1.0]),
        density=True,
        density_uncertainty="shape",
    )

    np.testing.assert_allclose(histogram.counts_scaled, [0.5])
    np.testing.assert_allclose(histogram.errs_scaled_fixed, [0.25])
    np.testing.assert_allclose(histogram.covariance_scaled, [[0.0]], atol=1.0e-15)
    np.testing.assert_allclose(histogram.errs_scaled, [0.0], atol=1.0e-15)
    assert histogram.integral_error() == pytest.approx(0.0, abs=1.0e-15)


# Check the density covariance has the expected null integral mode
def test_two_bin_shape_density_cov_rank_n_minus():
    histogram = hist.hobj(
        counts=np.asarray([3.0, 7.0]),
        errs=np.sqrt(np.asarray([3.0, 7.0])),
        bins=np.asarray([0.0, 1.0, 3.0]),
        cbins=np.asarray([0.5, 2.0]),
        density=True,
        density_uncertainty="shape",
    )

    covariance = histogram.covariance_scaled
    binwidth = histogram.binwidth
    np.testing.assert_allclose(covariance @ binwidth, np.zeros(2), atol=1.0e-15)
    assert np.linalg.matrix_rank(covariance, tol=1.0e-14) == 1
    assert histogram.integral_error() == pytest.approx(0.0, abs=1.0e-15)


# Check weighted MC density covariance against the explicit count Jacobian
def test_shape_density_cov_matches_count_jacobian():
    counts = np.asarray([2.0, 5.0, 7.0])
    errors = np.asarray([0.4, 1.5, 2.3])
    histogram = hist.hobj(
        counts=counts,
        errs=errors,
        bins=np.asarray([0.0, 1.0e8, 3.0e8, 6.0e8]),
        cbins=np.asarray([0.5e8, 2.0e8, 4.5e8]),
        density=True,
        density_uncertainty="shape",
    )

    total = np.sum(counts)
    binwidth = histogram.binwidth
    jacobian = np.diag(1.0 / (total * binwidth)) - np.outer(
        counts / (total**2 * binwidth),
        np.ones(len(counts)),
    )
    expected = jacobian @ np.diag(errors**2) @ jacobian.T

    assert np.any(expected != 0.0)
    np.testing.assert_allclose(histogram.covariance_scaled, expected, rtol=1.0e-12, atol=0.0)
    np.testing.assert_allclose(histogram.errs_scaled, np.sqrt(np.diag(expected)))


# Check density and cross-section histogram objects cannot be fused accidentally
def test_hist_fusion_rejects_mixed_norm():
    density_hist = one_bin_hobj(counts=1.0, errs=1.0, binscale=1.0)
    density_hist.density = True
    cross_section_hist = one_bin_hobj(counts=1.0, errs=1.0, binscale=1.0)

    with pytest.raises(ValueError, match="cannot mix density"):
        density_hist.fuse_independent_chunk(cross_section_hist)


# Check density chunks cannot silently mix uncertainty conventions
def test_hist_fusion_rejects_mixed_density_error():
    scaled = density_histogram(values=[0.2], weights=[1.0], density_uncertainty="scaled")["x"][
        "hdata"
    ]
    shape = density_histogram(values=[0.2], weights=[1.0], density_uncertainty="shape")["x"][
        "hdata"
    ]

    with pytest.raises(ValueError, match="matching density uncertainty"):
        scaled.fuse_independent_chunk(shape)


# Check density normalization excludes bins omitted from the comparison interval
def test_density_hist_normalizes_over_valid_bins():
    histogram = hist.hobj(
        counts=np.asarray([2.0, 8.0]),
        errs=np.asarray([1.0, 2.0]),
        bins=np.asarray([0.0, 1.0, 2.0]),
        cbins=np.asarray([0.5, 1.5]),
        binscale=1.0,
        valid=np.asarray([True, False]),
        density=True,
        density_uncertainty="shape",
    )

    np.testing.assert_allclose(histogram.counts_scaled, [1.0, 0.0])
    np.testing.assert_allclose(histogram.covariance_scaled, np.zeros((2, 2)), atol=1.0e-15)
    assert histogram.integral() == pytest.approx(1.0)


# Check independent process summation adds scaled cross sections and variances
def test_process_sum_vs_chunk_fusion():
    first = one_bin_hobj(counts=10.0, errs=2.0, binscale=2.0)
    second = one_bin_hobj(counts=3.0, errs=1.0, binscale=5.0)

    process_sum = hist.hobj.sum_independent_processes([first, second])
    chunk_fusion = hist.hobj.fuse_independent_chunks([first, second])

    np.testing.assert_allclose(process_sum.counts_scaled, [35.0])
    np.testing.assert_allclose(process_sum.errs_scaled, [np.hypot(4.0, 5.0)])
    assert process_sum.binscale == pytest.approx(1.0)
    assert chunk_fusion.counts_scaled[0] != pytest.approx(process_sum.counts_scaled[0])


# Build one ordinary MC record for process-stack tests
def process_record(label, histogram, color):
    return {
        "hdata": histogram,
        "hfunc": "hist",
        "label": label,
        "color": color,
        "style": dict(plot.hist_style_step),
        "obs": {
            "xlabel": "$x$",
            "ylabel": "$d\\sigma/dx$",
            "units": {"x": "1", "y": "pb"},
        },
    }


# Check filled stack records are cumulative while the total uses the process sum
def test_process_stack_sum():
    records = [
        process_record("A", one_bin_hobj(10.0, 2.0, 2.0), "red"),
        process_record("B", one_bin_hobj(3.0, 1.0, 5.0), "blue"),
    ]

    visual, total = plot.build_filled_process_stack_records(records)

    assert [record["label"] for record in visual] == ["B", "A"]
    np.testing.assert_allclose(visual[0]["hdata"].counts_scaled, [35.0])
    np.testing.assert_allclose(visual[1]["hdata"].counts_scaled, [20.0])
    np.testing.assert_allclose(total["hdata"].counts_scaled, [35.0])
    np.testing.assert_allclose(total["hdata"].errs_scaled, [np.hypot(4.0, 5.0)])
    assert total["label"] == "MC total"


# Check density stacking keeps cross-section fractions and normalizes only the total
@pytest.mark.parametrize("density_uncertainty", ["scaled", "shape"])
def test_density_process_stack_shared_final_norm(density_uncertainty):
    bins = np.array([0.0, 1.0, 3.0])
    centers = np.array([0.5, 2.0])
    first = hist.hobj(
        counts=np.array([2.0, 3.0]),
        errs=np.array([0.2, 0.3]),
        bins=bins,
        cbins=centers,
        binscale=np.array([4.0, 2.0]),
    )
    second = hist.hobj(
        counts=np.array([1.0, 2.0]),
        errs=np.array([0.1, 0.2]),
        bins=bins,
        cbins=centers,
        binscale=np.array([2.0, 5.0]),
    )
    records = [
        process_record("A", first, "red"),
        process_record("B", second, "blue"),
    ]

    visual, total = plot.build_filled_process_stack_records(
        records,
        density=True,
        density_uncertainty=density_uncertainty,
    )

    total_cross_section = 42.0
    np.testing.assert_allclose(total["hdata"].counts_scaled, [10.0 / 42.0, 16.0 / 42.0])
    np.testing.assert_allclose(visual[1]["hdata"].counts_scaled, [8.0 / 42.0, 6.0 / 42.0])
    direct_total = hist.hobj(
        counts=np.asarray([10.0, 32.0]),
        errs=np.asarray([np.hypot(0.8, 0.2), np.hypot(1.2, 2.0)]),
        bins=bins,
        cbins=centers,
        density=True,
        density_uncertainty=density_uncertainty,
    )
    np.testing.assert_allclose(total["hdata"].covariance_scaled, direct_total.covariance_scaled)
    assert total["hdata"].integral() == pytest.approx(1.0)
    assert visual[1]["hdata"].integral() == pytest.approx(20.0 / total_cross_section)
    if density_uncertainty == "shape":
        np.testing.assert_allclose(
            total["hdata"].covariance_scaled @ np.diff(bins),
            np.zeros(2),
            atol=1.0e-15,
        )
    assert total["obs"]["units"]["y"] == "1"
    assert records[0]["obs"]["units"]["y"] == "pb"


# Check every cumulative stack component uses the final common interval mask
def test_density_process_stack_final_valid_intervals():
    bins = np.asarray([0.0, 1.0, 2.0])
    centers = np.asarray([0.5, 1.5])
    first = hist.hobj(
        counts=np.asarray([2.0, 8.0]),
        errs=np.asarray([0.2, 0.8]),
        bins=bins,
        cbins=centers,
        valid=np.asarray([True, True]),
    )
    second = hist.hobj(
        counts=np.asarray([3.0, 4.0]),
        errs=np.asarray([0.3, 0.4]),
        bins=bins,
        cbins=centers,
        valid=np.asarray([True, False]),
    )

    visual, total = plot.build_filled_process_stack_records(
        [
            process_record("A", first, "red"),
            process_record("B", second, "blue"),
        ],
        density=True,
        density_uncertainty="shape",
    )

    lower = visual[1]["hdata"]
    np.testing.assert_array_equal(lower.valid, [True, False])
    np.testing.assert_allclose(lower.counts_scaled, [2.0 / 5.0, 0.0])
    np.testing.assert_allclose(total["hdata"].counts_scaled, [1.0, 0.0])
    assert total["hdata"].integral() == pytest.approx(1.0)


# Check foreground records are copied and drawn above every reference
def test_foreground_hist_records_raise_data_zorder():
    data = process_record("Data", one_bin_hobj(12.0, 2.0, 1.0), "black")
    mc = process_record("MC", one_bin_hobj(10.0, 1.0, 1.0), "red")
    foreground = plot.foreground_histogram_records([data], [mc])

    assert foreground[0]["style"]["zorder"] == 1

    points = data.copy()
    points["style"] = dict(plot.errorbar_style)
    foreground = plot.foreground_histogram_records([points], [mc])

    assert foreground[0]["style"]["zorder"] == plot.errorbar_style["zorder"]

    mc["style"]["zorder"] = 4
    foreground = plot.foreground_histogram_records([data], [mc])

    assert foreground[0] is not data
    assert foreground[0]["style"]["zorder"] > mc["style"]["zorder"]
    assert data["style"]["zorder"] == plot.hist_style_step["zorder"]


# Check process stacking rejects negative bins and normalized densities
def test_process_stacking_unsupported_hists():
    negative = process_record("negative", one_bin_hobj(-1.0, 1.0, 1.0), "red")
    density = process_record("density", one_bin_hobj(1.0, 1.0, 1.0), "blue")
    density["hdata"].density = True

    with pytest.raises(ValueError, match="nonnegative bins"):
        plot.build_filled_process_stack_records([negative])
    with pytest.raises(ValueError, match="does not accept density"):
        plot.build_filled_process_stack_records([density], density=True)


# Check explicit histogram methods reject empty and malformed inputs
def test_hist_combination_methods_validate_inputs():
    with pytest.raises(ValueError, match="at least one"):
        hist.hobj.fuse_independent_chunks([])
    with pytest.raises(ValueError, match="at least one"):
        hist.hobj.sum_independent_processes([])
    with pytest.raises(ValueError, match="at least one"):
        hist.hobj.stack_independent_processes([])

    valid = one_bin_hobj(1.0, 1.0, 1.0)
    negative_error = one_bin_hobj(1.0, -1.0, 1.0)
    nonfinite_scale = one_bin_hobj(1.0, 1.0, np.nan)
    with pytest.raises(ValueError, match="nonnegative errors"):
        valid.sum_independent_process(negative_error)
    with pytest.raises(ValueError, match="finite histogram values"):
        valid.fuse_independent_chunk(nonfinite_scale)


# Build one source-aware two-bin HEPData density input
def density_hepdata_input():
    values = np.asarray([2.0, 4.0])
    dataset = {
        "x": np.asarray([0.5, 1.5]),
        "bins": np.asarray([0.0, 1.0, 2.0]),
        "binwidth": np.ones(2),
        "xlim": np.asarray([0.0, 2.0]),
        "y": values,
        "scale": 1.0,
        "fitw": 1.0,
    }
    finalize_uncertainties(
        dataset,
        [
            uncertainty_source(
                "normalization",
                0.1 * values,
                category="systematic",
                correlation="collective",
                scope="test:normalization",
                effect="multiplicative",
            )
        ],
    )
    observable = {
        "x": {
            "xlabel": "$x$",
            "ylabel": "$dN/dx$",
            "units": {"x": "1", "y": "1"},
        }
    }
    return {"x": dataset}, observable


# Check density display errors follow the selected fixed-integral or shape convention
def test_density_error_marginals():
    hepdata, observable = density_hepdata_input()
    scaled = plot.histhepdata(
        hepdata=hepdata,
        obs=observable,
        density=True,
        density_uncertainty="scaled",
        MC_XS_SCALE=1.0,
    )["x"]
    shape = plot.histhepdata(
        hepdata=hepdata,
        obs=observable,
        density=True,
        density_uncertainty="shape",
        MC_XS_SCALE=1.0,
    )["x"]

    np.testing.assert_allclose(scaled["hdata"].counts_scaled, [1.0 / 3.0, 2.0 / 3.0])
    np.testing.assert_allclose(scaled["hdata"].errs_scaled, [1.0 / 30.0, 1.0 / 15.0])
    np.testing.assert_allclose(shape["hdata"].errs_scaled, [0.0, 0.0], atol=1.0e-15)
    for record in (scaled, shape):
        covariance = sum(
            (source_covariance(source) for source in record["uncertainties"]),
            start=np.zeros((2, 2)),
        )
        np.testing.assert_allclose(record["hdata"].errs_scaled, np.sqrt(np.diag(covariance)))


# Check HEPData records retain the transformed pure statistical uncertainty
def test_hepdata_statistical_error_plotting():
    values = np.asarray([2.0, 4.0])
    dataset = {
        "x": np.asarray([0.5, 1.5]),
        "bins": np.asarray([0.0, 1.0, 2.0]),
        "binwidth": np.ones(2),
        "xlim": np.asarray([0.0, 2.0]),
        "y": values,
        "scale": 1.0,
        "fitw": 1.0,
    }
    finalize_uncertainties(
        dataset,
        [
            uncertainty_source(
                "stat",
                np.asarray([0.2, 0.4]),
                category="statistical",
                correlation="uncorrelated",
            ),
            uncertainty_source(
                "syst",
                np.asarray([0.3, 0.1]),
                category="systematic",
                correlation="uncorrelated",
            ),
        ],
    )
    observable = {
        "x": {
            "xlabel": "$x$",
            "ylabel": "$dN/dx$",
            "units": {"x": "1", "y": "1"},
        }
    }

    record = plot.histhepdata(
        hepdata={"x": dataset},
        obs=observable,
        scale=2.0,
        MC_XS_SCALE=1.0,
    )["x"]

    np.testing.assert_allclose(record["stat_errs"], [[0.4, 0.8], [0.4, 0.8]])
    np.testing.assert_allclose(
        record["hdata"].errs_scaled,
        2.0 * np.sqrt(np.asarray([0.2, 0.4]) ** 2 + np.asarray([0.3, 0.1]) ** 2),
    )


# Check total errors define a covariance source before density normalization
def test_density_missing_error_source():
    dataset = {
        "x": np.asarray([0.5, 1.5]),
        "bins": np.asarray([0.0, 1.0, 2.0]),
        "binwidth": np.ones(2),
        "xlim": np.asarray([0.0, 2.0]),
        "y": np.asarray([2.0, 4.0]),
        "y_err": np.asarray([0.2, 0.4]),
        "scale": 1.0,
        "fitw": 1.0,
    }
    observable = {
        "x": {
            "xlabel": "$x$",
            "ylabel": "$dN/dx$",
            "units": {"x": "1", "y": "1"},
        }
    }

    record = plot.histhepdata(
        hepdata={"x": dataset},
        obs=observable,
        density=True,
        density_uncertainty="shape",
        MC_XS_SCALE=1.0,
    )["x"]
    covariance = source_covariance(record["uncertainties"][0])

    assert record["uncertainties"][0]["correlation"] == "covariance"
    np.testing.assert_allclose(covariance @ np.ones(2), np.zeros(2), atol=1.0e-15)
    np.testing.assert_allclose(record["hdata"].errs_scaled, np.sqrt(np.diag(covariance)))


# Check omitted data intervals do not enter density values or uncertainties
def test_density_invalid_intervals():
    values = np.asarray([1.0, 100.0, 3.0])
    dataset = {
        "x": np.asarray([0.5, 1.5, 3.0]),
        "bins": np.asarray([0.0, 1.0, 2.0, 4.0]),
        "binwidth": np.asarray([1.0, 1.0, 2.0]),
        "xlim": np.asarray([0.0, 4.0]),
        "y": values,
        "scale": 1.0,
        "fitw": 1.0,
        "valid": np.asarray([True, False, True]),
    }
    finalize_uncertainties(
        dataset,
        [
            uncertainty_source(
                "normalization",
                0.1 * values,
                category="systematic",
                correlation="collective",
                scope="test:normalization",
                effect="multiplicative",
            )
        ],
    )
    observable = {
        "x": {
            "xlabel": "$x$",
            "ylabel": "$dN/dx$",
            "units": {"x": "1", "y": "1"},
        }
    }

    record = plot.histhepdata(
        hepdata={"x": dataset},
        obs=observable,
        density=True,
        density_uncertainty="scaled",
        MC_XS_SCALE=1.0,
    )["x"]

    np.testing.assert_allclose(record["hdata"].counts_scaled, [1.0 / 7.0, 0.0, 3.0 / 7.0])
    np.testing.assert_allclose(record["hdata"].errs_scaled, [0.1 / 7.0, 0.0, 0.3 / 7.0])
    plotted_values, plotted_errors = plot._mask_invalid_plot_bins(
        record["hdata"],
        record["hdata"].counts_scaled,
        record["hdata"].errs_scaled,
    )
    np.testing.assert_allclose(plotted_values[[0, 2]], [1.0 / 7.0, 3.0 / 7.0])
    np.testing.assert_allclose(plotted_errors[[0, 2]], [0.1 / 7.0, 0.3 / 7.0])
    assert np.isnan(plotted_values[1])
    assert np.isnan(plotted_errors[1])
    assert record["hdata"].integral() == pytest.approx(1.0)


# Check every ratio uncertainty convention and the exact self-ratio behavior
def test_ratio_plot_error_modes_explicit():
    numerator = np.asarray([10.0])
    denominator = np.asarray([5.0])
    numerator_error = np.asarray([2.0])
    denominator_error = np.asarray([1.0])

    combined = plot.ratio_plot_errors(
        numerator,
        denominator,
        numerator_error,
        denominator_error,
        mode="combined",
        reference=False,
    )
    separate = plot.ratio_plot_errors(
        numerator,
        denominator,
        numerator_error,
        denominator_error,
        mode="separate",
        reference=False,
    )
    reference = plot.ratio_plot_errors(
        denominator,
        denominator,
        denominator_error,
        denominator_error,
        mode="separate",
        reference=True,
    )

    np.testing.assert_allclose(combined, [2.0 * np.sqrt(0.2**2 + 0.2**2)])
    np.testing.assert_allclose(separate, [0.4])
    np.testing.assert_allclose(reference, [0.2])
    for mode in ("combined", "numerator"):
        np.testing.assert_allclose(
            plot.ratio_plot_errors(
                denominator,
                denominator,
                denominator_error,
                denominator_error,
                mode=mode,
                reference=True,
            ),
            [0.2],
        )
    np.testing.assert_allclose(
        plot.ratio_plot_errors(
            denominator,
            denominator,
            denominator_error,
            denominator_error,
            mode="none",
            reference=True,
        ),
        [0.0],
    )

    ratio_numerator = hist.hobj(
        counts=np.asarray([10.0, 6.0]),
        errs=np.asarray([2.0, 1.0]),
        bins=np.asarray([0.0, 1.0, 2.0]),
        cbins=np.asarray([0.5, 1.5]),
        valid=np.asarray([True, True]),
    )
    ratio_reference = hist.hobj(
        counts=np.asarray([5.0, 3.0]),
        errs=np.asarray([1.0, 0.5]),
        bins=np.asarray([0.0, 1.0, 2.0]),
        cbins=np.asarray([0.5, 1.5]),
        valid=np.asarray([True, False]),
    )
    plotted_ratio, plotted_error = plot._mask_invalid_ratio_bins(
        ratio_numerator,
        ratio_reference,
        np.asarray([2.0, 2.0]),
        np.asarray([0.4, 0.3]),
    )
    assert plotted_ratio[0] == pytest.approx(2.0)
    assert plotted_error[0] == pytest.approx(0.4)
    assert np.isnan(plotted_ratio[1])
    assert np.isnan(plotted_error[1])

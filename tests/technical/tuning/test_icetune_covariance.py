# Tests for the icetune data covariance construction
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json

import numpy as np
import pytest
from core.io import steering
from core.io.hepdata_reader import load_dataset_reader, read
from core.stats import cov, hist, objective, uncertainty
from core.stats.uncertainty import source_covariance, uncertainty_source
from core.tune import likelihood as icetune_likelihood
from scipy import sparse

from tests.technical.support.iceplot import CDIR


# Check excluded bins cannot set the numerical scale or invalidate a finite fit
@pytest.mark.parametrize("excluded", [1.0e50, np.nan, np.inf])
def test_correlated_chi2_scale_ignores_excluded_bins(excluded):
    covariance = np.eye(2)
    params = {
        "mc_prediction": np.array([2.0, excluded]), "data_values": np.array([1.0, excluded]),
        "mc_stat_uncertainty": np.array([0.5, excluded]), "fit_weights": np.array([1.0, 0.0]),
        "data_total_covariance": covariance,
    }
    structure = cov.prepare_covariance_structure(
        cov.covariance_decomposition(covariance, np.empty((2, 0)))
    )
    for result in (objective.correlated_chi2_arrays(**params),
                   objective.structured_chi2_arrays(**params, covariance_structure=structure)):
        assert result["valid"]
        assert result["chi2"] == pytest.approx(1.0 / 1.25)
        assert result["rank"] == 1


# Reject nonfinite active data instead of silently fitting only the remaining bins
@pytest.mark.parametrize("field", ["mc_prediction", "data_values", "mc_stat_uncertainty"])
def test_correlated_chi2_nonfinite_active_values(field):
    params = {
        "mc_prediction": np.array([2.0, 2.0]), "data_values": np.array([1.0, 1.0]),
        "mc_stat_uncertainty": np.array([0.5, 0.5]), "fit_weights": np.ones(2),
        "data_total_covariance": np.eye(2),
    }
    params[field][1] = np.nan
    assert not objective.correlated_chi2_arrays(**params)["valid"]


# Compute one minimal histogram comparison item
def histogram_item(values, errors, *, name):
    bins = np.array([0.0, 1.0, 2.0])
    source = {
        "name": name,
        "category": "statistical",
        "correlation": "uncorrelated",
        "up": np.asarray(errors, dtype=float),
        "down": np.asarray(errors, dtype=float),
    }
    return {
        "hdata": hist.hobj(
            np.asarray(values, dtype=float),
            np.asarray(errors, dtype=float),
            bins,
            np.array([0.5, 1.5]),
        ),
        "fitw": 1.0,
        "uncertainties": [source],
    }


# Compute one unit-density item with its transformed statistical covariance
def density_histogram_item(values, errors, *, name):
    bins = np.array([0.0, 1.0, 2.0])
    values = np.asarray(values, dtype=float)
    source = uncertainty_source(
        name,
        errors,
        category="statistical",
        correlation="uncorrelated",
    )
    transformed = uncertainty.transform_source(
        source,
        values=values,
        binwidth=np.ones(2),
        binscale=1.0,
        density=True,
        density_uncertainty="shape",
    )
    covariance = source_covariance(transformed)
    return {
        "hdata": hist.hobj(
            values / np.sum(values),
            np.sqrt(np.diag(covariance)),
            bins,
            np.array([0.5, 1.5]),
        ),
        "fitw": 1.0,
        "uncertainties": [transformed],
    }


# Add a block and low-rank decomposition to one minimal covariance payload
def structured_payload(results, block, vectors):
    """Return one covariance payload with automatic independent components"""
    return {
        "layout": cov.build_layout(results["data"]),
        "data_total_covariance": block + vectors @ vectors.T,
        **cov.covariance_decomposition(block, vectors),
    }


# Check covariance chi-square rescales tiny values before forming variances
def test_chi2_mc_error_underflow():
    result = objective.correlated_chi2_arrays(
        mc_prediction=np.array([2.0e-200]),
        data_values=np.array([1.0e-200]),
        mc_stat_uncertainty=np.array([1.0e-200]),
        fit_weights=np.ones(1),
        data_total_covariance=np.zeros((1, 1)),
    )

    assert result["valid"]
    assert result["rank"] == 1
    assert result["chi2"] == pytest.approx(1.0)


# Check active measurements and explicitly identified theory references without quoted uncertainties
def test_all_hepdata_hists_error_sources():
    for card_path in sorted((CDIR / "icepack").rglob("dataset.json")):
        relative = card_path.relative_to(CDIR).as_posix()
        dataset, resolved = steering.load_dataset(relative, cdir=str(CDIR))
        if not dataset["active"] or dataset["type"] == "RAW_SCALAR":
            continue
        for subset in dataset["sets"]:
            if not subset.get("data", True):
                continue
            for histogram in subset["hist"]:
                filename = steering.resolve_data_reference(
                    histogram["file"],
                    datapath=dataset["datapath"],
                    dataset_path=resolved,
                    cdir=str(CDIR),
                )
                table = read(
                    dataset["reader"],
                    dataset_path=resolved,
                    cdir=str(CDIR),
                    filename=filename,
                    file_filter=histogram.get("file_filter"),
                    rebin_factor=histogram.get("rebin_factor"),
                    hist=histogram,
                    dataset=subset,
                )
                groups = table["uncertainties"]
                if not isinstance(groups, dict):
                    groups = {None: groups}
                for key, sources in groups.items():
                    if not sources:
                        assert isinstance(table.get("reference"), str) and table["reference"].strip(), relative
                        errors = table["y_err"] if key is None else table["y_err"][key]
                        assert np.count_nonzero(errors) == 0, relative
                    size = len(table["y"] if key is None else table["y"][key])
                    for source in sources:
                        assert source["category"] in {
                            "statistical",
                            "systematic",
                            "combined",
                        }
                        assert source_covariance(source).shape == (size, size)


# Check the STAR luminosity source spans different observables
def test_star_luminosity_scope_is_shared():
    reader = load_dataset_reader(
        "icepack/SOFTCEP/STAR_1792394/pipi/dataset.json",
        cdir=str(CDIR),
    )
    directory = CDIR / "HEPData" / "CEP" / "HEPData-ins1792394-v1-json"
    mass = reader.read(str(directory / "Figure8.json"))
    rapidity = reader.read(str(directory / "Figure10(left).json"))
    mass_lumi = next(source for source in mass["uncertainties"] if source["name"] == "luminosity")
    rapidity_lumi = next(
        source for source in rapidity["uncertainties"] if source["name"] == "luminosity"
    )

    assert mass_lumi["correlation"] == "collective"
    assert mass_lumi["scope"] == rapidity_lumi["scope"]
    assert mass_lumi["scope"] == "STAR_1792394:luminosity"


# Check one collective nuisance produces exactly one global outer product
def test_collective_source_builds_cross_obs_cov_once():
    first = histogram_item([2.0, 2.0], [1.0, 1.0], name="stat_a")
    second = histogram_item([2.0, 2.0], [1.0, 1.0], name="stat_b")
    first["uncertainties"].append(
        uncertainty_source(
            "luminosity",
            [0.2, 0.3],
            category="systematic",
            correlation="collective",
            scope="test:luminosity",
        )
    )
    second["uncertainties"].append(
        uncertainty_source(
            "luminosity",
            [0.4, 0.5],
            category="systematic",
            correlation="collective",
            scope="test:luminosity",
        )
    )
    data = [[{"a": first, "b": second}]]
    layout = cov.build_layout(data)
    data_stat_covariance, data_syst_covariance, data_count_stat_uncertainty, _ = (
        cov.assemble_data_covariances(data, layout)
    )

    luminosity = np.array([0.2, 0.3, 0.4, 0.5])
    assert data_stat_covariance == pytest.approx(np.eye(4))
    assert data_syst_covariance == pytest.approx(np.outer(luminosity, luminosity))
    assert data_count_stat_uncertainty == pytest.approx(np.ones(4))


# Check collective sources do not merge independent validation blocks
def test_collective_source_cov_independent_blocks(monkeypatch):
    first = histogram_item([2.0, 2.0], [1.0, 1.0], name="stat_a")
    second = histogram_item([2.0, 2.0], [1.0, 1.0], name="stat_b")
    for item, shift in ((first, [0.2, 0.3]), (second, [0.4, 0.5])):
        item["uncertainties"].append(
            uncertainty_source(
                "luminosity",
                shift,
                category="systematic",
                correlation="collective",
                scope="test:luminosity",
            )
        )
    data = [[{"a": first, "b": second}]]
    obs = [
        [
            {
                "a": {"bins": np.array([0.0, 1.0, 2.0])},
                "b": {"bins": np.array([0.0, 1.0, 2.0])},
            }
        ]
    ]
    dimensions = []
    original_eigvalsh = np.linalg.eigvalsh

    def tracked_eigvalsh(matrix, *args, **kwargs):
        dimensions.append(len(matrix))
        return original_eigvalsh(matrix, *args, **kwargs)

    monkeypatch.setattr(np.linalg, "eigvalsh", tracked_eigvalsh)
    payload = cov.build_data_covariance(
        data=data,
        obs=obs,
        mcdata_by_dataset=[None],
    )
    cov.validate_data_covariance(payload)

    assert max(dimensions) == 1
    assert payload["data_collective_vectors"].shape == (4, 1)
    assert len(np.unique(payload["data_component_labels"])) == 4


# Check collective subtraction tolerates floating point cancellation at data scale
def test_collective_cov_roundoff():
    vectors = np.array(
        [
            [1.8905338179353308e7, -5.2274844148074739e7],
            [-4.1306354339189343e7, -2.4414673826398557e8],
        ]
    )
    covariance = vectors @ vectors.T

    cov.validate_local_covariance(
        covariance,
        vectors,
        np.arange(2, dtype=np.int64),
        name="collective covariance",
    )


# Check the scaled subtraction tolerance still rejects physical negative variance
def test_collective_negative_variance():
    vectors = np.zeros((2, 0), dtype=float)
    covariance = np.diag([1.0, -1.0])

    with pytest.raises(ValueError, match="not positive semidefinite"):
        cov.validate_local_covariance(
            covariance,
            vectors,
            np.arange(2, dtype=np.int64),
            name="local covariance",
        )


# Check MC cross-correlation scaling excludes non-counting statistical sources
def test_count_error_excludes_collective():
    item = histogram_item([10.0, 20.0], [1.0, 2.0], name="counting")
    item["uncertainties"].append(
        uncertainty_source(
            "selection",
            np.array([3.0, 4.0]),
            category="statistical",
            correlation="collective",
            scope="selection",
        )
    )
    data = [[{"obs": item}]]
    data_stat_covariance, _, data_count_stat_uncertainty, _ = cov.assemble_data_covariances(
        data, cov.build_layout(data)
    )

    assert np.diag(data_stat_covariance) == pytest.approx([10.0, 20.0])
    assert data_count_stat_uncertainty == pytest.approx([1.0, 2.0])


# Check a pure normalization nuisance cancels from a unit-density shape
def test_density_transform_cancels_collective_norm():
    values = np.array([2.0, 3.0])
    source = uncertainty_source(
        "normalization",
        0.1 * values,
        category="systematic",
        correlation="collective",
        scope="test:normalization",
        effect="multiplicative",
    )
    transformed = uncertainty.transform_source(
        source,
        values=values,
        binwidth=np.ones(2),
        binscale=1.0,
        density=True,
        density_uncertainty="shape",
    )

    assert transformed["shift"] == pytest.approx(np.zeros(2), abs=1.0e-15)


# Check excluded measurement intervals do not enter density normalization
def test_density_jacobian_excludes_invalid_intervals():
    values = np.asarray([1.0, np.nan, 3.0])
    widths = np.asarray([1.0, 1.0, 2.0])
    valid = np.asarray([True, False, True])

    jacobian = uncertainty.density_jacobian(values, widths, valid=valid)
    scale = uncertainty.density_jacobian(values, widths, valid=valid, shape=False)
    density, errors = hist.normalize_bin_heights(
        values,
        np.asarray([0.1, np.nan, 0.3]),
        widths,
        valid=valid,
    )

    active_values = np.where(valid, values, 0.0)
    np.testing.assert_allclose(jacobian @ active_values, np.zeros(3), atol=1.0e-15)
    np.testing.assert_allclose((widths * valid) @ jacobian, np.zeros(3), atol=1.0e-15)
    np.testing.assert_allclose(jacobian[~valid], np.zeros((1, 3)))
    np.testing.assert_allclose(np.diag(scale), [1.0 / 7.0, 0.0, 1.0 / 7.0])
    np.testing.assert_allclose(density, [1.0 / 7.0, 0.0, 3.0 / 7.0])
    np.testing.assert_allclose(errors, [0.1 / 7.0, 0.0, 0.3 / 7.0])


# Check unit density rejects a histogram with no measured support
def test_density_transform_rejects_fully_masked_hist():
    values = np.asarray([1.0, 2.0])
    widths = np.asarray([0.5, 0.5])
    valid = np.asarray([False, False])

    with pytest.raises(ValueError, match="positive finite integral"):
        uncertainty.density_jacobian(values, widths, valid=valid)
    with pytest.raises(ValueError, match="positive finite integral"):
        uncertainty.density_jacobian(values, widths, valid=valid, shape=False)


# Check hyperbin assignments reject events outside each proton category
def test_hyperbin_assignment_masks_nonmembers():
    category_bins = {
        "proton_bin": (
            np.array([0.0, 1.0]),
            np.array([0.0, 1.0]),
            np.array([0.0, 1.0, 2.0]),
        )
    }
    mcdata = {
        "data": {
            "mass": np.array(
                [
                    [0.5, 0.5, 0.25],
                    [1.5, 0.5, 0.25],
                    [0.5, 0.5, 1.25],
                ]
            )
        },
        "event_ids": np.array([10, 11, 12]),
    }
    obs = {
        "mass": {
            "bins": category_bins,
            "symmetrized_fill": False,
        }
    }

    assignments = cov.subset_assignments(mcdata, obs)
    event_ids, local_bins = assignments["mass_[0.00-1.00]_[0.00-1.00]"]

    assert event_ids.tolist() == [10, 11, 12]
    assert local_bins.tolist() == [0, -1, 1]


# Check unit- and unequal-weight samples produce valid overlap correlations
@pytest.mark.parametrize(
    ("weights", "generator_weighted"),
    [
        (np.ones(4), False),
        (np.array([1.0, 2.0, 3.0, 4.0]), True),
    ],
)
def test_mc_sample_correlates_obs_projections(
    weights,
    generator_weighted,
):
    data = [
        [
            {
                "a": histogram_item([2.0, 2.0], [1.0, 1.0], name="stat_a"),
                "b": histogram_item([2.0, 2.0], [1.0, 1.0], name="stat_b"),
            }
        ]
    ]
    obs = [
        [
            {
                "a": {"bins": np.array([0.0, 1.0, 2.0])},
                "b": {"bins": np.array([0.0, 1.0, 2.0])},
            }
        ]
    ]
    mcdata = [
        [
            {
                "data": {
                    "a": np.array([0.25, 0.25, 1.25, 1.25]),
                    "b": np.array([0.25, 0.25, 1.25, 1.25]),
                },
                "event_ids": np.arange(4),
                "sample_event_ids": np.arange(4),
                "sample_weights": weights,
            }
        ]
    ]

    payload = cov.build_data_covariance(
        data=data,
        obs=obs,
        mcdata_by_dataset=mcdata,
        mc_correlation_datacards=[
            {
                "nevents": 4,
                "weighted": generator_weighted,
            }
        ],
    )

    assert payload["mc_cross_correlation"][0, 2] == pytest.approx(1.0)
    assert payload["mc_cross_correlation"][0, 3] == pytest.approx(0.0)
    assert payload["data_stat_correlation"][0, 2] == pytest.approx(1.0)
    assert payload["data_total_covariance"][0, 2] == pytest.approx(1.0)
    sample = payload["mc_correlation_samples"][0]
    assert sample["event_count_requested"] == 4
    assert sample["event_count_read"] == 4
    assert sample["generator_weighted"] is generator_weighted
    assert sample["event_weights_nonuniform"] is generator_weighted
    if generator_weighted:
        assert sample["effective_event_count"] < 4.0


# Check independently generated routed samples do not acquire false event overlap
def test_routed_mc_samples_independent_cross_cov():
    data = [
        [
            {"a": histogram_item([2.0, 2.0], [1.0, 1.0], name="stat_a")},
            {"b": histogram_item([2.0, 2.0], [1.0, 1.0], name="stat_b")},
        ]
    ]
    obs = [
        [
            {"a": {"bins": np.array([0.0, 1.0, 2.0])}},
            {"b": {"bins": np.array([0.0, 1.0, 2.0])}},
        ]
    ]
    first = {
        "data": {"a": np.array([0.25, 0.25, 1.25, 1.25])},
        "event_ids": np.arange(4),
        "sample_event_ids": np.arange(4),
        "sample_weights": np.ones(4),
    }
    second = {
        "data": {"b": np.array([0.25, 0.25, 1.25, 1.25])},
        "event_ids": np.arange(4),
        "sample_event_ids": np.arange(4),
        "sample_weights": np.ones(4),
    }
    grouped = [
        {
            "sample_groups": [
                {"mcdata": [first], "sample": "first", "set_indices": [0]},
                {"mcdata": [second], "sample": "second", "set_indices": [1]},
            ]
        }
    ]

    payload = cov.build_data_covariance(
        data=data,
        obs=obs,
        mcdata_by_dataset=grouped,
        mc_correlation_datacards=[{"nevents": 4, "weighted": False}],
    )

    np.testing.assert_allclose(payload["mc_cross_correlation"], np.zeros((4, 4)))
    assert [entry["sample"] for entry in payload["mc_correlation_samples"]] == [
        "first",
        "second",
    ]


# Check absolute and unit-density modes freeze their distinct covariance models
def test_xs_density_modes_build_expected_covariances():
    obs = [
        [
            {
                "a": {"bins": np.array([0.0, 1.0, 2.0])},
                "b": {"bins": np.array([0.0, 1.0, 2.0])},
            }
        ]
    ]
    mcdata = [
        [
            {
                "data": {
                    "a": np.array([0.25, 0.25, 1.25, 1.25]),
                    "b": np.array([0.25, 0.25, 1.25, 1.25]),
                },
                "event_ids": np.arange(4),
                "sample_event_ids": np.arange(4),
                "sample_weights": np.ones(4),
            }
        ]
    ]
    cross_section_data = [
        [
            {
                "a": histogram_item([2.0, 2.0], [1.0, 1.0], name="stat_a"),
                "b": histogram_item([2.0, 2.0], [1.0, 1.0], name="stat_b"),
            }
        ]
    ]
    density_data = [
        [
            {
                "a": density_histogram_item([2.0, 2.0], [1.0, 1.0], name="stat_a"),
                "b": density_histogram_item([2.0, 2.0], [1.0, 1.0], name="stat_b"),
            }
        ]
    ]

    cross_section = cov.build_data_covariance(
        data=cross_section_data,
        obs=obs,
        mcdata_by_dataset=mcdata,
    )
    density = cov.build_data_covariance(
        data=density_data,
        obs=obs,
        mcdata_by_dataset=mcdata,
    )

    assert {record["mc_correlation_model"] for record in cross_section["layout"]} == {"counting"}
    assert {record["mc_correlation_model"] for record in density["layout"]} == {"normalized"}
    assert cross_section["data_total_covariance"][0, 1] == pytest.approx(0.0)
    assert cross_section["data_total_covariance"][0, 2] == pytest.approx(1.0)
    assert cross_section["data_total_covariance"][0, 3] == pytest.approx(0.0)
    assert density["data_total_covariance"][0, 1] < 0.0
    assert density["data_total_covariance"][0, 2] > 0.0
    assert density["data_total_covariance"][0, 3] < 0.0
    for payload in (cross_section, density):
        cov.validate_data_covariance(payload)
        assert payload["data_total_covariance"] == pytest.approx(
            payload["data_stat_covariance"] + payload["data_syst_covariance"]
        )
        assert payload["data_stat_covariance"] == pytest.approx(
            payload["data_stat_correlation"]
            * np.outer(
                payload["data_stat_uncertainty"],
                payload["data_stat_uncertainty"],
            )
        )


# Check incompatible finite-sample cross terms are reduced to a valid matrix
def test_mc_cross_cov_psd():
    factor, covariance = cov.scale_mc_cross_covariance_to_psd(
        np.eye(2),
        np.array([[0.0, 2.0], [2.0, 0.0]]),
    )

    assert factor == pytest.approx(0.5, rel=1.0e-8)
    assert np.min(np.linalg.eigvalsh(covariance)) >= -1.0e-10


# Check independent variances cannot loosen the PSD condition of another component
@pytest.mark.parametrize("variance", [1.0e-8, 1.0e8])
def test_sparse_mc_cov_psd_scale(variance):
    statistical = sparse.diags([1.0, 1.0, variance], format="csr")
    cross = sparse.csr_matrix([[0.0, 2.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    factor = cov.sparse_mc_covariance_scale(statistical, cross, np.empty((3, 0)))
    corrected = (statistical + factor * cross).toarray()
    assert factor == pytest.approx(0.5, rel=1.0e-8)
    assert np.min(np.linalg.eigvalsh(corrected[:2, :2])) >= -1.0e-10


# Check diagonal covariance reproduces the independent-bin quadratic cost
def test_diagonal_cov_matches_independent_bin_chi2():
    data_item = histogram_item([2.0, 2.0], [1.0, 1.5], name="statistical")
    data_item["fitw"] = 1.7
    mc_item = histogram_item([3.0, 0.5], [0.2, 0.3], name="mc_statistical")
    results = {
        "data": [[{"x": data_item}]],
        "mc": [[{"x": mc_item}]],
    }
    independent, ndf = objective.chi2_cost(
        h_mc=mc_item["hdata"],
        h_data=data_item["hdata"],
        rho="quadratic",
    )
    payload = {
        "layout": cov.build_layout(results["data"]),
        "data_total_covariance": np.diag(np.array([1.0, 1.5]) ** 2),
    }
    joint = objective.correlated_chi2(results=results, payload=payload)

    assert joint["chi2"] == pytest.approx(1.7 * independent)
    assert joint["ndf"] == ndf


# Check fit weights form a symmetric positive quadratic objective
def test_fit_weight_quadratic_form():
    residual = np.array([0.7, -0.4, 3.0])
    covariance = np.array(
        [
            [1.2, 0.3, 0.2],
            [0.3, 0.9, -0.1],
            [0.2, -0.1, 1.5],
        ]
    )
    fitw = np.array([4.0, 0.25, 0.0])
    mc_error = np.array([0.2, 0.3, 0.4])
    result = objective.covariance_cost(
        residual=residual,
        covariance=covariance,
        fit_weights=fitw,
        mc_stat_uncertainty=mc_error,
        return_contributions=True,
    )
    active = fitw > 0.0
    weighted_residual = np.sqrt(fitw[active]) * residual[active]
    inverse_residual = np.linalg.solve(
        covariance[np.ix_(active, active)],
        weighted_residual,
    )
    expected = float(weighted_residual @ inverse_residual)
    sigma = 2.0 * np.linalg.norm(np.sqrt(fitw[active]) * inverse_residual * mc_error[active])

    assert result["valid"]
    assert result["cost"] == pytest.approx(expected)
    assert result["sigma_cost"] == pytest.approx(sigma)
    assert result["rank"] == 2
    assert np.sum(result["contributions"][active]) == pytest.approx(expected)
    assert np.isnan(result["contributions"][~active]).all()


# Check a common absolute Gaussian fit weight scales cost but not rank
def test_common_fit_weight_scaling():
    kwargs = {
        "residual": np.array([0.7, -0.4]),
        "covariance": np.array([[1.2, 0.3], [0.3, 0.9]]),
        "mc_stat_uncertainty": np.array([0.2, 0.3]),
    }
    unit = objective.covariance_cost(**kwargs, fit_weights=np.ones(2))
    scaled = objective.covariance_cost(**kwargs, fit_weights=np.full(2, 3.5))

    assert scaled["cost"] == pytest.approx(3.5 * unit["cost"])
    assert scaled["sigma_cost"] == pytest.approx(3.5 * unit["sigma_cost"])
    assert scaled["rank"] == unit["rank"]


# Check structured chi-square matches the original dense global eigensolve
def test_structured_chi2_matches_dense_cov(monkeypatch):
    data_a = histogram_item([2.0, 3.0], [1.0, 1.0], name="data_a")
    data_b = histogram_item([4.0, 2.0], [1.0, 1.0], name="data_b")
    data_a["fitw"] = 1.5
    data_b["fitw"] = 0.7
    results = {
        "data": [[{"a": data_a, "b": data_b}]],
        "mc": [
            [
                {
                    "a": histogram_item([2.4, 2.5], [0.2, 0.3], name="mc_a"),
                    "b": histogram_item([3.5, 2.3], [0.4, 0.2], name="mc_b"),
                }
            ]
        ],
    }
    block = np.array(
        [
            [1.0, 0.2, 0.0, 0.0],
            [0.2, 1.4, 0.0, 0.0],
            [0.0, 0.0, 1.2, -0.1],
            [0.0, 0.0, -0.1, 0.9],
        ]
    )
    vectors = np.array([[0.2], [0.3], [0.4], [0.25]])
    payload = structured_payload(results, block, vectors)
    dense = objective.correlated_chi2(
        results,
        {
            "layout": payload["layout"],
            "data_total_covariance": payload["data_total_covariance"],
        },
    )

    dimensions = []
    original_eigh = np.linalg.eigh

    def tracked_eigh(matrix, *args, **kwargs):
        dimensions.append(len(matrix))
        return original_eigh(matrix, *args, **kwargs)

    monkeypatch.setattr(np.linalg, "eigh", tracked_eigh)
    structured = objective.correlated_chi2(results, payload)

    assert max(dimensions) == 2
    assert structured["chi2"] == pytest.approx(dense["chi2"], rel=1.0e-12)
    assert structured["sigma_chi2"] == pytest.approx(dense["sigma_chi2"], rel=1.0e-12)
    assert structured["ndf"] == dense["ndf"]
    assert [record["objective_value"] for record in structured["observables"]] == pytest.approx(
        [record["objective_value"] for record in dense["observables"]]
    )


# Check unrelated physical blocks retain modes with very different numerical scales
def test_structured_cov_independent_spectral_scales():
    block = np.diag([1.0, 1.0e-16])
    structure = cov.prepare_covariance_structure(
        cov.covariance_decomposition(
            block,
            np.empty((2, 0), dtype=float),
        )
    )
    result = objective.structured_covariance_cost(
        residual=np.array([1.0, 1.0e-8]),
        structure=structure,
        fit_weights=np.ones(2),
        covariance_transform=np.ones(2),
        diagonal_uncertainty=np.zeros(2),
    )

    assert result is not None
    assert result["rank"] == 2
    assert result["cost"] == pytest.approx(2.0)


# Check normalized singular modes are truncated inside their physical block
def test_structured_cov_handles_local_singular_mode():
    block = np.array([[1.0, -1.0], [-1.0, 1.0]])
    vectors = np.array([[0.2], [-0.2]])
    structure = cov.prepare_covariance_structure(cov.covariance_decomposition(block, vectors))
    result = objective.structured_covariance_cost(
        residual=np.array([0.3, -0.3]),
        structure=structure,
        fit_weights=np.ones(2),
        covariance_transform=np.ones(2),
        diagonal_uncertainty=np.zeros(2),
        return_contributions=True,
    )
    dense = objective.covariance_cost(
        residual=np.array([0.3, -0.3]),
        covariance=block + vectors @ vectors.T,
        fit_weights=np.ones(2),
        return_contributions=True,
    )

    assert result is not None
    assert result["rank"] == 1
    assert result["cost"] == pytest.approx(dense["cost"])
    assert result["ndf_contributions"] == pytest.approx(dense["ndf_contributions"])


# Check structured ratio2 matches the original dense transformed covariance
@pytest.mark.parametrize("rho", ["quadratic", "huber"])
def test_structured_ratio2_matches_dense_cov(rho):
    data_a = histogram_item([2.0, 3.0], [0.3, 0.4], name="data_a")
    data_b = histogram_item([4.0, 2.0], [0.5, 0.3], name="data_b")
    data_a["fitw"] = 1.2
    data_b["fitw"] = 0.8
    results = {
        "data": [[{"a": data_a, "b": data_b}]],
        "mc": [
            [
                {
                    "a": histogram_item([2.4, 2.5], [0.2, 0.3], name="mc_a"),
                    "b": histogram_item([3.5, 2.3], [0.4, 0.2], name="mc_b"),
                }
            ]
        ],
    }
    block = np.array(
        [
            [0.09, 0.01, 0.0, 0.0],
            [0.01, 0.16, 0.0, 0.0],
            [0.0, 0.0, 0.25, -0.02],
            [0.0, 0.0, -0.02, 0.09],
        ]
    )
    vectors = np.array([[0.05], [0.08], [0.1], [0.06]])
    payload = structured_payload(results, block, vectors)
    dense = objective.correlated_ratio2(
        results,
        {
            "layout": payload["layout"],
            "data_total_covariance": payload["data_total_covariance"],
        },
        rho=rho,
    )
    structured = objective.correlated_ratio2(results, payload, rho=rho)

    assert structured["cost"] == pytest.approx(dense["cost"], rel=1.0e-12)
    assert structured["ndf"] == dense["ndf"]
    assert structured["rank"] == dense["rank"]
    assert sum(item["objective_value"] for item in structured["observables"]) == pytest.approx(
        structured["cost"]
    )


# Check observable diagnostics partition one joint precision-weighted objective
def test_joint_chi2_obs_partition_additive():
    data_a = histogram_item([2.0], [1.0], name="statistical")
    data_b = histogram_item([3.0], [1.0], name="statistical")
    data_a["fitw"] = 1.0
    data_b["fitw"] = 2.0
    results = {
        "data": [[{"a": data_a, "b": data_b}]],
        "mc": [
            [
                {
                    "a": histogram_item([2.5], [0.2], name="mc_statistical"),
                    "b": histogram_item([2.0], [0.3], name="mc_statistical"),
                }
            ]
        ],
    }
    payload = {
        "layout": cov.build_layout(results["data"]),
        "data_total_covariance": np.array([[1.0, 0.6], [0.6, 1.0]]),
    }

    joint = objective.correlated_chi2(results=results, payload=payload)

    assert sum(item["objective_value"] for item in joint["observables"]) == pytest.approx(
        joint["chi2"]
    )
    assert sum(item["ndf"] for item in joint["observables"]) == pytest.approx(joint["rank"])


# Check diagonal covariance reproduces the independent asinh ratio cost
def test_diagonal_cov_matches_independent_ratio2():
    data_item = histogram_item([2.0, 2.0], [1.0, 1.5], name="statistical")
    data_item["fitw"] = 1.7
    mc_item = histogram_item([3.0, 0.5], [0.2, 0.3], name="mc_statistical")
    results = {
        "data": [[{"x": data_item}]],
        "mc": [[{"x": mc_item}]],
    }
    independent, ndf = objective.ratio2_cost(
        h_mc=mc_item["hdata"],
        h_data=data_item["hdata"],
        rho="quadratic",
    )
    payload = {
        "layout": cov.build_layout(results["data"]),
        "data_total_covariance": np.diag(np.array([1.0, 1.5]) ** 2),
    }
    joint = objective.correlated_ratio2(results=results, payload=payload)

    assert joint["cost"] == pytest.approx(1.7 * independent)
    assert joint["ndf"] == ndf


# Check joint and independent ratio costs use the same MC validity mask
def test_correlated_ratio2_respects_mc_validity():
    data_item = histogram_item([2.0, 3.0], [0.4, 0.5], name="statistical")
    mc_item = histogram_item([2.5, 9.0], [0.2, 0.3], name="mc_statistical")
    mc_item["hdata"].valid[1] = False
    results = {"data": [[{"x": data_item}]], "mc": [[{"x": mc_item}]]}
    payload = {
        "layout": cov.build_layout(results["data"]),
        "data_total_covariance": np.diag(np.array([0.4, 0.5]) ** 2),
    }

    independent, ndf = objective.ratio2_cost(mc_item["hdata"], data_item["hdata"])
    joint = objective.correlated_ratio2(results=results, payload=payload)
    independent_chi2, chi2_ndf = objective.chi2_cost(mc_item["hdata"], data_item["hdata"])
    joint_chi2 = objective.correlated_chi2(results=results, payload=payload)

    assert joint["cost"] == pytest.approx(independent)
    assert joint["ndf"] == ndf == 1
    assert joint_chi2["chi2"] == pytest.approx(independent_chi2)
    assert joint_chi2["ndf"] == chi2_ndf == 1


# Check a residual outside the covariance support rejects the model
def test_cov_cost_rejects_exact_mode_disagreement():
    result = objective.covariance_cost(
        residual=np.array([1.0, -1.0]),
        covariance=np.ones((2, 2)),
        fit_weights=np.ones(2),
    )

    assert result["valid"] is False
    assert result["ndf"] == result["rank"] == 1


# Check the covariance cost is unchanged when MC and data are exchanged
def test_ratio2_mc_data_exchange():
    data_item = histogram_item([2.0, 2.0], [0.2, 0.2], name="data_statistical")
    mc_item = histogram_item([4.0, 1.0], [0.4, 0.1], name="mc_statistical")
    forward_results = {"data": [[{"x": data_item}]], "mc": [[{"x": mc_item}]]}
    reverse_results = {"data": [[{"x": mc_item}]], "mc": [[{"x": data_item}]]}
    forward_payload = {
        "layout": cov.build_layout(forward_results["data"]),
        "data_total_covariance": np.diag(np.array([0.2, 0.2]) ** 2),
    }
    reverse_payload = {
        "layout": cov.build_layout(reverse_results["data"]),
        "data_total_covariance": np.diag(np.array([0.4, 0.1]) ** 2),
    }

    forward = objective.correlated_ratio2(results=forward_results, payload=forward_payload)
    reverse = objective.correlated_ratio2(results=reverse_results, payload=reverse_payload)

    assert forward["valid"] is True
    assert reverse["valid"] is True
    assert reverse["ndf"] == forward["ndf"]
    assert reverse["cost"] == pytest.approx(forward["cost"])


# Check batched post-processing uses the direct generalized chi-square
def test_chi2_batch_vs_direct():
    covariance = np.array([[1.0, 0.6], [0.6, 2.0]])
    data = np.array([2.0, 3.0])
    predictions = np.array([[2.5, 2.5], [1.5, 4.0]])
    errors = np.array([[0.2, 0.3], [0.4, 0.1]])
    weights = np.array([1.0, 2.0])

    values, uncertainties = objective.correlated_chi2_batch(
        mc_prediction=predictions,
        data_values=data,
        mc_stat_uncertainty=errors,
        fit_weights=weights,
        data_total_covariance=covariance,
    )
    direct = [
        objective.correlated_chi2_arrays(
            mc_prediction=prediction,
            data_values=data,
            mc_stat_uncertainty=error,
            fit_weights=weights,
            data_total_covariance=covariance,
        )
        for prediction, error in zip(predictions, errors, strict=True)
    ]
    vectors = np.array([[0.5], [1.2]])
    block = covariance - vectors @ vectors.T
    structure = cov.prepare_covariance_structure(cov.covariance_decomposition(block, vectors))
    structured_values, structured_uncertainties = objective.correlated_chi2_batch(
        mc_prediction=predictions,
        data_values=data,
        mc_stat_uncertainty=errors,
        fit_weights=weights,
        data_total_covariance=covariance,
        covariance_structure=structure,
    )

    assert values == pytest.approx([item["chi2"] for item in direct])
    assert uncertainties == pytest.approx([item["sigma_chi2"] for item in direct])
    assert structured_values == pytest.approx(values)
    assert structured_uncertainties == pytest.approx(uncertainties)


# Check surrogate selections preserve covariance ordering and collective sources
def test_cov_structure_selection_matches_dense_cost():
    block = np.array(
        [
            [1.0, 0.2, 0.0, 0.0],
            [0.2, 1.1, 0.0, 0.0],
            [0.0, 0.0, 0.8, -0.1],
            [0.0, 0.0, -0.1, 1.3],
        ]
    )
    vectors = np.array([[0.2], [0.3], [0.4], [0.5]])
    covariance = block + vectors @ vectors.T
    selected = np.array([3, 0, 2])
    structure = cov.prepare_covariance_structure(cov.covariance_decomposition(block, vectors))
    selected_structure = cov.select_covariance_structure(structure, selected)
    mc = np.array([2.5, 1.2, 3.1])
    data = np.array([2.0, 1.0, 3.0])
    errors = np.array([0.2, 0.1, 0.3])
    weights = np.array([1.0, 2.0, 0.7])
    structured = objective.structured_chi2_arrays(
        mc_prediction=mc,
        data_values=data,
        mc_stat_uncertainty=errors,
        fit_weights=weights,
        data_total_covariance=covariance[np.ix_(selected, selected)],
        covariance_structure=selected_structure,
    )
    dense = objective.correlated_chi2_arrays(
        mc_prediction=mc,
        data_values=data,
        mc_stat_uncertainty=errors,
        fit_weights=weights,
        data_total_covariance=covariance[np.ix_(selected, selected)],
    )

    assert structured is not None
    assert structured["chi2"] == pytest.approx(dense["chi2"])
    assert structured["sigma_chi2"] == pytest.approx(dense["sigma_chi2"])


# Check full histogram manifests map exactly into fixed covariance order
def test_manifest_valid_cov_bins():
    layout = [
        {
            "dataset": 0,
            "subset": 0,
            "observable": "a",
            "bins": [0, 2],
            "start": 0,
            "stop": 2,
        },
        {
            "dataset": 0,
            "subset": 0,
            "observable": "b",
            "bins": [1],
            "start": 2,
            "stop": 3,
        },
    ]
    manifest = [
        {
            "dataset": 0,
            "subset": 0,
            "observable": "a",
            "start": 0,
            "stop": 3,
            "fitw": 2.0,
        },
        {
            "dataset": 0,
            "subset": 0,
            "observable": "b",
            "start": 3,
            "stop": 5,
            "fitw": 4.0,
        },
    ]

    assert cov.manifest_layout_indices(layout, manifest).tolist() == [0, 2, 4]
    assert cov.manifest_fit_weights(manifest).tolist() == [2.0, 2.0, 2.0, 4.0, 4.0]


# Check data covariance matrices are saved separately from readable metadata
def test_data_covariance_disk_format(tmp_path):
    data = [[{"x": histogram_item([2.0, 2.0], [1.0, 1.0], name="statistical")}]]
    payload = cov.build_data_covariance(data=data, obs=[{}], mcdata_by_dataset=[None])
    cov.validate_data_covariance(payload)
    npz_path, json_path = cov.save_data_covariance(
        payload,
        output_dir=str(tmp_path),
    )

    with np.load(npz_path) as arrays:
        assert "data_total_covariance" not in arrays.files
        np.testing.assert_allclose(cov.sparse_covariance_matrix(arrays, "stat_local").toarray(), np.eye(2))
        assert arrays["stat_local_values"].size == 2
    with open(json_path, encoding="utf-8") as source:
        metadata = json.load(source)
    assert metadata["matrix_file"] == "data_covariance.npz"
    assert len(metadata["matrix_sha256"]) == 64
    assert metadata["arrays"] == list(cov.DATA_COVARIANCE_ARRAYS)
    loaded = cov.load_data_covariance(
        npz_path=npz_path,
        json_path=json_path,
    )
    assert loaded["data_total_covariance"].tolist() == [
        [1.0, 0.0],
        [0.0, 1.0],
    ]


# Compare sparse collective chi-square values and derivatives with the complete covariance
def test_sparse_covariance_torch_derivatives():
    import torch

    data = [[{"a": histogram_item([2.0, 3.0], [1.0, 1.0], name="stat_a"),
              "b": histogram_item([4.0, 5.0], [1.0, 1.0], name="stat_b")}]]
    for name, shift in (("a", [0.2, 0.3]), ("b", [0.4, 0.5])):
        data[0][0][name]["uncertainties"].append(uncertainty_source(
            "luminosity", shift, category="systematic", correlation="collective", scope="shared"))
    payload = cov.build_data_covariance(data=data, obs=[{}], mcdata_by_dataset=[None])
    assert "data_total_covariance" not in payload
    structure = cov.prepare_covariance_structure(payload)
    shift = np.array([0.2, 0.3, 0.4, 0.5])
    covariance = np.eye(4) + np.outer(shift, shift)
    np.testing.assert_allclose(cov.covariance_selection(payload), covariance)
    point = torch.tensor([2.3, 2.8, 3.6, 5.2], dtype=torch.float64, requires_grad=True)

    # Differentiate the actual covariance objective including parameter dependent MC errors
    def evaluate(values, structured=True):
        arguments = dict(mc_prediction=values, data_values=np.arange(2.0, 6.0),
                         mc_stat_uncertainty=0.1 * values, fit_weights=np.ones(4))
        result = (objective.structured_chi2_arrays(**arguments, covariance_structure=structure) if structured else
                  objective.correlated_chi2_arrays(**arguments, data_total_covariance=covariance))
        assert result["valid"]
        return result["chi2"]

    structured, dense = evaluate(point), evaluate(point, False)
    torch.testing.assert_close(structured, dense)
    torch.testing.assert_close(torch.autograd.grad(structured, point)[0], torch.autograd.grad(dense, point)[0])
    assert torch.autograd.gradcheck(evaluate, (point,))
    assert torch.autograd.gradgradcheck(evaluate, (point,))


# Check new trial payloads contain only one objective convention
def test_likelihood_payload_objective_convention():
    payload = icetune_likelihood.build_joint_likelihood_payload(
        joint={"valid": True, "chi2": 12.0, "ndf": 6, "rank": 6, "bins": 6},
        covariance_mode="full",
    )

    assert payload["objective"] == {
        "name": "chi2",
        "value": 12.0,
        "ndf": 6.0,
        "reduced": 2.0,
    }
    assert "nll_total" not in payload
    assert "two_nll_total" not in payload
    assert "logL" not in payload


# Check joint histories retain additive diagnostics from the global solve
def test_joint_likelihood_payload_obs_impacts():
    observable = {
        "dataset": 0,
        "subset": 0,
        "observable": "mass",
        "ndf": 3.0,
        "valid": True,
        "objective_value": 3.0,
        "weight": 1.5,
    }
    payload = icetune_likelihood.build_joint_likelihood_payload(
        joint={
            "valid": True,
            "chi2": 12.0,
            "ndf": 6,
            "rank": 6,
            "bins": 6,
            "observables": [observable],
        },
        covariance_mode="full",
    )

    assert payload["observables"] == [observable]
    assert payload["components"]["observable_diagnostics"] == "joint_precision_partition_additive"

# Check the common icepack Monte Carlo precision criteria
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import pytest
from core.stats.validation import check_physics

from tests.physics.validation.test_icepacks import assert_event_kinematics
from tests.technical.support.hepmc import write_muon_events


# Build one minimal histogram report with configurable precision diagnostics
def precision_report(validation, **sample_updates):
    sample = {
        "label": "GRANIITTI",
        "empty": False,
        "mc_max_rel_uncertainty": 0.20,
        "mc_zero_bin_count": 0,
        "mc_valid_bin_count": 4,
        "mc_selected_events": 1000,
        "mc_ess_fraction": 0.25,
        "mc_effective_events": 250.0,
    }
    sample.update(sample_updates)
    return {
        "validation": validation,
        "sets": [
            {
                "name": "fiducial",
                "observables": [
                    {
                        "observable": "mass",
                        "kind": "histogram",
                        "samples": [sample],
                    }
                ],
            }
        ],
    }


# Accept finite binwise precision and effective population above the card bounds
def test_mc_precision_accepts_well_populated_hist():
    report = precision_report(
        {
            "max_mc_rel_uncertainty": 0.30,
            "min_mc_effective_events": 100,
            "min_mc_ess_fraction": 0.10,
        }
    )

    check_physics(report)


# Reject an observable whose worst populated bin is too imprecise
def test_mc_precision_rejects_large_relative_error():
    report = precision_report(
        {"max_mc_rel_uncertainty": 0.30},
        mc_max_rel_uncertainty=0.31,
    )

    with pytest.raises(ValueError, match="mc_max_rel_uncertainty=0.31 exceeds 0.3"):
        check_physics(report)


# Reject zero MC bins because their relative uncertainty is undefined
def test_mc_precision_rejects_zero_bins():
    report = precision_report(
        {"max_mc_rel_uncertainty": 0.30},
        mc_zero_bin_count=1,
    )

    with pytest.raises(ValueError, match="zero content and undefined relative uncertainty"):
        check_physics(report)


# Reject insufficient effective events and a low source ESS fraction
def test_mc_precision_rejects_source_population():
    report = precision_report(
        {
            "min_mc_effective_events": 100,
            "min_mc_ess_fraction": 0.10,
        },
        mc_effective_events=50.0,
        mc_ess_fraction=0.05,
    )

    with pytest.raises(ValueError, match="mc_effective_events=50 is below 100"):
        check_physics(report)


# Reject population bounds when a report cannot provide source diagnostics
def test_mc_missing_population():
    report = precision_report({"min_mc_effective_events": 100})
    sample = report["sets"][0]["observables"][0]["samples"][0]
    sample.pop("mc_effective_events")

    with pytest.raises(ValueError, match="mc_effective_events is unavailable"):
        check_physics(report)


# Construct an independent F/C comparison with an untested reference row
def sampling_report():
    report = precision_report(
        {
            "mc_reference": "F",
            "require_comparison": True,
            "require_differential_comparison": True,
            "require_fiducial_integral_comparison": True,
            "max_abs_integral_pull": 5.0,
            "max_chi2_ndf": 3.0,
            "max_shape_l1": 0.30,
            "max_integral_factor": 1.12,
            "min_mc_effective_events": 100,
        },
        comparison_status="reference",
    )
    samples = report["sets"][0]["observables"][0]["samples"]
    samples.append(
        {
            **samples[0],
            "label": "C",
            "comparison_status": "compared",
            "integral": 10.1,
            "data_integral": 10.0,
            "integral_pull": 0.2,
            "chi2_ndf": 1.1,
            "shape_l1": 0.04,
        }
    )
    return report


# A reference needs population checks but has no self-comparison chi-square or pull
def test_reference_no_self_comparison():
    check_physics(sampling_report())


# Each physical discrepancy must fail even when other comparison diagnostics pass
@pytest.mark.parametrize(
    ("field", "value", "failure"),
    [
        ("integral_pull", 5.1, "integral_pull"),
        ("integral_pull", -5.1, "integral_pull"),
        ("shape_l1", 0.31, "shape_l1"),
        ("integral", 12.0, "integral factor"),
        ("integral_pull", None, "integral_pull is unavailable"),
        ("shape_l1", -0.1, "shape_l1"),
        ("data_integral", None, "integral factor requires both integrals"),
        ("data_integral", float("nan"), "finite positive integrals"),
        ("integral", float("inf"), "finite positive integrals"),
    ],
)
def test_sampling_missing_or_biased(field, value, failure):
    report = sampling_report()
    report["sets"][0]["observables"][0]["samples"][1][field] = value
    with pytest.raises(ValueError, match=failure):
        check_physics(report)


# A chi-square discrepancy remains diagnostic even with an explicit bound
@pytest.mark.parametrize("value", [3.1, -1.0])
def test_sampling_chi2_does_not_raise(value):
    report = sampling_report()
    sample = report["sets"][0]["observables"][0]["samples"][1]
    sample["chi2_ndf"] = value
    check_physics(report)
    assert sample["chi2_ndf"] == pytest.approx(value)


# A noisy or empty reference cannot qualify the candidate as a physics agreement
def test_sampling_empty_reference():
    report = sampling_report()
    report["sets"][0]["observables"][0]["samples"][0]["mc_effective_events"] = 10
    with pytest.raises(ValueError, match="mc_effective_events"):
        check_physics(report)


# A single rate bin must use the requested pull tolerance instead of its squared value
def test_sampling_xs_integral_pull_bound():
    report = sampling_report()
    report["validation"]["require_differential_comparison"] = False
    observable = report["sets"][0]["observables"][0]
    observable["observable"] = "cross_section"
    observable["samples"][1].update(integral_pull=4.0, chi2_ndf=16.0, ndf=1)
    check_physics(report)
    observable["samples"][1]["integral_pull"] = 5.1
    with pytest.raises(ValueError, match="integral_pull"):
        check_physics(report)


# Every candidate must have comparison diagnostics, and a rate alone is not a shape check
@pytest.mark.parametrize("tag", ["cross_section", "rap"])
def test_sampling_requires_comparison_diff_obs(tag):
    report = sampling_report()
    observable = report["sets"][0]["observables"][0]
    candidate = dict(observable["samples"][1], label="C_flat")
    candidate.pop("comparison_status")
    observable["samples"].append(candidate)
    with pytest.raises(ValueError, match="MC comparison status is unavailable"):
        check_physics(report)
    candidate["comparison_status"] = "compared"
    observable["observable"] = tag
    candidate["mc_valid_bin_count"] = 1
    observable["samples"][1]["mc_valid_bin_count"] = 1
    with pytest.raises(ValueError, match="no finite differential comparison"):
        check_physics(report)


# Empty report collections cannot validate Monte Carlo precision
def test_mc_precision_rejects_missing_hists():
    report = precision_report({"max_mc_rel_uncertainty": 0.30})
    report["sets"][0]["observables"][0]["samples"] = []
    with pytest.raises(ValueError, match="no predictions"):
        check_physics(report)


# Require comparison diagnostics for every sample subject to a literature bound
@pytest.mark.parametrize("status", [None, "reference", "unavailable"])
def test_literature_comparison_status(status):
    report = precision_report({"max_chi2_ndf": 3.0})
    sample = report["sets"][0]["observables"][0]["samples"][0]
    if status is not None:
        sample["comparison_status"] = status
    with pytest.raises(ValueError, match="comparison|reference"):
        check_physics(report)


# Reject invalid uncertainties instead of interpreting them as high MC precision
@pytest.mark.parametrize("uncertainty", [-0.1, float("nan"), float("inf")])
def test_mc_precision_invalid_error(uncertainty):
    report = precision_report({"max_mc_rel_uncertainty": 0.3}, mc_max_rel_uncertainty=uncertainty)
    with pytest.raises(ValueError, match="mc_max_rel_uncertainty"):
        check_physics(report)


# Check conserving serialized events and reject absent or empty event input
def test_event_kinematics_requires_decoded_events(tmp_path):
    report = precision_report({})
    sample = report["sets"][0]["observables"][0]["samples"][0]
    with pytest.raises(AssertionError, match="no event input files"):
        assert_event_kinematics(report)
    path = tmp_path / "muons.hepmc3"
    sample["input_files"] = [str(path)]
    with pytest.raises(AssertionError):
        assert_event_kinematics(report)
    write_muon_events(path)
    assert_event_kinematics(report)
    sample["input_files"] = [write_muon_events(tmp_path / "empty.hepmc3", weights=())]
    with pytest.raises(AssertionError):
        assert_event_kinematics(report)

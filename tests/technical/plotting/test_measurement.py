# Statistical acceptance of real histogram comparisons with published measurements
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from core.iceplot import comparison_metrics
from core.io.serialize import load_json_file
from core.stats.hist import hobj
from core.stats.uncertainty import covariance_error_ratio
from core.stats.validation import check_measurements, check_physics, measurement_assessment


# Construct actual histogram predictions with known covariance and acceptance criteria
def comparison(mc=(10.0, 10.0), data=(10.0, 10.0), mc_error=0.1, data_error=1.0, density=False):
    histograms = [hobj(counts=np.asarray(values), errs=np.full(len(values), error),
                       bins=np.arange(len(values) + 1, dtype=float), density=density,
                       density_uncertainty="shape") for values, error in ((mc, mc_error), (data, data_error))]
    metrics = comparison_metrics(*histograms)
    criteria = load_json_file(Path(__file__).resolve().parents[3] / "icepack/_common/SETTINGS.json")
    return {"validation": criteria, "normalization": "unit_density" if density else "cross_section",
            "density_uncertainty": "shape" if density else None,
            "sets": [{"name": "measurement", "observables": [{"observable": "spectrum", "kind": "histogram",
                     "samples": [{"label": "prediction", **metrics}]}]}]}


# Report normalization and shape disagreement without interrupting physics checks
@pytest.mark.parametrize("values,accepted", [((10.0, 10.0), True), ((15.0, 15.0), False), ((5.0, 15.0), False)])
def test_agreement(values, accepted):
    report = comparison(mc=values)
    assessment = measurement_assessment(report)
    assert (assessment["status"] == "passed") is accepted
    check_measurements(report)
    check_physics(report)
    assert report["measurement"] == assessment
    sample = report["sets"][0]["observables"][0]["samples"][0]
    assert (sample["measurement"]["status"] == "passed") is accepted


# Reject statistical dilution even when adding large MC errors raises the p-value
def test_mc_noise_cannot_hide_disagreement():
    report = comparison(mc=(15.0, 15.0), mc_error=20.0)
    assessment = measurement_assessment(report)
    assert assessment["status"] == "failed"
    assert any("mc_data_error_ratio" in item for item in assessment["failures"])


# Resolve correlations when MC noise is small in each bin but large in a shape mode
def test_precision_checks_covariance_modes():
    data = np.array([[1.0, 0.99], [0.99, 1.0]])
    assert covariance_error_ratio(np.eye(2) * 0.01, data) == pytest.approx(1.0)
    angle = 0.3
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    for scale in (1e-180, 1.0, 1e180):
        assert covariance_error_ratio(scale * np.eye(2) * 0.01, scale * rotation @ data @ rotation.T) == pytest.approx(1.0)
    assert covariance_error_ratio(np.zeros((2, 2)), data) == pytest.approx(0.0)
    assert covariance_error_ratio(np.eye(2), np.zeros((2, 2))) is None


# Remove the normalization degree of freedom for actual density histograms
def test_density_has_only_shape_information():
    report = comparison(mc=(20.0, 20.0), density=True)
    assessment = measurement_assessment(report)
    sample = report["sets"][0]["observables"][0]["samples"][0]
    assert assessment["status"] == "passed"
    assert assessment["tests"] == 1
    assert sample["ndf"] == 1
    assert set(sample["measurement"]["pvalues"]) == {"chi2"}


# Reject absent comparisons and unpopulated measured intervals instead of accepting plots
@pytest.mark.parametrize("change", ["empty", "missing", "zero", "uncertainty"])
def test_missing_measurement_information_fails(change):
    report = comparison()
    sample = report["sets"][0]["observables"][0]["samples"][0]
    if change == "empty":
        report["sets"] = []
    elif change == "missing":
        sample.pop("comparison_status")
    elif change == "zero":
        sample["mc_zero_bin_count"] = 1
    else:
        sample["chi2"] = None
    assert measurement_assessment(report)["status"] == "failed"


# Correct for both spectra and rates without assuming that observables are independent
def test_probability_multiple_comparisons():
    report = comparison()
    report["sets"].append(deepcopy(report["sets"][0]))
    assessment = measurement_assessment(report)
    assert assessment["tests"] == 4
    assert assessment["min_pvalue"] == pytest.approx(report["validation"]["measurement"]["alpha"] / 4)


# Test extracted point measurements without inventing an integral over their x coordinates
def test_point_measurement_only_cov_statistic():
    report = comparison()
    report['sets'][0]['observables'][0]['kind'] = 'point'
    result = measurement_assessment(report)
    assert result['tests'] == 1
    assert result['status'] == 'passed'


# Prevent plotting overrides and omitted bins from bypassing the configured measurement
@pytest.mark.parametrize('change', ['normalization', 'mask', 'ndf'])
def test_incomplete_measurement_cannot_pass(change):
    report = comparison(density=change == 'normalization')
    sample = report['sets'][0]['observables'][0]['samples'][0]
    if change == 'normalization':
        report['reference_normalization'] = 'cross_section'
    elif change == 'mask':
        sample['mc_valid_bin_count'] -= 1
    else:
        sample['ndf'] = float('nan')
    assert measurement_assessment(report)['status'] == 'failed'


# Record unavailable numerical fields as failures so the report can still be written
@pytest.mark.parametrize("field", ["ndf", "mc_zero_bin_count", "data_valid_bin_count"])
def test_unavailable_statistics_are_reported(field):
    report = comparison()
    report["sets"][0]["observables"][0]["samples"][0][field] = None
    assert measurement_assessment(report)["status"] == "failed"


# Retain zero measured central values when a populated prediction is compatible with their errors
def test_zero_measured_bin_in_comparison():
    report = comparison(mc=(10.0, 0.1), data=(10.0, 0.0), mc_error=0.01)
    assessment = measurement_assessment(report)
    sample = report['sets'][0]['observables'][0]['samples'][0]
    assert sample['data_valid_bin_count'] == sample['mc_valid_bin_count'] == 2
    assert assessment['status'] == 'passed'


# Explicit data-only derived quantities do not become missing MC comparisons
def test_data_only_set():
    report = comparison()
    report['sets'].append(dict(name='derived', mc=False, observables=[dict(
        observable='suppression', kind='point', samples=[])]))
    assert measurement_assessment(report)['status'] == 'passed'
    report['sets'][-1]['mc'] = True
    assert measurement_assessment(report)['status'] == 'failed'


# Keep prediction panels outside the measured statistical family and retain data checks
def test_prediction_panel_with_measurement():
    report = comparison()
    prediction = deepcopy(report["sets"][0])
    prediction["name"] = "prediction"
    sample = prediction["observables"][0]["samples"][0]
    sample["comparison_status"] = "prediction"
    sample["chi2"] = None
    report["sets"].append(prediction)
    assessment = measurement_assessment(report)
    assert assessment["status"] == "passed"
    assert assessment["tests"] == 2
    assert "measurement" not in sample
    report["sets"][0]["observables"][0]["samples"][0].pop("comparison_status")
    assert measurement_assessment(report)["status"] == "failed"

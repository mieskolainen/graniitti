# Validate the sourced diffractive measurements used by the total cross section study
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
from dataclasses import replace

import numpy as np
import pytest

from tests.physics.studies.total_cross_sections.analysis.analyze import (
    DD_XI_005,
    SD_XI_005,
    CrossSectionPoint,
    Uncertainty,
    chi2_for_coupling,
    dd_measurements,
    fit_g3p,
    scan_columns,
    sd_measurements,
    total_measurements,
)


# Derive confidence intervals for a quadratic DD rate from delta chi2 = 1 and 4
def test_g3p_intervals_follow_chi2_quadratic_rate():
    coupling, uncertainty, reference = 0.16, 0.003, 0.2
    series = replace(dd_measurements()[0], fit_observable=DD_XI_005,
                     points=(CrossSectionPoint(900.0, coupling**2,
                                               (Uncertainty("stat", uncertainty, uncertainty),)),))
    columns = scan_columns()
    scan = np.zeros((1, max(vars(columns).values()) + 1))
    scan[0, columns.sqrts], scan[0, columns.dd] = 900.0, reference**2
    fit = fit_g3p(scan, columns, (series,), reference, quiet=True)
    step = fit.values[1] - fit.values[0]
    assert fit.best == pytest.approx(coupling, abs=step)
    for n, lower, upper in ((1, fit.low_1sigma, fit.up_1sigma),
                            (2, fit.low_2sigma, fit.up_2sigma)):
        assert lower == pytest.approx(math.sqrt(coupling**2 - n * uncertainty), abs=2 * step)
        assert upper == pytest.approx(math.sqrt(coupling**2 + n * uncertainty), abs=2 * step)


# Check exact power-law interpolation and prohibit unmeasured energy extrapolation
def test_scan_interpolation_matches_energy():
    from tests.physics.studies.total_cross_sections.analysis.analyze import matched_prediction

    columns = scan_columns()
    scan = np.zeros((2, 12))
    scan[:, columns.sqrts] = [100.0, 10000.0]
    scan[:, columns.sd] = [1.0, 100.0]
    series = sd_measurements()[0]
    point = replace(series.points[0], sqrts_gev=1000.0)
    assert matched_prediction(series, point, 0.2, scan, columns, 0.2) == pytest.approx((10.0, 0.0))
    with pytest.raises(ValueError, match='inside the ordered scan'):
        matched_prediction(series, replace(point, sqrts_gev=10.0), 0.2, scan, columns, 0.2)


# Propagate independent MC errors through logarithmic energy interpolation
@pytest.mark.parametrize('energy,value,error', [(100.0, 1.0, 0.1), (1000.0, 10.0, np.sqrt(0.5)), (10000.0, 100.0, 10.0)])
def test_scan_mc_error(energy, value, error):
    from tests.physics.studies.total_cross_sections.analysis.analyze import scan_prediction

    scan = np.array([[100.0, 1.0, 0.1], [10000.0, 100.0, 10.0]])
    prediction, uncertainty = scan_prediction(scan, energy, 1, 2)
    assert prediction == pytest.approx(value)
    assert uncertainty == pytest.approx(error)


# Test configured cross sections against real TOTEM data and reject a biased normalization
@pytest.mark.parametrize('scale,accepted', [(1.0, True), (2.0, False)])
def test_study_measurement_criteria(scale, accepted):
    from core.stats.validation import measurement_assessment

    from tests.physics.studies.total_cross_sections.analysis.analyze import measurement_report

    series = total_measurements()[0]
    columns = scan_columns()
    scan = np.zeros((1, 12))
    scan[0, 0] = series.points[0].sqrts_gev
    for point in series.points:
        scan[0, getattr(columns, point.channel)] = point.value_mb * scale
    report = measurement_report(scan, columns, (series,), ())
    assert (measurement_assessment(report)['status'] == 'passed') is accepted


# Include the coupling dependence of the SD and DD integration uncertainty
@pytest.mark.parametrize('observable,channel,rate,error', [(SD_XI_005, 'sd', 2.0, 0.3), (DD_XI_005, 'dd', 4.0, 0.6)])
def test_g3p_mc_error(observable, channel, rate, error):
    point = CrossSectionPoint(900.0, 1.5, (Uncertainty("stat", 0.2, 0.2),))
    series = replace(sd_measurements()[0], fit_observable=observable, points=(point,))
    columns = scan_columns()
    scan = np.zeros((1, 12))
    scan[0, columns.sqrts] = 900.0
    scan[0, getattr(columns, channel)] = rate
    scan[0, getattr(columns, channel + '_error')] = error
    assert chi2_for_coupling(0.1, scan, columns, (series,), 0.2) == pytest.approx(4.0)

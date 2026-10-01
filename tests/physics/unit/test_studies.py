# Check physical coordinates and histogram intervals in manual studies
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import base64
import json
import zlib

import numpy as np
import pytest
from core.stats.hist import hobj

from tests.physics.studies.inference_ampfit.inspect import integrals
from tests.physics.studies.sudakov.analyze import load_array


# Preserve complete mass-bin integrals and reject ambiguous partial-bin intervals
def test_mass_intervals_use_complete_bins():
    histogram = hobj(counts=np.array([2.0, 3.0]), errs=np.array([0.2, 0.3]),
                     bins=np.array([0.0, 1.0, 3.0]), cbins=np.array([0.5, 2.0]))
    row = integrals(histogram, [0.0, 3.0])[0]
    assert row['integral'] == pytest.approx(8.0)
    assert row['error'] == pytest.approx(np.hypot(0.2, 0.6))
    for edges in ([0.5, 3.0], [0.0, 2.0], [3.0, 0.0], [0.0, np.nan]):
        with pytest.raises(ValueError):
            integrals(histogram, edges)


# Decode both coordinate transformations on a non-square grid from file metadata
@pytest.mark.parametrize('log', [(False, False), (False, True), (True, False), (True, True)])
def test_sudakov_physical_axes(tmp_path, log):
    axes = [np.array([2., 4.]), np.array([0.1, 0.2, 0.4])]
    stored = [np.log(axis) if flag else axis for axis, flag in zip(axes, log, strict=True)]
    rows = [[x, y, axes[0][i] * axes[1][j], 0.]
            for i, x in enumerate(stored[0]) for j, y in enumerate(stored[1])]
    path = tmp_path / 'array'
    meta = dict(type='IArray2D', version=1, shape=[2, 3, 4],
                axes=['q2', 'x'], log=log, sqrts=123., pdf='input PDF',
                data=base64.b64encode(zlib.compress(np.asarray(rows, dtype='>f8').tobytes())).decode())
    path.write_text(json.dumps(meta))
    loaded, physical, values = load_array(path)
    assert loaded['sqrts'] == pytest.approx(meta['sqrts'])
    assert loaded['pdf'] == meta['pdf']
    for expected, actual in zip(axes, physical, strict=True):
        np.testing.assert_allclose(actual, expected)
    np.testing.assert_allclose(values, np.outer(axes[1], axes[0]))
    meta['type'] = 'unsupported'
    path.write_text(json.dumps(meta))
    with pytest.raises(ValueError, match='unsupported interpolation cache'):
        load_array(path)


# Keep screened cross sections separate from the appended integration-error columns
def test_screening_scan_uses_cross_sections(tmp_path):
    from tests.physics.studies.cross_section_scan.analysis.analyze import ScanColumns, load_combined_scan

    paths = [tmp_path / name for name in ('bare.csv', 'screened.csv')]
    rates = np.array([[100.0, 1.0, 2.0, 3.0, 10.0, 20.0, 30.0, 0.1, 0.2, 0.3]])
    for path, factor in zip(paths, (1.0, 0.5), strict=True):
        table = rates.copy()
        table[:, 4:7] *= factor
        np.savetxt(path, table, delimiter='\t', header='scan')
    scan = load_combined_scan(*paths)
    columns = ScanColumns()
    for bare, screened in ((columns.xs_el, columns.xs_el_screened),
                           (columns.xs_sd, columns.xs_sd_screened), (columns.xs_dd, columns.xs_dd_screened)):
        assert scan[0, screened] / scan[0, bare] == pytest.approx(0.5)

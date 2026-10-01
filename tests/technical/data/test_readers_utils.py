# Tests for the standalone readers utility functions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
import pytest
from core.analysis import obs as observables
from core.io import readers
from core.io.hepdata_reader import load_reader
from core.io.readers import generate_chunks, generate_partitions
from core.stats.hist import center2edgebins, rebin_histogram, rebin_histogram_groups
from pyHepMC3 import HepMC3 as h

from tests.technical.support.hepmc import write_muon_events


# Validate inclusive partition boundaries without overlap or lost events
@pytest.mark.parametrize(('total', 'size', 'expected'), [
    (5, 10, [(0, 4)]), (10, 10, [(0, 9)]),
    (25, 10, [(0, 9), (10, 19), (20, 24)]),
    (23, 10, [(0, 9), (10, 19), (20, 22)]),
    (100, 50, [(0, 49), (50, 99)]), (7, 10, [(0, 6)]), (8, 20, [(0, 7)]),
])
def test_generate_partitions(total, size, expected):
    assert list(generate_partitions(total, size)) == expected


# Reject empty input and zero partition sizes
@pytest.mark.parametrize(('total', 'size'), [(0, 10), (10, 0)])
def test_generate_partitions_invalid(total, size):
    with pytest.raises(ValueError):
        list(generate_partitions(total, size))


# Construct real HepMC input and two distinct invariant-mass histogram selections
@pytest.fixture
def worker(tmp_path):
    observable = readers.get_observables('default')['M']
    return {
        'chunk_range': [0, 2],
        'hepmc3file': write_muon_events(tmp_path / 'muons.hepmc3', energies=[1.0, 2.0, 3.0]),
        'obs': [{'M': {**observable, 'bins': np.array(bins)}} for bins in ([1., 3., 5., 7.], [1., 5., 7.])],
        'pid': [[13, -13]] * 2, 'cuts': ['core.analysis.cuts.default'] * 2,
        'xsmode': 'header', 'header_xsection': None, 'scales': [2.0, 3.0],
        'k': 1, 'label': 'MC', 'density': True, 'density_uncertainty': 'scaled', 'verbose': False,
    }


# Check the complete reader-to-histogram path against Poisson and multinomial covariance
@pytest.mark.parametrize('density', [False, True])
@pytest.mark.parametrize('density_uncertainty', ['scaled', 'shape'])
def test_parallel_wrapper_hist_cov(worker, density, density_uncertainty):
    worker.update(density=density, density_uncertainty=density_uncertainty)
    result = readers.parallel_wrapper(worker)
    for item, counts, scale in zip(result, ([1., 1., 1.], [2., 1.]), worker['scales'], strict=True):
        histogram = item['M']['hdata']
        counts = np.array(counts)
        widths = np.diff(histogram.bins)
        factor = (1.0 / 3.0 if density else scale) / widths
        covariance = np.diag(counts)
        if density and density_uncertainty == 'shape':
            covariance -= np.outer(counts, counts) / 3.0
        np.testing.assert_allclose(histogram.counts_scaled, counts * factor, rtol=1e-12)
        np.testing.assert_allclose(histogram.covariance_scaled, covariance * np.outer(factor, factor), atol=1e-15)
        assert histogram.integral() == pytest.approx(1.0 if density else 3.0 * scale)


# Reject a missing per-set scale after reading the same real event selections
def test_parallel_wrapper_requires_scale_per(worker):
    worker['scales'] = [1.0]
    with pytest.raises(ValueError, match='one value per dataset set'):
        readers.parallel_wrapper(worker)


# Validate sequence chunking including an empty sequence and a partial final chunk
@pytest.mark.parametrize(('items', 'size', 'expected'), [
    ([1, 2, 3, 4, 5, 6], 2, [[1, 2], [3, 4], [5, 6]]),
    ([1, 2, 3], 1, [[1], [2], [3]]), ([1, 2, 3], 5, [[1, 2, 3]]),
    ([1, 2, 3], 3, [[1, 2, 3]]), ([], 3, []),
])
def test_generate_chunks(items, size, expected):
    assert list(generate_chunks(items, size)) == expected


# Reject zero, negative and noninteger chunk sizes
@pytest.mark.parametrize(('size', 'error'), [(0, ValueError), (-1, ValueError), ('2', TypeError)])
def test_generate_chunks_invalid(size, error):
    with pytest.raises(error):
        list(generate_chunks([1, 2, 3], size))


# Check integral and differential rebinning against independent Poisson propagation
@pytest.mark.parametrize(('edges', 'values', 'differential', 'expected_edges', 'expected', 'variance'), [
    ([0, 1, 2, 3, 4, 5, 6], [100, 200, 300, 400, 500, 600], False,
     [0, 2, 4, 6], [300, 700, 1100], [300, 700, 1100]),
    ([0, 1, 2, 3, 4, 5, 6], [100, 200, 300, 400, 500, 600], True,
     [0, 2, 4, 6], [150, 350, 550], [75, 175, 275]),
    ([0, 2, 5, 9], [100, 200, 300], False, [0, 5, 9], [300, 300], [300, 300]),
    ([0, 2, 5, 9], [100, 200, 300], True, [0, 5, 9], [160, 300], [88, 300]),
])
def test_rebin_histogram(edges, values, differential, expected_edges, expected, variance):
    edges, values = np.array(edges), np.array(values)
    new_edges, new_values, new_errors = rebin_histogram(edges, values, np.sqrt(values), 2, differential=differential)
    np.testing.assert_allclose(new_edges, expected_edges)
    np.testing.assert_allclose(new_values, expected)
    np.testing.assert_allclose(new_errors, np.sqrt(variance))
    assert np.sum(new_values * (np.diff(new_edges) if differential else 1)) == pytest.approx(
        np.sum(values * (np.diff(edges) if differential else 1)))


# Check explicit groups preserve differential integrals over non-uniform bins
def test_rebin_histogram_groups():
    edges = np.asarray([0.0, 1.0, 3.0, 6.0, 10.0])
    values = np.asarray([2.0, 4.0, 6.0, 8.0])
    errors = np.asarray([0.2, 0.4, 0.6, 0.8])
    groups = np.asarray([0, 1, 3, 4])

    new_edges, new_values, new_errors = rebin_histogram_groups(
        edges,
        values,
        errors,
        groups,
        differential=True,
    )

    np.testing.assert_array_equal(new_edges, [0.0, 1.0, 6.0, 10.0])
    np.testing.assert_allclose(new_values, [2.0, 5.2, 8.0])
    np.testing.assert_allclose(
        new_errors,
        [0.2, np.hypot(0.4 * 2.0, 0.6 * 3.0) / 5.0, 0.8],
    )
    assert np.sum(new_values * np.diff(new_edges)) == pytest.approx(
        np.sum(values * np.diff(edges))
    )


# Check grouped rebinning rejects malformed physical inputs and partitions
@pytest.mark.parametrize(
    ("edges", "values", "errors", "groups", "error"),
    [
        ([0.0, 1.0, 2.0], [1.0, 2.0], [0.1, 0.2], [0, 1], ValueError),
        ([0.0, 1.0, 2.0], [1.0, 2.0], [0.1, 0.2], [0, 2, 1, 2], ValueError),
        ([0.0, 1.0, 2.0], [1.0, 2.0], [0.1, 0.2], [0, 1.5, 2], ValueError),
        ([0.0, 0.0, 2.0], [1.0, 2.0], [0.1, 0.2], [0, 2], ValueError),
        ([0.0, 1.0, 2.0], [1.0, np.nan], [0.1, 0.2], [0, 2], ValueError),
        ([0.0, 1.0, 2.0], [1.0, 2.0], [0.1, -0.2], [0, 2], ValueError),
    ],
)
def test_rebin_hist_groups_invalid_inputs(edges, values, errors, groups, error):
    with pytest.raises(error):
        rebin_histogram_groups(edges, values, errors, groups)


# Check grouped rebinning requires an explicit Boolean differential policy
def test_rebin_hist_groups_invalid_diff_policy():
    with pytest.raises(TypeError, match="differential must be boolean"):
        rebin_histogram_groups([0.0, 1.0], [1.0], [0.1], [0, 1], differential="true")


# Check fixed-factor rebinning requires a positive integer factor
@pytest.mark.parametrize("factor", [True, 0, -1, 1.5])
def test_rebin_hist_invalid_factor(factor):
    error = TypeError if isinstance(factor, (bool, float)) else ValueError
    with pytest.raises(error):
        rebin_histogram(
            np.asarray([0.0, 1.0]),
            np.asarray([1.0]),
            np.asarray([0.1]),
            factor,
        )


# Reject inconsistent histogram edge and content lengths
def test_rebin_histogram_invalid_input():
    bin_edges = np.array([0, 1, 2, 3])
    bin_contents = np.array([10, 20])

    with pytest.raises(ValueError):
        rebin_histogram(bin_edges, bin_contents, np.sqrt(bin_contents), 2)


# Validate center-to-edge conversion for uniform, irregular and empty bins
@pytest.mark.parametrize(('centers', 'edges'), [
    ([1, 2, 3, 4], [.5, 1.5, 2.5, 3.5, 4.5]),
    ([1, 3, 6, 10], [0, 2, 4.5, 8, 12]), ([2, 4], [1, 3, 5]), ([], [0]),
])
def test_center2edgebins(centers, edges):
    np.testing.assert_allclose(center2edgebins(np.array(centers)), edges, atol=1e-6)


# Keep independent histogram metadata when reading the same cached publication twice
def test_read_hepdata_isolates_cached_metadata():
    root = Path(__file__).resolve().parents[3]
    dataset_path = root / 'icepack/GAMMA/ATLAS_1377585/mumu/dataset.json'
    table = root / 'HEPData/CEP/HEPData-ins1377585-v1-json/Table5.json'
    reference = '../_common/reader.py'
    reader = load_reader(reference, dataset_path=str(dataset_path), cdir=str(root))
    shared = reader.read(str(table))
    dataset = {'hist': [{'file': str(table), 'obs': name, 'scale': scale, 'fitw': fitw}
                       for name, scale, fitw in [('first', 0.001, 1.0), ('second', 0.002, 5.0)]]}
    data, obs = readers.read_hepdata(dataset, str(table.parent), 'HEPDATA', {'first': {}, 'second': {}},
                                 cdir=str(root), reader=reference, dataset_path=str(dataset_path))
    assert [data[name]['scale'] for name in data] == pytest.approx([0.001, 0.002])
    assert [data[name]['fitw'] for name in data] == pytest.approx([1.0, 5.0])
    data['first']['y'][0] = -1.0
    np.testing.assert_array_equal(data['second']['y'], shared['y'])
    np.testing.assert_array_equal(reader.read(str(table))['y'], shared['y'])
    assert 'scale' not in shared
    assert not np.shares_memory(obs['first']['bins'], shared['bins'])


# Read equivalent external HepMC events using their declared momentum units
@pytest.mark.parametrize('unit', [h.Units.GEV, h.Units.MEV])
def test_hepmc_momentum_units_before_projection(tmp_path, unit):
    filename = write_muon_events(tmp_path / 'units.hepmc3', weights=[1.0], unit=unit)
    output = readers.read_hepmc3(
        hepmc3file=filename, obs=[{'M': {'func': observables.proj_1D_M}}],
        pid=[[13, -13]], cuts=['core.analysis.cuts.default'],
    )[0]
    np.testing.assert_allclose(output['data']['M'], [2.0])
    assert output['event_acceptance'] == pytest.approx(1.0)

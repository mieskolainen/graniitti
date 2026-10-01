# Test assigned MC samples through the native iceplot command line
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json

import pytest

from tests.technical.support import iceplot
from tests.technical.support.hepmc import write_muon_events
from tests.technical.support.output import PhysicsOutput


# Preserve distinct sample normalization, cuts and diagnostics in assigned MC-only panels
def test_assigned_samples(tmp_path):
    root = iceplot.CDIR
    paths = [tmp_path / f'{name}.hepmc3' for name in ('a', 'b')]
    write_muon_events(paths[0], energies=[0.8, 1.0, 1.2])
    write_muon_events(paths[1], weights=[1.0, 1.0], energies=[1.2, 1.4])
    card = {
        'active': True, 'type': 'MC_ONLY',
        'samples': [{'name': name, 'label': name, 'parameters': {},
                     'gencard': str(root / 'icepack/GAMMA/ATLAS_1377585/mumu/gencard.json')}
                    for name in ('a', 'b')],
        'plot': {'normalization': 'cross_section', 'ratio_uncertainty': 'combined',
                 'stack': False, 'data_style': 'hist'},
        'fit': {'normalization': 'cross_section'},
        'sets': [{'name': name, 'samples': [name], 'data': False, 'pid': [13, -13],
                  'cuts': str(root / 'icepack/_common/cuts_inclusive.py'),
                  'obs': str(root / 'python/src/core/analysis/observables/default.py'),
                  'hist': [{'obs': 'M', 'xmin': 0.0, 'xmax': 4.0, 'nbins': 8}]}
                 for name in ('a', 'b')],
    }
    dataset = tmp_path / 'dataset.json'
    dataset.write_text(json.dumps(card))
    output = PhysicsOutput(f'assigned_samples_{tmp_path.name}')
    iceplot.run_iceplot(
        hepmc3_tags=[str(path.with_suffix('')) for path in paths], plot_tag='assigned',
        labels=['a', 'b'], analysis=str(tmp_path), mc_hepdata_samples=['a', 'b'],
        mc_scales=[1.0, 2.0], unit='pb', output=output,
    )
    report = iceplot.load_iceplot_report(output.report('assigned'))
    for index, selected in enumerate(report['sets']):
        sample, = selected['observables'][0]['samples']
        assert sample['label'] == ('a', 'b')[index]
        assert sample['integral'] == pytest.approx((3.0, 6.0)[index])
        assert sample['mc_selected_events'] == (3, 2)[index]
        assert sample['mc_effective_events'] == pytest.approx((3.0, 2.0)[index])

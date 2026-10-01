# Unit tests for the shared iceplot test launcher
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import pickle

import numpy as np
import pytest

from tests.technical.support import iceplot
from tests.technical.support.hepmc import write_muon_events
from tests.technical.support.output import PhysicsOutput


# Run the actual launcher with automatic batching and cross-section mode selection
def test_run_iceplot_uses_automatic_batching(tmp_path):
    path = tmp_path / 'muons.hepmc3'
    write_muon_events(path, energies=[0.8, 1.0, 1.2])
    observables = tmp_path / 'mass.py'
    observables.write_text(
        '# Invariant mass projection for real event analysis\n#\n'
        '# (c) 2026 Mikael Mieskolainen\n'
        '# Licensed under the MIT License <http://opensource.org/licenses/MIT>.\n'
        'from core.analysis.observables.default import obs_M\n'
    )
    output = PhysicsOutput(f'automatic_batching_{tmp_path.name}')
    iceplot.run_iceplot(
        hepmc3_tags=[str(path.with_suffix(''))], plot_tag='muons', labels=['muons'],
        obs_module=str(observables), cuts=['core.analysis.cuts.default'], pid=[[13, -13]], unit='pb', output=output,
    )
    report = iceplot.load_iceplot_report(output.report('muons'))
    sample = report['sets'][0]['observables'][0]['samples'][0]
    assert sample['integral'] == pytest.approx(3.0)
    assert sample['mc_selected_events'] == 3
    assert sample['mc_effective_events'] == pytest.approx(3.0)


# Require worker histogram results to cross multiprocessing boundaries
def test_histogram_pickle():
    all_observables = iceplot.readers.get_observables("default")
    observables = {"M": all_observables["M"].copy()}
    result = iceplot.plots.histmc(
        mcdata={
            "data": {"M": np.asarray([0.5, 1.0, 1.5])},
            "weights": np.ones(3),
            "xsection_pb": 1.0,
        },
        obs=observables,
    )

    restored = pickle.loads(pickle.dumps(result))
    np.testing.assert_array_equal(restored["M"]["hdata"].counts, result["M"]["hdata"].counts)
    np.testing.assert_array_equal(restored["M"]["hdata"].errs, result["M"]["hdata"].errs)

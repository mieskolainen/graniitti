# Check absolute normalization and sample reuse through the real HepMC3 and Pythia APIs
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import pytest
from pyHepMC3 import HepMC3 as h
from pyHepMC3 import std

from tests.external.pythia.cases.z_merging import study
from tests.technical.drivers.test_pythia_lhe import DRIVER, write_lhe
from tests.technical.support.hepmc import read_events, write_muon_events


# Add LHE source weights to physical dimuon events using the official binding
def sample(path, weights, sigma, error):
    source = path.with_suffix('.source.hepmc3')
    write_muon_events(source, weights, xsections=[(sigma, error)] * len(weights))
    writer = h.WriterAscii(str(path))
    try:
        for event in read_events(source):
            event.add_attribute('graniitti_lhe_weight', h.DoubleAttribute(event.weights()[0]))
            writer.write_event(event)
    finally:
        writer.close()
    return path


# Add different multiplicities with unequal weights and retain physical zero-weight vetoes
@pytest.mark.parametrize('zero', [False, True])
def test_sum_multiplicities(tmp_path, zero):
    first = sample(tmp_path / 'first.hepmc3', [2.0, 1.0], 3.0, 0.2)
    second = sample(tmp_path / 'second.hepmc3', [0.0, 0.0] if zero else [0.0, 3.0],
                    0.0 if zero else 5.0, 0.0 if zero else 0.3)
    output = tmp_path / 'sum.hepmc3'
    expected = 3.0 if zero else 8.0
    assert study.combine([first, second], output) == pytest.approx(expected)
    events = read_events(output)
    assert len(events) == 4
    assert sum(event.weights()[0] for event in events) == pytest.approx(expected)
    assert events[-1].cross_section().xsec() == pytest.approx(expected)
    assert events[-1].cross_section().xsec_err() == pytest.approx(0.2 if zero else math.hypot(0.2, 0.3))
    assert events[2].weights()[0] == pytest.approx(0)


# Reject duplicate samples and output aliases before modifying any input
@pytest.mark.parametrize('alias', ['same', 'symlink', 'hardlink', 'duplicate'])
def test_input_aliases(tmp_path, alias):
    source = sample(tmp_path / 'first.hepmc3', [1.0], 3.0, 0.2)
    output = tmp_path / 'out.hepmc3'
    if alias == 'same':
        output = source
    elif alias == 'symlink':
        output.symlink_to(source)
    elif alias == 'hardlink':
        output.hardlink_to(source)
    before = source.read_bytes()
    with pytest.raises(ValueError, match='different file'):
        study.combine([source, source] if alias == 'duplicate' else [source], output)
    assert source.read_bytes() == before


# Reject damaged event graphs and invalid normalization without overwriting an existing output
@pytest.mark.parametrize('damage', ['empty', 'truncated', 'zero-normalization', 'nonfinite'])
def test_invalid_sample(tmp_path, damage):
    source = sample(tmp_path / 'source.hepmc3', [1.0, 1.0], 3.0, 0.2)
    if damage == 'empty':
        sample(source, [], 0.0, 0.0)
    elif damage == 'zero-normalization':
        sample(source, [0.0, 0.0], 3.0, 0.2)
    elif damage == 'nonfinite':
        write_muon_events(source, [float('nan')])
    else:
        text = source.read_text()
        source.write_text(text[:text.rfind('\nP ')])
    output = tmp_path / 'out.hepmc3'
    output.write_text('preserve this output')
    with pytest.raises(ValueError):
        study.combine([source], output)
    assert output.read_text() == 'preserve this output'


# Reuse complete samples but rerun when steering or output contents change
@pytest.mark.skipif(not DRIVER.is_file(), reason='Build the Pythia converter first')
def test_reuse_inputs(tmp_path, monkeypatch):
    monkeypatch.setattr(study, 'WORK', tmp_path)
    source, output, card = (tmp_path / name for name in ('source.lhe', 'out.hepmc3', 'shower.cmnd'))
    write_lhe(source, 'IPIP')
    card.write_text('PartonLevel:ISR = off\n')
    command = [DRIVER, source, output, '2', '58031', card, '5', 'fragment']
    options = dict(output=output, inputs=[DRIVER, source, card], reuse=True)
    study.run(command, 'convert', **options)
    stamp = output.stat().st_mtime_ns
    study.run(command, 'convert', **options)
    assert output.stat().st_mtime_ns == stamp
    card.write_text('PartonLevel:ISR = off\nPartonLevel:FSR = off\n')
    study.run(command, 'convert', **options)
    assert output.stat().st_mtime_ns != stamp
    output.write_text('incomplete output')
    study.run(command, 'convert', **options)
    assert len(read_events(output)) == 2


# A truncated merged sample must fail the source normalization check
def test_merging_check(tmp_path):
    path = sample(tmp_path / 'source.hepmc3', [1.0], 1.0, 0.1)
    event = read_events(path)[0]
    event.add_attribute('graniitti_merging_source_sum', h.DoubleAttribute(2.0))
    event.add_attribute('graniitti_merging_source_xsec', h.DoubleAttribute(2.0))
    event.add_attribute('graniitti_merging_weight', h.DoubleAttribute(1.0))
    event.add_attribute('graniitti_closure_pass', h.IntAttribute(1))
    output = tmp_path / 'partial.hepmc3'
    writer = h.WriterAscii(str(output))
    writer.write_event(event)
    writer.close()
    with pytest.raises(ValueError, match='incomplete source sample'):
        study.check(output)
    data = h.GenEventData()
    event.write_data(data)
    data.weights = std.vector_double([0.5])
    event.read_data(data)
    output = tmp_path / 'wrong.hepmc3'
    writer = h.WriterAscii(str(output))
    writer.write_event(event)
    writer.close()
    with pytest.raises(ValueError, match='inconsistent merging weight'):
        study.check(output)

# Check serialized shower conservation and reject inconsistent event records
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import pytest
from pyHepMC3 import HepMC3 as h

from tests.external.pythia.drivers.lhe_converter.check_matching import (
    read_hepmc_summary,
    read_lhe_events,
    validate_pair,
)
from tests.physics.unit.test_icepack_precision import precision_report
from tests.physics.validation.test_icepacks import assert_event_kinematics
from tests.technical.drivers.test_pythia_lhe import DRIVER, run_converter, write_lhe
from tests.technical.support.hepmc import read_events

pytestmark = pytest.mark.skipif(not DRIVER.is_file(), reason="Build the Pythia converter first")


# Generate real converter events with and without a requested lepton shower
@pytest.fixture(scope="module")
def converted(tmp_path_factory, request):
    folder = tmp_path_factory.mktemp("closure")
    source, output, card = (folder / name for name in ("source.lhe", "out.hepmc3", "shower.cmnd"))
    write_lhe(source, "direct")
    card.write_text(f"PartonLevel:ISR = off\nTimeShower:QEDshowerByL = {'on' if request.param else 'off'}\n")
    result = run_converter(source, output, card, converter_mode="fragment")
    assert result.returncode == 0, result.stdout + result.stderr
    return source, output, request.param


# Serialize the modified event graph through the official HepMC3 writer
def write_events(path, events, unit=h.Units.GEV):
    writer = h.WriterAscii(str(path))
    writer.set_precision(17)
    try:
        for event in events:
            event.set_units(unit, h.Units.MM)
            writer.write_event(event)
            assert not writer.failed()
    finally:
        writer.close()


# Recompute conservation and apply shower accuracy only to showered events
@pytest.mark.parametrize("converted", [False, True], indirect=True, ids=["direct", "shower"])
@pytest.mark.parametrize("unit", [h.Units.GEV, h.Units.MEV], ids=["GeV", "MeV"])
@pytest.mark.parametrize("fraction", [0.0, 0.5, 2.0])
def test_momentum_closure(tmp_path, converted, fraction, unit):
    source, original, shower = converted
    events = read_events(original)
    tolerance = float(events[0].attribute_as_string("graniitti_closure_tolerance"))
    delta = fraction * tolerance
    particle = next(p for p in events[0].particles() if p.pid() == 2212 and not p.end_vertex())
    p = particle.momentum()
    particle.set_momentum(h.FourVector(p.px(), p.py(), p.pz(), p.e() + delta))
    output = tmp_path / "changed.hepmc3"
    write_events(output, events, unit)
    report = precision_report({}, input_files=[str(output)])
    if fraction > 1.0 or (fraction > 0.0 and not shower):
        with pytest.raises(AssertionError):
            assert_event_kinematics(report)
        with pytest.raises(RuntimeError, match="closure too large"):
            validate_pair(source, output, len(events))
    else:
        assert_event_kinematics(report)
        validate_pair(source, output, len(events))
        summary = read_hepmc_summary(output, len(events), read_lhe_events(source, len(events)))
        assert summary["closure_max"] == pytest.approx(delta, abs=1e-10)


# Reject missing, nonfinite or failed shower diagnostics even for conserved momenta
@pytest.mark.parametrize("converted", [True], indirect=True, ids=["shower"])
@pytest.mark.parametrize("attribute,value", [
    ("tolerance", None), ("tolerance", "nan"), ("tolerance", "inf"), ("tolerance", "0"),
    ("tolerance", "-1"), ("max_abs", "nan"), ("max_abs", "-1"), ("max_abs", "1e30"), ("pass", "0"),
])
def test_invalid_shower_closure(tmp_path, converted, attribute, value):
    source, original, _ = converted
    events = read_events(original)
    name = f"graniitti_closure_{attribute}"
    events[0].remove_attribute(name)
    if value is not None:
        events[0].add_attribute(name, h.StringAttribute(value))
    output = tmp_path / "invalid.hepmc3"
    write_events(output, events)
    with pytest.raises((AssertionError, ValueError)):
        assert_event_kinematics(precision_report({}, input_files=[str(output)]))
    with pytest.raises(RuntimeError):
        validate_pair(source, output, len(events))

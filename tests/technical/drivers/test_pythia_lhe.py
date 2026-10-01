# Physics and input validation tests for the Pythia LHE converter
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import importlib.util
import math
import subprocess
from pathlib import Path

import pytest
from core.io import readers
from pyHepMC3 import HepMC3 as hepmc

from tests.technical.support.hepmc import read_events as read_hepmc_events

ROOT = Path(__file__).resolve().parents[3]
DRIVER = ROOT / "bin/pythia_lhe_hadronize"
pytestmark = pytest.mark.skipif(not DRIVER.is_file(), reason="Build the Pythia converter first")


# Retain merging vetoes and normalize the complete weighted source through the real shower
@pytest.mark.parametrize("resonance", [False, True])
def test_ckkwl_weights_and_cross_section(tmp_path, resonance):
    source, output, card = (tmp_path / name for name in ('source.lhe', 'merged.hepmc3', 'merging.cmnd'))
    write_lhe(source, 'full', count=32, resonance=resonance)
    size = 5 if resonance else 4
    source.write_text(source.read_text().replace(f'\n{size} 1 1 10 ', f'\n{size} 1 2 10 ', 16))
    card.write_text('PartonLevel:ISR = on\nMerging:doKTMerging = on\n'
                    'Merging:Process = pp>mu+mu-\nMerging:nJetMax = 1\n'
                    'Merging:nRequested = 0\nMerging:TMS = 1\n')
    result = run_converter(source, output, card, count=32, converter_mode='shower')
    assert result.returncode == 0, result.stdout + result.stderr
    events = read_hepmc_events(output)
    assert len(events) == 32
    raw, weights, factors = [], [], []
    for event in events:
        raw.append(float(event.attribute_as_string('graniitti_lhe_weight')))
        factors.append(float(event.attribute_as_string('graniitti_merging_weight')))
        weights.append(event.weights()[0])
        assert weights[-1] == pytest.approx(raw[-1] * factors[-1])
        assert float(event.attribute_as_string('graniitti_closure_pass')) > 0
    assert any(f <= 0 for f in factors) and any(f > 0 for f in factors)
    assert sum(raw) == pytest.approx(48)
    assert events[-1].cross_section().xsec() == pytest.approx(sum(weights) / sum(raw))
    summary = readers.read_hepmc3_weight_summary(str(output))
    assert readers.resolve_hepmc3_xsmode('auto', summary)[0] == 'header'
    assert readers.terminal_header_xsection(summary)[0] == pytest.approx(sum(weights) / sum(raw))
    ratio = sum(weights) / sum(raw)
    variance = len(raw) / (len(raw) - 1) * sum(
        (merged - ratio * source) ** 2 for source, merged in zip(raw, weights, strict=True)
    ) / sum(raw) ** 2
    assert events[-1].cross_section().xsec_err() == pytest.approx(math.sqrt(variance))
    assert events[-1].cross_section().xsec_err() > 0
    from tests.external.pythia.cases.z_merging.study import check
    assert check(output)["events"] == len(events)


# Keep a source event outside native remnant support as an explicitly vetoed merging record
def test_ckkwl_remnant_support_veto(tmp_path):
    source = Path(__file__).with_name('data') / 'beam_remnant_limit.lhe.txt'
    output, card = tmp_path / 'merged.hepmc3', tmp_path / 'merging.cmnd'
    card.write_text((ROOT / 'icepack/HARDPOM/z_merging/merging.cmnd').read_text())
    result = run_converter(source, output, card, count=1, converter_mode='shower')
    assert result.returncode == 0, result.stdout + result.stderr
    events = read_hepmc_events(output)
    assert len(events) == 1
    event = events[0]
    assert event.attribute_as_string('graniitti_merging_remnant_veto') == '1'
    assert event.weights()[0] == pytest.approx(0)
    assert event.cross_section().xsec() == pytest.approx(0)
    assert float(event.attribute_as_string('graniitti_lhe_weight')) > 0
    assert float(event.attribute_as_string('graniitti_closure_pass')) > 0


# Format a physical LHE particle with its generated invariant mass
def particle_row(pdg, status, px, py, pz, energy, colors=(0, 0)):
    mass = math.sqrt(max(0.0, energy * energy - px * px - py * py - pz * pz))
    mother = (0, 0) if status == -1 else (1, 2)
    return (
        f"{pdg} {status} {mother[0]} {mother[1]} {colors[0]} {colors[1]} "
        f"{px:.17g} {py:.17g} {pz:.17g} {energy:.17g} {mass:.17g} 0 9"
    )


# Write complete pp states or a reduced Pomeron hard process with exact source counters
def write_lhe(path, mode, count=2, flavour=21, jets=False, resonance=False):
    proton_mass = 0.93827208943
    beam_pz = math.sqrt(100.0**2 - proton_mass**2)
    leading_pz = math.sqrt(70.0**2 - proton_mass**2)
    full = mode in {"direct", "IPIP", "IPp"}
    rows = (
        [particle_row(2212, -1, 0, 0, beam_pz, 100), particle_row(2212, -1, 0, 0, -beam_pz, 100)]
        if full
        else [particle_row(2, -1, 0, 0, 20, 20, (501, 0)), particle_row(-2, -1, 0, 0, -20, 20, (0, 501))]
    )
    if full:
        rows.append(particle_row(2212, 1, 0, 0, leading_pz, 70))
    if mode in {"direct", "IPIP"}:
        rows.append(particle_row(2212, 1, 0, 0, -leading_pz, 70))
    if mode == "IPIP":
        rows += [particle_row(21, 1, 0, 0, 10, 10, (501, 502)), particle_row(21, 1, 0, 0, -10, 10, (502, 501))]
    elif mode == "IPp":
        # Supply standard Pythia quark, diquark and hadron identities for ordinary remnants
        states = {
            21: (2, 2101, 0.00216, 2 * proton_mass / 3, (501, 0), (0, 502), (502, 501)),
            2: (2101, 21, 2 * proton_mass / 3, 0.0, (0, 502), (502, 501), (501, 0)),
            3: (321, 2101, 0.493677, 2 * proton_mass / 3, (0, 0), (0, 501), (501, 0)),
            -5: (5122, 2, 5.61957, 0.00216, (0, 0), (501, 0), (0, 501)),
        }
        id1, id2, mq, md, c1, c2, spectator_color = states[flavour]
        mass, energy = 6.0, 70.0
        recoil_pz = -math.sqrt(energy * energy - mass * mass)
        quark_energy = (mass * mass + mq * mq - md * md) / (2 * mass)
        momentum = math.sqrt(quark_energy * quark_energy - mq * mq)
        rows += [
            particle_row(id1, 1, momentum, 0, recoil_pz * quark_energy / mass, energy * quark_energy / mass, c1),
            particle_row(
                id2,
                1,
                -momentum,
                0,
                recoil_pz * (mass - quark_energy) / mass,
                energy * (mass - quark_energy) / mass,
                c2,
            ),
            particle_row(flavour, 1, 0, 0, 10, 10, spectator_color),
        ]
        hard_pz = -leading_pz - recoil_pz - 10
        hard_mass = math.sqrt(50**2 - hard_pz**2)
        rows += [
            particle_row(-13, 1, hard_mass / 2, 0, hard_pz / 2, 25),
            particle_row(13, 1, -hard_mass / 2, 0, hard_pz / 2, 25),
        ]
    muon_energy = 30 if mode == "direct" else 20
    if mode != "IPp":
        rows += [
            particle_row(-13, 1, muon_energy, 0, 0, muon_energy),
            particle_row(13, 1, -muon_energy, 0, 0, muon_energy),
        ]
    side_mask = 1 if mode == "IPp" else 3
    xi = 1 - (70 + leading_pz) / (100 + beam_pz)
    side1 = f"side1_active=1 side1_lead_pdg=2212 side1_lead_pz={leading_pz:.17g} side1_lead_e=70 side1_lead_m={proton_mass} side1_xi={xi:.17g} side1_beta=0.6666666667"
    side2 = f"side2_active=1 side2_lead_pdg=2212 side2_lead_pz={-leading_pz:.17g} side2_lead_e=70 side2_lead_m={proton_mass} side2_xi={xi:.17g} side2_beta=0.6666666667"
    if mode == "full":
        side1 += " side1_rem_pdg=-2 side1_rem_pz=10 side1_rem_e=10"
        side2 += " side2_rem_pdg=2 side2_rem_pz=-10 side2_rem_e=10"
    if mode == "IPp":
        side2 = "side2_active=0"
    hard_rows = []
    if mode in {"IPp", "IPIP"}:
        # Cross the supplied remnant colors into the independent incoming hard partons
        if mode == "IPIP":
            ids, colors = (21, 21), ((502, 501), (501, 502))
        else:
            ids = (21 if flavour == 21 else -flavour, flavour)
            colors = (spectator_color[::-1], spectator_color)
        hard_e = sum(float(row.split()[9]) for row in rows[-2:])
        hard_pz = sum(float(row.split()[8]) for row in rows[-2:])
        hard_rows = [
            particle_row(ids[0], -1, 0, 0, (hard_e + hard_pz) / 2, (hard_e + hard_pz) / 2, colors[0]),
            particle_row(ids[1], -1, 0, 0, -(hard_e - hard_pz) / 2, (hard_e - hard_pz) / 2, colors[1]),
            *rows[-2:],
        ]
    if jets:
        assert mode == "IPIP"
        rows[5] = particle_row(21, 1, 0, 0, -10, 10, (503, 501))
        rows[-2] = particle_row(21, 1, 20, 0, 0, 20, (502, 504))
        rows[-1] = particle_row(21, 1, -20, 0, 0, 20, (504, 503))
        hard_rows[1] = particle_row(21, -1, 0, 0, -20, 20, (501, 503))
        hard_rows[-2:] = rows[-2:]
    if resonance:
        # Keep explicit Z ancestry while preserving the same daughter momenta
        for record in [rows, hard_rows] if hard_rows else [rows]:
            daughters = [row.split() for row in record[-2:]]
            momentum = [sum(float(row[i]) for row in daughters) for i in range(6, 10)]
            parent = len(record) - 1
            for row in daughters:
                row[2:4] = [str(parent), "0"]
            record[-2:] = [particle_row(23, 2, *momentum), *(" ".join(row) for row in daughters)]
    blocks = []
    id1, id2, x1, x2 = 2, -2, 0.2, 0.2
    if hard_rows:
        id1, id2 = (int(row.split()[0]) for row in hard_rows[:2])
        x1, x2 = (float(row.split()[9]) / 100 for row in hard_rows[:2])
    for index in range(count):
        comment = ""
        if mode != "standard":
            comment = (
                f"# graniitti_diff version=1 accepted_events={index + 1} attempted_events={2 * (index + 1)} "
                f"side_mask={side_mask} id1={id1} id2={id2} xhard1={x1:.17g} xhard2={x2:.17g} "
                f"beam1_e=100 beam2_e=100 {side1} {side2}\n"
            )
        comment += "".join(f"# graniitti_hard {row}\n" for row in hard_rows)
        blocks.append(f"<event>\n{len(rows)} 1 1 10 0.0073 0.118\n" + "\n".join(rows) + "\n" + comment + "</event>\n")
    path.write_text(
        '<LesHouchesEvents version="3.0">\n<init>\n2212 2212 100 100 0 0 0 0 3 1\n1 0 1 1\n</init>\n'
        + "".join(blocks)
        + "</LesHouchesEvents>\n"
    )


# Run the installed converter without building it or changing the runtime environment
def run_converter(source, output, card, count=2, converter_mode=None, cwd=ROOT):
    mode_args = [] if converter_mode is None else [converter_mode]
    return subprocess.run(
        [str(DRIVER), str(source), str(output), str(count), "58031", str(card), "5", *mode_args],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=45,
        check=False,
    )


# Read stable particles and weights from the converter's HepMC3 output
def read_events(path):
    return [{
        'particles': [(p.pid(), [p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().e()])
                      for p in event.particles() if p.status() == 1],
        'beams': [(p.pid(), [p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().e()])
                  for p in event.beams()],
        'weights': list(event.weights()),
        'attrs': {name: event.attribute_as_string(name).split() for name in event.attribute_names()},
        'statuses': [p.status() for p in event.particles()],
    } for event in read_hepmc_events(path)]


# Require complete fragmentation with finite momenta, no final partons and exact pp closure
@pytest.mark.parametrize(
    "mode,flavour",
    [("direct", 21), ("IPp", 21), ("IPp", 2), ("IPp", 3), ("IPp", -5), ("IPIP", 21), ("full", 21), ("standard", 21)],
)
def test_converter_physical_modes(tmp_path, mode, flavour):
    source, output, card = tmp_path / "source.lhe", tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    write_lhe(source, mode, flavour=flavour)
    card.write_text("PartonLevel:ISR = off\nParticleDecays:limitTau0 = on\nParticleDecays:tau0Max = 10\n")
    result = run_converter(source, output, card)
    assert result.returncode == 0, result.stdout + result.stderr
    events = read_events(output)
    assert len(events) == 2
    for event in events:
        assert event["weights"] == pytest.approx([1.0])
        total = [sum(p[axis] for _, p in event["particles"]) for axis in range(4)]
        assert total == pytest.approx([0, 0, 0, 200], abs=0.011)
        for pdg, momentum in event["particles"]:
            assert all(math.isfinite(value) for value in momentum)
            assert momentum[3] >= 0
            assert abs(pdg) not in {1, 2, 3, 4, 5, 21, 2101, 2203}
        if flavour == -5:
            assert all(abs(pdg) != 5122 for pdg, _ in event["particles"])
        if mode in {"IPp", "IPIP"}:
            assert len(event["particles"]) > 6


# Preserve a readable graph when isolated remnant roots join the beam vertex
def test_converter_remnant_vertices(tmp_path):
    source = Path(__file__).with_name("data") / "isolated_remnants.lhe.txt"
    output, card = tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    card.write_text("PartonLevel:ISR = on\nParticleDecays:limitTau0 = on\nParticleDecays:tau0Max = 10\n")
    result = run_converter(source, output, card, count=1, converter_mode="fragment")
    assert result.returncode == 0, result.stdout + result.stderr
    reader = hepmc.ReaderAscii(str(output))
    try:
        event = hepmc.GenEvent()
        assert reader.read_event(event) and not reader.failed()
        assert len(event.beams()) == 2
        assert all(vertex.particles_in() or vertex.particles_out() for vertex in event.vertices())
    finally:
        reader.close()


# Normalize converted variable LHE weights with the process cross section
def test_converter_weight_normalization(tmp_path):
    source, output, card = tmp_path / "source.lhe", tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    write_lhe(source, "direct")
    source.write_text(source.read_text().replace(" 1 1 10 0.0073", " 1 2 10 0.0073", 1))
    card.write_text("PartonLevel:ISR = off\n")
    result = run_converter(source, output, card)
    assert result.returncode == 0, result.stdout + result.stderr
    summary = readers.read_hepmc3_weight_summary(str(output))
    assert summary["is_nonconstant"]
    assert readers.resolve_hepmc3_xsmode("auto", summary)[0] == "header"
    combined = readers.combine_hepmc3_weight_summaries([summary, summary])
    assert readers.resolve_hepmc3_xsmode("auto", combined)[0] == "header"
    data = readers.read_hepmc3(str(output), obs=[{}], pid=[[-13, 13]], cuts=[None])[0]
    assert data["xsection_pb"] == pytest.approx(1.0)
    assert data["weights"] == pytest.approx([2.0, 1.0])


# Collapse a subthreshold remnant with tagged-proton recoil and unchanged hard muons
@pytest.mark.parametrize("fixture,recoil", [("low_mass_remnants", True), ("beam_remnant_limit", False)])
def test_converter_low_mass_remnants(tmp_path, fixture, recoil):
    source = Path(__file__).with_name("data") / f"{fixture}.lhe.txt"
    output, card = tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    card.write_text("PartonLevel:ISR = on\n")
    result = run_converter(source, output, card, count=1)
    assert result.returncode == 0, result.stdout + result.stderr
    events = list(read_hepmc_events(output))
    assert len(events) == 1
    event = events[0]
    assert event.attribute_as_string("graniitti_pythia_remnant_fallback") == "1"
    rows = source.read_text().split("<event>\n", 1)[1].split("</event>", 1)[0].splitlines()
    nup = int(rows[0].split()[0])
    assert event.weights()[0] == pytest.approx(float(rows[0].split()[2]), rel=1e-12, abs=0)
    final = [p for p in event.particles() if p.status() == 1]
    assert all(abs(p.pid()) not in {1, 2, 3, 4, 5, 21} for p in final)
    for row in rows[1:nup + 1]:
        fields = row.split()
        if abs(int(fields[0])) != 13:
            continue
        muon = next(p for p in final if p.pid() == int(fields[0]))
        p4 = muon.momentum()
        assert [p4.px(), p4.py(), p4.pz(), p4.e()] == pytest.approx([float(x) for x in fields[6:10]], abs=1e-10)
    tagged = [p for p in final if "graniitti_diff_side" in p.attribute_names()]
    assert len(tagged) == 2
    if recoil:
        assert any(abs(float(p.attribute_as_string("graniitti_recoil_de"))) > 1e-6 for p in tagged)
    for component in ["px", "py", "pz", "e"]:
        delta = sum(getattr(p.momentum(), component)() for p in final) - sum(
            getattr(p.momentum(), component)() for p in event.beams())
        assert abs(delta) < 1e-7


# Reject input and output aliases before truncating the source file
@pytest.mark.parametrize("alias", ["same", "symlink", "hardlink"])
def test_converter_preserves_input_file(tmp_path, alias):
    source, output, card = tmp_path / "source.lhe", tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    write_lhe(source, "direct")
    card.write_text("PartonLevel:ISR = off\n")
    if alias == "same":
        output = source
    elif alias == "symlink":
        output.symlink_to(source)
    else:
        output.hardlink_to(source)
    original = source.read_bytes()
    result = run_converter(source, output, card)
    assert result.returncode != 0 and "same file" in result.stderr
    assert source.read_bytes() == original


# Reject invalid source events instead of writing a successful partial sample
@pytest.mark.parametrize("damage", ["truncated", "bad-row", "bad-counters", "nonfinite-metadata"])
def test_converter_rejects_invalid_input(tmp_path, damage):
    source, output, card = tmp_path / "source.lhe", tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    write_lhe(source, "direct")
    text = source.read_text()
    if damage == "truncated":
        text = text.rsplit("</event>", 1)[0]
    elif damage == "bad-row":
        text = text.replace("2212 -1", "not_a_particle -1", 1)
    elif damage == "nonfinite-metadata":
        text = text.replace("side1_lead_e=70", "side1_lead_e=70 side1_lead_py=nan")
    else:
        text = text.replace("accepted_events=2", "accepted_events=1")
    source.write_text(text)
    card.write_text("PartonLevel:ISR = off\n")
    result = run_converter(source, output, card)
    assert result.returncode > 0, result.stdout + result.stderr


# Do not silently fall back to default physics when a requested settings file is absent
def test_converter_requires_requested_settings(tmp_path):
    source = tmp_path / "source.lhe"
    write_lhe(source, "direct")
    result = run_converter(source, tmp_path / "out.hepmc3", tmp_path / "absent.cmnd")
    assert result.returncode > 0, result.stdout + result.stderr


# Exercise both explicit hard-diffraction handoffs using the same complete source record
@pytest.mark.parametrize("mode", ["IPp", "IPIP"])
@pytest.mark.parametrize("converter_mode", ["fragment", "shower", "auto"])
@pytest.mark.parametrize("isr", ["on", "off"])
def test_converter_explicit_modes(tmp_path, mode, converter_mode, isr):
    source, output, card = tmp_path / "source.lhe", tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    count = 32 if converter_mode == "auto" else 2
    write_lhe(source, mode, count=count)
    card.write_text(f"PartonLevel:ISR = {isr}\nParticleDecays:limitTau0 = on\nParticleDecays:tau0Max = 10\n")
    result = run_converter(source, output, card, count=count, converter_mode=converter_mode)
    assert result.returncode == 0, result.stdout + result.stderr
    events = read_events(output)
    assert len(events) == count
    for event in events:
        expected = "isolated-shower" if converter_mode == "fragment" else "full-shower"
        assert event["attrs"]["graniitti_converter_mode"] == [expected]
        total = [sum(p[axis] for _, p in event["particles"]) for axis in range(4)]
        assert total == pytest.approx([0, 0, 0, 200], abs=0.011)
        assert all(abs(pdg) not in {1, 2, 3, 4, 5, 21, 2101, 2203} for pdg, _ in event["particles"])


# Preserve source normalization and signed weights for both full samples and prefixes
@pytest.mark.parametrize("mode", ["standard", "direct", "IPIP"])
@pytest.mark.parametrize("count", [1, 2])
def test_converter_source_normalization(tmp_path, mode, count):
    source, output, card = tmp_path / "source.lhe", tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    write_lhe(source, mode)
    data = source.read_text().replace("0 0 0 0 3 1\n1 0 1 1", "0 0 0 0 -4 1\n0.75 0.125 1 1")
    data = data.replace(" 1 1 10 0.0073", " 1 2.5 10 0.0073", 1)
    data = data.replace(" 1 1 10 0.0073", " 1 -1.0 10 0.0073", 1)
    source.write_text(data)
    card.write_text("PartonLevel:ISR = off\n")
    result = run_converter(source, output, card, count=count)
    assert result.returncode == 0, result.stdout + result.stderr
    events = read_events(output)
    assert [event["weights"] for event in events] == [[2.5], [-1.0]][:count]
    summary = readers.read_hepmc3_weight_summary(str(output))
    assert readers.resolve_hepmc3_xsmode("auto", summary)[0] == "header"
    for event in events:
        assert [float(value) for value in event["attrs"]["GenCrossSection"][:2]] == pytest.approx([0.75, 0.125])
        assert float(event["attrs"]["graniitti_scalup"][0]) == pytest.approx(10.0)


# Preserve signed weights, cumulative source counters and exact muons when fragmenting
def test_converter_weights_and_matching(tmp_path):
    source, output, card = tmp_path / "source.lhe", tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    write_lhe(source, "IPIP")
    data = source.read_text().replace("0 0 0 0 3 1", "0 0 0 0 -4 1")
    data = data.replace("0 0 0 0 -4 1\n1 0 1 1", "0 0 0 0 -4 2\n1 0 1 1\n2 0 1 2")
    data = data.replace("8 1 1 10", "8 1 -2.5 10")
    data = data.replace("8 1 -2.5 10", "8 2 -2.5 10", 1)
    data = data.replace("accepted_events=1 attempted_events=2", "accepted_events=11 attempted_events=22")
    data = data.replace("accepted_events=2 attempted_events=4", "accepted_events=12 attempted_events=24")
    source.write_text(data)
    card.write_text("PartonLevel:ISR = off\n")
    result = run_converter(source, output, card, converter_mode="fragment")
    assert result.returncode == 0, result.stdout + result.stderr
    events = read_events(output)
    assert [e["weights"] for e in events] == [[-2.5], [-2.5]]
    assert all(float(e["attrs"]["GenCrossSection"][0]) == pytest.approx(3.0) for e in events)
    spec = importlib.util.spec_from_file_location(
        "pythia_matching", ROOT / "tests/external/pythia/drivers/lhe_converter/check_matching.py"
    )
    matching = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(matching)
    summary = matching.read_hepmc_summary(output, 2, matching.read_lhe_events(source, 2))
    assert summary["accepted_counters"] == [11, 12]
    assert summary["attempted_counters"] == [22, 24]
    assert summary["protected_mass_delta_max"] < 1e-8


# Fail immediately when the requested HepMC destination cannot be opened
def test_converter_output_failure(tmp_path):
    source, card = tmp_path / "source.lhe", tmp_path / "shower.cmnd"
    write_lhe(source, "IPIP")
    card.write_text("PartonLevel:ISR = off\n")
    result = run_converter(source, tmp_path / "absent" / "out.hepmc3", card)
    assert result.returncode > 0, result.stdout + result.stderr


# Keep asymmetric massive nuclear beams and their exact final four momenta
@pytest.mark.parametrize("central", ["mumu", "gg"])
def test_converter_complete_nuclear_record(tmp_path, central):
    source, output, card = tmp_path / "source.lhe", tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    nucleus, nuclear_mass, proton_mass = 1000180400, 37.2155, 0.93827208943
    initial_pz = math.sqrt(100**2 - nuclear_mass**2) - math.sqrt(100**2 - proton_mass**2)
    nuclear_pz = math.sqrt(70**2 - nuclear_mass**2)
    proton_pz = -math.sqrt(70**2 - proton_mass**2)
    muon_pz = (initial_pz - nuclear_pz - proton_pz) / 2
    muon_px = math.sqrt(30**2 - muon_pz**2)
    rows = [
        particle_row(nucleus, -1, 0, 0, math.sqrt(100**2 - nuclear_mass**2), 100),
        particle_row(2212, -1, 0, 0, -math.sqrt(100**2 - proton_mass**2), 100),
        particle_row(nucleus, 1, 0, 0, nuclear_pz, 70),
        particle_row(2212, 1, 0, 0, proton_pz, 70),
        particle_row(-13, 1, muon_px, 0, muon_pz, 30),
        particle_row(13, 1, -muon_px, 0, muon_pz, 30),
    ]
    if central == "gg":
        rows[-2] = particle_row(21, 1, muon_px, 0, muon_pz, 30, (501, 502))
        rows[-1] = particle_row(21, 1, -muon_px, 0, muon_pz, 30, (502, 501))
    source.write_text(
        f'<LesHouchesEvents version="3.0">\n<init>\n{nucleus} 2212 100 100 0 0 0 0 3 1\n1 0 1 1\n</init>\n'
        + "<event>\n6 1 1 10 0.0073 0.118\n"
        + "\n".join(rows)
        + f"\n# graniitti_diff version=1 accepted_events=1 attempted_events=2 side_mask=3 beam1_e=100 beam2_e=100 "
        f"side1_active=1 side1_lead_pdg={nucleus} side1_lead_pz={nuclear_pz:.17g} side1_lead_e=70 "
        f"side1_lead_m={nuclear_mass} side2_active=1 side2_lead_pdg=2212 side2_lead_pz={proton_pz:.17g} "
        f"side2_lead_e=70 side2_lead_m={proton_mass}\n</event>\n</LesHouchesEvents>\n"
    )
    card.write_text("PartonLevel:ISR = off\n")
    result = run_converter(source, output, card, count=1)
    assert result.returncode == 0, result.stdout + result.stderr
    event = read_events(output)[0]
    total = [sum(p[axis] for _, p in event["particles"]) for axis in range(4)]
    assert total == pytest.approx([0, 0, initial_pz, 200], abs=1e-8)
    assert event["attrs"]["graniitti_closure_pass"] == ["1"]


# Reject a repeated or missing color tag before spending shower attempts on it
def test_converter_rejects_invalid_color(tmp_path):
    source, output, card = tmp_path / "source.lhe", tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    write_lhe(source, "IPIP")
    source.write_text(source.read_text().replace("21 1 1 2 501 502", "21 1 1 2 501 501"))
    card.write_text("PartonLevel:ISR = off\n")
    result = run_converter(source, output, card)
    assert result.returncode > 0, result.stdout + result.stderr


# Let color strings shared by hard jets and remnants absorb the full-shower recoil
def test_converter_colored_hard_shower(tmp_path):
    source, output, card = tmp_path / "source.lhe", tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    write_lhe(source, "IPIP", jets=True)
    card.write_text("PartonLevel:ISR = on\nParticleDecays:limitTau0 = on\n")
    result = run_converter(source, output, card, converter_mode="shower")
    assert result.returncode == 0, result.stdout + result.stderr
    events = read_events(output)
    assert len(events) == 2
    for event in events:
        assert event["attrs"]["graniitti_recoil_protected_particles"] == ["0"]
        total = [sum(p[axis] for _, p in event["particles"]) for axis in range(4)]
        assert total == pytest.approx([0, 0, 0, 200], abs=0.011)
        assert all(abs(pdg) not in {1, 2, 3, 4, 5, 21, 2101, 2203} for pdg, _ in event["particles"])


# Exercise Pythia beam-remnant construction for reduced records with ISR enabled
def test_converter_reduced_shower_isr(tmp_path):
    source, output, card = tmp_path / "source.lhe", tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    write_lhe(source, "full")
    card.write_text("PartonLevel:ISR = on\nParticleDecays:limitTau0 = on\n")
    result = run_converter(source, output, card)
    assert result.returncode == 0, result.stdout + result.stderr
    events = read_events(output)
    assert len(events) == 2
    for event in events:
        assert event["attrs"]["graniitti_converter_mode"] == ["full-shower"]
        total = [sum(p[axis] for _, p in event["particles"]) for axis in range(4)]
        assert total == pytest.approx([0, 0, 0, 200], abs=0.011)


# Write exact nonzero transverse recoil with optional forward strings and Durham jets
def write_cep_lhe(path, central, excitations, metadata=True):
    proton_mass, diquark_mass = 0.93827208943, 0.65
    beam_pz = math.sqrt(100**2 - proton_mass**2)
    rows = [particle_row(2212, -1, 0, 0, beam_pz, 100), particle_row(2212, -1, 0, 0, -beam_pz, 100)]
    forward_pz, sides = [], []
    for side, sign in enumerate([1, -1], start=1):
        excited = side <= excitations
        mass = 4.0 if excited else proton_mass
        px, pz = sign * 0.7, sign * math.sqrt(70**2 - mass**2 - 0.7**2)
        forward_pz.append(pz)
        if excited:
            energy = (mass**2 - diquark_mass**2) / (2 * mass)
            tag = 600 + side
            rows += [
                particle_row(2, 1, px * energy / mass, energy, pz * energy / mass, 70 * energy / mass, (tag, 0)),
                particle_row(
                    2101,
                    1,
                    px * (mass - energy) / mass,
                    -energy,
                    pz * (mass - energy) / mass,
                    70 * (mass - energy) / mass,
                    (0, tag),
                ),
            ]
            sides.append(
                f"side{side}_active=1 side{side}_excited=1 side{side}_fragmented=1 side{side}_string=1 "
                f"side{side}_system_pdg=90210 side{side}_system_px={px} side{side}_system_pz={pz:.17g} "
                f"side{side}_system_e=70 side{side}_system_m={mass}"
            )
        else:
            rows.append(particle_row(2212, 1, px, 0, pz, 70))
            sides.append(
                f"side{side}_active=1 side{side}_lead_pdg=2212 side{side}_lead_px={px} "
                f"side{side}_lead_pz={pz:.17g} side{side}_lead_e=70 side{side}_lead_m={mass}"
            )
    half_pz = -sum(forward_pz) / 2
    muon_mass = 0.1056583755 if central == "mumu" else 0.0
    pt = math.sqrt(30**2 - half_pz**2 - muon_mass**2)
    ids, colors = {
        "mumu": ((-13, 13), ((0, 0), (0, 0))),
        "gg": ((21, 21), ((501, 502), (502, 501))),
        "qq": ((2, -2), ((501, 0), (0, 501))),
    }[central]
    rows += [
        particle_row(ids[0], 1, pt, 0, half_pz, 30, colors[0]),
        particle_row(ids[1], 1, -pt, 0, half_pz, 30, colors[1]),
    ]
    comment = ""
    if metadata:
        comment = (
            "# graniitti_diff version=1 accepted_events=1 attempted_events=2 side_mask=3 "
            f"elastic_cep={int(central != 'mumu')} beam1_e=100 beam2_e=100 " + " ".join(sides) + "\n"
        )
    path.write_text(
        '<LesHouchesEvents version="3.0">\n<init>\n2212 2212 100 100 0 0 0 0 3 1\n1 0 1 1\n</init>\n'
        + f"<event>\n{len(rows)} 1 1 20 0.0073 0.118\n"
        + "\n".join(rows)
        + "\n"
        + comment
        + "</event>\n</LesHouchesEvents>\n"
    )
    return [
        (int(row.split()[0]), [float(x) for x in row.split()[6:10]])
        for row in rows[2:]
        if row.split()[4:6] == ["0", "0"]
    ]


# Fragment kT-EPA excitations and Durham QCD systems while retaining exact colorless recoil
@pytest.mark.parametrize(
    "central,excitations,metadata",
    [
        ("mumu", 1, True),
        ("mumu", 2, True),
        ("gg", 0, True),
        ("qq", 0, True),
        ("gg", 1, True),
        ("gg", 0, False),
        ("mumu", 1, False),
        ("mumu", 0, False),
    ],
)
def test_converter_epa_and_durham(tmp_path, central, excitations, metadata):
    source, output, card = tmp_path / "source.lhe", tmp_path / "out.hepmc3", tmp_path / "shower.cmnd"
    exact_particles = write_cep_lhe(source, central, excitations, metadata)
    card.write_text("PartonLevel:ISR = on\nParticleDecays:limitTau0 = on\nParticleDecays:tau0Max = 10\n")
    result = run_converter(source, output, card, count=1)
    assert result.returncode == 0, result.stdout + result.stderr
    events = read_events(output)
    assert len(events) == 1
    event = events[0]
    assert [pdg for pdg, _ in event["beams"]] == [2212, 2212]
    total = [sum(p[axis] for _, p in event["particles"]) for axis in range(4)]
    assert total == pytest.approx([0, 0, 0, 200], abs=1e-8)
    assert all(abs(pdg) not in {1, 2, 3, 4, 5, 21, 2101, 2203} for pdg, _ in event["particles"])
    for pdg, momentum in exact_particles:
        assert any(out_id == pdg and p == pytest.approx(momentum, abs=1e-8) for out_id, p in event["particles"])


# Read only the requested prefix and reproduce the same hadrons with the same seed
def test_converter_prefix_and_seed(tmp_path):
    source, card = tmp_path / "source.lhe", tmp_path / "shower.cmnd"
    write_lhe(source, "IPIP")
    first, second = source.read_text().split("</event>", maxsplit=1)
    source.write_text(first + "</event>" + second.replace("2212 -1", "invalid -1", 1))
    card.write_text("PartonLevel:ISR = off\nPartonLevel:FSR = off\n")
    events = []
    for index in range(2):
        output = tmp_path / f"out{index}.hepmc3"
        result = run_converter(source, output, card, count=1)
        assert result.returncode == 0, result.stdout + result.stderr
        converted = read_events(output)
        assert len(converted) == 1
        events.append(converted[0])
        # With FSR disabled there must be no Pythia shower branch in the event history
        assert 51 not in converted[0]["statuses"]
    assert len(events[0]["particles"]) == len(events[1]["particles"])
    for (id1, p1), (id2, p2) in zip(events[0]["particles"], events[1]["particles"], strict=True):
        assert id1 == id2
        assert p1 == pytest.approx(p2, abs=1e-12)


# Reject invalid steering even when the source could be copied without a shower
@pytest.mark.parametrize("mode", ["direct", "IPIP", "standard"])
def test_invalid_pythia_card(tmp_path, mode):
    source, output, card = (tmp_path / name for name in ("source.lhe", "out.hepmc3", "bad.cmnd"))
    write_lhe(source, mode)
    card.write_text("PartonLevel:ISR = off\nUnknownPhysics:setting = on\n")
    result = run_converter(source, output, card)
    assert result.returncode != 0
    assert not output.exists()


# Preserve steering files and use the compiled data path outside the repository root
def test_converter_paths(tmp_path):
    source, output, card = (tmp_path / name for name in ("source.lhe", "out.hepmc3", "shower.cmnd"))
    write_lhe(source, "IPIP")
    card.write_text("PartonLevel:ISR = off\n")
    before = card.read_bytes()
    result = run_converter(source, card, card, cwd=tmp_path)
    assert result.returncode != 0
    assert card.read_bytes() == before
    result = run_converter(source, output, card, cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(read_hepmc_events(output)) == 2


# Check low mass dimuon matching without requiring a Z resonance
@pytest.mark.parametrize("mode", ["direct", "IPIP"])
def test_matching_nonresonant_muons(tmp_path, mode):
    from tests.external.pythia.drivers.lhe_converter.check_matching import validate_pair
    source, output, card = (tmp_path / name for name in ("source.lhe", "out.hepmc3", "shower.cmnd"))
    write_lhe(source, mode)
    card.write_text("PartonLevel:ISR = off\n")
    result = run_converter(source, output, card)
    assert result.returncode == 0, result.stdout + result.stderr
    validate_pair(source, output, 2)
    with pytest.raises(RuntimeError, match="no Z peak"):
        validate_pair(source, output, 2, require_z_peak=True)


# An explicit all count consumes the complete LHE file instead of a short default prefix
def test_convert_all(tmp_path):
    source, output, card = (tmp_path / name for name in ("source.lhe", "out.hepmc3", "shower.cmnd"))
    write_lhe(source, "direct", count=21)
    card.write_text("PartonLevel:ISR = off\n")
    result = run_converter(source, output, card, count="all")
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(read_hepmc_events(output)) == 21


# Reject merging schemes whose signed subtraction weights this converter does not implement
@pytest.mark.parametrize("setting", ["doPTLundMerging", "doUMEPSTree", "doUNLOPSTree"])
def test_unsupported_merging(tmp_path, setting):
    source, output, card = (tmp_path / name for name in ("source.lhe", "out.hepmc3", "merging.cmnd"))
    write_lhe(source, "full")
    card.write_text(f"Merging:{setting} = on\n")
    result = run_converter(source, output, card, converter_mode="shower")
    assert result.returncode != 0
    assert not output.exists()


# A mixed full-shower source must never silently drop the events without diffraction metadata
def test_mixed_diffraction_metadata(tmp_path):
    source, output, card = (tmp_path / name for name in ("source.lhe", "out.hepmc3", "shower.cmnd"))
    write_lhe(source, "IPIP")
    lines = source.read_text().splitlines(keepends=True)
    first = next(i for i, line in enumerate(lines) if line.startswith("# graniitti_diff"))
    source.write_text("".join(line for i, line in enumerate(lines) if i != first))
    card.write_text("PartonLevel:ISR = off\n")
    result = run_converter(source, output, card, converter_mode="shower")
    assert result.returncode != 0
    assert not output.exists()


# Match the original hard muons when the requested QED shower adds photons to their final state
def test_matching_radiated_muons(tmp_path):
    from tests.external.pythia.drivers.lhe_converter.check_matching import validate_pair
    source, output, card = (tmp_path / name for name in ("source.lhe", "out.hepmc3", "shower.cmnd"))
    write_lhe(source, "direct", count=32)
    card.write_text("PartonLevel:ISR = off\nTimeShower:QEDshowerByL = on\n")
    result = run_converter(source, output, card, count=32, converter_mode="fragment")
    assert result.returncode == 0, result.stdout + result.stderr
    events = read_hepmc_events(output)
    assert len(events) == 32
    assert any(p.pid() == 22 and p.status() == 1 for event in events for p in event.particles())
    validate_pair(source, output, 32)

#!/usr/bin/env python3
# Validate GRANIITTI LHE to Pythia HepMC3 matching
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
import math
from itertools import islice
from pathlib import Path

from pyHepMC3 import HepMC3 as hepmc


# Compute invariant mass of one four-vector tuple
def mass(p4):
    px, py, pz, e = p4
    return math.sqrt(max(0.0, e * e - px * px - py * py - pz * pz))


# Compute the sum of two four-vector tuples
def add_p4(first, second):
    return tuple(first[i] + second[i] for i in range(4))


# Parse key-value fields from one GRANIITTI LHE metadata comment
def parse_lhe_metadata(text):
    values = {}
    for token in text.split():
        if "=" not in token:
            continue
        key, value = token.split("=", maxsplit=1)
        values[key] = value
    return values


# Compute the first max_events LHE weights and cumulative source counters
def read_lhe_events(path, max_events):
    events = []
    event = None
    with Path(path).open() as input_file:
        for line in input_file:
            text = line.strip()
            if text == "<event>" or text.startswith("<event "):
                event = {"weight": None, "accepted": None, "attempted": None, "muons": []}
                continue
            if text == "</event>":
                if event is None or event["weight"] is None:
                    raise RuntimeError(f"{path}: incomplete LHE event")
                events.append(event)
                event = None
                if len(events) >= max_events:
                    break
                continue
            if event is None or not text:
                continue
            if text.startswith("# graniitti_diff"):
                metadata = parse_lhe_metadata(text)
                event["accepted"] = int(metadata["accepted_events"])
                event["attempted"] = int(metadata["attempted_events"])
            elif text[0] not in "<#" and event["weight"] is None:
                event["weight"] = float(text.split()[2])
            elif text[0] not in "<#":
                fields = text.split()
                if len(fields) >= 13 and fields[1] == "1" and abs(int(fields[0])) == 13:
                    event["muons"].append((int(fields[0]), tuple(float(x) for x in fields[6:10])))
    if event is not None:
        raise RuntimeError(f"{path}: truncated LHE event")
    return events


# Read one event-level attribute as a float
def event_attr(event, name):
    if name not in event["attrs"]:
        raise RuntimeError(f"missing HepMC3 event attribute {name}")
    value = float(event["attrs"][name])
    if not math.isfinite(value):
        raise RuntimeError(f"nonfinite HepMC3 event attribute {name}")
    return value


# Read one event-level attribute as text
def event_text_attr(event, name):
    if name not in event["attrs"]:
        raise RuntimeError(f"missing HepMC3 event attribute {name}")
    return event["attrs"][name]


# Compute true when two floating-point metadata values agree
def metadata_values_match(expected, observed):
    scale = max(1.0, abs(expected), abs(observed))
    return abs(expected - observed) <= 1.0e-12 * scale


# Validate one side of the event-local Pythia beam configuration
def validate_pythia_beam_side(event, side):
    xi = event_attr(event, f"graniitti_diff_side{side}_xi")
    beta = event_attr(event, f"graniitti_diff_side{side}_beta")
    hard_x = event_attr(event, f"graniitti_diff_xhard{side}")
    original_energy = event_attr(event, f"graniitti_diff_beam{side}_e")
    pythia_id = round(event_attr(event, f"graniitti_pythia_beam{side}_id"))
    pythia_energy = event_attr(event, f"graniitti_pythia_beam{side}_e")
    pythia_x = event_attr(event, f"graniitti_pythia_x{side}")

    active = xi > 0.0
    physical_energy = xi * original_energy if active else original_energy
    expected_energy = physical_energy
    expected_x = beta if active else hard_x
    hard_energy_key = f"graniitti_pythia_hard{side}_e"
    if hard_energy_key in event["attrs"]:
        expected_x = event_attr(event, hard_energy_key) / expected_energy
    if not 0.0 < pythia_x < 1.0:
        raise RuntimeError(f"side {side}: Pythia hard fraction is outside (0, 1)")
    if active and pythia_id != 990:
        raise RuntimeError(f"side {side}: active DPDF beam id is {pythia_id}, not 990")
    if not active and pythia_id == 990:
        raise RuntimeError(f"side {side}: ordinary PDF side incorrectly uses beam id 990")
    if not metadata_values_match(expected_energy, pythia_energy):
        raise RuntimeError(f"side {side}: Pythia beam energy {pythia_energy} != {expected_energy}")
    if not metadata_values_match(expected_x, pythia_x):
        raise RuntimeError(f"side {side}: Pythia x {pythia_x} != {expected_x}")


# Validate event-local Pomeron beams for the full diffractive shower mode
def validate_pythia_beam_configuration(event):
    if event_text_attr(event, "graniitti_converter_mode") != "full-shower":
        return
    validate_pythia_beam_side(event, 1)
    validate_pythia_beam_side(event, 2)


# Read the complete event graph through the official HepMC3 library
def iter_hepmc_events(path):
    if not Path(path).is_file():
        raise RuntimeError(f"Missing HepMC3 input: {path}")
    reader = hepmc.ReaderAscii(str(path))
    try:
        while True:
            event = hepmc.GenEvent()
            if not reader.read_event(event):
                raise RuntimeError(f"{path}: invalid HepMC3 event graph")
            if reader.failed():
                if event.particles():
                    raise RuntimeError(f"{path}: truncated HepMC3 event")
                break
            event.set_units(hepmc.Units.GEV, event.length_unit())
            xs = event.cross_section()
            beams = list(event.beams())
            final = [p for p in event.particles() if p.end_vertex() is None]
            if not beams or not final:
                raise RuntimeError(f"{path}: event without incoming beams or final particles")
            residual = [math.fsum([p.momentum()[i] for p in beams] + [-p.momentum()[i] for p in final])
                        for i in range(4)]
            if not all(math.isfinite(value) for value in residual):
                raise RuntimeError(f"{path}: nonfinite momentum closure")
            closure = max(abs(value) for value in residual)
            yield {
                "closure": closure,
                "weight": float(event.weights()[0]) if event.weights() else None,
                "cross_section": None if xs is None else {
                    "value": xs.xsec(), "error": xs.xsec_err(),
                    "accepted": xs.get_accepted_events(), "attempted": xs.get_attempted_events(),
                },
                "attrs": {name: event.attribute_as_string(name) for name in event.attribute_names()},
                "particles": [{"pid": p.pid(), "status": p.status(),
                               "p4": (p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().e())}
                              for p in event.particles()],
            }
    finally:
        reader.close()


# Compute the opposite-sign stable dimuon mass closest to the Z pole
def selected_dimuon_mass(event):
    muons = [
        (particle["pid"], particle["p4"])
        for particle in event["particles"]
        if particle["status"] == 1 and abs(particle["pid"]) == 13
    ]

    best_mass = None
    best_score = None
    for i, first in enumerate(muons):
        for second in muons[i + 1 :]:
            if first[0] * second[0] > 0:
                continue
            m = mass(add_p4(first[1], second[1]))
            if m < 66.0 or m > 116.0:
                continue
            score = abs(m - 91.1876)
            if best_score is None or score < best_score:
                best_mass = m
                best_score = score
    return best_mass


# Compare source hard muons, including their intermediate copies before photon radiation
def copied_muon_mass_delta(event, source):
    isolated = event_text_attr(event, "graniitti_converter_mode") == "isolated-shower"
    remaining = [p for p in event["particles"] if abs(p["pid"]) == 13 and (p["status"] == 1 or isolated)]
    before = (0.0, 0.0, 0.0, 0.0)
    after = before
    for pdg, p4 in source["muons"]:
        match = next(
            (
                p
                for p in remaining
                if p["pid"] == pdg
                and all(abs(x - y) <= 1e-8 * max(1.0, abs(x)) for x, y in zip(p4, p["p4"], strict=True))
            ),
            None,
        )
        if match is None:
            raise RuntimeError("A copied hard muon changed momentum or disappeared")
        before = add_p4(before, p4)
        after = add_p4(after, match["p4"])
        remaining.remove(match)
    return abs(mass(after) - mass(before))


# Compute a summary of HepMC3 event weights and matching diagnostics
def read_hepmc_summary(path, max_events, lhe_events):
    weights = []
    closures = []
    protected_mass_deltas = []
    accepted_counters = []
    attempted_counters = []
    mass_bins = [0.0, 0.0, 0.0, 0.0, 0.0]
    bin_edges = [66.0, 76.0, 86.0, 96.0, 106.0, 116.0]

    for event in islice(iter_hepmc_events(path), max_events):
        if event["weight"] is None or not math.isfinite(event["weight"]):
            raise RuntimeError(f"{path}: event without a finite HepMC3 weight")
        weight = event["weight"]
        factor = float(event["attrs"].get("graniitti_merging_weight", "1"))
        if not math.isfinite(factor) or factor < 0:
            raise RuntimeError(f"{path}: invalid merging weight")
        if len(weights) >= len(lhe_events) or not weights_match(weight, lhe_events[len(weights)]["weight"] * factor):
            raise RuntimeError("HepMC3 weight does not match LHE XWGTUP times the merging weight")
        weights.append(weight)
        validate_pythia_beam_configuration(event)
        if event["cross_section"] is None:
            raise RuntimeError(f"{path}: event without GenCrossSection counters")
        accepted_counters.append(event["cross_section"]["accepted"])
        attempted_counters.append(event["cross_section"]["attempted"])

        recorded = event_attr(event, "graniitti_closure_max_abs")
        tolerance = event_attr(event, "graniitti_closure_tolerance")
        if tolerance <= 0.0 or recorded < 0.0 or recorded > tolerance:
            raise RuntimeError("Invalid Pythia momentum closure or tolerance")
        if event_text_attr(event, "graniitti_closure_pass") != "1":
            raise RuntimeError("Pythia momentum closure failed")
        if event_text_attr(event, "graniitti_converter_mode") not in {"isolated-shower", "full-shower"}:
            tolerance = min(tolerance, 1.0e-6)
        closure = event["closure"]
        if not math.isfinite(closure) or closure > tolerance:
            raise RuntimeError(f"closure too large: {closure:.6e}, tolerance: {tolerance:.6e}")
        closures.append(closure)
        if event_text_attr(event, "graniitti_converter_mode") == "full-shower":
            delta = abs(
                event_attr(event, "graniitti_protected_mass_after_recoil")
                - event_attr(event, "graniitti_protected_mass_before_recoil")
            )
        else:
            if len(weights) > len(lhe_events):
                raise RuntimeError("HepMC contains more events than the LHE source")
            delta = copied_muon_mass_delta(event, lhe_events[len(weights) - 1])
        protected_mass_deltas.append(delta)

        m = selected_dimuon_mass(event)
        if m is None:
            continue
        for i in range(len(mass_bins)):
            if bin_edges[i] <= m < bin_edges[i + 1]:
                mass_bins[i] += weight
                break

    return {
        "weights": weights,
        "accepted_counters": accepted_counters,
        "attempted_counters": attempted_counters,
        "closure_max": max(closures) if closures else float("nan"),
        "protected_mass_delta_max": (max(protected_mass_deltas) if protected_mass_deltas else float("nan")),
        "mass_bins": mass_bins,
    }


# Compute true when two floating-point weights agree at LHE precision
def weights_match(expected, observed):
    scale = max(1.0, abs(expected), abs(observed))
    return abs(expected - observed) <= 1.0e-12 * scale


# Validate one LHE and HepMC3 file pair
def validate_pair(lhe_path, hepmc_path, max_events, *, require_z_peak=False):
    if max_events <= 0:
        raise ValueError("max_events must be positive")
    lhe_events = read_lhe_events(lhe_path, max_events)
    lhe_weights = [event["weight"] for event in lhe_events]
    summary = read_hepmc_summary(hepmc_path, max_events, lhe_events)
    hepmc_weights = summary["weights"]

    if not lhe_weights or len(lhe_weights) != len(hepmc_weights):
        raise RuntimeError(f"event count mismatch: LHE={len(lhe_weights)} HepMC3={len(hepmc_weights)}")

    expected_accepted = [event["accepted"] for event in lhe_events]
    expected_attempted = [event["attempted"] for event in lhe_events]
    if summary["accepted_counters"] != expected_accepted:
        raise RuntimeError("HepMC3 accepted-event counters do not match LHE metadata")
    if summary["attempted_counters"] != expected_attempted:
        raise RuntimeError("HepMC3 attempted-event counters do not match LHE metadata")
    if any(
        attempted < accepted
        for accepted, attempted in zip(summary["accepted_counters"], summary["attempted_counters"], strict=True)
    ):
        raise RuntimeError("HepMC3 attempted-event counter is below accepted events")
    if any(
        current < previous
        for previous, current in zip(summary["attempted_counters"], summary["attempted_counters"][1:], strict=False)
    ):
        raise RuntimeError("HepMC3 attempted-event counters are not cumulative")
    if summary["protected_mass_delta_max"] > 1.0e-9:
        raise RuntimeError(
            f"protected hard-branch mass changed under recoil: {summary['protected_mass_delta_max']:.6e}"
        )

    bins = summary["mass_bins"]
    if require_z_peak and bins[2] <= max(bins[0], bins[1], bins[3], bins[4]):
        raise RuntimeError(f"weighted dimuon mass histogram has no Z peak: {bins}")

    print(f"{hepmc_path}")
    print(f"  events: {len(hepmc_weights)}")
    print(f"  weight range: {min(hepmc_weights):.6e} .. {max(hepmc_weights):.6e}")
    print(f"  weight mode: {'weighted' if not weights_match(max(hepmc_weights), min(hepmc_weights)) else 'unweighted'}")
    print(f"  closure max: {summary['closure_max']:.6e}")
    print(f"  protected mass delta max: {summary['protected_mass_delta_max']:.6e}")
    print(
        f"  source counters: accepted={summary['accepted_counters'][-1]} attempted={summary['attempted_counters'][-1]}"
    )
    print(f"  weighted mass bins [66,76,86,96,106,116]: {bins}")


# Run the command-line validation
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lhe", required=True, help="Input GRANIITTI LHE file")
    parser.add_argument("--hepmc", required=True, help="Converted HepMC3 file")
    parser.add_argument("--max-events", type=int, default=10000)
    parser.add_argument("--require-z-peak", action="store_true", help="Also require a populated Z mass peak")
    args = parser.parse_args()
    if args.max_events <= 0:
        parser.error("--max-events must be positive")
    validate_pair(args.lhe, args.hepmc, args.max_events, require_z_peak=args.require_z_peak)


if __name__ == "__main__":
    main()

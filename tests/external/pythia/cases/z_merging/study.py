# Run the SD CKKW-L icepack and compare absolute HepMC3 distributions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

from pyHepMC3 import HepMC3 as h
from pyHepMC3 import std

ROOT = Path(__file__).resolve().parents[5]
PACK = ROOT / "icepack/HARDPOM/z_merging"
WORK = ROOT / "tmp/z_merging"


# Hash steering, executables and outputs before deciding whether a sample can be reused
def digest(path):
    with Path(path).open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


# Run a generator or analysis command and reuse only complete, unchanged samples
def run(command, name, *, output=None, inputs=(), reuse=False, work=None):
    work = WORK if work is None else work
    command = [str(arg) for arg in command]
    record = work / f"{name}.command.json"
    state = {"command": command, "inputs": {str(p): digest(p) for p in inputs}}
    if reuse and output is not None and Path(output).is_file() and record.is_file():
        try:
            previous = json.loads(record.read_text())
            if previous == {**state, "output": digest(output)}:
                print(f"Reuse {name}", flush=True)
                return
        except (OSError, ValueError):
            pass
    print(name, flush=True)
    with (work / f"{name}.log").open("w") as log:
        subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    if output is not None:
        state["output"] = digest(output)
    record.write_text(json.dumps(state, indent=2))


# Iterate complete events with the official HepMC3 reader
def events(path):
    if not Path(path).is_file():
        raise ValueError(f"Missing HepMC3 input: {path}")
    reader = h.ReaderAscii(str(path))
    try:
        while not reader.failed():
            event = h.GenEvent()
            if not reader.read_event(event):
                raise ValueError(f"{path}: invalid HepMC3 event")
            if reader.failed():
                if event.particles():
                    raise ValueError(f"{path}: truncated HepMC3 event")
                break
            if any(p.pid() == 0 and p.status() == 0 for p in event.particles()):
                raise ValueError(f"{path}: incomplete HepMC3 particle")
            yield event
    finally:
        reader.close()


# Add independent multiplicities using their absolute cross sections and full weight sums
def combine(paths, output):
    paths, output = [Path(p) for p in paths], Path(output)
    if not paths:
        raise ValueError("At least one multiplicity is required")
    summaries = []
    for index, path in enumerate(paths):
        if (output.exists() and path.samefile(output)) or any(path.samefile(p) for p in paths[:index]):
            raise ValueError("Each input and the output must be a different file")
        total, count, nonzero, xs = 0.0, 0, 0, None
        for event in events(path):
            weights = list(event.weights())
            if len(weights) != 1 or not math.isfinite(weights[0]):
                raise ValueError(f"{path}: expected one finite nominal weight")
            total += weights[0]
            count += 1
            nonzero += abs(weights[0]) > 0
            xs = event.cross_section()
        if count == 0 or xs is None or not all(math.isfinite(v) and v >= 0 for v in (xs.xsec(), xs.xsec_err())):
            raise ValueError(f"{path}: empty sample or invalid terminal cross section")
        if not math.isfinite(total) or (total <= 0 and (nonzero > 0 or xs.xsec() > 0)):
            raise ValueError(f"{path}: nonpositive weight sum with a nonzero contribution")
        summaries.append((xs.xsec(), xs.xsec_err(), xs.xsec() / total if nonzero else 0.0, count))
    sigma = sum(s[0] for s in summaries)
    error = math.sqrt(sum(s[1] ** 2 for s in summaries))
    run_info = h.GenRunInfo()
    run_info.set_weight_names(std.vector_std_string(["Weight"]))
    writer = h.WriterAscii(str(output), run_info)
    writer.set_precision(17)
    count = 0
    try:
        if writer.failed():
            raise RuntimeError(f"Cannot write {output}")
        for path, (_, _, scale, expected) in zip(paths, summaries, strict=True):
            read = 0
            for event in events(path):
                # The Python weights accessor returns a copy, so update the serialized event data
                data = h.GenEventData()
                event.write_data(data)
                data.weights = std.vector_double([event.weights()[0] * scale])
                event.read_data(data)
                event.set_run_info(run_info)
                event.add_attribute("graniitti_sample_scale", h.DoubleAttribute(scale))
                event.set_event_number(count)
                count += 1
                read += 1
                header = h.GenCrossSection()
                header.set_cross_section(sigma, error, count, count)
                event.set_cross_section(header)
                writer.write_event(event)
                if writer.failed():
                    raise RuntimeError(f"Cannot write {output}")
            if read != expected:
                raise ValueError(f"{path}: event count changed during combination")
    finally:
        writer.close()
    return sigma


# Check the converter weight identity and the final normalization independently
def check(path):
    raw, merged, squares, zeros, count = 0.0, 0.0, 0.0, 0, 0
    remnant_weight = 0.0
    source_xs = None
    for event in events(path):
        if len(event.weights()) != 1:
            raise ValueError(f"{path}: expected one nominal merging weight")
        weight = event.weights()[0]
        source = float(event.attribute_as_string("graniitti_lhe_weight"))
        current_xs = float(event.attribute_as_string("graniitti_merging_source_xsec"))
        if source_xs is not None and not math.isclose(current_xs, source_xs, rel_tol=1e-12):
            raise ValueError(f"{path}: inconsistent source cross sections")
        source_xs = current_xs
        factor = float(event.attribute_as_string("graniitti_merging_weight"))
        if not all(math.isfinite(v) for v in (weight, source, source_xs, factor)):
            raise ValueError(f"{path}: nonfinite merging metadata")
        if factor < 0 or not math.isclose(weight, source * factor, rel_tol=1e-12, abs_tol=1e-15):
            raise ValueError(f"{path}: inconsistent merging weight")
        if event.attribute_as_string("graniitti_closure_pass") != "1":
            raise ValueError(f"{path}: four-momentum conservation failed")
        raw += source
        merged += weight
        squares += weight * weight
        zeros += factor <= 0
        if event.attribute_as_string("graniitti_merging_remnant_veto") == "1":
            remnant_weight += source
        count += 1
        xs = event.cross_section()
        if xs is None:
            raise ValueError(f"{path}: missing cross section")
        header = xs.xsec()
        source_sum = float(event.attribute_as_string("graniitti_merging_source_sum"))
    if count == 0 or not raw > 0:
        raise ValueError(f"{path}: empty sample or nonpositive source weight sum")
    if not math.isclose(raw, source_sum, rel_tol=1e-12):
        raise ValueError(f"{path}: incomplete source sample")
    expected = source_xs * merged / raw
    if not math.isclose(header, expected, rel_tol=1e-7):
        raise ValueError(f"{path}: incorrect merged cross section")
    return dict(
        events=count,
        zero_weights=zeros,
        source_pb=source_xs,
        merged_pb=header,
        effective_events=merged**2 / squares if squares > 0 else 0.0,
        remnant_veto_fraction=remnant_weight / raw,
    )


# Generate two multiplicities, vary the merging scale, and compare with native SD
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", default="z_merging", help="Output and log directory name")
    parser.add_argument("--reuse", action="store_true", help="Reuse existing generator and converter outputs")
    parser.add_argument("--events", type=int, default=None, help="Events per sample, including native Pythia")
    parser.add_argument("--native-events", type=int, default=None, help="Override the native Pythia event count")
    args = parser.parse_args()
    os.chdir(ROOT)
    tag = args.tag
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", tag) is None:
        parser.error("--tag must be a filename stem without directories")
    work = ROOT / "tmp" / tag
    work.mkdir(parents=True, exist_ok=True)
    config = json.loads((PACK / "study.json").read_text())
    dataset = json.loads((PACK / "dataset.json").read_text())
    count = dataset["validation"]["nevents"] if args.events is None else args.events
    native_count = args.native_events if args.native_events is not None else (
        config["native_events"] if args.events is None else count
    )
    if min(count, native_count) <= 0:
        parser.error("Event counts must be positive")
    seed = str(config["seed"])
    sources = [Path(f"output/{tag}_{n}.lhe") for n in range(2)]
    model_cards = sorted((ROOT / "modeldata").rglob("*.json"))
    for n, source in enumerate(sources):
        card = PACK / f"gencard_{n}.json"
        run(["bin/gr", "-i", card, "-o", source.stem, "-n", count], f"generate_{n}",
            output=source, inputs=[ROOT / "bin/gr", card, *model_cards], reuse=args.reuse, work=work)
    outputs, labels, diagnostics = [], [], {}
    for scale in config["merging_scales"]:
        parts = []
        for n, source in enumerate(sources):
            name = f"{tag}_t{scale}_{n}"
            card, output = work / f"{name}.cmnd", Path("output") / f"{name}.hepmc3"
            card.write_text(
                (PACK / "merging.cmnd").read_text() + f"\nMerging:TMS = {scale}\nMerging:nRequested = {n}\n"
            )
            run(
                [
                    "bin/pythia_lhe_hadronize",
                    str(source),
                    str(output),
                    str(count),
                    str(config["seed"] + n),
                    str(card),
                    "100",
                    "shower",
                ],
                name, output=output, inputs=[ROOT / "bin/pythia_lhe_hadronize", source, card], reuse=args.reuse, work=work,
            )
            diagnostics[name] = check(output)
            if diagnostics[name]["events"] != count:
                raise ValueError(f"{output}: incomplete or different sample, rerun without --reuse")
            parts.append(output)
        combined = Path(f"output/{tag}_t{scale}.hepmc3")
        combine(parts, combined)
        outputs.append(str(combined))
        labels.append(f"CKKW-L Qcut={scale} GeV")
    for name, command, label in [
        (
            f"{tag}_shower",
            [
                "bin/pythia_lhe_hadronize",
                str(sources[0]),
                f"output/{tag}_shower.hepmc3",
                str(count),
                seed,
                str(PACK / "shower.cmnd"),
                "100",
                "shower",
            ],
            "GRANIITTI Z + shower",
        ),
        (
            f"{tag}_native",
            [
                "bin/pythia_zmumu_hepmc3",
                str(PACK / "pythia.cmnd"),
                f"output/{tag}_native.hepmc3",
                str(native_count),
                seed,
                "--hard-muon-cuts",
                "off",
            ],
            "Native Pythia SD",
        ),
    ]:
        output = f"output/{name}.hepmc3"
        inputs = [ROOT / command[0], PACK / ("pythia.cmnd" if name.endswith("native") else "shower.cmnd")]
        if name.endswith("shower"):
            inputs.append(sources[0])
        run(command, name, output=output, inputs=inputs, reuse=args.reuse, work=work)
        expected = native_count if name.endswith("native") else count
        if sum(1 for _ in events(output)) != expected:
            raise ValueError(f"{output}: incomplete event sample")
        outputs.append(output)
        labels.append(label)
    (work / "weights.json").write_text(json.dumps(diagnostics, indent=2))
    run(
        [
            sys.executable,
            "-m",
            "core.iceplot",
            "--hepmc3",
            *outputs,
            "--mclabel",
            *labels,
            "--cuts",
            str(PACK / "cuts.py"),
            "--obs",
            str(PACK / "obs.py"),
            "--pid",
            "[[13,-13]]",
            "--unit",
            "pb",
            "--output",
            tag,
            "--report",
            str(work / "iceplot.json"),
            "--cores",
            "2",
        ],
        "iceplot", work=work,
    )
    print(f"Plots: figs/iceplot/{tag}, report: {work / 'iceplot.json'}", flush=True)


if __name__ == "__main__":
    main()

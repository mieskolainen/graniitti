# Validate prepared amplitude predictions against native evaluations on common events
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
import time
from pathlib import Path

import numpy as np
import torch
from pyHepMC3 import HepMC3 as hepmc3

from core.io.files import ensure_dir
from core.io.serialize import load_json_file, write_json_file
from core.numerics import array
from core.tune.drivers.graniitti.ampfit.amplitude import AmplitudeBank, evaluate_basis, input_card
from core.tune.drivers.graniitti.driver import GraniittiDriver
from core.tune.optimizers.ampfit.config import load_settings
from core.tune.tunesetup import load_tunesetup


# Copy selected events through the official HepMC3 reader and writer
def copy_events(source, destination, count, indices=None):
    indices = np.arange(count) if indices is None else np.asarray(indices)
    if (count < 1 or indices.shape != (count,) or not np.issubdtype(indices.dtype, np.integer)
            or indices[0] < 0 or np.any(np.diff(indices) <= 0)):
        raise ValueError("Amplitude validation requires increasing nonnegative event indices")
    reader = hepmc3.ReaderAscii(str(source))
    writer = hepmc3.WriterAscii(str(destination))
    copied = 0
    try:
        for index in range(int(indices[-1]) + 1):
            event = hepmc3.GenEvent()
            reader.read_event(event)
            if reader.failed():
                raise ValueError("Amplitude validation sample contains too few events")
            if index != indices[copied]:
                continue
            writer.write_event(event)
            if writer.failed():
                raise OSError("Cannot write amplitude validation events")
            copied += 1
    finally:
        reader.close()
        writer.close()
    return copied


# Evaluate explicit fit coordinates with the native program and retain all comparison outputs
def validate_points(*, driver, bank, run, points, directory, cdir, controls, events, timeout, data=None):
    directory = Path(directory).resolve()
    ensure_dir(directory, exist_ok=False)
    count = min(events, bank.amplitudes.shape[1])
    if count < 1 or not points:
        raise ValueError("Amplitude validation requires events and parameter points")
    source = Path(bank.tune).parent
    order = torch.argsort(bank.mass2).cpu().numpy()
    indices = np.sort(order[np.linspace(0, len(order) - 1, count, dtype=int)])
    copy_events(source / "events.hepmc3", directory / "events.hepmc3", count, indices=indices)
    positions = np.full(bank.amplitudes.shape[1], -1, dtype=int)
    positions[indices] = np.arange(count)
    samples = []
    for sample in bank.sample:
        selected = positions[np.asarray(sample["event_ids"], dtype=int)]
        accepted = selected >= 0
        samples.append({**sample, "event_ids": selected[accepted], "weights": np.asarray(sample["weights"])[accepted]})
    output = directory / "native.bin"
    cards = []
    configurations = [bank.initial | point for point in points]
    for index, parameters in enumerate(configurations):
        tune = directory / f"tune_{index}"
        driver.create_steering_card(param_space=parameters, tunename=str(tune), tune_default=str(bank.tune), cdir=cdir)
        card = directory / f"{index}.json"
        write_json_file(card, input_card(driver, run, tune))
        cards.append(card)
    evaluate_basis(cards=cards, events=directory / "events.hepmc3", output=output,
                   controls=controls, cdir=cdir, deadline=time.monotonic() + timeout,
                   log=directory / "native.json", indices=list(range(len(cards))))
    metadata = load_json_file(str(output) + ".json")
    if any(metadata["failures"]):
        raise ValueError("Native amplitude validation contains failed evaluations")
    native = np.fromfile(output, dtype=np.complex128).reshape(metadata["shape"])
    rows = []
    with torch.no_grad():
        for parameters, expected in zip(configurations, native, strict=True):
            actual = bank.amplitude(parameters)[indices].numpy()
            reference = bank.reference[indices].numpy()
            row = compare_amplitudes(actual, expected, reference)
            row["rates"] = compare_rates(actual, expected, reference, samples)
            if data is not None:
                row["histograms"] = histogram_statistics(bank.predict(parameters), data)
            rows.append(row)
    report = {"events": count, "event_indices": indices.tolist(), "parameters": configurations, "results": rows,
              "source": str(source), "screened": run.loopscreen == "true"}
    write_json_file(directory / "comparison.json", report)
    return report


# Quantify complex amplitude and intensity errors without hiding exact amplitude zeros
def compare_amplitudes(actual, expected, reference):
    error = np.sum(np.abs(actual - expected)**2, axis=1)
    native_intensity = np.sum(np.abs(expected)**2, axis=1)
    intensity = np.sum(np.abs(actual)**2, axis=1)
    nonzero = native_intensity > 0.0
    relative = np.divide(error, native_intensity, out=np.zeros_like(error), where=nonzero)
    relative[~nonzero & (error > 0.0)] = np.inf
    scale = np.maximum(native_intensity, reference)
    return {"max_relative_amplitude": float(np.sqrt(relative).max()),
        "max_source_relative_amplitude": float(np.sqrt(error / reference).max()),
        "max_scaled_intensity_error": float((np.abs(intensity - native_intensity) / scale).max()),
        "rms_relative_amplitude": float(np.sqrt(error.sum() / native_intensity.sum()))}


# Compare selected cross sections and effective sample sizes on common weighted events
def compare_rates(actual, expected, reference, samples):
    actual_ratio = np.sum(np.abs(actual)**2, axis=1) / reference
    native_ratio = np.sum(np.abs(expected)**2, axis=1) / reference
    rows = []
    for sample in samples:
        selected = np.asarray(sample["event_ids"], dtype=np.int64)
        accepted = selected < len(reference)
        selected = selected[accepted]
        weights = np.asarray(sample["weights"])[accepted]
        actual_weights = weights * actual_ratio[selected]
        native_weights = weights * native_ratio[selected]
        actual_sum, native_sum = actual_weights.sum(), native_weights.sum()
        actual_norm, native_norm = np.sum(actual_weights**2), np.sum(native_weights**2)
        rows.append({"selected_events": len(selected),
            "relative_rate_error": float((actual_sum - native_sum) / native_sum) if native_sum > 0.0 else None,
            "actual_effective_events": float(actual_sum**2 / actual_norm) if actual_norm > 0.0 else 0.0,
            "native_effective_events": float(native_sum**2 / native_norm) if native_norm > 0.0 else 0.0})
    return rows


# Measure MC coverage and effective counts in fitted bins using the complete prepared sample
def histogram_statistics(predictions, data):
    rows = []
    for subset, (prediction, measured) in enumerate(zip(predictions, data, strict=True)):
        for observable, entry in measured.items():
            mc, reference = prediction[observable]["hdata"], entry["hdata"]
            counts = np.asarray(array.to_numpy(mc.counts), dtype=float)
            errors = np.asarray(array.to_numpy(mc.errs), dtype=float)
            values = np.asarray(array.to_numpy(reference.counts_scaled), dtype=float)
            valid = mc.valid & reference.valid & np.isfinite(values) & (entry["fitw"] > 0.0)
            positive = valid & (values > 0.0)
            finite = np.isfinite(counts) & np.isfinite(errors)
            occupied = positive & finite & (counts > 0.0)
            # Raw weighted counts retain the effective sample size for density histograms
            effective = np.divide(counts, errors, out=np.zeros_like(counts), where=finite & (errors > 0.0))**2
            selected = effective[occupied]
            mc_error = np.asarray(array.to_numpy(mc.errs_scaled), dtype=float)
            data_error = np.asarray(array.to_numpy(reference.errs_scaled), dtype=float)
            comparable = occupied & np.isfinite(mc_error) & np.isfinite(data_error) & (data_error > 0.0)
            error_ratio = mc_error[comparable] / data_error[comparable]
            rows.append({"subset": subset, "observable": observable,
                "valid_bins": int(valid.sum()), "positive_data_bins": int(positive.sum()),
                "positive_data_bins_with_mc": int(occupied.sum()),
                "empty_mc_bins": np.flatnonzero(positive & finite & ~(counts > 0.0)).tolist(),
                "nonfinite_mc_bins": np.flatnonzero(valid & ~finite).tolist(),
                "min_effective_events_occupied_bin": float(selected.min()) if len(selected) else None,
                "median_effective_events_occupied_bin": float(np.median(selected)) if len(selected) else None,
                "mc_error_larger_than_data_bins": np.flatnonzero(comparable & (mc_error > data_error)).tolist(),
                "median_mc_to_data_error": float(np.median(error_ratio)) if len(error_ratio) else None})
    return rows


# Check one initialized sample using ordinary icetune steering and histogram input
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True)
    parser.add_argument("--tunesetup", required=True)
    parser.add_argument("--tune-default", help="Source tune for a Python tunesetup, read from a saved JSON tunesetup otherwise")
    parser.add_argument("--ampfit-config", required=True)
    parser.add_argument("--dataset", type=int, required=True)
    parser.add_argument("--sample", default="nominal")
    parser.add_argument("--points", required=True, help="JSON list of parameter dictionaries merged with the source point")
    parser.add_argument("--events", required=True, type=int)
    parser.add_argument("--timeout", required=True, type=float)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    root = Path.cwd()
    tunesetup = load_tunesetup(cdir=root, simdriver="GRANIITTI", name=args.tunesetup, tune_default=args.tune_default)
    controls = load_settings(args.ampfit_config)["bank"]
    driver = GraniittiDriver()
    directory = driver.amplitude_directory(cdir=str(root), run_name=args.run) / str(args.dataset) / args.sample
    closure = load_json_file(directory / "closure.json")
    error = max(float(closure[key]) for key in ("rms_relative_amplitude", "weighted_rms_relative_amplitude"))
    if not np.isfinite(error) or error > controls["closure_rtol"]:
        raise ValueError(f"Cannot validate a bank with failed source closure: {error:.6g}")
    driver.init_data(args.run, tunesetup.datacards, "default", str(root), pickle_dump=False)
    steering = {"tune_default": tunesetup.tune_default, "tunesetup_name": args.tunesetup}
    run = next(run for run in driver._run_cards(
        index=args.dataset, datacard=tunesetup.datacards[args.dataset], mc_steer=steering,
        tunename=tunesetup.tune_default, cdir=str(root)) if run.sample_name == args.sample)
    indices = driver._sample_set_indices(index=args.dataset, sample_name=args.sample)
    bank = AmplitudeBank(driver=driver, directory=directory, obs=[driver.obs[args.dataset][index] for index in indices],
                         pid=[driver.pid[args.dataset][index] for index in indices],
                         cuts=[driver.cuts[args.dataset][index] for index in indices], controls=controls)
    result = validate_points(driver=driver, bank=bank, run=run, points=load_json_file(args.points),
                             directory=args.output, cdir=str(root), controls=controls,
                             events=args.events, timeout=args.timeout,
                             data=[driver.data[args.dataset][index] for index in indices])
    for index, row in enumerate(result["results"]):
        summary = {key: value for key, value in row.items() if key != "histograms"}
        print(f"Point {index}: {summary}")


if __name__ == "__main__":
    main()

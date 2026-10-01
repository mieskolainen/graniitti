# Coherent event banks and shared icepack histogram predictions for icetune
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import errno
import json
import pickle
import shutil
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pyjson5
import torch
from torch.utils.checkpoint import checkpoint

from core.io import readers
from core.io.files import ensure_dir
from core.io.serialize import load_json_file, sha256_file, write_json_file
from core.plot import plot
from core.stats import hist as histogram
from core.tune.cache import file_records
from core.tune.drivers.graniitti.ampfit.coefficients import (
    AmplitudeSteering,
    ContinuumCoefficients,
    ResonanceCoefficients,
)
from core.tune.io import atomic_write_json
from core.tune.runtime import process as iceruntime


# Select source sample covariance while retaining the current predictions and occupancy
def source_mc(results):
    # Preserve predictions and occupancy while selecting the reference uncertainty
    def fixed(item):
        prediction = copy.copy(item["hdata"])
        prediction.stat_reference = item["source_statistics"]
        return {**item, "hdata": prediction}

    return {**results, "mc": [[{name: fixed(item) for name, item in subset.items()}
                               for subset in dataset] for dataset in results["mc"]]}


# Read the saved coefficient definitions identically during initialization and worker evaluation
def coefficients(driver, path, settings):
    if settings.get("interpolation") != "chebyshev":
        raise ValueError("Amplitude bank uses a different interpolation scheme and must be prepared again")
    common = dict(driver=driver, path=path,
                  **{key: settings[key] for key in ("model", "pid", "initial", "bounds", "nodes")})
    resonance = ResonanceCoefficients(**common, resonances=settings["resonances"], decay=settings["decay"])
    return settings, resonance, ContinuumCoefficients(**common)


# Materialize the effective source gencard using the driver's ordinary card field updates
def input_card(driver, run, tune):
    card = load_json_file(run.inputfile, loader=pyjson5.load)
    for override in run.overrides:
        field, _, value = override.partition("=")
        driver.update_target(field, json.loads(value), card)
    # Supplied events require no rejection sampling or maximum weight estimation
    card["GENERIC"].update(MODELPARAM=str(tune), INTEGRATOR=run.integrator, HIST=0, WEIGHTED=True)
    card["SCATTERING"].update(PROCESS=run.process_override or run.process, LOOPSCREEN=run.loopscreen == "true")
    return card


# Derive the prepared components from the ordinary source tune and active coordinates
def plan_bank(*, driver, run, tune, initial, bounds, pid, controls):
    card = input_card(driver, run, tune)
    process = card["SCATTERING"]["PROCESS"]
    model = process.split("[", 1)[0]
    if model not in {"MP", "XP", "GP"} or not process.startswith(f"{model}[RES+CON]<F>"):
        raise ValueError("ampfit requires a factorized MP, XP or GP RES+CON process")
    settings = dict(model=model, resonances=card["SCATTERING"]["RES"], pid=pid,
                    initial=initial, bounds=bounds, nodes=controls["nodes"], interpolation="chebyshev",
                    decay={key: controls[key] for key in ("decay_zero", "decay_steps")})
    return coefficients(driver, tune, settings)


# Record the native physics inputs independently of temporary directory names
def bank_request(*, settings, source, tune, cdir):
    source = {**source, "GENERIC": {**source["GENERIC"], "MODELPARAM": "tune"}}
    tune = Path(tune).resolve()
    return dict(settings=settings, source=source, tune=file_records(tune.rglob("*.json"), root=tune),
                binaries={name: sha256_file(Path(cdir) / "bin" / name) for name in ("gr", "ampfit")})


# Resolve one physics plan for complete-bank and event-sample reuse checks
class BankPlan:
    # Resolve coefficients and fingerprint the source physics once per requested bank
    def __init__(self, *, driver, run, tune, initial, bounds, pid, controls, cdir, nevents):
        self.driver, self.nevents = driver, nevents
        settings, self.resonance, self.continuum = plan_bank(
            driver=driver, run=run, tune=tune, initial=initial, bounds=bounds, pid=pid, controls=controls)
        self.request = bank_request(settings=settings, source=input_card(driver, run, tune), tune=tune, cdir=cdir)

    # Identify incompatible files and physics controls in one saved native bank
    def mismatch(self, directory):
        expected, nevents = self.request, self.nevents
        directory = Path(directory)
        missing = [name for name in ("request.json", "steering.json", "events.hepmc3", "amplitudes.bin.json")
                   if not (directory / name).is_file()]
        if missing:
            return ["missing " + name for name in missing]
        saved = load_json_file(directory / "request.json")
        metadata = load_json_file(directory / "amplitudes.bin.json")
        storage = Path(metadata.get("directory", directory))
        missing = [name + suffix for name in metadata.get("parts", ["amplitudes.bin"])
                   for suffix in ("", ".json", ".intensity", ".kinematics")
                   if not (storage / (name + suffix)).is_file()]
        if missing:
            return ["missing " + name for name in missing]
        differences = []
        if metadata["shape"][1] != nevents:
            differences.append(f"event count {metadata['shape'][1]} != {nevents}")
        if saved == expected:
            return differences
        # Changed exact coefficient bounds do not change the prepared native amplitudes
        expected = copy.deepcopy(expected)
        saved["settings"].pop("bounds")
        expected["settings"].pop("bounds")
        for field in ("binaries", "settings", "source", "tune"):
            old, new = saved[field], expected[field]
            if field == "tune":
                old, new = ({item["path"]: item for item in records} for records in (old, new))
            differences.extend(f"{field}.{key}" for key in sorted(old.keys() | new.keys()) if old.get(key) != new.get(key))
        if differences or saved != expected:
            return differences or ["request metadata"]
        _, old_resonance, old_continuum = coefficients(self.driver, directory / "tune", load_json_file(directory / "steering.json"))
        return [] if self.resonance.columns == old_resonance.columns and self.continuum.columns == old_continuum.columns else ["interpolation columns"]

    # Recover completed native shards after an interrupted merge without repeating event generation
    def recover(self, directory):
        directory = Path(directory)
        if (directory / "amplitudes.bin.json").is_file() or not (directory / "request.json").is_file():
            return
        if load_json_file(directory / "request.json") != self.request:
            return
        try:
            combine_bank_shards(directory=directory)
        except (FileNotFoundError, ValueError):
            return
        print(f"ampfit bank: recovered completed component shards in {directory}", flush=True)

    # Copy a compatible completed bank from an earlier run into the current initialization
    def reuse(self, directory, root):
        directory = Path(directory)
        closest = None
        for source in sorted(Path(root).glob("*/results/amplitude/*/*")):
            if source.name.endswith(("._old", ".partial")):
                continue
            if source.resolve() == directory.resolve():
                continue
            self.recover(source)
            differences = self.mismatch(source)
            if differences:
                if closest is None or len(differences) < len(closest[1]):
                    closest = source, differences
                continue
            if directory.exists():
                directory.rename(directory.with_name(f"{directory.name}.{time.time_ns()}._old"))
            temporary = directory.with_name(f"{directory.name}.{time.time_ns()}.partial")
            print(f"ampfit bank: copying {source} -> {directory}", flush=True)
            stage_bank(source, temporary)
            temporary.rename(directory)
            print(f"ampfit bank: copied {directory}", flush=True)
            return source
        if closest is None:
            print(f"ampfit bank: no previous bank directories found under {root}", flush=True)
        else:
            print(f"ampfit bank: closest previous bank {closest[0]} differs in {', '.join(closest[1])}", flush=True)
        return None

    # Reuse completed event sampling when native component preparation was interrupted
    def reusable_events(self, directory):
        directory = Path(directory)
        if not (directory / "request.json").is_file() or not (directory / "events.hepmc3").is_file():
            return False
        if load_json_file(directory / "request.json") != self.request:
            return False
        summary = readers.read_hepmc3_weight_summary(filename=str(directory / "events.hepmc3"))
        return summary["nevents"] == self.nevents


# Evaluate native amplitudes using the prepared source helicities when available
def evaluate_basis(*, cards, events, output, controls, cdir, deadline, log, indices, reference=None):
    command = [str(Path(cdir) / "bin/ampfit"), "--events", str(events),
               "--output", str(output), "--batch-events", str(controls["batch_events"]),
               "--closure-rtol", str(controls["closure_rtol"])]
    if reference is not None:
        command.extend(("--reference", str(reference)))
    for card in cards:
        command.extend(("--input", str(card)))
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise TimeoutError("ampfit initialization exceeded its time limit while preparing basis cards")
    if not iceruntime.execute_cmd(command, cwd=cdir, max_t=remaining, log_path=str(log)):
        raise RuntimeError(f"ampfit preparation failed, see {log}")
    metadata = load_json_file(str(output) + ".json")
    metadata["indices"] = indices
    write_json_file(str(output) + ".json", metadata)


# Combine native components in physical basis order without duplicating the source
def combine_basis(parts, output, *, materialize=True):
    metadata = [load_json_file(str(part) + ".json") for part in parts]
    first = metadata[0]
    kinematics = sha256_file(str(parts[0]) + ".kinematics")
    order = []
    for part_index, (part, entry) in enumerate(zip(parts, metadata, strict=True)):
        if entry["shape"][1:] != first["shape"][1:] or entry["layout"] != first["layout"]:
            raise ValueError("Native component ranges have incompatible dimensions")
        if sha256_file(str(part) + ".kinematics") != kinematics:
            raise ValueError("Native component ranges have different event kinematics")
        if len(entry["indices"]) != entry["shape"][0] or any(entry["failures"]):
            raise ValueError("Native component range is incomplete or contains failed amplitudes")
        order.extend((column, part_index, row) for row, column in enumerate(entry["indices"]))
    order.sort()
    if len({column for column, _, _ in order}) != len(order):
        raise ValueError("Native component ranges overlap")
    events, helicities = first["shape"][1:]
    ensure_dir(Path(output).parent)
    sizes = (("", events * helicities * np.dtype(np.complex128).itemsize),
             (".intensity", events * np.dtype(np.float64).itemsize))
    for suffix, size in sizes:
        if any(Path(str(part) + suffix).stat().st_size != entry["shape"][0] * size
               for part, entry in zip(parts, metadata, strict=True)):
            raise ValueError("Native component file size disagrees with its dimensions")
    merged = dict(first)
    for field in ("components", "failures"):
        merged[field] = [metadata[part][field][row] for _, part, row in order]
    merged["indices"] = [column for column, _, _ in order]
    merged["shape"] = [len(order), events, helicities]

    if not materialize:
        return merged
    required = len(order) * sum(size for _, size in sizes) + Path(str(parts[0]) + ".kinematics").stat().st_size
    try:
        for suffix, size in sizes:
            partial = Path(str(output) + suffix + ".partial")
            with partial.open("wb") as destination:
                for _, part_index, row in order:
                    # Preserve filesystem errno when a quota or storage limit interrupts a write
                    destination.write(np.memmap(str(parts[part_index]) + suffix, mode="r", dtype=np.uint8,
                                                offset=row * size, shape=(size,)))
        shutil.copyfile(str(parts[0]) + ".kinematics", str(output) + ".kinematics.partial")
        for suffix in ("", ".intensity", ".kinematics"):
            Path(str(output) + suffix + ".partial").replace(str(output) + suffix)
    except OSError as exc:
        exc.add_note(f"Amplitude bank merge into {output} requires {required} bytes in addition to the retained shards. "
                     "Check destination free space and disk quota before retrying")
        raise
    atomic_write_json(str(output) + ".json", merged)


# Prepare shared events and compute the source amplitude once before the basis jobs
def prepare_sample(*, driver, directory, run, tune, initial, bounds, pid, events, controls, cdir, deadline, shards=1):
    directory = Path(directory).resolve()
    ensure_dir(directory, exist_ok=False)
    shutil.copytree(tune, directory / "tune")
    shutil.copyfile(events, directory / "events.hepmc3")
    settings, resonance, continuum = plan_bank(driver=driver, run=run, tune=directory / "tune", initial=initial,
        bounds=bounds, pid=pid, controls=controls)
    write_json_file(directory / "steering.json", settings)
    source = input_card(driver, run, directory / "tune")
    write_json_file(directory / "request.json", bank_request(
        settings=settings, source=source, tune=directory / "tune", cdir=cdir))
    card = directory / "source.json"
    write_json_file(card, source)
    output = directory / "source.bin"
    evaluate_basis(cards=[card], events=events, output=output, controls=controls, cdir=cdir,
                   deadline=deadline, log=directory / "source.log", indices=[0])
    metadata = load_json_file(str(output) + ".json")
    if any(metadata["failures"]):
        raise ValueError("Source amplitude contains failed physical evaluations")
    metadata.update(shards=shards, columns=len(resonance.columns) + len(continuum.columns))
    write_json_file(str(output) + ".json", metadata)


# Interleave physical components across jobs using one local copy of the source tune
def write_basis_cards(*, driver, directory, temporary, run, cdir, shard=0, shards=1):
    directory, temporary = Path(directory).resolve(), Path(temporary).resolve()
    ensure_dir(temporary, exist_ok=False)
    tune = temporary / "tune"
    shutil.copytree(directory / "tune", tune)
    settings, resonance, continuum = coefficients(driver, tune, load_json_file(directory / "steering.json"))
    source = input_card(driver, run, tune)
    cards, indices = [], []
    columns = [*resonance.columns, *continuum.columns]
    for index in range(shard, len(columns), shards):
        column = columns[index]
        parameters = dict(column["parameters"])
        is_resonance = "resonance" in column
        if not is_resonance:
            parameters[f"REGGE|PARAM_CON:{settings['model']}"] = {
                **resonance.general["PARAM_CON"][settings["model"]], continuum.channel: [column["exchange_pair"]]}
        basis_tune = temporary / str(index)
        driver.create_steering_card(param_space=parameters, tunename=str(basis_tune), tune_default=str(tune), cdir=cdir)
        card = input_card(driver, run, basis_tune)
        card["SCATTERING"]["PROCESS"] = source["SCATTERING"]["PROCESS"].replace("[RES+CON]", "[RES]" if is_resonance else "[CON]", 1)
        card["SCATTERING"]["RES"] = [column["resonance"]] if is_resonance else []
        filename = temporary / f"{index}.json"
        write_json_file(filename, card)
        cards.append(filename)
        indices.append(index + 1)
    return cards, indices


# Evaluate component groups within the allocated CPU budget using a common source
def evaluate_components(*, cards, indices, reference, events, output, temporary, controls, cdir, deadline, processes):
    common = dict(events=events, controls=controls, cdir=cdir, deadline=deadline, reference=reference)
    count = min(processes, len(cards))
    if count == 0:
        ensure_dir(Path(output).parent)
        metadata = load_json_file(reference)
        metadata.update(shape=[0, *metadata["shape"][1:]], indices=[], components=[], failures=[])
        for suffix in ("", ".intensity"):
            Path(str(output) + suffix).write_bytes(b"")
        shutil.copyfile(str(reference).removesuffix(".json") + ".kinematics", str(output) + ".kinematics")
        write_json_file(str(output) + ".json", metadata)
    elif count == 1:
        evaluate_basis(cards=cards, indices=indices, output=output, log=temporary / "ampfit.json", **common)
    else:
        parts = [temporary / f"amplitudes_{index}.bin" for index in range(count)]
        with ThreadPoolExecutor(max_workers=count) as pool:
            futures = [pool.submit(evaluate_basis, cards=cards[index::count], indices=indices[index::count],
                                   output=part, log=temporary / f"ampfit_{index}.json", **common)
                       for index, part in enumerate(parts)]
            for future in futures:
                future.result()
        combine_basis(parts, output)


# Prepare complex amplitudes with the same source and basis stages as distributed jobs
def prepare_bank(*, driver, directory, temporary, run, tune, initial, bounds, pid, events, controls, cdir, deadline,
                 processes=1):
    prepare_sample(driver=driver, directory=directory, run=run, tune=tune, initial=initial, bounds=bounds,
                   pid=pid, events=events, controls=controls, cdir=cdir, deadline=deadline)
    prepare_bank_shard(driver=driver, directory=directory, temporary=temporary, run=run, controls=controls,
                       cdir=cdir, deadline=deadline, shard=0, shards=1, processes=processes)
    combine_bank_shards(directory=directory)


# Evaluate one basis shard on local events and publish its metadata after its binary files
def prepare_bank_shard(*, driver, directory, temporary, run, controls, cdir, deadline, shard, shards, processes=1):
    directory, temporary = Path(directory).resolve(), Path(temporary).resolve()
    reference = directory / "source.bin.json"
    source = load_json_file(reference)
    if source["shards"] != shards:
        raise ValueError("Amplitude shard count differs from the prepared sample")
    cards, indices = write_basis_cards(driver=driver, directory=directory, temporary=temporary, run=run, cdir=cdir,
                                      shard=shard, shards=shards)
    events = temporary / "events.hepmc3"
    _, count, helicities = source["shape"]
    native_bytes = len(cards) * count * (helicities * np.dtype(np.complex128).itemsize + np.dtype(np.float64).itemsize)
    native_bytes += count * 5 * np.dtype(np.float64).itemsize
    required = native_bytes * (2 if processes > 1 else 1) + (directory / "events.hepmc3").stat().st_size
    if required > shutil.disk_usage(temporary).free:
        raise OSError(errno.ENOSPC, f"Amplitude shard requires {required} free scratch bytes before evaluation")
    shutil.copyfile(directory / "events.hepmc3", events)
    output = temporary / "amplitudes.bin"
    evaluate_components(cards=cards, indices=indices, reference=reference, events=events, output=output,
                        temporary=temporary, controls=controls, cdir=cdir, deadline=deadline, processes=processes)
    part_dir = directory / "parts"
    ensure_dir(part_dir)
    for suffix in ("", ".intensity", ".kinematics", ".json"):
        part = part_dir / f"part_{shard}.bin{suffix}"
        partial = part.with_name(part.name + ".partial")
        shutil.copyfile(str(output) + suffix, partial)
        partial.replace(part)


# Combine the source and every basis shard in the original component order
def combine_bank_shards(*, directory):
    directory = Path(directory).resolve()
    source = directory / "source.bin"
    metadata = load_json_file(str(source) + ".json")
    parts = [source, *(directory / "parts" / f"part_{shard}.bin" for shard in range(metadata["shards"]))]
    missing = [str(part) + ".json" for part in parts if not Path(str(part) + ".json").is_file()]
    if missing:
        raise FileNotFoundError(f"Amplitude bank shards are incomplete: {missing}")
    indices = [index for part in parts for index in load_json_file(str(part) + ".json")["indices"]]
    if sorted(indices) != list(range(metadata["columns"] + 1)):
        raise ValueError("Amplitude bank basis components are incomplete")
    metadata = combine_basis(parts, directory / "amplitudes.bin", materialize=False)
    metadata["parts"] = [str(part.relative_to(directory)) for part in parts]
    atomic_write_json(directory / "amplitudes.bin.json", metadata)


# Stage small bank inputs while retaining immutable shared component paths
def stage_bank(source, destination, *, predictions=True):
    source, destination = Path(source), Path(destination)
    metadata = load_json_file(source / "amplitudes.bin.json")
    ignored = ["parts", "amplitudes.bin*", "*.partial", "*._old", "*.log", "source.bin*"]
    if not predictions:
        ignored.append("finalized.pkl")
    shutil.copytree(source, destination, ignore=shutil.ignore_patterns(*ignored))
    metadata["directory"] = str(Path(metadata.get("directory", source)).resolve())
    atomic_write_json(destination / "amplitudes.bin.json", metadata)


# Package only small bank inputs and explicit shared component references
def bank_files(directory):
    files = []
    for path in sorted(Path(directory).glob("*/*/amplitudes.bin.json")):
        if any(part.endswith(("._old", ".partial")) for part in path.relative_to(directory).parts):
            continue
        metadata = load_json_file(path)
        if "directory" not in metadata:
            atomic_write_json(path, {**metadata, "directory": str(path.parent.resolve())})
        files.extend(path.parent / name for name in ("request.json", "steering.json", "events.hepmc3", "amplitudes.bin.json"))
        files.extend((path.parent / "tune").rglob("*.json"))
    return files


# Load predictions only after the complete sample has been published
def load_finalized(directory):
    with (Path(directory) / "finalized.pkl").open("rb") as stream:
        result = pickle.load(stream)
    return result["mc"], result["mcdata"]


# Merge and validate one sample locally before publishing its completed predictions
def finalize_bank(*, directory, temporary, options, parameters, covariance, ready):
    directory, temporary = Path(directory), Path(temporary)
    completed = directory / "finalized.pkl"
    if ready and completed.is_file():
        print(f"ampfit finalize: reusing {directory}", flush=True)
        return
    if not ready:
        combine_bank_shards(directory=directory)
    print(f"ampfit finalize: closure and histograms {directory}", flush=True)
    with torch.no_grad():
        bank = AmplitudeBank(**{**options, "directory": directory})
        bank.validate_source(bank.amplitude(bank.initial), options["controls"]["closure_rtol"], directory / "closure.json")
        result = dict(mc=bank.predict(parameters), mcdata=bank.sample if covariance else None)
    # Publish completion only after all source closure and histogram checks succeed
    with completed.with_suffix(".pkl.partial").open("wb") as stream:
        pickle.dump(result, stream, protocol=pickle.HIGHEST_PROTOCOL)
    completed.with_suffix(".pkl.partial").replace(completed)
    print(f"ampfit finalize: completed {directory}", flush=True)


# Map distributed components in physical order without creating a second bank on disk
class Basis:
    # Validate and map each component shard once before any amplitude evaluation
    def __init__(self, directory, metadata):
        self.shape, self.dtype = tuple(metadata["shape"]), torch.complex128
        self.rows = [None] * self.shape[0]
        for name in metadata["parts"]:
            path = Path(metadata.get("directory", directory)) / name
            entry = load_json_file(str(path) + ".json")
            shape = tuple(entry["shape"])
            if (shape[1:] != self.shape[1:] or entry["layout"] != metadata["layout"]
                    or len(entry["indices"]) != shape[0] or any(entry["failures"])
                    or path.stat().st_size != int(np.prod(shape)) * np.dtype(np.complex128).itemsize):
                raise ValueError("Invalid amplitude bank shard dimensions or components")
            if not shape[0]:
                continue
            values = torch.from_numpy(np.memmap(path, mode="c", dtype=np.complex128, shape=shape))
            for row, column in enumerate(entry["indices"]):
                if not 0 <= column < len(self.rows) or self.rows[column] is not None:
                    raise ValueError("Amplitude bank shards have invalid or overlapping components")
                self.rows[column] = values[row]
        if any(row is None for row in self.rows):
            raise ValueError("Amplitude bank basis components are incomplete")

    # Read only the requested components and events from the mapped shards
    def __getitem__(self, index):
        components, events = index if isinstance(index, tuple) else (index, slice(None))
        rows = self.rows[components]
        if isinstance(components, int):
            return rows[events]
        if not rows:
            return torch.empty((0, *self.shape[1:]), dtype=self.dtype)[:, events]
        return torch.stack([row[events] for row in rows])


# Contract the resonance and continuum components with identical local and distributed physics
def contract_bank(resonance, continuum, basis, parameters, mass2, transfer, daughters, virtuality):
    steering = AmplitudeSteering(resonance.driver, resonance.cards, parameters)
    split = len(resonance.columns) + 1
    return (resonance.evaluate(steering, mass2, daughters, transfer, basis[1:split])
            + continuum.evaluate(steering, mass2, virtuality, transfer, basis[split:]))


# Hold physical complex amplitudes and the common selected event sample
class AmplitudeBank:
    # Validate the prepared amplitude bank before constructing any fit predictions
    def __init__(self, *, driver, directory, obs, pid, cuts, controls, density=False, scales=None, covariance_mode="diagonal", workers=None):
        directory = Path(directory)
        settings, resonance, continuum = coefficients(driver, directory / "tune", load_json_file(directory / "steering.json"))
        metadata = load_json_file(directory / "amplitudes.bin.json")
        storage = Path(metadata.get("directory", directory))
        path = storage / "amplitudes.bin"
        self.workers = workers
        shape = tuple(metadata["shape"])
        if len(shape) != 3 or min(shape) < 1 or metadata["layout"] != "component,event,helicity,re_im":
            raise ValueError("Invalid GRANIITTI amplitude bank dimensions or layout")
        if any(metadata["failures"]):
            raise ValueError("Amplitude bank contains failed physical evaluations")
        expected_bytes = int(np.prod(shape)) * np.dtype(np.complex128).itemsize
        if "parts" not in metadata and path.stat().st_size != expected_bytes:
            raise ValueError("Amplitude bank file size disagrees with its dimensions")
        # Copy on write maps retain shared file pages without exposing a read-only Tensor
        self.amplitudes = (Basis(directory, metadata) if "parts" in metadata else
                           torch.from_numpy(np.memmap(path, mode="c", dtype=np.complex128, shape=shape)))
        source = storage / "source.bin" if "parts" in metadata else path
        self.reference = torch.from_numpy(np.fromfile(str(source) + ".intensity", dtype=np.float64, count=shape[1]))
        if not torch.all(torch.isfinite(self.reference)) or torch.any(self.reference <= 0.0):
            raise ValueError("Reference event density must be finite and positive over the bank sample")
        self.sample = readers.read_hepmc3(str(directory / "events.hepmc3"), obs=obs, pid=pid, cuts=cuts)
        if any(len(sample["sample_weights"]) != shape[1] for sample in self.sample):
            raise ValueError("HepMC3 events and amplitude bank must describe the same sample")
        self.obs = obs
        self.covariance_mode = covariance_mode
        self.density = density
        self.scales = [1.0] * len(obs) if scales is None else scales
        self.source_statistics = None

        kinematics = torch.as_tensor(np.fromfile(str(source) + ".kinematics").reshape(-1, 5).copy())
        self.mass2 = kinematics[:, 0]
        self.transfer = kinematics[:, 1:3]
        if len(self.mass2) != self.amplitudes.shape[1]:
            raise ValueError("Amplitude bank kinematics and events have different lengths")
        self.daughters = [entry["mass"] for entry in metadata["daughters"]]
        # Extract the two meson vertex form factors using the smaller t or u virtuality
        self.virtuality = 2 * (self.daughters[0]**2 - kinematics[:, 3:].amax(1)).clamp_min(0)
        self.tune = directory / "tune"
        self.initial = settings["initial"]
        self.batch_events = controls["batch_events"]
        self.resonance, self.continuum = resonance, continuum
        if len(self.resonance.columns) + len(self.continuum.columns) + 1 != self.amplitudes.shape[0]:
            raise ValueError("Amplitude bank columns disagree with their physical steering")
        if controls["mc_stat"] == "source":
            reference = self.weighted_histograms(torch.ones_like(self.reference))
            self.source_statistics = [{name: histogram.numpy_output(item["hdata"]) for name, item in subset.items()}
                for subset in reference]

    # Evaluate a standard active parameter configuration through the shared coefficient decoder
    def coefficients(self, parameters):
        steering = AmplitudeSteering(self.resonance.driver, self.resonance.cards, parameters)
        return torch.cat((self.resonance.evaluate(steering, self.mass2, self.daughters, self.transfer),
                          self.continuum.evaluate(steering, self.mass2, self.virtuality, self.transfer)), dim=1)

    # Contract production rows before their common resonance and continuum factors
    def amplitude(self, parameters):
        if self.workers:
            from core.tune.drivers.graniitti.ampfit.distributed import remote_amplitude

            return remote_amplitude(parameters, self.workers)
        amplitudes = []
        for start in range(0, self.amplitudes.shape[1], self.batch_events):
            stop = start + self.batch_events
            amplitudes.append(checkpoint(self._batch, parameters, start, stop, use_reentrant=False)
                              if torch.is_grad_enabled() else self._batch(parameters, start, stop))
        return torch.cat(amplitudes)

    # Read one event batch inside the differentiable checkpoint so components can be released
    def _batch(self, parameters, start, stop):
        events = slice(start, stop)
        return contract_bank(self.resonance, self.continuum, self.amplitudes[:, events], parameters,
            self.mass2[events], self.transfer[events], self.daughters, self.virtuality[events])

    # Contract coherent basis amplitudes before summing physical external helicities
    def coherent(self, coefficients):
        amplitudes = []
        coefficients = coefficients.to(self.amplitudes.dtype)
        for start in range(0, self.amplitudes.shape[1], self.batch_events):
            stop = start + self.batch_events
            basis = self.amplitudes[1:, start:stop].to(coefficients.device)
            weights = coefficients[:, None] if coefficients.ndim == 1 else coefficients[start:stop].T
            amplitudes.append(torch.einsum("cnh,cn->nh", basis, weights))
        return torch.cat(amplitudes, dim=0)

    # Sum the coherent intensity over physical external helicities
    def intensity(self, coefficients):
        amplitude = self.coherent(coefficients)
        return amplitude.real.square().sum(1) + amplitude.imag.square().sum(1)

    # Require the reconstructed source amplitude to agree before publishing a fit bank
    def validate_source(self, amplitude, tolerance, output):
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ValueError("Amplitude closure tolerance must be finite and positive")
        with torch.no_grad():
            difference = amplitude - self.amplitudes[0]
            error = difference.abs().square().sum(1)
            relative = torch.sqrt(error / self.reference)
            weights = torch.as_tensor(self.sample[0]["sample_weights"], dtype=error.dtype).abs()
            result = {"max_relative_amplitude": float(relative.max()),
                      "rms_relative_amplitude": float(torch.sqrt(error.sum() / self.reference.sum())),
                      "weighted_rms_relative_amplitude": float(torch.sqrt((weights * relative.square()).sum() / weights.sum())),
                      "tolerance": tolerance}
        if not all(np.isfinite(value) for value in result.values()):
            raise ValueError("Reconstructed amplitude closure is nonfinite")
        write_json_file(output, result)
        error = max(result["rms_relative_amplitude"], result["weighted_rms_relative_amplitude"])
        if error > tolerance:
            raise ValueError(f"Amplitude bank fails RMS source closure: {error:.4g} > {tolerance:.4g}")
        return result

    # Fill the normal icepack observables, cuts, cross sections, MC errors and density normalization
    def predict(self, parameters):
        amplitude = self.amplitude(parameters)
        intensity = amplitude.real.square().sum(1) + amplitude.imag.square().sum(1)
        ratio = intensity / self.reference.to(intensity.device)
        predictions = self.weighted_histograms(ratio)
        if self.source_statistics is not None:
            for subset, reference in zip(predictions, self.source_statistics, strict=True):
                for name, item in subset.items():
                    item["source_statistics"] = reference[name]
        return predictions

    # Fill icepack histograms with event density ratios relative to the source sample
    def weighted_histograms(self, ratio):
        predictions = []
        for sample, obs, scale in zip(self.sample, self.obs, self.scales, strict=True):
            weights = torch.as_tensor(sample["weights"], dtype=ratio.dtype, device=ratio.device)
            selected = torch.as_tensor(sample["event_ids"], dtype=torch.int64, device=ratio.device)
            weights = weights * ratio[selected]
            cross_section = sample["sample_xsection_pb"] / np.sum(sample["sample_weights"]) * weights.sum()
            predictions.append(plot.histmc({**sample, "weights": weights, "xsection_pb": cross_section}, obs,
                    density=self.density, density_uncertainty="shape" if self.density else "scaled",
                    covariance_mode=self.covariance_mode, scale=scale, color=plot.colors(0), label="GRANIITTI"))
        return predictions

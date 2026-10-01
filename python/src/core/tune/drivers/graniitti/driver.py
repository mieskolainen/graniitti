# GRANIITTI simulation driver for icetune
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import copy
import hashlib
import json
import math
import multiprocessing
import os
import pathlib
import pickle
import re
import shutil
import signal
import tempfile
import time
import traceback
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from types import ModuleType

import numpy as np
import pyjson5 as json5

from core import resource
from core.io import jsonref, readers
from core.io import logger as log
from core.io import steering as dataset_steering
from core.io.files import ensure_dir
from core.io.serialize import json_safe, load_json_file, write_json_file
from core.numerics import array
from core.plot import loss, plot
from core.stats import cov, hist
from core.stats import objective as icecost
from core.tune import core as icetune
from core.tune import push as icetune_push
from core.tune.drivers.base import SimulatorDriver
from core.tune.drivers.graniitti import runtime as graniitti_runtime
from core.tune.drivers.graniitti.tunesetup import eikonal
from core.tune.drivers.graniitti.tunesetup.domains import angular_row_name, mp_spin_sectors, tp_con_names, tp_res_names
from core.tune.parameters import tools
from core.tune.runtime import process as iceruntime

logger = log.get_logger(__name__)

VGRID_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class _RunCard:
    """Resolved GRANIITTI paths and steering values for one dataset run"""

    gridfile: str
    inputfile: str
    integrator: str
    integrator_override: str | None
    output: str
    loopscreen: str
    process: str
    process_override: str | None
    sample_name: str
    sample_scale: float
    overrides: tuple[str, ...]
    beam_energy: tuple[float, float]
    xsmode: str


# Describe one complete optimizer coordinate group and its physical card values
@dataclass(frozen=True)
class _Coordinates:
    keys: tuple[str, ...]
    targets: tuple[str, ...]
    values: tuple[float, ...]
    encode: Callable
    amplitude: complex | None = None


# Decode complete coordinate groups into their physical card targets
def _decode_coordinates(param: dict, groups: list[_Coordinates]) -> dict:
    values = dict(param)
    for group in groups:
        for key in group.keys:
            values.pop(key, None)
        values.update(zip(group.targets, group.values, strict=True))
    return values


# Encode card defaults through the same coordinate groups used for decoding
def _encode_coordinate_defaults(original: dict, groups: list[_Coordinates]) -> dict:
    values = copy.deepcopy(original)
    for group in groups:
        missing = set(group.targets) - values.keys()
        if missing:
            raise ValueError(f"Missing original values for {sorted(missing)}")
        physical = [float(values.pop(key)) for key in group.targets]
        values.update(zip(group.keys, group.encode(physical), strict=True))
    return values


# Preserve structured GRANIITTI subprocess diagnostics across trial layers
class GraniittiCommandError(RuntimeError):
    # Initialize one simulator command error from its persistent failure payload
    def __init__(self, failure: dict):
        self.failure = copy.deepcopy(failure)
        super().__init__(str(failure["summary"]))


# Compute the time left in one complete GRANIITTI computation
def _remaining_time(deadline: float) -> float:
    remaining = float(deadline) - time.monotonic()
    if remaining <= 0.0:
        raise TimeoutError("GRANIITTI computation exceeded max_t")
    return remaining


# Replace the complete process tag before the final-state arrow
def _swap_process_tag(process: str, swap_process: str) -> str:
    if re.fullmatch(r"(?:MP|XP|GP|TP)\[(?:RES|CON|RES\+CON)\]<[FC]>", swap_process) is None:
        raise ValueError("GRANIITTI SWAP_PROCESS must be a full Pomeron process tag such as 'MP[RES+CON]<F>'")
    _, separator, final_state = process.partition("->")
    if not separator:
        raise ValueError("GRANIITTI SWAP_PROCESS requires a final-state arrow")
    return f"{swap_process} {separator}{final_state}"


# Compute one stable filename token for a beam energy
def _beam_energy_token(value: object) -> str:
    energy = float(value)
    if not math.isfinite(energy) or energy <= 0.0:
        raise ValueError("GRANIITTI beam energies must be finite and positive")
    return format(energy, ".12g").replace("+", "").replace("-", "m").replace(".", "p")


# Compute an input-card value after applying one exact dataset sample override
def _sample_steering_value(steering: dict, parameters: dict, path: str, *, default: object = None) -> object:
    if path in parameters:
        return parameters[path]
    value: object = steering
    for key in path.split("."):
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


# Apply an optional tunesetup density override without changing the source dataset
def _apply_force_density(dataset: dict, datacard: dict) -> dict:
    force_density = datacard.get("force_density", False)
    if not isinstance(force_density, bool):
        raise TypeError("GRANIITTI ICETUNE force_density must be boolean")
    if not force_density:
        return dataset
    output = copy.deepcopy(dataset)
    output["fit"]["normalization"] = "unit_density"
    return output


# Select one elastic prediction per kinematic configuration for the fitted eikonal model
def _apply_elastic_model(dataset: dict, datacard: dict) -> dict:
    model = datacard.get("elastic_model")
    if model is None:
        return dataset
    if model not in {"single", "double", "triple"}:
        raise ValueError("GRANIITTI icetune elastic_model must be single, double or triple")
    output = copy.deepcopy(dataset)
    # The triple model uses the single model samples as kinematic templates
    source = "single" if model == "triple" else model
    key = "GENERAL.json:PARAM_SOFT.active_model"
    output["samples"] = [sample for sample in output["samples"] if sample["parameters"].get(key) == source]
    if not output["samples"]:
        raise ValueError(f"GRANIITTI icetune elastic_model requires {source} model elastic samples")
    names = {sample["name"] for sample in output["samples"]}
    for entry in output["sets"]:
        selected = dataset_steering.set_sample_names(entry) or tuple(names)
        entry.pop("sample", None)
        entry["samples"] = [name for name in selected if name in names]
        if len(entry["samples"]) != 1:
            raise ValueError("GRANIITTI icetune elastic_model requires one elastic sample per data set")
    for sample in output["samples"]:
        sample["parameters"][key] = model
        sample["label"] = f"GRANIITTI X[EL] ({model})"
    return output


# Build a sample-specific grid identity including both beam energies
def _vgrid_card_identity(steering: dict, parameters: dict | None = None) -> tuple[str, tuple[float, float]]:
    parameters = {} if parameters is None else parameters
    energy = _sample_steering_value(steering, parameters, "SCATTERING.ENERGY")
    if not isinstance(energy, list) or len(energy) != 2:
        raise ValueError("GRANIITTI SCATTERING.ENERGY must contain two beam energies")
    beam_energy = (float(energy[0]), float(energy[1]))
    energy_label = "_".join(_beam_energy_token(value) for value in beam_energy)
    identity = steering if not parameters else {"parameters": parameters, "steering": steering}
    payload = json.dumps(identity, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    card_hash = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]
    return f"ENERGY_{energy_label}__CARD_{card_hash}", beam_energy


# Compute why serialized beam-energy metadata is incompatible
def _vgrid_energy_rebuild_reason(stored_energy: object, expected_energy: tuple[float, float]) -> str | None:
    if not isinstance(stored_energy, list) or len(stored_energy) != 2:
        return "BEAM_ENERGY metadata is invalid"
    try:
        energy_matches = all(math.isclose(float(stored), expected, rel_tol=1e-12, abs_tol=1e-12)
            for stored, expected in zip(stored_energy, expected_energy, strict=True))
    except (TypeError, ValueError):
        return "BEAM_ENERGY metadata is invalid"
    if not energy_matches:
        return f"beam energies {stored_energy} do not match {list(expected_energy)}"
    return None


# Compute why serialized integration-grid fields are incompatible
def _vgrid_field_rebuild_reason(griddata: dict, run: _RunCard) -> str | None:
    proposal_only = griddata.get("PROPOSAL_ONLY")
    if not isinstance(proposal_only, bool):
        return "PROPOSAL_ONLY metadata is invalid"
    expected_scalars = (("INTEGRATOR", run.integrator, str.upper),)
    for name, expected, normalizer in expected_scalars:
        stored = griddata.get(name)
        if normalizer(str(stored)) != expected:
            return f"{name}={stored} does not match {expected}"
    stored_loopscreen = str(griddata.get("LOOPSCREEN")).lower()
    if not proposal_only and stored_loopscreen != run.loopscreen.lower():
        return f"LOOPSCREEN={griddata.get('LOOPSCREEN')} does not match {run.loopscreen.lower()}"
    for name in (run.integrator, "STAT", "INTEGRATION", "MAX_WEIGHT"):
        if not isinstance(griddata.get(name), dict):
            return f"{name} metadata is invalid"
    return None


# Compute why one serialized integration grid must be rebuilt
def _vgrid_rebuild_reason(griddata: object, run: _RunCard) -> str | None:
    if not isinstance(griddata, dict):
        return "root value is not an object"
    if griddata.get("SCHEMA_VERSION") != VGRID_SCHEMA_VERSION:
        return "schema version is missing or incompatible"
    required = {"BEAM_ENERGY", "INTEGRATOR", "INTEGRATION", "MAX_WEIGHT", "LOOPSCREEN", "PROPOSAL_ONLY", "SQRTS",
        "STAT"}
    missing = sorted(required - set(griddata))
    if missing:
        return f"required fields are missing: {missing}"
    return _vgrid_energy_rebuild_reason(griddata.get("BEAM_ENERGY"), run.beam_energy) or _vgrid_field_rebuild_reason(
        griddata, run)


# Compute a human-readable subprocess termination reason
def _command_termination_reason(result: dict) -> str:
    returncode = result.get("returncode")
    status = str(result.get("status") or "failed")
    if isinstance(returncode, int) and returncode < 0:
        signum = -returncode
        try:
            signal_name = signal.Signals(signum).name
        except ValueError:
            signal_name = f"signal {signum}"
        return f"terminated by {signal_name} (signal {signum})"
    if isinstance(returncode, int):
        return f"exit code {returncode}"
    return status.replace("_", " ")


# Build one concise structured GRANIITTI command failure
def _graniitti_command_error(
    *, result: dict, stage: str, datacard_index: int, inputfile: str, gridfile: str, process: str, output: str
) -> GraniittiCommandError:
    root_cause = iceruntime.command_output_root_cause(result.get("output"))
    termination = _command_termination_reason(result)
    card_label = os.path.basename(inputfile)
    summary = f"GRANIITTI {stage.replace('_', ' ')} failed for datacard[{datacard_index}] {card_label} ({termination})"
    if root_cause:
        summary += f": {root_cause}"
    failure = {"argv": [str(value) for value in result.get("argv", [])], "datacard_index": int(datacard_index),
        "driver": "GRANIITTI", "duration_seconds": result.get("duration_seconds"),
        "error_type": "GraniittiCommandError", "gridfile": gridfile, "inputfile": inputfile, "output": output,
        "output_tail": iceruntime.command_output_tail(result.get("output")), "process": process,
        "returncode": result.get("returncode"), "root_cause": root_cause, "stage": stage,
        "status": result.get("status"), "summary": summary}
    return GraniittiCommandError(failure)


# Build one ordered GRANIITTI command from short-option pairs
def _gr_command(cdir: str, *options: tuple[str, object]) -> list[str]:
    command = [os.path.join(cdir, "bin/gr")]
    for name, value in options:
        command.extend((f"-{name}", str(value)))
    return command


def prepare_init_datacards(datacards: list[dict]) -> list[dict]:
    """
    Return bootstrap datacards for shared initialization.

    Initialization should preserve the logical xsmode/loopscreen choice of the
    original datacard and request proposal-only integration.
    """

    init_datacards = copy.deepcopy(datacards)
    for card in init_datacards:
        card["nevents"] = -1
    return init_datacards


# Configure the MC samples used for data statistical correlations
def prepare_mc_correlation_datacards(datacards: list[dict], *, event_count: int | None, weighting: str) -> list[dict]:
    if event_count is not None and event_count <= 0:
        raise ValueError("MC correlation event count must be positive")
    if weighting not in {"card", "weighted", "unweighted"}:
        raise ValueError(f"Unknown MC correlation weighting '{weighting}'")

    mc_correlation_datacards = copy.deepcopy(datacards)
    for card in mc_correlation_datacards:
        if event_count is not None:
            card["nevents"] = int(event_count)
        if int(card["nevents"]) <= 0:
            raise ValueError("MC correlation initialization requires events")
        if weighting != "card":
            card["weighted"] = weighting == "weighted"
    return mc_correlation_datacards


# Compute whether one GRANIITTI tune name is an explicit path
def _is_explicit_tune_path(tunename: str) -> bool:
    value = str(tunename)
    return os.path.isabs(value) or value.startswith(".") or os.sep in value


# Compute whether one name denotes a temporary GRANIITTI tune
def _is_temporary_tunename(tunename: str) -> bool:
    return str(tunename).startswith("TUNE_icetune")


# Resolve a GRANIITTI model tune or temporary tune directory
def _resolve_tune_dir(*, cdir: str, tunename: str) -> pathlib.Path:
    if _is_explicit_tune_path(tunename):
        return pathlib.Path(tunename) if os.path.isabs(tunename) else pathlib.Path(cdir) / tunename
    root = "tmp" if _is_temporary_tunename(tunename) else "modeldata"
    return pathlib.Path(cdir) / root / tunename


# Compute a flat output name that distinguishes explicit tune directories
def _tune_tag(*, cdir: str, tunename: str) -> str:
    if not _is_explicit_tune_path(tunename):
        return str(tunename)
    path = _resolve_tune_dir(cdir=cdir, tunename=tunename).resolve()
    digest = hashlib.sha256(str(path).encode()).hexdigest()[:12]
    return f"{path.name}_{digest}"


def resolve_gr_modelparam(*, cdir: str, tunename: str) -> str:
    """Return the MODELPARAM argument to pass to GRANIITTI."""

    if _is_temporary_tunename(tunename):
        return str(_resolve_tune_dir(cdir=cdir, tunename=tunename))
    return str(tunename)


# Remove one temporary GRANIITTI tuning card directory
def _cleanup_temporary_tune_dir(*, cdir: str, tunename: str) -> None:
    if not _is_temporary_tunename(tunename):
        return

    tune_dir = _resolve_tune_dir(cdir=cdir, tunename=tunename)
    shutil.rmtree(tune_dir, ignore_errors=True)
    for stage_dir in tune_dir.parent.glob(f".{tune_dir.name}.tmp.*"):
        shutil.rmtree(stage_dir, ignore_errors=True)


# Resolve file-backed dataset inputs which define cached data or bootstrap state
def _dataset_input_paths(
    *, cdir: pathlib.Path, dataset: dict, dataset_path: pathlib.Path, include_gencards: bool, strict_sets: bool
) -> set[pathlib.Path]:
    source = jsonref.JsonReader(json5.load)
    source.read(dataset_path)
    candidates = set(source.documents)

    # Add one resolved file-backed Python reference to the input set
    def add_python_reference(reference, package):
        resolved = dataset_steering.resolve_python_reference(reference, package=package, dataset_path=dataset_path, cdir=cdir)
        if resolved.endswith(".py"):
            path = pathlib.Path(resolved)
            candidates.add(path)
            if path.is_relative_to(cdir / "icepack"):
                candidates.update(path.parent.glob("*.py"))
                for parent in path.parents:
                    if parent == cdir:
                        break
                    candidates.update((parent / "_common").glob("*.py"))

    if include_gencards:
        for sample in dataset["samples"]:
            candidates.add(pathlib.Path(
                    dataset_steering.resolve_gencard_reference(sample["gencard"], dataset_path=dataset_path, cdir=cdir))
            )
    reader = dataset.get("reader")
    if reader:
        add_python_reference(reader, "core.io.hepdata_reader")
        candidates.update((cdir / "icepack/_common").glob("*.py"))
    entries = dataset["sets"] if strict_sets else dataset.get("sets", [])
    for entry in entries:
        for key, package in (("cuts", "core.analysis.cuts"), ("obs", "core.analysis.observables")):
            add_python_reference(entry[key], package)
        if not entry.get("data", True):
            continue
        histograms = entry["hist"] if strict_sets else entry.get("hist", [])
        for histogram in histograms:
            filename = pathlib.Path(dataset_steering.resolve_data_reference(
                histogram["file"], datapath=dataset["datapath"], dataset_path=dataset_path, cdir=cdir))
            candidates.add(filename)
            # Readers may also use normalization or correlation tables in the same HEPData entry
            if filename.suffix == ".json":
                candidates.update(filename.parent.glob("*.json"))
            for field in ("covariance", "flux"):
                if histogram.get(field):
                    candidates.add((filename.parent / histogram[field]).resolve())
    return candidates


# Build the penalty outputs for a failed simulation trial
def _penalty_trial_outputs(*, config: dict, trial_id: str, tunename: str, error: str, card_config: dict | None = None,
    failure: dict | None = None, error_type: str | None = None, stage: str | None = None,
    traceback_text: str | None = None,) -> dict:
    """Return a non-renderable penalty payload for a failed GRANIITTI trial evaluation."""

    penalty = float(np.inf)
    return {"config": copy.deepcopy(config), "card_config": copy.deepcopy(card_config), "cost_arr": None,
        "error": error, "error_type": error_type, "failure": copy.deepcopy(failure),
        "metrics": {"chi2": penalty, "gaussian": penalty, "ratio2": penalty, "wasserstein": penalty}, "ndf_arr": None, "results": None,
        "stage": stage, "traceback": traceback_text, "trial_id": trial_id, "tunename": tunename, "valid_arr": None,
        "weight_arr": None}


class GraniittiDriver(SimulatorDriver):
    BOOTSTRAP_SCHEMA_VERSION = 2
    CACHE_DIRS = ("eikonal", "sudakov", "vgrid")
    TRIAL_PREFIX = "TUNE_icetune"

    DEFAULT_TUNE = "TUNE0"
    AMPLITUDE_FIT = True

    # Compute the stable driver identifier
    @classmethod
    def driver_name(cls) -> str:
        return "GRANIITTI"


    # Recognize GRANIITTI hierarchical tune parameters
    @classmethod
    def matches_summary(cls, summary: dict) -> bool:
        config = summary.get("config")
        return isinstance(config, dict) and bool(config) and all("|" in str(key) for key in config)

    # Add GRANIITTI data covariance and MC correlation controls
    @classmethod
    def add_cli_arguments(cls, parser) -> None:
        parser.add_argument("--data_covariance_mode", choices=("full", "diagonal"), default="diagonal",
            help="use full data and MC covariance or their diagonals")
        parser.add_argument("--mc_correlation_events", type=int, default=None,
            help="MC correlation events per dataset, default uses each tune datacard")
        parser.add_argument("--mc_correlation_weighting", choices=("card", "weighted", "unweighted"), default="card",
            help="MC correlation sample weighting, default uses each tune datacard")

    # Validate GRANIITTI data covariance and MC correlation controls
    @classmethod
    def validate_cli_arguments(cls, parser, args) -> None:
        if args.algorithm == "ampfit":
            if not args.precompute:
                parser.error("ampfit requires shared amplitude initialization with --precompute true")
            if args.cost not in {"chi2", "gaussian", "ratio2"}:
                parser.error("ampfit requires a differentiable chi2, gaussian or ratio2 objective")
        if args.cost == "gaussian":
            if args.algorithm != "ampfit":
                parser.error("The gaussian objective currently requires ampfit")
            if args.cost_rho != "quadratic" or args.cost_avg != "sum":
                parser.error("The gaussian objective requires --cost_rho quadratic --cost_avg sum")
        if not args.precompute and args.data_covariance_mode == "full":
            parser.error("--data_covariance_mode full requires --precompute true")
        if args.data_covariance_mode == "full" and args.cost == "chi2" and args.cost_rho != "quadratic":
            parser.error("--data_covariance_mode full with --cost chi2 requires --cost_rho quadratic")
        if args.mc_correlation_events is not None and args.mc_correlation_events <= 0:
            parser.error("--mc_correlation_events must be positive")

    # Compute the default GRANIITTI external library root
    @classmethod
    def default_library_path(cls) -> str:
        return graniitti_runtime.default_library_path()

    # Build GRANIITTI run and MC correlation initialization steering
    def build_run_steering(self, args) -> dict:
        return super().build_run_steering(args) | {
            "ampfit": copy.deepcopy(args.ampfit_settings["bank"]) if args.algorithm == "ampfit" else None,
            "ampfit_reuse": getattr(args, "ampfit_reuse", True),
            "ampfit_reuse_root": getattr(args, "bank_shared_dir", None),
            "data_covariance_mode": args.data_covariance_mode, "mc_correlation_events": args.mc_correlation_events,
            "mc_correlation_weighting": args.mc_correlation_weighting}

    # Bind the source model directory while constructing parameter domains
    @classmethod
    def build_tunesetup(cls, *, config, cdir, tune_default):
        from core.tune.drivers.graniitti.tunesetup import build
        from core.tune.drivers.graniitti.tunesetup.domains import context
        with context(cdir=cdir, model_path=_resolve_tune_dir(cdir=cdir, tunename=tune_default),
                     settings=config["settings"]):
            return build(config)

    # Require every active amplitude coordinate to be represented before creating generator jobs
    def prepare_tunesetup(self, *, tunesetup, args) -> None:
        if getattr(args, "save_events", False) and any(int(card["nevents"]) <= 0 for card in tunesetup.datacards):
            raise ValueError("save_events requires positive event counts in the tuning datacards")
        if args.algorithm != "ampfit":
            return
        from core.tune.drivers.graniitti.ampfit.amplitude import BankPlan, plan_bank
        from core.tune.parameters.space import normalize_continuous_param_space, validate_initial_config

        bounds = normalize_continuous_param_space(tunesetup.param_space)
        controls = args.ampfit_settings["bank"]
        initial = self.get_initial_param(tunesetup.param_space, tunesetup.aux_param_space, args.cdir, args.tune_default)
        validate_initial_config(initial, tunesetup.param_space)
        initial.update(tunesetup.aux_param_space)
        self.init_data(args.run_name, tunesetup.datacards, args.obs_module, args.cdir, pickle_dump=False)
        if args.cost == "gaussian":
            self.validate_gaussian_data()
        source = pathlib.Path(args.cdir) / "runs" / "icetune" / args.run_name / "init" / "ampfit_preflight"
        if source.exists():
            source.rename(source.with_name(f"{source.name}.{time.time_ns()}._old"))
        self.create_steering_card(param_space=initial, tunename=str(source),
                                  cdir=args.cdir, tune_default=args.tune_default)
        unused = set(bounds)
        requests = []
        tunesetup.bank_samples = []
        for index, dataset in enumerate(self.datasets):
            if not dataset["active"]:
                continue
            if tunesetup.datacards[index]["nevents"] <= 0:
                raise ValueError("ampfit requires positive event counts in active tuning datacards")
            for run in self._run_cards(index=index, datacard=tunesetup.datacards[index], mc_steer=self.build_run_steering(args),
                                       tunename=args.tune_default, cdir=args.cdir):
                indices = self._sample_set_indices(index=index, sample_name=run.sample_name)
                pid = self.pid[index][indices[0]]
                if any(self.pid[index][set_index] != pid for set_index in indices):
                    raise ValueError("Amplitude observables routed to one sample must select the same final particles")
                _, resonance, continuum = plan_bank(
                    driver=self, run=run, tune=source, initial=initial, bounds=bounds, pid=pid, controls=controls)
                unused -= resonance.used | continuum.used
                requests.append((index, run, pid))
                tunesetup.bank_samples.append(dict(dataset=index, sample=run.sample_name,
                    nevents=tunesetup.datacards[index]["nevents"],
                    components=len(resonance.columns) + len(continuum.columns)))
        if unused:
            raise ValueError("Active parameters are not represented in the amplitude bank: " + ", ".join(sorted(unused)))
        if getattr(args, "preflight", False) and getattr(args, "ampfit_reuse", True) and not getattr(args, "init_force", False):
            directory = self.amplitude_directory(cdir=getattr(args, "bank_shared_dir", None) or args.cdir, run_name=args.run_name)
            print(f"ampfit bank search: {directory.parents[2]}", flush=True)
            for index, run, pid in requests:
                target = directory / str(index) / readers.clean_filename(run.sample_name)
                plan = BankPlan(driver=self, run=run, tune=source, initial=initial, bounds=bounds, pid=pid,
                                controls=controls, cdir=args.cdir, nevents=tunesetup.datacards[index]["nevents"])
                plan.recover(target)
                differences = plan.mismatch(target)
                if not differences:
                    print(f"ampfit bank: keeping compatible bank {target}", flush=True)
                    continue
                completed = target / "finalized.pkl"
                if completed.is_file():
                    completed.rename(completed.with_name(f"finalized.{time.time_ns()}.pkl._old"))
                print(f"ampfit bank: requested bank {target} differs in {', '.join(differences)}", flush=True)
                if plan.reuse(target, directory.parents[2]) is None:
                    print(f"ampfit bank: no compatible previous bank for dataset {index}, sample {run.sample_name}, preparation required", flush=True)

    # Build or apply the GRANIITTI worker runtime environment
    def runtime_environment(self, *, cdir: str, libdir: str, python_version: str, apply: bool = False
    ) -> dict[str, str]:
        return graniitti_runtime.environment(cdir=cdir, libdir=libdir, python_version=python_version, apply=apply)

    # Compute a temporary GRANIITTI tune name for one evaluation
    def trial_tunename(
        self, *, trial_id: str, node_id: str | None = None, pid: int | None = None, purpose: str = "trial") -> str:
        if purpose == "probe":
            trial_id = f"{trial_id}-{time.time_ns()}-{purpose}"
            purpose = "trial"
        return super().trial_tunename(trial_id=trial_id, node_id=node_id, pid=pid, purpose=purpose)

    # Initialize mutable per-process driver state
    def __init__(self):
        self.datasets = None
        self.dataset_paths = None
        self.data = None
        self.obs = None
        self.pid = None
        self.cuts = None
        self.data_covariance_payload = None
        self.wasserstein_replica_cache = None
        self.amplitude_banks = {}
        self.amplitude_path = None
        self.initialized = False

    # Require unit observable weights for a Gaussian likelihood before evaluating trials
    def validate_gaussian_data(self):
        for dataset in self.data:
            for subset in dataset:
                for item in subset.values():
                    weight = float(item["fitw"])
                    if weight > 0.0 and not np.isclose(weight, 1.0, rtol=0.0, atol=np.finfo(float).eps):
                        raise ValueError("The gaussian objective requires unit positive fit weights")

    # Prepare immutable cost state once before repeated trial evaluations
    def prepare_trial_runtime(self, param: dict) -> None:
        self.wasserstein_replica_cache = None
        if param.get("cost") == "gaussian" and self.data is not None:
            self.validate_gaussian_data()
        if param.get("cost") != "wasserstein" or self.data is None:
            return
        use_full_covariance = param.get("data_covariance_mode") == "full"
        covariance_payload = self.data_covariance_payload if use_full_covariance else None
        if use_full_covariance and covariance_payload is None:
            return
        self.wasserstein_replica_cache = icecost.build_wasserstein_replica_cache(
            data=self.data, covariance_payload=covariance_payload, rngseed=int(param.get("rngseed", 0)))

    # Compute unscreened initialization cards for reusable integration grids
    def prepare_init_datacards(self, datacards: list[dict]) -> list[dict]:
        return prepare_init_datacards(datacards)

    # Compute checksums for files which determine cached HEPData objects
    def _data_input_file_records(
        self, *, cdir: pathlib.Path, datacards: Iterable[dict], datasets: Iterable[dict] | None = None) -> list[dict]:
        from core.tune import cache as icetune_cache

        candidates: set[pathlib.Path] = set()
        supplied = list(datasets) if datasets is not None else None

        for index, card in enumerate(datacards):
            dataset_path = pathlib.Path(dataset_steering.resolve_dataset_reference(card["datacard"], cdir=cdir))
            dataset = (
                supplied[index] if supplied is not None else dataset_steering.load_dataset(str(dataset_path), cdir=cdir)[0]
            )
            # Add one file-backed reader, cut or observable definition
            candidates.update(_dataset_input_paths(
                    cdir=cdir, dataset=dataset, dataset_path=dataset_path, include_gencards=False, strict_sets=True))
        return icetune_cache.file_records(candidates, root=cdir)

    # Compute checksums for inputs which affect GRANIITTI initialization
    def _bootstrap_input_file_records(self, *, cdir: pathlib.Path, datacards: Iterable[dict], datasets: Iterable[dict],
        tune_default: str, tunesetup_path: str | None,) -> list[dict]:
        from core.tune import cache as icetune_cache

        candidates: set[pathlib.Path] = set()

        for relative in ["bin/gr", "VERSION.json", "modeldata"]:
            candidates.update(icetune_cache.files_below(cdir, relative))
        candidates.update(icetune_cache.files_below(cdir, str(_resolve_tune_dir(cdir=str(cdir), tunename=tune_default))))
        for card, dataset in zip(datacards, datasets, strict=False):
            if not isinstance(dataset, dict):
                continue
            dataset_path = pathlib.Path(dataset_steering.resolve_dataset_reference(card["datacard"], cdir=cdir))
            # Add one resolved file-backed Python reference to the fingerprint
            candidates.update(_dataset_input_paths(
                    cdir=cdir, dataset=dataset, dataset_path=dataset_path, include_gencards=True, strict_sets=False))
        if tunesetup_path:
            candidate = pathlib.Path(tunesetup_path).resolve()
            if candidate.is_file():
                candidates.add(candidate)
        return icetune_cache.file_records(candidates, root=cdir)

    # Collect the physics and runtime inputs independently of bank storage locations
    def bootstrap_inputs(
        self, *, cdir: str, datacards: list[dict], mc_steer: dict, runtime_sha256: str | None, tunesetup_path: str | None
    ) -> dict:
        from core.tune import cache as icetune_cache

        root = pathlib.Path(cdir).resolve()
        payload = {"cache_schema_version": self.BOOTSTRAP_SCHEMA_VERSION, "datacards": datacards,
            "input_files": self._bootstrap_input_file_records(cdir=root, datacards=datacards, datasets=self.datasets,
                tune_default=str(mc_steer["tune_default"]), tunesetup_path=tunesetup_path),
            "mc_steer": {key: value for key, value in mc_steer.items() if key != "ampfit_reuse_root"},
            "runtime_sha256": runtime_sha256}
        if mc_steer.get("data_covariance_mode") == "full":
            payload["data_covariance_schema"] = cov.DATA_COVARIANCE_SCHEMA
        if mc_steer.get("ampfit"):
            modules = ("tune/drivers/graniitti/driver", "tune/drivers/graniitti/ampfit/amplitude", "tune/drivers/graniitti/ampfit/coefficients",
                       "tune/parameters/tools", "numerics/array", "numerics/interp", "io/jsonref")
            payload["amplitude_files"] = icetune_cache.file_records([*(resource(f"{module}.py") for module in modules),
                 root / "bin/ampfit"], root=root)
        return copy.deepcopy(payload)

    # Compute the physics and runtime fingerprint for GRANIITTI initialization
    def bootstrap_fingerprint(
        self, *, cdir: str, datacards: list[dict], mc_steer: dict, runtime_sha256: str | None, tunesetup_path: str | None
    ) -> str:
        from core.tune import cache as icetune_cache

        return icetune_cache.json_fingerprint(self.bootstrap_inputs(cdir=cdir, datacards=datacards, mc_steer=mc_steer,
            runtime_sha256=runtime_sha256, tunesetup_path=tunesetup_path))

    # Compute regular files in the GRANIITTI reusable-cache directories
    def _bootstrap_cache_file_set(self, cdir: pathlib.Path) -> set[pathlib.Path]:
        from core.tune import cache as icetune_cache

        files = {path for directory in self.CACHE_DIRS for path in icetune_cache.files_below(cdir, directory)}
        return {path.resolve() for path in files if path.name != ".gitignore"}

    # Build or reuse one content-addressed GRANIITTI initialization cache
    def prepare_bootstrap(self, *, cdir: str, cache_base_url: str, datacards: list[dict], mc_steer: dict,
        runtime_sha256: str | None, tunesetup_path: str | None, init_force: bool, initialize: Callable[[], None],
        reusable_vgrids: Callable[[], Iterable[str | os.PathLike[str]]],
        reusable_outputs: Callable[[], Iterable[str | os.PathLike[str]]] | None = None,) -> dict:
        from core.tune import cache as icetune_cache

        if not cache_base_url:
            raise icetune_cache.PermanentConfigurationError("--bootstrap_cache_url is required for GRANIITTI init phase"
            )
        if mc_steer.get("ampfit") and not mc_steer.get("ampfit_reuse", True):
            cache_base_url = f"{cache_base_url.rstrip('/')}/rebuilds/{time.time_ns()}"
        root = pathlib.Path(cdir).resolve()
        inputs = self.bootstrap_inputs(cdir=str(root), datacards=datacards, mc_steer=mc_steer,
            runtime_sha256=runtime_sha256, tunesetup_path=tunesetup_path)
        fingerprint = icetune_cache.json_fingerprint(inputs)
        scratch_parent = root / "tmp"
        ensure_dir(scratch_parent)
        with tempfile.TemporaryDirectory(prefix="icetune-bootstrap-", dir=scratch_parent) as temporary:
            temporary_dir = pathlib.Path(temporary)
            if not init_force:
                existing = icetune_cache.reuse_content_archive(cache_base_url=cache_base_url, kind="graniitti",
                    fingerprint=fingerprint, temporary_dir=temporary_dir, schema_version=self.BOOTSTRAP_SCHEMA_VERSION,
                    label="GRANIITTI")
                if existing is not None:
                    return existing

            before = self._bootstrap_cache_file_set(root)
            initialize()
            after = self._bootstrap_cache_file_set(root)
            selected = {path for path in after - before if path.relative_to(root).parts[0] in {"eikonal", "sudakov"}}
            for label, outputs in (("vgrid", reusable_vgrids), ("initialization output", reusable_outputs)):
                for value in outputs() if outputs is not None else ():
                    path = pathlib.Path(value).resolve()
                    if not path.is_file():
                        raise icetune_cache.PermanentConfigurationError(f"Expected reusable {label} was not produced: {path}")
                    selected.add(path)

            records = icetune_cache.file_records(selected, root=root, allowed=graniitti_runtime.allowed_bootstrap_path)
            return icetune_cache.publish_content_archive(cache_base_url=cache_base_url, kind="graniitti",
                fingerprint=fingerprint, root=root, records=records, temporary_dir=temporary_dir,
                schema_version=self.BOOTSTRAP_SCHEMA_VERSION, label="GRANIITTI", manifest_extra={
                    "fingerprint_inputs": inputs, "reusable_vgrids": sorted(
                        record["path"] for record in records if record["path"].startswith("vgrid/"))})

    # Extract and verify one GRANIITTI bootstrap in worker-local scratch
    def stage_bootstrap(self, *, cdir: str, bootstrap: dict) -> None:
        graniitti_runtime.stage_bootstrap(cdir=cdir, bootstrap=bootstrap)
        self._bind_amplitude(bootstrap)

    # Resolve the bank location recorded by the verified initialization manifest
    def _bind_amplitude(self, bootstrap: dict) -> None:
        roots = {pathlib.PurePosixPath(*path.parts[:5]) for record in bootstrap.get("files", [])
                 if len((path := pathlib.PurePosixPath(record["path"])).parts) > 5
                 and path.parts[:2] == ("runs", "icetune") and path.parts[3:5] == ("results", "amplitude")}
        if len(roots) > 1:
            raise ValueError("An initialization manifest must contain one amplitude bank directory")
        self.amplitude_path = str(next(iter(roots))) if roots else None
        self.amplitude_banks.clear()

    # Load the data covariance pair from one staged GRANIITTI bootstrap
    def _load_bootstrap_covariance(self, *, bootstrap: dict, run_name: str, cdir: str, publish_current: bool) -> None:
        from core.tune import cache as icetune_cache

        covariance_paths = {pathlib.Path(record["path"]).name: os.path.join(cdir, record["path"])
            for record in bootstrap.get("files", [])
            if pathlib.Path(record.get("path", "")).name in {"data_covariance.npz", "data_covariance.json"}}
        expected = {"data_covariance.npz", "data_covariance.json"}
        if set(covariance_paths) != expected:
            raise icetune_cache.PermanentConfigurationError("GRANIITTI bootstrap has no data covariance pair")
        self.load_data_covariance(
            npz_path=covariance_paths["data_covariance.npz"], json_path=covariance_paths["data_covariance.json"])
        if publish_current:
            self.save_data_covariance(run_name=run_name, cdir=cdir)

    # Build GRANIITTI callbacks for a selected backend initialization
    def initialize_backend(self, *, args, tunesetup, mc_steer: dict) -> dict:
        if args.backend != "ray":
            raise ValueError(f'GRANIITTI does not support backend "{args.backend}"')
        init_datacards = self.prepare_init_datacards(tunesetup.datacards)
        use_full_covariance = mc_steer["data_covariance_mode"] == "full"

        # Compute the immutable physics and runtime initialization identity
        def fingerprint() -> str:
            return self.bootstrap_fingerprint(cdir=args.cdir, datacards=copy.deepcopy(tunesetup.datacards),
                mc_steer=copy.deepcopy(mc_steer), runtime_sha256=args.runtime_sha256,
                tunesetup_path=getattr(tunesetup, "__file__", None))

        # Compute the full covariance outputs stored in one reusable bootstrap
        def reusable_outputs() -> list:
            outputs = list(self.data_covariance_paths(run_name=args.run_name, cdir=args.cdir)) if use_full_covariance else []
            if getattr(args, "algorithm", None) == "ampfit":
                from core.tune.drivers.graniitti.ampfit.amplitude import bank_files

                outputs.extend(bank_files(self.amplitude_directory(cdir=args.cdir, run_name=args.run_name)))
            return outputs

        # Build or reuse the immutable GRANIITTI initialization archive
        def build_bootstrap() -> dict:
            return self.prepare_bootstrap(cdir=args.cdir, cache_base_url=args.bootstrap_cache_url,
                datacards=copy.deepcopy(tunesetup.datacards), mc_steer=copy.deepcopy(mc_steer),
                runtime_sha256=args.runtime_sha256, tunesetup_path=getattr(tunesetup, "__file__", None),
                init_force=args.init_force, initialize=lambda: self.initialize(
                    run_name=args.run_name, tunesetup=tunesetup, mc_steer=mc_steer,
                    obs_module=args.obs_module, cdir=args.cdir, init_force=args.init_force,
                    max_t=args.max_t, pickle_dump=args.pickle_dump, processes=args.cpu_per_trial, rngseed=args.rngseed,
                    bank_mode=({"phase": "combine", "jobs": int(args.bank_jobs),
                                "shared_root": getattr(args, "bank_shared_dir", None)}
                               if getattr(args, "bank_jobs", None) else (
                                   {"phase": "full", "shared_root": args.bank_shared_dir}
                                   if mc_steer.get("ampfit") and getattr(args, "bank_shared_dir", None) else None))),
                reusable_vgrids=lambda: self.reusable_vgrid_paths(
                    datacards=init_datacards, mc_steer=mc_steer, cdir=args.cdir), reusable_outputs=reusable_outputs)

        # Stage the archive and load its data covariance
        def stage_bootstrap(bootstrap: dict) -> None:
            self.stage_bootstrap(cdir=args.cdir, bootstrap=bootstrap)
            if use_full_covariance:
                self._load_bootstrap_covariance(
                    bootstrap=bootstrap, run_name=args.run_name, cdir=args.cdir, publish_current=True)
            else:
                self.data_covariance_payload = None
            self.prepare_trial_runtime({"data_covariance_mode": str(args.data_covariance_mode), "cost": str(args.cost),
                    "rngseed": int(args.rngseed)})

        return icetune.bootstrap_callbacks(args=args, tunesetup=tunesetup, mc_steer=mc_steer, simdriver=self,
            fingerprint_getter=fingerprint, bootstrap_builder=build_bootstrap, bootstrap_stager=stage_bootstrap)

    # Compute the shared location of amplitude banks in normal icetune initialization outputs
    def amplitude_directory(self, *, cdir: str, run_name: str) -> pathlib.Path:
        relative = self.amplitude_path or f"runs/icetune/{run_name}/results/amplitude"
        return pathlib.Path(cdir) / relative

    # Include initialized amplitude banks in the ordinary uploaded worker runtime
    def runtime_files(self, param: dict) -> list[str]:
        from core.tune.drivers.graniitti.ampfit.amplitude import bank_files

        files = super().runtime_files(param)
        if param.get("optimization", {}).get("optimizer") == "ampfit":
            root = pathlib.Path(param["cdir"])
            directory = self.amplitude_directory(cdir=param["cdir"], run_name=param["run_name"])
            if not directory.is_dir():
                raise ValueError("Amplitude initialization outputs are missing")
            files.extend(str(path.relative_to(root)) for path in bank_files(directory))
        return files

    # Compute the shared pickled HEPData cache path for one run
    def data_cache_path(self, *, run_name: str, cdir: str) -> str:
        return os.path.join(cdir, "runs", "icetune", run_name, "results", "data.pkl")

    # Compute the two persistent data covariance paths for one run
    def data_covariance_paths(self, *, run_name: str, cdir: str) -> tuple[str, str]:
        output_dir = os.path.join(cdir, "runs", "icetune", run_name, "results")
        return (os.path.join(output_dir, "data_covariance.npz"), os.path.join(output_dir, "data_covariance.json"))

    # Load one data covariance and its metadata
    def load_data_covariance(self, *, npz_path: str, json_path: str) -> None:
        self.data_covariance_payload = cov.load_data_covariance(npz_path=npz_path, json_path=json_path, data=self.data)

    # Save the currently loaded data covariance for one run
    def save_data_covariance(self, *, run_name: str, cdir: str) -> tuple[str, str]:
        if self.data_covariance_payload is None:
            raise RuntimeError("No data covariance is loaded")
        output_dir = os.path.join(cdir, "runs", "icetune", run_name, "results")
        return cov.save_data_covariance(self.data_covariance_payload, output_dir=output_dir)

    # Try to load previously initialized HEPData objects for this exact input setup
    def load_data_cache(self, *, run_name: str, datacards: list, obs_module: str, cdir: str) -> bool:
        target_path = self.data_cache_path(run_name=run_name, cdir=cdir)
        if not os.path.exists(target_path):
            return False
        input_files = self._data_input_file_records(cdir=pathlib.Path(cdir).resolve(), datacards=datacards)
        try:
            with open(target_path, "rb") as f:
                payload = pickle.load(f)
        except Exception as exc:
            logger.warning(".init_data: data cache load failed for %s (%s)", target_path, exc)
            return False

        if (not isinstance(payload, dict) or payload.get("replica_schema_version") != self.BOOTSTRAP_SCHEMA_VERSION
            or payload.get("datacards") != datacards or payload.get("obs_module") != obs_module
            or payload.get("input_files") != input_files):
            return False

        fields = ("datasets", "dataset_paths", "data", "obs", "pid", "cuts")
        if any(key not in payload for key in fields):
            return False
        for key in fields:
            setattr(self, key, payload[key])
        self.data_covariance_payload = payload.get("data_covariance_payload")
        self.initialized = True
        logger.info(".init_data: loaded data cache from %s", target_path)
        return True

    # Store initialized HEPData objects for later workers in the same run
    def write_data_cache(self, *, run_name: str, datacards: list, obs_module: str, cdir: str) -> None:
        store_dir = os.path.join(cdir, "runs", "icetune", run_name, "results")
        ensure_dir(store_dir)
        target_path = self.data_cache_path(run_name=run_name, cdir=cdir)
        payload = {"replica_schema_version": self.BOOTSTRAP_SCHEMA_VERSION, "datacards": copy.deepcopy(datacards),
            "input_files": self._data_input_file_records(
                cdir=pathlib.Path(cdir).resolve(), datacards=datacards, datasets=self.datasets),
            "datasets": self.datasets, "dataset_paths": self.dataset_paths, "data": self.data, "obs": self.obs,
            "obs_module": obs_module, "pid": self.pid, "cuts": self.cuts,
            "data_covariance_payload": self.data_covariance_payload, "unixtime": int(time.time())}

        tmp_path = icetune.default_tmp_path(target_path)
        with open(tmp_path, "wb") as f:
            logger.info(".init_data: saving data cache to %s", target_path)
            pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp_path, target_path)

    # Prepare shared simulator inputs and the requested amplitude bank phase
    def initialize(self, run_name: str, tunesetup: ModuleType, mc_steer: dict, obs_module: str, cdir: str,
        init_force: bool, max_t: int = 3600, pickle_dump: bool = True, processes: int = 1, rngseed: int = 0,
        bank_mode: dict | None = None,):
        """Run first init for VEGAS grids etc."""

        init_force = init_force or bool(mc_steer.get("ampfit") and not mc_steer.get("ampfit_reuse", True))
        logger.info(".initialize: tune_default=%s init_force=%s", mc_steer["tune_default"], init_force)

        init_datacards = prepare_init_datacards(tunesetup.datacards)

        self.init_data(
            run_name=run_name, datacards=tunesetup.datacards, obs_module=obs_module, cdir=cdir, pickle_dump=pickle_dump)

        # Ampfit integrates only the source proposal used to generate its fixed events
        if not mc_steer.get("ampfit"):
            self.compute(tunename=mc_steer["tune_default"], datacards=init_datacards, mc_steer=mc_steer, cdir=cdir,
                init_log_dir=os.path.join(cdir, "runs", "icetune", run_name, "init_logs"), max_t=max_t,
                init_force=init_force, processes=processes)
        bank_covariance = False
        if mc_steer.get("ampfit"):
            self.amplitude_path = None
            self.amplitude_banks.clear()
            # Distributed bank phases share the bank directory on the shared filesystem
            shared_root = bank_mode.get("shared_root") if bank_mode else None
            directory = self.amplitude_directory(cdir=shared_root or cdir, run_name=run_name)
            if directory.exists() and init_force and (bank_mode is None or bank_mode.get("phase") == "full"):
                directory.rename(directory.with_name(f"{directory.name}.{time.time_ns()}._old"))
            initial = self.get_initial_param(tunesetup.param_space, tunesetup.aux_param_space, cdir, mc_steer["tune_default"])
            initial.update(tunesetup.aux_param_space)
            from core.tune.parameters.space import normalize_continuous_param_space

            controls = mc_steer["ampfit"]
            bounds = normalize_continuous_param_space(tunesetup.param_space)
            temporary = pathlib.Path(cdir) / "runs" / "icetune" / run_name / "init" / f"ampfit.{time.time_ns()}"
            tune = self.trial_tunename(trial_id="source", pid=os.getpid(), purpose="ampfit")
            self.create_steering_card(param_space=initial, tunename=tune, cdir=cdir, tune_default=mc_steer["tune_default"])
            bank_covariance = (mc_steer["data_covariance_mode"] == "full"
                and all(card["weighted"] for card in tunesetup.datacards) and prepare_mc_correlation_datacards(
                    tunesetup.datacards, event_count=mc_steer.get("mc_correlation_events"),
                    weighting=mc_steer.get("mc_correlation_weighting", "card")) == tunesetup.datacards)
            amplitude = {"directory": directory, "prepare": True, "parameters": initial,
                         "bounds": bounds, "temporary": temporary, "controls": controls,
                         "reuse": not init_force and mc_steer.get("ampfit_reuse", True), "force": init_force,
                         "reuse_root": (pathlib.Path(mc_steer["ampfit_reuse_root"]) / "runs/icetune"
                                        if mc_steer.get("ampfit_reuse_root") else directory.parents[2])}
            if bank_mode is not None:
                amplitude.update(sample=bank_mode["phase"] == "sample", phase=bank_mode["phase"], covariance=bank_covariance)
                amplitude.update({key: bank_mode[key] for key in ("index", "shard", "jobs") if key in bank_mode})
            self.compute(tunename=tune, datacards=tunesetup.datacards, mc_steer=mc_steer,
                         cdir=cdir, max_t=max_t, processes=processes, rngseed=rngseed,
                         collect_data_covariance=bank_covariance and (bank_mode is None or bank_mode["phase"] in {"combine", "full"}),
                         data_covariance_run_name=run_name, amplitude=amplitude)
            # Mirror the combined shared bank into the local initialization outputs
            if bank_mode is not None and bank_mode.get("phase") in {"combine", "full"} and shared_root:
                local_directory = self.amplitude_directory(cdir=cdir, run_name=run_name)
                if local_directory.resolve() != directory.resolve():
                    if local_directory.exists():
                        local_directory.rename(local_directory.with_name(f"{local_directory.name}.{time.time_ns()}._old"))
                    from core.tune.drivers.graniitti.ampfit.amplitude import stage_bank

                    for bank in directory.glob("*/*/amplitudes.bin.json"):
                        relative = bank.parent.relative_to(directory)
                        if not any(part.endswith(("._old", ".partial")) for part in relative.parts):
                            stage_bank(bank.parent, local_directory / relative, predictions=False)
        # Per-sample jobs stop before the combined covariance and cache publication
        if bank_mode is not None and bank_mode.get("phase") in {"sample", "shard", "finalize"}:
            return None
        if mc_steer["data_covariance_mode"] == "full":
            if not bank_covariance:
                self.initialize_data_covariance(run_name=run_name, tunename=mc_steer["tune_default"],
                    datacards=tunesetup.datacards, mc_steer=mc_steer, cdir=cdir, max_t=max_t)
        else:
            self.data_covariance_payload = None
        if pickle_dump:
            self.write_data_cache(run_name=run_name, datacards=tunesetup.datacards, obs_module=obs_module, cdir=cdir)

        return self.get_initial_param(param_space=tunesetup.param_space, aux_param_space=tunesetup.aux_param_space,
            cdir=cdir, tune_default=mc_steer["tune_default"])

    # Load the verified shared INIT state without repeating grid generation
    def initialize_ray_bootstrap(self, *, args, tunesetup: ModuleType, mc_steer: dict) -> dict | None:
        self.init_data(run_name=args.run_name, datacards=tunesetup.datacards, obs_module=args.obs_module,
            cdir=args.cdir, pickle_dump=args.pickle_dump)
        from core.tune import cache as icetune_cache

        inputs = self.bootstrap_inputs(cdir=args.cdir, datacards=copy.deepcopy(tunesetup.datacards),
            mc_steer=copy.deepcopy(mc_steer), runtime_sha256=args.runtime_sha256,
            tunesetup_path=getattr(tunesetup, "__file__", None))
        bootstrap, points = icetune.load_ray_init_state(
            path=args.ray_init_state, fingerprint=icetune_cache.json_fingerprint(inputs), fingerprint_inputs=inputs)
        self._bind_amplitude(bootstrap)
        if mc_steer.get("ampfit") and self.amplitude_path is None:
            raise ValueError("The initialization manifest contains no amplitude banks")
        if mc_steer["data_covariance_mode"] == "full":
            self._load_bootstrap_covariance(bootstrap=bootstrap, run_name=args.run_name, cdir=args.cdir, publish_current=True)
        else:
            self.data_covariance_payload = None
        return points

    # Compute the physics and executable identity used to fence Ray restarts
    def physics_fingerprint(self, *, param: dict, tunesetup) -> str:
        return self.bootstrap_fingerprint(cdir=param["cdir"], datacards=param["datacards"], mc_steer=param["mc_steer"],
            runtime_sha256=param.get("runtime_sha256"), tunesetup_path=getattr(tunesetup, "__file__", None))

    # Generate the MC correlation sample and fix the data covariance
    def initialize_data_covariance(
        self, *, run_name: str, tunename: str, datacards: list, mc_steer: dict, cdir: str, max_t: int
    ) -> tuple[str, str]:
        mc_correlation_datacards = prepare_mc_correlation_datacards(datacards,
            event_count=mc_steer.get("mc_correlation_events"),
            weighting=mc_steer.get("mc_correlation_weighting", "card"))
        self.compute(tunename=tunename, datacards=mc_correlation_datacards, mc_steer=mc_steer, cdir=cdir, processes=1,
            max_t=max_t, collect_data_covariance=True, data_covariance_run_name=run_name)
        return self.data_covariance_paths(run_name=run_name, cdir=cdir)

    # Initialize HEPData comparisons with independent storage for each dataset and set
    def init_data(self, run_name, datacards: list, obs_module: str, cdir: str, pickle_dump: bool = True):
        if cdir is None:
            cdir = os.getcwd()
            logger.info(".init_data: cdir=%s", cdir)

        if not pickle_dump and self.load_data_cache(
            run_name=run_name, datacards=datacards, obs_module=obs_module, cdir=cdir):
            return

        self.datasets = [None] * len(datacards)
        self.dataset_paths = [None] * len(datacards)
        for i, card in enumerate(datacards):
            dataset, self.dataset_paths[i] = dataset_steering.load_dataset(card["datacard"], cdir=cdir)
            self.datasets[i] = _apply_force_density(_apply_elastic_model(dataset, card), card)

        if len(self.datasets) == 0:
            err_msg = __name__ + ".init_data: Error: did not find any input datasets -- check your input"
            logger.error(err_msg)
            raise Exception(err_msg)

        for name in ("data", "obs", "pid", "cuts"):
            setattr(self, name, [[None] * len(dataset["sets"]) if dataset["active"] else [] for dataset in self.datasets])
        for i, dataset in enumerate(self.datasets):
            if not dataset["active"]:
                continue
            density = dataset["fit"]["normalization"] == "unit_density"
            for j, dataset_set in enumerate(dataset["sets"]):
                all_obs, _ = dataset_steering.load_observables(
                    dataset_set["obs"], dataset_path=self.dataset_paths[i], cdir=cdir)
                if dataset_set.get("data", True):
                    hepdata, self.obs[i][j] = readers.read_hepdata(dataset=dataset_set,
                        datapath=dataset["datapath"], datatype=dataset["type"], all_obs=all_obs,
                        cdir=cdir, reader=dataset.get("reader"), dataset_path=self.dataset_paths[i])

                    self.data[i][j] = plot.histhepdata(hepdata=hepdata, obs=self.obs[i][j], density=density,
                        density_uncertainty="shape", label=dataset_set["name"])
                else:
                    self.obs[i][j] = dataset_steering.select_histogram_observables(dataset_set, all_obs)
                self.pid[i][j] = dataset_set["pid"]
                self.cuts[i][j] = dataset_steering.resolve_python_reference(dataset_set["cuts"],
                    package="core.analysis.cuts", dataset_path=self.dataset_paths[i], cdir=cdir)

        self.initialized = True

        if pickle_dump:
            self.write_data_cache(run_name=run_name, datacards=datacards, obs_module=obs_module, cdir=cdir)

    # Resolve one dataset and datacard into stable GRANIITTI run paths
    def _run_card(self, *, index: int, sample_index: int = 0, datacard: dict, mc_steer: dict, tunename: str, cdir: str
    ) -> _RunCard:
        sample = self.datasets[index]["samples"][sample_index]
        parameters = sample["parameters"]
        overrides = tuple(
            f"{path}={dataset_steering.serialize_override_value(value, context=f'Generator override {path!r}')}"
            for path, value in parameters.items())
        inputfile = dataset_steering.resolve_gencard_reference(
            sample["gencard"], dataset_path=self.dataset_paths[index], cdir=cdir)
        steering = load_json_file(inputfile, loader=json5.load)
        process = str(_sample_steering_value(steering, parameters, "SCATTERING.PROCESS"))
        loopscreen = "true" if datacard["loopscreen"] else "false"
        xsmode = str(datacard["xsmode"])
        swap_process = datacard.get("swap_process")
        integrator = datacard.get("integrator")
        if integrator is not None:
            integrator = str(integrator).upper()
            if integrator not in {"VEGAS", "NEUROJAC"}:
                raise ValueError("GRANIITTI ICETUNE integrator must be VEGAS or NEUROJAC")
        steering_integrator = str(_sample_steering_value(steering, parameters, "GENERIC.INTEGRATOR", default="VEGAS")
        ).upper()
        resolved_integrator = integrator or steering_integrator
        if resolved_integrator not in {"VEGAS", "NEUROJAC"}:
            raise ValueError("GRANIITTI gencard integrator must be VEGAS or NEUROJAC")
        card_identity, beam_energy = _vgrid_card_identity(steering, parameters)
        steering_output = str(_sample_steering_value(steering, parameters, "GENERIC.OUTPUT"))
        tunesetup_name = pathlib.Path(mc_steer["tunesetup_name"]).stem
        output = f"{steering_output}__{card_identity}__{tunesetup_name}__LOOPSCREEN_{loopscreen}"
        if swap_process is not None:
            output += f"__{swap_process}"
        if integrator is not None:
            output += f"__INTEGRATOR_{integrator}"
        output += "__" + _tune_tag(cdir=cdir, tunename=tunename if xsmode == "reset" else mc_steer["tune_default"])
        return _RunCard(gridfile=os.path.join(cdir, "vgrid", f"{output}.vgrid"), inputfile=inputfile,
            integrator=resolved_integrator, integrator_override=integrator, output=output, loopscreen=loopscreen,
            process=process,
            process_override=(_swap_process_tag(process, swap_process) if swap_process is not None else None),
            sample_name=sample["name"], sample_scale=float(sample.get("scale", 1.0)), overrides=overrides,
            beam_energy=beam_energy, xsmode=xsmode)

    # Resolve every generator sample belonging to one tuning dataset
    def _run_cards(self, *, index: int, datacard: dict, mc_steer: dict, tunename: str, cdir: str) -> list[_RunCard]:
        dataset = self.datasets[index]
        samples = dataset["samples"]
        if len(samples) > 1 and not all(dataset_steering.set_sample_names(entry) for entry in dataset.get("sets", [])):
            raise ValueError(f"GRANIITTI icetune datacard[{index}] with multiple samples requires routed sets")
        if any(len(dataset_steering.set_sample_names(entry)) > 1 for entry in dataset.get("sets", [])):
            raise ValueError(f"GRANIITTI icetune datacard[{index}] requires one MC sample per fitted set")
        return [self._run_card(index=index, sample_index=sample_index, datacard=datacard, mc_steer=mc_steer,
                tunename=tunename, cdir=cdir) for sample_index in range(len(samples))]

    # Compute the dataset set indices routed to one generator sample
    def _sample_set_indices(self, *, index: int, sample_name: str) -> list[int]:
        dataset_sets = self.datasets[index].get("sets", [])
        routed = any(dataset_steering.set_sample_names(entry) for entry in dataset_sets)
        if not routed:
            return list(range(len(dataset_sets)))
        return [
            set_index for set_index, entry in enumerate(dataset_sets) if sample_name in dataset_steering.set_sample_names(entry)
        ]

    # Execute one sample command with its steering overrides and failure diagnostics
    def _sample_command(self, cmd: list[str], *, run: _RunCard, cdir: str, deadline: float, stage: str,
        dataset_index: int, output: str, **logging,) -> None:
        for option, value in (("-p", run.process_override), ("-g", run.integrator_override)):
            if value is not None:
                cmd.extend((option, value))
        for override in run.overrides:
            cmd.extend(("--set", override))
        result = {}
        status = iceruntime.execute_cmd(cmd=cmd, cwd=cdir, max_t=_remaining_time(deadline), result_out=result, **logging
        )
        if status is False:
            error = _graniitti_command_error(result=result, stage=stage, datacard_index=dataset_index,
                inputfile=run.inputfile, gridfile=run.gridfile, process=run.process_override or run.process,
                output=output)
            logger.error("%s", error)
            raise error

    # Compute initialization-produced grids which are valid for every sampled trial
    def reusable_vgrid_paths(self, *, datacards: list, mc_steer: dict, cdir: str) -> list[str]:
        if not self.initialized:
            raise RuntimeError("GRANIITTI data must be initialized before resolving reusable vgrids")
        paths: list[str] = []
        if mc_steer.get("ampfit"):
            return paths
        for index, dataset in enumerate(self.datasets):
            if not dataset["active"] or str(datacards[index]["xsmode"]) == "reset":
                continue
            runs = self._run_cards(
                index=index, datacard=datacards[index], mc_steer=mc_steer, tunename=mc_steer["tune_default"], cdir=cdir)
            paths.extend(run.gridfile for run in runs)
        return paths

    # Generate and histogram one routed dataset sample
    def _compute_sample(self, *, dataset_index: int, set_indices: list[int], run: _RunCard, tunename: str,
        datacard: dict, cdir: str, init_log_dir: str | None, chunksize: int, processes: int, deadline: float,
        rng: np.random.Generator, init_force: bool, collect_data_covariance: bool, covariance_mode: str = "diagonal",
        amplitude: dict | None = None, event_dir: pathlib.Path | None = None,) -> tuple[list | None, list | None]:
        dataset = self.datasets[dataset_index]
        nevents = datacard["nevents"]
        weighted = "true" if datacard["weighted"] else "false"
        density = dataset["fit"]["normalization"] == "unit_density"
        local_obs = [self.obs[dataset_index][set_index] for set_index in set_indices]
        local_pid = [self.pid[dataset_index][set_index] for set_index in set_indices]
        local_cuts = [self.cuts[dataset_index][set_index] for set_index in set_indices]
        set_mc_scales = [
            run.sample_scale * float(dataset["sets"][set_index].get("mc_scale", 1.0)) for set_index in set_indices]
        saved_events = None
        shard_only = amplitude is not None and not amplitude.get("sample") and (
            amplitude.get("shard") is not None or amplitude.get("jobs") is not None)
        if amplitude is not None:
            from core.tune.drivers.graniitti.ampfit.amplitude import AmplitudeBank

            relative = pathlib.Path(str(dataset_index)) / readers.clean_filename(run.sample_name)
            bank_directory = pathlib.Path(amplitude["directory"]) / relative
            bank_options = dict(driver=self, directory=bank_directory, obs=local_obs, pid=local_pid,
                                cuts=local_cuts, density=density, scales=set_mc_scales,
                                controls=amplitude["controls"], covariance_mode=covariance_mode,
                                workers=amplitude.get("workers", {}).get(relative.as_posix()))
        if amplitude is not None and amplitude["prepare"] and nevents > 0:
            if amplitude.get("phase") == "combine":
                from core.tune.drivers.graniitti.ampfit.amplitude import load_finalized

                return load_finalized(bank_directory)
            from core.tune.drivers.graniitti.ampfit.amplitude import BankPlan

            plan = BankPlan(driver=self, run=run, tune=_resolve_tune_dir(cdir=cdir, tunename=tunename),
                         initial=amplitude["parameters"], bounds=amplitude["bounds"], pid=local_pid[0],
                         controls=amplitude["controls"], cdir=cdir, nevents=nevents)
            if amplitude.get("reuse", True) and not amplitude.get("force"):
                plan.recover(bank_directory)
            ready = not (amplitude.get("sample") and amplitude.get("force")) and not plan.mismatch(bank_directory)
            finalize = amplitude.get("phase") == "finalize" or (
                amplitude.get("phase") == "shard" and amplitude.get("jobs") == 1)
            if finalize and (ready or amplitude.get("phase") == "finalize"):
                from core.tune.drivers.graniitti.ampfit.amplitude import finalize_bank

                finalize_bank(directory=bank_directory, temporary=pathlib.Path(amplitude["temporary"]) / relative / "finalize",
                    options=bank_options, parameters=amplitude["parameters"],
                    covariance=amplitude["covariance"], ready=ready)
                return None, None
            if not ready and not shard_only and amplitude.get("reuse", True):
                source = plan.reuse(bank_directory, amplitude.get("reuse_root", pathlib.Path(amplitude["directory"]).parents[2]))
                ready = source is not None
                if ready:
                    logger.info(".compute: copied compatible amplitudes from previous run %s", source)
            if ready:
                if amplitude.get("sample") or amplitude.get("shard") is not None:
                    return None, None
                bank = AmplitudeBank(**bank_options)
                bank.validate_source(bank.amplitude(bank.initial), amplitude["controls"]["closure_rtol"], bank_directory / "closure.json")
                logger.info(".compute: reusing prepared amplitudes from %s", bank_directory)
                return bank.predict(amplitude["parameters"]), bank.sample if collect_data_covariance else None
            if bank_directory.exists() and not shard_only:
                reuse_events = not amplitude.get("force") and plan.reusable_events(bank_directory)
                previous = bank_directory.with_name(f"{bank_directory.name}.{time.time_ns()}._old")
                bank_directory.rename(previous)
                if reuse_events:
                    saved_events = previous / "events.hepmc3"
                    logger.info(".compute: reusing completed event sample from %s", saved_events)
        if amplitude is not None and not amplitude["prepare"]:
            key = str(bank_directory)
            bank = self.amplitude_banks.get(key)
            if bank is None:
                bank = AmplitudeBank(**bank_options)
            self.amplitude_banks[key] = bank
            return bank.predict(amplitude["parameters"]), None
        init_log_path = None
        if init_log_dir is not None:
            ensure_dir(init_log_dir)
            init_log_path = os.path.join(init_log_dir, f"{readers.clean_filename(run.output)}.json")

        # Basis shards and the combine step reuse the shared event sample
        command_stage = None
        rebuild_reason = None
        init_loopscreen = "false" if run.xsmode == "sample" else run.loopscreen
        modelparam = resolve_gr_modelparam(cdir=cdir, tunename=tunename)
        command_log_stage = None
        if not shard_only and (init_force or not os.path.exists(run.gridfile)):
            command_stage = "vgrid_initialization"
            command_log_stage = "vgrid_init"
        elif not shard_only:
            try:
                griddata = load_json_file(run.gridfile, loader=json5.load)
                rebuild_reason = _vgrid_rebuild_reason(griddata, run)
            except Exception as exc:
                rebuild_reason = f"grid cannot be read: {exc}"
            if rebuild_reason is not None:
                logger.warning(".compute: incompatible vgrid file, re-computing (%s)", rebuild_reason)
                command_stage = "vgrid_reinitialization"
                command_log_stage = "vgrid_reinit_incompatible"

        if command_stage is not None:
            logger.info(".compute: running vgrid initialization to %s", run.gridfile)
            cmd = _gr_command(cdir, ("h", 0), ("c", processes), ("w", weighted), ("n", -1), ("i", run.inputfile),
                ("o", run.output), ("m", modelparam), ("l", init_loopscreen))
            self._sample_command(cmd, run=run, cdir=cdir, deadline=deadline, stage=command_stage,
                dataset_index=dataset_index, output=run.output, log_path=init_log_path, log_metadata={
                    "datacard_index": dataset_index, "gridfile": run.gridfile, "inputfile": run.inputfile,
                    "integrator": run.integrator_override, "output": run.output, "sample": run.sample_name,
                    "stage": command_log_stage, "tunename": tunename, "vgrid_rebuild_reason": rebuild_reason,
                    "xsmode": run.xsmode})

        if nevents <= 0:
            return None, None

        hepmc3output = f"{run.output}_{_tune_tag(cdir=cdir, tunename=tunename)}"
        rngseed = int(rng.integers(1, int(1e6)))
        cmd = _gr_command(cdir, ("h", 0), ("c", processes), ("w", weighted), ("l", run.loopscreen), ("n", nevents),
            ("m", modelparam), ("d", run.gridfile), ("i", run.inputfile), ("o", hepmc3output), ("r", rngseed))
        if saved_events is None and not shard_only:
            self._sample_command(cmd, run=run, cdir=cdir, deadline=deadline, stage="event_generation",
                dataset_index=dataset_index, output=hepmc3output)

        outputfile = os.path.join(cdir, f"output/{hepmc3output}")
        hepmc3file = outputfile + ".hepmc3" if saved_events is None else str(saved_events)
        read_mode = "header" if run.xsmode == "reset" else run.xsmode
        mcdata = None

        if amplitude is not None:
            from core.tune.drivers.graniitti.ampfit.amplitude import prepare_bank, prepare_bank_shard, prepare_sample

            bank_shard = amplitude.get("shard")
            bank_jobs = amplitude.get("jobs")
            preparation = dict(driver=self, directory=bank_directory, run=run, controls=amplitude["controls"],
                               cdir=cdir, deadline=deadline)
            source = dict(tune=_resolve_tune_dir(cdir=cdir, tunename=tunename), initial=amplitude["parameters"],
                          bounds=amplitude["bounds"], pid=local_pid[0], events=hepmc3file)
            if amplitude.get("sample"):
                # Publish the shared event sample for the parallel basis shards
                prepare_sample(**preparation, **source, shards=bank_jobs or 1)
                return None, None
            if bank_shard is not None:
                # Evaluate only this basis shard on the shared event sample
                prepare_bank_shard(**preparation, temporary=pathlib.Path(amplitude["temporary"]) / relative,
                                   shard=bank_shard, shards=bank_jobs, processes=processes)
                if bank_jobs == 1:
                    from core.tune.drivers.graniitti.ampfit.amplitude import finalize_bank

                    finalize_bank(directory=bank_directory, temporary=pathlib.Path(amplitude["temporary"]) / relative / "finalize",
                        options=bank_options, parameters=amplitude["parameters"],
                        covariance=amplitude["covariance"], ready=False)
                return None, None
            prepare_bank(**preparation, **source, temporary=pathlib.Path(amplitude["temporary"]) / relative,
                         processes=processes)
            bank = AmplitudeBank(**bank_options)
            bank.validate_source(bank.amplitude(bank.initial), amplitude["controls"]["closure_rtol"], bank_directory / "closure.json")
            mc = bank.predict(amplitude["parameters"])
            mcdata = bank.sample
        elif processes > 1:
            weight_summary = readers.read_hepmc3_weight_summary(filename=hepmc3file)
            event_count = weight_summary["nevents"]
            if event_count is None or event_count == 0:
                err_msg = __name__ + f".compute: Error, no events found in file: {hepmc3file}"
                logger.error(err_msg)
                raise Exception(err_msg)
            read_mode, _ = readers.resolve_hepmc3_xsmode(requested_xsmode=read_mode, weight_summary=weight_summary)
            same_param = {"chunk_range": None, "hepmc3file": hepmc3file, "obs": local_obs, "pid": local_pid,
                "cuts": local_cuts, "xsmode": read_mode,
                "header_xsection": (readers.terminal_header_xsection(weight_summary) if read_mode == "header" else None),
                "scales": set_mc_scales, "k": 0, "label": f"GRANIITTI {graniitti_runtime.version(cdir)}",
                "density": density, "density_uncertainty": "shape" if density else "scaled",
                "covariance_mode": covariance_mode, "verbose": False}
            paramlist = [{**copy.deepcopy(same_param), "chunk_range": chunk_range}
                for chunk_range in readers.generate_partitions(totalsize=event_count, chunksize=chunksize)]
            with multiprocessing.Pool(processes=processes) as pool:
                hists = pool.map_async(readers.parallel_wrapper, paramlist).get(timeout=_remaining_time(deadline))
            mc = plot.fuse_worker_chunk_outputs(hists)
        else:
            mcdata = readers.read_hepmc3(
                hepmc3file=hepmc3file, obs=local_obs, pid=local_pid, cuts=local_cuts, xsmode=read_mode)
            mc = [plot.histmc(mcdata=mcdata[local_index], obs=local_obs[local_index], density=density,
                    density_uncertainty="shape" if density else "scaled", covariance_mode=covariance_mode,
                    scale=set_mc_scales[local_index], color=plot.colors(0),
                    label=f"GRANIITTI {graniitti_runtime.version(cdir)}") for local_index in range(len(set_indices))]

        if event_dir is not None:
            ensure_dir(event_dir, exist_ok=False)
            shutil.copy2(run.inputfile, event_dir / "gencard.json")
            shutil.copy2(run.gridfile, event_dir / "vgrid.json")
            (event_dir / "sample.json").write_text(json.dumps(json_safe({
                "dataset_index": dataset_index, "dataset": dataset, "datacard": datacard,
                "run": asdict(run), "command": cmd, "rngseed": rngseed,
                "weighted": datacard["weighted"], "nevents": nevents}), sort_keys=True, indent=2))
            shutil.move(hepmc3file, event_dir / "events.hepmc3")
        elif saved_events is None:
            pathlib.Path(hepmc3file).unlink(missing_ok=True)
        if run.xsmode == "reset":
            pathlib.Path(run.gridfile).unlink(missing_ok=True)
        return mc, mcdata if collect_data_covariance else None

    # Prepare independent amplitude samples within the allocated native CPU budget
    def _compute_samples(self, jobs: list[dict], processes: int):
        if jobs and jobs[0]["amplitude"] is not None and jobs[0]["amplitude"].get("shard") is not None:
            # Give each bank the full subprocess budget, largest event sample first
            for job in sorted(jobs, key=lambda job: job["datacard"]["nevents"], reverse=True):
                yield job, self._compute_sample(**{**job, "processes": processes})
            return
        parallel = jobs and jobs[0]["amplitude"] is not None and jobs[0]["amplitude"]["prepare"]
        if not parallel or processes <= 1:
            for job in jobs:
                yield job, self._compute_sample(**job)
            return
        cpus, remainder = divmod(processes, min(processes, len(jobs)))
        jobs = [{**job, "processes": cpus + int(index < remainder)} for index, job in enumerate(jobs)]
        with ThreadPoolExecutor(max_workers=min(processes, len(jobs))) as pool:
            futures = [pool.submit(self._compute_sample, **job) for job in jobs]
            try:
                for job, future in zip(jobs, futures, strict=True):
                    yield job, future.result()
            finally:
                for future in futures:
                    future.cancel()

    # Compute routed MC samples and their shared histogram predictions
    def compute(self, tunename: str, datacards: list, mc_steer: dict, cdir: str | None = None,
        init_log_dir: str | None = None, chunksize: int = 10000, processes: int = 8, max_t: int = 3600,
        deadline: float | None = None, rngseed: int | None = None, init_force: bool = False,
        collect_data_covariance: bool = False, data_covariance_run_name: str | None = None,
        amplitude: dict | None = None, event_dir: pathlib.Path | None = None,) -> list:
        """Compute MC sample and compare with HEPData over all datacards input."""

        if not self.initialized:
            raise Exception(__name__ + ".compute: Error -- data is not initialized!")

        if cdir is None:
            cdir = os.getcwd()
            logger.info(".compute: cdir=%s", cdir)

        logger.info(
            '.compute: executing tunename="%s" multiprocessing.cpu_count=%s', tunename, multiprocessing.cpu_count())

        if deadline is None:
            deadline = time.monotonic() + float(max_t)
        rng = np.random.default_rng(rngseed)

        mc_all = [[] for _ in self.datasets]
        mc_correlation_data_by_dataset = [None] * len(self.datasets)
        jobs = []
        for i, dataset in enumerate(self.datasets):
            if not dataset["active"]:
                continue
            mc_all[i] = [None] * len(dataset.get("sets", []))
            mc_correlation_data_by_dataset[i] = {"sample_groups": []}
            runs = self._run_cards(index=i, datacard=datacards[i], mc_steer=mc_steer, tunename=tunename, cdir=cdir)
            for run in runs:
                set_indices = self._sample_set_indices(index=i, sample_name=run.sample_name)
                jobs.append(dict(event_dir=None if event_dir is None else event_dir / f"sample_{len(jobs):04d}",
                    dataset_index=i, set_indices=set_indices, run=run, tunename=tunename, datacard=datacards[i],
                    cdir=cdir, init_log_dir=init_log_dir, chunksize=chunksize, processes=processes, deadline=deadline,
                    rng=np.random.default_rng(rng.integers(0, 2**63)) if amplitude is not None else rng,
                    init_force=init_force, collect_data_covariance=collect_data_covariance,
                    covariance_mode=mc_steer.get("data_covariance_mode", "diagonal"), amplitude=amplitude))
        # Select one bank after assigning independent sample seeds
        bank_index = amplitude.get("index") if amplitude is not None else None
        if bank_index is not None:
            if not 0 <= bank_index < len(jobs):
                raise ValueError("Amplitude bank index is outside the active sample list")
            jobs = [jobs[bank_index]]
        for job, (sample_mc, sample_mcdata) in self._compute_samples(jobs, processes):
            i, set_indices = job["dataset_index"], job["set_indices"]
            if sample_mc is not None:
                for local_index, set_index in enumerate(set_indices):
                    mc_all[i][set_index] = sample_mc[local_index]
            if sample_mcdata is not None:
                mc_correlation_data_by_dataset[i]["sample_groups"].append(
                    {"mcdata": sample_mcdata, "sample": job["run"].sample_name, "set_indices": set_indices})
        bank_only = amplitude is not None and amplitude.get("phase") in {"sample", "shard", "finalize"}
        for i, dataset in enumerate(self.datasets):
            if not dataset["active"]:
                continue
            if datacards[i]["nevents"] > 0 and not bank_only:
                if any(value is None for value in mc_all[i]):
                    raise RuntimeError(f"GRANIITTI icetune dataset[{i}] has unrouted MC sets")
            else:
                mc_all[i] = None

        _cleanup_temporary_tune_dir(cdir=cdir, tunename=tunename)

        if collect_data_covariance:
            self.data_covariance_payload = cov.build_data_covariance(data=self.data, obs=self.obs,
                mcdata_by_dataset=mc_correlation_data_by_dataset, mc_correlation_datacards=datacards)
            cov.validate_data_covariance(self.data_covariance_payload, data=self.data)
            run_name = data_covariance_run_name
            if run_name is None:
                raise ValueError("collect_data_covariance requires data_covariance_run_name")
            self.save_data_covariance(run_name=run_name, cdir=cdir)

        return mc_all

    def cleanup_trial_outputs(self, tunename: str, datacards: list, mc_steer: dict, cdir: str | None = None) -> None:
        """Best-effort cleanup for per-trial outputs."""

        if cdir is None:
            cdir = os.getcwd()

        for i in range(len(self.datasets)):
            if not self.datasets[i]["active"]:
                continue

            try:
                runs = self._run_cards(index=i, datacard=datacards[i], mc_steer=mc_steer, tunename=tunename, cdir=cdir)
            except Exception:
                continue

            for run in runs:
                if datacards[i]["nevents"] > 0:
                    hepmc3output = f"{run.output}_{_tune_tag(cdir=cdir, tunename=tunename)}"
                    hepmc3file = os.path.join(cdir, f"output/{hepmc3output}.hepmc3")
                    pathlib.Path(hepmc3file).unlink(missing_ok=True)

                if run.xsmode == "reset":
                    pathlib.Path(run.gridfile).unlink(missing_ok=True)

        _cleanup_temporary_tune_dir(cdir=cdir, tunename=tunename)

    def get_initial_param(self, param_space: dict, aux_param_space: dict, cdir: str, tune_default: str = "TUNE0"
    ) -> dict:
        """Read default values for tunable parameters from the JSON cards."""

        p = dict.fromkeys(param_space.keys(), 0.0)
        for key in param_space:
            if key.endswith("FF_transfer.LambdaInv2"):
                p[key] = 1.0
        for key in aux_param_space:
            p[key] = aux_param_space[key]

        tunename = self.trial_tunename(trial_id="probe", pid=os.getpid(), purpose="probe")
        orig_param_values = self.create_steering_card(
            tunename=tunename, param_space=p, cdir=cdir, tune_default=tune_default)
        _cleanup_temporary_tune_dir(cdir=cdir, tunename=tunename)

        for key in aux_param_space:
            del orig_param_values[key]

        return orig_param_values

    def create_steering_card(
        self, param_space: dict, tunename: str, cdir: str | None = None, tune_default: str = "TUNE0") -> dict:
        """
        Create a new set of steering cards by copy and update parameter json files.
        """

        if cdir is None:
            cdir = os.getcwd()
            logger.info(".create_steering_card: cdir=%s", cdir)

        if tunename != tune_default:
            src = _resolve_tune_dir(cdir=cdir, tunename=tune_default)
            dst = _resolve_tune_dir(cdir=cdir, tunename=tunename)
            ensure_dir(dst.parent)
            _cleanup_temporary_tune_dir(cdir=cdir, tunename=tunename)
            stage = icetune.default_tmp_stage_path(dst)

            try:
                ensure_dir(stage, exist_ok=False)
                for item in src.iterdir():
                    if item.is_file():
                        shutil.copyfile(item, stage / item.name)
                    elif item.is_dir():
                        shutil.copytree(item, stage / item.name)
                os.rename(stage, dst)
            except Exception:
                shutil.rmtree(stage, ignore_errors=True)
                raise
        else:
            logger.info(".create_steering_card: tunename == tune_default (%s), not creating a new one", tune_default)

        path = str(_resolve_tune_dir(cdir=cdir, tunename=tunename))
        prepared_param, physical_keys, groups = self._prepare_card_params(param_space, path=path)
        eikonal_groups, phase_groups, complex_groups, projective_groups = groups
        orig_gen_values = self.update_general_param(path=path, param=prepared_param)
        orig_res_values = self.update_resonance_param(path=path, param=prepared_param)
        orig_br_values = self.update_branching_param(path=path, param=prepared_param)
        orig_prod_values = self.update_production_param(path=path, param=prepared_param)

        orig_values = {**orig_gen_values, **orig_res_values, **orig_br_values, **orig_prod_values}
        orig_values = self._direction_defaults(orig_values, projective_groups)
        orig_values = _encode_coordinate_defaults(orig_values, phase_groups)
        orig_values = _encode_coordinate_defaults(orig_values, complex_groups)
        orig_values = _encode_coordinate_defaults(orig_values, eikonal_groups)
        return self._restore_physical_optimizer_keys(orig_values, physical_keys)

    # Validate and return the selected GRANIITTI push card families
    def _push_card_selectors(self, options: dict) -> set[str]:
        allowed = {"GENERAL.json", "RES*.json", "DECAYS.json", "CON_MP.json", "CON_XP.json", "CON_GP.json",
            "CON_TP.json"}
        unknown_options = set(options) - {"cards"}
        if unknown_options:
            raise ValueError(f"Unknown GRANIITTI push options: {sorted(unknown_options)}")
        cards = options.get("cards", sorted(allowed))
        if not isinstance(cards, list) or not cards or any(not isinstance(card, str) for card in cards):
            raise ValueError('GRANIITTI push option "cards" must be a non-empty string array')
        unknown_cards = set(cards) - allowed
        if unknown_cards:
            raise ValueError(f"Unknown GRANIITTI push cards: {sorted(unknown_cards)}")
        return set(cards)

    # Compute the GRANIITTI card family containing one physical parameter
    def _push_parameter_card(self, parameter: str) -> str:
        parameter_class = str(parameter).partition("|")[0]
        if parameter_class == "RES":
            return "RES*.json"
        if parameter_class == "DECAY":
            return "DECAYS.json"
        if parameter_class in {"CON_MP", "CON_XP", "CON_GP", "CON_TP"}:
            return f"{parameter_class}.json"
        return "GENERAL.json"

    # Publish optimized parameters only into the selected GRANIITTI cards
    def push_parameters(self, *, summary: dict, target_path: str, cdir: str, options: dict,
        confirm: Callable[[list[tuple[str, object, object]]], bool],) -> list[tuple[str, object, object]]:
        config, target = self._resolve_push_inputs(summary, target_path, cdir)
        if not target.is_dir():
            raise ValueError(f'GRANIITTI push target "{target}" is not a tune directory')
        cards = self._push_card_selectors(options)
        selected_config = {key: value for key, value in config.items() if self._push_parameter_card(key) in cards}
        if not selected_config:
            raise ValueError(f"GRANIITTI push selected no parameters for cards {sorted(cards)}")

        tunename = self.trial_tunename(trial_id="probe", pid=os.getpid(), purpose="probe")
        stage = _resolve_tune_dir(cdir=cdir, tunename=tunename)
        try:
            original_config = self.create_steering_card(
                param_space=selected_config, tunename=tunename, cdir=cdir, tune_default=str(target))
            old_card = self._decoded_card_config(config=original_config, path=str(target))
            new_card = self._decoded_card_config(config=selected_config, path=str(stage))
            expected = summary.get("card_config")
            if isinstance(expected, dict):
                expected_values = icetune_push.flatten_card_config(expected)
                generated_values = icetune_push.flatten_card_config(new_card)
                if not generated_values.keys() <= expected_values.keys() or any(not math.isclose(
                        float(expected_values[key]), float(generated_values[key]), rel_tol=1.0e-10, abs_tol=1.0e-12)
                    for key in generated_values):
                    raise ValueError("GRANIITTI summary card_config does not match its optimizer config")

            rendered: dict[pathlib.Path, str] = {}
            new_fields: list[tuple[pathlib.Path, tuple]] = []
            for staged_path in sorted(stage.rglob("*.json")):
                relative = staged_path.relative_to(stage)
                # Include source cards changed through references from the selected parameters
                source_path = target / relative
                if not source_path.is_file():
                    raise ValueError(f'GRANIITTI push target is missing "{relative}"')
                desired = json5.loads(staged_path.read_text(encoding="utf-8"))
                source = source_path.read_text(encoding="utf-8")
                allowed_additions = set()
                if relative.name in {"CON_MP.json", "CON_XP.json", "CON_GP.json", "CON_TP.json"}:
                    parameter_class = relative.stem
                    for parameter in selected_config:
                        parts = tools.substring_split(s=parameter, delimiters=["|", ":"])
                        if (len(parts) == 4 and parts[0] == parameter_class
                            and parts[3].split(".", 1)[0] in {"FF_transfer", "FF_offshell"}):
                            allowed_additions.add((parts[1], parts[2], parts[3].split(".", 1)[0]))
                updated, changed = icetune_push.render_json5_scalars(
                    source, desired, allowed_object_additions=allowed_additions)
                new_fields.extend((relative, path) for path in changed if path in allowed_additions)
                if changed:
                    rendered[source_path] = updated
            rows = icetune_push.card_config_rows(old_card, new_card)
            for relative, path in new_fields:
                field = "".join(f"[{json.dumps(component)}]" for component in path)
                print(f"New steering field: {relative}{field}")
            if not confirm(rows):
                raise icetune_push.PushCancelled
            icetune_push.atomic_write_texts(rendered)
            return rows
        finally:
            _cleanup_temporary_tune_dir(cdir=cdir, tunename=tunename)

    # Update one scalar or indexed value in a card block
    def update_target(self, parname: str, new_value: float, target: dict) -> float | str | int:
        if "@" in parname:
            raise Exception(__name__ + f".update_target: Unsupported parameter decorator in {parname}")

        tokens = jsonref.field_parts(parname)

        container = target
        try:
            for token in tokens[:-1]:
                container = container[token]
            final = tokens[-1]
            old_value = container[final]
        except (KeyError, IndexError, TypeError):
            logger.error(".update_target: unknown parameter %s", parname)
            logger.debug(".update_target: target=%s", target)
            raise Exception(__name__ + f".update_target: Unknown parameter {parname}") from None
        if isinstance(old_value, bool) and isinstance(new_value, (int, float)):
            new_value = bool(new_value)
        container[final] = new_value
        return old_value

    # Update one power transfer form factor through its inverse scale
    def _update_transfer_inv2(self, target: dict, new_value: float) -> float:
        current = target["FF_transfer"]
        if not isinstance(current, dict):
            raise ValueError("Transfer form factor is missing")
        if current.get("type") == "none":
            raise ValueError("Inactive transfer form factors are not icetune parameters")
        if current.get("type") != "power":
            raise ValueError(f"Transfer inverse scale requires a power form, got {current.get('type')}")
        old_value = 1.0 / float(current["Lambda2"])

        inv2 = float(new_value)
        if not math.isfinite(inv2) or inv2 <= 0.0:
            raise ValueError(f"Transfer inverse scale must be finite and positive, got {inv2}")
        current["Lambda2"] = 1.0 / inv2
        return old_value

    # Compute one minimal continuum subchannel form factor object
    def _continuum_ff_template(self, ff_type: str, kernels: int = 1) -> dict:
        if ff_type == "exp":
            return {"type": ff_type, "norm": "pole", "b": 1.0}
        if ff_type == "logexp":
            return {"type": ff_type, "norm": "pole", "b": 1.0, "Lambda2": 1.0}
        if ff_type == "power":
            return {"type": ff_type, "norm": "pole", "Lambda2": 1.0, "n": 1.0}
        if ff_type == "orear":
            return {"type": ff_type, "norm": "pole", "b": 1.0, "a": 0.25}
        if ff_type == "gkernel" and kernels >= 1:
            term = {"a": 1.0, "p": 1.0, "nu": 0.0, "mu2": 0.0}
            return {"type": ff_type, "norm": "pole", "terms": [copy.deepcopy(term) for _ in range(kernels)]}
        raise ValueError(f"Unsupported continuum FF_offshell type {ff_type}")

    # Update one value inside a named card block
    def update_block(self, parname: str, blockname: str, new_value: float, json_data: dict) -> float | str | int:
        target = json_data
        for key in blockname.split("."):
            target = target[key]
        return self.update_target(parname=parname, new_value=new_value, target=target)

    # Compute one charge-sector block from a central continuum card
    def production_sector_block(self, mother_block: dict, channel: str) -> dict:
        if "/" in channel:
            family, sector = channel.split("/", 1)
            if family not in mother_block:
                raise Exception(__name__ + f".production_sector_block: Unknown production family {family}")
            family_block = mother_block[family]
            sectors = [name for name in ("same", "opposite", "self") if name in family_block]
            if "self" in sectors and len(sectors) != 1:
                raise Exception(
                    __name__ + f".production_sector_block: Production family {family} mixes self with same or opposite")
            if sector not in sectors:
                raise Exception(__name__ + f".production_sector_block: Unknown production sector {family}/{sector}")
            return family_block[sector]

        if channel not in mother_block:
            raise Exception(__name__ + f".production_sector_block: Unknown production family {channel}")

        family_block = mother_block[channel]
        sectors = [sector for sector in ("same", "opposite", "self") if sector in family_block]
        if len(sectors) != 1:
            raise Exception(__name__
                + f".production_sector_block: Production family {channel} has sectors {sectors}; use family/sector in CON model parameter keys"
            )
        return family_block[sectors[0]]

    # Split one resonance tuning key into model and parameter names
    def _split_resonance_param(self, par: str) -> tuple[str | None, str]:
        remainder = par.partition(":")[2]
        parts = remainder.split(":", 1)
        if len(parts) == 2 and parts[0] in {"XP", "GP", "MP", "TP"}:
            return parts[0], parts[1]
        return None, remainder

    # Compute the resonance card name from one tuning key
    def _resonance_name_from_key(self, par: str) -> str:
        if not par.startswith("RES|"):
            raise Exception(__name__ + f"._resonance_name_from_key: Not a resonance parameter key: {par}")
        remainder = par.partition("|")[2]
        return remainder.split(":", 1)[0]

    # Compute the active eta-style inline block for central resonance steering
    def _resonance_block(self, json_data: dict, model: str) -> dict:
        model_block = json_data["PARAM_RES"]["MODELS"][model]
        keys = [key for key in model_block if re.fullmatch(r"\[-?\d+,-?\d+\]", key)]
        if len(keys) != 1:
            raise Exception(__name__ + f"._resonance_block: {model} central steering expects one active channel")
        return model_block[keys[0]]

    # Resolve the card block targeted by one resonance parameter
    def _resonance_target(self, json_data: dict, model: str | None, parname: str) -> dict:
        if model is None:
            return json_data["PARAM_RES"]
        if model == "TP" and parname == "BW":
            raise ValueError("TP resonances do not use BW steering")

        models = json_data["PARAM_RES"]["MODELS"]
        if model in {"XP", "GP", "MP"}:
            if parname == "g" or parname.startswith(("g[", "g_ls", "helicity")) or parname in {"basis", "Lambda"}:
                return self._resonance_block(json_data=json_data, model=model)
            if model == "MP" and (parname.startswith("polarization.") or parname == "basis"):
                return self._resonance_block(json_data=json_data, model="MP")
        return models[model]

    # Fix the common phase of one decay shape through its first active row
    def _canonicalize_relative_shape(self, block: dict) -> None:
        basis = block.get("basis")
        if basis not in {"alpha_ls", "helicity"}:
            return
        rows = block[basis]
        reference = next((index for index, row in enumerate(rows) if float(row[2]) > 0.0), None)
        if reference is None:
            raise Exception(__name__ + "._canonicalize_relative_shape: zero decay shape")
        reference_magnitude = float(rows[reference][2])
        if not math.isfinite(reference_magnitude) or reference_magnitude <= 0.0:
            raise Exception(__name__ + "._canonicalize_relative_shape: zero reference")
        phase_shift = float(rows[reference][3]) if len(rows[reference]) > 3 else 0.0
        for row in rows:
            relative_phase = float(row[3]) if len(row) > 3 else 0.0
            relative_phase = tools.canonicalize_phase(relative_phase - phase_shift)
            if abs(relative_phase) <= 1.0e-12:
                relative_phase = 0.0
            if len(row) > 3:
                row[3] = relative_phase
            elif relative_phase != 0.0:
                row.append(relative_phase)

    # Materialize both MP polarization coordinate systems for every mode
    def _canonicalize_mp_polarization(self, json_data: dict) -> None:
        block = self._resonance_block(json_data=json_data, model="MP")
        polarization = block["polarization"]
        mode = polarization["mode"]
        basis = block["basis"]
        if basis not in {"auto_min_L", "auto_min_S", "auto_equal_ls", "auto_equal_helicity", "g_ls", "helicity"}:
            raise ValueError(f"Unsupported MP production basis: {basis}")
        field = "g" if basis.startswith("auto_") else basis
        if field not in block:
            raise ValueError(f"MP basis {basis} requires {field}")
        for inactive in {"g", "g_ls", "helicity"} - {field}:
            block.pop(inactive, None)
        if "Lambda" not in block or not math.isfinite(block["Lambda"]) or block["Lambda"] <= 0.0:
            raise ValueError("MP production requires finite positive Lambda")
        unknown = polarization.keys() - {"mode", "a_Jz", "rho_mag", "rho_phase", "random_rho"}
        if unknown:
            raise ValueError(f"Unknown MP polarization fields: {sorted(unknown)}")
        spin_x2 = int(json_data["PARAM_RES"]["spinX2"])
        if spin_x2 < 0 or spin_x2 % 2 != 0:
            raise Exception(__name__ + "._canonicalize_mp_polarization: MP requires integer resonance spin")
        j = spin_x2 // 2
        size = spin_x2 + 1

        if mode not in {"none", "a_Jz", "rho"}:
            raise Exception(__name__ + f'._canonicalize_mp_polarization: Unknown polarization mode "{mode}"')
        if "a_Jz" not in polarization:
            polarization["a_Jz"] = (
                [[0, 1.0, 0.0]] if j == 0 else [[jz, 1.0 if jz == 0 else 0.0, 0.0] for jz in range(-j, 1)])
        if "rho_mag" not in polarization:
            polarization["rho_mag"] = (np.eye(size, dtype=float) / size).tolist()
        if "rho_phase" not in polarization:
            polarization["rho_phase"] = np.zeros((size, size), dtype=float).tolist()
        polarization.setdefault("random_rho", False)

    # Compute one coherent MP direction angle key
    def _ajz_angle_key(self, res: str, k: int, geometry: str) -> str:
        base = f"RES|{res}:MP:polarization.a_Jz"
        if geometry == "projective":
            return tools.ajzp_angle_key(base, k)
        if geometry == "sphere":
            return tools.spherical_angle_key(base, k)
        raise ValueError(f'Unknown coherent MP direction geometry "{geometry}"')

    # Compute whether one parameter is a coherent MP direction angle
    def _is_ajz_angle_key(self, key: str) -> bool:
        if tools.is_ajzp_angle_key(key):
            return True
        if not tools.is_spherical_angle_key(key):
            return False
        base, _ = tools.parse_spherical_angle_key(key)
        return base.endswith(":MP:polarization.a_Jz")

    # Identify a relative coherent spin phase after common phase decoding
    def _is_ajz_phase_key(self, key: str) -> bool:
        return re.fullmatch(r"RES\|[^:]+:MP:polarization\.phi_Jz\[\d+\]", key) is not None

    # Compute validated coherent a_Jz rows from the active RES card block
    def _ajzp_card_rows(self, res: str, j: int, json_data: dict) -> list[tuple[int, float, float]]:
        rows = []
        seen = set()
        res_block = self._resonance_block(json_data=json_data, model="MP")
        for entry in res_block["polarization"].get("a_Jz", []):
            if len(entry) != 3:
                raise Exception(__name__ + f'._ajzp_card_rows: Bad a_Jz entry for "{res}"')
            if not isinstance(entry[0], int) or isinstance(entry[0], bool):
                raise Exception(__name__ + f'._ajzp_card_rows: Non-integer Jz entry for "{res}"')
            jz, mag, phase = int(entry[0]), float(entry[1]), float(entry[2])
            if jz < -j or jz > j:
                raise Exception(__name__ + f'._ajzp_card_rows: Jz outside [-J,J] for "{res}"')
            if jz in seen:
                raise Exception(__name__ + f'._ajzp_card_rows: Duplicate Jz entry for "{res}"')
            if not math.isfinite(mag) or not math.isfinite(phase) or mag < 0.0:
                raise Exception(__name__ + f'._ajzp_card_rows: Bad a_Jz magnitude or phase for "{res}"')
            seen.add(jz)
            rows.append((jz, mag, phase))
        return rows

    # Compute a complete in-memory coherent a_Jz map with parity partners filled
    def _complete_ajzp_values(self, res: str, j: int, json_data: dict) -> dict:
        a_map = {jz: (mag, phase) for jz, mag, phase in self._ajzp_card_rows(res, j, json_data)}

        if not a_map:
            sector_labels = list(range(j + 1))
            reference = sector_labels[0]
            a_map = self._ajzp_values_from_parameters(j=j, sector_labels=sector_labels,
                sectors=[1.0 if label == reference else 0.0 for label in sector_labels])

        for m in range(1, j + 1):
            if -m in a_map and m not in a_map:
                mag, phase = a_map[-m]
                a_map[m] = mag, tools.canonicalize_phase(phase + (math.pi if m % 2 else 0.0))
            elif -m not in a_map and m in a_map:
                mag, phase = a_map[m]
                a_map[-m] = mag, tools.canonicalize_phase(phase + (math.pi if m % 2 else 0.0))

        # Sparse omitted projections are zero, as in the native resonance reader
        for m in range(-j, j + 1):
            a_map.setdefault(m, (0.0, 0.0))

        return a_map

    # Compute a normalized in-memory coherent a_Jz map from @AJZP signed sectors
    def _ajzp_values_from_parameters(self, j: int, sector_labels: list[int], sectors: list[float]
    ) -> dict[int, tuple[float, float]]:
        a_map = {jz: (0.0, 0.0) for jz in range(-j, j + 1)}
        for m, sector in zip(sector_labels, sectors, strict=False):
            coefficient = sector if m == 0 else sector / math.sqrt(2.0)
            reference = (abs(coefficient), tools.real_phase(coefficient))
            if m == 0:
                a_map[0] = reference
                continue
            a_map[-m] = reference
            a_map[m] = reference[0], tools.canonicalize_phase(reference[1] + (math.pi if m % 2 else 0.0))
        return a_map

    # Compute default coherent MP angles using the spin-basis reference sector
    def _default_ajz_values(self, res: str, j: int, json_data: dict, geometry: str, *, complex_phases=False, frame="CS") -> dict:
        a_map = self._complete_ajzp_values(res=res, j=j, json_data=json_data)

        sector_labels = mp_spin_sectors(j, frame)
        reference = sector_labels[0]
        ordered_labels = [reference, *[label for label in sector_labels if label != reference]]
        sectors = []
        for m in ordered_labels:
            sign = 1.0 if complex_phases or math.cos(a_map[-m][1]) >= 0.0 else -1.0
            scale = 1.0 if m == 0 else math.sqrt(2.0)
            sectors.append(sign * scale * a_map[-m][0])

        default = {}
        angles = (tools.projective_angles_from_vector(sectors) if geometry == "projective"
            else tools.spherical_angles_from_vector(sectors))
        for k, angle in enumerate(angles):
            default[self._ajz_angle_key(res=res, k=k, geometry=geometry)] = angle
        return default

    # Decode normalized coherent spin sectors without discarding derivatives at zero
    def _ajz_sectors(self, *, res: str, j: int, geometry: str, param: dict, json_data: dict) -> dict:
        labels = mp_spin_sectors(j, param.get("REGGE|MP_FRAME", "CS"))
        reference = labels[0]
        ordered = [reference, *[m for m in labels if m != reference]]
        angles = [param[self._ajz_angle_key(res=res, k=k, geometry=geometry)] for k in range(len(labels) - 1)]
        values = (tools.projective_vector_from_angles(angles)
            if geometry == "projective" else tools.spherical_vector_from_angles(angles))
        prefix = f"RES|{res}:MP:polarization.phi_Jz"
        phased = any(key.startswith(prefix) for key in param)
        phase0 = self._complete_ajzp_values(res, j, json_data)[reference][1] if phased else 0.0
        sectors = {}
        for m, value in zip(ordered, values, strict=True):
            phase = array.asarray(phase0 + param.get(f"{prefix}[{m}]", 0.0), dtype=float)
            xp = array.namespace(phase)
            sectors[m] = value * (xp.cos(phase) + 1j * xp.sin(phase))
        return sectors

    # Apply grouped coherent MP direction parameters while preserving explicit rows
    def _update_ajzp_param(self, res: str, res_params: list, param: dict, json_data: dict, orig_values: dict) -> None:
        ajzp_params = [par for par in res_params if self._is_ajz_angle_key(par)]
        phase_params = [par for par in res_params if self._is_ajz_phase_key(par)]
        if not ajzp_params:
            if phase_params:
                raise ValueError(f'Coherent spin phases require direction angles for "{res}"')
            return
        if self._resonance_block(json_data=json_data, model="MP")["polarization"]["mode"] != "a_Jz":
            raise ValueError("MP spin coordinates require polarization.mode a_Jz")
        geometries = {"projective" if tools.is_ajzp_angle_key(par) else "sphere" for par in ajzp_params}
        if len(geometries) != 1:
            raise Exception(__name__ + f'._update_ajzp_param: Mixed direction geometries for "{res}"')
        geometry = next(iter(geometries))

        direct_ajz = []
        for par in res_params:
            if self._is_ajz_angle_key(par):
                continue
            model, parname = self._split_resonance_param(par)
            if model not in (None, "MP"):
                continue
            if parname.startswith("polarization.a_Jz"):
                direct_ajz.append(par)

        if direct_ajz:
            raise Exception(
                __name__ + f'._update_ajzp_param: Cannot mix direct a_Jz edits and direction angles for "{res}"')

        spin_x2 = int(json_data["PARAM_RES"]["spinX2"])
        if spin_x2 % 2 != 0:
            raise Exception(__name__
                + f'._update_ajzp_param: Coherent direction angles require integer spin, got spinX2={spin_x2} for "{res}"'
            )
        j = spin_x2 // 2
        frame = param.get("REGGE|MP_FRAME", "CS")
        sector_labels = mp_spin_sectors(j, frame)

        expected_angles = {self._ajz_angle_key(res=res, k=k, geometry=geometry) for k in range(len(sector_labels) - 1)}
        found = set()
        for par in ajzp_params:
            if par not in expected_angles:
                _, parname = self._split_resonance_param(par)
                raise Exception(__name__
                    + f'._update_ajzp_param: Unexpected direction parameter "{parname}" for J={j} resonance "{res}"')
            found.add(par)

        if found != expected_angles:
            missing = sorted(expected_angles - found)
            raise Exception(
                __name__ + f'._update_ajzp_param: Incomplete direction angle set for "{res}" (missing={missing})')
        reference = sector_labels[0]
        expected_phases = {f"RES|{res}:MP:polarization.phi_Jz[{m}]" for m in sector_labels if m != reference}
        if phase_params and set(phase_params) != expected_phases:
            raise ValueError(f'Coherent spin phases must cover every nonreference sector for "{res}"')
        defaults = self._default_ajz_values(
            res=res, j=j, json_data=json_data, geometry=geometry, complex_phases=bool(phase_params), frame=frame)
        for par in sorted(found):
            orig_values[par] = defaults[par]
        original = self._complete_ajzp_values(res, j, json_data)
        for m in sector_labels:
            key = f"RES|{res}:MP:polarization.phi_Jz[{m}]"
            if key in phase_params:
                orig_values[key] = tools.canonicalize_phase(original[m][1] - original[reference][1])
        sectors = self._ajz_sectors(res=res, j=j, geometry=geometry, param=param, json_data=json_data)
        a_map = {m: tools.decode_cartesian(value.real, value.imag)
            for m, value in ((m, (-1 if m > 0 and m % 2 else 1) * sectors.get(abs(m), 0j) / math.sqrt(1.0 if m == 0 else 2.0)) for m in range(-j, j + 1))
        }
        if phase_params:
            for m, (magnitude, _phase) in a_map.items():
                if magnitude <= tools.COMPLEX_EPS:
                    key = f"RES|{res}:MP:polarization.phi_Jz[{abs(m)}]"
                    a_map[m] = magnitude, tools.canonicalize_phase(original[reference][1] + param.get(key, 0.0))
        a_jz = [[m, a_map[m][0], a_map[m][1]] for m in range(-j, 1)]

        polarization = self._resonance_block(json_data=json_data, model="MP")["polarization"]
        polarization["mode"] = "a_Jz"
        polarization["a_Jz"] = a_jz

    # Compute the base name and integer indices of one bracketed parameter
    def _indexed_parameter(self, parname: str) -> tuple[str, list[int]] | None:
        if "[" not in parname or "]" not in parname:
            return None
        indices = tools.string_to_index(s=tools.find_str_between(s=parname, start="[", end="]"), delim=",")
        return parname.partition("[")[0], indices

    # Compute magnitude and phase columns for one angular coupling table
    def _angular_columns(self, key: str) -> tuple[int, int]:
        parts = tools.substring_split(s=key, delimiters=["|", ":"])
        gp_nonphoton = len(parts) > 1 and parts[0] == "CON_GP" and int(parts[1]) != 22
        magnitude = 3 if gp_nonphoton else 2
        return magnitude, magnitude + 1

    # Read a parameter card from disk or an immutable decoder snapshot
    def _parameter_card(self, path, filename):
        return copy.deepcopy(path[filename]) if isinstance(path, dict) else load_json_file(
            os.path.join(path, filename), loader=json5.load)

    # Resolve one coupling table in a resonance, decay or continuum card
    def _coupling_table(self, base_key: str, path: str) -> dict:
        parts = tools.substring_split(s=base_key, delimiters=["|", ":"])
        cls, field = parts[0], parts[-1]
        if cls == "RES" and len(parts) in {4, 5}:
            data = self._parameter_card(path, f"RES/{parts[1]}.json")
            target = (data["PARAM_RES"]["MODELS"]["TP"][parts[3]] if len(parts) == 5 and parts[2] == "TP"
                else self._resonance_target(json_data=data, model=parts[2], parname=field))
        elif cls == "DECAY" and len(parts) == 4:
            data = self._parameter_card(path, "DECAYS.json")
            target = data[parts[1]][parts[2]]
        elif cls in {"CON_MP", "CON_XP", "CON_GP", "CON_TP"} and len(parts) == 4:
            data = self._parameter_card(path, f"{cls}.json")
            target = (data[parts[1]][parts[2]] if cls == "CON_TP"
                else self.production_sector_block(mother_block=data[parts[1]], channel=parts[2]))
        else:
            raise ValueError(f'Unsupported coupling table "{base_key}"')
        rows = copy.deepcopy(target.get(field))
        table = {"basis": field, "rows": rows}
        if field == "g_tensor":
            table["names"] = tp_res_names(data["PARAM_RES"]) if cls == "RES" else tp_con_names(rows)
        return table

    # Compute the physical angular table selected by one parameter name
    def _angular_card_rows(self, *, base_key: str, path: str) -> tuple[str, list, int]:
        parts = tools.substring_split(s=base_key, delimiters=["|", ":"])
        if len(parts) != 4:
            raise ValueError(f'Invalid physical angular parameter "{base_key}"')
        field_match = re.fullmatch(r"(g_ls|alpha_ls|helicity)\([^()]+\)", parts[-1])
        if field_match is None:
            raise ValueError(f'Invalid physical angular parameter "{base_key}"')
        field = field_match.group(1)
        if parts[0] not in {"RES", "DECAY", "CON_MP", "CON_XP", "CON_GP"}:
            raise ValueError(f'Unsupported physical angular parameter "{base_key}"')
        rows = self._coupling_table(f"{base_key.rsplit(':', 1)[0]}:{field}", path)["rows"]
        if not isinstance(rows, list) or not rows:
            raise ValueError(f'Physical angular parameter has no {field} table: "{base_key}"')
        magnitude_column, _ = self._angular_columns(base_key)
        row_size = magnitude_column + 2
        if any(not isinstance(row, list) or len(row) != row_size for row in rows):
            raise ValueError(f'Physical angular parameter has invalid {field} rows: "{base_key}"')
        return field, rows, magnitude_column

    # Convert physical coupling names to private steering card coordinates
    def _physical_params(self, param: dict, *, path: str) -> tuple[dict, dict]:
        suffixes = (tools.PHASE_U_SUFFIX, tools.PHASE_V_SUFFIX, tools.PROJECTIVE_VECTOR_SUFFIX,
            tools.SPHERICAL_VECTOR_SUFFIX, tools.COMPLEX_RE_SUFFIX, tools.COMPLEX_IM_SUFFIX, tools.MAGNITUDE_SUFFIX,
            tools.RAW_PHASE_SUFFIX)
        transformed = {}
        exact = {}
        components = {}
        table_cache = {}
        physical = re.compile(r"^(.*:)(g_ls|alpha_ls|helicity)\([^()]+\)$")

        for key, value in param.items():
            name = str(key)
            suffix = next((item for item in suffixes if name.endswith(item)), None)
            base_key = name if suffix is None else name[: -len(suffix)]
            field = base_key.rsplit(":", 1)[-1]
            tensor = field == "g_tensor" or re.fullmatch(r"g_tensor\(.*\)", field) is not None
            match = physical.fullmatch(base_key)
            if not tensor and match is None:
                transformed[key] = value
                continue
            if tensor:
                if suffix in {tools.PROJECTIVE_VECTOR_SUFFIX, tools.SPHERICAL_VECTOR_SUFFIX}:
                    transformed[key] = value
                    continue
                if suffix is not None:
                    raise ValueError(f'Tensor coupling has no supported component marker: "{key}"')
                normalized = self._tensor_coupling_card_key(base_key, path=path)
            else:
                if suffix is None:
                    raise ValueError(f'Physical angular parameter requires a component marker: "{key}"')
                table_key = f"{match.group(1)}{match.group(2)}"
                if table_key not in table_cache:
                    table_cache[table_key] = self._angular_card_rows(base_key=base_key, path=path)
                field, rows, magnitude_column = table_cache[table_key]
                row_names = [angular_row_name(field, row) for row in rows]
                if len(row_names) != len(set(row_names)):
                    raise ValueError(f'Duplicate physical angular rows for "{base_key}"')
                requested = base_key.rsplit(":", 1)[1]
                raw_base = f"{match.group(1)}{field}"
                phase_column = magnitude_column + 1
                for index, row_name in enumerate(row_names):
                    public_base = f"{match.group(1)}{row_name}"
                    components[f"{raw_base}[{index},{magnitude_column}]"] = tools.magnitude_key(public_base)
                    components[f"{raw_base}[{index},{phase_column}]"] = tools.raw_phase_key(public_base)
                equal_m_suffix = ",m=-2,-1,0)"
                if requested.endswith(equal_m_suffix):
                    if suffix != tools.MAGNITUDE_SUFFIX or field != "helicity":
                        raise ValueError(f'Equal m coupling requires a helicity magnitude: "{key}"')
                    prefix = requested[: -len(equal_m_suffix)]
                    selected_names = [f"{prefix},{m})" for m in (-2, -1, 0)]
                    if any(row_name not in row_names for row_name in selected_names):
                        raise ValueError(f'Unknown equal m helicity rows for "{base_key}"')
                    for row_name in selected_names:
                        row = row_names.index(row_name)
                        raw_mag = f"{raw_base}[{row},{magnitude_column}]"
                        if raw_mag in transformed:
                            raise ValueError(f'Duplicate physical coupling target "{key}"')
                        transformed[raw_mag] = value
                        exact[raw_mag] = key
                        components[raw_mag] = key
                    continue
                if requested not in row_names:
                    raise ValueError(f'Unknown physical angular row "{base_key}"')
                row = row_names.index(requested)
                raw_mag = f"{raw_base}[{row},{magnitude_column}]"
                raw_phase = f"{raw_base}[{row},{phase_column}]"
                normalized = {tools.MAGNITUDE_SUFFIX: raw_mag, tools.RAW_PHASE_SUFFIX: tools.raw_phase_key(raw_phase),
                    tools.PHASE_U_SUFFIX: tools.phase_u_key(raw_phase),
                    tools.PHASE_V_SUFFIX: tools.phase_v_key(raw_phase),
                    tools.COMPLEX_RE_SUFFIX: tools.complex_re_key(f"{raw_base}[{row}]"),
                    tools.COMPLEX_IM_SUFFIX: tools.complex_im_key(f"{raw_base}[{row}]")}.get(suffix)
                if normalized is None:
                    if row == 0:
                        raise ValueError(f'Direction reference row has no optimizer angle: "{base_key}"')
                    normalized = (tools.spherical_angle_key(raw_base, row - 1)
                        if suffix == tools.SPHERICAL_VECTOR_SUFFIX else tools.projective_angle_key(raw_base, row - 1))
            if normalized in transformed:
                raise ValueError(f'Duplicate physical coupling target "{key}"')
            transformed[normalized] = value
            exact[normalized] = key

        return transformed, {"components": components, "exact": exact}

    # Restore public coupling names after optimizer coordinate encoding
    def _restore_physical_optimizer_keys(self, values: dict, mapping: dict) -> dict:
        exact = mapping.get("exact", {})
        return {exact.get(key, key): value for key, value in values.items()}

    # Restore physical names in decoded steering card coordinates
    def _restore_physical_card_keys(self, values: dict, mapping: dict) -> dict:
        names = mapping.get("exact", {}) | mapping.get("components", {})
        return {names.get(key, key): value for key, value in values.items()}

    # Resolve one physical Tensor Pomeron coefficient to its card index
    def _tensor_coupling_card_key(self, key: str, *, path: str) -> str:
        parts = tools.substring_split(s=key, delimiters=["|", ":"])
        supported = (
            parts[0] == "RES" and len(parts) == 5 and parts[2] == "TP" or parts[0] == "CON_TP" and len(parts) == 4)
        if not supported:
            raise ValueError(f'Unsupported physical TP coupling "{key}"')
        table = self._coupling_table(f"{key.rsplit(':', 1)[0]}:g_tensor", path)
        couplings, names = table["rows"], table["names"]
        if not isinstance(couplings, list) or len(couplings) != len(names):
            raise ValueError(f'Invalid TP g_tensor table for "{key}"')
        if parts[-1] not in names:
            raise ValueError(f'Unknown physical TP coupling "{key}"')
        return f"{key.rsplit(':', 1)[0]}:g_tensor[{names.index(parts[-1])}]"

    # Compute true for a scalar phase that may use raw or Cayley coordinates
    def _is_scalar_phase_target(self, par: str) -> bool:
        cls = par.partition("|")[0]
        parname = tools.substring_split(s=par, delimiters=["|", ":"])[-1]
        if cls == "DECAY":
            return parname in {"zeta.MP", "zeta.XP", "zeta.GP", "zeta.TP"}
        if cls == "RES":
            model, parname = self._split_resonance_param(par)
            if self._is_ajz_phase_key(par):
                return True
            if model in {"MP", "XP", "GP", "TP"} and parname == "phi":
                return True
            if model not in {"MP", "XP", "GP"}:
                return False
        elif cls not in {"CON_MP", "CON_XP", "CON_GP"}:
            return False
        if parname == "g[1]":
            return True
        target = self._indexed_parameter(parname)
        _, phase_column = self._angular_columns(par)
        return bool(target is not None and target[0] in {"g_ls", "helicity"} and len(target[1]) == 2
            and target[1][1] == phase_column)

    # Compute true for normalized signed projective vector base keys
    def _is_direction_target(self, base_key: str) -> bool:
        cls = base_key.partition("|")[0]
        if cls == "DECAY":
            parname = tools.substring_split(s=base_key, delimiters=["|", ":"])[-1]
            return parname in {"alpha_ls", "helicity"}
        if cls == "RES":
            parts = tools.substring_split(s=base_key, delimiters=["|", ":"])
            return (len(parts) == 4 and parts[2] in {"XP", "GP"} and parts[3] == "g_ls") or (
                len(parts) == 5 and parts[2] == "TP" and parts[4] == "g_tensor")
        if cls in {"CON_MP", "CON_XP", "CON_GP"}:
            parts = tools.substring_split(s=base_key, delimiters=["|", ":"])
            return len(parts) == 4 and parts[3] in {"helicity", "g_ls"}
        return False

    # Compute raw magnitude and phase keys for one projective vector row
    def _projective_vector_row_keys(self, base_key: str, row: int) -> tuple[str, str]:
        if base_key.endswith(":g_tensor"):
            return f"{base_key}[{row}]", ""
        magnitude, phase = self._angular_columns(base_key)
        return f"{base_key}[{row},{magnitude}]", f"{base_key}[{row},{phase}]"

    # Compute true for an observable complex channel coupling
    def _is_channel_coupling_target(self, par: str) -> bool:
        cls = par.partition("|")[0]
        if cls == "RES":
            model, parname = self._split_resonance_param(par)
            if model not in {"MP", "XP", "GP"}:
                return False
        elif cls in {"CON_MP", "CON_XP", "CON_GP"}:
            parname = tools.substring_split(s=par, delimiters=["|", ":"])[-1]
        else:
            return False
        target = self._indexed_parameter(parname)
        return parname == "g" or bool(target is not None and target[0] in {"helicity", "g_ls"} and len(target[1]) == 1)

    # Convert an encoded coupling base key to raw magnitude and phase keys
    def _raw_channel_coupling_keys(self, base_key: str) -> tuple[str, str]:
        if base_key.endswith(":g"):
            return f"{base_key}[0]", f"{base_key}[1]"
        if "[" not in base_key or "]" not in base_key:
            raise Exception(__name__ + f'._raw_channel_coupling_keys: Missing row index in "{base_key}"')
        prefix, row_token = base_key.rsplit("[", 1)
        row = row_token.rstrip("]")
        magnitude_column, phase_column = self._angular_columns(base_key)
        return (f"{prefix}[{row},{magnitude_column}]", f"{prefix}[{row},{phase_column}]")

    # Collect complete canonical eikonal coordinate groups
    def _collect_eikonal_canonical_groups(self, param: dict) -> list[_Coordinates]:
        grouped = {}
        definitions = ((eikonal.is_descending_coupling_key, eikonal.descending_coupling_target, "descending"),
            (eikonal.is_symmetric_coupling_key, eikonal.symmetric_coupling_target, "symmetric"),
            (eikonal.is_ordered_dpow_key, eikonal.ordered_dpow_target, "DPOW"))
        for key, value in param.items():
            definition = next((item for item in definitions if item[0](key)), None)
            if definition is None:
                continue
            _, target_key, coord_type = definition
            raw = target_key(key)
            target = self._indexed_parameter(raw)
            if target is None or len(target[1]) != 2 or not target[0].startswith("SOFT|"):
                raise ValueError(f'Invalid {coord_type} target "{raw}"')
            base, (row, column) = target
            matrix = ":EXCHANGE." in base and base.endswith(".g")
            valid = {"descending": matrix and row == column, "symmetric": matrix and row < column,
                "DPOW": ":FF." in base and base.endswith(".param") and column in {0, 1}}[coord_type]
            if not valid:
                raise ValueError(f'Invalid {coord_type} target "{raw}"')
            mirror = f"{base}[{column},{row}]" if coord_type == "symmetric" else None
            if raw in param or mirror in param:
                raise ValueError(f'Mixed raw and {coord_type} target "{raw}"')
            group_key = (coord_type, raw if mirror else base, row if coord_type == "DPOW" else None)
            grouped.setdefault(group_key, []).append((column, key, raw, value, mirror))

        groups = []
        for (coord_type, base, _), entries in grouped.items():
            entries.sort(key=lambda entry: entry[0])
            if coord_type == "symmetric":
                _, key, raw, value, mirror = entries[0]
                groups.append(_Coordinates((key,), (raw, mirror), (value, value), eikonal.encode_symmetric_coupling))
                continue
            if coord_type == "descending":
                subblock = base.split("|", 1)[1].split(":", 1)[0]
                prefix, separator, model = subblock.partition(".")
                if prefix != "MODEL" or not separator:
                    raise ValueError(f'Invalid model target "{base}"')
                count = {"single": 1, "double": 2, "triple": 3}.get(model)
                decode, encode = (eikonal.decode_descending_couplings, eikonal.encode_descending_couplings)
            else:
                count = 2
                limits = {entry[1].rsplit(eikonal.ORDERED_DPOW_SUFFIX, 1)[1] for entry in entries}
                if len(limits) != 1:
                    raise ValueError(f'Inconsistent DPOW limits for "{base}"')
                limit = float(limits.pop())

                # Decode the two ordered DPOW scales with the saved domain limit
                def decode(values, limit=limit):
                    return eikonal.decode_ordered_dpow_scales(*values, limit=limit)

                # Encode the ordered scales with the same saved domain limit
                def encode(values, limit=limit):
                    return eikonal.encode_ordered_dpow_scales(*values, limit=limit)

            if count is None or [entry[0] for entry in entries] != list(range(count)):
                raise ValueError(f'Incomplete {coord_type} group "{base}"')
            groups.append(_Coordinates(tuple(entry[1] for entry in entries), tuple(entry[2] for entry in entries),
                    tuple(decode([entry[3] for entry in entries])), encode))
        return groups

    # Decode canonical eikonal optimizer coordinates into flat card targets
    def _prepare_eikonal_canonical_params(self, param: dict) -> tuple[dict, list[_Coordinates]]:
        groups = self._collect_eikonal_canonical_groups(param)
        return _decode_coordinates(param, groups), groups

    # Collect complete phase or complex coordinate groups through one suffix decoder
    def _collect_coordinate_groups(self, param: dict, coord_type: str) -> list[_Coordinates]:
        method = "_collect_coordinate_groups"
        is_phase = coord_type == "phase"
        if is_phase:
            is_encoded = tools.is_phase_encoded_key
            base_from_key = tools.encoded_phase_base_key
            is_target = self._is_scalar_phase_target
            components = {tools.PHASE_U_SUFFIX: "u", tools.PHASE_V_SUFFIX: "v"}
            target_label, pair_label = "encoded phase target", "encoded phase"
        elif coord_type == "complex":
            is_encoded = tools.is_complex_encoded_key
            base_from_key = tools.encoded_complex_base_key
            is_target = self._is_channel_coupling_target
            components = {tools.COMPLEX_RE_SUFFIX: "re", tools.COMPLEX_IM_SUFFIX: "im"}
            target_label, pair_label = "complex coupling target", "complex coupling"
        else:
            raise ValueError(f'Unknown coordinate group type "{coord_type}"')

        encoded = {}
        for key, value in param.items():
            if not is_encoded(key):
                continue
            base_key = base_from_key(key)
            if not is_target(base_key):
                raise Exception(f'{__name__}.{method}: Unsupported {target_label} "{base_key}"')
            field = next(name for suffix, name in components.items() if key.endswith(suffix))
            encoded.setdefault(base_key, {})[field] = value
        required = set(components.values())
        for base_key, values in encoded.items():
            if required - values.keys():
                raise Exception(f'{__name__}.{method}: Incomplete {pair_label} pair for "{base_key}"')

        groups = []
        for base_key, values in encoded.items():
            if is_phase:
                targets = (base_key,)
                keys = (tools.phase_u_key(base_key), tools.phase_v_key(base_key))
                decoded = (tools.decode_phase_components(values["u"], values["v"]),)

                # Encode the canonical scalar phase in its Cayley coordinates
                def encode(x):
                    return tools.encode_phase(tools.canonicalize_phase(x[0]))
            else:
                targets = self._raw_channel_coupling_keys(base_key)
                keys = (tools.complex_re_key(base_key), tools.complex_im_key(base_key))
                decoded = tools.decode_cartesian(values["re"], values["im"])

                # Encode the physical magnitude and phase in Cartesian coordinates
                def encode(x):
                    return tools.encode_polar(*x)

            if not param.keys().isdisjoint(targets):
                label = "raw and encoded phase" if is_phase else "raw and complex-encoded channel"
                raise Exception(f'{__name__}.{method}: Cannot mix {label} inputs for "{base_key}"')
            amplitude = None if is_phase else values["re"] + 1j * values["im"]
            groups.append(_Coordinates(keys, targets, tuple(decoded), encode, amplitude))

        if is_phase:
            for key, value in param.items():
                if not tools.is_raw_phase_key(key):
                    continue
                base_key = tools.raw_phase_base_key(key)
                if not self._is_scalar_phase_target(base_key):
                    raise Exception(f'{__name__}.{method}: Unsupported marked phase target "{base_key}"')
                if base_key in param or base_key in encoded:
                    raise Exception(f'{__name__}.{method}: Cannot mix phase encodings for "{base_key}"')
                groups.append(_Coordinates((key,), (base_key,), (tools.canonicalize_phase(value),),
                        lambda x: (tools.canonicalize_phase(x[0]),)))
        return groups

    # Collect complete signed hyperspherical projective angle groups
    def _direction_groups(self, param: dict, *, path: str) -> dict:
        angle_groups = {}

        for key, value in param.items():
            if not tools.is_projective_norm_key(key):
                continue
            base_key, reference = tools.parse_projective_norm_key(key)
            if not self._is_direction_target(base_key):
                raise Exception(__name__ + f'._direction_groups: Unsupported projective norm target "{base_key}"')
            entry = angle_groups.setdefault(base_key, {"angles": {}})
            if "norm_key" in entry:
                raise Exception(__name__ + f'._direction_groups: Duplicate projective norm for "{base_key}"')
            entry["norm_key"] = key
            entry["norm"] = value
            entry["reference"] = reference

        for key, value in param.items():
            if not tools.is_projective_coefficient_key(key):
                continue
            base_key, reference, component = tools.parse_projective_coefficient_key(key)
            if not self._is_direction_target(base_key):
                raise Exception(
                    __name__ + f'._direction_groups: Unsupported Cartesian projective coefficient target "{base_key}"')
            entry = angle_groups.setdefault(base_key, {"angles": {}})
            entry.setdefault("geometry", "projective")
            if "norm_key" in entry:
                raise Exception(
                    __name__ + f'._direction_groups: Mixed polar norm and Cartesian coefficient for "{base_key}"')
            if entry.setdefault("reference", reference) != reference:
                raise Exception(
                    __name__ + f'._direction_groups: Mixed Cartesian coefficient references for "{base_key}"')
            components = entry.setdefault("cartesian", {})
            if component in components:
                raise Exception(
                    __name__ + f'._direction_groups: Duplicate Cartesian coefficient component for "{base_key}"')
            components[component] = (key, value)

        for key, value in param.items():
            projective = tools.is_projective_angle_key(key)
            spherical = tools.is_spherical_angle_key(key)
            if not projective and not spherical:
                continue
            if spherical and self._is_ajz_angle_key(key):
                continue

            geometry = "projective" if projective else "sphere"
            base_key, angle_index = (
                tools.parse_projective_angle_key(key) if projective else tools.parse_spherical_angle_key(key))
            if not self._is_direction_target(base_key):
                raise Exception(__name__ + f'._direction_groups: Unsupported projective target "{base_key}"')

            entry = angle_groups.setdefault(base_key, {"angles": {}})
            if entry.setdefault("geometry", geometry) != geometry:
                raise Exception(__name__ + f'._direction_groups: Mixed direction geometries for "{base_key}"')
            if angle_index in entry["angles"]:
                raise Exception(__name__ + f'._direction_groups: Duplicate angle index {angle_index} for "{base_key}"')
            angle = value
            entry["angles"][angle_index] = (key, angle)

        for base_key, entry in angle_groups.items():
            if "cartesian" in entry:
                if set(entry["cartesian"]) != {"re", "im"}:
                    raise Exception(__name__ + f'._direction_groups: Incomplete Cartesian coefficient for "{base_key}"')
                re_key, re_value = entry["cartesian"]["re"]
                im_key, im_value = entry["cartesian"]["im"]
                entry["cartesian_keys"] = [re_key, im_key]
                entry["norm"], entry["residue_phase"] = tools.decode_cartesian(re_value, im_value)
            phase_key = tools.projective_phase_base_key(base_key)
            if "cartesian" in entry:
                if phase_key is None:
                    if base_key.startswith(("CON_MP|", "CON_XP|", "CON_GP|")):
                        entry["row_phase"] = True
                    else:
                        raise Exception(__name__
                            + f'._direction_groups: Cartesian coefficient has no separate common phase for "{base_key}"'
                        )
                elif phase_key in param:
                    raise Exception(__name__ + f'._direction_groups: Mixed polar and Cartesian phases for "{base_key}"')
                else:
                    entry["coupled_phase_key"] = phase_key
            elif entry.get("geometry") == "projective" and phase_key in param:
                entry["coupled_phase_key"] = phase_key
            found = set(entry["angles"])
            table = self._projective_card_table(base_key=base_key, path=path)
            field = str(table["basis"])
            common_phase = "coupled_phase_key" in entry or entry.get("row_phase", False)
            row_names = (
                list(table["names"]) if field == "g_tensor" else [angular_row_name(field, row) for row in table["rows"]]
            )
            row_count = len(table["rows"])
            if "reference" in entry:
                reference_name = f"{field}({entry['reference']})"
                if reference_name not in row_names:
                    raise Exception(__name__
                        + f'._direction_groups: Unknown projective reference "{reference_name}" for "{base_key}"')
                reference_row = row_names.index(reference_name)
                indexed = []
                for index, item in entry["angles"].items():
                    row = index + 1 if isinstance(index, int) else row_names.index(f"{field}({index})")
                    indexed.append((row, item))
                indexed.sort(key=lambda item: item[0])
                selected_rows = [reference_row, *[row for row, _ in indexed]]
                if not selected_rows or len(selected_rows) != len(set(selected_rows)):
                    raise Exception(__name__ + f'._direction_groups: Invalid production row selection for "{base_key}"')
                if field == "g_tensor":
                    if entry.get("geometry") == "projective" and not common_phase:
                        entry["orientation"] = 1.0 if float(table["rows"][reference_row]) >= 0.0 else -1.0
                else:
                    _, phase_column = self._angular_columns(base_key)
                    reference_phase = tools.canonicalize_phase(float(table["rows"][reference_row][phase_column]))
                    if abs(math.sin(reference_phase)) > 1.0e-8 and not entry.get("row_phase", False):
                        raise Exception(__name__ + f'._direction_groups: Non-real projective reference for "{base_key}"'
                        )
                    if entry.get("geometry") == "projective" and not common_phase:
                        entry["orientation"] = 1.0 if math.cos(reference_phase) >= 0.0 else -1.0
                ordered = [item for _, item in indexed]
            else:
                expected = set(range(max(0, row_count - 1)))
                if found != expected:
                    missing = sorted(expected - found)
                    extra = sorted(found - expected)
                    raise Exception(__name__
                        + f'._direction_groups: Incomplete projective angle set for "{base_key}" (missing={missing}, extra={extra}, rows={row_count})'
                    )
                selected_rows = list(range(row_count))
                ordered = [entry["angles"][k] for k in range(row_count - 1)]
            entry["ordered_angle_keys"] = [item[0] for item in ordered]
            row_keys = [self._projective_vector_row_keys(base_key, row) for row in selected_rows]
            entry["mag_keys"] = [keys[0] for keys in row_keys]
            entry["phase_keys"] = [keys[1] for keys in row_keys if keys[1]]
            entry["real"] = field == "g_tensor"
            entry["complex_keys"] = [key for row in selected_rows
                for key in (tools.complex_re_key(f"{base_key}[{row}]"), tools.complex_im_key(f"{base_key}[{row}]"))]
            conflicts = (("direct magnitude inputs", entry["mag_keys"]), ("Cartesian inputs", entry["complex_keys"]),
                ("direct phases", entry["phase_keys"]))
            for label, keys in conflicts:
                if not param.keys().isdisjoint(keys):
                    raise Exception(__name__ + f'._direction_groups: Cannot mix angles and {label} for "{base_key}"')
            angles = [item[1] for item in ordered]
            direction = (tools.hemisphere_vector_from_angles(angles)
                if entry.get("geometry") == "projective" and common_phase
                else tools.projective_vector_from_angles(angles) if entry.get("geometry") == "projective"
                else tools.spherical_vector_from_angles(angles))
            completion_weights = [value for value in table.get("completion_weights", [1.0] * row_count)]
            selected_weights = [completion_weights[row] for row in selected_rows]
            completed_norm = array.sqrt(
                sum(weight * value * value for weight, value in zip(selected_weights, direction, strict=True)))
            if completed_norm <= tools.COMPLEX_EPS:
                raise Exception(__name__ + f'._direction_groups: Zero completed norm for "{base_key}"')
            if table.get("completion_weights") is not None:
                direction = [value / completed_norm for value in direction]
                entry["completion_weights"] = selected_weights
            entry["decoded_vector"] = [entry.get("orientation", 1.0) * value for value in direction]
            # Preserve signed rows and Cartesian couplings through zero for amplitude autograd
            if "cartesian" in entry:
                residue = entry["cartesian"]["re"][1] + 1j * entry["cartesian"]["im"][1]
            else:
                phase = param.get(entry.get("coupled_phase_key"), 0.0)
                phase = array.asarray(phase, dtype=float)
                xp = array.namespace(phase)
                residue = entry.get("norm", 1.0) * (xp.cos(phase) + 1j * xp.sin(phase))
            entry["amplitudes"] = [residue * value for value in entry["decoded_vector"]]

        return angle_groups

    # Decode marked raw and Cayley phase inputs into card-space phases
    def _prepare_scalar_phase_params(self, param: dict) -> tuple[dict, list[_Coordinates]]:
        groups = self._collect_coordinate_groups(param, "phase")
        transformed = _decode_coordinates(param, groups)
        for key, value in transformed.items():
            if self._is_scalar_phase_target(key):
                transformed[key] = tools.canonicalize_phase(value)
        return transformed, groups

    # Decode Cartesian channel inputs into card-space magnitude and phase
    def _prepare_channel_coupling_params(self, param: dict) -> tuple[dict, list[_Coordinates]]:
        groups = self._collect_coordinate_groups(param, "complex")
        return _decode_coordinates(param, groups), groups

    # Decode hyperspherical angles into signed magnitude and phase rows
    def _direction_params(self, param: dict, *, path: str) -> tuple[dict, dict]:
        transformed = dict(param)
        angle_groups = self._direction_groups(param, path=path)

        for entry in angle_groups.values():
            if "norm_key" in entry:
                transformed.pop(entry["norm_key"], None)
            for key in entry.get("cartesian_keys", []):
                transformed.pop(key, None)
            for angle_key in entry["ordered_angle_keys"]:
                transformed.pop(angle_key, None)
            if "residue_phase" in entry and "coupled_phase_key" in entry:
                transformed[entry["coupled_phase_key"]] = entry["residue_phase"]
            if entry["real"]:
                for real_key, coefficient in zip(entry["mag_keys"], entry["decoded_vector"], strict=False):
                    transformed[real_key] = entry.get("norm", 1.0) * coefficient
            else:
                for mag_key, phase_key, coefficient in zip(
                    entry["mag_keys"], entry["phase_keys"], entry["decoded_vector"], strict=False):
                    transformed[mag_key] = entry.get("norm", 1.0) * abs(coefficient)
                    phase = tools.real_phase(coefficient)
                    if entry.get("row_phase", False):
                        phase += entry["residue_phase"]
                    transformed[phase_key] = tools.canonicalize_phase(phase) if transformed[mag_key] > 0.0 else 0.0

        return transformed, angle_groups

    # Encode original real coupling rows as hyperspherical angles
    def _direction_defaults(self, orig_values: dict, angle_groups: dict) -> dict:
        encoded_values = copy.deepcopy(orig_values)

        for base_key, entry in angle_groups.items():
            coefficients = []
            if entry["real"]:
                for real_key in entry["mag_keys"]:
                    if real_key not in encoded_values:
                        raise Exception(__name__ + f'._direction_defaults: Missing original value for "{base_key}"')
                    coefficients.append(float(encoded_values.pop(real_key)))
            else:
                complex_coefficients = []
                for mag_key, phase_key in zip(entry["mag_keys"], entry["phase_keys"], strict=False):
                    if mag_key not in encoded_values or phase_key not in encoded_values:
                        raise Exception(__name__ + f'._direction_defaults: Missing original value for "{base_key}"')
                    magnitude = float(encoded_values.pop(mag_key))
                    phase = tools.canonicalize_phase(encoded_values.pop(phase_key))
                    complex_coefficients.append(magnitude * complex(math.cos(phase), math.sin(phase)))
                common_row_phase = 0.0
                if entry.get("row_phase", False):
                    reference = next((value for value in complex_coefficients if abs(value) > tools.COMPLEX_EPS), 0.0j)
                    common_row_phase = tools.canonicalize_phase(math.atan2(reference.imag, reference.real))
                    rotation = complex(math.cos(common_row_phase), -math.sin(common_row_phase))
                    complex_coefficients = [value * rotation for value in complex_coefficients]
                if any(abs(value.imag) > 1.0e-8 * max(1.0, abs(value)) for value in complex_coefficients):
                    raise Exception(__name__ + f'._direction_defaults: "{base_key}" is not a real residue')
                coefficients.extend(value.real for value in complex_coefficients)

            if "coupled_phase_key" in entry:
                phase_key = entry["coupled_phase_key"]
                if phase_key not in encoded_values:
                    raise Exception(__name__ + f'._direction_defaults: Missing common phase for "{base_key}"')
                sign_reference = next((value for value in coefficients if abs(value) > tools.COMPLEX_EPS), 0.0)
                if sign_reference < 0.0:
                    coefficients = [-value for value in coefficients]
                    encoded_values[phase_key] = tools.canonicalize_phase(float(encoded_values[phase_key]) + math.pi)
            elif entry.get("row_phase", False):
                sign_reference = next((value for value in coefficients if abs(value) > tools.COMPLEX_EPS), 0.0)
                if sign_reference < 0.0:
                    coefficients = [-value for value in coefficients]
                    common_row_phase = tools.canonicalize_phase(common_row_phase + math.pi)

            angles = (tools.projective_angles_from_vector(coefficients) if entry.get("geometry") == "projective"
                else tools.spherical_angles_from_vector(coefficients))
            weights = entry.get("completion_weights", [1.0] * len(coefficients))
            norm = math.sqrt(sum(weight * value * value for weight, value in zip(weights, coefficients, strict=True)))
            if "cartesian_keys" in entry:
                phase = (common_row_phase if entry.get("row_phase", False)
                    else float(encoded_values.pop(entry["coupled_phase_key"])))
                re_value, im_value = tools.encode_polar(norm, phase)
                encoded_values[entry["cartesian_keys"][0]] = re_value
                encoded_values[entry["cartesian_keys"][1]] = im_value
            elif "norm_key" in entry:
                encoded_values[entry["norm_key"]] = norm
            for angle_key, angle in zip(entry["ordered_angle_keys"], angles, strict=False):
                encoded_values[angle_key] = angle

        return encoded_values

    # Read one decoded normalized vector exactly as stored in a generated tune card
    def _projective_card_table(self, *, base_key: str, path: str) -> dict:
        parts = tools.substring_split(s=base_key, delimiters=["|", ":"])
        supported = self._is_direction_target(base_key)
        if not supported:
            raise Exception(__name__ + f'._projective_card_table: Unsupported projective base "{base_key}"')
        table = self._coupling_table(base_key, path)
        if parts[0] in {"CON_MP", "CON_XP", "CON_GP"}:
            labels = [angular_row_name(parts[3], row).partition("(")[2][:-1] for row in table["rows"]]
            if (weights := tools.continuum_completion_weights(base_key, labels)) is not None:
                table["completion_weights"] = weights
        return table

    # Decode optimizer coordinates in the same order for steering cards and summaries
    def _prepare_card_params(self, param: dict, *, path: str) -> tuple[dict, dict, tuple]:
        prepared, physical = self._physical_params(param, path=path)
        prepared, eikonal = self._prepare_eikonal_canonical_params(prepared)
        prepared, phase = self._prepare_scalar_phase_params(prepared)
        prepared, coupling = self._prepare_channel_coupling_params(prepared)
        prepared, direction = self._direction_params(prepared, path=path)
        return prepared, physical, (eikonal, phase, coupling, direction)

    # Decode a rho population into the parity related diagonal entries
    def _rho_diagonal(self, param: dict, base: str):
        entries = sorted((tools.parse_simplex_theta_key(key)[1], value) for key, value in param.items()
                         if tools.is_simplex_theta_key(key) and tools.parse_simplex_theta_key(key)[0] == base)
        probabilities = tools.simplex_probabilities_from_angles([value for _, value in entries])
        side = [0.5 * value for value in probabilities[1:]]
        return [*reversed(side), probabilities[0], *side]

    # Compute named physical values through the same decoders used to write steering cards
    def _physical_card_parameters(self, prepared: dict, physical_keys: dict, *, path: str) -> dict:
        parameters = self._restore_physical_card_keys({key: value for key, value in prepared.items()
             if not self._is_ajz_angle_key(key) and not self._is_ajz_phase_key(key) and not tools.is_simplex_theta_key(key)},
            physical_keys)
        spin_res = {self._resonance_name_from_key(key) for key in prepared if self._is_ajz_angle_key(key)}
        for res in sorted(spin_res):
            data = self._parameter_card(path, f"RES/{res}.json")
            general = self._parameter_card(path, "GENERAL.json")
            param = {**prepared, "REGGE|MP_FRAME": prepared.get("REGGE|MP_FRAME", general["PARAM_REGGE"]["MP_FRAME"])}
            self._update_ajzp_param(res, [key for key in param if key.startswith(f"RES|{res}:")], param, data, {})
            for m, magnitude, phase in self._resonance_block(json_data=data, model="MP")["polarization"]["a_Jz"]:
                base = f"RES|{res}:MP:polarization.a_Jz({m})"
                parameters[tools.magnitude_key(base)] = magnitude
                parameters[tools.raw_phase_key(base)] = phase
        simplex = {tools.parse_simplex_theta_key(key)[0] for key in prepared if tools.is_simplex_theta_key(key)}
        for base in sorted(simplex):
            diagonal = self._rho_diagonal(prepared, base)
            j = len(diagonal) // 2
            for index, value in enumerate(diagonal):
                parameters[f"{base.rsplit('.', 1)[0]}.rho({index - j},{index - j})@MAG"] = value
        return parameters

    # Bind the named steering decoder for generic autograd covariance and impact calculations
    def parameter_transform(self, names, reference, *, cdir: str, metadata: dict):
        from core.stats.transform import ParameterTransform

        tune = metadata.get("mc_steer", {}).get("tune_default")
        if not tune:
            raise ValueError("Physical parameter plots require the saved mc_steer.tune_default")
        path = _resolve_tune_dir(cdir=cdir, tunename=tune)
        if not path.is_dir():
            raise FileNotFoundError(f"Physical parameter source tune is missing: {path}")
        cards = {str(filename.relative_to(path)): load_json_file(filename, loader=json5.load)
                 for filename in path.rglob("*.json")}

        # Preserve the complete tensor graph through the standard parameter decoders
        def decode(config):
            prepared, keys, _ = self._prepare_card_params(config, path=cards)
            return self._physical_card_parameters(prepared, keys, path=cards)

        physical = decode(dict(zip(names, reference, strict=True)))
        periods = {key: 2.0 * math.pi for key in physical
                   if tools.is_raw_phase_key(key) or self._is_scalar_phase_target(key)}
        return ParameterTransform(names, decode, reference, periods)

    # Compute every optimized value after decoding into generated-card coordinates
    def _decoded_card_config(self, *, config: dict, path: str) -> dict:
        prepared, physical_keys, groups = self._prepare_card_params(config, path=path)
        projective_groups = groups[-1]

        parameters = self._physical_card_parameters(prepared, physical_keys, path=path)
        tables = {base_key: self._projective_card_table(base_key=base_key, path=path)
            for base_key in sorted(projective_groups)}
        for base_key, table in tables.items():
            table["representation"] = ("signed_real_projective"
                if projective_groups[base_key].get("geometry") == "projective" else "signed_real_spherical")
            if table["basis"] != "g_tensor":
                table["phase_encoding"] = "0 for nonnegative coefficients and -pi for negative coefficients"

        ajzp_bases = set()
        for key in config:
            if tools.is_ajzp_angle_key(key):
                ajzp_bases.add(tools.parse_ajzp_angle_key(key)[0])
            elif self._is_ajz_angle_key(key):
                ajzp_bases.add(tools.parse_spherical_angle_key(key)[0])
        for base_key in sorted(ajzp_bases):
            resonance = self._resonance_name_from_key(base_key)
            json_data = load_json_file(os.path.join(path, "RES", f"{resonance}.json"), loader=json5.load)
            block = self._resonance_block(json_data=json_data, model="MP")
            spin_x2 = int(json_data["PARAM_RES"]["spinX2"])
            j = spin_x2 // 2
            sector_labels = mp_spin_sectors(j, config.get("REGGE|MP_FRAME", load_json_file(
                os.path.join(path, "GENERAL.json"), loader=json5.load)["PARAM_REGGE"]["MP_FRAME"]))
            tables[base_key] = {"basis": "a_Jz", "representation": ("complex_spin"
                    if any(key.startswith(f"{base_key.rsplit('.', 1)[0]}.phi_Jz[") for key in config) else
                    "signed_real_projective" if any(
                        tools.is_ajzp_angle_key(key) and tools.parse_ajzp_angle_key(key)[0] == base_key
                        for key in config) else "signed_real_spherical"), "reference_sector": sector_labels[0],
                "rows": copy.deepcopy(block["polarization"]["a_Jz"])}

        simplex_bases = {tools.parse_simplex_theta_key(key)[0] for key in config if tools.is_simplex_theta_key(key)}
        for base_key in sorted(simplex_bases):
            resonance = self._resonance_name_from_key(base_key)
            json_data = load_json_file(os.path.join(path, "RES", f"{resonance}.json"), loader=json5.load)
            block = self._resonance_block(json_data=json_data, model="MP")
            polarization = block["polarization"]
            tables[base_key] = {"basis": "rho", "rho_mag": copy.deepcopy(polarization["rho_mag"]),
                "rho_phase": copy.deepcopy(polarization["rho_phase"])}

        return {"schema_version": 1, "parameters": parameters, "tables": tables}

    def update_general_param(self, path: str, param: dict) -> dict:
        filename = os.path.join(path, "GENERAL.json")
        json_data = load_json_file(filename, loader=json5.load)

        orig_values = {}
        for par in param:
            cls = par.partition("|")[0]
            if cls in {"RES", "DECAY", "CON_MP", "CON_XP", "CON_GP", "CON_TP"}:
                continue

            if ":" not in par:
                parname = par.partition("|")[2]
                orig_values[par] = self.update_block(
                    parname=parname, blockname=f"PARAM_{cls}", new_value=param[par], json_data=json_data)
            else:
                subname = tools.find_str_between(s=par, start="|", end=":")
                parname = par.partition(":")[2]
                block = json_data[f"PARAM_{cls}"]
                old_value = self.update_block(parname=parname, blockname=subname, new_value=param[par], json_data=block)
                if par not in orig_values:
                    orig_values[par] = old_value

        write_json_file(filename, json_data, indent=4, reference_root=path)
        return orig_values

    # Update resonance fields within their selected production model
    def update_resonance_param(self, path: str, param: dict) -> dict:
        general = load_json_file(os.path.join(path, "GENERAL.json"), loader=json5.load)
        resonances = set()
        for par in param:
            if "RES|" in par:
                resonances.add(self._resonance_name_from_key(par))

        orig_values = {}
        for res in resonances:
            filename = os.path.join(path, "RES", f"{res}.json")
            json_data = load_json_file(filename, loader=json5.load)

            res_params = [j for j in param if j.startswith(f"RES|{res}:")]
            for par in res_params:
                if self._is_ajz_angle_key(par) or self._is_ajz_phase_key(par):
                    continue
                if tools.is_simplex_theta_key(par):
                    continue
                model, parname = self._split_resonance_param(par)
                if model == "TP" and ":" in parname:
                    channel, parname = parname.split(":", 1)
                    target = json_data["PARAM_RES"]["MODELS"]["TP"][channel]
                else:
                    target = self._resonance_target(json_data=json_data, model=model, parname=parname)
                if parname == "FF_transfer.LambdaInv2":
                    orig_values[par] = self._update_transfer_inv2(target, param[par])
                else:
                    orig_values[par] = self.update_target(parname=parname, new_value=param[par], target=target)

            if any(self._split_resonance_param(par)[0] == "MP" for par in res_params):
                self._canonicalize_mp_polarization(json_data=json_data)

            self._update_ajzp_param(res=res, res_params=res_params,
                param={**param, "REGGE|MP_FRAME": param.get("REGGE|MP_FRAME", general["PARAM_REGGE"]["MP_FRAME"])},
                json_data=json_data, orig_values=orig_values)

            simplex_entries = []
            expected_simplex_base = f"RES|{res}:MP:polarization.rho_population"
            for par in res_params:
                if not tools.is_simplex_theta_key(par):
                    continue
                base_key, angle_index = tools.parse_simplex_theta_key(par)
                if base_key != expected_simplex_base:
                    raise Exception(__name__ + f'update_resonance_param: Unsupported simplex target "{base_key}"')
                simplex_entries.append((angle_index, par))
            if simplex_entries:
                simplex_entries.sort()
                if [index for index, _ in simplex_entries] != list(range(len(simplex_entries))):
                    raise Exception(__name__ + f'update_resonance_param: Incomplete rho simplex for "{res}"')
                spin_x2 = int(json_data["PARAM_RES"]["spinX2"])
                if spin_x2 < 2 or spin_x2 % 2 != 0:
                    raise Exception(
                        __name__ + f'update_resonance_param: rho simplex requires positive integer spin for "{res}"')
                j = spin_x2 // 2
                if len(simplex_entries) != j:
                    raise Exception(__name__ + f'update_resonance_param: rho simplex size mismatch for "{res}"')
                res_block = self._resonance_block(json_data=json_data, model="MP")
                polarization = res_block["polarization"]
                rho_mag = polarization["rho_mag"]
                default_probabilities = [float(rho_mag[j][j])]
                default_probabilities.extend(2.0 * float(rho_mag[j - m][j - m]) for m in range(1, j + 1))
                default_angles = tools.simplex_angles_from_probabilities(default_probabilities)
                for (_, key), angle in zip(simplex_entries, default_angles, strict=False):
                    orig_values[key] = angle
                size = 2 * j + 1
                diagonal = self._rho_diagonal(param, expected_simplex_base)
                polarization["rho_mag"] = np.diag(diagonal).tolist()
                polarization["rho_phase"] = np.zeros((size, size), dtype=float).tolist()

            write_json_file(filename, json_data, indent=4, reference_root=path)

        return orig_values

    # Update one particle-indexed decay or production channel card
    def _update_channel_card(
        self, *, path: str, param: dict, parameter_class: str, filename: str, select_target: Callable[[dict, str], dict]
    ) -> dict:
        prefix = f"{parameter_class}|"
        particle_ids = {
            tools.substring_split(s=key, delimiters=["|", ":"])[1] for key in param if key.startswith(prefix)}
        card_path = os.path.join(path, filename)
        json_data = load_json_file(card_path, loader=json5.load)
        orig_values = {}
        for particle_id in particle_ids:
            block = copy.deepcopy(json_data[particle_id])
            relative_shape_targets = {}
            selected = [item for item in param if item.startswith(f"{prefix}{particle_id}:")]
            if parameter_class.startswith("CON_"):
                for key in selected:
                    parts = tools.substring_split(s=key, delimiters=["|", ":"])
                    if parts[-1] != "FF_offshell.type":
                        continue
                    pair = parts[-2]
                    target = block[pair]
                    current = target["FF_offshell"]
                    term_indices = [int(index) for item in selected
                        if item.startswith(f"{prefix}{particle_id}:{pair}:FF_offshell.terms[")
                        for index in re.findall(r"terms\[(\d+)\]", item)]
                    kernels = max(term_indices, default=0) + 1
                    rebuild = current.get("type") != param[key] or (
                        param[key] == "gkernel" and len(current.get("terms", [])) != kernels)
                    if rebuild:
                        orig_values[key] = current.get("type")
                        target["FF_offshell"] = self._continuum_ff_template(str(param[key]), kernels)
                        if current.get("type") in {"exp", "logexp"} and param[key] == "orear":
                            # Match the logarithmic slope at the meson pole
                            target["FF_offshell"]["b"] = 2 * target["FF_offshell"]["a"] * current["b"]

            for key in selected:
                parts = tools.substring_split(s=key, delimiters=["|", ":"])
                if parameter_class.startswith("CON_") and parts[-1].startswith(
                    ("FF_transfer.", "FF_offshell.", "pveto.", "reggeize.")):
                    pair = parts[-2]
                    if pair not in block:
                        raise Exception(__name__ + f"._update_channel_card: Unknown continuum pair {pair}")
                    target = block[pair]
                    if parts[-1] == "FF_transfer.LambdaInv2":
                        orig_values[key] = self._update_transfer_inv2(target, param[key])
                        continue
                else:
                    target = select_target(block, parts[-2])
                if parameter_class == "DECAY" and target.get("basis") in {"alpha_ls", "helicity"}:
                    relative_shape_targets[id(target)] = target
                old_value = self.update_target(parname=parts[-1], new_value=param[key], target=target)
                if key not in orig_values:
                    orig_values[key] = old_value
            for target in relative_shape_targets.values():
                self._canonicalize_relative_shape(target)
            json_data[particle_id] = block
        write_json_file(card_path, json_data, indent=4, reference_root=path)
        return orig_values

    # Apply decay-channel parameter updates to DECAYS.json
    def update_branching_param(self, path: str, param: dict) -> dict:
        return self._update_channel_card(path=path, param=param, parameter_class="DECAY", filename="DECAYS.json",
            select_target=lambda block, channel: block[channel])

    # Apply central continuum parameter updates to the four model cards
    def update_production_param(self, path: str, param: dict) -> dict:
        values = {}
        for model in ("MP", "XP", "GP"):
            parameter_class = f"CON_{model}"
            values.update(self._update_channel_card(path=path, param=param, parameter_class=parameter_class,
                    filename=f"{parameter_class}.json", select_target=lambda block, channel: (
                        block if channel == "[*]" else self.production_sector_block(mother_block=block, channel=channel)
                    )))
        values.update(self._update_channel_card(path=path, param=param, parameter_class="CON_TP",
                filename="CON_TP.json", select_target=lambda block, channel: block[channel]))
        return values

    # Evaluate the shared histogram objective without detaching amplitude derivatives
    def trial_costs(self, results: dict, param: dict) -> dict:
        use_full_covariance = param.get("data_covariance_mode") == "full"
        covariance_payload = None
        if use_full_covariance:
            if self.data_covariance_payload is None:
                npz_path, json_path = self.data_covariance_paths(run_name=param["run_name"], cdir=param["cdir"])
                self.load_data_covariance(npz_path=npz_path, json_path=json_path)
            covariance_payload = self.data_covariance_payload

        cost_results = results
        if (param["mc_steer"].get("ampfit") or {}).get("mc_stat") == "source":
            from core.tune.drivers.graniitti.ampfit.amplitude import source_mc

            cost_results = source_mc(results)
        return icecost.evaluate_cost_bundle(results=cost_results, selected_cost=str(param["cost"]),
            cost_rho=str(param["cost_rho"]), cost_avg=str(param["cost_avg"]), covariance_payload=covariance_payload,
            rngseed=int(param.get("rngseed", 0)), wasserstein_cache=self.wasserstein_replica_cache)

    # Evaluate a worker trial including its steering cards and output bundle
    def evaluate_trial_outputs(self, *, config: dict, param: dict, trial_id: str, tunename: str) -> dict:
        """Evaluate one GRANIITTI trial and return a renderable output bundle."""

        param_space = {**config, **param["aux_param_space"]}
        card_config = None
        stage = "steering_card"
        deadline = time.monotonic() + float(param["max_t"])
        event_dir = None

        try:
            self.create_steering_card(param_space=param_space, tunename=tunename, cdir=param["cdir"],
                                      tune_default=param["mc_steer"]["tune_default"])
            card_config = self._decoded_card_config(
                config=config, path=str(_resolve_tune_dir(cdir=param["cdir"], tunename=tunename)))
            if param.get("save_events"):
                scratch = pathlib.Path(param["cdir"], "tmp")
                ensure_dir(scratch)
                event_dir = pathlib.Path(tempfile.mkdtemp(prefix="icetune-events-", dir=scratch))
                shutil.copytree(_resolve_tune_dir(cdir=param["cdir"], tunename=tunename), event_dir / "modeldata")
                shutil.copy2(pathlib.Path(param["cdir"]) / "VERSION.json", event_dir / "VERSION.json")

            stage = "simulation"
            amplitude = None
            if param.get("optimization", {}).get("optimizer") == "ampfit":
                import torch

                names = sorted(config)
                theta = torch.tensor([config[name] for name in names], dtype=torch.float64, requires_grad=True)
                amplitude = {"directory": self.amplitude_directory(cdir=param["cdir"], run_name=param["run_name"]),
                             "prepare": False, "parameters": dict(zip(names, theta, strict=True)) | param["aux_param_space"],
                             "controls": param["mc_steer"]["ampfit"], "workers": param.get("ampfit_workers", {})}
            mc = self.compute(tunename=tunename, datacards=param["datacards"], mc_steer=param["mc_steer"],
                cdir=param["cdir"], chunksize=param["chunksize"], processes=param["processes"], max_t=param["max_t"],
                deadline=deadline, rngseed=int(param.get("rngseed", 0)), amplitude=amplitude, event_dir=event_dir)
            stage = "cost_evaluation"
            _remaining_time(deadline)
            results = {"mc": mc, "data": self.data, "obs": self.obs, "datasets": self.datasets}

            costs = self.trial_costs(results, param)
            _remaining_time(deadline)
            gradient = None
            if amplitude is not None:
                gradient = {}
                if names:
                    objective = costs["metrics"][param["cost"]]
                    if not np.isfinite(array.to_numpy(objective)):
                        raise ValueError("Amplitude objective is nonfinite at the proposed parameters")
                    derivative, = torch.autograd.grad(objective, theta)
                    if not torch.all(torch.isfinite(derivative)):
                        raise ValueError("Amplitude objective has a nonfinite autograd gradient")
                    gradient = dict(zip(names, derivative.detach().tolist(), strict=True))
                costs, results = hist.numpy_output((costs, results))

            return {"event_samples": None if event_dir is None else str(event_dir),
                "card_config": copy.deepcopy(card_config), "config": copy.deepcopy(config),
                "cost_arr": copy.deepcopy(costs["cost_arr"]), "error": None, "gradient": gradient,
                "likelihood": copy.deepcopy(costs["likelihood"]), "metrics": copy.deepcopy(costs["metrics"]),
                "ndf_arr": copy.deepcopy(costs["ndf_arr"]), "results": results, "trial_id": trial_id,
                "tunename": tunename, "valid_arr": copy.deepcopy(costs["valid_arr"]),
                "weight_arr": copy.deepcopy(costs["weight_arr"])}
        except Exception as exc:
            if event_dir is not None:
                shutil.rmtree(event_dir, ignore_errors=True)
            logger.exception(".evaluate_trial_outputs: trial_id=%s failed during %s step", trial_id, stage)
            return _penalty_trial_outputs(config=config, trial_id=trial_id, tunename=tunename, card_config=card_config,
                error=str(exc), error_type=exc.__class__.__name__,
                failure=exc.failure if isinstance(exc, GraniittiCommandError) else None, stage=stage,
                traceback_text=traceback.format_exc())

    # Render one GRANIITTI trial figure tree into a caller-provided directory
    def render_trial_figures_to_dir(
        self, *, outputs: dict, param: dict, summary_payload: dict, output_dir: str, summary_file: str | None = None
    ) -> dict:
        if outputs.get("results") is None or outputs.get("valid_arr") is None:
            raise RuntimeError(f'Cannot render figures for non-renderable trial "{outputs.get("trial_id")}"')

        return loss.visualize_losses(results=outputs["results"], run_name=param["run_name"], cdir=param["cdir"],
            tunename=outputs["tunename"], valid_arr=outputs["valid_arr"]["chi2"], likelihood=outputs.get("likelihood"),
            summary=summary_payload, output_dir=output_dir, summary_file=summary_file)

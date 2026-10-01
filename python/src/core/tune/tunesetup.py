# Load and save explicit JSON tuning definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path
from types import SimpleNamespace

import pyjson5
from ray import tune

from core.io.files import ensure_dir
from core.io.serialize import load_json_file
from core.tune.drivers.registry import get_driver_type
from core.tune.parameters.space import normalize_param_space


# Resolve a JSON tuning definition with an optional JSON Pointer fragment
def resolve_tunesetup(*, cdir, simdriver, name):
    root = Path(cdir).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"icetune working directory does not exist: {root}")
    source = Path(str(name).partition("#")[0]).expanduser()
    source = (source if source.is_absolute() else root / source).resolve()
    if source.suffix != ".json":
        raise ValueError(f"Tuning setup must be a JSON file: {name}")
    if not source.is_file():
        raise FileNotFoundError(f"Tuning setup for {simdriver} does not exist: {source}")
    return source


# Construct parameter domains from explicit bounds using the shared representation
def parameter_space(parameters):
    return {key: tune.randint(spec["lower"], spec["upper"] + 1) if spec["type"] == "int"
            else tune.uniform(spec["lower"], spec["upper"])
            for key, spec in normalize_param_space(parameters).items()}


# Build a campaign selection or restore its complete resolved definition
def load_tunesetup(*, cdir, simdriver, name, tune_default=None):
    path = resolve_tunesetup(cdir=cdir, simdriver=simdriver, name=name)
    fragment = str(name).partition("#")[2]
    source = str(path) + ("#" + fragment if fragment else "")
    payload = load_json_file(source, loader=pyjson5.load)
    resolved = payload.get("schema_version") == 1
    if str(payload.get("simdriver" if resolved else "driver", "")).upper() != simdriver.upper():
        raise ValueError(f"Tuning setup driver disagrees with {simdriver}: {source}")
    driver = get_driver_type(simdriver)
    baseline = payload["tune_default"] if resolved else payload.get("fit", {}).get("tune_default", driver.DEFAULT_TUNE)
    if resolved and tune_default is not None and baseline != tune_default:
        raise ValueError("The saved tuning setup uses a different baseline model tune")
    baseline = baseline if tune_default is None else tune_default
    if resolved:
        cards, parameters, fixed = payload["datacards"], parameter_space(payload["param_space"]), payload["aux_param_space"]
    else:
        cards, parameters, fixed = driver.build_tunesetup(config=payload, cdir=cdir, tune_default=baseline)
        explicit = parameter_space(payload.get("parameters", {}))
        constants = payload.get("fixed", {})
        if explicit.keys() & constants.keys():
            raise ValueError("A parameter cannot be both varied and fixed")
        parameters.update(explicit)
        for key in explicit:
            fixed.pop(key, None)
        fixed.update(constants)
        for key in constants:
            parameters.pop(key, None)
    study = SimpleNamespace(**payload)
    study.__file__, study.source, study.tune_default = str(path), payload.get("source", source), baseline
    study.datacards, study.param_space, study.aux_param_space = cards, parameters, fixed
    return study


# Save the complete sampled and fixed tuning definition before any simulation
def save_tunesetup(*, tunesetup, path: str | Path, simdriver: str, tune_default: str) -> None:
    payload = {
        "schema_version": 1,
        "source": getattr(tunesetup, "source", tunesetup.__file__),
        "simdriver": simdriver.upper(),
        "tune_default": tune_default,
        "datacards": tunesetup.datacards,
        "param_space": normalize_param_space(tunesetup.param_space),
        "aux_param_space": tunesetup.aux_param_space,
        "optimizer": getattr(tunesetup, "optimizer", {}),
    }
    if hasattr(tunesetup, "bank_samples"):
        payload["bank_samples"] = tunesetup.bank_samples
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    path = Path(path)
    if path.is_file():
        if json.loads(path.read_text(encoding="utf-8")) != json.loads(encoded):
            raise ValueError(f"Run already has a different tuning definition: {path}")
        return
    ensure_dir(path.parent)
    path.write_text(encoded, encoding="utf-8")

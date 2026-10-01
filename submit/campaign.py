#!/usr/bin/env python3
# Read and resolve icetune campaigns
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import pathlib
import re

import yaml

from submit import CAMPAIGN_DIR

SUPPORTED_SIMDRIVERS = frozenset({"GRANIITTI", "PANDORA"})
SUPPORTED_ALGORITHMS = {"ampfit", "ax", "basic", "bayesopt", "hebo", "hyperopt", "icebo", "lbfgs",
        "optuna", "scikit", }
COMMON_REQUIRED_ENV = frozenset(["ALGORITHM", "COST", "COST_AVG", "CPU_PER_TRIAL", "DATA_COVARIANCE_MODE", "MAX_T", "MAX_CONCURRENT", "NUM_TRIALS", "FULL_OUTPUT", "PLOT", "PROPOSAL_BATCH_FRACTION", "PROPOSAL_BATCH_SIZE", "RAND_TRIALS", "RNGSEED", "RUN_NAME", "SIMDRIVER", "TUNESETUP", "TUNE_DEFAULT"])

RAY_PORTS = tuple(["RAY_CLIENT_SERVER_PORT", "RAY_DASHBOARD_AGENT_GRPC_PORT", "RAY_DASHBOARD_AGENT_LISTEN_PORT", "RAY_MAX_WORKER_PORT", "RAY_METRICS_EXPORT_PORT", "RAY_MIN_WORKER_PORT", "RAY_NODE_MANAGER_PORT", "RAY_OBJECT_MANAGER_PORT", "RAY_RUNTIME_ENV_AGENT_PORT"])

RAY_MEMORY = tuple(["RAY_INIT_REQUEST_DISK", "RAY_INIT_REQUEST_MEMORY", "RAY_REQUEST_DISK", "RAY_HEAD_REQUEST_MEMORY", "RAY_REQUEST_MEMORY", "RAY_STEER_REQUEST_MEMORY"])

RAY_POSITIVE = tuple(["RAY_WORKER_CPU", "RAY_HEAD_CPU", "RAY_HISTORY_INTERVAL_S", "RAY_INIT_MAX_RUNTIME_S", "RAY_INIT_CPU", "RAY_INIT_JOBS", "RAY_LXPLUS_MAX_STARTING", "RAY_LXPLUS_NODE_START_TIMEOUT_S", "RAY_LXPLUS_POLL_INTERVAL_S", "RAY_LXPLUS_SHUTDOWN_TIMEOUT_S", "RAY_MAX_RUNTIME_S", "RAY_WORKERS", "RAY_PORT", "RAY_STARTUP_TIMEOUT_S", "RAY_WORKER_ADMISSION_TIMEOUT_S", "RAY_WORKER_POLL_INTERVAL_S", "RAY_WORKER_FAILURE_LIMIT", "RAY_WORKER_RECOVERY_TIMEOUT_S", "RAY_RESTORE_RETRY_DELAY_S", "RAY_RESOURCE_STATUS_INTERVAL_S", "RAY_STATUS_INTERVAL_S", "RAY_RAYLET_START_TIMEOUT_S", "RAY_CONDOR_COMMAND_TIMEOUT_S", "RAY_LXPLUS_DEATH_TIMEOUT_S", "RAY_WORKER_SHUTDOWN_TIMEOUT_S"])

RAY_NONNEGATIVE = tuple(["GPU_PER_TRIAL", "RAY_WORKER_GPU", "RAY_HEAD_GPU", "RAY_INIT_GPU", "RAY_LXPLUS_JOIN_RETRY_LIMIT", "RAY_WORKER_REPLACEMENT_LIMIT", "RAY_RESTORE_RETRY_LIMIT", "RAY_TRIAL_RETRY_LIMIT"])
RAY_REQUIRED_ENV = frozenset((*RAY_POSITIVE, *RAY_NONNEGATIVE, *RAY_PORTS, *RAY_MEMORY, 'RAY_LXPLUS_START_INTERVAL_S'))

CAMPAIGN_FIELDS = {"optimizer": {name.lower(): name for name in ("ALGORITHM", "HEBO_CONFIG", "ICEBO_CONFIG",
            "AMPFIT_CONFIG", "LBFGS_CONFIG", "NUM_TRIALS", "RAND_TRIALS", "RNGSEED", "PROPOSAL_BATCH_SIZE",
            "PROPOSAL_BATCH_FRACTION", )}, "fit": {name.lower(): name
        for name in ("COST", "COST_AVG", "DATA_COVARIANCE_MODE", "PRECOMPUTE", "PLOT", "FULL_OUTPUT", "SAVE_EVENTS", "TUNE_DEFAULT", "AMPFIT_REUSE")
    }, "ray": {name.removeprefix("RAY_").lower(): name
        for name in (*RAY_REQUIRED_ENV, "CPU_PER_TRIAL", "MAX_T", "MAX_CONCURRENT")}, "runtime": {
        "shared_output_dir": "SHARED_OUTPUT_DIR", "pfa_dir": "PANDORA_PFA_DIR", "key4hep_setup": "PANDORA_KEY4HEP_SETUP",
        "key4hep_release": "PANDORA_KEY4HEP_RELEASE", "lcg_view": "ICETUNE_LCG_VIEW",
        "conda_env": "ICETUNE_CONDA_ENV", "env_setup": "ICETUNE_ENV_SETUP", }, }


# Compute one scalar YAML value in the launcher environment representation
def environment_value(value: object) -> str:
    if isinstance(value, bool):
        return str(int(value))
    if value is None or isinstance(value, (dict, list, tuple)):
        raise ValueError(f"Campaign environment values must be scalar, got {value!r}")
    rendered = str(value)
    if any(character in rendered for character in ('"', "\n", "\r")) or any(
        character.isspace() for character in rendered):
        raise ValueError(f"Campaign environment value cannot contain quotes or whitespace: {value!r}")
    return rendered


# Compute the selected simulator driver name
def simdriver_name(environment: dict[str, str]) -> str:
    name = str(environment["SIMDRIVER"]).upper()
    if name not in SUPPORTED_SIMDRIVERS:
        raise ValueError(f"Unsupported SIMDRIVER={environment['SIMDRIVER']}")
    return name


# Compute the shared output directory setting for the selected simulator driver
def shared_output_dir_key(environment: dict[str, str]) -> str:
    return f"{simdriver_name(environment)}_SHARED_OUTPUT_DIR"


# Compute the shared output directory used by the selected simulator driver
def shared_output_dir(environment: dict[str, str]) -> pathlib.Path:
    key = shared_output_dir_key(environment)
    if not environment.get(key):
        raise ValueError(f"{simdriver_name(environment)} Ray campaign requires {key}")
    return pathlib.Path(environment[key])


# Parse repeatable section.field=value campaign overrides and reject duplicates
def parse_overrides(items: list[str] | None) -> dict[str, str]:
    overrides = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"Invalid --set value {item!r}, expected section.field=value")
        key, value = item.split("=", 1)
        section, separator, field = key.partition(".")
        if not separator or field not in CAMPAIGN_FIELDS.get(section, {}):
            raise ValueError(f"Unknown campaign field {key!r}, expected section.field")
        key = CAMPAIGN_FIELDS[section][field]
        if key in overrides:
            raise ValueError(f"Duplicate --set override for {key}")
        overrides[key] = environment_value(value)
    return overrides


# Map one campaign section to launcher environment values
def section_environment(mapping: object, *, section: str, driver: str) -> dict[str, str]:
    if not isinstance(mapping, dict):
        raise ValueError(f"Campaign {section} must be a mapping")
    output = {}
    for raw_key, value in mapping.items():
        if raw_key not in CAMPAIGN_FIELDS[section]:
            raise ValueError(f"Unknown campaign field {section}.{raw_key}")
        key = CAMPAIGN_FIELDS[section][raw_key]
        if key == "SHARED_OUTPUT_DIR":
            key = f"{driver}_SHARED_OUTPUT_DIR"
        output[key] = environment_value(value)
    return output


# Validate the shared and backend default sections
def _validate_defaults(defaults: object) -> None:
    expected_defaults = set(CAMPAIGN_FIELDS)
    if not isinstance(defaults, dict) or set(defaults) != expected_defaults:
        raise ValueError(f"Campaign defaults must contain exactly {sorted(expected_defaults)}")
    for section in expected_defaults:
        section_environment(defaults[section], section=section, driver="GRANIITTI")


# Validate one named campaign entry
def _validate_campaign(name: object, campaign: object) -> None:
    if re.fullmatch(r"[a-z0-9][a-z0-9_-]*", str(name)) is None:
        raise ValueError(f"Invalid campaign name: {name!r}")
    if not isinstance(campaign, dict) or not {"driver", "tunesetup"} <= set(campaign):
        raise ValueError(f'Campaign "{name}" must select a driver and tunesetup')
    description = campaign.get("description")
    if not isinstance(description, str) or not description.strip():
        raise ValueError(f'Campaign "{name}" must contain a non-empty description')
    driver = simdriver_name({"SIMDRIVER": campaign["driver"]})
    allowed_campaign_sections = {"description", "driver", "tunesetup", "run_name", *CAMPAIGN_FIELDS}
    unknown = set(campaign) - allowed_campaign_sections
    if unknown:
        raise ValueError(f'Campaign "{name}" has unknown sections: {sorted(unknown)}')
    for section, mapping in campaign.items():
        if section not in CAMPAIGN_FIELDS:
            continue
        section_environment(mapping, section=section, driver=driver)


# Load and validate the shared campaign catalog
def load_campaign_catalog(path: pathlib.Path) -> dict:
    path = pathlib.Path(path).resolve()
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise ValueError(f"Unsupported ICETUNE campaign catalog schema: {path}")
    unknown = set(payload) - {"schema_version", "defaults", "campaigns"}
    if unknown:
        raise ValueError(f"Unknown campaign catalog keys: {sorted(unknown)}")
    _validate_defaults(payload.get("defaults"))
    campaigns = payload.get("campaigns")
    if not isinstance(campaigns, dict) or not campaigns:
        raise ValueError('Campaign catalog "campaigns" must be a non-empty mapping')
    for name, campaign in campaigns.items():
        _validate_campaign(name, campaign)
    return payload


# Compute the available simultaneous Ray trial count
def derived_max_concurrent(environment: dict[str, str]) -> int:
    workers = int(environment["RAY_WORKERS"])
    capacity = workers * (int(environment["RAY_WORKER_CPU"]) // int(environment["CPU_PER_TRIAL"]))
    gpu_per_trial = int(environment["GPU_PER_TRIAL"])
    if gpu_per_trial > 0:
        gpu_capacity = workers * (int(environment["RAY_WORKER_GPU"]) // gpu_per_trial)
        capacity = min(capacity, gpu_capacity)
    return max(1, capacity)


# Resolve defaults, campaign settings and explicit command line values
def resolve(catalog: dict, *, campaign_name: str, catalog_path=None, run_name: str | None = None,
    repo_dir: pathlib.Path | None = None,
    workers: int | None = None, worker_cpu: int | None = None, max_concurrent: int | None = None,
    num_trials: int | None = None, rand_trials: int | None = None, rngseed: int | None = None,
    environment_overrides: dict[str, str] | None = None, ) -> dict[str, str]:
    named_overrides = {"MAX_CONCURRENT": max_concurrent, "NUM_TRIALS": num_trials, "RAND_TRIALS": rand_trials,
        "RNGSEED": rngseed, "RAY_WORKER_CPU": worker_cpu, "RAY_WORKERS": workers, "RUN_NAME": run_name}
    campaigns = catalog["campaigns"]
    if campaign_name not in campaigns:
        choices = ", ".join(sorted(campaigns))
        raise ValueError(f'Unknown campaign "{campaign_name}"; available campaigns: {choices}')
    campaign = campaigns[campaign_name]
    directory = pathlib.Path(catalog_path or CAMPAIGN_DIR / "campaigns.yml").resolve().parent
    driver = simdriver_name({"SIMDRIVER": campaign["driver"]})
    environment = {"SIMDRIVER": driver, "TUNESETUP": str((directory / campaign["tunesetup"]).resolve()),
        "RUN_NAME": environment_value(campaign.get("run_name", campaign_name)), }
    for section in CAMPAIGN_FIELDS:
        values = {**catalog["defaults"][section], **campaign.get(section, {})}
        environment.update(section_environment(values, section=section, driver=driver))

    for key, value in tuple(environment.items()):
        if key.endswith("_CONFIG") and not value.startswith("package:"):
            environment[key] = str((directory / value).resolve())

    catalog_has_max_concurrent = "MAX_CONCURRENT" in environment
    named = {key: value for key, value in (named_overrides or {}).items() if value is not None}
    explicit = {f"{driver}_SHARED_OUTPUT_DIR" if key == "SHARED_OUTPUT_DIR" else key: environment_value(value)
        for key, value in (environment_overrides or {}).items()}
    overridden = set(named) | set(explicit)
    allowed = {key for fields in CAMPAIGN_FIELDS.values() for key in fields.values()}
    unknown = sorted(overridden - (allowed | set(environment) | {shared_output_dir_key(environment)}))
    if unknown:
        raise ValueError(f"Unknown campaign environment override(s): {', '.join(unknown)}")
    conflicting = sorted(set(named) & set(explicit))
    if conflicting:
        raise ValueError("Dedicated option and --set both override: " + ", ".join(conflicting))
    environment.update({key: environment_value(value) for key, value in named.items()})
    environment.update(explicit)

    # Resolve shared storage relative to the simulator submission directory
    root = pathlib.Path(repo_dir or pathlib.Path.cwd()).expanduser().resolve()
    key = shared_output_dir_key(environment)
    environment[key] = str((root / pathlib.Path(environment.get(key, ".")).expanduser()).resolve())
    if driver == "PANDORA":
        environment["PANDORA_PFA_DIR"] = str(
            (root / pathlib.Path(environment.get("PANDORA_PFA_DIR", "../PandoraPFA")).expanduser()).resolve())

    if not catalog_has_max_concurrent and "MAX_CONCURRENT" not in overridden:
        environment["MAX_CONCURRENT"] = str(derived_max_concurrent(environment))
    validate_environment(environment)
    return environment


# Require a complete environment section
def _require_keys(environment: dict[str, str], keys, *, label: str) -> None:
    missing = sorted(keys - set(environment))
    if missing:
        raise ValueError(f"Campaign is missing required {label} keys: {', '.join(missing)}")


# Require positive numeric environment values
def _require_positive(environment: dict[str, str], keys, *, number_type=int) -> None:
    for key in keys:
        if number_type(environment[key]) <= 0:
            raise ValueError(f"Campaign requires {key} > 0")


# Require non-negative integer environment values
def _require_nonnegative(environment: dict[str, str], keys) -> None:
    for key in keys:
        if int(environment[key]) < 0:
            raise ValueError(f"Campaign requires {key} >= 0")


# Validate that adaptive optimizers have separate warm-up and total budgets
def validate_trial_budget(environment: dict[str, str]) -> None:
    missing = [name for name in ("NUM_TRIALS", "RAND_TRIALS") if name not in environment]
    if missing:
        raise ValueError(f"Campaign must define {', '.join(missing)} explicitly")
    num_trials = int(environment["NUM_TRIALS"])
    rand_trials = int(environment["RAND_TRIALS"])
    if num_trials <= 0 or rand_trials < 0 or rand_trials > num_trials:
        raise ValueError(f"Invalid trial budget: NUM_TRIALS={num_trials}, RAND_TRIALS={rand_trials}")
    if environment.get("ALGORITHM", "").lower() in {"icebo", "hebo"} and rand_trials >= num_trials:
        raise ValueError(
            "Bayesian-optimization campaign requires NUM_TRIALS > RAND_TRIALS so adaptive proposals are generated")


# Validate scalar settings shared by execution backends
def validate_common_environment(environment: dict[str, str]) -> None:
    _require_keys(environment, COMMON_REQUIRED_ENV, label="environment")
    simdriver_name(environment)
    if environment["DATA_COVARIANCE_MODE"] not in {"diagonal", "full"}:
        raise ValueError("Campaign requires DATA_COVARIANCE_MODE to be diagonal or full")
    algorithm = environment["ALGORITHM"].lower()
    if environment.get("SAVE_EVENTS", "0") == "1" and (
        environment["SIMDRIVER"].upper() != "GRANIITTI" or algorithm in {"ampfit", "lbfgs"}
    ):
        raise ValueError("save_events requires GRANIITTI Ray trials generating fresh events, not ampfit or lbfgs")
    if algorithm in {"icebo", "hebo", "ampfit", "lbfgs"} and not environment.get(algorithm.upper() + "_CONFIG"):
        raise ValueError(f"{algorithm} campaigns require {algorithm.upper()}_CONFIG")
    if algorithm == "ampfit" and (
        environment["SIMDRIVER"].upper() != "GRANIITTI" or environment.get("PRECOMPUTE", "1") != "1"):
        raise ValueError("ampfit campaigns require GRANIITTI and shared initialization")
    if environment["SIMDRIVER"].upper() == "PANDORA" and environment["DATA_COVARIANCE_MODE"] != "diagonal":
        raise ValueError("Pandora campaigns require DATA_COVARIANCE_MODE=diagonal")
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", environment["RUN_NAME"]) is None:
        raise ValueError(f"Unsafe RUN_NAME={environment['RUN_NAME']!r}")
    _require_positive(environment, ("CPU_PER_TRIAL", "MAX_CONCURRENT", "MAX_T"))
    _require_nonnegative(environment, ("PROPOSAL_BATCH_SIZE", "RNGSEED"))
    fraction = float(environment["PROPOSAL_BATCH_FRACTION"])
    if not 0.0 < fraction <= 1.0:
        raise ValueError("Campaign requires 0 < PROPOSAL_BATCH_FRACTION <= 1")
    binary_keys = {"PRECOMPUTE", "PLOT", "FULL_OUTPUT", "SAVE_EVENTS", "AMPFIT_REUSE"}
    for key in binary_keys:
        if key not in environment:
            continue
        if environment[key] not in {"0", "1"}:
            raise ValueError(f"Campaign requires {key} to be 0 or 1")
    validate_trial_budget(environment)
    if environment.get("PRECOMPUTE", "1") == "0" and environment["DATA_COVARIANCE_MODE"] == "full":
        raise ValueError("Full MC covariance requires PRECOMPUTE=1")


# Validate Ray resource, port, and allocation settings
def _validate_ray_environment(environment: dict[str, str]) -> None:
    _require_positive(environment, RAY_POSITIVE, )
    _require_nonnegative(environment, RAY_NONNEGATIVE, )
    if float(environment["RAY_LXPLUS_START_INTERVAL_S"]) <= 0.0:
        raise ValueError("Campaign requires RAY_LXPLUS_START_INTERVAL_S > 0")
    if int(environment["RAY_PORT"]) > 65535:
        raise ValueError("Campaign requires RAY_PORT <= 65535")
    optional_ports = RAY_PORTS
    for key in optional_ports:
        if not 0 <= int(environment[key]) <= 65535:
            raise ValueError(f"Campaign requires 0 <= {key} <= 65535")
    minimum_worker_port = int(environment["RAY_MIN_WORKER_PORT"])
    maximum_worker_port = int(environment["RAY_MAX_WORKER_PORT"])
    if (minimum_worker_port == 0) != (maximum_worker_port == 0):
        raise ValueError("Set both RAY_MIN_WORKER_PORT and RAY_MAX_WORKER_PORT, or neither")
    if minimum_worker_port and minimum_worker_port > maximum_worker_port:
        raise ValueError("RAY_MIN_WORKER_PORT exceeds RAY_MAX_WORKER_PORT")
    workers = int(environment["RAY_WORKERS"])
    worker_cpu = int(environment["RAY_WORKER_CPU"])
    cpu_per_trial = int(environment["CPU_PER_TRIAL"])
    if cpu_per_trial > worker_cpu:
        raise ValueError("CPU_PER_TRIAL exceeds RAY_WORKER_CPU")
    capacity = workers * (worker_cpu // cpu_per_trial)
    gpu_per_trial = int(environment["GPU_PER_TRIAL"])
    worker_gpu = int(environment["RAY_WORKER_GPU"])
    if gpu_per_trial > 0:
        if gpu_per_trial > worker_gpu:
            raise ValueError("GPU_PER_TRIAL exceeds RAY_WORKER_GPU")
        capacity = min(capacity, workers * (worker_gpu // gpu_per_trial))
    if int(environment["MAX_CONCURRENT"]) > capacity:
        raise ValueError("MAX_CONCURRENT exceeds the available Ray trial slots")
    for key in RAY_MEMORY:
        if re.fullmatch(r"[1-9][0-9]*(KB|MB|GB|TB)", environment[key].upper()) is None:
            raise ValueError(f"{key} must use KB, MB, GB, or TB units")
    shared_output_dir(environment)


# Validate campaign and Ray environment settings
def validate_environment(environment: dict[str, str]) -> None:
    validate_common_environment(environment)
    algorithm = environment["ALGORITHM"].lower()
    if algorithm not in SUPPORTED_ALGORITHMS:
        raise ValueError(f"Unsupported ray optimizer {algorithm!r}")
    _require_keys(environment, RAY_REQUIRED_ENV, label="ray")
    _validate_ray_environment(environment)

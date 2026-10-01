#!/usr/bin/env python
#
# icetune: distributed blackbox optimization
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
import copy
import json
import multiprocessing
import os
import pathlib
import re
import signal
import sys
import time

from core import __AUTHOR__, __RELEASE__, __version__, resource
from core.analysis import obs
from core.io import cli
from core.io import logger as log
from core.tune import core as icetune_core
from core.tune import push as icetune_push
from core.tune import short_id
from core.tune.cache import PermanentConfigurationError
from core.tune.drivers.base import SimulatorDriver
from core.tune.drivers.registry import create_driver, default_driver_name, get_driver_type, infer_driver_name
from core.tune.parameters import space as parameter_space
from core.tune.parameters import tools
from core.tune.runtime import process as iceruntime
from core.tune.tunesetup import load_tunesetup, save_tunesetup

logger = log.get_logger(__name__)
DEFAULT_ICEBO_CONFIG = resource("tune/settings/icebo.json")
BACKEND_ALGORITHMS = {
    "ray": {
        "ampfit",
        "ax",
        "basic",
        "bayesopt",
        "hebo",
        "hyperopt",
        "icebo",
        "lbfgs",
        "optuna",
        "scikit",
    }
}


# Resolve optimizer settings once for preflight, distributed snapshots and run history
def resolve_optimizer_settings(args, parser: argparse.ArgumentParser) -> None:
    from importlib import import_module

    loaders = {"hebo": "load_hebo_config", "icebo": "load_icebo_config", "ampfit": "load_settings", "lbfgs": "load_settings"}
    for name, loader in loaders.items():
        setattr(args, name + "_settings", None)
        if args.algorithm != name:
            continue
        value = getattr(args, name + "_config")
        source = resource(value.removeprefix("package:")) if value.startswith("package:") else pathlib.Path(value)
        source = (pathlib.Path(args.cdir) / source).resolve()
        try:
            module = import_module(f"core.tune.optimizers.{name}.config")
            setattr(args, name + "_settings", getattr(module, loader)(source))
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        setattr(args, name + "_config", str(source))


# Parse and validate the optimizer and driver settings
def parse_arguments():
    """Parse and normalize the icetune command line"""

    selector = argparse.ArgumentParser(add_help=False)
    cli.add_value(selector, "simdriver", None)
    selected, _ = selector.parse_known_args()
    selected_driver_name = (
        str(selected.simdriver).upper() if selected.simdriver else default_driver_name()
    )
    selected_driver_type = SimulatorDriver if not selected.simdriver and any(flag in sys.argv for flag in ("-h", "--help")) else get_driver_type(selected_driver_name)
    parser = argparse.ArgumentParser(
        description=f"%(prog)s {__version__} {__RELEASE__} [{__AUTHOR__}]",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="command:\n  history SOURCE  Print completed trial costs",
    )

    # System
    cli.add_value(parser, "cdir", os.getcwd(), help="Main absolute directory")
    cli.add_value(
        parser,
        "libdir",
        selected_driver_type.default_library_path(),
        help="Auxiliary library directory",
    )
    cli.add_value(parser, "backend", "ray", help="Tuning backend")

    # Ray backend
    cli.add_value(parser, "address", None, help="Ray node address or local")
    cli.add_value(parser, "ray_worker_cpu", None, value_type=int, help="CPU capacity of the local worker allocation")
    cli.add_value(parser, "ray_worker_gpu", 0, value_type=int, help="GPU capacity of the local worker allocation")
    cli.add_value(parser, "ray_worker_poll_interval_s", 60, value_type=int, help="Ray worker poll interval s")
    cli.add_value(parser, "ray_worker_failure_limit", 2, value_type=int, help="Ray worker failure limit")
    cli.add_value(parser, "ray_worker_recovery_timeout_s", 600, value_type=int, help="Ray worker recovery timeout s")
    cli.add_value(parser, "ray_restore_retry_limit", 2, value_type=int, help="Ray restore retry limit")
    cli.add_value(parser, "ray_trial_retry_limit", 3, value_type=int, help="Ray trial retry limit")
    cli.add_value(parser, "ray_restore_retry_delay_s", 10, value_type=int, help="Ray restore retry delay s")
    cli.add_value(parser, "ray_resource_status_interval_s", 5, value_type=int, help="Ray resource status interval s")
    cli.add_value(parser, "ray_status_interval_s", 60, value_type=int, help="Ray status interval s")
    cli.add_value(parser, "ray_worker_admission_timeout_s", 600, value_type=int,
                     help="Ray worker admission probe timeout in seconds")
    cli.add_value(parser, "ray_temp_dir", None, help="Local Ray session directory")
    cli.add_value(parser, "ray_storage_path", None, help="Persistent Tune storage path")
    cli.add_value(parser, "ray_init_state", None, help="Completed initialization state")
    cli.add_value(parser, "ray_runtime_uri", None, help="Prepared runtime archive URI accessible to every Ray node")
    cli.add_flag(
        parser,
        "ray_upload_runtime",
        help="Upload the runtime to Ray nodes instead of requiring a shared checkout",
    )

    # Backend initialization
    cli.add_value(parser, "precompute", True, value_type=cli.parse_bool, help="Precompute shared simulator grids and MC covariance")
    cli.add_value(parser, "coord_dir", None, help="Initialization state directory")
    cli.add_value(parser, "history_interval_s", 60.0, help="History and plot refresh interval")
    cli.add_value(
        parser,
        "proposal_batch_size",
        0,
        help="Maximum joint acquisition block, zero selects it automatically",
    )
    cli.add_value(
        parser,
        "proposal_batch_fraction",
        0.1,
        help="Automatic adaptive proposal fraction",
    )
    cli.add_value(
        parser,
        "phase",
        None,
        choices=("bank", "init"),
        help="Explicit backend initialization phase",
    )
    cli.add_value(parser, "bootstrap_cache_url", None, help="Initialization cache URL")
    cli.add_value(parser, "runtime_sha256", None, help="Static runtime archive checksum")
    cli.add_value(parser, "bank_index", None, value_type=int, help="Active amplitude sample index")
    cli.add_value(parser, "bank_shard", None, value_type=int, help="Amplitude bank shard index")
    cli.add_value(parser, "bank_jobs", None, value_type=int, help="Number of amplitude bank shards")
    cli.add_value(parser, "bank_shared_dir", None, help="Shared filesystem root for distributed bank phases")
    cli.add_flag(parser, "bank_sample", help="Generate the shared amplitude bank event sample")
    cli.add_flag(parser, "bank_finalize", help="Finalize one sample after its amplitude shards complete")

    # Computational resources
    cli.add_value(parser, "cpu_per_trial", -1, help="CPUs per trial; -1 selects automatically")
    cli.add_value(parser, "gpu_per_trial", 0, help="GPUs per trial")
    cli.add_value(parser, "max_concurrent_trials", 4, help="Maximum simultaneous trials")
    cli.add_value(parser, "wait_workers", None, value_type=int,
                     help="Required Ray worker allocations, excluding the dedicated head")
    cli.add_value(parser, "wait_workers_timeout_s", 300.0, help="Ray worker allocation wait timeout")
    cli.add_value(parser, "chunksize", 10000, help="Events read per process")

    # icetune steering input
    cli.add_value(parser, "obs_module", "default", help="Fallback observable module")
    cli.add_value(parser, "tunesetup", "tunesetup.json", help="JSON tuning definition, optionally file.json#/pointer")
    cli.add_value(parser, "simdriver", None, help="Simulator driver")
    cli.add_value(parser, "plot-brand", None, help="Plot producer name")
    cli.add_value(
        parser,
        "push",
        None,
        nargs="+",
        metavar="VALUE",
        help="Publish parameters: INPUT TARGET_PATH [DRIVER_OPTIONS_JSON]",
    )

    # MC generator steering
    cli.add_value(parser, "tune_default", selected_driver_type.DEFAULT_TUNE, help="Default MC model tune")
    cli.add_value(parser, "rngseed", 1234, help="Random seed")
    cli.add_value(parser, "run_name", None, help="Run name")

    # icetune technical controls
    cli.add_flag(parser, "init_force", help="Force initialization of precomputed MC arrays")
    cli.add_value(parser, "ampfit_reuse", True, value_type=cli.parse_bool, help="Reuse compatible ampfit banks from the current or previous runs")
    cli.add_flag(parser, "restore", help="Restore the selected run")
    cli.add_flag(parser, "no_initial_point", help="Disable initial parameter values")
    cli.add_flag(parser, "no_fit", help="Disable fitting")
    cli.add_flag(parser, "preflight", help="Validate the tuning card and initial fit point, then exit")
    cli.add_value(parser, "save_tunesetup", None, help="Save the resolved tuning definition as JSON")
    cli.add_value(
        parser,
        "plot",
        True,
        value_type=cli.parse_bool,
        help="Enable plotting <true|false>",
    )
    cli.add_flag(parser, "pickle_dump", help="Write full trial pickle payloads")
    cli.add_flag(parser, "save_events", help="Retain complete GRANIITTI trial event samples")
    cli.add_value(parser, "events_dir", None, help="Shared destination for completed trial event samples")
    cli.add_value(parser, "pickle_dir", None, help="Full trial pickle output directory")
    cli.add_value(parser, "verbose", 1, help="Verbosity level")

    # Optimization algorithm
    cli.add_value(parser, "algorithm", "icebo", help="Optimization algorithm")
    cli.add_value(parser, "ampfit_config", str(resource("tune/settings/ampfit.json")), help="Bounded Torch L-BFGS configuration")
    cli.add_value(
        parser, "hebo_config", str(resource("tune/settings/hebo.json")), help="HEBO JSON configuration",
    )
    cli.add_value(
        parser,
        "icebo_config",
        str(DEFAULT_ICEBO_CONFIG),
        help="Authoritative ICEBO JSON configuration",
    )

    # Cost function
    cli.add_value(parser, "cost", "chi2", help="Objective function")
    cli.add_value(parser, "cost_rho", "quadratic", help="Objective residual shape")
    cli.add_value(parser, "cost_avg", "global-mean", help="Objective combiner")

    # Sampling
    cli.add_value(parser, "rand_trials", 1000, help="Random exploration trials")
    cli.add_value(parser, "num_trials", 100000, help="Target successful trials")
    cli.add_value(parser, "max_t", 3600, help="Maximum complete-trial runtime")
    cli.add_value(
        parser,
        "surrogate_fit_percentile",
        1.0,
        help="Lowest-cost trial fraction used for surrogate fitting",
    )

    cli.add_value(parser, "lbfgs_config", str(resource("tune/settings/lbfgs.json")), help="Distributed L-BFGS configuration")

    selected_driver_type.add_cli_arguments(parser)
    args = parser.parse_args()
    if args.push is not None and len(args.push) not in {2, 3}:
        parser.error("--push requires INPUT TARGET_PATH [DRIVER_OPTIONS_JSON]")
    args.backend = args.backend.lower()
    args.algorithm = args.algorithm.lower()
    if args.backend not in BACKEND_ALGORITHMS:
        parser.error(f'unsupported backend "{args.backend}"')
    if args.algorithm not in BACKEND_ALGORITHMS[args.backend]:
        parser.error(f'unsupported {args.backend} optimizer "{args.algorithm}"')
    resolve_optimizer_settings(args, parser)
    if args.algorithm == "ampfit" and not selected_driver_type.AMPLITUDE_FIT:
        parser.error("ampfit requires a driver with amplitude fitting support")
    if args.algorithm == "ampfit" and args.no_initial_point:
        parser.error("ampfit relative starts require initial parameter values")
    if not 0.0 < args.surrogate_fit_percentile <= 1.0:
        parser.error("--surrogate_fit_percentile must be in the interval (0, 1]")
    if args.proposal_batch_size < 0:
        parser.error("--proposal_batch_size must be non-negative")
    if not 0.0 < args.proposal_batch_fraction <= 1.0:
        parser.error("--proposal_batch_fraction must be in the interval (0, 1]")
    if args.max_concurrent_trials <= 0:
        parser.error("--max_concurrent_trials must be positive")
    if args.num_trials <= 0 or not 0 <= args.rand_trials <= args.num_trials:
        parser.error("require --num_trials > 0 and 0 <= --rand_trials <= --num_trials")
    if args.algorithm in {"hebo", "icebo"} and args.rand_trials >= args.num_trials:
        parser.error("adaptive optimization requires --num_trials > --rand_trials")
    if args.rngseed < 0 or args.gpu_per_trial < 0:
        parser.error("--rngseed and --gpu_per_trial must be non-negative")
    if args.max_t <= 0 or args.chunksize <= 0:
        parser.error("--max_t and --chunksize must be positive")
    selected_driver_type.validate_cli_arguments(parser, args)
    if not args.precompute and (args.phase in {"bank", "init"} or args.ray_init_state or args.init_force):
        parser.error("--precompute false cannot be combined with INIT state, phase init or init_force")
    if args.phase != "bank" and (args.bank_sample or args.bank_finalize or args.bank_shard is not None or args.bank_index is not None):
        parser.error("--bank_sample, --bank_finalize, --bank_index and --bank_shard require --phase bank")
    if args.bank_index is not None and args.bank_index < 0:
        parser.error("--bank_index must be nonnegative")
    if args.bank_jobs is not None and (args.bank_jobs < 1 or args.algorithm != "ampfit" or args.phase not in {"bank", "init"}):
        parser.error("--bank_jobs requires ampfit, --phase bank or init, and a positive shard count")
    if args.phase == "bank":
        if args.algorithm != "ampfit":
            parser.error("--phase bank requires the ampfit algorithm")
        if args.bank_sample and args.bank_finalize:
            parser.error("--bank_sample and --bank_finalize are mutually exclusive")
        if args.bank_finalize and (args.bank_index is None or args.bank_jobs is None):
            parser.error("--bank_finalize requires --bank_index and --bank_jobs")
        if args.bank_sample or args.bank_finalize:
            if args.bank_shard is not None:
                parser.error("--bank_sample and --bank_finalize cannot be combined with --bank_shard")
        elif args.bank_shard is None or args.bank_jobs is None:
            parser.error("--phase bank requires --bank_shard and --bank_jobs")
        elif not 0 <= args.bank_shard < args.bank_jobs:
            parser.error("--bank_shard must satisfy 0 <= shard < --bank_jobs")

    if args.run_name is None:
        timestamp = iceruntime.get_current_time()
        tune_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", pathlib.Path(args.tunesetup).stem)
        args.run_name = f"tunesetup_{tune_name}_algorithm_{args.algorithm}_{timestamp}"
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.run_name) is None:
        parser.error("--run_name must be one safe path component")
    if args.phase in {"bank", "init"} and args.coord_dir is None:
        args.coord_dir = os.path.join(args.cdir, "runs", "icetune", args.run_name, "ray", "init")

    args.auto_cpu_per_trial = args.cpu_per_trial <= 0
    if args.auto_cpu_per_trial:
        args.cpu_per_trial = max(1, multiprocessing.cpu_count() // args.max_concurrent_trials)
    if args.backend == "ray" and args.address is None:
        args.address = "local"
    if args.ray_worker_cpu is not None and args.ray_worker_cpu <= 0:
        parser.error("--ray_worker_cpu must be positive")
    if args.ray_worker_gpu < 0:
        parser.error("--ray_worker_gpu must be non-negative")
    if args.ray_worker_poll_interval_s <= 0:
        parser.error("--ray_worker_poll_interval_s must be positive")
    if args.ray_worker_failure_limit <= 0:
        parser.error("--ray_worker_failure_limit must be positive")
    if args.ray_worker_recovery_timeout_s <= 0:
        parser.error("--ray_worker_recovery_timeout_s must be positive")
    if args.ray_restore_retry_limit < 0:
        parser.error("--ray_restore_retry_limit must be non-negative")
    if args.ray_trial_retry_limit < 0:
        parser.error("--ray_trial_retry_limit must be non-negative")
    if args.ray_restore_retry_delay_s <= 0:
        parser.error("--ray_restore_retry_delay_s must be positive")
    if args.ray_resource_status_interval_s <= 0:
        parser.error("--ray_resource_status_interval_s must be positive")
    if args.ray_status_interval_s <= 0:
        parser.error("--ray_status_interval_s must be positive")
    if args.ray_worker_admission_timeout_s <= 0:
        parser.error("--ray_worker_admission_timeout_s must be positive")
    if args.backend == "ray" and args.ray_temp_dir is None:
        args.ray_temp_dir = os.path.join(args.cdir, "tmp", "ray")
    args.ray_temp_dir = cli.absolute_path(args.ray_temp_dir)
    args.ray_storage_path = cli.absolute_path(args.ray_storage_path)
    args.pickle_dir = cli.absolute_path(args.pickle_dir)
    args.events_dir = cli.absolute_path(args.events_dir)
    if args.save_events and (selected_driver_name.upper() != "GRANIITTI" or args.backend != "ray"
                            or args.algorithm in {"ampfit", "lbfgs"}):
        parser.error("--save_events requires GRANIITTI Ray trials generating fresh events, not ampfit or lbfgs")
    if args.simdriver is not None or args.push is None:
        args.simdriver = selected_driver_name
    return args


# Publish one summary through its selected driver and print the physical changes
def run_push(args):
    """Run the driver-generic parameter push operation."""

    input_path = pathlib.Path(args.push[0]).expanduser()
    if not input_path.is_absolute():
        input_path = pathlib.Path(args.cdir) / input_path
    summary = icetune_push.load_summary(input_path)
    options = icetune_push.parse_driver_options(args.push[2]) if len(args.push) == 3 else {}
    print(f'Push parameter source: {summary["_push_parameter_source"]} in "{input_path}"')
    print(f"Push driver options: {json.dumps(options, sort_keys=True)}")
    args.simdriver = str(args.simdriver).upper() if args.simdriver else infer_driver_name(summary)
    try:
        create_driver(args.simdriver).push_parameters(
            summary=summary,
            target_path=args.push[1],
            cdir=args.cdir,
            options=options,
            confirm=icetune_push.confirm_change_table,
        )
    except icetune_push.PushCancelled:
        print("Push cancelled; target files were not changed.")
        return False
    print("Push completed.")
    return True


def initialize_backend_context(*, args, tunesetup, simdriver, mc_steer):
    """
    Initialize simulator state and backend-specific bootstrap helpers.
    """

    if args.phase == "init":
        return simdriver.initialize_backend(
            args=args,
            tunesetup=tunesetup,
            mc_steer=mc_steer,
        )

    if args.phase == "bank":
        bank_mode = ({"phase": "sample", "jobs": args.bank_jobs or 1} if args.bank_sample
                     else {"phase": "shard", "shard": args.bank_shard, "jobs": args.bank_jobs})
        if args.bank_finalize:
            bank_mode = {"phase": "finalize", "jobs": args.bank_jobs}
        if args.bank_index is not None:
            bank_mode["index"] = args.bank_index
        if args.bank_shared_dir:
            bank_mode["shared_root"] = args.bank_shared_dir
        simdriver.initialize(
            run_name=args.run_name,
            tunesetup=tunesetup,
            mc_steer=mc_steer,
            obs_module=args.obs_module,
            cdir=args.cdir,
            init_force=args.init_force,
            max_t=args.max_t,
            pickle_dump=False,
            processes=args.cpu_per_trial,
            rngseed=args.rngseed,
            bank_mode=bank_mode,
        )
        return {"bootstrapper": None, "bootstrap_fingerprint_getter": None,
                "initial_point_getter": None, "initial_points": None}

    if not args.precompute:
        simdriver.init_data(
            run_name=args.run_name, datacards=tunesetup.datacards, obs_module=args.obs_module,
            cdir=args.cdir, pickle_dump=args.pickle_dump,
        )
        initial_points = simdriver.get_initial_param(
            param_space=tunesetup.param_space, aux_param_space=tunesetup.aux_param_space,
            cdir=args.cdir, tune_default=mc_steer["tune_default"],
        )
    elif args.backend == "ray" and args.ray_init_state:
        initial_points = simdriver.initialize_ray_bootstrap(
            args=args,
            tunesetup=tunesetup,
            mc_steer=mc_steer,
        )
    else:
        initial_points = simdriver.initialize(
            run_name=args.run_name,
            tunesetup=tunesetup,
            mc_steer=mc_steer,
            obs_module=args.obs_module,
            cdir=args.cdir,
            init_force=args.init_force,
            max_t=args.max_t,
            pickle_dump=args.pickle_dump,
            processes=args.cpu_per_trial,
            rngseed=args.rngseed,
        )

    return {
        "bootstrapper": None,
        "initial_point_getter": None,
        "initial_points": initial_points,
    }


# Validate the driver default point before initialization or cluster startup
def preflight_initial_config(*, args, tunesetup, simdriver, mc_steer) -> dict | None:
    if args.no_initial_point:
        print_preflight_parameter_table(param_space=tunesetup.param_space, initial=None)
        return None
    initial = simdriver.get_initial_param(
        param_space=tunesetup.param_space,
        aux_param_space=tunesetup.aux_param_space,
        cdir=args.cdir,
        tune_default=mc_steer["tune_default"],
    )
    print_preflight_parameter_table(param_space=tunesetup.param_space, initial=initial)
    try:
        return parameter_space.validate_initial_config(initial, tunesetup.param_space)
    except (KeyError, TypeError, ValueError) as exc:
        raise PermanentConfigurationError(str(exc)) from exc


# Format one preflight table scalar without truncating significant bounds
def _preflight_scalar(value) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, int):
        return str(value)
    try:
        return f"{float(value):.12g}"
    except (TypeError, ValueError):
        return str(value)


# Print every active optimizer coordinate with its initial value and bounds
def print_preflight_parameter_table(*, param_space: dict, initial: dict | None) -> None:
    bounds = parameter_space.normalize_param_space(param_space)
    headers = ("#", "parameter", "initial", "lower", "upper", "type")
    rows = [
        (
            str(index),
            name,
            _preflight_scalar(None if initial is None else initial.get(name)),
            _preflight_scalar(spec["lower"]),
            _preflight_scalar(spec["upper"]),
            str(spec["type"]),
        )
        for index, (name, spec) in enumerate(sorted(bounds.items()), start=1)
    ]
    widths = [
        max([len(headers[column]), *(len(row[column]) for row in rows)])
        for column in range(len(headers))
    ]

    def render(row) -> str:
        return "  ".join(value.ljust(widths[column]) for column, value in enumerate(row)).rstrip()

    print(f"ICETUNE tunable parameters ({len(rows)})")
    print(render(headers))
    print(render(tuple("-" * width for width in widths)))
    for row in rows:
        print(render(row))


def build_common_param(*, args, tunesetup, mc_steer, python_version):
    """Build the backend-independent trial parameter payload."""

    bounds = parameter_space.normalize_param_space(tunesetup.param_space)
    parameter_metadata = icetune_core.build_parameter_space(bounds)
    optimization = icetune_core.build_optimization_metadata(args)
    parameter_topology = tools.build_parameter_topology(tunesetup.param_space.keys())

    return {
        "aux_param_space": copy.deepcopy(tunesetup.aux_param_space),
        "datacards": copy.deepcopy(tunesetup.datacards),
        "mc_steer": copy.deepcopy(mc_steer),
        "cost": copy.deepcopy(args.cost),
        "cost_definition": icetune_core.optimizer_cost_definition(args.cost),
        "cost_rho": copy.deepcopy(args.cost_rho),
        "cost_avg": copy.deepcopy(args.cost_avg),
        "data_covariance_mode": copy.deepcopy(getattr(args, "data_covariance_mode", "diagonal")),
        "obs_module": copy.deepcopy(args.obs_module),
        "cdir": copy.deepcopy(args.cdir),
        "libdir": copy.deepcopy(args.libdir),
        "PYTHON_VERSION": python_version,
        "backend": copy.deepcopy(args.backend),
        "parameter_space": copy.deepcopy(parameter_metadata),
        "parameter_topology": copy.deepcopy(parameter_topology),
        "optimization": copy.deepcopy(optimization),
        "simdriver": copy.deepcopy(args.simdriver),
        "plot_brand": copy.deepcopy(args.plot_brand),
        "chunksize": copy.deepcopy(args.chunksize),
        "processes": copy.deepcopy(args.cpu_per_trial),
        "rngseed": int(args.rngseed),
        "runtime_sha256": copy.deepcopy(args.runtime_sha256),
        "ray_upload_runtime": bool(args.ray_upload_runtime),
        "ray_worker_admission_timeout_s": args.ray_worker_admission_timeout_s,
        "ray_worker_poll_interval_s": args.ray_worker_poll_interval_s,
        "ray_worker_failure_limit": args.ray_worker_failure_limit,
        "ray_worker_recovery_timeout_s": args.ray_worker_recovery_timeout_s,
        "ray_restore_retry_limit": args.ray_restore_retry_limit,
        "ray_trial_retry_limit": args.ray_trial_retry_limit,
        "ray_restore_retry_delay_s": args.ray_restore_retry_delay_s,
        "ray_resource_status_interval_s": args.ray_resource_status_interval_s,
        "ray_status_interval_s": args.ray_status_interval_s,
        "plot": args.plot,
        "pickle_dump": args.pickle_dump,
        "save_events": args.save_events and args.phase == "fit",
        "events_dir": args.events_dir,
        "pickle_dir": copy.deepcopy(args.pickle_dir),
        "run_name": copy.deepcopy(args.run_name),
        "max_t": copy.deepcopy(args.max_t),
    }


# Run the selected distributed backend through a lazy dispatcher
def run_backend(*, args, tunesetup, simdriver, initial_points, param, storage_path, start_time):
    # Run Ray without importing its runtime before backend selection
    def run_ray() -> None:
        from core.tune.runtime.ray import configure_ray_runtime

        configure_ray_runtime()
        from core.tune.backends import ray as icetune_ray

        icetune_ray.run_ray_backend(
            args=args,
            tunesetup=tunesetup,
            simdriver=simdriver,
            initial_points=initial_points,
            param=param,
            storage_path=storage_path,
            start_time=start_time,
        )

    runners = {"ray": run_ray}
    try:
        runner = runners[args.backend]
    except KeyError:
        raise ValueError(f'Unsupported backend "{args.backend}"') from None
    runner()


# Run history inspection or the tuning command
def main():
    if len(sys.argv) > 1 and sys.argv[1] == "history":
        from core.tune.history import main as history
        return history(sys.argv[2:])
    wall_start = time.perf_counter()

    python_version = f"{sys.version_info.major}.{sys.version_info.minor}"

    # Set the environment variable to disable strict metric checking
    os.environ["TUNE_DISABLE_STRICT_METRIC_CHECKING"] = "1"

    args = parse_arguments()
    log.configure(level=log.level_from_verbosity(args.verbose), force=True)
    iceruntime.configure_numerical_threads(args.cpu_per_trial)

    def signal_handler(sig, frame):
        del frame
        logger.warning(".icetune: received signal %s, stopping cleanly", sig)
        raise SystemExit(128 + int(sig))

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    logger.info(".icetune: argv=%s", " ".join(sys.argv))
    if args.push is not None:
        logger.debug(".icetune: args=%s", args)
        run_push(args)
        elapsed_seconds = time.perf_counter() - wall_start
        logger.info(iceruntime.format_done_message("icetune", elapsed_seconds))
        return

    simdriver = create_driver(args.simdriver)

    if not args.preflight:
        logger.info(".icetune: IP=%s", iceruntime.get_ip())
    if args.auto_cpu_per_trial:
        logger.info(".icetune: Set automatic cpu_per_trial=%s", args.cpu_per_trial)
    logger.debug(".icetune: args=%s", args)

    # Fail before initialization or trial execution if the worker Python runtime is incomplete
    if not args.preflight:
        obs.preflight_numba_runtime()
        logger.info(".icetune: Numba runtime preflight completed")

    iceruntime.set_random_seeds(args.rngseed)

    mc_steer = simdriver.build_run_steering(args)
    saved_tunesetup = pathlib.Path(args.cdir).resolve() / "runs" / "icetune" / args.run_name / "tunesetup.json"
    tunesetup_input = str(saved_tunesetup) if args.restore and saved_tunesetup.is_file() else args.tunesetup
    tunesetup = load_tunesetup(
        cdir=args.cdir,
        simdriver=args.simdriver,
        name=tunesetup_input, tune_default=args.tune_default,
    )
    if getattr(tunesetup, "tune_default", args.tune_default) != args.tune_default:
        raise ValueError("The saved tunesetup uses a different baseline model tune")
    simdriver.prepare_tunesetup(tunesetup=tunesetup, args=args)
    print(f"icetune baseline tune: {args.tune_default}")
    print(f"icetune datasets: {json.dumps(tunesetup.datacards, indent=2)}")
    print(f"icetune fixed parameters: {json.dumps(tunesetup.aux_param_space, indent=2)}")
    initial = preflight_initial_config(
        args=args,
        tunesetup=tunesetup,
        simdriver=simdriver,
        mc_steer=mc_steer,
    )
    if initial is not None:
        logger.info(".icetune: Initial fit configuration is inside all parameter bounds")
    if args.save_tunesetup or not args.preflight:
        tunesetup.optimizer = {args.algorithm: getattr(args, args.algorithm + "_settings", None)}
        definition_path = pathlib.Path(args.save_tunesetup) if args.save_tunesetup else saved_tunesetup
        save_tunesetup(tunesetup=tunesetup, path=definition_path, simdriver=args.simdriver, tune_default=args.tune_default)
        print(f"icetune resolved tunesetup: {definition_path}")
        tunesetup = load_tunesetup(cdir=args.cdir, simdriver=args.simdriver, name=str(definition_path.resolve()))
    if args.preflight:
        logger.info(".icetune: Preflight completed, exit")
        return
    backend_context = initialize_backend_context(
        args=args,
        tunesetup=tunesetup,
        simdriver=simdriver,
        mc_steer=mc_steer,
    )
    if args.phase == "bank":
        if args.bank_finalize:
            logger.info(".icetune: Amplitude sample %s finalized", args.bank_index)
        elif args.bank_sample:
            logger.info(".icetune: Amplitude bank event sample completed")
        else:
            logger.info(".icetune: Amplitude bank shard %s/%s completed", args.bank_shard, args.bank_jobs)
        return
    param = build_common_param(
        args=args,
        tunesetup=tunesetup,
        mc_steer=mc_steer,
        python_version=python_version,
    )

    if args.no_fit:
        logger.info(".icetune: No-fit mode chosen, exit.")
        sys.exit(0)
    simdriver.prepare_trial_runtime(param)

    logger.info('.icetune: Launching backend "%s" ...', args.backend)
    start = time.time()
    storage_path = args.ray_storage_path or os.path.join(args.cdir, "runs", "icetune")
    logger.info('.icetune: Results storage path "%s"', storage_path)

    if args.phase == "init":
        from core.tune import init as icetune_init

        results = icetune_init.run_initialization(
            args=args,
            bootstrapper=backend_context["bootstrapper"],
            fingerprint_getter=backend_context["bootstrap_fingerprint_getter"],
            initial_point_getter=backend_context["initial_point_getter"],
        )
        logger.info(
            ".icetune: %s initialization completed: %s",
            args.backend,
            short_id(results["bootstrap"]["fingerprint"]),
        )

    else:
        run_backend(
            args=args,
            tunesetup=tunesetup,
            simdriver=simdriver,
            initial_points=backend_context["initial_points"],
            param=param,
            storage_path=storage_path,
            start_time=start,
        )

    elapsed_seconds = time.perf_counter() - wall_start
    logger.info(iceruntime.format_done_message("icetune", elapsed_seconds))


if __name__ == "__main__":
    try:
        main()
    except PermanentConfigurationError as exc:
        logger.error(".icetune: permanent configuration error: %s", exc)
        sys.exit(64)

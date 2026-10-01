# Test icetune Ray and lbfgs optimization backends
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import errno
import hashlib
import json
import os
import pickle
import runpy
import socket
import sys
import time
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from core import resource

ROOT = Path(__file__).resolve().parents[3]

ray = pytest.importorskip("ray")
torch = pytest.importorskip("torch")
from core.tune import likelihood as icetune_likelihood
from core.tune.backends import ray as icetune_ray
from core.tune.optimizers.hebo.config import load_hebo_config
from core.tune.optimizers.icebo.config import load_icebo_config, override_icebo_config
from core.tune.optimizers.lbfgs import optimizer as icetune_lbfgs
from core.tune.optimizers.lbfgs.config import load_settings as load_lbfgs_settings
from core.tune.parameters import space as parameter_space
from core.tune.runtime import ray as ray_runtime
from ray import tune


class _ToyTunesetup:
    param_space = {
        "y": tune.uniform(-2.0, 2.0),
        "x": tune.uniform(-2.0, 2.0),
    }


class _ChoiceTunesetup:
    param_space = {"x": tune.choice([0.0, 1.0])}


class _ToyDriver:
    # Keep the analytic driver free of external input files
    def shared_runtime_files(self, param):
        return []

    # Render only in a worker process for the full lbfgs backend test
    def render_trial_figures_to_dir(self, *, outputs, param, summary_payload, output_dir, summary_file=None):
        assert os.getpid() != param["head_pid"]
        Path(output_dir, "fit.txt").write_text(str(outputs["metrics"]["loss"]), encoding="utf-8")
        Path(summary_file).write_text(json.dumps(summary_payload), encoding="utf-8")
        return summary_payload

    # Compute a synthetic driver-owned trial name
    def trial_tunename(self, *, trial_id, node_id=None, pid=None, purpose="trial"):
        return f"test-{purpose}-{node_id}-{pid}-{trial_id}"

    # Accept worker runtime environment requests
    def runtime_environment(self, **kwargs):
        return {"LD_LIBRARY_PATH": "", "PYTHONPATH": ""}

    # Compute a deterministic quadratic objective for one synthetic trial
    def evaluate_trial_outputs(self, *, config, param, trial_id, tunename):
        loss = (float(config["x"]) - 0.4) ** 2 + (float(config["y"]) + 0.2) ** 2
        return {
            "config": dict(config),
            "metrics": {"loss": float(loss)},
            "results": {"ok": True},
            "trial_id": trial_id,
            "tunename": tunename,
        }

    # Accept cleanup calls from core.icetune helpers
    def cleanup_trial_outputs(self, **kwargs):
        return None


class _RayFigureDriver:
    # Keep the analytic driver free of external input files
    def shared_runtime_files(self, param):
        return []

    # Compute a synthetic driver-owned trial name
    def trial_tunename(self, *, trial_id, node_id=None, pid=None, purpose="trial"):
        return f"test-{purpose}-{node_id}-{pid}-{trial_id}"

    # Compute a deterministic one-dimensional objective for a Ray figure test
    def evaluate_trial_outputs(self, *, config, param, trial_id, tunename):
        loss = (float(config["x"]) - 0.2) ** 2
        return {
            "config": dict(config),
            "metrics": {"loss": float(loss)},
            "results": {"ok": True},
            "trial_id": trial_id,
            "tunename": tunename,
        }

    # Render a minimal non-empty figure payload for publication tests
    def render_trial_figures_to_dir(
        self, *, outputs, param, summary_payload, output_dir, summary_file=None
    ):
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        Path(output_dir, "output.txt").write_text(summary_payload["trial_id"], encoding="utf-8")
        if summary_file is not None:
            Path(summary_file).write_text(json.dumps(summary_payload), encoding="utf-8")
        return dict(summary_payload)

    # Accept cleanup calls from core.icetune helpers
    def cleanup_trial_outputs(self, **kwargs):
        return None


class _ImmediateRemoteMethod:
    # Wrap a local callable with the Ray actor .remote() call shape
    def __init__(self, func):
        self.func = func

    # Execute the wrapped callable synchronously
    def remote(self, *args, **kwargs):
        return self.func(*args, **kwargs)


# Verify the render selector marks only strict finite improvements
def test_ray_best_state_selects_worker_renders():
    state = icetune_ray.BestState()

    assert state.finalize_trial_result("trial-1", 4.0)["locally_improving"]
    assert not state.finalize_trial_result("trial-2", 5.0)["locally_improving"]
    assert state.finalize_trial_result("trial-3", 2.0)["locally_improving"]
    snapshot = state.snapshot()
    assert snapshot["best_cost"] == 2.0
    assert snapshot["best_trial_id"] == "trial-3"


# Check render selection never reserves or blocks a simulation slot
def test_ray_best_state_selects_pending_hists():
    state = icetune_ray.BestState(render_figures=True)
    assert state.finalize_trial_result("trial-1", 10.0)["store_render"]
    assert state.finalize_trial_result("trial-2", 20.0)["store_render"]
    state.published_cost = 10.0
    assert not state.finalize_trial_result("trial-2", 20.0)["store_render"]
    assert state.finalize_trial_result("trial-3", 5.0)["store_render"]


# Check Ray transfers durable trial outputs without touching head-private storage
def test_ray_worker_transfers_written_trial_outputs(tmp_path, monkeypatch):
    publications = []
    outputs = {
        "metrics": {"loss": 1.0},
        "trial_id": "trial-1",
        "tunename": "TUNE_trial-1",
    }
    state = SimpleNamespace(
        finalize_trial_result=_ImmediateRemoteMethod(
            lambda trial_id, cost: {"locally_improving": False}
        ),
        publish_trial_outputs=_ImmediateRemoteMethod(
            lambda **kwargs: publications.append(kwargs) or {"output_token": 1}
        ),
    )
    monkeypatch.setattr(icetune_ray.ray, "get", lambda value: value)
    monkeypatch.setattr(
        icetune_ray.icetune_main,
        "evaluate_trial_outputs",
        lambda **kwargs: outputs,
    )
    monkeypatch.setattr(
        icetune_ray.icetune_main,
        "cleanup_trial_outputs",
        lambda **kwargs: None,
    )
    monkeypatch.setattr(
        icetune_ray.icetune_main,
        "maybe_dump_trial_payload",
        lambda **kwargs: {
            "filename": "TUNE_icetune_trial-1.pkl",
            "relative_path": "results/TUNE_icetune_trial-1.pkl",
            "sha256": "a" * 64,
            "size": 1,
        },
    )
    monkeypatch.setattr(
        icetune_ray,
        "collect_ray_trial_outputs",
        lambda **kwargs: {
            "figures": None,
            "pickle": b"x" if kwargs["descriptor"] is not None else None,
        },
    )

    param = {
        "cost": "loss",
        "pickle_dump": False,
        "plot": False,
        "run_name": "ray-transfer",
    }
    trainable = SimpleNamespace(
        config={},
        deadline=time.monotonic() + 60.0,
        started_at=time.time(),
        setup_error=None,
        global_state=state,
        initial_points=None,
        logdir=str(tmp_path),
        output_param=param,
        param=param,
        save=lambda: pytest.fail("worker checkpoints must stay disabled"),
        simdriver=object(),
        trial_id="trial-1",
    )
    icetune_ray.CFunc.step(trainable)
    assert publications == []

    param["pickle_dump"] = True
    result = icetune_ray.CFunc.step(trainable)
    assert len(publications) == 1
    assert publications[0]["transfer"]["pickle"] == b"x"
    assert result["trial_pickle"]["relative_path"].startswith("results/")


# Build the minimal argparse-like object needed by the lbfgs backend
def _args(*, num_trials=80):
    return SimpleNamespace(
        algorithm="lbfgs",
        lbfgs_settings=load_lbfgs_settings(resource("tune/settings/lbfgs.json")),
        cpu_per_trial=1,
        gpu_per_trial=0,
        max_concurrent_trials=2,
        no_initial_point=False,
        num_trials=num_trials,
        rngseed=123,
    )


# Build one Ray searcher over the scalar toy bound
def _search(args, initial_points=None):
    bounds = {"x": {"type": "uniform", "lower": 0.0, "upper": 1.0}}
    search = icetune_ray.PercentileSearch(args=args, bounds=bounds, initial_points=initial_points)
    search.param_space = {"x": tune.uniform(0.0, 1.0)}
    return search


# Build shared asynchronous search arguments with explicit overrides
def _search_args(algorithm, **overrides):
    values = {
        "cost": "loss",
        "max_concurrent_trials": 3,
        "no_initial_point": True,
        "num_trials": 6,
        "proposal_batch_fraction": 1.0,
        "proposal_batch_size": 2,
        "rand_trials": 0,
        "rngseed": 123,
        "surrogate_fit_percentile": 0.95,
    }
    if algorithm == "hebo":
        settings = load_hebo_config()
        settings["gp"].update(num_epochs=2, optimizer="adam")
        settings["acquisition"].update(population=12, generations=3)
        values["hebo_settings"] = settings
    values.update(overrides)
    return SimpleNamespace(algorithm=algorithm, **values)


# Build the minimal trial parameter payload needed by the toy driver
def _param(tmp_path):
    return {
        "cdir": str(tmp_path),
        "cost": "loss",
        "datacards": [],
        "mc_steer": {},
        "max_t": 60.0,
        "pickle_dump": False,
        "plot": False,
        "run_name": "pytest-lbfgs",
    }


# Build the minimal trial parameter payload needed by the Ray figure driver
def _ray_param(tmp_path):
    return {
        "cdir": str(tmp_path),
        "cost": "loss",
        "datacards": [],
        "libdir": None,
        "mc_steer": {},
        "max_t": 60.0,
        "node_id": "pytest",
        "pickle_dump": False,
        "plot": True,
        "PYTHON_VERSION": "3.12",
        "run_name": "pytest-ray-default",
    }


# Start local Ray for tests, using the same startup patch module as icetune
def _ensure_ray(*, local_mode=True):
    from core.tune.backends import ray as icetune_ray  # noqa: F401

    if not ray.is_initialized():
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
                probe.bind(("127.0.0.1", 0))
        except OSError as exc:
            if exc.errno not in {errno.EACCES, errno.EAFNOSUPPORT, errno.EPERM}:
                raise
            pytest.skip(f"Local Ray requires TCP socket support: {exc}")
        ray.init(num_cpus=2, include_dashboard=False, local_mode=local_mode)


def test_cli_accepts_lbfgs_algorithm(monkeypatch):
    module = runpy.run_module("core.icetune")
    monkeypatch.setattr(sys, "argv", ["icetune", "--algorithm", "lbfgs", "--address", "local"])

    args = module["parse_arguments"]()

    assert args.algorithm == "lbfgs"
    assert args.history_interval_s == pytest.approx(60.0)
    assert args.lbfgs_settings["grad_step_rel"] == pytest.approx(1e-3)
    assert args.surrogate_fit_percentile == pytest.approx(1.0)
    assert args.wait_workers_timeout_s == pytest.approx(300.0)
    assert args.preflight is False


# Check initial-point validation reports the parameter value and fit bounds
def test_initial_config_bounds():
    bounds = {"coupling": tune.uniform(0.0, 2.0)}

    with pytest.raises(
        ValueError,
        match=r'Initial fit configuration is invalid: Parameter "coupling".*2\.5.*\[0\.0, 2\.0\]',
    ):
        parameter_space.validate_initial_config({"coupling": 2.5}, bounds)


# Check icetune promotes an invalid starting point to a permanent configuration error
def test_preflight_invalid_initial_config():
    module = runpy.run_module("core.icetune")
    args = SimpleNamespace(cdir=str(ROOT), no_initial_point=False)
    tunesetup = SimpleNamespace(
        aux_param_space={},
        param_space={"coupling": tune.uniform(0.0, 2.0)},
    )
    driver = SimpleNamespace(get_initial_param=lambda **kwargs: {"coupling": 2.5})

    with pytest.raises(
        module["PermanentConfigurationError"],
        match=r'"coupling".*2\.5.*\[0\.0, 2\.0\]',
    ):
        module["preflight_initial_config"](
            args=args,
            tunesetup=tunesetup,
            simdriver=driver,
            mc_steer={"tune_default": "TUNE0"},
        )


# Reject unsafe paths, invalid budgets, and unknown optimizers before initialization
@pytest.mark.parametrize(
    "flags",
    (
        ("--run_name", "../../escape"),
        ("--max_concurrent_trials", "0"),
        ("--num_trials", "20", "--rand_trials", "20"),
        ("--algorithm", "broken"),
        ("--ray_worker_cpu", "0"),
        ("--ray_worker_gpu", "-1"),
        ("--ray_worker_admission_timeout_s", "0"),
        ("--ray_worker_poll_interval_s", "0"),
        ("--ray_worker_failure_limit", "0"),
        ("--ray_worker_recovery_timeout_s", "0"),
        ("--ray_restore_retry_limit", "-1"),
        ("--ray_trial_retry_limit", "-1"),
        ("--ray_restore_retry_delay_s", "0"),
        ("--ray_resource_status_interval_s", "0"),
        ("--ray_status_interval_s", "0"),
    ),
)
def test_cli_invalid_campaign_controls(monkeypatch, flags):
    module = runpy.run_module("core.icetune")
    monkeypatch.setattr(sys, "argv", ["icetune", *flags])
    with pytest.raises(SystemExit):
        module["parse_arguments"]()


# Check local and external Ray initialization receive different resource options
@pytest.mark.parametrize("gpu", [0, 2])
def test_ray_init_args_runtime_specific(tmp_path, gpu):
    from core.tune.backends import ray as icetune_ray

    local_args = SimpleNamespace(
        address="local",
        ray_worker_cpu=3,
        ray_worker_gpu=gpu,
        ray_temp_dir=str(tmp_path / "ray"),
    )
    local = icetune_ray.ray_init_arguments(args=local_args, runtime_env="unused")
    assert local == {
        "include_dashboard": False,
        "_node_ip_address": "127.0.0.1",
        "num_cpus": 3,
        "num_gpus": gpu,
        "_temp_dir": str((tmp_path / "ray").resolve()),
    }

    external_args = SimpleNamespace(address="head.example:60010")
    external = icetune_ray.ray_init_arguments(
        args=external_args,
        runtime_env="runtime-env",
    )
    assert external == {
        "address": "head.example:60010",
        "runtime_env": "runtime-env",
    }


# Check Ray Tune stages enough placement groups to fill the effective trial budget
@pytest.mark.parametrize("initial", (None, "auto", " AUTO "))
def test_tune_pending_trials_follows_concurrency(monkeypatch, initial):
    if initial is None:
        monkeypatch.delenv(icetune_ray.RAY_TUNE_PENDING_ENV, raising=False)
    else:
        monkeypatch.setenv(icetune_ray.RAY_TUNE_PENDING_ENV, initial)
    args = SimpleNamespace(max_concurrent_trials=50, num_trials=30)

    assert icetune_ray.configure_tune_pending_trials(args) == 30
    assert os.environ[icetune_ray.RAY_TUNE_PENDING_ENV] == "30"


# Check an explicit positive Ray Tune staging limit remains authoritative
def test_tune_pending_trials_respects_explicit_limit(monkeypatch):
    monkeypatch.setenv(icetune_ray.RAY_TUNE_PENDING_ENV, "7")
    args = SimpleNamespace(max_concurrent_trials=50, num_trials=3000)

    assert icetune_ray.configure_tune_pending_trials(args) == 7
    assert os.environ[icetune_ray.RAY_TUNE_PENDING_ENV] == "7"


# Reject Ray Tune staging limits which would stall or misconfigure the controller
@pytest.mark.parametrize("value", ("0", "-1", "broken"))
def test_tune_pending_trials_invalid_limit(monkeypatch, value):
    monkeypatch.setenv(icetune_ray.RAY_TUNE_PENDING_ENV, value)
    args = SimpleNamespace(max_concurrent_trials=50, num_trials=3000)

    with pytest.raises(ValueError, match=icetune_ray.RAY_TUNE_PENDING_ENV):
        icetune_ray.configure_tune_pending_trials(args)


# Check driver status exposes only live per-node Ray resource markers
def test_ray_resource_state_live_node_markers(monkeypatch):
    monkeypatch.setattr(icetune_ray.ray, "is_initialized", lambda: True)
    monkeypatch.setattr(
        icetune_ray.ray,
        "available_resources",
        lambda: {"CPU": 8.0, "icetune_worker_252171_1": 1.0},
    )
    monkeypatch.setattr(
        icetune_ray.ray,
        "cluster_resources",
        lambda: {
            "CPU": 8.0,
            "icetune_head": 1.0,
            "icetune_worker_252171_1": 1.0,
        },
    )
    monkeypatch.setattr(
        icetune_ray.ray,
        "nodes",
        lambda: [
            {
                "Alive": True,
                "NodeID": "head-id",
                "NodeManagerAddress": "192.0.2.1",
                "Resources": {"icetune_head": 1.0},
            },
            {
                "Alive": True,
                "NodeID": "worker-id",
                "NodeManagerAddress": "192.0.2.2",
                "Resources": {"CPU": 8.0, "icetune_worker_252171_1": 1.0},
            },
            {
                "Alive": False,
                "NodeID": "dead-id",
                "NodeManagerAddress": "192.0.2.3",
                "Resources": {"icetune_worker_252171_2": 1.0},
            },
        ],
    )

    state = icetune_ray.ray_resource_state()

    assert state["cluster"]["icetune_head"] == 1.0
    assert state["nodes"]["worker-id"]["resources"] == {
        "CPU": 8.0,
        "icetune_worker_252171_1": 1.0,
    }
    assert "dead-id" not in state["nodes"]


# Check the zero CPU head cannot satisfy the first Ray trial capacity gate
def test_ray_trial_capacity_waits_worker_cpu(monkeypatch):
    nodes = iter(
        (
            [{"Alive": True, "Resources": {"icetune_head": 1.0}}],
            [
                {"Alive": True, "Resources": {"icetune_head": 1.0}},
                {"Alive": True, "Resources": {"CPU": 8.0}},
            ],
        )
    )
    resources = iter(
        (
            {"available": {"icetune_head": 1.0}, "cluster": {"icetune_head": 1.0}},
            {
                "available": {"CPU": 8.0, "icetune_head": 1.0},
                "cluster": {"CPU": 8.0, "icetune_head": 1.0},
            },
        )
    )
    monkeypatch.setattr(icetune_ray.ray, "nodes", lambda: next(nodes))
    monkeypatch.setattr(icetune_ray, "ray_resource_state", lambda: next(resources))
    monkeypatch.setattr(icetune_ray.time, "sleep", lambda seconds: None)
    args = SimpleNamespace(
        cpu_per_trial=8,
        gpu_per_trial=0,
        wait_workers=1,
        wait_workers_timeout_s=10,
    )

    assert icetune_ray.wait_ray_trial_capacity(args)["cluster"]["CPU"] == 8.0


# Check separate Ray allocations on one host cannot combine their CPUs for a trial
def test_ray_same_host_capacity():
    from ray.cluster_utils import Cluster

    ray.shutdown()
    cluster = Cluster()
    try:
        cluster.add_node(num_cpus=0, resources={"icetune_head": 1}, include_dashboard=False)
        for index in range(2):
            cluster.add_node(num_cpus=2, resources={f"icetune_worker_12_{index + 1}": 1})
        ray.init(address=cluster.address)
        cluster.wait_for_nodes()
        live = [node for node in ray.nodes() if node["Alive"]]
        assert len(live) == 3
        assert len({node["NodeManagerAddress"] for node in live}) == 1
        args = SimpleNamespace(cpu_per_trial=2, gpu_per_trial=0,
                               wait_workers=3, wait_workers_timeout_s=0)
        with pytest.raises(RuntimeError, match="workers=2"):
            icetune_ray.wait_ray_trial_capacity(args)
        args.wait_workers = 2
        assert icetune_ray.wait_ray_trial_capacity(args)["cluster"]["CPU"] == 4
        args.cpu_per_trial = 3
        with pytest.raises(RuntimeError, match="no complete trial resource bundle"):
            icetune_ray.wait_ray_trial_capacity(args)

        # Failure of one allocation must not mark the other allocation on this host
        workers = icetune_ray.RayWorkers(False)
        for node in live:
            for marker in node["Resources"]:
                if marker.startswith("icetune_worker_"):
                    workers.node_ids[marker] = node["NodeID"]
        workers.failure(f"{workers.node_ids['icetune_worker_12_1']}: trial worker lost")
        assert set(workers.failed) == {"icetune_worker_12_1"}
    finally:
        ray.shutdown()
        cluster.shutdown()


# Check uploaded Ray runtimes replace only checkout paths with relative paths
def test_portable_ray_environment_uploaded_checkout(tmp_path):
    checkout = tmp_path / "graniitti"
    external = tmp_path / "external"
    variables = {
        "LD_LIBRARY_PATH": str(external / "lib"),
        "PYTHONPATH": os.pathsep.join(
            (str(checkout / "python"), str(checkout), str(external / "python"))
        ),
    }

    output = icetune_ray.portable_ray_environment(variables, cdir=str(checkout))

    assert output["PYTHONPATH"].split(os.pathsep) == [
        "python",
        ".",
        str(external / "python"),
    ]
    assert output["LD_LIBRARY_PATH"] == str(external / "lib")
    assert output["RAY_CHDIR_TO_TRIAL_DIR"] == "0"


# Check uploaded Ray trials ignore a shared driver working directory
def test_ray_trial_uses_uploaded_worker_root(tmp_path, monkeypatch):
    head_root = tmp_path / "head-runtime"
    local_root = tmp_path / "ray-runtime"
    external = tmp_path / "external"
    calls = []
    driver = SimpleNamespace(
        runtime_environment=lambda **kwargs: calls.append(kwargs),
        shared_runtime_files=lambda param: [],
        dataset_paths=[str(head_root / "icepack" / "sample" / "dataset.json")],
        cuts=[
            [
                str(head_root / "icepack" / "sample" / "cuts.py"),
                str(external / "cuts.py"),
            ]
        ],
    )
    trainable = SimpleNamespace(logdir=str(tmp_path / "trial"))
    param = {
        "cdir": str(head_root),
        "datacards": [
            {"datacard": str(head_root / "icepack" / "sample" / "dataset.json")}
        ],
        "libdir": None,
        "PYTHON_VERSION": "3.12",
        "ray_upload_runtime": True,
        "max_t": 60.0,
    }
    monkeypatch.setattr(ray_runtime, "ray_worker_root", lambda: local_root)
    monkeypatch.setattr(
        icetune_ray.obs,
        "preflight_numba_runtime",
        lambda: calls.append("numba_preflight"),
    )

    icetune_ray.CFunc.setup(
        trainable,
        config={"x": 0.5},
        simdriver=driver,
        param=param,
        global_state=object(),
    )

    assert trainable.param["cdir"] == str(local_root)
    assert trainable.param["datacards"][0]["datacard"] == str(
        local_root / "icepack" / "sample" / "dataset.json"
    )
    assert driver.dataset_paths == [str(local_root / "icepack" / "sample" / "dataset.json")]
    assert driver.cuts == [
        [str(local_root / "icepack" / "sample" / "cuts.py"), str(external / "cuts.py")]
    ]
    assert calls[0]["cdir"] == str(local_root)
    assert calls[1] == "numba_preflight"
    assert param["cdir"] == str(head_root)
    assert param["datacards"][0]["datacard"] == str(
        head_root / "icepack" / "sample" / "dataset.json"
    )


# Check Ray storage accepts only the exact immutable campaign identity
def test_ray_campaign_identity_fences_restore_state(tmp_path):
    args = _search_args("hebo", run_name="run", no_initial_point=False)
    inputs = dict(
        args=args,
        tunesetup=_ToyTunesetup,
        initial_points={"x": 0.5, "y": -0.5},
        param=_param(tmp_path),
        experiment_dir=str(tmp_path / "run"),
    )
    stored = icetune_ray.ensure_ray_campaign_identity(**inputs)
    assert icetune_ray.ensure_ray_campaign_identity(**inputs) == stored

    changed = SimpleNamespace(**vars(args))
    changed.rngseed += 1
    with pytest.raises(RuntimeError, match="different physics or optimizer campaign"):
        icetune_ray.ensure_ray_campaign_identity(**{**inputs, "args": changed})

    legacy = tmp_path / "legacy"
    legacy.mkdir()
    (legacy / "tuner.pkl").write_bytes(b"ray-state")
    with pytest.raises(RuntimeError, match="no icetune campaign identity"):
        icetune_ray.ensure_ray_campaign_identity(**{**inputs, "experiment_dir": str(legacy)})


# Check a changed runtime resets only a provably empty Ray Tune checkpoint
def test_ray_empty_state_reset(tmp_path):
    args = _search_args("hebo", run_name="run", no_initial_point=False)
    param = _param(tmp_path)
    inputs = dict(
        args=args,
        tunesetup=_ToyTunesetup,
        initial_points={"x": 0.5, "y": -0.5},
        param=param,
        experiment_dir=str(tmp_path / "run"),
    )
    icetune_ray.ensure_ray_campaign_identity(**inputs, physics_fingerprint="old-runtime")
    state = Path(inputs["experiment_dir"])
    (state / "tuner.pkl").write_bytes(b"pending-ray-state")
    (state / "status.json").write_text(
        json.dumps(
            {
                "backend": "ray",
                "best_result": None,
                "completed": 0,
                "failed": 0,
            }
        ),
        encoding="utf-8",
    )
    (state / "experiment_state-000.json").write_text(
        json.dumps(
            {
                "trial_data": [
                    [
                        json.dumps({"status": "PENDING", "trial_id": "pending-1"}),
                        json.dumps(
                            {
                                "last_result": {
                                    "config": {"x": 0.5},
                                    "trial_id": "pending-1",
                                    "training_iteration": 0,
                                },
                                "num_failures": 0,
                            }
                        ),
                    ]
                ]
            }
        ),
        encoding="utf-8",
    )
    (state / "ray_init.json").write_text("current-init\n", encoding="utf-8")
    from core.tune.drivers.graniitti.driver import GraniittiDriver

    bank = GraniittiDriver().amplitude_directory(cdir=str(tmp_path), run_name="run")
    # Use this test's state root to exercise retention beside the Tune checkpoint
    bank = state / bank.relative_to(tmp_path / "runs/icetune/run")
    bank.mkdir(parents=True)
    (bank / "current.bin").write_bytes(b"current bank")
    (state / "submit/logs").mkdir(parents=True)
    (state / "submit/job.condor").write_text("queue 1\n")

    current = icetune_ray.ensure_ray_campaign_identity(
        **inputs,
        physics_fingerprint="current-runtime",
    )

    assert current["identity"]["physics_fingerprint"] == "current-runtime"
    assert not (state / "tuner.pkl").exists()
    assert (state / "ray_init.json").read_text(encoding="utf-8") == "current-init\n"
    assert list(tmp_path.glob("run._old-*"))
    assert (bank / "current.bin").read_bytes() == b"current bank"
    assert (state / "submit/job.condor").is_file()

    missing = tmp_path / "missing"
    missing.mkdir()
    (missing / "tuner.pkl").write_bytes(b"pending-ray-state")
    (missing / "status.json").write_text(
        json.dumps(
            {
                "backend": "ray",
                "best_result": None,
                "completed": 0,
                "failed": 0,
            }
        ),
        encoding="utf-8",
    )
    (missing / "experiment_state-000.json").write_text(
        json.dumps({"trial_data": []}),
        encoding="utf-8",
    )
    installed = icetune_ray.ensure_ray_campaign_identity(
        **{**inputs, "experiment_dir": str(missing)},
        physics_fingerprint="current-runtime",
    )
    assert installed["identity"]["physics_fingerprint"] == "current-runtime"


# Check a changed runtime cannot mix with any completed Ray evaluation
def test_ray_nonempty_state_mismatch(tmp_path):
    args = _search_args("hebo", run_name="run", no_initial_point=False)
    inputs = dict(
        args=args,
        tunesetup=_ToyTunesetup,
        initial_points={"x": 0.5, "y": -0.5},
        param=_param(tmp_path),
        experiment_dir=str(tmp_path / "run"),
    )
    icetune_ray.ensure_ray_campaign_identity(**inputs, physics_fingerprint="old-runtime")
    state = Path(inputs["experiment_dir"])
    (state / "tuner.pkl").write_bytes(b"completed-ray-state")
    (state / "status.json").write_text(
        json.dumps(
            {
                "backend": "ray",
                "best_result": None,
                "completed": 0,
                "failed": 0,
            }
        ),
        encoding="utf-8",
    )
    (state / "experiment_state-000.json").write_text(
        json.dumps(
            {
                "trial_data": [
                    [
                        json.dumps({"status": "TERMINATED"}),
                        json.dumps(
                            {
                                "last_result": {"loss": 1.0},
                                "num_failures": 0,
                            }
                        ),
                    ]
                ]
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="different physics or optimizer campaign"):
        icetune_ray.ensure_ray_campaign_identity(
            **inputs,
            physics_fingerprint="current-runtime",
        )


# Check restored campaigns replace worker sync with current head callbacks
@pytest.mark.parametrize("retry_limit", [0, 4])
def test_ray_restore_output_transport(retry_limit):
    restored_search = _search(_search_args("basic", num_trials=2))
    old_callback = object()
    internal = SimpleNamespace(
        _run_config=SimpleNamespace(
            callbacks=[old_callback],
            checkpoint_config=None,
            sync_config=None,
        ),
        _tune_config=SimpleNamespace(search_alg=restored_search),
    )
    tuner = SimpleNamespace(_local_tuner=internal)
    history = SimpleNamespace(search_alg=None)
    callbacks = [object(), history]

    output = icetune_ray.rebind_restored_ray_tuner(
        tuner=tuner,
        callbacks=callbacks,
        history_plot=history,
        retry_limit=retry_limit,
    )

    assert output is restored_search
    assert history.search_alg is restored_search
    assert internal._run_config.callbacks is callbacks
    assert not internal._run_config.sync_config.sync_artifacts
    assert not internal._run_config.checkpoint_config.checkpoint_at_end
    assert internal._run_config.failure_config.max_failures == retry_limit


# Check real history rendering leaves Tune free to accept results and propose trials
@pytest.mark.parametrize("publish_failure", [False, True])
def test_ray_history_plot_callback_completion_timer(tmp_path, publish_failure):
    from core.tune import core as icetune_main

    output = tmp_path / "figs" / "icetune" / "run" / "cost_evolution.png"
    args = _search_args("basic", num_trials=2)
    search = _search(args)
    search.completed_records = [
        {
            "completed_at_datetime": "2026-08-29 00:00:00 CEST",
            "completed_at_unix": 1787954400.0,
            "config": {"x": 0.4},
            "metrics": {"loss": 1.0},
            "node_id": "ray-worker",
            "search_payload": {"kind": "basic"},
            "trial_id": "trial-000000",
        }
    ]
    param = {
        **_param(tmp_path),
        "optimization": {"backend": "ray", "cost": "loss", "optimizer": "basic"},
        "parameter_space": [{"name": "x", "lower": 0.0, "upper": 1.0}],
        "parameter_topology": {"schema_version": 1, "groups": []},
        "plot_brand": "test",
        "simdriver": "GRANIITTI",
    }
    callback = icetune_ray.HistoryPlotCallback(
        experiment_dir=str(tmp_path / "ray" / "run"),
        cdir=str(tmp_path),
        cost="loss",
        run_name="run",
        interval_s=300.0,
        search_alg=search,
        param=param,
        num_trials=2,
        max_concurrent_trials=2,
    )
    token = icetune_main.acquire_publish_lock(cdir=str(tmp_path), run_name="run")
    try:
        # Hold real publication busy while requiring the completion callback to finish
        with ThreadPoolExecutor(max_workers=1) as calls:
            try:
                calls.submit(callback.on_trial_complete, iteration=1, trials=[], trial=object()).result(timeout=5.0)
                pending = callback.plot_future
                assert pending is not None and not pending.done()
                restored = pickle.loads(pickle.dumps(callback))
                assert restored.plot_pool is None and restored.plot_future is None

                search.next_index = 1
                assert search.suggest("ray-next") is not None
                search.on_trial_complete(
                    "ray-next",
                    result={
                        "loss": 0.5,
                        "completed_at_unix": 1787954401.0,
                        "completed_at_datetime": "2026-08-29 00:00:01 CEST",
                    },
                )
                callback.on_trial_complete(iteration=2, trials=[], trial=object())
                assert callback.plot_future is pending
                assert json.loads(callback.status_path.read_text())["completed"] == 2
                if publish_failure:
                    output.parent.mkdir(parents=True, exist_ok=True)
                    output.mkdir()
            finally:
                icetune_main.release_publish_lock(cdir=str(tmp_path), run_name="run", token=token)

        if publish_failure:
            with pytest.raises(IsADirectoryError):
                pending.result(timeout=60.0)
            assert callback.plot() is None
            assert callback.plot_pool is None
            output.rename(output.with_suffix(".png._old"))
        else:
            assert pending.result(timeout=60.0) == output
            assert callback.plot() == output
        assert callback.plot_future is None
        callback.on_step_end(iteration=3, trials=[])
        assert callback.plot_future is None
        assert callback.plot(force=True) == output

        history = json.loads(callback.history_path.read_text())
        assert [record["metrics"]["loss"] for record in history["trials"]] == [1.0, 0.5]
        assert json.loads(callback.status_path.read_text())["done"]
    finally:
        callback.close()
    assert callback.plot_pool is None


# Check Ray publishes readable live state before the first completion
def test_ray_running_queued_trials(tmp_path):
    args = _search_args("basic", max_concurrent_trials=3, num_trials=5)
    search = _search(args)
    assert search.suggest("ray-a") is not None
    assert search.suggest("ray-b") is not None
    param = {
        **_param(tmp_path),
        "optimization": {"backend": "ray", "cost": "loss", "optimizer": "basic"},
        "parameter_space": [{"name": "x", "lower": 0.0, "upper": 1.0}],
        "parameter_topology": {"schema_version": 1, "groups": []},
    }
    callback = icetune_ray.HistoryPlotCallback(
        experiment_dir=str(tmp_path / "ray" / "run"),
        cdir=str(tmp_path),
        cost="loss",
        run_name="run",
        interval_s=300.0,
        search_alg=search,
        param=param,
        num_trials=5,
        max_concurrent_trials=3,
    )
    trials = [
        SimpleNamespace(trial_id="ray-a", status="RUNNING"),
        SimpleNamespace(trial_id="ray-b", status="PENDING"),
    ]

    callback.on_step_end(iteration=1, trials=trials)

    history = json.loads(callback.history_path.read_text(encoding="utf-8"))
    status = json.loads(callback.status_path.read_text(encoding="utf-8"))
    assert history["trials"] == []
    assert status["pending"] == 3
    assert status["queued"] == 2
    assert status["running"] == 1


# Check Ray writes full failure diagnostics outside status and history
def test_ray_failure_log_separate_complete(tmp_path):
    args = _search_args("basic", num_trials=1)
    search = _search(args)
    assert search.suggest("ray-a") is not None
    search.on_trial_complete("ray-a", error=True)
    param = {
        **_param(tmp_path),
        "optimization": {"backend": "ray", "cost": "loss", "optimizer": "basic"},
        "parameter_space": [{"name": "x", "lower": 0.0, "upper": 1.0}],
        "parameter_topology": {"schema_version": 1, "groups": []},
    }
    callback = icetune_ray.HistoryPlotCallback(
        experiment_dir=str(tmp_path / "ray" / "run"),
        cdir=str(tmp_path),
        cost="loss",
        run_name="run",
        interval_s=300.0,
        search_alg=search,
        param=param,
        num_trials=1,
        max_concurrent_trials=1,
    )
    traceback_text = "Traceback (most recent call last):\nGRANIITTI event generation failed"
    trial = SimpleNamespace(
        config={"x": 0.5},
        get_error=lambda: traceback_text,
        num_failures=1,
        status="ERROR",
        trial_id="ray-a",
    )

    callback.on_trial_recover(iteration=1, trials=[trial], trial=trial)
    retry_path = callback.failure_dir / "trial-000000-retry-1.json"
    retry_bytes = retry_path.read_bytes()
    assert json.loads(retry_bytes)["error"] == traceback_text
    callback.on_trial_recover(iteration=1, trials=[trial], trial=trial)
    assert retry_path.read_bytes() == retry_bytes
    callback.on_trial_error(iteration=1, trials=[trial], trial=trial)

    failure = json.loads((callback.failure_dir / "trial-000000.json").read_text())
    assert failure["error"] == traceback_text
    assert failure["ray_trial_id"] == "ray-a"
    assert failure["trial_id"] == "trial-000000"
    assert traceback_text not in callback.status_path.read_text(encoding="utf-8")
    assert traceback_text not in callback.history_path.read_text(encoding="utf-8")


# Check the head output state receives worker figures and replica pickles
def test_ray_head_state_publishes_worker_outputs(tmp_path):
    from core.tune import core as icetune_main

    driver = _RayFigureDriver()
    param = _ray_param(tmp_path)
    param["pickle_dump"] = True
    experiment_dir = tmp_path / "ray" / param["run_name"]
    _ensure_ray(local_mode=False)
    head_actor = icetune_ray.GlobalState.remote(
        experiment_dir=str(experiment_dir), cdir=str(tmp_path),
        run_name=param["run_name"], cost=param["cost"],
    )
    callback = icetune_ray.TrialOutputCallback(
        experiment_dir=str(experiment_dir),
        param=param,
        global_state=head_actor,
    )

    # Render one worker result into its own Tune trial directory
    def trial_output(trial_id, x, *, initial, improving, commit=True):
        trial_dir = tmp_path / "trials" / trial_id
        output_param = {**param, "cdir": str(trial_dir)}
        outputs = icetune_main.evaluate_trial_outputs(
            config={"x": x},
            param=param,
            simdriver=driver,
            trial_id=trial_id,
        )
        icetune_ray._render_ray_trial_figures(
            outputs=outputs,
            param=output_param,
            simdriver=driver,
            initial=initial,
            improving=improving,
        )
        relative = Path("results", icetune_main.trial_pickle_filename(outputs))
        descriptor = icetune_main.maybe_dump_trial_payload(
            outputs=outputs,
            param=param,
            destination=trial_dir / relative,
        )
        descriptor["relative_path"] = relative.as_posix()
        metrics = {
            **outputs["metrics"],
            "is_initial": initial,
            "rendered_best": improving,
            "rendered_initial": initial,
            "trial_id": trial_id,
            "trial_pickle": descriptor,
        }
        transfer = icetune_ray.collect_ray_trial_outputs(
            trial_dir=trial_dir,
            run_name=param["run_name"],
            descriptor=descriptor,
            rendered_initial=initial,
            rendered_best=improving,
        )
        staged = ray.get(head_actor.publish_trial_outputs.remote(
            trial_id=trial_id,
            cost=outputs["metrics"][param["cost"]],
            transfer=transfer,
            descriptor=descriptor,
            rendered_initial=initial,
            rendered_best=improving,
        )
        )
        metrics["ray_output_token"] = staged["output_token"]
        if commit:
            callback.on_trial_complete(
                iteration=1,
                trials=[],
                trial=SimpleNamespace(last_result=metrics),
            )
        return metrics

    trial_output("ray-orphan", 0.4, initial=False, improving=True, commit=False)
    figure_root = tmp_path / "figs" / "icetune" / param["run_name"]
    assert not (figure_root / "summary.json").exists()
    assert not list((experiment_dir / "results").glob("TUNE_icetune_*.pkl"))
    callback.on_trial_error(
        iteration=1,
        trials=[],
        trial=SimpleNamespace(trial_id="ray-orphan", last_result={}),
    )

    best_metrics = trial_output(
        "ray-best", 0.2, initial=False, improving=True
    )
    trial_output(
        "ray-initial", 0.9, initial=True, improving=False
    )

    callback.flush(wait=True)
    assert json.loads((figure_root / "summary.json").read_text())["trial_id"] == "ray-best"
    assert json.loads((figure_root / "init" / "summary.json").read_text())[
        "trial_id"
    ] == "ray-initial"
    pickle_files = sorted((experiment_dir / "results").glob("TUNE_icetune_*.pkl"))
    assert len(pickle_files) == 2
    with open(pickle_files[0], "rb") as handle:
        assert pickle.load(handle)["replica_schema_version"] == 1
    callback.validate(
        best_result=SimpleNamespace(metrics=best_metrics),
        require_initial=True,
    )
    ray.kill(head_actor)
    ray.shutdown()


# Check a failed promotion retains its disk-backed transfer for an idempotent retry
def test_ray_head_output_commit_retries_local_stage(tmp_path, monkeypatch):
    experiment_dir = tmp_path / "ray" / "run"
    state = icetune_ray.BestState(
        experiment_dir=str(experiment_dir),
        cdir=str(tmp_path),
        run_name="run",
        cost="loss",
    )
    pickle_payload = b"complete worker pickle"
    descriptor = {
        "filename": "TUNE_icetune_trial-retry.pkl",
        "relative_path": "results/TUNE_icetune_trial-retry.pkl",
        "sha256": hashlib.sha256(pickle_payload).hexdigest(),
        "size": len(pickle_payload),
    }
    figure_source = tmp_path / "worker-figures"
    figure_source.mkdir()
    (figure_source / "summary.json").write_text(
        json.dumps({"metrics": {"loss": 1.0}, "trial_id": "trial-retry"}),
        encoding="utf-8",
    )
    figure_payload = icetune_ray.pack_ray_figure_tree(figure_source)
    staged = state.publish_trial_outputs(
        trial_id="trial-retry",
        cost=1.0,
        transfer={"figures": figure_payload, "pickle": pickle_payload},
        descriptor=descriptor,
        rendered_initial=False,
        rendered_best=True,
    )
    token = staged["output_token"]
    pending = state.pending_outputs["trial-retry"]
    stage = Path(pending["stage_dir"])
    assert stage.is_dir()
    assert stage.stat().st_mode & 0o777 == 0o700
    assert pending["pickle_path"] == str(stage / "trial.pkl")
    assert pending["figure_path"] == str(stage / "figure-tree" / "figures")
    assert all(not isinstance(value, bytes) for value in pending.values())
    repeated = state.publish_trial_outputs(
        trial_id="trial-retry",
        cost=1.0,
        transfer={"figures": figure_payload, "pickle": pickle_payload},
        descriptor=descriptor,
        rendered_initial=False,
        rendered_best=True,
    )
    assert repeated["output_token"] == token
    assert list(state.output_stage_root.glob("pending-*")) == [stage]

    # Simulate one canonical figure publication failure after the pickle succeeds
    def fail_figure_publish(**kwargs):
        del kwargs
        raise OSError("temporary figure failure")

    publish_figures = state._publish_figures
    monkeypatch.setattr(state, "_publish_figures", fail_figure_publish)
    with pytest.raises(OSError, match="temporary figure failure"):
        state.commit_trial_outputs("trial-retry", token)
    assert state.pending_outputs["trial-retry"]["token"] == token
    assert stage.is_dir()
    assert (experiment_dir / "results" / descriptor["filename"]).is_file()

    monkeypatch.setattr(state, "_publish_figures", publish_figures)
    published = state.commit_trial_outputs("trial-retry", token)
    assert published == {"published_best": True, "published_initial": False}
    assert not stage.exists()
    assert "trial-retry" not in state.pending_outputs
    assert json.loads(
        (tmp_path / "figs" / "icetune" / "run" / "summary.json").read_text(
            encoding="utf-8"
        )
    )["trial_id"] == "trial-retry"
    assert state.commit_trial_outputs("trial-retry", token) == published


# Check a reconstructed head actor restores pending stages and committed receipts
def test_ray_head_output_recovery(tmp_path):
    experiment_dir = tmp_path / "ray" / "run"
    settings = {
        "experiment_dir": str(experiment_dir),
        "cdir": str(tmp_path),
        "run_name": "run",
        "cost": "loss",
    }
    pickle_payload = b"restart-safe worker pickle"
    descriptor = {
        "filename": "TUNE_icetune_trial-restart.pkl",
        "relative_path": "results/TUNE_icetune_trial-restart.pkl",
        "sha256": hashlib.sha256(pickle_payload).hexdigest(),
        "size": len(pickle_payload),
    }
    first = icetune_ray.BestState(**settings)
    staged = first.publish_trial_outputs(
        trial_id="trial-restart",
        cost=2.0,
        transfer={"figures": None, "pickle": pickle_payload},
        descriptor=descriptor,
        rendered_initial=False,
        rendered_best=False,
    )
    token = staged["output_token"]

    reconstructed = icetune_ray.BestState(**settings)
    assert reconstructed.pending_outputs["trial-restart"]["token"] == token
    published = reconstructed.commit_trial_outputs("trial-restart", token)
    assert published == {"published_best": False, "published_initial": False}
    target = experiment_dir / "results" / descriptor["filename"]
    assert target.read_bytes() == pickle_payload

    after_commit = icetune_ray.BestState(**settings)
    assert after_commit.commit_trial_outputs("trial-restart", token) == published
    assert after_commit.pending_outputs == {}


# Check the cluster preflight reports the exact inconsistent node fields
def test_ray_node_config_mismatch():
    from core.tune.backends import ray as icetune_ray

    valid = {
        "cdir": True,
        "entrypoint": True,
        "hostname": "worker-a",
        "libdir": True,
        "python_version": "3.12",
        "ray_version": "2.54.1",
        "storage": True,
        "storage_writable": True,
    }
    icetune_ray.validate_ray_node_contracts(
        [valid],
        expected_python="3.12",
        expected_ray="2.54.1",
    )

    invalid = dict(valid, hostname="worker-b", storage_writable=False, ray_version="2.53.0")
    with pytest.raises(RuntimeError, match="worker-b: storage_writable, ray_version=2.53.0"):
        icetune_ray.validate_ray_node_contracts(
            [invalid],
            expected_python="3.12",
            expected_ray="2.54.1",
        )


# Verify every Ray search algorithm uses the shared asynchronous searcher
@pytest.mark.parametrize(
    "algorithm",
    [
        "basic",
        "hyperopt",
        "optuna",
        "icebo",
        "hebo",
        "bayesopt",
        "ax",
        "scikit",
    ],
)
def test_ray_searchers_shared_async_searcher(algorithm):
    from core.tune.backends import ray as icetune_ray

    args = _search_args(algorithm, surrogate_fit_percentile=0.5)

    search_alg = icetune_ray.set_search_algo(
        args=args,
        tunesetup=_ToyTunesetup(),
        initial_points=None,
    )

    assert isinstance(search_alg, icetune_ray.PercentileSearch)


# Verify random warmup fills every free Ray slot without feedback
def test_ray_random_warmup_fills_all_available_slots():
    args = _search_args(
        "hebo",
        max_concurrent_trials=8,
        num_trials=20,
        rand_trials=20,
    )
    search = _search(args)

    configs = [search.suggest(f"ray-{index}") for index in range(8)]

    assert all(config is not None for config in configs)
    assert len(search.live_records) == 8
    assert all(
        record["search_payload"]["kind"] == "random"
        for record in search.live_records.values()
    )
    assert search.suggest("ray-overflow") is None


# Preserve physical parameters and likelihood uncertainty through Ray result flattening
def test_ray_search_objective_error(tmp_path):
    from core.tune import history as icetune_history
    from core.tune.search import SearchState
    from ray.tune.experiment import Trial
    from ray.tune.utils import flatten_dict

    args = _search_args("basic")
    search = _search(args)
    search.live_records["ray-feedback"] = {
        "config": {"x": 0.4},
        "search_payload": {"kind": "basic"},
        "trial_id": "trial-000000",
    }
    likelihood = {
        "objective": {"name": "loss"},
        "mc_uncertainty": {"sigma_objective": 0.125},
        "mc_statistics": [{"dataset": 0, "subset": 0, "observable": "M",
                           "mc_empty_bin_count": 2, "mc_empty_bins": [1, 4], "mc_valid_bin_count": 8}],
    }
    base = "CON_GP|990:[211,211]/opposite:helicity"
    card_config = {
        "schema_version": 1,
        "parameters": {f"{base}(0,0,0)@MAG": 0.8},
        "tables": {base: {"basis": "helicity", "rows": [[0, 0, 0, 0.8, 0.0]]}},
    }
    result = {
        "card_config": card_config,
        "completed_at_datetime": "2026-08-29T00:00:00Z",
        "completed_at_unix": 1787961600.0,
        "likelihood": likelihood,
        "loss": 1.0,
        "optimizer_metrics": {"loss": 1.0, "objective_error": 0.2},
        "theta": {"names": ["x"], "values": [0.4]},
    }
    callback = icetune_ray.HistoryPlotCallback(
        experiment_dir=str(tmp_path / "ray"),
        cdir=str(tmp_path),
        cost="loss",
        run_name="run",
        interval_s=300.0,
        search_alg=search,
        param={
            "optimization": {"backend": "ray", "cost": "loss"},
            "parameter_space": [{"name": "x", "lower": 0.0, "upper": 1.0}],
            "parameter_topology": {"schema_version": 1, "groups": []},
        },
        num_trials=args.num_trials,
        max_concurrent_trials=args.max_concurrent_trials,
    )
    tune.register_trainable("icetune-history", icetune_ray.CFunc)
    trial = Trial("icetune-history", trial_id="ray-feedback")
    callback.on_trial_result(iteration=1, trials=[trial], trial=trial, result=result)
    search.on_trial_complete(trial.trial_id, result=flatten_dict(result))

    record = search.completed_records[0]
    assert record["card_config"] == card_config
    assert record["metrics"] == {"loss": 1.0, "objective_error": 0.2}
    assert record["likelihood"] == likelihood
    assert record["theta"] == result["theta"]
    assert "result" not in record
    state = SearchState(args=args, bounds=search.bounds, initial_points=None)
    assert state._record_cost_error(record) == pytest.approx(0.2)
    likelihood_only = {**record, "metrics": {"loss": 1.0}}
    assert state._record_cost_error(likelihood_only) == pytest.approx(0.125)

    callback._write_state(force_history=True)
    history, _, points = icetune_history.load_history(callback.history_path)
    assert history["trials"][0]["likelihood"]["mc_statistics"] == likelihood["mc_statistics"]
    _, names, physical = icetune_history.parameter_evolution(history, points, physical=True)
    _, raw_names, raw = icetune_history.parameter_evolution(history, points, physical=False)
    assert f"{base}(0,0,0)@MAG" in names
    assert physical[f"{base}(0,0,0)@MAG"] == pytest.approx([0.8])
    assert names == [f"{base}(0,0,0)@MAG"]
    assert raw_names == ["x"]
    assert raw["x"] == pytest.approx([0.4])


# Verify Ray refills from HEBO without waiting for a feedback block
def test_ray_hebo_refills_while_next_fit_runs(monkeypatch):
    fitted = []
    args = _search_args("hebo")
    search = _search(args)
    search.completed_records = [
        {
            "config": {"x": 0.5},
            "metrics": {"loss": 1.0},
            "search_payload": {"kind": "random"},
            "trial_id": "trial-000000",
        }
    ]
    search.next_index = 1
    calls = []

    # Compute an already completed background fit with unique proposals
    def submit(request):
        calls.append((request["next_index"], request["count"]))
        future = Future()
        future.set_result(icetune_ray._fit_search_state(request))
        fitted.extend(future.result()["proposals"])
        return future

    monkeypatch.setattr(search, "_submit_fit", submit)

    assert search.suggest("ray-a") == fitted[0][0]
    assert search.suggest("ray-b") == fitted[1][0]
    assert search.suggest("ray-c") == fitted[2][0]
    assert search.suggest("ray-blocked") is None
    search.on_trial_complete("ray-a", result={"loss": 1.0})
    assert search.suggest("ray-refill") == fitted[3][0]
    assert "ray-c" in search.live_records
    assert calls == [(1, 3), (4, 2)]


# Verify one completed outcome does not discard unissued HEBO recommendations
def test_ray_hebo_prefetch_after_feedback(monkeypatch):
    fitted = []
    args = _search_args("hebo", max_concurrent_trials=3, proposal_batch_size=2)
    search = _search(args)
    search.completed_records = [
        {
            "config": {"x": 0.5},
            "metrics": {"loss": 1.0},
            "search_payload": {"kind": "random"},
            "trial_id": "trial-000000",
        }
    ]
    search.next_index = 1
    calls = []

    # Compute one ready proposal reserve
    def submit(request):
        calls.append((request["next_index"], request["count"]))
        future = Future()
        future.set_result(icetune_ray._fit_search_state(request))
        fitted.extend(future.result()["proposals"])
        return future

    monkeypatch.setattr(search, "_submit_fit", submit)

    assert search.suggest("ray-first") == fitted[0][0]
    search.on_trial_complete("ray-first", result={"loss": 0.5})
    assert search.suggest("ray-next") == fitted[1][0]
    assert calls == [(1, 3)]


# Verify Ray prefetch cannot cross from independent warmup into HEBO
def test_ray_hebo_prefetch_waits_warmup_boundary():
    args = _search_args("hebo", num_trials=10, proposal_batch_size=3, rand_trials=5)
    search = _search(args)
    search.next_index = 4
    search.live_records = {
        "ray-a": {"config": {"x": 0.0}, "trial_id": "trial-000002"},
        "ray-b": {"config": {"x": 1.0}, "trial_id": "trial-000003"},
    }
    proposals = search._new_proposal_batch()

    assert len(proposals) == 1
    assert proposals[0][1]["kind"] == "random"
    assert search.fit_future is None


# Verify Ray starts HEBO with unfinished random configurations marked pending
def test_ray_hebo_starts_with_pending_warmup(monkeypatch):
    fitted = []
    args = _search_args("hebo", num_trials=10, rand_trials=5)
    search = _search(args)
    search.next_index = 5
    search.completed_records = [
        {
            "config": {"x": 0.1 * index},
            "metrics": {"loss": 1.0},
            "search_payload": {"kind": "random"},
            "trial_id": f"trial-{index:06d}",
        }
        for index in range(3)
    ]
    search.live_records = {
        f"ray-{index}": {
            "config": {"x": 0.1 * index},
            "search_payload": {"kind": "random"},
            "trial_id": f"trial-{index:06d}",
        }
        for index in range(3, 5)
    }
    calls = []

    # Compute pending aware adaptive proposals at the issued warmup boundary
    def submit(request):
        calls.append(
            (
                request["next_index"],
                request["count"],
                copy.deepcopy(request["pending_configs"]),
            )
        )
        future = Future()
        future.set_result(icetune_ray._fit_search_state(request))
        fitted.extend(future.result()["proposals"])
        return future

    monkeypatch.setattr(search, "_submit_fit", submit)

    proposals = search._new_proposal_batch()

    assert proposals
    assert all(payload["kind"] == "hebo" for _, payload in proposals)
    assert calls[0][:2] == (5, 3)
    assert calls[0][2] == [
        {"x": pytest.approx(0.3)},
        {"x": pytest.approx(0.4)},
    ]
    assert set(search.live_records) == {"ray-3", "ray-4"}


# Verify cached HEBO proposals refill Ray while the next fit is unfinished
def test_ray_reserve_refills_during_background_fit(monkeypatch):
    fitted = []
    args = _search_args("hebo", num_trials=10, rand_trials=5)
    search = _search(args)
    search.next_index = 5
    search.completed_records = [
        {
            "config": {"x": 0.1 * index},
            "metrics": {"loss": 1.0},
            "search_payload": {"kind": "random"},
            "trial_id": f"trial-{index:06d}",
        }
        for index in range(3)
    ]
    search.live_records = {
        f"ray-{index}": {
            "config": {"x": 0.1 * index},
            "search_payload": {"kind": "random"},
            "trial_id": f"trial-{index:06d}",
        }
        for index in range(3, 5)
    }
    futures = []

    # Complete the first reserve fit and keep the next fit pending
    def submit(request):
        future = Future()
        if not futures:
            future.set_result(icetune_ray._fit_search_state(request))
            fitted.extend(future.result()["proposals"])
        futures.append(future)
        return future

    monkeypatch.setattr(search, "_submit_fit", submit)

    assert search.suggest("ray-5") == fitted[0][0]
    search.on_trial_complete("ray-3", result={"loss": 1.0})
    assert search.suggest("ray-6") == fitted[1][0]
    search.on_trial_complete("ray-4", result={"loss": 1.0})
    assert search.fit_future is futures[1]
    assert not search.fit_future.done()
    assert search.suggest("ray-7") == fitted[2][0]


# Verify Ray replaces missing warmup outcomes before entering HEBO
def test_ray_hebo_missing_warmup_is_reissued():
    args = _search_args("hebo", rand_trials=2)
    search = _search(args)
    search.next_index = 2

    first = search.suggest("ray-random-0")
    assert first == search._uniform_config(2)
    assert search.live_records["ray-random-0"]["search_payload"]["kind"] == "random"
    second = search.suggest("ray-random-1")
    assert second == search._uniform_config(3)
    assert search.live_records["ray-random-1"]["search_payload"]["kind"] == "random"
    assert search.suggest("ray-wait") is None

    search.on_trial_complete("ray-random-0", error=True)
    assert search.suggest("ray-still-waiting") is None


# Verify a HEBO fit runs while Ray trials remain live
def test_ray_hebo_fit_runs_while_trials_live(monkeypatch):
    fitted = []
    args = _search_args("hebo", max_concurrent_trials=3, num_trials=8)
    search = _search(args)
    search.completed_records = [
        {
            "config": {"x": 0.5},
            "metrics": {"loss": 1.0},
            "search_payload": {"kind": "random"},
            "trial_id": "trial-000000",
        }
    ]
    search.next_index = 1
    requests = []

    # Complete the first fit and leave the refill fit running
    def submit(request):
        requests.append(copy.deepcopy(request))
        future = Future()
        if len(requests) == 1:
            future.set_result(icetune_ray._fit_search_state(request))
            fitted.extend(future.result()["proposals"])
        return future

    monkeypatch.setattr(search, "_submit_fit", submit)
    for trial_id in ("ray-a", "ray-b", "ray-c"):
        assert search.suggest(trial_id) is not None
    assert search.fit_future is not None
    assert not search.fit_future.done()
    assert len(search.live_records) == 3
    search.fit_future.set_result(icetune_ray._fit_search_state(requests[-1]))
    search.on_trial_complete("ray-a", error=True)
    assert search.suggest("ray-next") is not None
    assert search.live_records["ray-next"]["search_payload"]["kind"] == "hebo"
    assert [(request["next_index"], request["count"]) for request in requests] == [
        (1, 3),
        (4, 3),
    ]


# Verify an active HEBO fit restarts safely after a Ray searcher checkpoint
def test_ray_hebo_background_fit_checkpoint(tmp_path, monkeypatch):
    fitted = []
    args = _search_args("hebo", max_concurrent_trials=3, num_trials=8)
    search = _search(args)
    search.completed_records = [
        {
            "config": {"x": 0.5},
            "metrics": {"loss": 1.0},
            "search_payload": {"kind": "random"},
            "trial_id": "trial-000000",
        }
    ]
    search.next_index = 1
    calls = []

    # Complete the initial fit and leave its replacement pending
    def submit(request):
        calls.append((request["next_index"], request["count"]))
        future = Future()
        if len(calls) == 1:
            future.set_result(icetune_ray._fit_search_state(request))
            fitted.extend(future.result()["proposals"])
        return future

    monkeypatch.setattr(search, "_submit_fit", submit)

    for trial_id in ("ray-a", "ray-b", "ray-c"):
        assert search.suggest(trial_id) is not None
    assert search.fit_request is not None
    del search.__dict__["_submit_fit"]
    checkpoint = tmp_path / "search.pkl"
    search.save(str(checkpoint))

    restored = _search(args)
    restored.restore(str(checkpoint))

    # Restart the serialized request and return its proposals
    def restore_submit(request):
        calls.append((request["next_index"], request["count"]))
        future = Future()
        future.set_result(icetune_ray._fit_search_state(request))
        fitted.extend(future.result()["proposals"])
        return future

    monkeypatch.setattr(restored, "_submit_fit", restore_submit)
    restored.on_trial_complete("ray-a", result={"loss": 0.5})
    assert restored.suggest("ray-next") == fitted[3][0]
    assert calls == [(1, 3), (4, 3), (4, 3)]


# Verify the custom Ray searcher retries a failed initial point
def test_ray_searcher_retries_failed_initial_point():
    args = SimpleNamespace(
        algorithm="bayesopt",
        cost="loss",
        rand_trials=10,
        rngseed=123,
        surrogate_fit_percentile=0.95,
    )
    search = _search(args, {"x": 0.9})

    first = search.suggest("ray-first")
    second = search.suggest("ray-second")
    search.on_trial_complete("ray-first", error=True)
    retry = search.suggest("ray-retry")
    search.on_trial_complete("ray-retry", result={"loss": 1.0})
    following = search.suggest("ray-following")

    assert first == {"x": pytest.approx(0.9)}
    assert second is None
    assert retry == first
    assert following != first


# Verify every Ray optimizer refills a terminal slot with stragglers still live
@pytest.mark.parametrize(
    "algorithm",
    ["basic", "icebo", "hebo", "hyperopt", "optuna", "bayesopt", "ax", "scikit"],
)
def test_ray_searcher_refills_stragglers(algorithm):
    args = _search_args(algorithm, max_concurrent_trials=2, rand_trials=5)
    search = _search(args)
    first = search.suggest("ray-a")
    second = search.suggest("ray-b")
    assert first is not None and second is not None and first != second
    assert search.suggest("ray-blocked") is None
    search.on_trial_complete("ray-a", result={"loss": (first["x"] - 0.37)**2})
    third = search.suggest("ray-refill")
    assert third is not None and third != second
    assert "ray-b" in search.live_records
    search.on_trial_complete("ray-b", error=True)
    assert search.suggest("ray-next") is not None
    assert len(search.live_records) == 2


# Reproduce native proposals with the submitted seed
@pytest.mark.parametrize("algorithm", ["bayesopt", "ax", "scikit"])
def test_native_optimizer_seed(algorithm):
    # Compute proposals from the same completed scalar objective measurements
    def proposals(seed):
        search = _search(_search_args(algorithm, rngseed=seed, surrogate_fit_percentile=1.0))
        search.completed_records = [
            {"trial_id": f"trial-{index}", "config": {"x": x}, "metrics": {"loss": (x - 0.37) ** 2}}
            for index, x in enumerate((0.1, 0.5, 0.9))
        ]
        return [config["x"] for config, _ in search._ask_native_batch(2)]

    first = proposals(713)
    np.testing.assert_allclose(first, proposals(713), rtol=1e-10, atol=1e-12)
    assert not np.allclose(first, proposals(719), rtol=1e-10, atol=1e-12)


# Run the feedback handoff through each shared Ray optimizer
@pytest.mark.parametrize("algorithm", ["hebo", "optuna", "hyperopt", "icebo", "bayesopt", "ax", "scikit"])
def test_real_ray_optimizer_feedback_handoff(algorithm):
    args = _search_args(algorithm, rand_trials=2, surrogate_fit_percentile=1.0)
    if algorithm == "icebo":
        settings = load_icebo_config(resource("tune/settings/icebo.json"))
        args.icebo_settings = override_icebo_config(
            settings,
            {
                "acquisition.gradient_steps": 2,
                "acquisition.mc_samples": 16,
                "acquisition.optimization_batch_size": 16,
                "acquisition.raw_samples": 32,
                "acquisition.restarts": 2,
                "gp.fit_restarts": 1,
                "gp.fit_steps": 3,
            },
        )
    search = _search(args)

    # Finish one live trial through the deterministic toy objective
    def finish(trial_id):
        x = float(search.live_records[trial_id]["config"]["x"])
        search.on_trial_complete(trial_id, result={"loss": (x - 0.37) ** 2})

    # Wait only in the test harness for a real background optimizer fit
    def suggest_when_ready(trial_id):
        while True:
            suggestion = search.suggest(trial_id)
            if suggestion is not None:
                return suggestion
            assert search.fit_future is not None
            search.fit_future.result(timeout=300.0)

    assert search.suggest("ray-warm-0") is not None
    assert search.suggest("ray-warm-1") is not None
    assert search.suggest("ray-warm-2") is None
    assert search.suggest("ray-full") is None

    finish("ray-warm-1")
    assert suggest_when_ready("ray-adaptive-pending") is not None
    payload = search.live_records["ray-adaptive-pending"]["search_payload"]
    assert payload["kind"] == algorithm
    assert "ray-warm-0" in search.live_records
    if algorithm == "icebo":
        assert payload["acquisition"] == "qlognei"
    finish("ray-warm-0")

    next_trial = 0
    while search.next_index < args.num_trials or search.live_records:
        while search.next_index < args.num_trials:
            trial_id = f"ray-adaptive-{next_trial}"
            suggestion = search.suggest(trial_id)
            if suggestion is None:
                break
            next_trial += 1
        if not search.live_records and search.next_index < args.num_trials:
            trial_id = f"ray-adaptive-{next_trial}"
            assert suggest_when_ready(trial_id) is not None
            next_trial += 1
        assert search.live_records
        finish(next(iter(search.live_records)))
    assert len(search.completed_records) == args.num_trials
    assert search.suggest("ray-finished") == icetune_ray.Searcher.FINISHED


def test_parameter_unpack_is_deterministic():
    bounds = icetune_lbfgs.normalize_param_space(_ToyTunesetup.param_space)
    names = sorted(bounds)
    decoded = icetune_lbfgs.unpack_config([0.25, -0.5], names)

    assert names == ["x", "y"]
    assert list(decoded) == ["x", "y"]
    assert decoded == pytest.approx({"x": 0.25, "y": -0.5})


def test_uniform_bounds_extracted_ray_tune_space():
    bounds = icetune_lbfgs.normalize_param_space(_ToyTunesetup.param_space)

    assert bounds["x"]["lower"] == pytest.approx(-2.0)
    assert bounds["x"]["upper"] == pytest.approx(2.0)
    assert bounds["y"]["lower"] == pytest.approx(-2.0)
    assert bounds["y"]["upper"] == pytest.approx(2.0)


def test_unsupported_space_fails_clearly():
    with pytest.raises(ValueError, match="continuous uniform"):
        icetune_lbfgs.normalize_param_space(_ChoiceTunesetup.param_space)


def test_finite_diff_gradient_matches_quadratic():
    bounds = {
        "x": {"lower": -2.0, "upper": 2.0},
        "y": {"lower": -2.0, "upper": 2.0},
    }
    names = ["x", "y"]
    x0 = np.array([0.3, -0.4])
    stencils = icetune_lbfgs.build_gradient_stencils(
        x0,
        bounds,
        names,
        rel_step=1e-4,
        abs_step=1e-6,
    )

    values_by_point = {}
    for stencil in stencils:
        for point in stencil.points:
            shifted = np.array(x0, copy=True)
            shifted[stencil.index] = point
            values_by_point[(stencil.index, point)] = (shifted[0] - 0.4) ** 2 + (
                shifted[1] + 0.2
            ) ** 2

    grad = icetune_lbfgs.finite_difference_gradient_from_values(stencils, values_by_point)

    assert grad == pytest.approx(np.array([-0.2, -0.4]), abs=1e-8)


def test_bounded_stencil_stays_inside_limits():
    bounds = {"x": {"lower": 0.0, "upper": 1.0}}
    stencils = icetune_lbfgs.build_gradient_stencils(
        np.array([0.0]),
        bounds,
        ["x"],
        rel_step=1e-3,
        abs_step=1e-6,
    )

    assert stencils[0].points[0] >= 0.0
    assert all(0.0 <= point <= 1.0 for point in stencils[0].points)
    assert len(stencils[0].points) == 3


# Verify equal bounded probes share one Ray task while retaining request diagnostics
def test_lbfgs_batch_deduplicates_fresh_configs(tmp_path, monkeypatch):
    calls = []

    # Compute one immediate synthetic simulator result
    def evaluate(config, param, simdriver, trial_id, seed):
        del param, simdriver, seed
        calls.append(float(config["x"]))
        cost = float(config["x"]) ** 2
        return {
            "config": dict(config),
            "cost": cost,
            "metrics": {"loss": cost},
            "success": True,
            **icetune_lbfgs.trial_times(time.time(), time.time()),
            "trial_id": trial_id,
        }

    remote = SimpleNamespace(options=lambda **_: SimpleNamespace(remote=evaluate))
    monkeypatch.setattr(icetune_lbfgs, "_evaluate_config_remote", remote)
    monkeypatch.setattr(icetune_lbfgs.ray, "get", lambda values: values)
    objective = icetune_lbfgs.DistributedLBFGSObjective(
        args=_args(),
        bounds={"x": {"lower": 0.0, "upper": 1.0}},
        names=["x"],
        param=_param(tmp_path),
        simdriver=object(),
        output_dir=str(tmp_path / "lbfgs"),
    )
    requests = [
        {"config": {"x": value}, "label": label, "seed": 7, "trial_id": label}
        for value, label in ((0.25, "first"), (0.75, "second"), (0.25, "duplicate"))
    ]

    records = objective._evaluate_many(requests)
    assert calls == [0.25, 0.75]
    assert [record["label"] for record in records] == ["first", "second", "duplicate"]
    assert [record["cached"] for record in records] == [False, False, True]
    assert Path(objective.publish_history()).is_file()
    history = json.loads(Path(objective.history_path).read_text())
    assert len(history["trials"]) == 2
    assert [record["label"] for record in history["requests"]] == [
        "first",
        "second",
        "duplicate",
    ]


@pytest.mark.parametrize("budget", [1, 80])
def test_ray_lbfgs_quadratic(tmp_path, budget):
    _ensure_ray(local_mode=False)
    try:
        summary = icetune_lbfgs.run_lbfgs_backend(
            args=_args(num_trials=budget),
            tunesetup=_ToyTunesetup,
            simdriver=_ToyDriver(),
            initial_points={"x": 1.5, "y": -1.5},
            param={**_param(tmp_path), "plot": True, "head_pid": os.getpid()},
            start_time=time.time(),
        )
    finally:
        ray.shutdown()

    assert summary["num_unique_points"] <= budget
    if budget == 1:
        assert not summary["diagnostics"]["success"]
        return
    assert summary["best_config"]["x"] == pytest.approx(0.4, abs=5e-2)
    assert summary["best_config"]["y"] == pytest.approx(-0.2, abs=5e-2)
    assert summary["best_metrics"]["loss"] < 1e-3
    assert summary["best_fit"]["schema_version"] == 1
    assert summary["best_fit"]["source"] == "realized"
    assert summary["best_fit"]["objective"]["name"] == "loss"
    assert summary["best_fit"]["parameters"]["x"]["uncertainty"] is None
    assert (tmp_path / "runs" / "icetune" / "pytest-lbfgs" / "summary.json").is_file()
    history_path = tmp_path / "runs" / "icetune" / "pytest-lbfgs" / "history.json"
    assert Path(summary["history_path"]) == history_path
    assert Path(summary["history_plot_path"]).is_file()
    assert summary["published_summary"]["config"] == summary["best_config"]
    assert (tmp_path / "figs" / "icetune" / "pytest-lbfgs" / "fit.txt").is_file()


# Check final figures cross the worker boundary and uploaded paths are rebased
def test_lbfgs_worker_figure_transfer(tmp_path, monkeypatch):
    head = tmp_path / "head"
    worker = tmp_path / "worker"
    worker.mkdir()
    param = {**_param(head), "ray_upload_runtime": True}
    driver = _RayFigureDriver()
    driver.dataset_paths = [str(head / "datasets")]
    driver.cuts = {"path": str(head / "cuts")}
    calls = []
    render = icetune_lbfgs._render_best_remote._function

    class Remote:
        # Record resources reserved for the final simulator call
        def options(self, **kwargs):
            calls.append(kwargs)
            return self

        # Run the worker body against a separate runtime directory
        def remote(self, *args):
            return render(*args)

    monkeypatch.setattr(icetune_lbfgs, "_render_best_remote", Remote())
    monkeypatch.setattr(icetune_lbfgs.ray, "get", lambda value: value)
    monkeypatch.setattr(ray_runtime, "ray_worker_root", lambda: worker)
    monkeypatch.setattr(icetune_lbfgs.icetune_main, "publish_trial_rerun",
                        lambda **kwargs: pytest.fail("The head must not rerun the simulator"))
    summary = icetune_lbfgs._publish_best(
        args=_args(), config={"x": 0.4}, param=param, simdriver=driver, metrics={"loss": 0.0}
    )
    assert calls == [{"num_cpus": 1, "num_gpus": 0}]
    assert summary["config"] == {"x": 0.4}
    assert driver.dataset_paths == [str(worker / "datasets")]
    assert driver.cuts == {"path": str(worker / "cuts")}
    assert param["cdir"] == str(head)
    target = head / "figs" / "icetune" / param["run_name"]
    assert (target / "output.txt").read_text() == "lbfgs-final-best"
    assert json.loads((target / "summary.json").read_text())["config"] == {"x": 0.4}


# Check worker likelihood and parameter identity survive head history publication
@pytest.mark.parametrize("pickle_dump", [False, True])
def test_lbfgs_worker_history_likelihood(tmp_path, monkeypatch, pickle_dump):
    class LikelihoodDriver(_ToyDriver):
        # Supply the canonical Gaussian objective produced by a histogram driver
        def evaluate_trial_outputs(self, **kwargs):
            outputs = super().evaluate_trial_outputs(**kwargs)
            outputs["likelihood"] = icetune_likelihood.build_likelihood_payload(
                cost=[[{"obs": outputs["metrics"]["loss"]}]], ndf=[[{"obs": 1}]],
                weight=[[{"obs": 1.0}]], valid=[[{"obs": 1}]],
            )
            outputs["likelihood"]["mc_uncertainty"]["sigma_objective"] = 0.125
            return outputs

    bounds = {name: {"lower": -2.0, "upper": 2.0} for name in ("x", "y")}
    param = {
        **_param(tmp_path), "parameter_space": icetune_lbfgs.icetune_main.build_parameter_space(bounds),
        "pickle_dump": pickle_dump,
    }
    driver = LikelihoodDriver()
    evaluate = icetune_lbfgs._evaluate_config_remote._function

    class Remote:
        # Accept the worker resource request without starting a Ray cluster
        def options(self, **kwargs):
            return self

        # Exercise the actual worker evaluation through the head result collector
        def remote(self, *args):
            return evaluate(*args)

    monkeypatch.setattr(icetune_lbfgs, "_evaluate_config_remote", Remote())
    monkeypatch.setattr(icetune_lbfgs.ray, "get", lambda values: values)
    objective = icetune_lbfgs.DistributedLBFGSObjective(
        args=_args(), bounds=bounds, names=["x", "y"], param=param,
        simdriver=driver, output_dir=str(tmp_path / "lbfgs"),
    )
    record = objective._evaluate_many([
        {"config": {"x": 1.4, "y": -0.2}, "trial_id": "likelihood-trial", "seed": 42, "label": "objective"}
    ])[0]
    assert record["success"]
    assert record["theta_hash"]
    assert "pickle_payload" not in record
    stored_path = (
        icetune_lbfgs.icetune_main.trial_pickle_dir(param)
        / icetune_lbfgs.icetune_main.trial_pickle_filename(record)
    )
    assert stored_path.is_file() is pickle_dump
    if pickle_dump:
        with stored_path.open("rb") as stream:
            stored = pickle.load(stream)
        assert stored["results"] == {"ok": True}
        assert stored["likelihood"] == record["likelihood"]
        assert stored["theta_hash"] == record["theta_hash"]
        assert stored["param"]["cdir"] == str(tmp_path)
    objective.publish_history()
    history = icetune_likelihood.load_likelihood_history(objective.history_path)
    np.testing.assert_allclose(history["chi2"], [1.0])
    np.testing.assert_allclose(history["X"], [[1.4, -0.2]])
    assert history["valid"].tolist() == [True]
    assert history["theta_hash"] == [record["theta_hash"]]
    assert history["mc_uncertainty"][0]["sigma_two_nll"] == pytest.approx(0.125)


# Check final output validation reports missing figures without rerunning a simulator
def test_ray_output_validation_rejects_missing_best(tmp_path):
    param = _ray_param(tmp_path)
    callback = icetune_ray.TrialOutputCallback(
        experiment_dir=str(tmp_path / "ray" / param["run_name"]),
        param=param,
        global_state=None,
    )
    best = SimpleNamespace(metrics={"loss": 0.0, "trial_id": "missing"})

    with pytest.raises(RuntimeError, match="final best figure"):
        callback.validate(best_result=best, require_initial=False)


# Inject the raw OOM exception through Ray's real actor manager and check native cleanup
def test_ray_oom_actor_cleanup(monkeypatch):
    from ray.air.execution._internal.actor_manager import RayActorManager
    from ray.air.execution._internal.tracked_actor import TrackedActor
    from ray.air.execution._internal.tracked_actor_task import TrackedActorTask
    from ray.exceptions import OutOfMemoryError, RayActorError

    monkeypatch.setattr(RayActorManager, "_actor_task_failed", RayActorManager._actor_task_failed)
    icetune_ray.install_ray_oom_handler()
    handler = RayActorManager._actor_task_failed
    icetune_ray.install_ray_oom_handler()
    assert RayActorManager._actor_task_failed is handler
    freed, errors = [], []
    manager = RayActorManager(SimpleNamespace(free_resources=lambda **kw: freed.append(kw)))
    actor = TrackedActor(7, on_error=lambda actor, error: errors.append(error))
    task = TrackedActorTask(actor, on_error=lambda actor, error: errors.append(error))
    manager._live_actors_to_ray_actors_resources[actor] = (None, "trial-cpus")
    oom = OutOfMemoryError("node running low on memory")
    manager._actor_task_failed(task, oom)
    assert manager.num_live_actors == 0
    assert freed == [{"acquired_resource": "trial-cpus"}]
    assert len(errors) == 2
    assert all(isinstance(error, RayActorError) and error.__cause__ is oom for error in errors)
    with pytest.raises(RuntimeError, match="unexpected exception"):
        manager._actor_task_failed(task, ValueError("invalid card"))


# Check campaign recovery is bounded and does not retry invalid input
@pytest.mark.parametrize("error,recoveries", [(ray.exceptions.OutOfMemoryError("OOM"), 2), (ValueError("card"), 0)])
def test_ray_campaign_recovery_limit(monkeypatch, error, recoveries):
    attempts, restores = [], []
    monkeypatch.setattr(icetune_ray.time, "sleep", lambda seconds: None)

    # Fail the real recovery boundary without executing a physics trial
    def fit():
        attempts.append(1)
        raise error

    tuner = SimpleNamespace(fit=fit)

    # Preserve one saved campaign across every restore attempt
    def restore():
        restores.append(1)
        return tuner

    with pytest.raises(type(error)):
        icetune_ray.fit_ray_tuner(tuner=tuner, restore=restore)
    assert len(restores) == recoveries
    assert len(attempts) == recoveries + 1


# Check a recovered campaign returns its results once without issuing another restore
def test_ray_campaign_recovery_success(monkeypatch):
    delays = []
    monkeypatch.setattr(icetune_ray.time, "sleep", delays.append)
    results = object()

    # Simulate a raw infrastructure failure escaping Tune
    def fail():
        raise ray.exceptions.OutOfMemoryError("OOM")

    restored = SimpleNamespace(fit=lambda: results)
    assert icetune_ray.fit_ray_tuner(
        tuner=SimpleNamespace(fit=fail), restore=lambda: restored, limit=1, delay_s=7
    ) is results
    assert delays == [7]


# Check worker loss, replacement identity and exhausted capacity using the real health state
def test_ray_worker_health_loss_and_rejoin(monkeypatch):
    marker = "icetune_worker_12_1"
    workers = icetune_ray.RayWorkers(False, poll_s=15, failure_limit=3, recovery_s=45)
    workers.node_ids[marker] = "old-node"
    clock = [100.0]
    monkeypatch.setattr(icetune_ray.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(ray, "is_initialized", lambda: True)
    monkeypatch.setattr(ray, "nodes", lambda: [])
    assert workers.poll()[marker]["count"] == 1
    clock[0] += 14
    assert workers.poll()[marker]["count"] == 1
    clock[0] += 1
    assert workers.poll()[marker]["count"] == 2
    workers.require_capacity()
    assert workers.empty_since is None
    clock[0] += 15
    assert workers.poll()[marker]["count"] == 3
    workers.require_capacity()
    clock[0] += 44
    workers.require_capacity()
    clock[0] += 2
    with pytest.raises(icetune_ray.RayClusterError):
        workers.require_capacity()
    workers.failed.clear()
    workers.require_capacity()
    assert workers.empty_since is None


# Check a successful real probe clears a timeout and stops repeated admission checks
def test_ray_worker_admission_recovery():
    ray.shutdown()
    marker = "icetune_worker_12_1"
    ray.init(num_cpus=2, include_dashboard=False, resources={marker: 1})
    workers = icetune_ray.RayWorkers(False, timeout_s=300)
    try:
        assert workers.poll() == {}
        actor, future, _ = workers.pending[marker]
        assert ray.get(future, timeout=60) == workers.node_ids[marker]

        # Keep a pending probe alive past the former fixed 120 second limit
        workers.pending[marker] = (actor, ray.ObjectRef.from_random(), time.monotonic() - 121)
        workers.checked_at = time.monotonic()
        assert workers.poll() == {}
        assert marker in workers.pending

        # Exercise the configured timeout with the same unresolved Ray object
        actor, future, _ = workers.pending[marker]
        workers.pending[marker] = (actor, future, time.monotonic() - workers.timeout_s - 1)
        failure = workers.poll()[marker]
        assert failure["count"] == 1
        assert "Ray worker admission timed out" in failure["error"]
        assert not workers.pending

        # Admit the same node again and consume the completed probe
        workers.checked_at = 0.0
        workers.poll()
        assert ray.get(workers.pending[marker][1], timeout=60) == workers.node_ids[marker]
        assert workers.poll() == {}
        assert not workers.pending

        # A later scan must leave the recovered node admitted without another probe
        workers.checked_at = 0.0
        assert workers.poll() == {}
        assert not workers.pending
        workers.require_capacity()
        assert workers.empty_since is None
    finally:
        workers.close()
        ray.shutdown()


# Inject a raw memory kill into Tune and verify retry and restore retain completed trials
def test_ray_native_oom_retry_and_restore(tmp_path, monkeypatch):
    _ensure_ray(local_mode=False)

    class MemoryTrial(tune.Trainable):
        # Expose the actor handle only so the test can emulate Ray's memory monitor
        def step(self):
            return {"loss": self.config["x"], "kill_actor": ray.get_runtime_context().current_actor}

    original_get = ray.get
    killed, evaluated = [], []

    # Kill one real trial actor and deliver the same raw exception as the memory monitor
    def get(*args, **kwargs):
        value = original_get(*args, **kwargs)
        if isinstance(value, dict) and "kill_actor" in value:
            actor = value.pop("kill_actor")
            evaluated.append(value["loss"])
            if not killed:
                killed.append(True)
                ray.kill(actor, no_restart=True)
                raise ray.exceptions.OutOfMemoryError("injected worker memory kill")
        return value

    monkeypatch.setattr(ray, "get", get)
    icetune_ray.install_ray_oom_handler()
    try:
        tuner = tune.Tuner(MemoryTrial, param_space={"x": tune.grid_search([1.0, 2.0])},
                           tune_config=tune.TuneConfig(max_concurrent_trials=1),
                           run_config=tune.RunConfig(storage_path=str(tmp_path), name="oom", verbose=0,
                                                     stop={"training_iteration": 1},
                                                     checkpoint_config=tune.CheckpointConfig(checkpoint_at_end=False, checkpoint_frequency=0),
                                                     failure_config=tune.FailureConfig(max_failures=3)))
        results = tuner.fit()
        assert sorted(result.metrics["loss"] for result in results) == [1.0, 2.0]
        assert all(result.error is None for result in results)
        assert len(evaluated) == 3
        restored = tune.Tuner.restore(str(tmp_path / "oom"), trainable=MemoryTrial).fit()
        assert sorted(result.metrics["loss"] for result in restored) == [1.0, 2.0]
        assert len(evaluated) == 3
    finally:
        ray.shutdown()


# Exercise real Tune acceptance while a separate renderer waits for later trials
@pytest.mark.parametrize("pickle_dump", [False, True])
def test_ray_rendering_releases_trial_slots(tmp_path, pickle_dump):
    _ensure_ray(local_mode=False)
    count = 3

    class WaitingDriver(_RayFigureDriver):
        # Apply the same worker setup API as the simulator drivers
        def runtime_environment(self, **kwargs):
            return {}

        # Record each completed computation before any plotting starts
        def evaluate_trial_outputs(self, **kwargs):
            outputs = super().evaluate_trial_outputs(**kwargs)
            Path(kwargs["param"]["cdir"], "tmp", outputs["trial_id"]).touch()
            return outputs

        # Require later computations to finish while the first render is active
        def render_trial_figures_to_dir(self, **kwargs):
            assert os.getpid() != head_pid
            deadline = time.monotonic() + 30.0
            while len(list((tmp_path / "tmp").glob("????????"))) < count:
                assert time.monotonic() < deadline, "Rendering blocked later trial computations"
                time.sleep(0.05)
            return super().render_trial_figures_to_dir(**kwargs)

    head_pid = os.getpid()
    param = {**_ray_param(tmp_path), "pickle_dump": pickle_dump,
             "optimization": {"backend": "ray", "cost": "loss", "optimizer": "basic"},
             "parameter_space": [{"name": "x", "lower": 0.0, "upper": 1.0}],
             "parameter_topology": {"schema_version": 1, "groups": []}}
    (tmp_path / "tmp").mkdir()
    experiment = tmp_path / "ray" / "run"
    actor = icetune_ray.GlobalState.remote(experiment_dir=str(experiment), cdir=str(tmp_path),
                                          run_name=param["run_name"], cost="loss", render_figures=True,
                                          keep_pickles=pickle_dump)
    driver = WaitingDriver()
    outputs = icetune_ray.TrialOutputCallback(experiment_dir=str(experiment), param=param,
                                            global_state=actor, simdriver=driver)
    search = _search(_search_args("basic", num_trials=count, max_concurrent_trials=1), {"x": 0.9})
    history = icetune_ray.HistoryPlotCallback(experiment_dir=str(experiment), cdir=str(tmp_path),
                                            run_name=param["run_name"], cost="loss", interval_s=300,
                                            search_alg=search, param=param, num_trials=count, max_concurrent_trials=1)
    trainable = tune.with_resources(tune.with_parameters(icetune_ray.CFunc, simdriver=driver, param=param,
                                                        global_state=actor, initial_points={"x": 0.9}), {"cpu": 1})
    try:
        results = tune.Tuner(trainable, tune_config=tune.TuneConfig(search_alg=search, num_samples=count),
                             run_config=tune.RunConfig(storage_path=str(experiment.parent), name=experiment.name,
                                                       callbacks=[outputs, history], verbose=0,
                                                       stop={"training_iteration": 1},
                                                       checkpoint_config=tune.CheckpointConfig(checkpoint_at_end=False),
                                                       failure_config=tune.FailureConfig(max_failures=3))).fit()
        best = results.get_best_result(metric="loss", mode="min")
        outputs.validate(best_result=best, require_initial=True)
        history.close()
        records = json.loads(history.history_path.read_text())["trials"]
        assert len(records) == count
        assert all(row["elapsed_seconds"] > 0.0 for row in records)
        assert len(list((experiment / "results").glob("*.pkl"))) == (count if pickle_dump else 0)
        restored = icetune_ray.BestState(experiment_dir=str(experiment), cdir=str(tmp_path),
                                        run_name=param["run_name"], cost="loss", render_figures=True,
                                        keep_pickles=pickle_dump)
        for result in results:
            if result.metrics.get("ray_output_token"):
                restored.recover_trial_outputs(result.metrics)
        assert restored.snapshot()["published_cost"] == pytest.approx(best.metrics["loss"])
    finally:
        history.close()
        ray.shutdown()


# Check actual worker timeouts become failed history entries with measured runtime
@pytest.mark.parametrize("algorithm", ["basic", "lbfgs"])
def test_worker_timeout_records_runtime(tmp_path, algorithm):
    class SlowDriver(_ToyDriver):
        # Exceed the shared worker budget inside the actual driver API
        def evaluate_trial_outputs(self, **kwargs):
            time.sleep(2.0)
            Path(kwargs["param"]["cdir"], "completed").touch()
            return super().evaluate_trial_outputs(**kwargs)

    param = {**_param(tmp_path), "max_t": 0.1}
    config = {"x": 0.5, "y": 0.0}
    if algorithm == "lbfgs":
        result = icetune_lbfgs._evaluate_config_remote._function(config, param, SlowDriver(), "trial-timeout", 123)
        assert not result["success"]
        assert "max_t" in result["error"]
    else:
        _ensure_ray(local_mode=False)
        search = _search(_search_args("basic", num_trials=1))
        trainable = tune.with_parameters(icetune_ray.CFunc, simdriver=SlowDriver(), param=param, global_state=None)
        try:
            results = tune.Tuner(trainable, tune_config=tune.TuneConfig(search_alg=search, num_samples=1),
                                 run_config=tune.RunConfig(storage_path=str(tmp_path), name="timeout", verbose=0,
                                                           stop={"training_iteration": 1},
                                                           failure_config=tune.FailureConfig(max_failures=3))).fit()
            result = results[0].metrics
            assert result["trial_failed"] and results[0].error is None
            assert "max_t" in result["trial_error"]
            assert len(search.failed_records) == 1
            assert search.failed_records[0]["elapsed_seconds"] == pytest.approx(result["elapsed_seconds"])
        finally:
            ray.shutdown()
    assert result["elapsed_seconds"] >= param["max_t"]
    assert result["elapsed_seconds"] == pytest.approx(
        result["completed_at_unix"] - result["started_at_unix"]
    )
    assert not (tmp_path / "completed").exists()


# Preserve shared Pandora inputs even beneath excluded checkout directories
def test_pandora_uploaded_inputs(tmp_path, monkeypatch):
    from core.io.files import file_identity
    from core.tune.drivers.pandora.driver import PandoraDriver
    from core.tune.drivers.pandora.runtime import stage_file

    source, worker = tmp_path / "source", tmp_path / "worker"
    input_path = source / "output/input.root"
    input_path.parent.mkdir(parents=True)
    input_path.write_bytes(b"simulation input")
    identity = file_identity(input_path)
    worker.mkdir()
    monkeypatch.chdir(worker)
    param = ray_runtime.worker_runtime({"cdir": str(source), "ray_upload_runtime": True,
        "datacards": [{"inputs": [{"path": str(input_path), **identity}]}]}, PandoraDriver())
    assert param["cdir"] == str(worker)
    path = param["datacards"][0]["inputs"][0]["path"]
    assert path == str(input_path)
    assert stage_file(path, worker / "cache/input.root", identity).read_bytes() == input_path.read_bytes()


# Reject invalid L-BFGS controls before any worker allocation
@pytest.mark.parametrize("key,value", [("maxiter", 0), ("grad_step_rel", float("inf")),
                                      ("grad_step_abs", 0), ("line_search_decay", 1.0),
                                      ("retry_failed_search", 1)])
def test_lbfgs_invalid_settings(tmp_path, key, value):
    settings = load_lbfgs_settings(resource("tune/settings/lbfgs.json"))
    settings[key] = value
    path = tmp_path / "settings.json"
    path.write_text(json.dumps(settings))
    with pytest.raises(ValueError):
        load_lbfgs_settings(path)


# A repeatedly failing initial render must not block a successful best candidate
@pytest.mark.parametrize("published_cost", [None, 0.5])
def test_failed_initial_render(tmp_path, published_cost):
    callback = icetune_ray.TrialOutputCallback(experiment_dir=str(tmp_path / "run"),
        param=dict(cdir=str(tmp_path), run_name="plots", cost="loss"), global_state=None)
    callback.render_records = {
        "initial": dict(token="initial", trial_id="initial", cost=2.0, initial=True),
        "best": dict(token="best", trial_id="best", cost=1.0, initial=False)}
    callback.output_jobs["render"]["failures"]["initial"] = 3
    callback.rendered_cost = published_cost
    selected = callback._next_render()
    if published_cost is None:
        assert selected["trial_id"] == "best" and selected["render_best"] and not selected["render_initial"]
    else:
        assert selected is None

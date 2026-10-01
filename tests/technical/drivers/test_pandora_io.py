# Pandora runtime I/O recovery and Ray trial failures
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import ctypes.util
import errno
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
from core.tune import io
from core.tune.drivers.pandora import runtime
from core.tune.drivers.pandora.driver import PandoraDriver
from core.tune.tunesetup import load_tunesetup

from submit import campaign_source
from submit.lxplus import cluster, nodes


# Load the production Pandora catalog before exercising runtime I/O failures
@pytest.fixture
def pandora_aux(pandora_inputs):
    return load_tunesetup(cdir=Path(__file__).resolve().parents[3], simdriver="PANDORA",
                         name=campaign_source("tune-pandora-v0")).aux_param_space


# Freeze real files and expose a kernel EIO through the external setup dependency
@pytest.fixture
def frozen(tmp_path):
    if not Path("/proc/self/mem").is_file():
        pytest.skip("Kernel EIO test requires procfs")
    source = tmp_path / "source"
    source.mkdir()
    setup = source / "setup.sh"
    setup.write_text("export PANDORA_TEST=1\n")
    (source / "settings.xml").write_text("<pandora/>\n")
    manifest, archive = runtime.freeze_runtime(
        source, ["settings.xml"], tmp_path / "frozen",
        dict(setup="", dependencies={str(setup): runtime.file_identity(setup)}),
    )
    manifest["archive"] = str(archive)
    saved = setup.with_suffix(".sh._old")
    setup.rename(saved)
    setup.symlink_to("/proc/self/mem")
    return manifest, setup, saved


# Retry actual failed reads and keep the frozen checksum check after recovery
@pytest.mark.parametrize("recovery", ["valid", "changed", "absent"])
def test_dependency_read(frozen, tmp_path, monkeypatch, recovery):
    manifest, setup, saved = frozen
    delays = []

    # Restore the dependency during the retry delay without replacing the file reader
    def recover(delay):
        delays.append(delay)
        if recovery != "absent":
            if recovery == "changed":
                saved.write_text("changed setup\n")
            setup.rename(setup.with_suffix(".sh.failed"))
            saved.replace(setup)

    monkeypatch.setattr(io.time, "sleep", recover)
    if recovery == "valid":
        root = runtime.stage_runtime(manifest, tmp_path, tmp_path / "cache")
        assert runtime.file_identity(root / "setup.sh") == manifest["files"]["setup.sh"]
    elif recovery == "changed":
        with pytest.raises(ValueError, match="dependency changed"):
            runtime.stage_runtime(manifest, tmp_path, tmp_path / "cache")
    else:
        with pytest.raises(OSError) as caught:
            runtime.stage_runtime(manifest, tmp_path, tmp_path / "cache")
        assert caught.value.errno == errno.EIO
    assert len(delays) == (io.IO_RETRIES - 1 if recovery == "absent" else 1)


# Preserve exhausted dependency I/O failures through the real Pandora and Ray APIs
def test_trial_io(frozen, pandora_aux, tmp_path, monkeypatch):
    pytest.importorskip("ray")
    from core.tune.backends.ray import CFunc
    from ray import tune

    manifest, _, _ = frozen
    monkeypatch.setattr("ray.tune.trainable.trainable.DEFAULT_STORAGE_PATH", str(tmp_path))
    monkeypatch.setattr(io.time, "sleep", lambda delay: None)
    param = dict(cdir=str(tmp_path), run_name="io", cost="pflow", max_t=60,
                 aux_param_space=pandora_aux, datacards=[dict(weight=1.0, runtime=manifest)])
    trial = tune.with_parameters(CFunc, simdriver=PandoraDriver(), param=param, global_state=None)(
        config={},
    )
    try:
        with pytest.raises(OSError) as caught:
            trial.step()
        assert caught.value.errno == errno.EIO
    finally:
        trial.stop()


# Retry storage outages while rejecting invalid paths, permissions and full disks
@pytest.mark.parametrize("code,retry", [(errno.EIO, True), (errno.ESTALE, True), (errno.ETIMEDOUT, True),
                                       (errno.ENOENT, False), (errno.EACCES, False), (errno.EPERM, False),
                                       (errno.ENOSPC, False)])
def test_ray_io_errors(code, retry):
    pytest.importorskip("ray")
    from core.tune.backends.ray import ray_infrastructure_error

    error = OSError(code, "runtime dependency")
    wrapped = RuntimeError("trial failed")
    wrapped.__cause__ = error
    assert ray_infrastructure_error(error) is retry
    assert ray_infrastructure_error(wrapped) is retry


# Reject unreadable or changed dependencies before the library command or Ray can start
@pytest.mark.parametrize("failure", ["eio", "changed", "missing"])
def test_worker_dependency(frozen, tmp_path, monkeypatch, failure):
    manifest, setup, saved = frozen
    monkeypatch.setattr(io.time, "sleep", lambda delay: None)
    if failure != "eio":
        setup.rename(setup.with_suffix(".sh.failed"))
        if failure == "changed":
            setup.write_text("changed setup\n")
    executed = tmp_path / "library-check"
    check = dict(command=[sys.executable, "-c", "from pathlib import Path; Path(__import__('sys').argv[1]).touch()",
                          str(executed)], env=runtime.command_environment(inherit=False), timeout=10)
    result = nodes.start_ray_worker_node(cpus=1, temp_dir=str(tmp_path / "irt-0123456789ab"), gpus=0,
        head_address="127.0.0.1:10000", resource_marker="icetune_worker_12_1",
        dependencies=manifest["dependencies"], worker_check=check)
    assert result["status"] == "failed" and result["retry"] is False
    assert str(setup) in result["error"]
    assert ("Input/output error" if failure == "eio" else
            "dependency changed" if failure == "changed" else "No such file") in result["error"]
    assert not executed.exists() and not (tmp_path / "irt-0123456789ab").exists()


# A real dynamic loader error must reject the allocation without starting Ray or retrying admission
@pytest.mark.parametrize("failure", ["library", "timeout"])
def test_worker_library_failure(tmp_path, failure):
    setup = tmp_path / "setup.sh"
    setup.write_text("# verified dependency\n")
    command = [sys.executable, "-c", "import ctypes, sys; ctypes.CDLL(sys.argv[1])", str(tmp_path / "missing.so")]
    if failure == "timeout":
        command = [sys.executable, "-c", "import time; time.sleep(60)"]
    check = dict(command=command, env=runtime.command_environment(inherit=False), timeout=0.1 if failure == "timeout" else 10)
    result = nodes.start_ray_worker_node(cpus=1, temp_dir=str(tmp_path / "irt-0123456789ab"), gpus=0,
        head_address="127.0.0.1:10000", resource_marker="icetune_worker_12_1",
        dependencies={str(setup): runtime.file_identity(setup)}, worker_check=check)
    assert result["status"] == "failed" and result["retry"] is False
    assert ("missing.so" if failure == "library" else "timed out") in result["error"]
    assert not (tmp_path / "irt-0123456789ab").exists()
    state = dict(pending=[], seen=set(), joined=set(), failed=set(), attempts={}, errors={})
    cluster.apply_ray_worker_starts(results={"worker": result}, head_address="127.0.0.1:10000", retry_limit=3, **state)
    assert state["failed"] == {"worker"} and not state["joined"] and not state["pending"]


# Read real dependencies and load a system shared library in an isolated child environment
def test_worker_runtime(tmp_path, monkeypatch):
    library = ctypes.util.find_library("m")
    assert library
    setup = tmp_path / "setup.sh"
    setup.write_text("# verified dependency\n")
    marker = tmp_path / "loaded"
    monkeypatch.setenv("PYTHONPATH", str(tmp_path / "wrong-python"))
    monkeypatch.setenv("LD_LIBRARY_PATH", str(tmp_path / "wrong-libraries"))
    command = [sys.executable, "-c", """import ctypes, os, pathlib, sys
assert 'PYTHONPATH' not in os.environ and 'LD_LIBRARY_PATH' not in os.environ
ctypes.CDLL(sys.argv[1], mode=ctypes.RTLD_GLOBAL)
pathlib.Path(sys.argv[2]).write_text('loaded')
""", library, str(marker)]
    nodes.check_worker_runtime({str(setup): runtime.file_identity(setup)},
        dict(command=command, env=runtime.command_environment(inherit=False), timeout=10))
    assert marker.read_text() == "loaded"


# Preserve the Ray node identity of an actual kernel I/O error for worker replacement
def test_trial_io_node(frozen, pandora_aux, tmp_path, monkeypatch):
    ray = pytest.importorskip("ray")
    from core.tune.backends.ray import CFunc, HistoryPlotCallback, RayProbe
    from ray import tune

    manifest, setup, saved = frozen
    monkeypatch.setattr("ray.tune.trainable.trainable.DEFAULT_STORAGE_PATH", str(tmp_path))
    monkeypatch.setattr(io.time, "sleep", lambda delay: None)
    ray.shutdown()
    ray.init(num_cpus=1, include_dashboard=False)
    trial = None
    callback = HistoryPlotCallback(experiment_dir=str(tmp_path / "results"), cdir=str(tmp_path), cost="pflow",
        run_name="io", interval_s=60, search_alg=None, num_trials=1, max_concurrent_trials=1,
        param={"aux_param_space": {"runtime_dependencies": manifest["dependencies"]}})
    try:
        node_id = ray.get_runtime_context().get_node_id()
        callback.workers.node_ids["icetune_worker_12_1"] = node_id
        param = dict(cdir=str(tmp_path), run_name="io", cost="pflow", max_t=60,
                     aux_param_space=pandora_aux, datacards=[dict(weight=1.0, runtime=manifest)])
        trial = tune.with_parameters(CFunc, simdriver=PandoraDriver(), param=param, global_state=None)(config={})
        with pytest.raises(OSError) as caught:
            trial.step()
        assert caught.value.errno == errno.EIO
        from traceback import format_exception
        error = "".join(format_exception(caught.value))
        callback._worker_failure(SimpleNamespace(last_result={"trial_error": error}))
        assert callback.workers.failed["icetune_worker_12_1"]["count"] == 1
        assert callback.workers.dependencies == manifest["dependencies"]
        with pytest.raises(OSError):
            RayProbe().ready(False, callback.workers.dependencies)
        setup.rename(setup.with_suffix(".sh.failed"))
        saved.replace(setup)
        assert RayProbe().ready(False, callback.workers.dependencies) == node_id
    finally:
        if trial is not None:
            trial.stop()
        callback.workers.close()
        ray.shutdown()

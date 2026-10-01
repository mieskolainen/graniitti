# Test command execution, logging and runtime environment paths
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import signal
import sys
import time

import pytest
from core.tune.drivers.graniitti import runtime as graniitti_runtime
from core.tune.runtime import process as iceruntime


# Prefer GRANIITTI_IO_PATH when constructing external library paths
def test_io_path_environment(monkeypatch):
    monkeypatch.setenv("GRANIITTI_IO_PATH", "/tmp/graniitti-io")

    assert graniitti_runtime.default_library_path() == "/tmp/graniitti-io"


# Fall back to HOME/local when GRANIITTI_IO_PATH is not exported
def test_io_path_falls_back_to_home(monkeypatch):
    monkeypatch.delenv("GRANIITTI_IO_PATH", raising=False)
    monkeypatch.setenv("HOME", "/tmp/graniitti-home")

    assert graniitti_runtime.default_library_path() == "/tmp/graniitti-home/local"


# Keep icepack imports available alongside libraries selected by GRANIITTI_IO_PATH
def test_setpaths_io_path(monkeypatch):
    monkeypatch.setenv("GRANIITTI_IO_PATH", "/tmp/graniitti-io")
    monkeypatch.setenv("LD_LIBRARY_PATH", "")
    monkeypatch.setenv("PYTHONPATH", "")

    variables = graniitti_runtime.environment(
        cdir="/tmp/graniitti-repo",
        libdir=None,
        python_version="3.12",
    )

    ld_library_path = variables["LD_LIBRARY_PATH"]
    pythonpath = variables["PYTHONPATH"]
    assert ld_library_path.startswith("/tmp/graniitti-io/HEPMC3/lib:")
    assert "/tmp/graniitti-io/LHAPDF/lib64" in ld_library_path
    assert pythonpath.startswith("/tmp/graniitti-repo:")
    assert "/tmp/graniitti-io/HEPMC3/lib/python3.12/site-packages" in pythonpath


# Preserve child failure details in the structured command result
def test_execute_cmd_preserves_child_failure():
    result = {}
    status = iceruntime.execute_cmd(
        [
            sys.executable,
            "-c",
            "import sys; print('graniitti-child-output'); sys.exit(4)",
        ],
        max_t=5,
        result_out=result,
    )

    assert status is False
    assert result["status"] == "failed"
    assert result["returncode"] == 4
    assert "graniitti-child-output" in result["output"]


# Extract the GRANIITTI exception instead of a destructor line
def test_command_error_prefers_simulator_exception():
    output = """
banner
Exception catched: MHelicityConfig::ProcessHelicityStructure: Kinematic coupling problem
~MGraniitti [DONE]
"""

    assert (
        iceruntime.command_output_root_cause(output)
        == "MHelicityConfig::ProcessHelicityStructure: Kinematic coupling problem"
    )
    assert "Exception catched:" in iceruntime.command_output_tail(output)


# Preserve the algorithm failure when Gaudi reports subsequent initialization failures
@pytest.mark.parametrize("cause", [
    "Tracking ERROR Cannot open calibration input",
    "FileNotFoundError: calibration input is absent",
    "PandoraApi::ReadSettings return STATUS_CODE_NOT_FOUND",
    None,
])
def test_command_error_before_gaudi_shutdown(cause):
    shutdown = (
        "EventLoopMgr ERROR Unable to initialize Algorithm: k4FWCore__Sequencer\n"
        "ServiceManager ERROR Unable to initialize Service: EventLoopMgr\n"
        "ApplicationMgr ERROR Application Manager Terminated with error code 1"
    )
    output = f"{cause or ''}\nDDMarlinPandora INFO MaxErrorCount: 10\nDDMarlinPandora INFO Init processor\n{shutdown}"
    assert iceruntime.command_output_root_cause(output) == (cause or shutdown.splitlines()[-1])


# Check that external simulators can be pinned to a scratch-local project tree
def test_execute_cmd_requested_working_dir(tmp_path):
    marker = tmp_path / "cwd.txt"

    status = iceruntime.execute_cmd(
        [
            sys.executable,
            "-c",
            "import pathlib; pathlib.Path('cwd.txt').write_text(str(pathlib.Path.cwd()))",
        ],
        max_t=5,
        cwd=tmp_path,
    )

    assert status is True
    assert marker.read_text(encoding="utf-8") == str(tmp_path)


def test_execute_cmd_writes_log_file_success(tmp_path):
    log_path = tmp_path / "success.json"

    status = iceruntime.execute_cmd(
        [
            sys.executable,
            "-c",
            "print('graniitti-success-output')",
        ],
        max_t=5,
        log_path=str(log_path),
        log_metadata={"stage": "unit_success"},
    )

    assert status is True
    payload = json.loads(log_path.read_text(encoding="utf-8"))
    assert payload["status"] == "ok"
    assert payload["returncode"] == 0
    assert payload["metadata"]["stage"] == "unit_success"
    assert "graniitti-success-output" in payload["output"]
    assert not (tmp_path / "success.txt").exists()


def test_execute_cmd_writes_log_file_failure(tmp_path):
    log_path = tmp_path / "failure.json"

    status = iceruntime.execute_cmd(
        [
            sys.executable,
            "-c",
            "import sys; print('graniitti-failure-output'); sys.exit(7)",
        ],
        max_t=5,
        log_path=str(log_path),
        log_metadata={"stage": "unit_failure"},
    )

    assert status is False
    payload = json.loads(log_path.read_text(encoding="utf-8"))
    assert payload["status"] == "failed"
    assert payload["returncode"] == 7
    assert payload["metadata"]["stage"] == "unit_failure"
    assert "graniitti-failure-output" in payload["output"]
    assert not (tmp_path / "failure.txt").exists()


def test_execute_cmd_writes_log_file_timeout(tmp_path):
    log_path = tmp_path / "timeout.json"

    status = iceruntime.execute_cmd(
        [
            sys.executable,
            "-c",
            "import time; print('graniitti-timeout-output', flush=True); time.sleep(2)",
        ],
        max_t=0.1,
        log_path=str(log_path),
        log_metadata={"stage": "unit_timeout"},
    )

    assert status is False
    payload = json.loads(log_path.read_text(encoding="utf-8"))
    assert payload["status"] == "timeout"
    assert payload["returncode"] is None
    assert payload["metadata"]["stage"] == "unit_timeout"
    assert "graniitti-timeout-output" in payload["output"]
    assert not (tmp_path / "timeout.txt").exists()


# Check the child can observe its running log before completion
def test_running_log_before_completion(tmp_path):
    log_path = tmp_path / "running.json"

    status = iceruntime.execute_cmd(
        [
            sys.executable,
            "-c",
            "import json, sys; payload = json.load(open(sys.argv[1])); "
            "assert payload['status'] == 'running'; "
            "assert payload['returncode'] is None; assert payload['output'] == ''; print('done')",
            str(log_path),
        ],
        max_t=5,
        log_path=str(log_path),
        log_metadata={"stage": "unit_running"},
    )

    assert status is True
    payload = json.loads(log_path.read_text(encoding="utf-8"))
    assert payload["status"] == "ok"
    assert payload["metadata"]["stage"] == "unit_running"
    assert "done" in payload["output"]
    assert not (tmp_path / "running.txt").exists()


# Check timeout stops descendants that would otherwise keep simulating
@pytest.mark.parametrize("whole_trial", [False, True])
def test_execute_cmd_timeout_stops_descendants(tmp_path, whole_trial):
    marker = tmp_path / "surviving-child.txt"
    child = "import pathlib, sys, time; time.sleep(0.8); pathlib.Path(sys.argv[1]).write_text('alive')"
    parent = (
        "import subprocess, sys, time; "
        f"subprocess.Popen([sys.executable, '-c', {child!r}, sys.argv[1]]); "
        "print('child-started', flush=True); time.sleep(20)"
    )
    cmd = [sys.executable, "-c", parent, str(marker)]
    if whole_trial:
        handler = signal.getsignal(signal.SIGALRM)
        with pytest.raises(iceruntime.TrialTimeout), iceruntime.trial_limit(time.monotonic() + 0.3):
            iceruntime.run_process(cmd)
        assert signal.getsignal(signal.SIGALRM) == handler
        assert signal.getitimer(signal.ITIMER_REAL)[0] <= 0.0
    else:
        result = {}
        assert not iceruntime.execute_cmd(cmd, max_t=0.3, result_out=result)
        assert result["status"] == "timeout"
        assert "child-started" in result["output"]
    time.sleep(1.0)
    assert not marker.exists()

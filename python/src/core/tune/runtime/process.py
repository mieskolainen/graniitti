# Process and runtime utilities shared by the Python command-line tools
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import contextlib
import copy
import os
import random
import re
import shlex
import signal
import socket
import subprocess
import sys
import threading
import time
from dataclasses import asdict, dataclass
from datetime import datetime

import numpy as np

from core.io.files import ensure_dir
from core.io.serialize import write_json_file


# Identify an exhausted trial budget separately from infrastructure failures
class TrialTimeout(TimeoutError):
    """The complete worker computation exceeded its wall time budget."""


# Bound one worker attempt and let subprocess cleanup handle interrupted commands
@contextlib.contextmanager
def trial_limit(deadline: float):
    remaining = float(deadline) - time.monotonic()
    if remaining <= 0.0:
        raise TrialTimeout("Trial exceeded max_t")
    if threading.current_thread() is not threading.main_thread():
        raise RuntimeError("Trial time limits require the worker main thread")

    # Interrupt Python waits and propagate the deadline through command cleanup
    def expired(signum, frame):
        raise TrialTimeout("Trial exceeded max_t")

    handler = signal.getsignal(signal.SIGALRM)
    previous = signal.getitimer(signal.ITIMER_REAL)
    started = time.monotonic()
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, min(remaining, previous[0]) if previous[0] > 0.0 else remaining)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, handler)
        if previous[0] > 0.0:
            signal.setitimer(signal.ITIMER_REAL, max(1e-6, previous[0] - (time.monotonic() - started)), previous[1])


@dataclass(frozen=True)
class CommandResult:
    """Immutable result from one external command"""

    argv: list[str]
    command: str
    duration_seconds: float
    end_unix: float
    max_t: float
    metadata: dict | None
    output: str
    returncode: int | None
    start_unix: float
    status: str

    # Convert the immutable result into a caller-owned JSON payload
    def to_dict(self) -> dict:
        return asdict(self)


# Convert subprocess output into text without losing undecodable bytes
def _normalize_output(output) -> str:
    if output is None:
        return ""
    if isinstance(output, bytes):
        return output.decode(errors="replace")
    return str(output)


# Build one immutable command result
def _command_result(cmd, max_t, output, start, end, status, returncode, metadata) -> CommandResult:
    return CommandResult(
        argv=[str(value) for value in cmd],
        command=shlex.join(map(str, cmd)),
        duration_seconds=end - start,
        end_unix=end,
        max_t=max_t,
        metadata=copy.deepcopy(metadata),
        output=_normalize_output(output),
        returncode=returncode,
        start_unix=start,
        status=status,
    )


# Publish one command result to optional memory and disk destinations
def _publish_result(result: CommandResult, result_out: dict | None, log_path: str | None) -> None:
    payload = result.to_dict()
    if result_out is not None:
        result_out.clear()
        result_out.update(copy.deepcopy(payload))
    if log_path is not None:
        ensure_dir(os.path.dirname(log_path))
        write_json_file(log_path, payload, indent=2, sort_keys=True, newline=True)


# Compute a bounded tail suitable for persistent failure diagnostics
def command_output_tail(output, *, max_lines: int = 20, max_chars: int = 4000) -> str:
    lines = _normalize_output(output).strip().splitlines()
    tail = "\n".join(lines[-max(1, int(max_lines)) :])
    return tail[-max_chars:] if len(tail) > max_chars else tail


# Extract the most informative single-line cause from command output
def command_output_root_cause(output) -> str | None:
    lines = [line.strip() for line in _normalize_output(output).splitlines() if line.strip()]
    shutdown = re.compile(
        r"(?:Unable to initialize (?:Algorithm|Service):|Application Manager Terminated with error code)",
        re.IGNORECASE,
    )
    patterns = (
        re.compile(r"Exception\s+(?:catched|caught):\s*(.+)", re.IGNORECASE),
        re.compile(
            r'.*(?:fatal error|segmentation fault|illegal instruction|version [`\'"].+not found).*', re.IGNORECASE
        ),
        re.compile(r".*(?:\b\w*(?:Error|Exception):|\bSTATUS_CODE_(?!SUCCESS\b)[A-Z_]+).*"),
        re.compile(r".*\b(?:ERROR|FATAL)\b.*"),
        re.compile(r".*(?:\berror\b|\bfailed\b|not found|\bkilled\b|\baborted\b).*", re.IGNORECASE),
    )
    for candidates in ([line for line in lines if not shutdown.search(line)], lines):
        for pattern in patterns:
            for line in reversed(candidates):
                match = pattern.search(line)
                if match is not None:
                    return match.group(1).strip() if match.lastindex else line
    ignored = re.compile(r"^(?:~\S+\s+\[DONE\]|[-.`]+)$")
    return next((line for line in reversed(lines) if ignored.fullmatch(line) is None), None)


# Run a command with process group cleanup on timeout or interruption
def run_process(cmd, *, cwd=None, stdout=subprocess.PIPE, timeout=None, env=None) -> subprocess.CompletedProcess:
    with subprocess.Popen(
        cmd, cwd=cwd, stdout=stdout, stderr=subprocess.STDOUT, text=True, env=env, start_new_session=True
    ) as process:
        try:
            output, _ = process.communicate(timeout=timeout)
        except BaseException as exc:
            # Stop descendants as well as the command before draining its output
            with contextlib.suppress(ProcessLookupError):
                os.killpg(process.pid, signal.SIGKILL)
            output, _ = process.communicate()
            if isinstance(exc, subprocess.TimeoutExpired):
                exc.stdout = output
            raise
    return subprocess.CompletedProcess(cmd, process.returncode, stdout=output)


# Execute one external command and publish structured diagnostics
def execute_cmd(
    cmd: list,
    max_t: float = 3600,
    log_path: str | None = None,
    log_metadata: dict | None = None,
    cwd: str | os.PathLike[str] | None = None,
    result_out: dict | None = None,
) -> bool:
    command = shlex.join(map(str, cmd))
    print(f"{__name__}.execute_cmd: argv = {cmd} (max_t = {max_t})")
    print(f"{__name__}.execute_cmd: command = {command}")
    start = time.time()
    _publish_result(_command_result(cmd, max_t, "", start, start, "running", None, log_metadata), result_out, log_path)

    try:
        process = run_process(cmd, cwd=None if cwd is None else os.fspath(cwd), timeout=max_t)
        status = "ok" if process.returncode == 0 else "failed"
        output, returncode = process.stdout, process.returncode
        message = "Command succeeded" if returncode == 0 else f"Command failed with exit code {returncode}"
    except subprocess.TimeoutExpired as exc:
        status, output, returncode = "timeout", exc.stdout, None
        message = f"Command timed out after {exc.timeout} seconds"
    except FileNotFoundError:
        status, output, returncode = "not_found", "", None
        message = "Command not found; check the command name or path"
    except OSError as exc:
        status, output, returncode = "os_error", str(exc), None
        message = f"OS error occurred: {exc}"
    except Exception as exc:
        status, output, returncode = "exception", str(exc), None
        message = f"Unexpected error occurred: {exc}"

    end = time.time()
    result = _command_result(cmd, max_t, output, start, end, status, returncode, log_metadata)
    _publish_result(result, result_out, log_path)
    print(f"{__name__}.execute_cmd: {message} after {end - start:.1f} sec")
    if status != "ok" and result.output.strip():
        print(f"{__name__}.execute_cmd: Command output follows:\n{result.output.rstrip()}")
    return status == "ok"


# Compute the best available local IP address
def get_ip() -> str:
    with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as sock:
        try:
            sock.connect(("10.255.255.255", 1))
            return sock.getsockname()[0]
        except OSError:
            return "127.0.0.1"


# Seed Python, NumPy and optional PyTorch random number generators
def set_random_seeds(seed: int = 42, use_torch: bool = False) -> None:
    print(f"Setting random seed: {seed}")
    random.seed(seed)
    np.random.seed(seed)
    if use_torch:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


# Limit native numerical libraries to the CPUs assigned to one trial
def configure_numerical_threads(thread_count: int) -> int:
    count = max(1, int(thread_count))
    for name in (
        "BLIS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[name] = str(count)
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(limits=count, user_api="blas")
    except ImportError:
        pass
    # Update an already loaded Torch runtime as well as future library imports
    if "torch" in sys.modules:
        sys.modules["torch"].set_num_threads(count)
    return count


# Compute a filesystem-safe current timestamp
def get_current_time() -> str:
    return f"{datetime.now()}".replace(":", "-").replace(" ", "--").split(".")[0]


# Compute one local ISO timestamp for training progress records
def training_timestamp() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


# Format one command-line tool completion message
def format_done_message(tool_name: str, elapsed_seconds: float) -> str:
    return f"[{tool_name}: done in {elapsed_seconds:.1f} sec]"

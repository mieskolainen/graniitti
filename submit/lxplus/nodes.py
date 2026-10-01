# Start Ray nodes inside CERN batch allocations
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import hashlib
import json
import os
import pathlib
import re
import shutil
import socket
import stat
import subprocess
import time
from collections.abc import Iterator
from contextlib import contextmanager, suppress

from core.io.files import check_dependencies, ensure_dir

from submit.lxplus import common


# Compute the original CERN schedd recorded in the steering job ClassAd
def condor_schedd() -> str:
    schedd = os.environ.get("RAY_LXPLUS_SCHEDD", "").strip()
    if not schedd:
        job_ad = pathlib.Path(os.environ.get("_CONDOR_JOB_AD", ""))
        if not job_ad.is_file():
            raise RuntimeError("Ray lxplus steering must run as a submitted HTCondor job")
        match = None
        for line in job_ad.read_text(encoding="utf-8").splitlines():
            match = re.match(r'^\s*GlobalJobId\s*=\s*"([^#"]+)#', line)
            if match is not None:
                break
        if match is None:
            raise RuntimeError(f"HTCondor job ad has no GlobalJobId: {job_ad}")
        schedd = match.group(1)
    if re.fullmatch(r"[A-Za-z0-9_.-]+", schedd) is None:
        raise ValueError(f"Invalid HTCondor schedd name: {schedd!r}")
    return schedd


# Compute the submitted steering cluster ID from its HTCondor job ad
def condor_cluster_id() -> int:
    job_ad = pathlib.Path(os.environ.get("_CONDOR_JOB_AD", ""))
    if not job_ad.is_file():
        raise RuntimeError("Ray lxplus steering job ad is missing")
    for line in job_ad.read_text(encoding="utf-8").splitlines():
        match = re.fullmatch(r"\s*ClusterId\s*=\s*([0-9]+)\s*", line)
        if match is not None:
            return int(match.group(1))
    raise RuntimeError(f"HTCondor job ad has no ClusterId: {job_ad}")


# Compute the worker array process identity from its local HTCondor job ad
def condor_job_identity() -> dict[str, int | str]:
    job_ad = pathlib.Path(os.environ.get("_CONDOR_JOB_AD", ""))
    if not job_ad.is_file():
        raise RuntimeError("Ray lxplus worker job ad is missing")
    values = {}
    for line in job_ad.read_text(encoding="utf-8").splitlines():
        match = re.fullmatch(r"\s*(ClusterId|ProcId)\s*=\s*([0-9]+)\s*", line)
        if match is not None:
            values[match.group(1)] = int(match.group(2))
    if set(values) != {"ClusterId", "ProcId"}:
        raise RuntimeError(f"Ray lxplus worker job identity is incomplete: {job_ad}")
    return {"cluster_id": values["ClusterId"], "job_id": f"{values['ClusterId']}.{values['ProcId']}",
        "proc_id": values["ProcId"], }


# Compute the explicit HTCondor address options for the original CERN schedd
def condor_target() -> tuple[str, list[str]]:
    schedd = condor_schedd()
    address = os.environ.get("RAY_LXPLUS_SCHEDD_ADDR", "").strip()
    if re.fullmatch(r"<[^<>\s]+>", address) is None:
        raise ValueError(f"Invalid HTCondor schedd address: {address!r}")
    return schedd, ["-addr", address]


# Compute the exact physical worker name from the local HTCondor machine ad
def condor_machine() -> str:
    machine_ad = pathlib.Path(os.environ.get("_CONDOR_MACHINE_AD", ""))
    if not machine_ad.is_file():
        raise RuntimeError(f"HTCondor machine ad is missing: {machine_ad}")
    machine = ""
    for line in machine_ad.read_text(encoding="utf-8").splitlines():
        match = re.fullmatch(r'\s*Machine\s*=\s*"([^"]+)"\s*', line)
        if match is not None:
            machine = match.group(1)
            break
    if re.fullmatch(r"[A-Za-z0-9_.-]+", machine) is None:
        raise RuntimeError(f"HTCondor machine ad has no valid Machine: {machine_ad}")
    return machine


# Compute one HTCondor requirement for physical worker placement
def machine_requirements(machines: list[str]) -> str:
    unique = list(dict.fromkeys(machines))
    if not unique:
        raise ValueError("Ray lxplus worker placement needs at least one excluded machine")
    for machine in unique:
        if re.fullmatch(r"[A-Za-z0-9_.-]+", machine) is None:
            raise ValueError(f"Invalid HTCondor machine name: {machine!r}")
    return " && ".join(f'(TARGET.Machine =!= "{machine}")' for machine in unique)


# Check that the submitted steering process can write to the original schedd
def validate_condor_write(address: str) -> None:
    completed = subprocess.run(
        ["condor_ping", "-address", address, "-quiet", "WRITE"], check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"Cannot submit HTCondor allocations for Ray workers to the original CERN schedd: {detail}")


# Compute one currently unused CERN lxbatch Dask port
def dask_scheduler_port() -> int:
    node_ip = socket.gethostbyname(socket.getfqdn())
    for port in range(10000, 10010):
        handle = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            handle.bind((node_ip, port))
        except OSError:
            pass
        else:
            return port
        finally:
            handle.close()
    raise RuntimeError("Dask scheduler needs one free port in the CERN lxbatch range 10000:10009")


# Open one host shared Ray port coordination file
def open_port_file(path: pathlib.Path) -> int:
    flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags, 0o666)
    metadata = os.fstat(descriptor)
    if metadata.st_uid == os.geteuid() and stat.S_IMODE(metadata.st_mode) != 0o666:
        os.fchmod(descriptor, 0o666)
    return descriptor


# Compute worker ports including spare capacity for Python runtime changes
def ray_worker_ports(cpus: int) -> int:
    # Allow idle runtime workers and actor replacement to overlap
    return max(32, int(cpus) + 16)


# Compute the safe Ray allocation count sharing one physical host port window
def ray_host_slots(cpus: int) -> int:
    cpus = int(cpus)
    if cpus <= 0:
        raise ValueError("Ray lxplus CPUs per allocation must be positive")
    worker_ports = ray_worker_ports(cpus)
    worker_total = 7 + worker_ports
    head_total = 9 + worker_ports
    if head_total > 91:
        return 0
    return min(91 // worker_total, 1 + (91 - head_total) // worker_total)


# Compute one short path spelling which resolves to separate local storage on each host
def ray_cluster_tmp(schedd: str, cluster_id: int) -> pathlib.Path:
    token = hashlib.sha256(f"{schedd}#{int(cluster_id)}".encode()).hexdigest()[:12]
    return pathlib.Path("/tmp") / f"irt-{token}"


# Create or validate one private node-local parent for the Ray session tree
def prepare_ray_tmp(path: pathlib.Path) -> None:
    with suppress(FileExistsError):
        ensure_dir(path, exist_ok=False, mode=0o700)
    metadata = os.lstat(path)
    if not stat.S_ISDIR(metadata.st_mode) or metadata.st_uid != os.getuid():
        raise RuntimeError(f"Ray temporary path is not a private owned directory: {path}")
    if stat.S_IMODE(metadata.st_mode) != 0o700:
        path.chmod(0o700)


# Compute the routable address already used by this Dask allocation
def dask_node_ip() -> str:
    import ipaddress
    from urllib.parse import urlsplit

    from distributed import get_worker

    address = str(get_worker().address)
    host = urlsplit(address).hostname
    if host is None:
        raise RuntimeError(f"Dask worker has no network address: {address!r}")
    try:
        node_ip = socket.gethostbyname(host)
        parsed = ipaddress.IPv4Address(node_ip)
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"Dask worker address is not reachable by Ray: {address!r}") from exc
    if parsed.is_link_local or parsed.is_loopback or parsed.is_multicast or parsed.is_unspecified:
        raise RuntimeError(f"Dask worker address is not routable by Ray: {node_ip}")
    return node_ip


# Reserve host ports across Condor PID and mount namespaces
@contextmanager
def reserve_ray_ports(node_ip: str, count: int, cluster: str) -> Iterator[tuple[list[int], bool]]:
    import fcntl

    root = pathlib.Path(os.environ.get("RAY_LXPLUS_PORT_DIR", "/tmp"))
    ensure_dir(root)
    lock = open_port_file(root / "icetune-ray-ports.lock")
    held = []
    try:
        fcntl.flock(lock, fcntl.LOCK_EX)
        ports = []
        for port in range(10010, 10101):
            descriptor = open_port_file(root / f"icetune-ray-port-{port}.lock")
            lease = socket.socket(socket.AF_UNIX, socket.SOCK_DGRAM)
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                # Abstract sockets share the host network namespace even with private /tmp
                lease.bind(f"\0icetune-ray-port-{node_ip}-{port}")
                with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as handle:
                    handle.bind(("", port))
            except OSError:
                lease.close()
                os.close(descriptor)
                continue
            ports.append(port)
            held.extend((descriptor, lease.detach()))
            if len(ports) == count:
                break
        if len(ports) != count:
            raise RuntimeError(f"Ray needs {count} free reserved ports in CERN range 10010:10100, "
                f"but only {len(ports)} are available on {node_ip}")
        monitor = open_port_file(root / f"icetune-ray-monitor-{cluster}.lock")
        try:
            fcntl.flock(monitor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            os.close(monitor)
            first = False
        else:
            held.append(monitor)
            first = True
        yield ports, first
    except BaseException:
        for descriptor in held:
            os.close(descriptor)
        raise
    finally:
        os.close(lock)


# Bound Ray memory by the requested allocation and the detected cgroup limit
def ray_memory(requested: str) -> dict[str, int]:
    from dask.utils import parse_bytes
    from distributed.system import MEMORY_LIMIT

    requested = parse_bytes(requested)
    total = min(requested, MEMORY_LIMIT)
    if total < 512 * 1024**2:
        raise ValueError(f"Ray allocation memory is too small: {total} bytes")
    store = max(80 * 1024**2, min(total // 8, 1024**3))
    common.write_log(f"Ray allocation memory: requested={requested} detected={MEMORY_LIMIT} "
        f"budget={total} object_store={store} bytes")
    return {"memory": total * 3 // 4, "object_store_memory": store}


# Read bounded Ray process logs even after the processes have exited
def ray_log_tails(root: pathlib.Path) -> dict[str, str]:
    logs = {}
    for name in ("gcs_server.err", "gcs_server.out", "raylet.err", "raylet.out", "dashboard.err", "dashboard.log",
        "runtime_env_agent.err", "runtime_env_agent.log", ):
        try:
            with (root / name).open("rb") as handle:
                handle.seek(max(0, os.fstat(handle.fileno()).st_size - 32768))
                logs[name] = handle.read(32768).decode("utf-8", errors="replace")
        except FileNotFoundError:
            continue
        except OSError as exc:
            logs[name] = str(exc)
    return logs


# Include Ray process logs in the exception returned before allocation cleanup
def ray_startup_detail(command: list[str], ray_tmp: pathlib.Path, stdout, stderr) -> str:
    root = ray_tmp / "session_latest" / "logs"
    if not root.is_dir():
        sessions = sorted(ray_tmp.glob("session_*/logs"))
        if not sessions:
            # Ray workers can inherit the head's session path from GCS
            sessions = sorted((ray_tmp.parent.parent / "icetune_head" / "ray").glob("session_*/logs"))
        if sessions:
            root = sessions[-1]
    output = [value.decode("utf-8", errors="replace") if isinstance(value, bytes) else str(value or "")
        for value in (stdout, stderr)]
    return "\n".join(
        (*output, json.dumps({"command": command, "ray_log_dir": str(root), "logs": ray_log_tails(root)}, indent=2)))


# Keep unusable software outside the Ray trial pool
class WorkerRuntimeError(RuntimeError):
    """Worker dependencies or shared libraries are unavailable"""


# Verify external files and load the simulator libraries before starting Ray
def check_worker_runtime(dependencies, check):
    try:
        check_dependencies(dependencies or {})
        if check:
            environment = dict(check["env"])
            environment.update({key: os.environ[key] for key in ("HOME", "USER", "LOGNAME", "TMPDIR") if key in os.environ})
            if os.environ.get("_CONDOR_SCRATCH_DIR"):
                environment["TMPDIR"] = str(common.condor_scratch_dir())
            result = subprocess.run(check["command"], env=environment, capture_output=True, text=True,
                                    timeout=check["timeout"], check=False)
            if result.returncode:
                raise RuntimeError(f"Library check exited {result.returncode}:\n{result.stdout}\n{result.stderr}")
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
        raise WorkerRuntimeError(f"Worker runtime unavailable on {socket.getfqdn()}: {exc}") from exc


# Start one Ray daemon inside a dask-lxplus worker slot
def start_ray_node(*, role: str, cpus: int, temp_dir: str, gpus: int = 0, head_address: str | None = None,
    resource_marker: str | None = None, conda_prefix: str | None = None, ray_bin: str | None = None,
    request_memory: str | None = None, dependencies: dict | None = None, worker_check: dict | None = None, ) -> dict:
    import socket

    if role not in {"head", "worker"}:
        raise ValueError(f"Unknown Ray node role {role!r}")
    if role == "worker" and not head_address:
        raise ValueError("Ray worker requires the head address")
    if role == "worker" and re.fullmatch(r"icetune_worker_[0-9]+_[0-9]+", str(resource_marker or "")) is None:
        raise ValueError("Ray worker requires its Condor resource marker")
    if role == "head" and resource_marker is not None:
        raise ValueError("Ray head uses the fixed icetune_head resource marker")
    if int(cpus) < 0 or (role == "worker" and int(cpus) == 0) or int(gpus) < 0:
        raise ValueError("Ray workers require positive CPUs and all node resources must be non-negative")
    if role == "worker":
        check_worker_runtime(dependencies, worker_check)
    node_ip = dask_node_ip()
    service_count = 9 if role == "head" else 7
    worker_count = ray_worker_ports(cpus)
    memory_key = "RAY_HEAD_REQUEST_MEMORY" if role == "head" else "RAY_REQUEST_MEMORY"
    memory = ray_memory(request_memory or os.environ[memory_key])
    requested_tmp = pathlib.Path(temp_dir)
    if not requested_tmp.is_absolute() or re.fullmatch(r"irt-[0-9a-f]{12}", requested_tmp.name) is None:
        raise ValueError(f"Invalid Ray lxplus temporary path: {requested_tmp}")
    cluster_token = requested_tmp.name.removeprefix("irt-")
    node_name = "icetune_head" if role == "head" else str(resource_marker)
    node_tmp = requested_tmp / node_name
    ray_tmp = node_tmp / "ray"
    spill_dir = common.condor_scratch_dir() / "ray-spill" / cluster_token

    with reserve_ray_ports(node_ip, service_count + worker_count, cluster_token) as reservation:
        ports, first_cluster_node = reservation
        services = ports[:service_count]
        worker_ports = ports[service_count:]
        prepare_ray_tmp(requested_tmp)
        prepare_ray_tmp(node_tmp)
        ensure_dir(spill_dir)
        environment = os.environ.copy()
        environment.update(common.condor_scratch_environment())
        environment.pop("RAY_ADDRESS", None)
        environment.pop("PYTHONHOME", None)
        environment.pop("PYTHONPATH", None)
        prefix_value = conda_prefix or environment.get("ICETUNE_CONDA_PREFIX")
        if not prefix_value:
            raise ValueError("Ray node startup received no Graniitti Conda prefix")
        prefix_path = pathlib.Path(prefix_value)
        environment["PATH"] = f"{prefix_path / 'bin'}:{environment.get('PATH', '')}"
        environment["CONDA_DEFAULT_ENV"] = str(prefix_path)
        environment["CONDA_PREFIX"] = str(prefix_path)
        environment["ICETUNE_CONDA_PREFIX"] = str(prefix_path)
        environment["RAY_raylet_start_wait_time_s"] = environment.get("RAY_RAYLET_START_TIMEOUT_S", "600")
        environment["RAY_TMPDIR"] = str(node_tmp)
        requested_ray = ray_bin or environment.get("ICETUNE_RAY_BIN") or "ray"
        ray_executable = shutil.which(str(requested_ray), path=environment["PATH"])
        if ray_executable is None:
            raise ValueError(f'Graniitti Ray executable is unavailable on this allocation: "{requested_ray}"')
        environment["ICETUNE_RAY_BIN"] = ray_executable

        command = [ray_executable, "start", f"--node-ip-address={node_ip}", f"--num-cpus={int(cpus)}",
            f"--memory={memory['memory']}", f"--object-store-memory={memory['object_store_memory']}",
            f"--node-manager-port={services[0]}", f"--object-manager-port={services[1]}",
            f"--runtime-env-agent-port={services[2]}", f"--dashboard-agent-listen-port={services[3]}",
            f"--dashboard-agent-grpc-port={services[4]}", f"--metrics-export-port={services[5]}",
            f"--object-spilling-directory={spill_dir}",
            f"--worker-port-list={','.join(str(port) for port in worker_ports)}", ]
        if not first_cluster_node:
            command.append("--include-log-monitor=False")
        # Reserve head GPUs for the fitting process, outside Ray trial scheduling
        command.append(f"--num-gpus={0 if role == 'head' else int(gpus)}")
        result = {"head_address": head_address, "hostname": socket.getfqdn(), "node_ip": node_ip,
            "ray_tmp": str(ray_tmp), "resource_marker": "icetune_head" if role == "head" else resource_marker,
            "role": role, "spill_dir": str(spill_dir), }
        if role == "head":
            result["head_address"] = f"{node_ip}:{services[6]}"
            command.extend(("--head", f"--port={services[6]}", f"--ray-client-server-port={services[7]}",
                    f"--dashboard-port={services[8]}", "--include-dashboard=False",
                    f"--resources={json.dumps({'icetune_head': 1}, separators=(',', ':'))}", f"--temp-dir={ray_tmp}", ))
        else:
            command.extend((f"--ray-client-server-port={services[6]}", f"--address={head_address}",
                    f"--resources={json.dumps({resource_marker: 1}, separators=(',', ':'))}", ))

        timeout = float(os.environ.get("RAY_LXPLUS_NODE_START_TIMEOUT_S", "1200"))
        started_at = time.monotonic()
        try:
            completed = subprocess.run(
                command, check=False, capture_output=True, env=environment, text=True, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            detail = ray_startup_detail(command, ray_tmp, exc.stdout, exc.stderr)
            raise RuntimeError(f"Ray {role} startup timed out after {timeout:.0f} seconds on "
                f"{node_ip}: {detail or 'no launcher output'}") from exc
        if completed.returncode != 0:
            detail = ray_startup_detail(command, ray_tmp, completed.stdout, completed.stderr)
            raise RuntimeError(f"Ray {role} startup failed on {node_ip}: {detail}")
        result["startup_seconds"] = time.monotonic() - started_at
    return result


# Start one Ray worker and return a serializable startup result
def start_ray_worker_node(*, cpus: int, temp_dir: str, gpus: int, head_address: str, resource_marker: str,
    conda_prefix: str | None = None, ray_bin: str | None = None, request_memory: str | None = None,
    dependencies: dict | None = None, worker_check: dict | None = None, ) -> dict:
    try:
        node = start_ray_node(role="worker", cpus=cpus, temp_dir=temp_dir, gpus=gpus, head_address=head_address,
            resource_marker=resource_marker, conda_prefix=conda_prefix, ray_bin=ray_bin,
            request_memory=request_memory, dependencies=dependencies, worker_check=worker_check, )
    except Exception as exc:
        return {"error": str(exc), "status": "failed", "retry": not isinstance(exc, WorkerRuntimeError)}
    return {"node": node, "status": "started"}


# Collect bounded Ray log tails from processes belonging to this allocation
def ray_worker_logs(marker: str) -> dict:
    import psutil

    logs = {}
    for process in psutil.process_iter(["cmdline"]):
        command = process.info["cmdline"] or []
        if not command or not command[0].endswith("/raylet") or marker not in " ".join(command):
            continue
        for argument in command:
            if not argument.startswith("--log_dir="):
                continue
            root = pathlib.Path(argument.split("=", 1)[1])
            logs.update(ray_log_tails(root))
    return {"hostname": condor_machine(), "job": condor_job_identity(), "logs": logs}

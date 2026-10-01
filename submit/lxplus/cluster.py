# Acquire and maintain CERN Dask allocations for Ray
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import json
import os
import pathlib
import re
import subprocess
import sys
import time

from core.io.files import ensure_dir

from submit.lxplus import common, nodes


# Compute a CERN job class which owns one fixed size HTCondor array
def cern_array_job(job_cls):
    class CernArrayJob(job_cls):
        # Initialize one array submission with unique Dask worker names
        def __init__(self, *args, array_size: int, array_name: str, head_cpu: int, head_gpu: int,
            head_memory: str, worker_gpu: int, worker_memory: str, name=None, **kwargs, ):
            self.array_size = int(array_size)
            if self.array_size < 0:
                raise ValueError("Ray lxplus worker array size must be non-negative")
            if re.fullmatch(r"[A-Za-z0-9_.-]+", array_name) is None:
                raise ValueError(f"Invalid Ray lxplus array name: {array_name!r}")
            # Apply after dask-lxplus appends its default range so Dask cannot take Ray ports
            kwargs["worker_extra_args"] = [*(kwargs.get("worker_extra_args") or []), "--worker-port", "10000:10009"]
            super().__init__(*args, name=f"{name}-$(ProcId)", **kwargs)
            self.job_header_dict.pop("batch_name", None)
            for key in tuple(self.job_header_dict):
                normalized = str(key).lstrip("+").lower()
                if normalized in {"my.sendcredential", "sendcredential", "should_transfer_files",
                    "transfer_output_files", "when_to_transfer_output", }:
                    self.job_header_dict.pop(key)
            self.job_header_dict["JobBatchName"] = f'"{array_name}"'
            # Use Condor scratch so the forwarded credential never lands in the checkout
            self.job_header_dict["MY.SendCredential"] = "True"
            self.job_header_dict["should_transfer_files"] = "YES"
            self.job_header_dict["when_to_transfer_output"] = "ON_EXIT"
            self.job_header_dict["transfer_output_files"] = '""'
            self.head_cpu = int(head_cpu)
            self.head_gpu = int(head_gpu)
            self.head_memory = head_memory
            self.worker_gpu = int(worker_gpu)
            self.worker_memory = worker_memory

        # Render one dedicated head followed by the Ray worker array
        def job_script(self) -> str:
            script = super().job_script()
            if re.search(r"(?im)^\s*\+?(?:MY\.)?SendCredential\s*=\s*True\s*$", script) is None:
                raise RuntimeError("dask-lxplus worker description must forward the CERN credential")
            transfer_modes = [value.strip('"').upper()
                for value in re.findall(r"(?im)^\s*should_transfer_files\s*=\s*([^\s#]+)", script)]
            if not transfer_modes or any(value != "YES" for value in transfer_modes):
                raise RuntimeError("dask-lxplus worker description must use an isolated transfer sandbox")
            returned = re.findall(r"(?im)^\s*transfer_output_files\s*=\s*([^\r\n#]*)", script)
            if not returned or any(value.strip().strip('"') for value in returned):
                raise RuntimeError("dask-lxplus worker description must return no scratch files")
            match = re.search(r"(?m)^Queue\s*$", script)
            if match is None or script[match.end() :].strip():
                raise RuntimeError("dask-lxplus worker description has no terminal Queue statement")
            prefix = script[: match.start()]
            cpu = re.search(r"(?im)^[ \t]*(request_?cpus[ \t]*=)(.*)$", prefix)
            if cpu is None:
                raise RuntimeError("dask-lxplus worker description has no CPU request")
            worker_cpu_line = f"{cpu.group(1)}{cpu.group(2)}\n"
            prefix = f"{prefix[: cpu.start()]}{cpu.group(1)} {self.head_cpu}{prefix[cpu.end() :]}"
            memory = re.search(r"(?im)^\s*(request_?memory\s*=).*$", prefix)
            if memory is None:
                raise RuntimeError("dask-lxplus worker description has no memory request")
            head = f"{memory.group(1)} {self.head_memory}"
            prefix = f"{prefix[: memory.start()]}{head}{prefix[memory.end() :]}"
            # Set both queue groups explicitly so head GPUs cannot carry into workers
            gpu = re.search(r"(?im)^[ \t]*request_?gpus[ \t]*=.*$", prefix)
            head_gpu = f"request_gpus = {self.head_gpu}"
            if gpu is None:
                prefix += f"{head_gpu}\n"
            else:
                prefix = f"{prefix[: gpu.start()]}{head_gpu}{prefix[gpu.end() :]}"
            worker_gpu_line = f"request_gpus = {self.worker_gpu}\n"
            head_queue = f"{prefix}Queue 1\n"
            if self.array_size == 0:
                return head_queue
            return (f"{head_queue}{worker_cpu_line}{memory.group(1)} {self.worker_memory}\n" f"{worker_gpu_line}"
                "on_exit_remove = False\non_exit_hold = NumJobStarts >= 3\n" f"Queue {self.array_size}\n")

        # Compute the cluster ID so one condor_rm removes the complete array
        def _job_id_from_submit_output(self, output: str) -> str:
            job_id = super()._job_id_from_submit_output(output)
            match = re.fullmatch(r"(\d+)(?:\.0)?", job_id)
            if match is None:
                raise RuntimeError(f"Invalid dask-lxplus array cluster ID: {job_id!r}")
            cluster_id = match.group(1)
            if self.array_size:
                common.write_log(f"Ray lxplus head: {cluster_id}.0; worker array: {cluster_id}.1 "
                    f"through {cluster_id}.{self.array_size}")
            else:
                common.write_log(f"Ray lxplus head: {cluster_id}.0; worker array is empty")
            return cluster_id

    return CernArrayJob


# Compute newly connected and placement validated Dask slots
def new_ray_slots(*, client, known: set[str], steering_machine: str) -> list[str]:
    workers = sorted(set(client.scheduler_info()["workers"]) - known)
    if not workers:
        return []
    placement = client.run(nodes.condor_machine, workers=workers)
    if set(placement) != set(workers):
        raise RuntimeError("dask-lxplus did not return every HTCondor allocation placement")
    machines = [placement[worker] for worker in workers]
    if steering_machine in machines:
        raise RuntimeError(f"dask-lxplus placed a worker allocation on the steering machine {steering_machine}")
    for worker, machine in zip(workers, machines, strict=True):
        common.write_log(
            f"HTCondor allocation ready: Dask control={worker}, host={machine}; Ray has not started on it yet")
    return workers


# Submit one Ray slot array without waiting for Condor matchmaking
def submit_ray_slots(*, cluster, workers: int, cpus: int, steering_machine: str) -> None:
    try:
        spec = cluster.new_spec
        directives = spec["options"]["job_extra_directives"]
    except (AttributeError, KeyError, TypeError) as exc:
        raise RuntimeError("dask-lxplus worker specification has no HTCondor directives") from exc
    if not isinstance(spec, dict) or not isinstance(directives, dict):
        raise RuntimeError("dask-lxplus worker specification must use dictionaries")

    # Map every array process name back to its single Dask job owner
    group = [f"-{index}" for index in range(workers + 1)]
    if spec.get("group") not in (None, group):
        raise RuntimeError("dask-lxplus worker specification has an incompatible process group")
    spec["group"] = group

    # Reconcile inside one Dask loop call so the array submission failure propagates
    async def submit_slots() -> None:
        cluster.scale(jobs=1)
        await cluster._correct_state()

    if nodes.ray_host_slots(cpus) <= 0:
        raise ValueError("Ray lxplus CPU request exceeds the CERN port range")
    directives["requirements"] = nodes.machine_requirements([steering_machine])
    directives.pop("rank", None)
    common.write_log(
        f"Submitting one dedicated Ray head and an array of {workers} Ray workers with fragmented-slot matchmaking")
    cluster.sync(submit_slots)
    common.write_log(f"Submitted one Ray head and {workers} Ray worker array processes")


# Wait only for the first submitted HTCondor worker
def wait_ray_slots(*, client, timeout: float, steering_machine: str) -> list[str]:
    common.write_log(f"Waiting up to {timeout:.0f} seconds for the first HTCondor allocation; "
        "Ray starts after the dedicated head allocation is ready")
    deadline = time.monotonic() + timeout
    known = set()
    workers = []
    client.wait_for_workers(1, timeout=timeout)
    while time.monotonic() < deadline:
        workers.extend(new_ray_slots(client=client, known=known, steering_machine=steering_machine))
        known.update(workers)
        identities = client.run(nodes.condor_job_identity, workers=workers)
        heads = [worker for worker in workers if identities[worker]["proc_id"] == 0]
        if heads:
            head = heads[0]
            return [head, *(worker for worker in workers if worker != head)]
        time.sleep(min(1.0, max(0.0, deadline - time.monotonic())))
    raise TimeoutError(f"Ray head HTCondor allocation process 0 did not start within {timeout:.0f} seconds")


# Compute the unique Ray custom resource owned by one Condor worker process
def ray_worker_marker(identity: dict) -> str:
    try:
        cluster_id = int(identity["cluster_id"])
        proc_id = int(identity["proc_id"])
        job_id = str(identity["job_id"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Ray worker Condor identity is invalid") from exc
    if cluster_id <= 0 or proc_id <= 0 or job_id != f"{cluster_id}.{proc_id}":
        raise ValueError(f"Ray worker Condor identity is invalid: {identity}")
    return f"icetune_worker_{cluster_id}_{proc_id}"


# Submit independent Ray worker starts without waiting for their completion
def submit_ray_worker_starts(*, client, workers: list[str], markers: dict[str, str], starting: dict, cpus: int,
    temp_dir: str, gpus: int, head_address: str, conda_prefix: str | None = None, ray_bin: str | None = None,
    dependencies: dict | None = None, worker_check: dict | None = None,
) -> dict[str, dict]:
    if not workers:
        return {}
    missing_markers = sorted(set(workers) - set(markers))
    if missing_markers:
        raise RuntimeError(f"Ray workers have no Condor resource marker: {missing_markers}")
    results = {}
    for worker in workers:
        if worker in starting:
            raise RuntimeError(f"Ray worker startup is already active: {worker}")
        try:
            starting[worker] = client.submit(nodes.start_ray_worker_node, cpus=cpus, temp_dir=temp_dir, gpus=gpus,
                head_address=head_address, resource_marker=markers[worker], conda_prefix=conda_prefix,
                ray_bin=ray_bin, request_memory=os.environ.get("RAY_REQUEST_MEMORY"), dependencies=dependencies,
                worker_check=worker_check, workers=[worker],
                allow_other_workers=False, pure=False, )
        except Exception as exc:
            results[worker] = {"error": f"Dask startup submission failed: {exc}", "status": "failed"}
    return results


# Collect only completed Ray worker starts so slow allocations do not block admission
def collect_ray_worker_starts(*, starting: dict) -> dict[str, dict]:
    results = {}
    for worker, future in tuple(starting.items()):
        try:
            if not future.done():
                continue
        except Exception as exc:
            results[worker] = {"error": f"Dask startup status failed: {exc}", "status": "failed"}
            starting.pop(worker, None)
            continue
        starting.pop(worker, None)
        try:
            results[worker] = future.result()
        except Exception as exc:
            results[worker] = {"error": f"Dask startup task failed: {exc}", "status": "failed"}
    return results


# Validate and log one completed Ray worker startup outcome
def normalize_ray_worker_start(*, worker: str, result: dict | None, head_address: str) -> dict:
    if not isinstance(result, dict):
        result = {"error": "no matching startup result", "status": "failed"}
    status = str(result.get("status", "failed"))
    if status == "started":
        node = result.get("node")
        if not isinstance(node, dict) or not node.get("hostname") or not node.get("spill_dir"):
            result = {"error": "Ray worker startup returned incomplete node metadata", "status": "failed"}
            status = "failed"
    if status == "started":
        common.write_log(f"Ray worker launcher completed: host={result['node']['hostname']}, "
            f"Dask control={worker}, Ray head={head_address}, " f"local spill={result['node']['spill_dir']}, "
            f"startup={float(result['node'].get('startup_seconds', 0.0)):.1f}s")
    else:
        result["status"] = "failed"
        result.setdefault("error", "unknown error")
        common.write_log(f"Ray worker {worker} startup failed: {result['error']}", file=sys.stderr)
    return result


# Apply completed Ray worker starts to retry and live node state
def apply_ray_worker_starts(*, results: dict[str, dict], pending: list[str], seen: set[str], joined: set[str],
    failed: set[str], attempts: dict[str, int], errors: dict[str, str], head_address: str, retry_limit: int, ) -> None:
    for worker, raw_result in results.items():
        result = normalize_ray_worker_start(worker=worker, result=raw_result, head_address=head_address)
        status = str(result.get("status"))
        if status == "started":
            seen.add(worker)
            joined.add(worker)
            continue
        count = attempts.get(worker, 0) + 1
        attempts[worker] = count
        errors[worker] = str(result.get("error", "unknown error"))
        if count <= retry_limit and result.get("retry", True):
            pending.append(worker)
            common.write_log(
                f"Ray worker {worker} will retry startup {count}/{retry_limit}: {errors[worker]}", file=sys.stderr)
            continue
        seen.add(worker)
        failed.add(worker)
        common.write_log(f"Ray worker {worker} rejected before joining Ray: {errors[worker]}", file=sys.stderr)


# Compute the next steering delay without rounding fractional launch intervals
def ray_start_delay(*, pending, starting, max_starting, launch_state) -> float:
    if pending and len(starting) < max_starting:
        return min(1.0, max(0.0, launch_state["next_at"] - time.monotonic()))
    return 1.0


# Advance one centrally paced wave of independent Ray worker launches
def advance_ray_workers(*, client, head_key: str, pending: list[str], seen: set[str], joined: set[str],
    failed: set[str], starting: dict, attempts: dict[str, int], errors: dict[str, str], markers: dict[str, str],
    launch_state: dict[str, float], cpus: int, gpus: int, temp_dir: str, head_address: str, steering_machine: str,
    start_interval_s: float, max_starting: int, retry_limit: int, conda_prefix: str | None = None,
    ray_bin: str | None = None, dependencies: dict | None = None, worker_check: dict | None = None, ) -> None:
    if start_interval_s <= 0.0:
        raise ValueError("Ray worker start interval must be positive")
    if max_starting <= 0:
        raise ValueError("Ray simultaneous worker start limit must be positive")
    live = set(client.scheduler_info()["workers"])
    if head_key not in live:
        raise RuntimeError("The dedicated Ray head HTCondor allocation disappeared")
    lost = ((seen | set(starting)) - {head_key}) - live
    for worker in sorted(lost):
        future = starting.pop(worker, None)
        if future is not None:
            try:
                future.cancel()
            except Exception as exc:
                common.write_log(f"Ray worker startup cancellation deferred for {worker}: {exc}", file=sys.stderr)
        seen.discard(worker)
        joined.discard(worker)
        failed.discard(worker)
        attempts.pop(worker, None)
        errors.pop(worker, None)
        markers.pop(worker, None)
        common.write_log(f"Ray worker allocation disappeared and may rejoin: Dask control={worker}", file=sys.stderr)
    pending[:] = list(dict.fromkeys(worker for worker in pending if worker in live))
    state = dict(pending=pending, seen=seen, joined=joined, failed=failed, attempts=attempts, errors=errors,
        head_address=head_address, retry_limit=retry_limit, )
    joined_before = len(joined)
    completed = collect_ray_worker_starts(starting=starting)
    apply_ray_worker_starts(results=completed, **state)
    pending.extend(
        new_ray_slots(client=client, known=seen | set(pending) | set(starting), steering_machine=steering_machine))
    now = time.monotonic()
    next_at = launch_state.setdefault("next_at", now)
    selected = []
    if pending and len(starting) < max_starting and now >= next_at:
        selected.append(pending.pop(0))
        launch_state["next_at"] = max(next_at + start_interval_s, now + start_interval_s)
    unknown = [worker for worker in selected if worker not in markers]
    if unknown:
        batch_identities = client.run(nodes.condor_job_identity, workers=unknown)
        if set(batch_identities) != set(unknown):
            raise RuntimeError("Ray worker returned incomplete Condor metadata")
        for worker in unknown:
            markers[worker] = ray_worker_marker(batch_identities[worker])
    immediate = submit_ray_worker_starts(client=client, workers=selected, markers=markers, starting=starting,
        cpus=cpus, temp_dir=temp_dir, gpus=gpus, head_address=head_address, conda_prefix=conda_prefix,
        ray_bin=ray_bin, dependencies=dependencies, worker_check=worker_check, )
    apply_ray_worker_starts(results=immediate, **state)
    completed = collect_ray_worker_starts(starting=starting)
    apply_ray_worker_starts(results=completed, **state)
    if selected:
        common.write_log(f"Ray startup queue: {len(joined) - 1} joined, {len(starting)} starting, "
            f"{len(pending)} waiting, limit={max_starting}")
    if len(joined) > joined_before:
        common.write_log(f"Ray launch progress: dedicated head and {len(joined) - 1} workers started")


# Requeue only a failed worker process and exclude its previous physical host
def replace_ray_worker(*, marker: str, machine: str, exhausted: bool) -> None:
    match = re.fullmatch(r"icetune_worker_([1-9][0-9]*)_([1-9][0-9]*)", marker)
    if match is None:
        raise ValueError(f"Invalid Ray worker replacement marker: {marker}")
    job = ".".join(match.groups())
    _, options = nodes.condor_target()
    timeout_s = float(os.environ.get("RAY_CONDOR_COMMAND_TIMEOUT_S", "30"))
    completed = subprocess.run(["condor_q", *options, job, "-json", "-attributes", "Requirements"], check=True,
        capture_output=True, text=True, timeout=timeout_s, )
    ads = json.loads(completed.stdout)
    if len(ads) != 1 or not ads[0].get("Requirements"):
        raise RuntimeError(f"Cannot read worker requirements for {job}")
    requirements = f"({ads[0]['Requirements']}) && {nodes.machine_requirements([machine])}"
    subprocess.run(["condor_hold", *options, "-reason", "icetune Ray worker failed", job], check=True,
        capture_output=True, text=True, timeout=timeout_s, )
    if exhausted:
        return
    subprocess.run(["condor_qedit", *options, job, "Requirements", requirements], check=True, capture_output=True,
        text=True, timeout=timeout_s, )
    subprocess.run(["condor_release", *options, job], check=True, capture_output=True, text=True, timeout=timeout_s)


# Process each unhealthy Ray node once and bound allocation replacements across rejoins
def recover_ray_workers(*, client, resources: dict | None, markers: dict, recovery: dict, log_dir: pathlib.Path
) -> None:
    failure_limit = int(os.environ.get("RAY_WORKER_FAILURE_LIMIT", "2"))
    replacement_limit = int(os.environ.get("RAY_WORKER_REPLACEMENT_LIMIT", "2"))
    failures = dict(((resources or {}).get("resources") or {}).get("unhealthy_workers", {}))
    failures.update({marker: record for marker, record in recovery.items() if not record["done"]})
    for marker, failure in failures.items():
        if int(failure.get("count", 0)) < failure_limit:
            continue
        node_id = str(failure["node_id"])
        previous = recovery.get(marker, {})
        if previous.get("done") and (previous.get("node_id") == node_id or previous.get("exhausted")):
            continue
        try:
            if previous and not previous["done"]:
                record = previous
            else:
                worker = next((worker for worker, value in markers.items() if value == marker), None)
                if worker is None:
                    continue
                detail = client.run(nodes.ray_worker_logs, marker=marker, workers=[worker])[worker]
                count = int(previous.get("attempts", 0)) + 1
                record = {**detail, **failure, "attempts": count, "exhausted": count > replacement_limit, "done": False}
                ensure_dir(log_dir)
                (log_dir / f"{marker}-{node_id}.json").write_text(json.dumps(record, indent=2) + "\n")
                recovery[marker] = record
            replace_ray_worker(marker=marker, machine=record["hostname"], exhausted=record["exhausted"])
            record["done"] = True
            common.write_log(f"Ray worker {marker}: {'held after retry limit' if record['exhausted'] else 'requeued'} "
                f"on failure {record['attempts']}, previous host={record['hostname']}")
        except Exception as exc:
            common.write_log(f"Ray worker replacement deferred for {marker}: {exc}", file=sys.stderr)

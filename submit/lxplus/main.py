# Coordinate the CERN icetune cluster and output return
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import hashlib
import json
import os
import pathlib
import re
import signal
import socket
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

from core.io.files import ensure_dir

from submit import campaign
from submit import runtime as staging
from submit.lxplus import cluster as allocations
from submit.lxplus import common, nodes, outputs


# Write one atomic process record for the head local icetune driver
def write_runtime_pid(*, environment: dict[str, str], pid: int) -> pathlib.Path:
    birth = common.process_birth(pid)
    if birth is None or os.getpgid(pid) != pid:
        raise RuntimeError(f"icetune driver process group is invalid: {pid}")
    path = common.runtime_pid_path(environment)
    payload = (json.dumps({"birth": birth, "pid": pid}, sort_keys=True) + "\n").encode()
    outputs.atomic_write(path, payload)
    return path


# Request a scoped interrupt of the head local icetune driver process group
def request_runtime_stop(*, environment: dict[str, str]) -> dict:
    path = common.runtime_pid_path(environment)
    if not path.is_file():
        return {"requested": False, "reason": "driver process record is absent"}
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
        pid = int(record["pid"])
        birth = str(record["birth"])
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Invalid icetune driver process record: {path}") from exc
    if pid <= 0 or common.process_birth(pid) != birth:
        return {"requested": False, "reason": "driver process already exited"}
    try:
        if os.getpgid(pid) != pid:
            raise RuntimeError(f"icetune driver process group changed: {pid}")
        os.killpg(pid, signal.SIGINT)
    except ProcessLookupError:
        return {"requested": False, "reason": "driver process already exited"}
    return {"pid": pid, "requested": True}


# Run icetune once on the Ray head slot and compute its compact output tree
def run_icetune_runtime(*, environment: dict[str, str]) -> dict:
    import os
    import pathlib
    import subprocess
    import tarfile
    from urllib.request import urlopen

    work_root = common.runtime_work_root(environment)
    ensure_dir(work_root)
    # Keep this record after completion or failure because Dask may recompute lost tasks
    started = work_root / "icetune-driver.started"
    try:
        with started.open("x", encoding="ascii") as handle:
            handle.write(f"{os.getpid()}\n")
            handle.flush()
            os.fsync(handle.fileno())
    except FileExistsError as exc:
        raise RuntimeError(
            f"icetune driver already started in this head allocation: {started}. Refusing to repeat the campaign"
        ) from exc
    with urlopen(environment["RAY_RUNTIME_URI"]) as source, tarfile.open(fileobj=source, mode="r|gz") as archive:
        archive.extractall(work_root, filter="data")
    runtime = work_root / "runtime"
    for relative in ("eikonal", "figs", "output", "runs/icetune", "sudakov", "tmp", "vgrid"):
        ensure_dir(runtime / relative)
    storage = runtime / "runs" / "icetune"
    job_environment = os.environ.copy()
    # Keep the GPU assignment of this allocation, not the steering process
    job_environment.update({key: value for key, value in environment.items() if key not in (
                "CUDA_VISIBLE_DEVICES", "NVIDIA_VISIBLE_DEVICES", "GPU_DEVICE_ORDINAL", "ROCR_VISIBLE_DEVICES",
                "HIP_VISIBLE_DEVICES", )})
    job_environment.update(common.condor_scratch_environment())
    job_environment.pop("GRANIITTI_ENV", None)
    job_environment.pop("PYTHONHOME", None)
    job_environment.pop("PYTHONPATH", None)
    conda_prefix = pathlib.Path(job_environment["ICETUNE_CONDA_PREFIX"])
    job_environment["PATH"] = f"{conda_prefix / 'bin'}:{job_environment.get('PATH', '')}"
    job_environment.update({"RAY_REPO_DIR": str(runtime), "RAY_STORAGE_PATH": str(storage), "RAY_UPLOAD_RUNTIME": "1"})
    process = subprocess.Popen(["bash", str(runtime / "submit/shell" / "run_ray.sh")],
        cwd=runtime, env=job_environment, start_new_session=True, )
    pid_path = None
    try:
        pid_path = write_runtime_pid(environment=environment, pid=process.pid)
        returncode = process.wait()
    except BaseException:
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
                process.wait(timeout=10.0)
            except ProcessLookupError:
                pass
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        raise
    finally:
        if pid_path is not None:
            pid_path.unlink(missing_ok=True)

    snapshot = outputs.snapshot_final_outputs(environment=environment)
    return {"checkpoint": snapshot["checkpoint"], "history_figures": snapshot["history_figures"],
        "outputs": snapshot["outputs"], "recovery": snapshot["recovery"], "returncode": returncode, }


# Close the Dask client and allocation cluster without skipping either cleanup
def close_dask_cluster(*, client, cluster, reason: str) -> None:
    for name, resource in (("client", client), ("cluster", cluster)):
        if resource is not None:
            try:
                resource.close()
            except Exception as exc:
                common.write_log(f"Dask {name} cleanup after {reason}: {exc}", file=sys.stderr)


# Stop the Ray driver, return its final state and release every Dask allocation
def graceful_ray_shutdown(*, client, cluster, future, head_key: str | None, environment: dict[str, str],
    shared_root: pathlib.Path, run_name: str, campaign_fingerprint: str, figure_revision: dict | None,
    live_revision: dict[str, tuple[int, int]] | None, timeout_s: float, ) -> None:
    common.write_log("Graceful Ray lxplus shutdown requested")
    deadline = time.monotonic() + max(1.0, float(timeout_s))
    final_returned = False
    if client is not None and future is not None and head_key is not None:
        try:
            if not future.done():
                responses = client.run(request_runtime_stop, environment=environment, workers=[head_key])
                response = responses.get(head_key, {})
                if not response.get("requested"):
                    common.write_log(
                        f"Ray driver interrupt was not needed or available: {response.get('reason', 'unknown reason')}",
                        file=sys.stderr, )
                else:
                    common.write_log(f"Ray driver interrupt sent to process {response['pid']}")
            result = future.result(timeout=max(1.0, float(timeout_s) - 30.0))
            outputs.publish_runtime_result(
                result=result, shared_root=shared_root, run_name=run_name, campaign_fingerprint=campaign_fingerprint)
            final_returned = True
            if isinstance(result, dict) and result.get("checkpoint") is not None:
                common.write_log("Ray driver exited and returned its final checkpoint")
            else:
                common.write_log("Ray driver exited and returned its final outputs")
        except Exception as exc:
            common.write_log(f"Ray driver graceful stop deferred: {exc}", file=sys.stderr)
    if not final_returned and client is not None and head_key is not None:
        try:
            remaining = deadline - time.monotonic()
            if remaining <= 5.0:
                raise TimeoutError("no time remains for the last live output return")
            outputs.return_live_outputs(client=client, head_key=head_key, environment=environment,
                shared_root=shared_root, run_name=run_name, campaign_fingerprint=campaign_fingerprint,
                figure_revision=figure_revision, live_revision=live_revision, timeout_s=min(30.0, remaining - 5.0), )
            common.write_log("Returned the last available live Ray outputs before shutdown")
        except Exception as exc:
            common.write_log(f"Final live Ray output return deferred: {exc}", file=sys.stderr)
        common.write_log(
            "No final checkpoint was created because driver termination and final packing were not confirmed",
            file=sys.stderr, )
    close_dask_cluster(client=client, cluster=cluster, reason="graceful shutdown")
    common.write_log("Ray and Dask allocations closed")


# Restore verified initialization and previous outputs before starting any allocations
def prepare_runtime(environment: dict, stage_path: pathlib.Path) -> None:
    shared_root = campaign.shared_output_dir(environment).resolve()
    run_name, campaign_fingerprint = environment["RUN_NAME"], environment["RAY_CAMPAIGN_FINGERPRINT"]
    driver = campaign.simdriver_name(environment)
    runtime_archive = pathlib.Path(environment["RAY_INIT_RUNTIME_ARCHIVE"])
    runtime_sha256 = environment["RAY_INIT_RUNTIME_SHA256"]
    common.write_log(f"Staging Ray lxplus runtime from {runtime_archive}")
    runtime = staging.extract_lxplus_runtime(
        archive=runtime_archive, target=stage_path / "runtime", checksum=runtime_sha256)
    previous = outputs.stage_previous_outputs(
        shared_root=shared_root, runtime=runtime, run_name=run_name, campaign_fingerprint=campaign_fingerprint)
    if environment.get("PRECOMPUTE", "1") == "1":
        coord_dir = pathlib.Path(environment["RAY_INIT_COORD_DIR"])
        expected_coord = shared_root / "runs" / "icetune" / run_name / "ray" / "init" / campaign_fingerprint[:16]
        if coord_dir.resolve() != expected_coord.resolve():
            raise ValueError(f"Ray INIT coordination path differs from this runtime: {coord_dir}")
        fingerprint = staging.stage_ray_bootstrap(
            runtime=runtime, coord_dir=coord_dir, run_name=run_name, simdriver=driver)
        if environment["DATA_COVARIANCE_MODE"] == "full":
            result_dir = runtime / "runs" / "icetune" / run_name / "results"
            outputs.restore_outputs(payload=outputs.pack_runtime_outputs(figure_dir=None, result_files=[
                        result_dir / name for name in ("data_covariance.npz", "data_covariance.json")],
                    run_name=run_name, ), shared_root=shared_root, run_name=run_name,
                campaign_fingerprint=campaign_fingerprint, )
        common.write_log(f"Staged current verified Ray INIT bootstrap {fingerprint[:16]}")
    common.write_log("Packing Ray lxplus runtime")
    archive = staging.pack_runtime(runtime)
    environment["RAY_RUNTIME_URI"] = staging.publish_ray_runtime(
        archive=archive, directory=shared_root / "runs" / "icetune" / run_name / "ray" / "runtime")
    common.write_log(f"Packed Ray lxplus runtime: {archive.stat().st_size / (1024 * 1024):.1f} MiB")
    if previous:
        common.write_log(f"Ray lxplus restore: {previous[0]}")


# Advance worker admission while returning completed live output snapshots
def monitor(launch: dict, publication: dict, future, *, poll_interval: float, shutdown_timeout: float):
    campaign_figures = outputs.campaign_figure_path(publication["shared_root"], publication["run_name"],
                                                   publication["campaign_fingerprint"])
    with ThreadPoolExecutor(max_workers=1, thread_name_prefix="icetune-return") as publisher:
        next_return_at = time.monotonic()
        recovery = {}
        transfer = None
        while not future.done():
            allocations.advance_ray_workers(**launch)
            now = time.monotonic()
            if transfer is None and now >= next_return_at:
                startup_failures = {launch["markers"][worker]: {"node_id": hashlib.sha256(worker.encode()).hexdigest(),
                        "count": int(publication["environment"]["RAY_WORKER_FAILURE_LIMIT"]), "error": launch["errors"][worker], }
                    for worker in launch["failed"] if worker in launch["markers"]}
                allocations.recover_ray_workers(client=publication["client"],
                    resources={"resources": {"unhealthy_workers": startup_failures}}, markers=launch["markers"],
                    recovery=recovery, log_dir=common.condor_scratch_dir() / "ray-worker-failures", )
                transfer = publisher.submit(outputs.return_live_outputs, **publication, timeout_s=shutdown_timeout)
            if transfer is not None and transfer.done():
                try:
                    publication["figure_revision"], publication["live_revision"], figure_updates, _ray_resources = transfer.result()
                    allocations.recover_ray_workers(client=publication["client"], resources=_ray_resources,
                        markers=launch["markers"], recovery=recovery,
                        log_dir=common.condor_scratch_dir() / "ray-worker-failures", )
                    if "best" in figure_updates:
                        common.write_log(f"Ray lxplus best trial figures updated under {campaign_figures}")
                    if "history" in figure_updates:
                        common.write_log(
                            f"Ray lxplus history figures updated under {campaign_figures}"
                        )
                except Exception as exc:
                    common.write_log(f"Ray lxplus incremental output return deferred: {exc}", file=sys.stderr)
                transfer = None
                next_return_at = time.monotonic() + poll_interval
            time.sleep(allocations.ray_start_delay(
                    pending=launch["pending"], starting=launch["starting"], max_starting=launch["max_starting"], launch_state=launch["launch_state"]))


# Run the existing Ray campaign after dask-lxplus has acquired batch slots
def main() -> int:
    try:
        from dask_lxplus import CernCluster
        from distributed import Client
    except ImportError as exc:
        raise RuntimeError("Ray lxplus mode requires dask-lxplus and distributed from a CERN LCG view") from exc

    run_name = os.environ["RUN_NAME"]
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", run_name) is None:
        raise ValueError(f"Invalid Ray lxplus run name: {run_name!r}")
    workers = int(os.environ["RAY_WORKERS"])
    total_allocations = workers + 1
    cpus = int(os.environ["RAY_WORKER_CPU"])
    gpus = int(os.environ["RAY_WORKER_GPU"])
    max_runtime = int(os.environ["RAY_MAX_RUNTIME_S"])
    timeout = float(os.environ.get("RAY_LXPLUS_STARTUP_TIMEOUT_S", os.environ.get("RAY_STARTUP_TIMEOUT_S", "10800")))
    start_interval = float(os.environ.get("RAY_LXPLUS_START_INTERVAL_S", "1.5"))
    max_starting = int(os.environ.get("RAY_LXPLUS_MAX_STARTING", "16"))
    join_retry_limit = int(os.environ.get("RAY_LXPLUS_JOIN_RETRY_LIMIT", "3"))
    node_start_timeout = float(os.environ.get("RAY_LXPLUS_NODE_START_TIMEOUT_S", "1200"))
    shutdown_timeout = float(os.environ.get("RAY_LXPLUS_SHUTDOWN_TIMEOUT_S", "120"))
    if min(timeout, start_interval, node_start_timeout, shutdown_timeout) <= 0.0:
        raise ValueError("Ray lxplus startup, node startup and shutdown intervals must be positive")
    if max_starting <= 0:
        raise ValueError("Ray lxplus simultaneous worker start limit must be positive")
    if join_retry_limit < 0:
        raise ValueError("Ray lxplus worker startup retry limit must be nonnegative")
    shared_root = campaign.shared_output_dir(os.environ).resolve()
    worker_log_dir = pathlib.Path(os.environ["RAY_SUBMIT_DIR"]) / "logs"
    ensure_dir(worker_log_dir)
    conda_prefix = os.environ["ICETUNE_CONDA_PREFIX"]
    ray_bin = os.environ.get("ICETUNE_RAY_BIN") or str(pathlib.Path(conda_prefix) / "bin" / "ray")
    if nodes.ray_host_slots(cpus) <= 0:
        raise ValueError("Ray lxplus CPU request exceeds the CERN port range")
    if max_runtime > 604800:
        raise ValueError("Ray lxplus maximum runtime is 604800 seconds")
    schedd, condor_options = nodes.condor_target()
    steering_machine = nodes.condor_machine()
    ray_tmp = nodes.ray_cluster_tmp(schedd, nodes.condor_cluster_id())
    campaign_fingerprint = os.environ["RAY_CAMPAIGN_FINGERPRINT"]
    campaign_figures = outputs.campaign_figure_path(shared_root, run_name, campaign_fingerprint)
    address = condor_options[1]
    common.write_log(f"Ray lxplus steering: schedd={schedd} address={address} "
        f"machine={steering_machine} workers={workers} worker_cpu={cpus} "
        f"node-local_ray_parent={ray_tmp} start_interval={start_interval:g}s " f"max_starting={max_starting} "
        f"campaign={campaign_fingerprint[:16]}")
    nodes.validate_condor_write(address)
    common.write_log("Ray lxplus HTCondor submission is ready")

    # Let Condor terminate Ray daemons before removing their allocation
    job_options = {"MY.MaxRuntime": str(max_runtime), "MY.SendCredential": "True", "MY.WantOS": '"el9"',
        "job_max_vacate_time": os.environ["RAY_WORKER_SHUTDOWN_TIMEOUT_S"], "kill_sig": "15",
        "should_transfer_files": "YES", "transfer_output_files": '""', "want_graceful_removal": "True",
        "when_to_transfer_output": "ON_EXIT", }
    array_name = f"icetune-ray-{run_name}-workers"
    cluster_options = {"array_name": array_name, "array_size": workers, "batch_name": array_name,
        "container_runtime": "none", "cores": cpus, "death_timeout": os.environ["RAY_LXPLUS_DEATH_TIMEOUT_S"],
        "disk": os.environ["RAY_REQUEST_DISK"], "job_cls": allocations.cern_array_job(CernCluster.job_cls),
        "job_directives_skip": ["batch_name"], "job_extra_directives": job_options,
        # dask-lxplus forwards the LCG environment validated by STEER
        "job_script_prologue": [], "lcg": True, "log_directory": str(worker_log_dir),
        "memory": os.environ["RAY_REQUEST_MEMORY"], "head_memory": os.environ["RAY_HEAD_REQUEST_MEMORY"],
        "head_cpu": int(os.environ["RAY_HEAD_CPU"]), "head_gpu": int(os.environ["RAY_HEAD_GPU"]),
        "worker_gpu": gpus, "worker_memory": os.environ["RAY_REQUEST_MEMORY"], "nanny": False, "processes": 1,
        "submit_command_extra": list(condor_options), "cancel_command_extra": list(condor_options), }

    environment = os.environ.copy()
    with tempfile.TemporaryDirectory(prefix=f"icetune-ray-{run_name}-", dir=common.condor_scratch_dir()) as stage:
        scheduler_host = socket.gethostname()
        scheduler_port = nodes.dask_scheduler_port()
        cluster_options["scheduler_options"] = {"dashboard_address": None, "host": scheduler_host,
            "port": scheduler_port, }
        cluster = None
        client = None
        future = None
        head_key = None
        publication = {"figure_revision": None, "live_revision": None}
        checkpoint_returned = False
        try:
            prepare_runtime(environment, pathlib.Path(stage))
            definition = json.loads((pathlib.Path(stage) / "runtime" / environment["TUNESETUP"]).read_text())
            dependencies = definition["aux_param_space"].get("runtime_dependencies", {})
            cluster = CernCluster(**cluster_options)
            client = Client(cluster)
            common.write_log(f"Ray lxplus Dask scheduler started: host={scheduler_host} port={scheduler_port}")
            allocations.submit_ray_slots(cluster=cluster, workers=workers, cpus=cpus, steering_machine=steering_machine)
            slots = allocations.wait_ray_slots(client=client, timeout=timeout, steering_machine=steering_machine)
            head_key = slots[0]
            head_future = client.submit(nodes.start_ray_node, role="head", cpus=0, temp_dir=str(ray_tmp), gpus=0,
                conda_prefix=conda_prefix, ray_bin=ray_bin, request_memory=environment["RAY_HEAD_REQUEST_MEMORY"],
                workers=[head_key], allow_other_workers=False, pure=False, )
            head_result = head_future.result()
            if not isinstance(head_result, dict) or not head_result.get("head_address"):
                raise RuntimeError("dask-lxplus did not start the requested Ray head")
            head_address = head_result["head_address"]
            common.write_log(f"Ray head started: host={head_result['hostname']}, "
                f"Dask control={head_key}, Ray address={head_address}, "
                f"host-local session={head_result['ray_tmp']}, " f"local spill={head_result['spill_dir']}, "
                f"startup={float(head_result.get('startup_seconds', 0.0)):.1f}s")
            joined = {head_key}
            seen = {head_key}
            failed = set()
            starting = {}
            attempts = {}
            errors = {}
            markers = {}
            launch_state = {"next_at": time.monotonic()}
            pending = list(slots[1:])
            environment["RAY_ADDRESS"] = head_address
            poll_interval = max(1.0, float(environment.get("RAY_LXPLUS_POLL_INTERVAL_S", "60")))
            common.write_log(f"Changed Ray status, history and figures return every {poll_interval:g} seconds at most")
            launch = dict(client=client, head_key=head_key, pending=pending, seen=seen, joined=joined,
                failed=failed, starting=starting, attempts=attempts, errors=errors, markers=markers,
                launch_state=launch_state, cpus=cpus, gpus=gpus, temp_dir=str(ray_tmp), head_address=head_address,
                steering_machine=steering_machine, start_interval_s=start_interval, max_starting=max_starting,
                retry_limit=join_retry_limit, conda_prefix=conda_prefix, ray_bin=ray_bin, dependencies=dependencies,
                worker_check=definition["aux_param_space"].get("worker_check"), )
            first_worker_deadline = time.monotonic() + timeout
            while len(joined) == 1:
                allocations.advance_ray_workers(**launch)
                if len(joined) > 1:
                    break
                if len(seen) == total_allocations and not pending and not starting:
                    detail = next(iter(errors.values()), "no worker launcher succeeded")
                    raise RuntimeError(f"No Ray worker allocation could start: {detail}")
                if time.monotonic() >= first_worker_deadline:
                    raise TimeoutError(f"No Ray worker allocation started within {timeout:.0f} seconds")
                time.sleep(allocations.ray_start_delay(
                        pending=pending, starting=starting, max_starting=max_starting, launch_state=launch_state))
            environment["RAY_WAIT_WORKERS"] = "1"
            common.write_log(f"Starting icetune after the first Ray worker launcher completed; "
                f"{len(joined) - 1}/{workers} requested workers started and later " "HTCondor allocations may join")
            future = client.submit(run_icetune_runtime, environment=environment,
                workers=[head_key], allow_other_workers=False, pure=False, retries=0, )
            publication = dict(client=client, head_key=head_key, environment=environment, shared_root=shared_root,
                               run_name=run_name, campaign_fingerprint=campaign_fingerprint,
                               figure_revision=None, live_revision=None)
            monitor(launch, publication, future, poll_interval=poll_interval, shutdown_timeout=shutdown_timeout)
            if len(seen) == total_allocations and failed:
                allocation_word = "allocation" if len(failed) == 1 else "allocations"
                common.write_log(f"Ray worker launchers completed for {len(joined) - 1}/{workers} workers. "
                    f"Startup failed for {len(failed)} requested {allocation_word}", file=sys.stderr, )
            elif len(joined) == total_allocations:
                common.write_log(f"Ray launchers completed for all {workers} requested workers")
            result = future.result()
            returncode = outputs.publish_runtime_result(
                result=result, shared_root=shared_root, run_name=run_name, campaign_fingerprint=campaign_fingerprint)
            checkpoint_returned = bool(isinstance(result, dict) and result.get("checkpoint") is not None)
        except KeyboardInterrupt:
            graceful_ray_shutdown(client=client, cluster=cluster, future=future, head_key=head_key,
                environment=environment, shared_root=shared_root, run_name=run_name,
                campaign_fingerprint=campaign_fingerprint, figure_revision=publication["figure_revision"],
                live_revision=publication["live_revision"], timeout_s=shutdown_timeout, )
            return 0
        except BaseException:
            close_dask_cluster(client=client, cluster=cluster, reason="failure")
            raise
        else:
            try:
                client.close()
            finally:
                cluster.close()
    if checkpoint_returned:
        common.write_log(
            f"Ray lxplus checkpoint returned to {outputs.checkpoint_path(shared_root, run_name, campaign_fingerprint)}"
        )
    common.write_log(f"Ray lxplus figures returned to {campaign_figures}")
    return returncode


# Serialize steering functions independently of a shared checkout on Dask allocations
def register_worker_modules():
    import cloudpickle

    for name, module in tuple(sys.modules.items()):
        if name.startswith("submit") or name == "core.io.files":
            cloudpickle.register_pickle_by_value(module)


if __name__ == "__main__":
    register_worker_modules()
    try:
        raise SystemExit(main())
    except (KeyError, OSError, RuntimeError, ValueError) as exc:
        common.write_log(f"[icetune lxplus] {exc}", file=sys.stderr)
        raise SystemExit(common.steer_exit_code(exc)) from exc

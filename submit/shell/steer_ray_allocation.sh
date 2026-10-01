#!/bin/bash
#
# Start one Ray worker allocation per physical host, plus a dedicated head

set -euo pipefail

SCHEDULER="${1:-}"
RAY_REPO_DIR="${RAY_REPO_DIR:-${PWD}}"
RAY_WORKERS="${RAY_WORKERS:-1}"

[[ "${SCHEDULER}" =~ ^(pbs|slurm)$ ]] || {
    echo "Usage: ${0} <pbs|slurm>" >&2
    exit 64
}

RAY_REPO_DIR="$(cd "${RAY_REPO_DIR}" && pwd -P)"
RAY_LAUNCHER="$(dirname "${BASH_SOURCE[0]}")/launch_ray_node.sh"
RUNNER="$(dirname "${BASH_SOURCE[0]}")/run_ray.sh"
source "$(dirname "${BASH_SOURCE[0]}")/ray_runtime.sh"
prepare_icetune_ray_environment "${RAY_REPO_DIR}"
prepare_icetune_ray_scratch "${RAY_REPO_DIR}" "${RUN_NAME:-cluster}"

# Compute the unique hosts allocated by the active scheduler
allocated_hosts() {
    if [[ "${SCHEDULER}" == "pbs" ]]; then
        [[ -n "${PBS_NODEFILE:-}" && -f "${PBS_NODEFILE}" ]] || {
            echo "[${0}] PBS_NODEFILE is missing" >&2
            return 64
        }
        sort -u "${PBS_NODEFILE}"
    else
        [[ -n "${SLURM_NODELIST:-}" ]] || {
            echo "[${0}] SLURM_NODELIST is missing" >&2
            return 64
        }
        scontrol show hostname "${SLURM_NODELIST}" | sort -u
    fi
}

# Run one short Ray control command on an allocated worker host
run_remote() {
    local host="$1"
    shift
    if [[ "${SCHEDULER}" == "slurm" ]]; then
        exec srun \
            --exclusive \
            --nodes=1 \
            --ntasks=1 \
            --cpus-per-task="${RAY_WORKER_CPU}" \
            --nodelist="${host}" \
            "$@"
    elif command -v pbsdsh >/dev/null 2>&1; then
        exec pbsdsh -h "${host}" "$@"
    else
        exec ssh -o BatchMode=yes -o StrictHostKeyChecking=no "${host}" "$@"
    fi
}

mapfile -t hosts < <(allocated_hosts)
required_hosts=$((RAY_WORKERS + 1))
(( ${#hosts[@]} >= required_hosts )) || {
    echo "[${0}] Requested one Ray head and ${RAY_WORKERS} workers but scheduler allocated ${#hosts[@]} hosts" >&2
    exit 64
}

head_short="$(hostname -s)"
worker_hosts=()
for host in "${hosts[@]}"; do
    [[ "${host%%.*}" == "${head_short}" ]] && continue
    worker_hosts+=("${host}")
    (( ${#worker_hosts[@]} >= RAY_WORKERS )) && break
done
(( ${#worker_hosts[@]} == RAY_WORKERS )) || {
    echo "[${0}] Could not select ${RAY_WORKERS} worker hosts beside the dedicated head" >&2
    exit 64
}

head_ip="$(hostname -I | awk '{print $1}')"
[[ -n "${head_ip}" ]] || {
    echo "[${0}] Could not determine Ray head IP" >&2
    exit 64
}

remote_step_pids=()
head_pid=""

# Stop the allocation-local Ray daemons when tuning exits
cleanup_ray() {
    for pid in "${remote_step_pids[@]}"; do
        kill "${pid}" >/dev/null 2>&1 || true
        wait "${pid}" >/dev/null 2>&1 || true
    done
    if [[ -n "${head_pid}" ]]; then
        kill "${head_pid}" >/dev/null 2>&1 || true
        wait "${head_pid}" >/dev/null 2>&1 || true
    fi
}
trap cleanup_ray EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

echo "[${0}] Starting Ray head on ${head_short} (${head_ip})"
head_address_file="${RAY_TMP_DIR}/ray_current_cluster"
if [[ -f "${head_address_file}" ]]; then
    mv "${head_address_file}" "${head_address_file}._old-${BASHPID}"
fi
RAY_NODE_IP="${head_ip}" bash "${RAY_LAUNCHER}" head "" 1 &
head_pid=$!
export RAY_ADDRESS="${head_ip}:${RAY_PORT:-60010}"

# Wait for this allocation's head before connecting workers or the fitting process
deadline=$((SECONDS + RAY_STARTUP_TIMEOUT_S))
until [[ -s "${head_address_file}" ]] && \
    timeout "$((deadline - SECONDS > 0 ? deadline - SECONDS : 1))" \
    ray health-check --address="${RAY_ADDRESS}" >/dev/null 2>&1; do
    if ! kill -0 "${head_pid}" 2>/dev/null || (( SECONDS >= deadline )); then
        echo "[${0}] Ray head allocation did not become ready" >&2
        exit 1
    fi
    sleep 1
done

for host in "${worker_hosts[@]}"; do
    echo "[${0}] Starting Ray worker on ${host}"
    run_remote "${host}" bash "${RAY_LAUNCHER}" worker "${head_ip}" 1 &
    remote_step_pids+=("$!")
done

export RAY_WAIT_WORKERS="${RAY_WORKERS}"
bash "${RUNNER}"

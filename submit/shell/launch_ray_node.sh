#!/bin/bash
#
# Start one Ray head or worker daemon with the resources of its batch allocation

set -euo pipefail

ROLE="${1:-}"
HEAD_ADDRESS="${2:-}"
BLOCK="${3:-0}"
RAY_REPO_DIR="${RAY_REPO_DIR:-${PWD}}"
RAY_PORT="${RAY_PORT:-60010}"
RAY_WORKER_CPU="${RAY_WORKER_CPU:-1}"
RAY_WORKER_GPU="${RAY_WORKER_GPU:-0}"

[[ "${ROLE}" =~ ^(head|worker)$ ]] || {
    echo "Usage: ${0} <head|worker> [head-address]" >&2
    exit 64
}
if [[ "${ROLE}" == "worker" && -z "${HEAD_ADDRESS}" ]]; then
    echo "[${0}] Worker role requires the Ray head address" >&2
    exit 64
fi

RAY_REPO_DIR="$(cd "${RAY_REPO_DIR}" && pwd -P)"
source "$(dirname "${BASH_SOURCE[0]}")/ray_runtime.sh"
prepare_icetune_ray_environment "${RAY_REPO_DIR}"
prepare_icetune_ray_scratch "${RAY_REPO_DIR}" "${RUN_NAME:-cluster}"

# Keep Ray state on node-local scratch and use a routable address
RAY_NODE_IP="${RAY_NODE_IP:-$(hostname -I | awk '{print $1}')}"
[[ -n "${RAY_NODE_IP}" ]] || {
    echo "[${0}] Could not determine a non-loopback node address" >&2
    exit 64
}

if [[ "${ROLE}" == "head" ]]; then
    resource_args=(
        "--num-cpus=0"
        "--num-gpus=0"
        '--resources={"icetune_head":1}'
    )
else
    resource_args=("--num-cpus=${RAY_WORKER_CPU}" "--num-gpus=${RAY_WORKER_GPU}")
fi
resource_args+=("--object-spilling-directory=${RAY_SPILL_DIR}")
block_args=()
if [[ "${BLOCK}" == "1" ]]; then
    block_args+=("--block")
fi

# Add one explicitly configured Ray network option
append_network_arg() {
    local variable="$1"
    local option="$2"
    local value=""

    if [[ -v "${variable}" ]]; then
        value="${!variable}"
    fi

    if [[ -n "${value}" && "${value}" != "0" ]]; then
        network_args+=("${option}=${value}")
    fi
}

network_args=()
append_network_arg RAY_NODE_MANAGER_PORT --node-manager-port
append_network_arg RAY_OBJECT_MANAGER_PORT --object-manager-port
append_network_arg RAY_RUNTIME_ENV_AGENT_PORT --runtime-env-agent-port
append_network_arg RAY_DASHBOARD_AGENT_LISTEN_PORT --dashboard-agent-listen-port
append_network_arg RAY_DASHBOARD_AGENT_GRPC_PORT --dashboard-agent-grpc-port
append_network_arg RAY_METRICS_EXPORT_PORT --metrics-export-port
if [[ -n "${RAY_WORKER_PORT_LIST:-}" ]]; then
    network_args+=("--worker-port-list=${RAY_WORKER_PORT_LIST}")
elif [[ "${RAY_MIN_WORKER_PORT:-0}" != "0" || "${RAY_MAX_WORKER_PORT:-0}" != "0" ]]; then
    [[ -n "${RAY_MIN_WORKER_PORT:-}" && -n "${RAY_MAX_WORKER_PORT:-}" ]] || {
        echo "[${0}] Set both RAY_MIN_WORKER_PORT and RAY_MAX_WORKER_PORT" >&2
        exit 64
    }
    network_args+=(
        "--min-worker-port=${RAY_MIN_WORKER_PORT}"
        "--max-worker-port=${RAY_MAX_WORKER_PORT}"
    )
fi

# Ray --block owns its child processes and cleans up only this allocation on exit
if [[ "${ROLE}" == "head" ]]; then
    head_args=()
    if [[ "${RAY_CLIENT_SERVER_PORT:-0}" != "0" ]]; then
        head_args+=("--ray-client-server-port=${RAY_CLIENT_SERVER_PORT}")
    fi
    echo "[${0}] Starting Ray head at ${RAY_NODE_IP}:${RAY_PORT}"
    exec ray start \
        --head \
        --node-ip-address="${RAY_NODE_IP}" \
        --port="${RAY_PORT}" \
        --include-dashboard=False \
        --temp-dir="${RAY_TMP_DIR}" \
        "${resource_args[@]}" \
        "${network_args[@]}" \
        "${head_args[@]}" \
        "${block_args[@]}"
else
    if [[ "${HEAD_ADDRESS}" == *:* ]]; then
        RAY_CLUSTER_ADDRESS="${HEAD_ADDRESS}"
    else
        RAY_CLUSTER_ADDRESS="${HEAD_ADDRESS}:${RAY_PORT}"
    fi
    echo "[${0}] Joining Ray head ${RAY_CLUSTER_ADDRESS} from ${RAY_NODE_IP}"
    exec ray start \
        --address="${RAY_CLUSTER_ADDRESS}" \
        --node-ip-address="${RAY_NODE_IP}" \
        "${resource_args[@]}" \
        "${network_args[@]}" \
        "${block_args[@]}"
fi

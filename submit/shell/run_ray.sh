#!/bin/bash
#
# Run one catalog-resolved icetune Ray campaign

set -euo pipefail

RAY_REPO_DIR="${RAY_REPO_DIR:-${PWD}}"
RAY_REPO_DIR="$(cd "${RAY_REPO_DIR}" && pwd -P)"
source "$(dirname "${BASH_SOURCE[0]}")/ray_runtime.sh"
prepare_icetune_ray_environment "${RAY_REPO_DIR}"

RAY_ADDRESS="${RAY_ADDRESS:-local}"
RAY_UPLOAD_RUNTIME="${RAY_UPLOAD_RUNTIME:-0}"
RAY_STARTUP_TIMEOUT_S="${RAY_STARTUP_TIMEOUT_S:-10800}"
RAY_WAIT_WORKERS="${RAY_WAIT_WORKERS:-${RAY_WORKERS}}"
RAY_STORAGE_PATH="${RAY_STORAGE_PATH:-${RAY_REPO_DIR}/runs/icetune}"

case "${RAY_REPO_DIR}" in
    /eos|/eos/*)
        echo "[${0}] Ray requires a non-EOS repository path: ${RAY_REPO_DIR}" >&2
        exit 64
        ;;
esac
[[ "${RAY_STORAGE_PATH}" == /* ]] || {
    echo "[${0}] RAY_STORAGE_PATH must be absolute: ${RAY_STORAGE_PATH}" >&2
    exit 64
}
case "${RAY_STORAGE_PATH}" in
    /eos|/eos/*)
        echo "[${0}] Ray Tune storage must not use EOS: ${RAY_STORAGE_PATH}" >&2
        exit 64
        ;;
esac
if [[ "${RAY_ADDRESS,,}" == "local" && "${RAY_WORKERS}" != "1" ]]; then
    echo "[${0}] Local Ray requires RAY_WORKERS=1" >&2
    exit 64
fi

prepare_icetune_ray_scratch "${RAY_REPO_DIR}" "${RUN_NAME}"
mkdir -p "${RAY_STORAGE_PATH}"

echo "[${0}] campaign=${RUN_NAME} driver=${SIMDRIVER} address=${RAY_ADDRESS}"
export RAY_ADDRESS RAY_STORAGE_PATH RAY_UPLOAD_RUNTIME RAY_WAIT_WORKERS RAY_STARTUP_TIMEOUT_S
exec python -m submit.runtime --phase fit --cdir "${RAY_REPO_DIR}"

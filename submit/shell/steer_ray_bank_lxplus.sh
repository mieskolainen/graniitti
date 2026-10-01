#!/bin/bash
#
# Evaluate an ampfit basis shard or finalize its complete sample on lxplus

set -euo pipefail

# Exit with the permanent-configuration status
die_config() {
    echo "[${0}] Configuration error: $*" >&2
    exit 64
}

for required in RUN_NAME RAY_INIT_RUNTIME_ARCHIVE RAY_INIT_RUNTIME_SHA256 RAY_INIT_JOBS; do
    [[ -n "${!required:-}" ]] || die_config "${required} is required"
done
[[ "${RAY_INIT_RUNTIME_SHA256}" =~ ^[0-9a-f]{64}$ ]] ||
    die_config "RAY_INIT_RUNTIME_SHA256 must be a lowercase SHA-256"
[[ "${RAY_INIT_JOBS}" =~ ^[0-9]+$ ]] || die_config "RAY_INIT_JOBS must be a positive integer"

SCRATCH_PARENT="${LOCAL_SCRATCH_ROOT:-${TMPDIR:-${_CONDOR_SCRATCH_DIR:-}}}"
[[ -n "${SCRATCH_PARENT}" ]] || die_config "HTCondor scratch directory is missing"
SCRATCH_RESOLVED="$(readlink -f -- "${SCRATCH_PARENT}")" ||
    die_config "scratch directory does not exist: ${SCRATCH_PARENT}"
case "${SCRATCH_RESOLVED}" in
    /eos|/eos/*) die_config "scratch directory resolves to EOS: ${SCRATCH_RESOLVED}" ;;
esac

# DAG arguments select the sample and its shard or finalization
SHARD_INDEX="${2:-}"
[[ "${SHARD_INDEX}" == finalize || "${SHARD_INDEX}" =~ ^[0-9]+$ ]] || die_config "a Condor shard index or finalize is required"

SAFE_RUN_NAME="${RUN_NAME//[^A-Za-z0-9_.-]/_}"
JOB_SCRATCH="$(mktemp -d "${SCRATCH_PARENT%/}/icetune-${SAFE_RUN_NAME}-ray-bank.XXXXXX")"
LOCAL_GRDEV="${JOB_SCRATCH}/grdev"
mkdir -p \
    "${LOCAL_GRDEV}" \
    "${JOB_SCRATCH}/tmp" \
    "${JOB_SCRATCH}/cache/matplotlib" \
    "${JOB_SCRATCH}/cache/numba" \
    "${JOB_SCRATCH}/cache/pycache"
export TMPDIR="${JOB_SCRATCH}/tmp"
export XDG_CACHE_HOME="${JOB_SCRATCH}/cache"
export MPLCONFIGDIR="${JOB_SCRATCH}/cache/matplotlib"
export NUMBA_CACHE_DIR="${JOB_SCRATCH}/cache/numba"
export PYTHONPYCACHEPREFIX="${JOB_SCRATCH}/cache/pycache"

echo "${RAY_INIT_RUNTIME_SHA256}  ${RAY_INIT_RUNTIME_ARCHIVE}" |
    sha256sum --check --status || die_config "runtime checksum mismatch"
tar --zstd -xf "${RAY_INIT_RUNTIME_ARCHIVE}" -C "${LOCAL_GRDEV}"
cd "${LOCAL_GRDEV}"
# Conda activation hooks may read unset optional variables
set +u
source install/setconda_lxplus.sh
set -u
export PYTHONPATH="${PWD}/python${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p eikonal sudakov vgrid output tmp

if [[ "${SHARD_INDEX}" == finalize ]]; then
    exec python -m submit.runtime --phase finalize --cdir "${LOCAL_GRDEV}" --bank_index "$1" --bank_jobs "$3"
fi

echo "[${0}] phase=bank run=${RUN_NAME} bank=$1 shard=${SHARD_INDEX}/$3 scratch=${JOB_SCRATCH} coord=${RAY_INIT_COORD_DIR}"
exec python -m submit.runtime --phase bank --cdir "${LOCAL_GRDEV}" \
    --bank_index "$1" --bank_shard "${SHARD_INDEX}" --bank_jobs "$3"

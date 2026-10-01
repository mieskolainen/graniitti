#!/bin/bash
#
# Shared portable environment and scratch setup for icetune Ray processes

source "${RAY_REPO_DIR:-${PWD}}/tests/environment.sh"

# Activate an optional environment and validate the Ray Python runtime
prepare_icetune_ray_environment() {
    local repo_dir="$1"
    local ray_python=""
    local restore_nounset=0

    export RAY_raylet_start_wait_time_s="${RAY_RAYLET_START_TIMEOUT_S:-600}"
    cd "${repo_dir}"
    if [[ "$-" == *u* ]]; then
        restore_nounset=1
        set +u
    fi
    if [[ -n "${ICETUNE_ENV_SETUP:-}" ]]; then
        source "${ICETUNE_ENV_SETUP}"
    else
        activate_conda_environment "${ICETUNE_CONDA_PREFIX:-${ICETUNE_CONDA_ENV:?Set runtime.conda_env in the campaign steering}}" || return
    fi

    source install/setenv.sh
    export PYTHONPATH="${PWD}/python${PYTHONPATH:+:${PYTHONPATH}}"
    if [[ -n "${ICETUNE_CONDA_PREFIX:-${CONDA_PREFIX:-}}" &&
          -x "${ICETUNE_CONDA_PREFIX:-${CONDA_PREFIX}}/bin/python" ]]; then
        ray_python="${ICETUNE_CONDA_PREFIX:-${CONDA_PREFIX}}/bin/python"
        export PATH="${ICETUNE_CONDA_PREFIX:-${CONDA_PREFIX}}/bin:${PATH}"
    else
        ray_python="$(command -v python)"
    fi
    if (( restore_nounset )); then
        set -u
    fi
    env -u PYTHONHOME -u PYTHONPATH "${ray_python}" -c "import ray" >/dev/null 2>&1 || {
        echo "[${0}] Python cannot import Ray; activate the ICETUNE environment or set ICETUNE_ENV_SETUP/ICETUNE_CONDA_ENV" >&2
        return 64
    }
}

# Put Ray sessions and process caches on node-local scratch
prepare_icetune_ray_scratch() {
    local repo_dir="$1"
    local run_name="${2:-cluster}"
    local cache_dir runtime_id safe_run scratch_id scratch_parent session_parent

    scratch_parent="${RAY_TMP_PARENT:-${TMPDIR:-/tmp}}"
    session_parent="${RAY_SESSION_TMP_PARENT:-/tmp}"
    runtime_id="${SLURM_JOB_ID:-${PBS_JOBID:-${CONDOR_CLUSTER_ID:-manual}}}"
    safe_run="${run_name//[^A-Za-z0-9_.-]/_}"
    read -r scratch_id _ < <(
        printf "%s" "${safe_run}:${runtime_id}" | cksum
    )
    # The spelling is common across nodes, but each path is host-local
    export RAY_TMPDIR="${session_parent%/}/r-${scratch_id:0:5}"
    export RAY_TMP_DIR="${RAY_TMPDIR}/ray"
    cache_dir="${scratch_parent%/}/icetune-cache/${safe_run}/${runtime_id}/$(hostname -s)"
    export RAY_SPILL_DIR="${scratch_parent%/}/icetune-spill/${safe_run}/${runtime_id}/$(hostname -s)"
    export XDG_CACHE_HOME="${cache_dir}"
    export MPLCONFIGDIR="${cache_dir}/matplotlib"
    export NUMBA_CACHE_DIR="${cache_dir}/numba"
    export PYTHONPYCACHEPREFIX="${cache_dir}/pycache"
    export RAY_ENABLE_UV_RUN_RUNTIME_ENV=0
    [[ ! -L "${RAY_TMPDIR}" ]] || {
        echo "[${0}] Ray session scratch must not be a symbolic link: ${RAY_TMPDIR}" >&2
        return 64
    }
    mkdir -p \
        "${MPLCONFIGDIR}" \
        "${NUMBA_CACHE_DIR}" \
        "${PYTHONPYCACHEPREFIX}" \
        "${RAY_SPILL_DIR}"
    install -d -m 700 "${RAY_TMPDIR}"
    [[ "$(stat -c %u "${RAY_TMPDIR}")" == "$(id -u)" ]] || {
        echo "[${0}] Ray session scratch is not owned by this job: ${RAY_TMPDIR}" >&2
        return 64
    }
}

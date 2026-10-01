#!/bin/bash
# Load the CERN LCG view and start the icetune Ray lxplus steering
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

set -euo pipefail

[[ -n "${_CONDOR_SCRATCH_DIR:-}" && -d "${_CONDOR_SCRATCH_DIR}" ]] || {
    echo "[${0}] HTCondor scratch directory is missing: ${_CONDOR_SCRATCH_DIR:-unset}" >&2
    exit 64
}
export TMPDIR="${_CONDOR_SCRATCH_DIR}"
export TMP="${_CONDOR_SCRATCH_DIR}"
export TEMP="${_CONDOR_SCRATCH_DIR}"
export XDG_CACHE_HOME="${_CONDOR_SCRATCH_DIR}/cache"
export MPLCONFIGDIR="${XDG_CACHE_HOME}/matplotlib"
export NUMBA_CACHE_DIR="${XDG_CACHE_HOME}/numba"
export PYTHONPYCACHEPREFIX="${XDG_CACHE_HOME}/pycache"
export RAY_raylet_start_wait_time_s="${RAY_RAYLET_START_TIMEOUT_S:-600}"
mkdir -p "${MPLCONFIGDIR}" "${NUMBA_CACHE_DIR}" "${PYTHONPYCACHEPREFIX}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAY_REPO_DIR="${RAY_REPO_DIR:-${PWD}}"
ICETUNE_LCG_VIEW="${ICETUNE_LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_110/x86_64-el9-gcc13-opt/setup.sh}"

[[ -f "${ICETUNE_LCG_VIEW}" ]] || {
    echo "[${0}] CERN LCG view is missing: ${ICETUNE_LCG_VIEW}" >&2
    exit 64
}

cd "${RAY_REPO_DIR}"
set +u
source tests/environment.sh
export GRANIITTI_ENV=
source install/setconda_lxplus.sh
ICETUNE_CONDA_PREFIX="${CONDA_PREFIX}"
CONDA_PYTHON="${ICETUNE_CONDA_PREFIX}/bin/python"
ICETUNE_RAY_BIN="${ICETUNE_CONDA_PREFIX}/bin/ray"
if [[ ! -x "${CONDA_PYTHON}" || ! -x "${ICETUNE_RAY_BIN}" ]]; then
    echo "[${0}] Graniitti conda environment has no Python or Ray executable: ${ICETUNE_CONDA_PREFIX}" >&2
    exit 64
fi
export ICETUNE_CONDA_PREFIX ICETUNE_RAY_BIN

source "${ICETUNE_LCG_VIEW}"
LCG_PYTHON="$(command -v python)"
LCG_PYTHONHOME="${PYTHONHOME:-}"
LCG_PYTHONPATH="${PYTHONPATH:-}"
set -u

# Require dask-lxplus and its interpreter from the configured CVMFS view
PYTHONHOME="${LCG_PYTHONHOME}" PYTHONPATH="${LCG_PYTHONPATH}" "${LCG_PYTHON}" -c '
import dask_lxplus
import distributed
import sys

if not sys.executable.startswith("/cvmfs/"):
    raise RuntimeError(f"Dask Python must come from CVMFS, got {sys.executable}")
for module in (dask_lxplus, distributed):
    if not module.__file__.startswith("/cvmfs/"):
        raise RuntimeError(f"{module.__name__} must come from CVMFS, got {module.__file__}")
'
echo "[${0}] CVMFS Dask environment is ready"

# Require Ray from the active Graniitti conda environment
env -u PYTHONHOME -u PYTHONPATH "${CONDA_PYTHON}" -c '
import ray
import sys

if not ray.__file__.startswith(f"{sys.prefix}/"):
    raise RuntimeError(f"Ray must come from {sys.prefix}, got {ray.__file__}")
'
echo "[${0}] Graniitti Ray environment is ready"

export PYTHONPATH="${RAY_REPO_DIR}:${RAY_REPO_DIR}/python/src:${LCG_PYTHONPATH}"
exec "${LCG_PYTHON}" -m submit.lxplus.main

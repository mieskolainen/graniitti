# Activate the selected Conda environment for CERN launchers
# 
# Source from the repository root, using Conda on PATH
# or set first
# export CONDA_EXE=/path/to/bin/conda
# 
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

source tests/environment.sh || return
_icetune_conda_env="${ICETUNE_CONDA_PREFIX:-${ICETUNE_CONDA_ENV:-graniitti}}"
if [[ -n "${CONDA_EXE:-}" ]]; then
    if [[ ! -x "${CONDA_EXE}" ]]; then
        echo "CONDA_EXE is not executable: ${CONDA_EXE}" >&2
        return 64
    fi
    _graniitti_conda_hook="$("${CONDA_EXE}" shell.bash hook)" || return
    eval "${_graniitti_conda_hook}" || return
    unset _graniitti_conda_hook
elif ! command -v conda >/dev/null 2>&1 && ! conda_environment_is_active "${_icetune_conda_env}"; then
    echo "Conda is unavailable: initialize Conda or set CONDA_EXE to its executable" >&2
    return 64
fi
activate_conda_environment "${_icetune_conda_env}" || return
unset _icetune_conda_env
source install/setenv.sh

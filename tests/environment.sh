#!/usr/bin/env bash

# Compute whether the current process already uses the requested conda environment
conda_environment_is_active() {
  local requested_env="$1"
  local default_env="${CONDA_DEFAULT_ENV:-}"
  local conda_prefix="${CONDA_PREFIX:-}"
  local python_prefix=""

  requested_env="${requested_env%/}"
  default_env="${default_env%/}"
  conda_prefix="${conda_prefix%/}"
  if [[ "${requested_env}" == */* ]]; then
    [[ "${conda_prefix}" == "${requested_env}" && "$(command -v python)" == "${requested_env}/bin/python" ]]
    return
  fi
  if [[ "${default_env}" == "${requested_env}" || "${conda_prefix}" == "${requested_env}" ||
        "${default_env##*/}" == "${requested_env##*/}" ||
        "${conda_prefix##*/}" == "${requested_env##*/}" ]]; then
    return 0
  fi

  if command -v python >/dev/null 2>&1; then
    python_prefix="$(python -c 'import sys; print(sys.prefix)' 2>/dev/null || true)"
  fi
  python_prefix="${python_prefix%/}"
  [[ "${python_prefix}" == "${requested_env}" ||
     "${python_prefix##*/}" == "${requested_env##*/}" ]]
}


# Compute whether the current process already uses the graniitti conda environment
graniitti_environment_is_active() {
  conda_environment_is_active graniitti
}


# Activate a conda environment only when the current interpreter is outside it
activate_conda_environment() {
  local requested_env="$1"
  if conda_environment_is_active "${requested_env}"; then
    return
  fi

  if ! command -v conda >/dev/null 2>&1; then
    echo "Cannot activate '${requested_env}' because conda is unavailable" >&2
    return 64
  fi
  if [[ "$(type -t conda)" != "function" ]]; then
    eval "$(conda shell.bash hook)"
  fi
  conda activate "${requested_env}"
}


# Run a command directly in an active environment or resolve graniitti through conda
run_in_graniitti() {
  if graniitti_environment_is_active; then
    "$@"
    return
  fi

  if ! command -v conda >/dev/null 2>&1; then
    echo "The graniitti environment is inactive and conda is unavailable" >&2
    return 64
  fi
  conda run --no-capture-output -n graniitti "$@"
}

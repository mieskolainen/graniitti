#!/usr/bin/env bash
#
# Shared generation and analysis controls for manual studies
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>

# Resolve and validate the requested study action
study_resolve_action() {
  if [[ -z "${STUDY_ACTION:-}" ]]; then
    STUDY_ACTION="both"
  fi

  case "${STUDY_ACTION}" in
    both|analyze)
      ;;
    *)
      echo "STUDY_ACTION must be 'both' or 'analyze'" >&2
      return 64
      ;;
  esac

  export STUDY_ACTION
}

# Initialize one study from any working directory
study_initialize() {
  local activation_status=0
  local common_dir
  local default_events="${1:-10000}"
  local restore_nounset=0

  common_dir="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
  STUDY_REPO_ROOT="$(CDPATH= cd -- "${common_dir}/../../.." && pwd)"

  cd -- "${STUDY_REPO_ROOT}"
  source "${STUDY_REPO_ROOT}/tests/environment.sh"

  study_resolve_action

  NEVENTS="${NEVENTS:-${default_events}}"
  if [[ ! "${NEVENTS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "NEVENTS must be a positive integer" >&2
    return 64
  fi
  WEIGHTED="${WEIGHTED:-true}"
  LOOPSCREEN="${LOOPSCREEN:-false}"
  if [[ ! "${WEIGHTED}" =~ ^(0|1|true|false)$ ]]; then
    echo "WEIGHTED must be 0, 1, true or false" >&2
    return 64
  fi
  if [[ ! "${LOOPSCREEN}" =~ ^(0|1|true|false)$ ]]; then
    echo "LOOPSCREEN must be 0, 1, true or false" >&2
    return 64
  fi

  export NEVENTS WEIGHTED LOOPSCREEN STUDY_REPO_ROOT

  if [[ "$-" == *u* ]]; then
    restore_nounset=1
    set +u
  fi

  activate_conda_environment graniitti || activation_status=$?

  if [[ "${activation_status}" == "0" && -z "${GRANIITTI_ENV:-}" ]]; then
    source "${STUDY_REPO_ROOT}/install/setenv.sh"
  fi

  if [[ "${restore_nounset}" == "1" ]]; then
    set -u
  fi

  return "${activation_status}"
}


# Run gr with the common event generation controls
study_gr() {
  local argument
  local -a arguments=()

  while [[ "$#" -gt 0 ]]; do
    argument="$1"
    case "${argument}" in
      -n|-w|-l|--NEVENTS|--WEIGHTED|--LOOPSCREEN)
        if [[ "$#" -lt 2 ]]; then
          echo "${argument} requires a value" >&2
          return 64
        fi
        shift 2
        ;;
      --NEVENTS=*|--WEIGHTED=*|--LOOPSCREEN=*)
        shift
        ;;
      *)
        arguments+=("${argument}")
        shift
        ;;
    esac
  done

  "${STUDY_REPO_ROOT}/bin/gr" "${arguments[@]}" \
    --NEVENTS "${NEVENTS}" --WEIGHTED "${WEIGHTED}" --LOOPSCREEN "${LOOPSCREEN}"
}


# Run the histogram analyzer with the common event record limit
study_analyze() {
  local argument
  local -a arguments=()

  while [[ "$#" -gt 0 ]]; do
    argument="$1"
    case "${argument}" in
      -t)
        local screening="bare"
        [[ "${LOOPSCREEN}" =~ ^(1|true)$ ]] && screening="screened"
        arguments+=("-t" "$2 [${SCREENING_LABEL:-${screening}}]")
        shift 2
        ;;
      -X|--maximum)
        if [[ "$#" -lt 2 ]]; then
          echo "${argument} requires a value" >&2
          return 64
        fi
        shift 2
        ;;
      --maximum=*)
        shift
        ;;
      *)
        arguments+=("${argument}")
        shift
        ;;
    esac
  done

  "${STUDY_REPO_ROOT}/bin/analyze" "${arguments[@]}" --maximum "${NEVENTS}"
}


# Run the harmonic analyzer with the common event record limit
study_fitharmonic() {
  local argument
  local -a arguments=()

  while [[ "$#" -gt 0 ]]; do
    argument="$1"
    case "${argument}" in
      -X|--maximum)
        if [[ "$#" -lt 2 ]]; then
          echo "${argument} requires a value" >&2
          return 64
        fi
        shift 2
        ;;
      --maximum=*)
        shift
        ;;
      *)
        arguments+=("${argument}")
        shift
        ;;
    esac
  done

  "${STUDY_REPO_ROOT}/bin/fitharmonic" "${arguments[@]}" \
    --maximum "${NEVENTS}"
}


# Compute whether event generation is enabled
study_generates() {
  [[ "${STUDY_ACTION}" == "both" ]]
}


# End one unavailable optional study without failing an automatic run
study_skip() {
  echo "[study: skip] $*"
  exit 0
}


# Require an explicit environment opt in for an external workflow
study_require_opt_in() {
  local variable="$1"
  local description="$2"

  if [[ "${!variable:-0}" != "1" ]]; then
    study_skip "${description} requires ${variable}=1"
  fi
}


# Skip an optional study when one required command is unavailable
study_require_command() {
  local command="$1"
  local description="$2"

  if ! command -v "${command}" >/dev/null 2>&1; then
    study_skip "${description} command is unavailable: ${command}"
  fi
}


# Skip an optional study when any required input file is unavailable
study_require_files() {
  local description="$1"
  shift
  local path

  for path in "$@"; do
    if [[ ! -f "${path}" ]]; then
      study_skip "${description} input is unavailable: ${path}"
    fi
  done
}


# Require generated outputs when running an analysis only study
study_require_outputs() {
  local description="$1"
  shift
  local path

  for path in "$@"; do
    if [[ ! -f "${path}" ]]; then
      echo "${description} output is unavailable: ${path}" >&2
      return 66
    fi
  done
}

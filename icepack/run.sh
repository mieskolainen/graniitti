#!/usr/bin/env bash
#
# Generate and compare one icepack entry
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

set -euo pipefail

# Print the command-line interface
usage() {
  printf '%s\n' \
    "Usage: bash icepack/run.sh [options] <icepack-entry>..." \
    "Select a folder recursively or a dataset.json file for exactly one entry" \
    "" \
    "Options:" \
    "  --NEVENTS N                 Override dataset or card event count" \
    "  --LOOPSCREEN BOOL           Override soft survival screening (ELASTIC always uses 1)" \
    "  --WEIGHTED BOOL             Weighted event generation" \
    "  --DENSITY BOOL              Force iceplot density normalization" \
    "  --INTEGRATOR NAME           VEGAS, FLAT or NEUROJAC" \
    "  --output-dir PATH           Output root containing one folder per icepack" \
    "  --report-dir PATH           JSON report directory" \
    "  --skip-generate             Use existing generated samples" \
    "  --all                       Run every active icepack entry" \
    "  --list                      List selected entries without running them" \
    "  -h, --help                  Show this help"
}

# Resolve one input relative to the icepack root
resolve_entry_input() {
  local input="${1%/}"
  local icepack_root="$2"

  input="${input#./}"
  if [[ "${input}" == "*" || "${input}" == "icepack" || "${input}" == "icepack/*" ]]; then
    printf '%s\n' "${icepack_root}"
    return
  fi
  if [[ "${input}" == "${icepack_root}" || "${input}" == "${icepack_root}/"* ]]; then
    realpath -m "${input}"
    return
  fi

  input="${input#icepack/}"
  realpath -m "${icepack_root}/${input}"
}

# Compute whether an icepack path contains a helper component
is_helper_path() {
  local relative="$1"
  local component
  local -a components

  IFS="/" read -r -a components <<< "${relative}"
  for component in "${components[@]}"; do
    if [[ "${component}" == _* || "${component}" == *._old* ]]; then
      return 0
    fi
  done
  return 1
}

# Save a replay command with the project environment and original argument boundaries
save_command() {
  local filename="$1"
  shift
  {
    printf '#!/usr/bin/env bash\nset -eo pipefail\n'
    printf 'cd %q\n' "${REPO_ROOT}"
    printf '%s\n' 'source tests/environment.sh' 'activate_conda_environment graniitti' 'source install/setenv.sh'
    printf '%q ' "$@"
    printf '\n'
  } > "${filename}"
}

# Keep terminal progress active while flushing output to the console and log
run_logged() {
  local filename="$1"
  local command
  shift
  printf -v command '%q ' "$@"
  SHELL="${BASH}" run_in_graniitti script \
    --quiet --return --flush --command "${command}" "${filename}" </dev/null
}

# Reserve a UTC timestamp shared by the plots, run log, replay script and generation records
new_run_log() {
  local directory="$1" stamp name index=1
  stamp="$(date -u +%Y-%m-%d_%H-%M-%S_UTC)"
  name="${stamp}"
  while ! (set -o noclobber; : > "${directory}/${name}.log") 2>/dev/null; do
    [[ -e "${directory}/${name}.log" ]] || return 1
    index=$((index + 1))
    name="${stamp}_${index}"
  done
  printf '%s\n' "${directory}/${name}.log"
}

# Process and compare one complete icepack entry
run_dataset() {
  local dataset_dir="$1"
  local icepack_root="$2"
  local nevents="$3"
  local loopscreen="$4"
  local weighted="$5"
  local integrator="$6"
  local density_force="$7"
  local skip_generate="$8"
  local output_dir="$9"
  local report_dir="${10}"
  local run_log="${11}"
  local plot_name="plots/$(basename "${run_log}" .log)"
  local plot_dir="${output_dir}/${plot_name}"
  local dataset="${dataset_dir}/dataset.json"
  local entry="${dataset_dir#"${icepack_root}/"}"
  local sample_tag="${entry//\//__}"
  local report="${report_dir}/${sample_tag}.json"
  local generation_dir sample_dir sample_count stack density output label gencard scale assignment
  local fragment_engine fragment_card fragment_mode fragment_seed fragment_attempts override_count
  local cursor=0 sample_index override_index
  local -a plan_fields=() sample_outputs=() sample_labels=() sample_scales=() sample_assignments=()
  local -a gr_args=() fragment_args=() ice_args=()

  mapfile -d '' -t plan_fields < <(
    run_in_graniitti python -m core.io.steering generation-plan \
      --dataset "${dataset}" \
      --cdir "${REPO_ROOT}" \
      --output-prefix "${sample_tag}"
  )
  if [[ "${#plan_fields[@]}" -lt 5 ]]; then
    echo "Could not resolve generator samples from ${dataset}" >&2
    return 65
  fi

  sample_count="${plan_fields[cursor++]}"
  stack="${plan_fields[cursor++]}"
  density="${plan_fields[cursor++]}"
  # Always include eikonal screening for elastic scattering
  if [[ "${entry}" == ELASTIC/* ]]; then
    loopscreen=1
  elif [[ -z "${loopscreen}" ]]; then
    loopscreen="${plan_fields[cursor]}"
  fi
  cursor=$((cursor + 1))
  if [[ -z "${nevents}" ]]; then
    nevents="${plan_fields[cursor]}"
  fi
  cursor=$((cursor + 1))
  if [[ "${skip_generate}" != "1" ]]; then
    mkdir -p "${output_dir}/generation"
    generation_dir="${output_dir}/generation/$(basename "${run_log}" .log)"
    mkdir "${generation_dir}"
    cp -- "${dataset}" "${generation_dir}/dataset.json"
    echo "Generation records: ${generation_dir}"
  fi
  for ((sample_index = 0; sample_index < sample_count; ++sample_index)); do
    output="${plan_fields[cursor++]}"
    label="${plan_fields[cursor++]}"
    gencard="${plan_fields[cursor++]}"
    scale="${plan_fields[cursor++]}"
    assignment="${plan_fields[cursor++]}"
    fragment_engine="${plan_fields[cursor++]}"
    fragment_card="${plan_fields[cursor++]}"
    fragment_mode="${plan_fields[cursor++]}"
    fragment_seed="${plan_fields[cursor++]}"
    fragment_attempts="${plan_fields[cursor++]}"
    override_count="${plan_fields[cursor++]}"

    local output_format="hepmc3"
    if [[ "${fragment_engine}" == "pythia" ]]; then
      output_format="lhe"
    fi
    gr_args=(-i "${gencard}" -f "${output_format}" -o "${output}")
    if [[ -n "${nevents}" ]]; then
      gr_args+=(--NEVENTS "${nevents}")
    fi
    if [[ -n "${loopscreen}" ]]; then
      gr_args+=(--LOOPSCREEN "${loopscreen}")
    fi
    if [[ -n "${weighted}" ]]; then
      gr_args+=(--WEIGHTED "${weighted}")
    fi
    if [[ -n "${integrator}" ]]; then
      gr_args+=(--INTEGRATOR "${integrator}")
    fi
    for ((override_index = 0; override_index < override_count; ++override_index)); do
      gr_args+=(--set "${plan_fields[cursor++]}")
    done

    if [[ "${skip_generate}" != "1" ]]; then
      if [[ "${fragment_engine}" == "pythia" ]]; then
        if [[ -z "${nevents}" ]]; then
          echo "Pythia fragmentation requires validation.nevents or --NEVENTS in ${dataset}" >&2
          return 65
        fi
        if [[ ! -x ./bin/pythia_lhe_hadronize ]]; then
          bash tests/external/pythia/drivers/lhe_converter/build.sh
        fi
      fi
      sample_dir="${generation_dir}/${output}"
      mkdir -p "${sample_dir}"
      cp -- "${gencard}" "${sample_dir}/gencard.json"
      gr_args[1]="${sample_dir}/gencard.json"
      save_command "${sample_dir}/command.sh" ./bin/gr "${gr_args[@]}"
      run_logged "${sample_dir}/generator.log" ./bin/gr "${gr_args[@]}"
      if [[ "${fragment_engine}" == "pythia" ]]; then
        cp -- "${fragment_card}" "${sample_dir}/fragmentation.cmnd"
        fragment_args=(
          "output/${output}.lhe" "output/${output}.hepmc3" "${nevents}"
          "${fragment_seed}" "${sample_dir}/fragmentation.cmnd" "${fragment_attempts}"
          "${fragment_mode}"
        )
        save_command "${sample_dir}/fragmentation.sh" ./bin/pythia_lhe_hadronize "${fragment_args[@]}"
        run_logged "${sample_dir}/fragmentation.log" ./bin/pythia_lhe_hadronize "${fragment_args[@]}"
      fi
    fi
    sample_outputs+=("output/${output}.hepmc3")
    sample_labels+=("${label}")
    sample_scales+=("${scale}")
    sample_assignments+=("${assignment}")
  done

  if [[ "${cursor}" -ne "${#plan_fields[@]}" ]]; then
    echo "Invalid generator plan field count for ${dataset}" >&2
    return 65
  fi

  ice_args=(
    --hepmc3 "${sample_outputs[@]}"
    --mclabel "${sample_labels[@]}"
    --mcscale "${sample_scales[@]}"
    --analysis "${dataset_dir}"
    --output "${plot_name}"
    --output-dir "${output_dir}"
    --report "${report}"
    --validate
  )
  if [[ "${stack}" == "1" ]]; then
    ice_args+=(--stack)
  fi
  case "${density_force:-${density}}" in
    1 | true) ice_args+=(--density) ;;
    0 | false) ice_args+=(--no-density) ;;
  esac
  if [[ -n "${sample_assignments[0]}" ]]; then
    ice_args+=(--mc-hepdata-sample "${sample_assignments[@]}")
  fi
  save_command "${run_log%.log}.sh" python -m core.iceplot "${ice_args[@]}"
  run_in_graniitti python -m core.iceplot "${ice_args[@]}"
  printf 'Plots: %s\n' "${plot_dir}"
}

NEVENTS_VALUE="${NEVENTS:-}"
LOOPSCREEN_VALUE="${LOOPSCREEN:-}"
WEIGHTED_VALUE="${WEIGHTED:-1}"
DENSITY_VALUE="${DENSITY:-}"
INTEGRATOR_VALUE="${INTEGRATOR:-}"
OUTPUT_DIR_VALUE="${ICEPACK_OUTPUT_DIR:-figs/icepack}"
REPORT_DIR_VALUE="${ICEPACK_REPORT_DIR:-}"
SKIP_GENERATE=0
LIST_DATASETS=0
ENTRY_INPUTS=()

while [[ "$#" -gt 0 ]]; do
  option="${1%%=*}"
  case "${option}" in
    --NEVENTS | --LOOPSCREEN | --WEIGHTED | --DENSITY | --INTEGRATOR | --output-dir | --report-dir)
      if [[ "$1" == *=* ]]; then
        value="${1#*=}"
      else
        if [[ "$#" -lt 2 ]]; then
          echo "Missing value for $1" >&2
          exit 64
        fi
        value="$2"
        shift
      fi
      case "${option}" in
        --NEVENTS) NEVENTS_VALUE="${value}" ;;
        --LOOPSCREEN) LOOPSCREEN_VALUE="${value}" ;;
        --WEIGHTED) WEIGHTED_VALUE="${value}" ;;
        --DENSITY) DENSITY_VALUE="${value}" ;;
        --INTEGRATOR) INTEGRATOR_VALUE="${value}" ;;
        --output-dir) OUTPUT_DIR_VALUE="${value}" ;;
        --report-dir) REPORT_DIR_VALUE="${value}" ;;
      esac
      ;;
    *)
      case "$1" in
        --skip-generate) SKIP_GENERATE=1 ;;
        --all) ENTRY_INPUTS+=("*") ;;
        --list) LIST_DATASETS=1 ;;
        -h | --help) usage; exit 0 ;;
        -*) echo "Unknown option: $1" >&2; usage >&2; exit 64 ;;
        *) ENTRY_INPUTS+=("$1") ;;
      esac
      ;;
  esac
  shift
done

if [[ "${#ENTRY_INPUTS[@]}" -eq 0 ]]; then
  usage >&2
  exit 64
fi

if [[ -n "${NEVENTS_VALUE}" && ! "${NEVENTS_VALUE}" =~ ^[1-9][0-9]*$ ]]; then
  echo "NEVENTS must be a positive integer: ${NEVENTS_VALUE}" >&2
  exit 64
fi
for option in LOOPSCREEN WEIGHTED DENSITY; do
  variable="${option}_VALUE"
  if [[ -n "${!variable}" && ! "${!variable}" =~ ^(0|1|true|false)$ ]]; then
    echo "${option} must be 0, 1, true or false: ${!variable}" >&2
    exit 64
  fi
done
if [[ -n "${INTEGRATOR_VALUE}" && ! "${INTEGRATOR_VALUE}" =~ ^(VEGAS|FLAT|NEUROJAC)$ ]]; then
  echo "INTEGRATOR must be VEGAS, FLAT or NEUROJAC: ${INTEGRATOR_VALUE}" >&2
  exit 64
fi

SCRIPT_DIR="$(realpath -e "$(dirname "${BASH_SOURCE[0]}")")"
REPO_ROOT="$(realpath -e "${SCRIPT_DIR}/..")"
ICEPACK_ROOT="$(realpath -e "${REPO_ROOT}/icepack")"
OUTPUT_DIR_VALUE="$(realpath -m "${OUTPUT_DIR_VALUE}")"
if [[ -n "${REPORT_DIR_VALUE}" ]]; then
  REPORT_DIR_VALUE="$(realpath -m "${REPORT_DIR_VALUE}")"
fi

shopt -s globstar nullglob
declare -A SEEN_DATASET_DIRS=()
DATASET_DIRS=()
for entry_input in "${ENTRY_INPUTS[@]}"; do
  candidate_dir="$(resolve_entry_input "${entry_input}" "${ICEPACK_ROOT}")"
  case "${candidate_dir}" in
    "${ICEPACK_ROOT}" | "${ICEPACK_ROOT}"/*) ;;
    *) continue ;;
  esac
  if [[ -f "${candidate_dir}" && "${candidate_dir}" == */dataset.json ]]; then
    dataset_cards=("${candidate_dir}")
  elif [[ -d "${candidate_dir}" ]]; then
    dataset_cards=("${candidate_dir}"/**/dataset.json)
  else
    continue
  fi

  for dataset_card in "${dataset_cards[@]}"; do
    dataset_dir="${dataset_card%/dataset.json}"
    relative_dir="${dataset_dir#"${ICEPACK_ROOT}/"}"
    if is_helper_path "${relative_dir}"; then
      continue
    fi
    if [[ -z "${SEEN_DATASET_DIRS["${dataset_dir}"]+set}" ]]; then
      SEEN_DATASET_DIRS["${dataset_dir}"]=1
      DATASET_DIRS+=("${dataset_dir}")
    fi
  done
done
shopt -u globstar nullglob

if [[ "${#DATASET_DIRS[@]}" -eq 0 ]]; then
  echo "No complete icepack entries found" >&2
  exit 66
fi
cd "${REPO_ROOT}"
source tests/environment.sh
set +u
activate_conda_environment graniitti
set -u
if [[ -z "${GRANIITTI_ENV:-}" ]]; then
  set +u
  source install/setenv.sh >/dev/null
  set -u
fi

active_datasets="$(run_in_graniitti python -m core.io.steering list-datasets \
  --cdir "${REPO_ROOT}" --dataset "${DATASET_DIRS[@]/%//dataset.json}")"
if [[ -z "${active_datasets}" ]]; then
  echo 'No active icepack entries selected' >&2
  exit 66
fi
mapfile -t DATASET_DIRS <<< "${active_datasets}"
if [[ "${LIST_DATASETS}" == "1" ]]; then
  printf '%s\n' "${DATASET_DIRS[@]#"${ICEPACK_ROOT}/"}"
  exit 0
fi

run_index=0
failed_entries=()
for dataset_dir in "${DATASET_DIRS[@]}"; do
  run_index=$((run_index + 1))
  entry="${dataset_dir#"${ICEPACK_ROOT}/"}"
  printf '\nicepack/run.sh: [%d/%d] %s\n' \
    "${run_index}" "${#DATASET_DIRS[@]}" "${entry}"
  output_dir="${OUTPUT_DIR_VALUE}/${entry//\//__}"
  report_dir="${REPORT_DIR_VALUE:-${output_dir}/reports}"
  log_dir="${output_dir}/logs"
  mkdir -p "${log_dir}" "${report_dir}"
  run_log="$(new_run_log "${log_dir}")"
  printf 'Run log: %s\nReport: %s\n' "${run_log}" "${report_dir}/${entry//\//__}.json"
  # Keep errexit inside each dataset while allowing subsequent datasets to run
  set +e
  (
    set -e
    run_dataset \
    "${dataset_dir}" \
    "${ICEPACK_ROOT}" \
    "${NEVENTS_VALUE}" \
    "${LOOPSCREEN_VALUE}" \
    "${WEIGHTED_VALUE}" \
    "${INTEGRATOR_VALUE}" \
    "${DENSITY_VALUE}" \
    "${SKIP_GENERATE}" \
    "${output_dir}" \
    "${report_dir}" \
    "${run_log}"
  ) 2>&1 | tee "${run_log}"
  run_status=$?
  set -e
  if [[ "${run_status}" -ne 0 ]]; then
    failed_entries+=("${entry}")
    printf 'FAILED: %s (exit %d), log: %s\n' "${entry}" "${run_status}" "${run_log}" >&2
  fi
done
if [[ "${#failed_entries[@]}" -gt 0 ]]; then
  printf '\nFailed icepacks (%d/%d):\n' "${#failed_entries[@]}" "${#DATASET_DIRS[@]}" >&2
  printf '  %s\n' "${failed_entries[@]}" >&2
  exit 1
fi

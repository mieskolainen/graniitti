#!/usr/bin/env bash
set -euo pipefail

# Generate or consume a GRANIITTI LHE file and pass it through Pythia8
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
cd "$REPO_ROOT"
if [[ "${1:-}" == --help || "${1:-}" == -h ]]; then
  echo "Usage: bash tests/external/pythia/drivers/lhe_converter/run.sh CARD_OR_LHE TAG [GRANIITTI_OPTIONS]"
  exit 0
fi
if [[ $# -lt 2 ]]; then
  echo "Usage: bash tests/external/pythia/drivers/lhe_converter/run.sh CARD_OR_LHE TAG [GRANIITTI_OPTIONS]" >&2
  exit 2
fi

INPUT="$1"
TAG="$2"
shift 2
if [[ "$INPUT" == *.lhe ]]; then
  NEVENTS="${NEVENTS:-all}"
else
  NEVENTS="${NEVENTS:-20}"
fi
SEED="${SEED:-12345}"
LOOPSCREEN="${LOOPSCREEN:-}"
WEIGHTED="${WEIGHTED:-}"
PYTHIA_CMND="${PYTHIA_CMND-tests/external/pythia/drivers/lhe_converter/shower.cmnd}"
RESHOWER_ATTEMPTS="${RESHOWER_ATTEMPTS:-500}"
CONVERTER="${CONVERTER:-bin/pythia_lhe_hadronize}"

# Require an optional GRANIITTI binary switch to be 0 or 1
validate_binary_switch() {
  local name="$1"
  local value="$2"

  if [[ -n "$value" && ! "$value" =~ ^[01]$ ]]; then
    echo "${name} must be 0 or 1: ${value}" >&2
    exit 2
  fi
}

validate_binary_switch "LOOPSCREEN" "$LOOPSCREEN"
validate_binary_switch "WEIGHTED" "$WEIGHTED"

if [[ ! "$TAG" =~ ^[a-zA-Z0-9][a-zA-Z0-9_.-]*$ ]]; then
  echo "TAG must be a filename stem without directories" >&2
  exit 2
fi
if [[ ! "$NEVENTS" =~ ^[1-9][0-9]*$ && !( "$INPUT" == *.lhe && "$NEVENTS" == all ) ]] ||
   [[ ! "$SEED" =~ ^[0-9]+$ || ${#SEED} -gt 9 ]]; then
  echo "NEVENTS must be positive (or all for LHE input), SEED must be between 0 and 900000000" >&2
  exit 2
fi
if (( 10#$SEED > 900000000 )); then
  echo "SEED must be between 0 and 900000000" >&2
  exit 2
fi
if [[ ! "$RESHOWER_ATTEMPTS" =~ ^[1-9][0-9]*$ || ! "${CONVERTER_MODE:-auto}" =~ ^(auto|fragment|shower)$ ]]; then
  echo "RESHOWER_ATTEMPTS must be positive, CONVERTER_MODE must be auto, fragment or shower" >&2
  exit 2
fi
if [[ ! -s "$INPUT" || ( -n "$PYTHIA_CMND" && ! -s "$PYTHIA_CMND" ) ]]; then
  echo "Missing or empty input card, LHE file or PYTHIA_CMND" >&2
  exit 2
fi

if [[ "$INPUT" == *.lhe && ( $# -ne 0 || -n "$LOOPSCREEN" || -n "$WEIGHTED" ) ]]; then
  echo "GRANIITTI options, LOOPSCREEN and WEIGHTED require a generator card" >&2
  exit 2
fi
set +u
source tests/environment.sh
activate_conda_environment graniitti
source install/setenv.sh
set -u

if [[ ! -x "$CONVERTER" ]]; then
  OUTPUT="$CONVERTER" bash tests/external/pythia/drivers/lhe_converter/build.sh
fi
mkdir -p output
if [[ "$INPUT" == *.lhe ]]; then
  LHE_FILE="$INPUT"
else
  GENERATOR_OPTIONS=(--set 'SCATTERING.BEAMFRAG="diquark"' "$@")
  if [[ -n "$LOOPSCREEN" ]]; then
    GENERATOR_OPTIONS+=(--LOOPSCREEN "$LOOPSCREEN")
  fi
  if [[ -n "$WEIGHTED" ]]; then
    GENERATOR_OPTIONS+=(--WEIGHTED "$WEIGHTED")
  fi
  LHE_FILE="output/${TAG}.lhe"
  run_in_graniitti \
    ./bin/gr -i "$INPUT" -f lhe -o "$TAG" -n "$NEVENTS" -r "$SEED" "${GENERATOR_OPTIONS[@]}"
fi

HEPMC_OUT="output/${TAG}.hepmc3"
run_in_graniitti \
  "$CONVERTER" "$LHE_FILE" "$HEPMC_OUT" "$NEVENTS" "$SEED" \
    "$PYTHIA_CMND" "$RESHOWER_ATTEMPTS" "${CONVERTER_MODE:-auto}"

echo "LHE input: ${LHE_FILE}"
echo "Pythia output: ${HEPMC_OUT}"

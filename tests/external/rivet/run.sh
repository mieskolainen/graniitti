#!/usr/bin/env bash
set -euo pipefail
# Generate events and run one Rivet analysis with separate outputs for each case
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO_ROOT"
CASE="${1:-}"
ACTION="${2:-all}"
if [[ "$CASE" == -h || "$CASE" == --help || "$ACTION" == -h || "$ACTION" == --help ]]; then
  echo "Usage: bash tests/external/rivet/run.sh cms|generic|star [all|generate|analyze]"
  exit 0
fi
if [[ $# -gt 2 || ! "$ACTION" =~ ^(all|generate|analyze)$ ]]; then
  echo "Expected all, generate or analyze" >&2
  exit 2
fi
case "$CASE" in
  cms|generic)
    CARD=icepack/GAMMA/integrated/CMS_2012/gencard.json
    ANALYSIS=CMS_2011_I954992
    COUNT=100000
    if [[ "$CASE" == generic ]]; then
      ANALYSIS=MC_FSPARTICLES
      COUNT=1000
    fi
    ;;
  star)
    CARD=icepack/SOFTCEP/STAR_1792394/pipi/gencard.json
    ANALYSIS=STAR_2020_I1792394
    COUNT=50000
    ;;
  *) echo "Expected cms, generic or star" >&2; exit 2 ;;
esac
TAG="${TAG:-rivet_$CASE}"
if [[ ! "$TAG" =~ ^[a-zA-Z0-9][a-zA-Z0-9_.-]*$ ]]; then
  echo "TAG must be a filename stem without directories" >&2
  exit 2
fi
INPUT="output/$TAG.hepmc3"
RUN_DIR="runs/tests/rivet/$TAG"

# Activate Rivet only in a subshell so GRANIITTI keeps its own libraries
setup_rivet() {
  local setup="${RIVET_ENV:-$REPO_ROOT/.local/rivet/rivetenv.sh}"
  if [[ -f "$setup" ]]; then
    set +u
    source "$setup"
    set -u
  elif [[ -n "${RIVET_ENV:-}" ]]; then
    echo "Cannot read RIVET_ENV: $setup" >&2
    return 1
  fi
  command -v rivet >/dev/null && command -v rivet-mkhtml >/dev/null || {
    echo "Rivet is unavailable, activate it or set RIVET_ENV to its setup script" >&2
    return 1
  }
}

if [[ "$ACTION" != generate ]]; then
  (setup_rivet; rivet --show-analysis "$ANALYSIS" >/dev/null)
fi
if [[ "$ACTION" != analyze ]]; then
  NEVENTS="${NEVENTS:-$COUNT}"
  if [[ ! "$NEVENTS" =~ ^[1-9][0-9]*$ && ( "$ACTION" != generate || "$NEVENTS" != -1 ) ]]; then
    echo "NEVENTS must be positive, or -1 with generate to save a proposal" >&2
    exit 2
  fi
  OPTIONS=()
  [[ -z "${CORES:-}" ]] || OPTIONS+=(--CORES "$CORES")
  [[ -z "${VGRID:-}" ]] || OPTIONS+=(-d "$VGRID")
  for name in LOOPSCREEN WEIGHTED; do
    value="${!name:-}"
    if [[ -n "$value" ]]; then
      [[ "$value" =~ ^[01]$ ]] || { echo "$name must be 0 or 1" >&2; exit 2; }
      OPTIONS+=("--$name" "$value")
    fi
  done
  (
    set +u
    source tests/environment.sh
    activate_conda_environment graniitti
    source install/setenv.sh
    [[ -x bin/gr ]] || { echo "Build GRANIITTI first with cmake --build build -j4" >&2; exit 1; }
    ./bin/gr -i "$CARD" -n "$NEVENTS" -r "${SEED:-12345}" -o "$TAG" -f hepmc3 -h 0 "${OPTIONS[@]}"
  )
fi
if [[ "$ACTION" != generate ]]; then
  [[ -s "$INPUT" ]] || { echo "Missing or empty input: $INPUT" >&2; exit 1; }
  mkdir -p "$RUN_DIR"
  (
    setup_rivet
    rivet "$INPUT" --analysis="$ANALYSIS" -o "$RUN_DIR/Rivet.yoda"
    rivet-mkhtml "$RUN_DIR/Rivet.yoda:GRANIITTI" -o "$RUN_DIR/plots"
  )
fi

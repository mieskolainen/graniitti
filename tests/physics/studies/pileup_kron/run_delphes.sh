#!/usr/bin/env bash
set -euo pipefail

# Run a Delphes PU study and report Delphes jet baselines
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
cd "$REPO_ROOT"
set +u
source tests/environment.sh
activate_conda_environment graniitti
source install/setenv.sh
set -u

NEVENTS="${NEVENTS:-300}"
PU="${PU:-50}"
SEED_HARD="${SEED_HARD:-31001}"
SEED_PU="${SEED_PU:-91001}"
REGENERATE="${REGENERATE:-0}"
OUTPUT_DIR="${OUTPUT_DIR:-output/graph_pu${PU}}"
HARD_CMND="${HARD_CMND:-tests/physics/studies/pileup_kron/pythia_hard_qcd.cmnd}"
PU_CMND="${PU_CMND:-tests/physics/studies/pileup_kron/pythia_minbias.cmnd}"
DELPHES_DIR="${DELPHES_DIR:?Set DELPHES_DIR to your Delphes installation directory}"
DELPHES_SOURCE_DIR="${DELPHES_SOURCE_DIR:-${DELPHES_DIR}}"
DELPHES_BIN_DIR="${DELPHES_BIN_DIR:-${DELPHES_DIR}}"
DELPHES_HEPMC3="${DELPHES_HEPMC3:-}"
DELPHES_HEPMC2PILEUP="${DELPHES_HEPMC2PILEUP:-}"
DELPHES_CARD_TEMPLATE="${DELPHES_CARD_TEMPLATE:-${DELPHES_SOURCE_DIR}/cards/CMS_PhaseII/CMS_PhaseII_200PU_v03_nodtf.tcl}"
DELPHES_CARD="${DELPHES_CARD:-${OUTPUT_DIR}/delphes_card_CMS_PhaseII_PU${PU}_graph.tcl}"
DELPHES_ROOT="${DELPHES_ROOT:-${OUTPUT_DIR}/delphes_pu${PU}_${NEVENTS}.root}"
MATCH_RADIUS="${MATCH_RADIUS:-0.30}"
REFERENCE="${REFERENCE:-GenHardJet}"
GEN_PT_MIN="${GEN_PT_MIN:-30.0}"
RECO_PT_MIN="${RECO_PT_MIN:-15.0}"
ETA_MAX="${ETA_MAX:-4.7}"
MAX_EVENTS="${MAX_EVENTS:-0}"
EVENT_START="${EVENT_START:-0}"
EVENT_INDICES="${EVENT_INDICES:-}"
EVENT_INDICES_FILE="${EVENT_INDICES_FILE:-}"
CHECK_DELPHES_RUNTIME="${CHECK_DELPHES_RUNTIME:-1}"
INCLUDE_GRAPH_KRON="${INCLUDE_GRAPH_KRON:-1}"
GRAPH_KRON_COLLECTION="${GRAPH_KRON_COLLECTION:-JetGraphKron}"
GRAPH_KRON_BASE_AXIS="${GRAPH_KRON_BASE_AXIS:-raw}"
GRAPH_KRON_CANDIDATE_PT_MIN="${GRAPH_KRON_CANDIDATE_PT_MIN:-0.50}"
GRAPH_KRON_CLUSTER_MODE="${GRAPH_KRON_CLUSTER_MODE:-posterior}"
GRAPH_KRON_MAX_CANDIDATES="${GRAPH_KRON_MAX_CANDIDATES:-1200}"
GRAPH_KRON_GRAPH_RADIUS="${GRAPH_KRON_GRAPH_RADIUS:-0.45}"
GRAPH_KRON_GRAPH_SIGMA="${GRAPH_KRON_GRAPH_SIGMA:-0.25}"
GRAPH_KRON_RADIATION_RADIUS="${GRAPH_KRON_RADIATION_RADIUS:-1.20}"
GRAPH_KRON_RADIATION_CORE="${GRAPH_KRON_RADIATION_CORE:-0.04}"
GRAPH_KRON_RADIATION_ANGULAR_POWER="${GRAPH_KRON_RADIATION_ANGULAR_POWER:-1.0}"
GRAPH_KRON_RADIATION_PT_POWER="${GRAPH_KRON_RADIATION_PT_POWER:-1.0}"
GRAPH_KRON_RADIATION_EDGE_STRENGTH="${GRAPH_KRON_RADIATION_EDGE_STRENGTH:-0.08}"
GRAPH_KRON_RADIATION_EDGE_MAX_NEIGHBOURS="${GRAPH_KRON_RADIATION_EDGE_MAX_NEIGHBOURS:-16}"
GRAPH_KRON_RADIATION_UNARY_STRENGTH="${GRAPH_KRON_RADIATION_UNARY_STRENGTH:-2.0}"
GRAPH_KRON_ANCHOR_STRENGTH="${GRAPH_KRON_ANCHOR_STRENGTH:-25.0}"
GRAPH_KRON_LOCAL_SUPPORT_STRENGTH="${GRAPH_KRON_LOCAL_SUPPORT_STRENGTH:-1.5}"
GRAPH_KRON_FISHER_MODE="${GRAPH_KRON_FISHER_MODE:-sample}"
GRAPH_KRON_FISHER_STRENGTH="${GRAPH_KRON_FISHER_STRENGTH:-2.0}"
GRAPH_KRON_FISHER_SHRINKAGE="${GRAPH_KRON_FISHER_SHRINKAGE:-0.20}"
GRAPH_KRON_FISHER_MIN_CLASS_COUNT="${GRAPH_KRON_FISHER_MIN_CLASS_COUNT:-25}"
GRAPH_KRON_FISHER_EVENT_PRIOR_WEIGHT="${GRAPH_KRON_FISHER_EVENT_PRIOR_WEIGHT:-0.35}"
GRAPH_KRON_ORPHAN_NEUTRAL_STRENGTH="${GRAPH_KRON_ORPHAN_NEUTRAL_STRENGTH:-0.0}"
GRAPH_KRON_ORPHAN_NEUTRAL_PT_SCALE="${GRAPH_KRON_ORPHAN_NEUTRAL_PT_SCALE:-15.0}"
GRAPH_KRON_ORPHAN_NEUTRAL_PILEUP_SCALE="${GRAPH_KRON_ORPHAN_NEUTRAL_PILEUP_SCALE:-1.0}"
GRAPH_KRON_SINK_STRENGTH="${GRAPH_KRON_SINK_STRENGTH:-0.20}"
GRAPH_KRON_ITERATIONS="${GRAPH_KRON_ITERATIONS:-20}"
GRAPH_KRON_PU_TERMINALS="${GRAPH_KRON_PU_TERMINALS:-12}"
GRAPH_KRON_WEIGHT_MIN="${GRAPH_KRON_WEIGHT_MIN:-0.03}"
GRAPH_KRON_JET_RADIUS="${GRAPH_KRON_JET_RADIUS:-0.40}"
GRAPH_KRON_CLUSTER_PT_MIN="${GRAPH_KRON_CLUSTER_PT_MIN:-5.0}"
GRAPH_KRON_OUTPUT_SCALE="${GRAPH_KRON_OUTPUT_SCALE:-0.45}"
GRAPH_KRON_NEUTRAL_PILEUP_SCALE="${GRAPH_KRON_NEUTRAL_PILEUP_SCALE:-1.0}"
GRAPH_KRON_SUBTRACTED_HARD_STRENGTH="${GRAPH_KRON_SUBTRACTED_HARD_STRENGTH:-0.0}"
GRAPH_KRON_MIN_JET_HARD_FRACTION="${GRAPH_KRON_MIN_JET_HARD_FRACTION:-0.70}"
GRAPH_KRON_MIN_JET_HARD_PT="${GRAPH_KRON_MIN_JET_HARD_PT:-18.0}"
GRAPH_KRON_RESIDUAL_PU_FRACTION="${GRAPH_KRON_RESIDUAL_PU_FRACTION:-0.35}"
GRAPH_KRON_RECOVERY_MIN_HARD_FRACTION="${GRAPH_KRON_RECOVERY_MIN_HARD_FRACTION:-0.20}"
GRAPH_KRON_RECOVERY_MIN_HARD_PT="${GRAPH_KRON_RECOVERY_MIN_HARD_PT:-15.0}"
GRAPH_KRON_RECOVERY_PT_MIN="${GRAPH_KRON_RECOVERY_PT_MIN:-15.0}"
GRAPH_KRON_RECOVERY_DEDUPE_RADIUS="${GRAPH_KRON_RECOVERY_DEDUPE_RADIUS:-0.28}"
GRAPH_KRON_RECOVERY_MAX_ADD="${GRAPH_KRON_RECOVERY_MAX_ADD:-1}"
GRAPH_KRON_RECOVERY_RESIDUAL_FRACTION="${GRAPH_KRON_RECOVERY_RESIDUAL_FRACTION:-0.35}"
GRAPH_KRON_RECOVERY_SCORE_PILEUP_PENALTY="${GRAPH_KRON_RECOVERY_SCORE_PILEUP_PENALTY:-0.50}"
GRAPH_KRON_RECOVERY_LOG_ODDS_MIN="${GRAPH_KRON_RECOVERY_LOG_ODDS_MIN:-2.0}"
GRAPH_KRON_RECOVERY_LIKELIHOOD_WEIGHT="${GRAPH_KRON_RECOVERY_LIKELIHOOD_WEIGHT:-1.0}"
GRAPH_KRON_RECOVERY_OUTPUT_SCALE="${GRAPH_KRON_RECOVERY_OUTPUT_SCALE:-0.45}"
GRAPH_KRON_LATENT_RECOVERY_RAW_PT_MIN="${GRAPH_KRON_LATENT_RECOVERY_RAW_PT_MIN:-70.0}"
GRAPH_KRON_LATENT_RECOVERY_MIN_HARD_PT="${GRAPH_KRON_LATENT_RECOVERY_MIN_HARD_PT:-1.0}"
GRAPH_KRON_LATENT_RECOVERY_MAX_ADD="${GRAPH_KRON_LATENT_RECOVERY_MAX_ADD:-1}"
GRAPH_KRON_LATENT_RECOVERY_OUTPUT_SCALE="${GRAPH_KRON_LATENT_RECOVERY_OUTPUT_SCALE:-0.30}"

mkdir -p "$OUTPUT_DIR"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-grdev}"

# Compute an absolute path without requiring the target to exist
abs_path() {
  local path="$1"
  local directory
  local basename
  directory="$(cd "$(dirname "$path")" && pwd)"
  basename="$(basename "$path")"
  printf '%s/%s\n' "$directory" "$basename"
}

# Require that a file exists before continuing
require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "$path" ]]; then
    echo "${label} not found: ${path}" >&2
    exit 1
  fi
}

# Require that an executable exists before continuing
require_executable() {
  local path="$1"
  local label="$2"
  if [[ ! -x "$path" ]]; then
    echo "${label} not executable: ${path}" >&2
    exit 1
  fi
}

# Rename an existing generated file before regenerating it
backup_existing_file() {
  local path="$1"
  local backup
  local index
  if [[ ! -e "$path" ]]; then
    return
  fi
  backup="${path}._old"
  index=1
  while [[ -e "$backup" ]]; do
    backup="${path}._old${index}"
    index=$((index + 1))
  done
  mv "$path" "$backup"
  echo "Renamed existing file: ${path} -> ${backup}"
}

# Source the Delphes runtime environment if the installation provides it
source_delphes_env() {
  if [[ -f "${DELPHES_DIR}/DelphesEnv.sh" ]]; then
    set +u
    source "${DELPHES_DIR}/DelphesEnv.sh"
    set -u
  fi
  if [[ -d "${DELPHES_DIR}/lib" ]]; then
    export LD_LIBRARY_PATH="${DELPHES_DIR}/lib:${LD_LIBRARY_PATH:-}"
  fi
  if [[ -d "${DELPHES_DIR}" ]]; then
    export LD_LIBRARY_PATH="${DELPHES_DIR}:${LD_LIBRARY_PATH:-}"
  fi
  if [[ -d "${DELPHES_BIN_DIR}" ]]; then
    export LD_LIBRARY_PATH="${DELPHES_BIN_DIR}:${LD_LIBRARY_PATH:-}"
  fi
}

# Compute the Delphes binary path for an in-source or installed layout
delphes_binary() {
  local name="$1"
  if [[ -x "${DELPHES_BIN_DIR}/${name}" ]]; then
    printf '%s/%s\n' "$DELPHES_BIN_DIR" "$name"
    return
  fi
  if [[ -x "${DELPHES_DIR}/bin/${name}" ]]; then
    printf '%s/bin/%s\n' "$DELPHES_DIR" "$name"
    return
  fi
  printf '%s/%s\n' "$DELPHES_DIR" "$name"
}

# Fail early if a Delphes binary has unresolved shared libraries
check_runtime_dependencies() {
  local binary="$1"
  local missing
  missing="$(ldd "$binary" | awk '/not found/ {print $1}')"
  if [[ -n "$missing" ]]; then
    echo "Missing shared libraries for ${binary}:" >&2
    echo "$missing" | sed 's/^/  /' >&2
    echo "Set DELPHES_DIR to a Delphes build matching the available ROOT runtime" >&2
    exit 2
  fi
}

# Run one C++ executable inside the graniitti conda environment
run_graniitti_binary() {
  run_in_graniitti "$@"
}

# Compute true when an existing Delphes ROOT file has graph-study inputs and event count
root_has_graph_inputs() {
  local path="$1"
  local expected_events="$2"
  if [[ ! -s "$path" ]]; then
    return 1
  fi
  run_in_graniitti python -c "import sys, uproot; tree = uproot.open(sys.argv[1], handler=uproot.source.file.MultithreadedFileSource)['Delphes']; expected = int(sys.argv[2]); required = ['Vertex.Z', 'EFlowTrackAll.PT', 'EFlowTrackAll.Charge', 'EFlowTrackAll.Z', 'EFlowTrackAll.DZ', 'EFlowPhoton.ET', 'EFlowNeutralHadron.ET', 'Particle.PT', 'Particle.Eta', 'Particle.Phi', 'Particle.Mass', 'Particle.PID', 'Particle.Status', 'Particle.IsPU']; ok = all(name in tree for name in required) and tree.num_entries == expected; sys.exit(0 if ok else 1)" "$path" "$expected_events" >/dev/null 2>&1
}

source_delphes_env

if [[ -z "$DELPHES_HEPMC3" ]]; then
  DELPHES_HEPMC3="$(delphes_binary DelphesHepMC3)"
fi
if [[ -z "$DELPHES_HEPMC2PILEUP" ]]; then
  DELPHES_HEPMC2PILEUP="$(delphes_binary hepmc2pileup)"
fi

require_file "$DELPHES_CARD_TEMPLATE" "Delphes card template"
require_executable "$DELPHES_HEPMC3" "DelphesHepMC3"
require_executable "$DELPHES_HEPMC2PILEUP" "hepmc2pileup"

bash tests/physics/studies/pileup_kron/build_pythia_driver.sh
bash tests/physics/studies/pileup_kron/build_hepmc_converter.sh

HARD_HEPMC="${OUTPUT_DIR}/hard_qcd_${NEVENTS}.hepmc3"
PU_TOTAL=$((NEVENTS * PU))
PU_HEPMC3="${OUTPUT_DIR}/minbias_${PU_TOTAL}.hepmc3"
PU_HEPMC2="${OUTPUT_DIR}/minbias_${PU_TOTAL}.hepmc2"
PILEUP_FILE="${OUTPUT_DIR}/minbias_${PU_TOTAL}.pileup"

# Match generated outputs to the cards, seeds and executable versions used
GENERATION_KEY="$(
    { printf '%s\n' "$NEVENTS" "$PU" "$SEED_HARD" "$SEED_PU" "$DELPHES_ROOT";
      sha256sum "$HARD_CMND" "$PU_CMND" "$DELPHES_CARD_TEMPLATE" \
        "$DELPHES_HEPMC3" "$DELPHES_HEPMC2PILEUP" \
        tests/physics/studies/pileup_kron/{pythia_hepmc3.cc,hepmc3_to_hepmc2.cc,prepare_delphes_card.py};
    } | sha256sum | cut -d ' ' -f 1
)"
GENERATION_RECORD="${DELPHES_ROOT}.sha256"
if [[ ! -f "$GENERATION_RECORD" || "$(cat "$GENERATION_RECORD")" != "$GENERATION_KEY" ]]; then
    REGENERATE=1
fi

if [[ "$REGENERATE" == "1" || ! -s "$HARD_HEPMC" ]]; then
  if [[ "$REGENERATE" == "1" ]]; then
    backup_existing_file "$HARD_HEPMC"
  fi
  run_graniitti_binary \
    env -u PYTHIA8DATA tests/physics/studies/pileup_kron/pythia_hepmc3 \
    "$HARD_CMND" "$HARD_HEPMC" "$NEVENTS" "$SEED_HARD"
fi

if [[ "$REGENERATE" == "1" || ! -s "$PU_HEPMC3" ]]; then
  if [[ "$REGENERATE" == "1" ]]; then
    backup_existing_file "$PU_HEPMC3"
  fi
  run_graniitti_binary \
    env -u PYTHIA8DATA tests/physics/studies/pileup_kron/pythia_hepmc3 \
    "$PU_CMND" "$PU_HEPMC3" "$PU_TOTAL" "$SEED_PU"
fi

if [[ "$REGENERATE" == "1" || ! -s "$PU_HEPMC2" ]]; then
  if [[ "$REGENERATE" == "1" ]]; then
    backup_existing_file "$PU_HEPMC2"
  fi
  run_graniitti_binary tests/physics/studies/pileup_kron/hepmc3_to_hepmc2 "$PU_HEPMC3" "$PU_HEPMC2"
fi

if [[ "$REGENERATE" == "1" || ! -s "$PILEUP_FILE" ]]; then
  if [[ "$REGENERATE" == "1" ]]; then
    backup_existing_file "$PILEUP_FILE"
  fi
  if [[ "$CHECK_DELPHES_RUNTIME" == "1" ]]; then
    check_runtime_dependencies "$DELPHES_HEPMC2PILEUP"
  fi
  run_graniitti_binary "$DELPHES_HEPMC2PILEUP" "$PILEUP_FILE" "$PU_HEPMC2"
fi

if [[ "$REGENERATE" == "1" ]]; then
  backup_existing_file "$DELPHES_CARD"
fi
python tests/physics/studies/pileup_kron/prepare_delphes_card.py \
  --template "$DELPHES_CARD_TEMPLATE" \
  --output "$DELPHES_CARD" \
  --pileup-file "$(abs_path "$PILEUP_FILE")" \
  --mean-pu "$PU" \
  --fixed-pu

if [[ -s "$DELPHES_ROOT" ]] && ! root_has_graph_inputs "$DELPHES_ROOT" "$NEVENTS"; then
  backup_existing_file "$DELPHES_ROOT"
fi

if [[ "$REGENERATE" == "1" || ! -s "$DELPHES_ROOT" ]]; then
  if [[ "$REGENERATE" == "1" ]]; then
    backup_existing_file "$DELPHES_ROOT"
  fi
  if [[ "$CHECK_DELPHES_RUNTIME" == "1" ]]; then
    check_runtime_dependencies "$DELPHES_HEPMC3"
  fi
  run_graniitti_binary "$DELPHES_HEPMC3" "$DELPHES_CARD" "$DELPHES_ROOT" "$HARD_HEPMC"
fi

printf '%s\n' "$GENERATION_KEY" > "$GENERATION_RECORD"

ANALYZE_ARGS=(
  python tests/physics/studies/pileup_kron/analyze_delphes.py
  --input "$DELPHES_ROOT"
  --output-dir "$OUTPUT_DIR"
  --reference "$REFERENCE"
  --collections Jet JetPUPPI
  --match-radius "$MATCH_RADIUS"
  --gen-pt-min "$GEN_PT_MIN"
  --reco-pt-min "$RECO_PT_MIN"
  --eta-max "$ETA_MAX"
  --event-start "$EVENT_START"
  --max-events "$MAX_EVENTS"
)

if [[ -n "$EVENT_INDICES" ]]; then
  ANALYZE_ARGS+=(--event-indices "$EVENT_INDICES")
fi

if [[ -n "$EVENT_INDICES_FILE" ]]; then
  ANALYZE_ARGS+=(--event-indices-file "$EVENT_INDICES_FILE")
fi

if [[ "$INCLUDE_GRAPH_KRON" == "1" ]]; then
  ANALYZE_ARGS+=(
    --include-graph-kron
    --graph-kron-collection "$GRAPH_KRON_COLLECTION"
    --graph-kron-base-axis "$GRAPH_KRON_BASE_AXIS"
    --graph-kron-cluster-mode "$GRAPH_KRON_CLUSTER_MODE"
    --graph-kron-candidate-pt-min "$GRAPH_KRON_CANDIDATE_PT_MIN"
    --graph-kron-max-candidates "$GRAPH_KRON_MAX_CANDIDATES"
    --graph-kron-graph-radius "$GRAPH_KRON_GRAPH_RADIUS"
    --graph-kron-graph-sigma "$GRAPH_KRON_GRAPH_SIGMA"
    --graph-kron-radiation-radius "$GRAPH_KRON_RADIATION_RADIUS"
    --graph-kron-radiation-core "$GRAPH_KRON_RADIATION_CORE"
    --graph-kron-radiation-angular-power "$GRAPH_KRON_RADIATION_ANGULAR_POWER"
    --graph-kron-radiation-pt-power "$GRAPH_KRON_RADIATION_PT_POWER"
    --graph-kron-radiation-edge-strength "$GRAPH_KRON_RADIATION_EDGE_STRENGTH"
    --graph-kron-radiation-edge-max-neighbours "$GRAPH_KRON_RADIATION_EDGE_MAX_NEIGHBOURS"
    --graph-kron-radiation-unary-strength "$GRAPH_KRON_RADIATION_UNARY_STRENGTH"
    --graph-kron-anchor-strength "$GRAPH_KRON_ANCHOR_STRENGTH"
    --graph-kron-local-support-strength "$GRAPH_KRON_LOCAL_SUPPORT_STRENGTH"
    --graph-kron-fisher-mode "$GRAPH_KRON_FISHER_MODE"
    --graph-kron-fisher-strength "$GRAPH_KRON_FISHER_STRENGTH"
    --graph-kron-fisher-shrinkage "$GRAPH_KRON_FISHER_SHRINKAGE"
    --graph-kron-fisher-min-class-count "$GRAPH_KRON_FISHER_MIN_CLASS_COUNT"
    --graph-kron-fisher-event-prior-weight "$GRAPH_KRON_FISHER_EVENT_PRIOR_WEIGHT"
    --graph-kron-orphan-neutral-strength "$GRAPH_KRON_ORPHAN_NEUTRAL_STRENGTH"
    --graph-kron-orphan-neutral-pt-scale "$GRAPH_KRON_ORPHAN_NEUTRAL_PT_SCALE"
    --graph-kron-orphan-neutral-pileup-scale "$GRAPH_KRON_ORPHAN_NEUTRAL_PILEUP_SCALE"
    --graph-kron-sink-strength "$GRAPH_KRON_SINK_STRENGTH"
    --graph-kron-iterations "$GRAPH_KRON_ITERATIONS"
    --graph-kron-pu-terminals "$GRAPH_KRON_PU_TERMINALS"
    --graph-kron-weight-min "$GRAPH_KRON_WEIGHT_MIN"
    --graph-kron-jet-radius "$GRAPH_KRON_JET_RADIUS"
    --graph-kron-cluster-pt-min "$GRAPH_KRON_CLUSTER_PT_MIN"
    --graph-kron-output-scale "$GRAPH_KRON_OUTPUT_SCALE"
    --graph-kron-neutral-pileup-scale "$GRAPH_KRON_NEUTRAL_PILEUP_SCALE"
    --graph-kron-subtracted-hard-strength "$GRAPH_KRON_SUBTRACTED_HARD_STRENGTH"
    --graph-kron-min-jet-hard-fraction "$GRAPH_KRON_MIN_JET_HARD_FRACTION"
    --graph-kron-min-jet-hard-pt "$GRAPH_KRON_MIN_JET_HARD_PT"
    --graph-kron-residual-pu-fraction "$GRAPH_KRON_RESIDUAL_PU_FRACTION"
    --graph-kron-recovery-min-hard-fraction "$GRAPH_KRON_RECOVERY_MIN_HARD_FRACTION"
    --graph-kron-recovery-min-hard-pt "$GRAPH_KRON_RECOVERY_MIN_HARD_PT"
    --graph-kron-recovery-pt-min "$GRAPH_KRON_RECOVERY_PT_MIN"
    --graph-kron-recovery-dedupe-radius "$GRAPH_KRON_RECOVERY_DEDUPE_RADIUS"
    --graph-kron-recovery-max-add "$GRAPH_KRON_RECOVERY_MAX_ADD"
    --graph-kron-recovery-residual-fraction "$GRAPH_KRON_RECOVERY_RESIDUAL_FRACTION"
    --graph-kron-recovery-score-pileup-penalty "$GRAPH_KRON_RECOVERY_SCORE_PILEUP_PENALTY"
    --graph-kron-recovery-log-odds-min "$GRAPH_KRON_RECOVERY_LOG_ODDS_MIN"
    --graph-kron-recovery-likelihood-weight "$GRAPH_KRON_RECOVERY_LIKELIHOOD_WEIGHT"
    --graph-kron-recovery-output-scale "$GRAPH_KRON_RECOVERY_OUTPUT_SCALE"
    --graph-kron-latent-recovery-raw-pt-min "$GRAPH_KRON_LATENT_RECOVERY_RAW_PT_MIN"
    --graph-kron-latent-recovery-min-hard-pt "$GRAPH_KRON_LATENT_RECOVERY_MIN_HARD_PT"
    --graph-kron-latent-recovery-max-add "$GRAPH_KRON_LATENT_RECOVERY_MAX_ADD"
    --graph-kron-latent-recovery-output-scale "$GRAPH_KRON_LATENT_RECOVERY_OUTPUT_SCALE"
  )
fi

run_in_graniitti "${ANALYZE_ARGS[@]}"

echo "Delphes card:      ${DELPHES_CARD}"
echo "Delphes ROOT:      ${DELPHES_ROOT}"
echo "Delphes metrics:   ${OUTPUT_DIR}/delphes_metrics_summary.md"
echo "Matched jets:      ${OUTPUT_DIR}/delphes_matched_jets.csv"

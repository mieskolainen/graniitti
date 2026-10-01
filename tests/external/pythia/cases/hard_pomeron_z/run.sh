#!/usr/bin/env bash
set -euo pipefail

# Generate the 13 TeV hard Pomeron Z workflow samples
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
cd "$REPO_ROOT"
# Print command line usage for the hard Pomeron Z case
usage() {
  cat <<'EOF'
Usage: bash tests/external/pythia/cases/hard_pomeron_z/run.sh [options]

Options:
  --hard-muon-cuts on|off  Toggle the native Pythia hard process dimuon filter
                            (default: on)
  -h, --help                Show this help

Environment:
  NEVENTS=N                 Number of generated events
  LOOPSCREEN=0|1               Pomeron loop screening for GRANIITTI samples
  WEIGHTED=0|1              Weighted GRANIITTI event generation
EOF
}

# Require a GRANIITTI binary switch to be 0 or 1
validate_binary_switch() {
  local name="$1"
  local value="$2"

  if [[ ! "$value" =~ ^[01]$ ]]; then
    echo "${name} must be 0 or 1: ${value}" >&2
    exit 2
  fi
}

HARD_MUON_CUTS="${HARD_MUON_CUTS:-on}"
HARD_MUON_CUTS_SEEN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --hard-muon-cuts)
      if [[ "$HARD_MUON_CUTS_SEEN" -eq 1 ]]; then
        echo "--hard-muon-cuts was specified more than once" >&2
        exit 2
      fi
      if [[ $# -lt 2 ]]; then
        echo "--hard-muon-cuts requires on or off" >&2
        usage >&2
        exit 2
      fi
      HARD_MUON_CUTS="$2"
      HARD_MUON_CUTS_SEEN=1
      shift 2
      ;;
    --hard-muon-cuts=*)
      if [[ "$HARD_MUON_CUTS_SEEN" -eq 1 ]]; then
        echo "--hard-muon-cuts was specified more than once" >&2
        exit 2
      fi
      HARD_MUON_CUTS="${1#*=}"
      HARD_MUON_CUTS_SEEN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ "$HARD_MUON_CUTS" != "on" && "$HARD_MUON_CUTS" != "off" ]]; then
  echo "--hard-muon-cuts must be either on or off" >&2
  exit 2
fi

set +u
source tests/environment.sh
activate_conda_environment graniitti
source install/setenv.sh
set -u

NEVENTS="${NEVENTS:-10000}"
SEED="${SEED:-12345}"
LOOPSCREEN="${LOOPSCREEN:-0}"
WEIGHTED="${WEIGHTED:-1}"
PYTHIA_CMND="${PYTHIA_CMND-tests/external/pythia/drivers/lhe_converter/shower.cmnd}"
RESHOWER_ATTEMPTS="${RESHOWER_ATTEMPTS:-500}"

validate_binary_switch "LOOPSCREEN" "$LOOPSCREEN"
validate_binary_switch "WEIGHTED" "$WEIGHTED"

IPP_CARD="${IPP_CARD:-icepack/HARDPOM/z_with_pythia/gencard_IPp_Z_mumu.json}"
IPIP_CARD="${IPIP_CARD:-icepack/HARDPOM/z_with_pythia/gencard_IPIP_Z_mumu.json}"
DY_CMND="${DY_CMND:-icepack/HARDPOM/z_with_pythia/pythia_inclusive.cmnd}"
PYTHIA_HARD_DIFF_CMND="${PYTHIA_HARD_DIFF_CMND:-icepack/HARDPOM/z_with_pythia/pythia_hard_diffraction.cmnd}"
PYTHIA_HARD_DIFF_NOMPI_CMND="${PYTHIA_HARD_DIFF_NOMPI_CMND:-icepack/HARDPOM/z_with_pythia/pythia_hard_diffraction_no_mpi.cmnd}"

IPP_TAG="${IPP_TAG:-hard_pomeron_z_graniitti_IPp}"
IPIP_TAG="${IPIP_TAG:-hard_pomeron_z_graniitti_IPIP}"
DY_TAG="${DY_TAG:-hard_pomeron_z_pythia_inclusive}"
PYTHIA_HARD_DIFF_TAG="${PYTHIA_HARD_DIFF_TAG:-hard_pomeron_z_pythia_hard_diffraction}"
PYTHIA_HARD_DIFF_NOMPI_TAG="${PYTHIA_HARD_DIFF_NOMPI_TAG:-hard_pomeron_z_pythia_hard_diffraction_no_mpi}"

mkdir -p output

# Print a compact section header for workflow logs
section() {
  printf '\n==== %s ====\n' "$1"
}

# Generate one GRANIITTI LHE sample and convert it to HepMC3
generate_graniitti_lhe_pythia() {
  local card="$1"
  local tag="$2"
  NEVENTS="$NEVENTS" SEED="$SEED" LOOPSCREEN="$LOOPSCREEN" WEIGHTED="$WEIGHTED" \
    PYTHIA_CMND="$PYTHIA_CMND" RESHOWER_ATTEMPTS="$RESHOWER_ATTEMPTS" \
    bash tests/external/pythia/drivers/lhe_converter/run.sh "$card" "$tag"
}

# Generate gamma*/Z to muons directly with Pythia
generate_pythia_zmumu() {
  local cmnd="$1"
  local tag="$2"
  local title="$3"
  local hepmc_file="output/${tag}.hepmc3"

  section "Pythia ${title} ${tag}"
  run_in_graniitti \
    bin/pythia_zmumu_hepmc3 \
      "$cmnd" "$hepmc_file" "$NEVENTS" "$SEED" \
      --hard-muon-cuts "$HARD_MUON_CUTS"

  echo "Pythia output:  ${hepmc_file}"
}

section "Build Pythia helpers"
echo "Native Pythia hard process muon cuts: ${HARD_MUON_CUTS}"
echo "GRANIITTI LOOPSCREEN: ${LOOPSCREEN}"
echo "GRANIITTI WEIGHTED: ${WEIGHTED}"
if [[ ! -x bin/pythia_lhe_hadronize ]]; then
  bash tests/external/pythia/drivers/lhe_converter/build.sh
fi
if [[ ! -x bin/pythia_zmumu_hepmc3 ]]; then
  bash tests/external/pythia/drivers/pythia_drell_yan/build.sh
fi

generate_graniitti_lhe_pythia "$IPP_CARD" "$IPP_TAG"
generate_graniitti_lhe_pythia "$IPIP_CARD" "$IPIP_TAG"
generate_pythia_zmumu "$DY_CMND" "$DY_TAG" "inclusive DY"
generate_pythia_zmumu \
  "$PYTHIA_HARD_DIFF_CMND" "$PYTHIA_HARD_DIFF_TAG" "hard diffractive DY"
generate_pythia_zmumu \
  "$PYTHIA_HARD_DIFF_NOMPI_CMND" "$PYTHIA_HARD_DIFF_NOMPI_TAG" \
  "hard diffractive DY without MPI"

section "Generated HepMC3 files"
printf '  %s\n' \
  "output/${IPP_TAG}.hepmc3" \
  "output/${IPIP_TAG}.hepmc3" \
  "output/${DY_TAG}.hepmc3" \
  "output/${PYTHIA_HARD_DIFF_TAG}.hepmc3" \
  "output/${PYTHIA_HARD_DIFF_NOMPI_TAG}.hepmc3"

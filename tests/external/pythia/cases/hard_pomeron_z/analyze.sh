#!/usr/bin/env bash
set -euo pipefail

# Analyze the hard Pomeron Z samples with direct physics observables
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
cd "$REPO_ROOT"
set +u
source tests/environment.sh
activate_conda_environment graniitti
source install/setenv.sh
set -u

IPP_TAG="${IPP_TAG:-hard_pomeron_z_graniitti_IPp}"
IPIP_TAG="${IPIP_TAG:-hard_pomeron_z_graniitti_IPIP}"
DY_TAG="${DY_TAG:-hard_pomeron_z_pythia_inclusive}"
PYTHIA_HARD_DIFF_TAG="${PYTHIA_HARD_DIFF_TAG:-hard_pomeron_z_pythia_hard_diffraction}"
PYTHIA_HARD_DIFF_NOMPI_TAG="${PYTHIA_HARD_DIFF_NOMPI_TAG:-hard_pomeron_z_pythia_hard_diffraction_no_mpi}"
OUTPUT_TAG="${OUTPUT_TAG:-hard_pomeron_z}"
CORES="${CORES:-0}"
LEGEND="${LEGEND:-outside-right}"
LEGEND_NCOL="${LEGEND_NCOL:-1}"
LEGEND_FONTSIZE="${LEGEND_FONTSIZE:-6}"
LEGEND_BBOX="${LEGEND_BBOX:-0.70 0.58}"
MPLCONFIGDIR="${MPLCONFIGDIR:-${REPO_ROOT}/tmp/matplotlib}"
XDG_CACHE_HOME="${XDG_CACHE_HOME:-${REPO_ROOT}/tmp/cache}"

export MPLCONFIGDIR
export XDG_CACHE_HOME
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

SAMPLES=(
  "$IPP_TAG"
  "$IPIP_TAG"
  "$DY_TAG"
  "$PYTHIA_HARD_DIFF_TAG"
  "$PYTHIA_HARD_DIFF_NOMPI_TAG"
)
LABELS=(
  'GRANIITTI $IPp[Z]$'
  'GRANIITTI $IPIP[Z]$'
  'Pythia 8 inclusive $Z/\gamma^{*}$'
  'Pythia 8 hard SD with MPI'
  'Pythia 8 hard SD without MPI'
)
CUTS=()
for _ in "${SAMPLES[@]}"; do
  CUTS+=(icepack/HARDPOM/z_with_pythia/cuts.py)
done

# Require one generated HepMC3 file before starting the analysis
require_hepmc3() {
  local sample="$1"
  local file="output/${sample}.hepmc3"
  if [[ ! -s "$file" ]]; then
    echo "Missing or empty HepMC3 file: ${file}" >&2
    exit 1
  fi
}

for sample in "${SAMPLES[@]}"; do
  require_hepmc3 "$sample"
done

ICE_ARGS=(
  -m core.iceplot
  --hepmc3 "${SAMPLES[@]}"
  --obs icepack/HARDPOM/z_with_pythia/obs.py
  --cuts "${CUTS[@]}"
  --pid '[[13,-13]]'
  --mclabel "${LABELS[@]}"
  --unit pb
  --output "$OUTPUT_TAG"
  --title '(13 TeV)'
  --title_loc left
  --legend "$LEGEND"
  --legend-ncol "$LEGEND_NCOL"
  --legend-fontsize "$LEGEND_FONTSIZE"
  --cores "$CORES"
)

if [[ -n "$LEGEND_BBOX" ]]; then
  read -r LEGEND_BBOX_X LEGEND_BBOX_Y <<< "$LEGEND_BBOX"
  ICE_ARGS+=(--legend-bbox "$LEGEND_BBOX_X" "$LEGEND_BBOX_Y")
fi

ICE_ARGS+=("$@")

run_in_graniitti python "${ICE_ARGS[@]}"

PLOT_OUTPUT_TAG="$OUTPUT_TAG"
for argument in "$@"; do
  if [[ "$argument" == "--density" ]]; then
    PLOT_OUTPUT_TAG="${OUTPUT_TAG}__[density]"
  fi
done
echo "iceplot plots: figs/iceplot/${PLOT_OUTPUT_TAG}"

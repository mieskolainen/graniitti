#!/usr/bin/env bash
set -euo pipefail

# Compare quark and gluon jets after particle level anti-kT clustering
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
cd "$REPO_ROOT"
set +u
source tests/environment.sh
activate_conda_environment graniitti
source install/setenv.sh
set -u

GG_TAG="${GG_TAG:-durham_qcd_gg_to_gg}"
UUBAR_TAG="${UUBAR_TAG:-durham_qcd_gg_to_uubar}"
OUTPUT_TAG="${OUTPUT_TAG:-durham_qcd_jets}"
MPLCONFIGDIR="${MPLCONFIGDIR:-${REPO_ROOT}/tmp/matplotlib}"
XDG_CACHE_HOME="${XDG_CACHE_HOME:-${REPO_ROOT}/tmp/cache}"

export MPLCONFIGDIR
export XDG_CACHE_HOME
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

SAMPLES=("$GG_TAG" "$UUBAR_TAG")
LABELS=('GRANIITTI+Pythia $gg \to gg$ (Durham)' 'GRANIITTI+Pythia $gg \to u\bar{u}$ (Durham)')
CUTS=()
for _ in "${SAMPLES[@]}"; do
  CUTS+=(icepack/DURHAM/dijet_with_pythia/cuts.py)
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

run_in_graniitti python -m core.iceplot \
  --hepmc3 "${SAMPLES[@]}" \
  --obs icepack/DURHAM/dijet_with_pythia/obs.py \
  --cuts "${CUTS[@]}" \
  --pid '[[0]]' \
  --mclabel "${LABELS[@]}" \
  --unit pb \
  --output "$OUTPUT_TAG" \
  --title '(13 TeV)' \
  --title_loc left \
  "$@"

echo "iceplot plots: figs/iceplot/${OUTPUT_TAG}"

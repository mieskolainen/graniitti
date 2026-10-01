#!/usr/bin/env bash
set -euo pipefail

# Plot dimuon and forward excitation particle level observables
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
cd "$REPO_ROOT"
set +u
source tests/environment.sh
activate_conda_environment graniitti
source install/setenv.sh
set -u

TAG="${TAG:-yy_mumu_dissociative}"
OUTPUT_TAG="${OUTPUT_TAG:-yy_mumu_dissociative}"
CORES="${CORES:-0}"
MPLCONFIGDIR="${MPLCONFIGDIR:-${REPO_ROOT}/tmp/matplotlib}"
XDG_CACHE_HOME="${XDG_CACHE_HOME:-${REPO_ROOT}/tmp/cache}"

export MPLCONFIGDIR
export XDG_CACHE_HOME
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

if [[ ! -s "output/${TAG}.hepmc3" ]]; then
  echo "Missing or empty HepMC3 file: output/${TAG}.hepmc3" >&2
  exit 1
fi

run_in_graniitti python -m core.iceplot \
  --hepmc3 "$TAG" \
  --obs icepack/GAMMA/mumu_with_pythia/obs.py \
  --cuts icepack/GAMMA/mumu_with_pythia/cuts.py \
  --pid '[[13,-13]]' \
  --mclabel 'GRANIITTI $\gamma\gamma \to \mu^+\mu^-$ with one excited proton' \
  --unit pb \
  --output "$OUTPUT_TAG" \
  --title '(13 TeV)' \
  --title_loc left \
  --cores "$CORES" \
  "$@"

echo "iceplot plots: figs/iceplot/${OUTPUT_TAG}"

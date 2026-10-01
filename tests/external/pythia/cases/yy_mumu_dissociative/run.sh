#!/usr/bin/env bash
set -euo pipefail

# Generate gamma gamma to muons with one string fragmented forward excitation
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
cd "$REPO_ROOT"
set +u
source tests/environment.sh
activate_conda_environment graniitti
source install/setenv.sh
set -u

export CONVERTER_MODE="${CONVERTER_MODE:-fragment}"

TAG="${TAG:-yy_mumu_dissociative}"
export LOOPSCREEN="${LOOPSCREEN:-0}"
export WEIGHTED="${WEIGHTED:-1}"

bash tests/external/pythia/drivers/lhe_converter/run.sh \
  icepack/GAMMA/mumu_with_pythia/gencard.json \
  "$TAG" \
  --set 'SCATTERING.NSTARS=1' \
  --set 'SCATTERING.BEAMFRAG="diquark"' \
  "$@"

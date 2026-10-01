#!/usr/bin/env bash
set -euo pipefail

# Generate elastic Durham gg to gg and gg to u ubar samples
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
cd "$REPO_ROOT"
set +u
source tests/environment.sh
activate_conda_environment graniitti
source install/setenv.sh
set -u

export CONVERTER_MODE="${CONVERTER_MODE:-fragment}"

CARD="${CARD:-icepack/DURHAM/dijet_with_pythia/gencard.json}"
GG_TAG="${GG_TAG:-durham_qcd_gg_to_gg}"
UUBAR_TAG="${UUBAR_TAG:-durham_qcd_gg_to_uubar}"
export LOOPSCREEN="${LOOPSCREEN:-1}"
export WEIGHTED="${WEIGHTED:-1}"

bash tests/external/pythia/drivers/lhe_converter/run.sh \
  "$CARD" \
  "$GG_TAG" \
  --set 'SCATTERING.PROCESS="gg[QCD]<F> -> g g"' \
  --set 'SCATTERING.NSTARS=0' \
  "$@"

bash tests/external/pythia/drivers/lhe_converter/run.sh \
  "$CARD" \
  "$UUBAR_TAG" \
  --set 'SCATTERING.PROCESS="gg[QCD]<F> -> u u~"' \
  --set 'SCATTERING.NSTARS=0' \
  "$@"

#!/usr/bin/env bash
set -eo pipefail
# Run the SD merging icepack from the repository root
cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.."
source tests/environment.sh
activate_conda_environment graniitti
source install/setenv.sh
[[ -x bin/gr ]] || { echo "Build GRANIITTI first with cmake --build build -j4" >&2; exit 1; }
[[ -x bin/pythia_lhe_hadronize ]] || bash tests/external/pythia/drivers/lhe_converter/build.sh
[[ -x bin/pythia_zmumu_hepmc3 ]] || bash tests/external/pythia/drivers/pythia_drell_yan/build.sh
export MPLCONFIGDIR="$PWD/tmp/matplotlib"
export XDG_CACHE_HOME="$PWD/tmp/cache"
python tests/external/pythia/cases/z_merging/study.py "$@"

#!/usr/bin/env bash
set -euo pipefail
# Compile one Pythia driver with absolute library and XML data paths
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
cd "$REPO_ROOT"
set +u
source tests/environment.sh
activate_conda_environment graniitti
source install/setenv.sh
set -u
PYTHIA_DIR="$(cd "${PYTHIA_DIR:-../pythia8317}" && pwd)"
if [[ ! -f "$PYTHIA_DIR/lib/libpythia8.a" && ! -f "$PYTHIA_DIR/lib/libpythia8.so" ]]; then
  echo "Build Pythia first with bash tests/external/pythia/build_pythia.sh" >&2
  exit 1
fi
mkdir -p "$(dirname "$2")"
"${CXX:-g++}" -std=c++17 -O2 "$1" -o "$2" \
  "-DPYTHIA_XML_DIR=\"$PYTHIA_DIR/share/Pythia8/xmldoc\"" \
  -I"$PYTHIA_DIR/include" -I"$HEPMC3SYS/include" \
  -L"$PYTHIA_DIR/lib" -L"$HEPMC3SYS/lib" -L"$HEPMC3SYS/lib64" \
  -Wl,-rpath,"$PYTHIA_DIR/lib" -Wl,-rpath,"$HEPMC3SYS/lib" -Wl,-rpath,"$HEPMC3SYS/lib64" \
  -lpythia8 -lHepMC3 -ldl -pthread
echo "Built $2"

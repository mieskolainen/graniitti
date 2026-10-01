#!/usr/bin/env bash
set -euo pipefail

# Compile the generic Pythia-to-HepMC3 driver
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
cd "$REPO_ROOT"
set +u
source tests/environment.sh
activate_conda_environment graniitti
source install/setenv.sh
set -u

PYTHIA_DIR="${PYTHIA_DIR:-../pythia8317}"
CXX="${CXX:-g++}"
OUTPUT="${OUTPUT:-tests/physics/studies/pileup_kron/pythia_hepmc3}"

if [[ ! -f "${PYTHIA_DIR}/lib/libpythia8.a" && ! -f "${PYTHIA_DIR}/lib/libpythia8.so" ]]; then
  echo "Pythia library not found under ${PYTHIA_DIR}/lib" >&2
  echo "Run tests/external/pythia/build_pythia.sh first" >&2
  exit 1
fi

"$CXX" -std=c++17 -O2 \
  tests/physics/studies/pileup_kron/pythia_hepmc3.cc \
  -o "$OUTPUT" \
  -I"${PYTHIA_DIR}/include" \
  -I"${HEPMC3SYS}/include" \
  -L"${PYTHIA_DIR}/lib" \
  -L"${HEPMC3SYS}/lib" \
  -L"${HEPMC3SYS}/lib64" \
  -Wl,-rpath,"${PYTHIA_DIR}/lib" \
  -Wl,-rpath,"${HEPMC3SYS}/lib" \
  -Wl,-rpath,"${HEPMC3SYS}/lib64" \
  -lpythia8 -lHepMC3 -lHepMC3search -ldl -pthread

echo "Built ${OUTPUT}"

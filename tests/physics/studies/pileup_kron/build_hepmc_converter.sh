#!/usr/bin/env bash
set -euo pipefail

# Compile the HepMC3-to-HepMC2 converter used by Delphes pile-up preparation
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
cd "$REPO_ROOT"
set +u
source tests/environment.sh
activate_conda_environment graniitti
source install/setenv.sh
set -u

CXX="${CXX:-g++}"
OUTPUT="${OUTPUT:-tests/physics/studies/pileup_kron/hepmc3_to_hepmc2}"

"$CXX" -std=c++17 -O2 \
  tests/physics/studies/pileup_kron/hepmc3_to_hepmc2.cc \
  -o "$OUTPUT" \
  -Iinclude \
  -I"${HEPMC3SYS}/include" \
  -L"${HEPMC3SYS}/lib" \
  -L"${HEPMC3SYS}/lib64" \
  -Wl,-rpath,"${HEPMC3SYS}/lib" \
  -Wl,-rpath,"${HEPMC3SYS}/lib64" \
  -lHepMC3

echo "Built ${OUTPUT}"

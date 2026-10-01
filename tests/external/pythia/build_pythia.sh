#!/usr/bin/env bash
set -euo pipefail

# Build the local Pythia8 tree against the GRANIITTI HepMC3 and LHAPDF stack
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO_ROOT"
JOBS="${JOBS:-4}"
if [[ ! "$JOBS" =~ ^[1-4]$ ]]; then
  echo 'JOBS must be between 1 and 4' >&2
  exit 64
fi

set +u
source tests/environment.sh
activate_conda_environment graniitti
source install/setenv.sh
set -u

PYTHIA_DIR="$(cd "${PYTHIA_DIR:-../pythia8317}" && pwd)"

if [[ ! -x "${PYTHIA_DIR}/configure" ]]; then
  echo "Pythia configure script not found under ${PYTHIA_DIR}" >&2
  exit 1
fi

cd "$PYTHIA_DIR"
./configure \
  --prefix="$PYTHIA_DIR" \
  --with-hepmc3="$HEPMC3SYS" \
  --with-lhapdf6="$LHAPDFSYS" \
  --with-gzip

make -j "$JOBS"
make install

echo "Pythia built under ${PYTHIA_DIR}"

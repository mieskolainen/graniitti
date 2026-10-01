#!/usr/bin/env bash
set -euo pipefail
# Build this Pythia driver with the shared compiler and library settings
DRIVER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$DRIVER_DIR/../build.sh" "$DRIVER_DIR/pythia_zmumu_hepmc3.cc" "${OUTPUT:-bin/pythia_zmumu_hepmc3}"

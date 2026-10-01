#!/usr/bin/env bash
set -euo pipefail
# Run the cms Rivet comparison through the common driver
exec bash "$(dirname "${BASH_SOURCE[0]}")/../run.sh" cms "$@"

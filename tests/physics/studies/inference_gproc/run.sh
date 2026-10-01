#!/usr/bin/env bash
# Run Gaussian process input and target uncertainty studies
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize

if ! study_generates; then
    study_skip "Gaussian process study has no separate analysis stage"
fi
exec python "${STUDY_LAUNCHER_DIR}/gaussian_process.py"

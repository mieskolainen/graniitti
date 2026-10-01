#!/usr/bin/env bash
#
# Submit a GRANIITTI test suite
#
# Run with: bash tests/condor/submit.sh
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

set -euo pipefail

REPO_ROOT="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
exec bash "${REPO_ROOT}/tests/condor/run_job.sh" "${REPO_ROOT}" submit "$@"

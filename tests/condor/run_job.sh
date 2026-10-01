#!/usr/bin/env bash
#
# Enter the GRANIITTI environment for submission, execution or reporting
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

set -euo pipefail

REPO_ROOT="${1:?missing repository root}"
shift
MODULE="${1:?missing Condor module}"
shift
case "${MODULE}" in
  submit|job_runner|verify|icepacks) ;;
  *) echo "Unknown Condor module: ${MODULE}" >&2; exit 64 ;;
esac

if [[ ! -f "${REPO_ROOT}/install/setenv.sh" ]]; then
  echo "Repository root is not accessible: ${REPO_ROOT}" >&2
  exit 64
fi
cd "${REPO_ROOT}"

source tests/environment.sh
set +u
activate_conda_environment graniitti
if [[ -z "${GRANIITTI_ENV:-}" ]]; then
  source install/setenv.sh
fi
set -u

exec python -m "tests.condor.${MODULE}" "$@"

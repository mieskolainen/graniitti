#!/usr/bin/env bash
set -euo pipefail

# Run the reset Delphes baseline workflow
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize 300
study_require_opt_in RUN_DELPHES "pileup Kron Delphes study"
if [[ -z "${DELPHES_DIR:-}" ]]
then
	study_skip "pileup Kron Delphes study requires DELPHES_DIR"
fi
study_require_files "pileup Kron Delphes study" \
	"${DELPHES_DIR}/cards/CMS_PhaseII/CMS_PhaseII_200PU_v03_nodtf.tcl"
if ! study_generates
then
	study_skip "pileup Kron Delphes study has no separate analysis stage"
fi
exec bash tests/physics/studies/pileup_kron/run_delphes.sh "$@"

#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Compare single and double channel eikonal models across the pp and ppbar energy scans
#
# Run with: bash tests/physics/studies/multichannel/run.sh
# Plot existing outputs with: STUDY_ACTION=analyze bash tests/physics/studies/multichannel/run.sh

read -r -a PP_ENERGIES <<< "${SQRTS:-200 7000 13000 60000}"
PPBAR_ENERGIES=(546 1960)

SCAN_ROOT="${STUDY_REPO_ROOT}/tmp/multichannel"
for MODEL in single double
do
    EIKONAL_MODEL="$MODEL" SCAN_DIR="$SCAN_ROOT/$MODEL" \
        FIG_DIR="${STUDY_REPO_ROOT}/figs/multichannel" \
        PP_ENERGIES="$(IFS=,; echo "${PP_ENERGIES[*]}")" \
        PPBAR_ENERGIES="$(IFS=,; echo "${PPBAR_ENERGIES[*]}")" \
        bash "${STUDY_LAUNCHER_DIR}/../elastic/run.sh"
done

for BEAM in pp ppbar
do
    BEAM2=2212
    ENERGIES=("${PP_ENERGIES[@]}")
    if [[ "$BEAM" == ppbar ]]; then
        BEAM2=-2212
        ENERGIES=("${PPBAR_ENERGIES[@]}")
    fi
    python -m tests.physics.studies.multichannel.analysis.analyze \
        --scan-dir "$SCAN_ROOT" --models single double --sqrts "${ENERGIES[@]}" \
        --beam2 "$BEAM2" --fig-dir "${STUDY_REPO_ROOT}/figs/multichannel/$BEAM"
done

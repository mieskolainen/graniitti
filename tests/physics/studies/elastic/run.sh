#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Elastic pp scattering at different energies
#
# Generate the scan and plots with: bash ./tests/physics/studies/elastic/run.sh
# Select eikonal model with: EIKONAL_MODEL=double bash ./tests/physics/studies/elastic/run.sh

# Screening loop must be on, set with -l true !

EIKONAL_MODEL=${EIKONAL_MODEL:-single}

case "$EIKONAL_MODEL" in
single|double|triple)
	;;
*)
	echo "run_elastic: EIKONAL_MODEL must be one of: single, double, triple"
	return 1 2>/dev/null || exit 1
	;;
esac

MODEL_OVERRIDE="GENERAL.json:PARAM_SOFT.active_model=\"${EIKONAL_MODEL}\""
PPBAR_BEAM_OVERRIDE='SCATTERING.BEAM=["p+","p-"]'

echo "run_elastic: EIKONAL_MODEL=${EIKONAL_MODEL}"

SCAN_DIR="${SCAN_DIR:-${STUDY_REPO_ROOT}/tmp/elastic/${EIKONAL_MODEL}}"
PP_ENERGIES="${PP_ENERGIES:-200,500,7000,13000,60000}"
PPBAR_ENERGIES="${PPBAR_ENERGIES:-62,546,1960}"

if study_generates
then
    mkdir -p "${SCAN_DIR}/pp" "${SCAN_DIR}/ppbar"
    (
        cd "${SCAN_DIR}/ppbar"
        "${STUDY_REPO_ROOT}/bin/xscan" -i "${STUDY_LAUNCHER_DIR}/gencard_el.json" -e "${PPBAR_ENERGIES}" -l true \
            --set "$MODEL_OVERRIDE" --set "$PPBAR_BEAM_OVERRIDE" 2>&1 | tee scan.log
    )
    (
        cd "${SCAN_DIR}/pp"
        "${STUDY_REPO_ROOT}/bin/xscan" -i "${STUDY_LAUNCHER_DIR}/gencard_el.json" -e "${PP_ENERGIES}" -l true \
            --set "$MODEL_OVERRIDE" 2>&1 | tee scan.log
    )
fi

for BEAM in pp ppbar
do
    ENERGIES="$PP_ENERGIES"
    [[ "$BEAM" == ppbar ]] && ENERGIES="$PPBAR_ENERGIES"
    for SPACE in momentum impact
    do
        python -m "tests.physics.studies.elastic.analysis.analyze_${SPACE}" \
            --model "$EIKONAL_MODEL" --scan-dir "$SCAN_DIR" --beam "$BEAM" \
            --sqrts ${ENERGIES//,/ } --fig-dir "${FIG_DIR:-${STUDY_REPO_ROOT}/figs/elastic}/$BEAM"
    done
done

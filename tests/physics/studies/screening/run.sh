#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Simulation K+K- with screening loop off/on
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/screening/run.sh

if study_generates
then

# Generate
LOOPSCREEN=false study_gr -i ./tests/physics/studies/screening/gencard_ALICECUTS.json -p "GP[CON]<F> -> K+ K-" \
-w true -o "continuum" -f "hepmc3"

LOOPSCREEN=true study_gr -i ./tests/physics/studies/screening/gencard_ALICECUTS.json -p "GP[CON]<F> -> K+ K-" \
-w true -o "continuum_screened" -f "hepmc3"

fi

# Analyze

SCREENING_LABEL="bare and screened" study_analyze \
-i "continuum, continuum_screened" \
-g "321, 321" \
-n "2, 2" \
-l "K^{+}K^{-}, K^{+}K^{-} (S^{2})" \
-M "95, 0.0, 3.0" \
-Y "95, -2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-R false \
-t "#sqrt{s} = 13 TeV, |#eta| < 0.9, p_{T} > 0.15 GeV"
#-X 1000

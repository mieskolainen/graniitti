#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Simulation pi+pi- and analysis with ALICE cuts
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/alice_single/run.sh

study_require_opt_in RUN_EXTERNAL_DATA "ALICE data comparison study"
study_require_files "ALICE data comparison study" ALICE7_2pi.csv



if study_generates
then

# Generate
study_gr -i ./tests/physics/studies/alice_single/gencard_ALICECUTS.json -p "GP[RES+CON]<F> -> pi+ pi-" -e 7000 -o "ALICE7_2pi" -f "hepmc3"

fi


# Analyze
study_analyze \
-i "ALICE7_2pi, ALICE7_2pi.csv" \
-g "211, 211" \
-n "2, 2" \
-l "GRANIITTI, ALICE (arb.norm)" \
-t "#sqrt{s} = 7 TeV, #pi^{+}#pi^{-}, |#eta| < 0.9 #wedge p_{T} > 0.15 GeV" \
-M "95, 0.0, 2.5" \
-Y "95, -1.25, 1.25" \
-P "95, 0.0, 2.0" \
-u ub

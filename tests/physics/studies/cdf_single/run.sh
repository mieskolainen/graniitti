#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Simulation pi+pi- and analysis with CDF cuts
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/cdf_single/run.sh



if study_generates
then

# Generate
study_gr -i ./tests/physics/studies/cdf_single/gencard_CDF14_2pi.json

fi

# Analyze
study_analyze \
-i CDF14_2pi \
-g 211 \
-n 2 \
-l 'GRANIITTI #pi^{+}#pi^{-}' \
-t '#sqrt{s} = 1.96 TeV, |#eta| < 1.3, p_{T} > 0.4 GeV, |Y_{x}| < 1.0' \
-M "95, 0.0, 2.5" \
-Y "95,-1.5, 1.5" \
-P "95, 0.0, 1.5" \
-u ub

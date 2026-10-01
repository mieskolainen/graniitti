#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# LUXpdf and inclusive lepton pair production via gamma-gamma fusion
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/lux/run.sh

if study_generates
then

# Generate
study_gr -i ./tests/physics/studies/lux/gencard_LUX.json

fi

# Scale factor 3 x for three lepton flavors (we generate a sample only for one flavor)
SCALE=3

# Analyze
study_analyze \
-i LUX \
-g 11 \
-n 2 \
-l "l^{+}l^{-} / 13 TeV" \
-M "95, 0.0, 5000" \
-Y "95,-1.5, 1.5" \
-P "95, 0.0, 2.0" \
-u fb \
-S $SCALE #-X 1000

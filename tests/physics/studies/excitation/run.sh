#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Simulation of full pi+pi- spectrum with forward proton excitation (0, 1 or 2)
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/excitation/run.sh

if study_generates
then

# Generate
study_gr -i ./tests/physics/studies/excitation/gencard_ALICECUTS.json -p "GP[RES+CON]<F> -> pi+ pi-" \
-s 0 -o "2pi_excite_0" -f "hepmc3"

study_gr -i ./tests/physics/studies/excitation/gencard_ALICECUTS.json -p "GP[RES+CON]<F> -> pi+ pi-" \
-s 1 -o "2pi_excite_1" -f "hepmc3"

study_gr -i ./tests/physics/studies/excitation/gencard_ALICECUTS.json -p "GP[RES+CON]<F> -> pi+ pi-" \
-s 2 -o "2pi_excite_2" -f "hepmc3"

fi

# Analyze

study_analyze \
-i "2pi_excite_0, 2pi_excite_1, 2pi_excite_2" \
-g "211, 211, 211" \
-n "2, 2, 2" \
-l "p + #pi^{+}#pi^{-} + p, X + #pi^{+}#pi^{-} + p, X + #pi^{+}#pi^{-} + Y" \
-M "95, 0.0, 3.0" \
-Y "95,-1.5, 1.5" \
-P "95, 0.0, 3.0" \
-u ub \
-R false \
-t '#sqrt{s} = 13 TeV, |#eta| < 0.9, p_{T} > 0.15 GeV'
#-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

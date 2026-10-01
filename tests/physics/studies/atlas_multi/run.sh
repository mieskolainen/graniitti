#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize 100000
#
# Simulation pi+pi- with ATLAS roman pot cuts
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/atlas_multi/run.sh

if study_generates
then

# Full
study_gr -i ./tests/physics/studies/atlas_multi/gencard_ATLAS.json \
 -p "GP[RES+CON]<F> -> pi+ pi-" -o "ATLAS" -f "hepmc3"

# Forward transverse dot product < 0
study_gr -i ./tests/physics/studies/atlas_multi/gencard_ATLAS_NEG.json \
 -p "GP[RES+CON]<F> -> pi+ pi-" -o "ATLAS_NEG" -f "hepmc3"

# Forward transverse dot product > 0
study_gr -i ./tests/physics/studies/atlas_multi/gencard_ATLAS_POS.json \
 -p "GP[RES+CON]<F> -> pi+ pi-" -o "ATLAS_POS" -f "hepmc3"

fi
# Analyze

study_analyze \
-i "ATLAS, ATLAS_NEG, ATLAS_POS" \
-g "211,211,211" \
-n "2,2,2" \
-l "#pi^{+}#pi^{-}, #pi^{+}#pi^{-} [#delta#phi > #pi/2],  #pi^{+}#pi^{-} [#delta#phi < #pi/2]" \
-M "95, 0.0, 3.0" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-R false \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV, |t| > 0.035 GeV^{2}'

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

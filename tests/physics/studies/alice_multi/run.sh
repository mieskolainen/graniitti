#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Simulation pi+pi-, K+K-, ppbar and analysis with ALICE cuts
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/alice_multi/run.sh

if study_generates
then

# Generate
study_gr -i ./tests/physics/studies/alice_multi/gencard_ALICECUTS.json -p "GP[RES+CON]<F> -> pi+ pi-" -o "ALICE_2pi"   -f "hepmc3"
study_gr -i ./tests/physics/studies/alice_multi/gencard_ALICECUTS.json -p "GP[RES+CON]<F> -> K+ K- @RES{f0_980:1,phi_1020:1,f2_1270:1,f0_1500:1,f2_1525:1,f0_1710:1,f2_1950:1}" -o "ALICE_2K"    -f "hepmc3"
study_gr -i ./tests/physics/studies/alice_multi/gencard_ALICECUTS.json -p "GP[CON]<F> -> p+ p-" -o "ALICE_ppbar" -f "hepmc3"

fi

# Analyze

study_analyze \
-i "ALICE_2pi, ALICE_2K, ALICE_ppbar" \
-g "211, 321, 2212" \
-n "2, 2, 2" \
-l "#pi^{+}#pi^{-}, #it{K}^{+}#it{K}^{-}, #it{p}#bar{#it{p}}" \
-M "95, 0.0, 3.0" \
-Y "95,-1.5, 1.5" \
-P "95, 0.0, 2.0" \
-u ub \
-R false \
-t '#sqrt{s} = 13 TeV, |#eta| < 0.9, p_{T} > 0.15 GeV'

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

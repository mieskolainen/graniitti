#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Simulation and analysis of Durham QCD chic0
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/durham_chic0/run.sh

if study_generates
then

# Include the physical pion decay amplitude and branching fraction

study_gr -i ./tests/physics/studies/durham_chic0/gencard_gg2chic0.json -p "gg[chic(0)]<F> -> pi+ pi-" \
-o "chic0_MSTW2008lo68cl" -f "hepmc3" -q "MSTW2008lo68cl"

study_gr -i ./tests/physics/studies/durham_chic0/gencard_gg2chic0.json -p "gg[chic(0)]<F> -> pi+ pi-" \
-o "chic0_CT10nlo"        -f "hepmc3" -q "CT10nlo"

study_gr -i ./tests/physics/studies/durham_chic0/gencard_gg2chic0.json -p "gg[chic(0)]<F> -> pi+ pi-" \
-o "chic0_MMHT2014lo68cl" -f "hepmc3" -q "MMHT2014lo68cl"

fi

# Analyze

study_analyze \
-i "chic0_MSTW2008lo68cl, chic0_CT10nlo, chic0_MMHT2014lo68cl" \
-g "211, 211, 211" \
-n "2, 2, 2" \
-l "#chi_{c0} (MSTW2008lo68cl), #chi_{c0} (CT10nlo), #chi_{c0} (MMHT2014lo68cl)" \
-M "95, 3.2, 3.6" \
-Y "95,-1.0, 1.0" \
-P "95, 0.0, 3.0" \
-u pb \
-t '#sqrt{s} = 13 TeV, |#eta| < 0.9, p_{T} > 0.15 GeV | #chi_{c0}#rightarrow#pi^{+}#pi^{-}'

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

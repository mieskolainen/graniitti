#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Simulation and analysis of Durham QCD MMbar continuum pi+pi-, K+K-, etaeta, eta'eta'
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/durham_mmbar/run.sh

if study_generates
then

# Generate
study_gr -i ./tests/physics/studies/durham_mmbar/gencard_gg2MMbar.json \
-p "gg[MM]<F> -> pi+ pi-" -o "gg2pipi"     -f "hepmc3"

study_gr -i ./tests/physics/studies/durham_mmbar/gencard_gg2MMbar.json \
-p "gg[MM]<F> -> K+ K-" -o "gg2KK"       -f "hepmc3"

study_gr -i ./tests/physics/studies/durham_mmbar/gencard_gg2MMbar.json \
-p "gg[MM]<F> -> eta0 eta0" -o "gg2etaeta"   -f "hepmc3"

study_gr -i ./tests/physics/studies/durham_mmbar/gencard_gg2MMbar.json \
-p "gg[MM]<F> -> eta'(958)0 eta'(958)0" -o "gg2etapetap" -f "hepmc3"

fi

# Analyze

study_analyze \
-i "gg2pipi, gg2KK, gg2etaeta, gg2etapetap" \
-g "211, 321, 221, 331" \
-n "2, 2, 2, 2" \
-l "#pi^{+}#pi^{-}, #it{K}^{+}#it{K}^{-}, #eta#eta, #eta'#eta'" \
-M "95, 3.0, 15.0" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 6.0" \
-u nb \
-t '#sqrt{s} = 7 TeV, |#eta| < 1.8, p_{T} > 2 GeV (NNPDF31_lo_as_0118)'

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Simulation and analysis with different tensor Pomeron (scalar, pseudoscalar) couplings
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/tensor0/run.sh

if study_generates
then

study_gr -i ./tests/physics/studies/tensor0/gencard_CMSCUTS.json \
-p "TP[RES]<F> -> pi+ pi- @RES{f0_980:1}  @R[f0_980]{g0:1.0, g1:0.0}" -o "f0_0" -f "hepmc3"

study_gr -i ./tests/physics/studies/tensor0/gencard_CMSCUTS.json \
-p "TP[RES]<F> -> pi+ pi- @RES{f0_980:1}  @R[f0_980]{g0:0.0, g1:1.0}" -o "f0_1" -f "hepmc3"

study_gr -i ./tests/physics/studies/tensor0/gencard_CMSCUTS.json \
-p "TP[RES]<F> -> 22 22   @RES{eta:1}     @R[eta]{g0:1.0, g1:0.0}" -o "eta_0" -f "hepmc3"

study_gr -i ./tests/physics/studies/tensor0/gencard_CMSCUTS.json \
-p "TP[RES]<F> -> 22 22   @RES{eta:1}     @R[eta]{g0:0.0, g1:1.0}" -o "eta_1" -f "hepmc3"

study_gr -i ./tests/physics/studies/tensor0/gencard_CMSCUTS.json \
-p "TP[RES]<F> -> pi+ pi- @RES{rho_770:1}" -o "rho" -f "hepmc3"

fi

# Analyze

study_analyze \
-i "f0_0, f0_1, eta_0, eta_1, rho" \
-g "211, 211, 22, 22, 211" \
-n "2, 2, 2, 2, 2" \
-l "f_{0}: <00>, f_{0}: <22>, #eta: <11>, #eta: <33>, #rho" \
-M "95, 0.3, 1.3" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV'

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

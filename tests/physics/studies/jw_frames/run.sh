#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize 100000
#
# Simulation and analysis with different rest frame definitions for the spin polarization density
# 
# Compare the decay angles in each polarization frame
# Diagonal spin densities give flat azimuth only before the finite acceptance cuts
# These plots complement the native complex amplitude covariance tests
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/jw_frames/run.sh

if study_generates
then

study_gr -i ./tests/physics/studies/jw_frames/gencard_FULL13.json -p "MP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:0, JZ1:0, JZ2:1} @MP_FRAME:CM" \
-o "f2_JZ2_CM" -f "hepmc3"

study_gr -i ./tests/physics/studies/jw_frames/gencard_FULL13.json -p "MP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:0, JZ1:0, JZ2:1} @MP_FRAME:HX" \
-o "f2_JZ2_HX" -f "hepmc3"

study_gr -i ./tests/physics/studies/jw_frames/gencard_FULL13.json -p "MP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:0, JZ1:0, JZ2:1} @MP_FRAME:CS" \
-o "f2_JZ2_CS" -f "hepmc3"

fi

# Analyze

study_analyze \
-i "f2_JZ2_CM, f2_JZ2_HX, f2_JZ2_CS" \
-g "211, 211, 211" \
-n "2, 2, 2" \
-l "f_{2}: #lambda=#pm2 (CM), f_{2}: #lambda=#pm2 (HX), f_{2}: #lambda=#pm2 (CS)" \
-M "95, 0.5, 1.6" \
-Y "95,-9.0, 9.0" \
-P "95, 0.0, 2.0" \
-u ub \
-R false \
-t '#sqrt{s} = 13 TeV, |#eta| < 9.0, p_{T} > 0 GeV'

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

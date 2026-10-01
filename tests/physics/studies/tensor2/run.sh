#!/usr/bin/env bash
set -euo pipefail

STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize

#
# Simulation and analysis with different tensor Pomeron spin-2 couplings
#

# Run with: NEVENTS=10000 bash ./tests/physics/studies/tensor2/run.sh

if study_generates; then

  study_gr -i ./tests/physics/studies/tensor2/gencard_CMSCUTS.json \
  -p "TP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{g0:1.0, g1:0.0, g2:0.0, g3:0.0, g4:0.0, g5:0.0, g6:0.0}" -o "f2_0" -f "hepmc3"

  study_gr -i ./tests/physics/studies/tensor2/gencard_CMSCUTS.json \
  -p "TP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{g0:0.0, g1:1.0, g2:0.0, g3:0.0, g4:0.0, g5:0.0, g6:0.0}" -o "f2_1" -f "hepmc3"

  study_gr -i ./tests/physics/studies/tensor2/gencard_CMSCUTS.json \
  -p "TP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{g0:0.0, g1:0.0, g2:1.0, g3:0.0, g4:0.0, g5:0.0, g6:0.0}" -o "f2_2" -f "hepmc3"

  study_gr -i ./tests/physics/studies/tensor2/gencard_CMSCUTS.json \
  -p "TP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{g0:0.0, g1:0.0, g2:0.0, g3:1.0, g4:0.0, g5:0.0, g6:0.0}" -o "f2_3" -f "hepmc3"

  study_gr -i ./tests/physics/studies/tensor2/gencard_CMSCUTS.json \
  -p "TP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{g0:0.0, g1:0.0, g2:0.0, g3:0.0, g4:1.0, g5:0.0, g6:0.0}" -o "f2_4" -f "hepmc3"

  study_gr -i ./tests/physics/studies/tensor2/gencard_CMSCUTS.json \
  -p "TP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{g0:0.0, g1:0.0, g2:0.0, g3:0.0, g4:0.0, g5:1.0, g6:0.0}" -o "f2_5" -f "hepmc3"

  study_gr -i ./tests/physics/studies/tensor2/gencard_CMSCUTS.json \
  -p "TP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{g0:0.0, g1:0.0, g2:0.0, g3:0.0, g4:0.0, g5:0.0, g6:1.0}" -o "f2_6" -f "hepmc3"

fi

# Analyze

study_analyze \
-i "f2_0, f2_1, f2_2, f2_3, f2_4, f2_5, f2_6" \
-g "211, 211, 211, 211, 211, 211, 211" \
-n "2, 2, 2, 2, 2, 2, 2" \
-l "f_{2}: <02>, f_{2}: <20> - <22>, f_{2}: <20> + <22>, f_{2}: <24>, f_{2}: <42>, f_{2}: <44>, f_{2}: <64>" \
-M "95, 0.75, 1.75" \
-Y "95, -2.5, 2.5" \
-P "95,  0.0, 2.0" \
-u ub \
-R false \
-t "#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV"

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

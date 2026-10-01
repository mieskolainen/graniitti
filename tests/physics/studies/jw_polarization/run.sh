#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize 100000
#
# Simulation and analysis with different diagonal spin-density matrix elements (Jacob-Wick amplitudes)
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/jw_polarization/run.sh

if study_generates
then

study_gr -i ./tests/physics/studies/jw_polarization/gencard_CMSCUTS.json -p "MP[RES]<F> -> pi+ pi- @RES{f0_980:1}  " \
-o "f0_980" -f "hepmc3"

study_gr -i ./tests/physics/studies/jw_polarization/gencard_CMSCUTS.json -p "MP[RES]<F> -> pi+ pi- @RES{rho_770:1} @R[rho_770]{JZ0:1, JZ1:0} " \
-o "rho_JZ0" -f "hepmc3"

study_gr -i ./tests/physics/studies/jw_polarization/gencard_CMSCUTS.json -p "MP[RES]<F> -> pi+ pi- @RES{rho_770:1} @R[rho_770]{JZ0:0, JZ1:1} " \
-o "rho_JZ1" -f "hepmc3"

study_gr -i ./tests/physics/studies/jw_polarization/gencard_CMSCUTS.json -p "MP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:1, JZ1:0, JZ2:0} " \
-o "f2_JZ0" -f "hepmc3"

study_gr -i ./tests/physics/studies/jw_polarization/gencard_CMSCUTS.json -p "MP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:0, JZ1:1, JZ2:0} " \
-o "f2_JZ1" -f "hepmc3"

study_gr -i ./tests/physics/studies/jw_polarization/gencard_CMSCUTS.json -p "MP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:0, JZ1:0, JZ2:1} " \
-o "f2_JZ2" -f "hepmc3"

fi

# Analyze

study_analyze \
-i "f0_980, rho_JZ0, rho_JZ1, f2_JZ0, f2_JZ1, f2_JZ2" \
-g "211, 211, 211, 211, 211, 211" \
-n "2, 2, 2, 2, 2, 2" \
-l "f_{0}, #rho: J_{z}=0, #rho: J_{z}=#pm1, f_{2}: J_{z}=0, f_{2}: J_{z}=#pm1, f_{2}: J_{z}=#pm2" \
-M "95, 0.5, 1.6" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-R false \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV'
#-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

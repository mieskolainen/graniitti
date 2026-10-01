#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize 1000000
#
# Simulation Pomeron-Odderon and Gamma-Pomeron -> phi(1020)
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/odderon_photoprod/run.sh

if study_generates
then

# Generate

study_gr -i ./tests/physics/studies/odderon_photoprod/gencard_CMSCUTS.json \
-p "GP[RES]<F> -> K+ K- @RES{phi_1020:1}" -o "photo_phi"        -f "hepmc3"

study_gr -i ./tests/physics/studies/odderon_photoprod/gencard_CMSCUTSPOTS.json \
-p "GP[RES]<F> -> K+ K- @RES{phi_1020:1}" -o "photo_phi_pots"   -f "hepmc3"

study_gr -i ./tests/physics/studies/odderon_photoprod/gencard_CMSCUTS.json \
-p "GP[RES]<F> -> K+ K- @RES{phi_1020_odd:1}" -o "odderon_phi"      -f "hepmc3"

study_gr -i ./tests/physics/studies/odderon_photoprod/gencard_CMSCUTSPOTS.json \
-p "GP[RES]<F> -> K+ K- @RES{phi_1020_odd:1}" -o "odderon_phi_pots" -f "hepmc3"

fi

# Analyze

study_analyze \
-i "photo_phi, photo_phi_pots, odderon_phi, odderon_phi_pots" \
-g "321, 321, 321, 321" \
-n "2, 2, 2, 2" \
-l "#gammaP #rightarrow #phi, #gammaP #rightarrow #phi (RP), OP #rightarrow #phi, OP #rightarrow #phi (RP)"  \
-M "95, 0.99, 1.05" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u nb \
-t '#phi #rightarrow #it{K}^{+}#it{K}^{-} | #sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV (RP: |t_{1,2}| > 0.05 GeV^{2})'

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

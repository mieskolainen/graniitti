#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Simulation of rhorho, phiphi and with the same with tensor pomeron
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/rhorho_phiphi/run.sh

if study_generates
then

study_gr -i ./tests/physics/studies/rhorho_phiphi/gencard_CMSCUTS.json -p "MP[CON]<F> -> rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}" \
-o "rhorho" -f "hepmc3"

study_gr -i ./tests/physics/studies/rhorho_phiphi/gencard_CMSCUTS.json -p "MP[CON]<F> -> phi(1020)0 > {K+ K-} phi(1020)0 > {K+ K-}" \
-o "phiphi" -f "hepmc3"

study_gr -i ./tests/physics/studies/rhorho_phiphi/gencard_CMSCUTS.json -p "TP[CON]<F> -> rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}" \
-o "rhorho_tensor" -f "hepmc3"

study_gr -i ./tests/physics/studies/rhorho_phiphi/gencard_CMSCUTS.json -p "TP[CON]<F> -> phi(1020)0 > {K+ K-} phi(1020)0 > {K+ K-}" \
-o "phiphi_tensor" -f "hepmc3"

fi

# Analyze

study_analyze \
-i "rhorho, phiphi, rhorho_tensor, phiphi_tensor" \
-g "113, 333, 113, 333" \
-n "2, 2, 2, 2" \
-l "#rho^{0}#rho^{0} #rightarrow #pi^{+}#pi^{-}, #phi#phi #rightarrow K^{+}K^{-}, #rho^{0}#rho^{0} #rightarrow #pi^{+}#pi^{-} (TP), #phi#phi #rightarrow K^{+}K^{-} (TP)" \
-M "95, 0.0, 7.0" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV'

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

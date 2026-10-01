#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Simulation pi+pi-, K+K-, ppbar and analysis with Tensor Pomeron
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/tensor_spectrum/run.sh

if study_generates
then

study_gr -i ./tests/physics/studies/tensor_spectrum/gencard_CMSCUTS.json \
-p "TP[RES+CON]<F> -> pi+ pi- @RES{f0_500:1,rho_770:1,f0_980:1,f2_1270:1,f0_1500:1,f2_1525:1,f0_1710:1,f2_1950:1}" \
-o "tensor_pipi" -f "hepmc3"

study_gr -i ./tests/physics/studies/tensor_spectrum/gencard_CMSCUTS.json \
-p "TP[RES+CON]<F> -> K+ K-   @RES{f0_980:1,phi_1020:1,f2_1270:1,f0_1500:1,f2_1525:1,f0_1710:1,f2_1950:1}" \
-o "tensor_KK" -f "hepmc3"

study_gr -i ./tests/physics/studies/tensor_spectrum/gencard_CMSCUTS.json \
-p "TP[CON]<F> -> p+ p-" \
-o "tensor_ppbar" -f "hepmc3"

fi

# Analyze

study_analyze \
-i "tensor_pipi,tensor_KK,tensor_ppbar" \
-g "211, 321, 2212" \
-n "2, 2, 2" \
-l "#pi^{+}#pi^{-} (TP), K^{+}K^{-} (TP), p#bar{p} (TP)" \
-M "95, 0.0, 2.5" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t "#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV"

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

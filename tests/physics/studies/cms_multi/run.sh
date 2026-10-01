#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Simulation pi+pi-, K+K-, ppbar and analysis with CMS cuts
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/cms_multi/run.sh

if study_generates
then

# Generate
study_gr -i tests/physics/studies/cms_multi/gencard_CMSCUTS.json -p "GP[RES+CON]<F> -> pi+ pi- @RES{f0_980:1}" -o "CMS19_2pi"   -f "hepmc3"
study_gr -i tests/physics/studies/cms_multi/gencard_CMSCUTS.json -p "GP[RES+CON]<F> -> K+ K- @RES{f0_980:1}" -o "CMS19_2K"    -f "hepmc3"
study_gr -i tests/physics/studies/cms_multi/gencard_CMSCUTS.json -p "GP[CON]<F> -> p+ p-" -o "CMS19_ppbar" -f "hepmc3"

fi

# Analyze

study_analyze \
-i "CMS19_2pi, CMS19_2K, CMS19_ppbar" \
-g "211, 321, 2212" \
-n "2, 2, 2" \
-l "#pi^{+}#pi^{-}, #it{K}^{+}#it{K}^{-}, #it{p}#bar{#it{p}}" \
-M "95, 0.0, 3.0" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t "#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV"

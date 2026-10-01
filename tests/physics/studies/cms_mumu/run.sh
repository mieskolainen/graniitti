#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize
#
# Simulation yy -> mu+mu- (low mass and high mass domains) and analysis
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/cms_mumu/run.sh



if study_generates
then

# Simulate mu+mu-
study_gr -i ./tests/physics/studies/cms_mumu/gencard_CMS_mumu_lowmass.json \
-p "yy[EPA]<F> -> mu+ mu-" -o "yy_EPA_lo" -f "hepmc3"

study_gr -i ./tests/physics/studies/cms_mumu/gencard_CMS_mumu_lowmass.json \
-p "yy[QED]<F> -> mu+ mu-" -o "yy_QED_lo" -f "hepmc3"

study_gr -i ./tests/physics/studies/cms_mumu/gencard_CMS_mumu_himass.json  \
-p "yy[EPA]<F> -> mu+ mu-" -o "yy_EPA_hi" -f "hepmc3"

study_gr -i ./tests/physics/studies/cms_mumu/gencard_CMS_mumu_himass.json  \
-p "yy[QED]<F> -> mu+ mu-" -o "yy_QED_hi" -f "hepmc3"

fi

# Analyze lowmass
study_analyze \
-i "yy_EPA_lo, yy_QED_lo" \
-g "13, 13" \
-n "2, 2" \
-l "#gamma#gamma #rightarrow #mu^{+}#mu^{-} (kt-EPA), #gamma#gamma #rightarrow #mu^{+}#mu^{-} (QED)" \
-t "#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV, 0.5 < M < 5 GeV" \
-M "95, 0.0, 6.0" \
-Y "95,-3.0, 3.0" \
-P "95, 0.0, 2.0" \
-u nb

# Analyze highmass
study_analyze \
-i "yy_EPA_hi, yy_QED_hi" \
-g "13, 13" \
-n "2, 2" \
-l "#gamma#gamma #rightarrow #mu^{+}#mu^{-} (kt-EPA), #gamma#gamma #rightarrow #mu^{+}#mu^{-} (QED)" \
-t "#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV, M > 5 GeV" \
-M "95, 4.0, 20.0" \
-Y "95,-3.0, 3.0" \
-P "95, 0.0, 2.0" \
-u nb

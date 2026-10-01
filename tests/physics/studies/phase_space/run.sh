#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize 100000
#
# Simulation of 2pi, 4pi, 6pi with different phase space constructions
# With these cuts, <F> and <C> should match one to one!
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/phase_space/run.sh

if study_generates
then

# Generate
study_gr -i ./tests/physics/studies/phase_space/gencard_four_body.json -p "GP[CON]<F> -> pi+ pi-" -o "F2" -f "hepmc3"
study_gr -i ./tests/physics/studies/phase_space/gencard_four_body.json -p "GP[CON]<C> -> pi+ pi-" -o "C2" -f "hepmc3"

study_gr -i ./tests/physics/studies/phase_space/gencard_four_body.json -p "GP[CON]<F> -> pi+ pi- pi+ pi-" -o "F4" -f "hepmc3"
study_gr -i ./tests/physics/studies/phase_space/gencard_four_body.json -p "GP[CON]<C> -> pi+ pi- pi+ pi-" -o "C4" -f "hepmc3"

study_gr -i ./tests/physics/studies/phase_space/gencard_four_body.json -p "GP[CON]<F> -> pi+ pi- pi+ pi- pi+ pi-" -o "F6" -f "hepmc3"
#study_gr -i ./tests/physics/studies/phase_space/gencard_four_body.json -p "GP[CON]<C> -> pi+ pi- pi+ pi- pi+ pi-" -w true -l false -n $NEVENTS -o "C6" -f "hepmc3"


fi

# Analyze

study_analyze \
-i "F2, C2, F4, C4" \
-g "211, 211, 211, 211" \
-n "2, 2, 4, 4" \
-l "2#pi <F>, 2#pi <C>, 4#pi <F>, 4#pi <C>" \
-M "95, 0.0, 4.0" \
-Y "95,-1.5, 1.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 0.9, p_{T} > 0.15 GeV'

study_analyze \
-i "F6" \
-g "211" \
-n "6" \
-l "6#pi <F>" \
-M "95, 0.0, 4.0" \
-Y "95,-1.5, 1.5" \
-P "95, 0.0, 2.0" \
-u nb \
-t '#sqrt{s} = 13 TeV, |#eta| < 0.9, p_{T} > 0.15 GeV'

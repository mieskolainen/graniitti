#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize 6000
#
# Optimal Transport tests
#

if study_generates
then

# Generate
study_gr -i ./tests/physics/studies/optimal_transport/gencard_test.json -p "GP[CON]<F> -> pi+ pi-"     -o continuum1    -f hepmc3 -r 222
study_gr -i ./tests/physics/studies/optimal_transport/gencard_test.json -p "GP[CON]<F> -> pi+ pi-"     -o continuum2    -f hepmc3 -r 333
study_gr -i ./tests/physics/studies/optimal_transport/gencard_test.json -p "GP[RES+CON]<F> -> pi+ pi-" -o fullspectrum1 -f hepmc3 -r 555
study_gr -i ./tests/physics/studies/optimal_transport/gencard_test.json -p "GP[RES+CON]<F> -> pi+ pi-" -o fullspectrum2 -f hepmc3 -r 777
fi

# Analyze
LAMBDA=0.01
ITER=1500

#./bin/ot -i ./output/fullspectrum1.hepmc3,./output/fullspectrum1.hepmc3 -a $LAMBDA -r $ITER
#./bin/ot -i ./output/fullspectrum1.hepmc3,./output/fullspectrum2.hepmc3 -a $LAMBDA -r $ITER

./bin/ot -i ./output/fullspectrum1.hepmc3,./output/continuum1.hepmc3    -a $LAMBDA -r $ITER
#./bin/ot -i ./output/fullspectrum1.hepmc3,./output/continuum2.hepmc3    -a $LAMBDA -r $ITER

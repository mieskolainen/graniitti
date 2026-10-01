#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize 100000
#
# Simulation of kt-epa, DZ-epa and Durham flux only processes (matrix element = 1)
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/flux/run.sh

if study_generates
then

# Generate
study_gr -i ./tests/physics/studies/flux/gencard_fluxcuts.json -p "gg[FLUX]<F> -> 22 22" \
-o "pure_gg_flux" -f "hepmc3"

study_gr -i ./tests/physics/studies/flux/gencard_fluxcuts.json -p "yy[FLUX]<F> -> 22 22" \
-o "pure_yy_flux" -f "hepmc3"

study_gr -i ./tests/physics/studies/flux/gencard_fluxcuts.json -p "yy_DZ[FLUX]<P> -> 22 22" \
-o "pure_yy_dz_flux" -f "hepmc3"

fi

# Analyze

study_analyze \
-i "pure_gg_flux, pure_yy_flux, pure_yy_dz_flux" \
-g "22, 22, 22" \
-n "2, 2, 2" \
-l "Durham flux, kt-EPA flux , DZ EPA flux" \
-M "95, 50.0, 900.0" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 3.0" \
-u nb \
-R false \
-t '|A|^{2} = 1.0, #sqrt{s} = 14 TeV, Y_{X} = [-2.5, 2.5]'
#-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

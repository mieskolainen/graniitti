#!/bin/sh
#
# Elastic scattering at different energies
#
# Run with: source ./tests/run_elastic/run.sh

# Screening loop must be on!

./bin/xscan -i ./tests/run_minbias/el.json -e 62,500,546,1800,7000,8000,13000,60000 -l true

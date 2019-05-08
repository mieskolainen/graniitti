#!/bin/sh
#
# GLOBAL test of the program, output to be inspected by human
#
# Run with: source ./tests/megarun.sh > megarun.out

make -j4 TEST=TRUE


# single process and plotting
yes | source ./tests/run_cdf_single/run.sh
yes | source ./tests/run_alice_single/run.sh

# multiple processes and plotting
yes | source ./tests/run_alice_multi/run.sh
yes | source ./tests/run_cms_multi/run.sh

# harmonic fitting
yes | source ./tests/run_cms_harmonic/run.sh
yes | source ./tests/run_alice_harmonic/run.sh

# eikonal/elastic scattering
yes | source ./tests/run_elastic/run.sh

# luxpdf
yes | source ./tests/run_luxtest/run.sh

# minimum bias
./bin/minbias 7000,13000 1000

# catch2 driven tests
./bin/testbench0
./bin/testbench1
./bin/testbench2

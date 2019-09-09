#!/bin/sh
#
# GLOBAL benchmark test of the program, output Quality Assurance done by human
#
# Run with: source ./tests/megarun.sh > megarun.out

make -j4 TEST=TRUE ROOT14=TRUE

# For different ROOT installation versions
#make -j4 TEST=TRUE

read -p "megarun: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then
CMD=y
else
CMD=n
fi

# single process and plotting
yes $CMD | source ./tests/run_screening/run.sh
yes $CMD | source ./tests/run_cdf_single/run.sh
yes $CMD | source ./tests/run_alice_single/run.sh

# multiple processes and plotting
yes $CMD | source ./tests/run_alice_multi/run.sh
yes $CMD | source ./tests/run_cms_multi/run.sh

# spherical harmonic expansion
yes $CMD | source ./tests/run_cms_harmonic/run.sh
yes $CMD | source ./tests/run_alice_harmonic/run.sh

# eikonal/elastic scattering
yes $CMD | source ./tests/run_elastic/run.sh

# kt-EPA vs full QED yy in low-mass and high-mass domain
yes $CMD | source ./tests/run_cms_mumu_low_high/run.sh

# luxpdf gamma-gamma
yes $CMD | source ./tests/run_luxtest/run.sh

# tensor Pomeron
yes $CMD | source ./tests/run_tensor0_multi/run.sh
yes $CMD | source ./tests/run_tensor2_multi/run.sh
yes $CMD | source ./tests/run_tensor_spectrum/run.sh

# minimum bias
./bin/minbias "7000,13000" 1000

# catch2 driven tests
./bin/testbench0
./bin/testbench1
./bin/testbench2
./bin/testbench3

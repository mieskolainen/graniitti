#!/bin/sh
#
# GLOBAL benchmark test of the program, output Quality Assurance done by human
#
# Run with: source ./tests/megarun.sh > megarun.out

make -j4 TEST=TRUE ROOT14=TRUE

# For different ROOT installation versions
#make -j4 TEST=TRUE

read -p "plotrun: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then
CMD=y
else
CMD=n
fi

# Not set
#if [ -z "$EVENTS" ]
#then
EVENTS=1000000
#fi

# Processes and plotting
yes $CMD | source ./tests/run_screening/run.sh
yes $CMD | source ./tests/run_alice_single/run.sh
yes $CMD | source ./tests/run_cms_multi/run.sh
yes $CMD | source ./tests/run_excitation/run.sh

# Spin density generation and control
yes $CMD | source ./tests/run_JW_polarizations/run.sh
yes $CMD | source ./tests/run_JW_frames/run.sh

# Tensor Pomeron
yes $CMD | source ./tests/run_tensor0_multi/run.sh
yes $CMD | source ./tests/run_tensor2_multi/run.sh
yes $CMD | source ./tests/run_tensor_spectrum/run.sh

# Spherical harmonic expansion
yes $CMD | source ./tests/run_cms_harmonic/run.sh
yes $CMD | source ./tests/run_alice_harmonic/run.sh

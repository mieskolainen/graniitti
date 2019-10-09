#!/bin/sh
#
# Generate simulations and analysis plots
#
# Run with: time xxx.sh
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
#fi

EVENTS=1000000
yes $CMD | source ./tests/run_atlas_multi/run.sh


EVENTS=2000000
yes $CMD | source ./tests/run_excitation/run.sh


EVENTS=2000000
yes $CMD | source ./tests/run_screening/run.sh
yes $CMD | source ./tests/run_alice_multi/run.sh
yes $CMD | source ./tests/run_JW_polarization/run.sh
yes $CMD | source ./tests/run_JW_frames/run.sh


# Tensor Pomeron
EVENTS=1000000
#yes $CMD | source ./tests/run_tensor0_multi/run.sh
yes $CMD | source ./tests/run_tensor2_multi/run.sh
#yes $CMD | source ./tests/run_tensor_spectrum/run.sh


# Spherical harmonic expansion (slowest)
EVENTS=1000000
yes $CMD | source ./tests/run_cms_harmonic/run.sh
#yes $CMD | source ./tests/run_alice_harmonic/run.sh


echo "plotrun.sh [Done]"

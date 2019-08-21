#!/bin/sh
#
# Simulation and analysis with different tensor Pomeron spin-0 couplings
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
source ./tests/run_tensor0_multi/combogr.sh

fi

# Analyze
source ./tests/run_tensor0_multi/comboanalyze.sh

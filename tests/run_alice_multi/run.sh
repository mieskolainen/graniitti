#!/bin/sh
#
# Simulation and analysis with ALICE cuts
#
# Run with: source ./tests/run_xxx/run.sh

read -p "alicerun: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
source ./tests/run_alice_multi/combogr.sh

fi

# Analyze
source ./tests/run_alice_multi/comboanalyze.sh

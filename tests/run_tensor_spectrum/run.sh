#!/bin/sh
#
# Simulation and analysis with Tensor Pomeron
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
source ./tests/run_tensor_spectrum/combogr.sh

fi

# Analyze
source ./tests/run_tensor_spectrum/comboanalyze.sh

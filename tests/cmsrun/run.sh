#!/bin/sh
#
# Run with: source ./tests/alicerun/run.sh

read -p "cmsrun: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
source ./tests/cmsrun/combogr.sh

fi

# Analyze
source ./tests/cmsrun/comboanalyze.sh

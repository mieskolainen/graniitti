#!/bin/sh
#
# Run with: source ./tests/alicerun/run.sh

read -p "alicerun: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
source ./tests/alicerun/combogr.sh

fi

# Analyze
source ./tests/alicerun/comboanalyze.sh

#!/bin/sh
#
# Simulation and analysis with CMS cuts
#
# Run with: source ./tests/run_xxx/run.sh

read -p "cmsrun: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
source ./tests/run_cms_multi/combogr.sh

fi

# Analyze
source ./tests/run_cms_multi/comboanalyze.sh

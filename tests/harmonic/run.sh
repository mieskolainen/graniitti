#!/bin/sh
#
# Run with: source ./tests/harmonic/run.sh

read -p "Harmonic expansion test: Generate events (or only analyze)? [y/n]" -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
./bin/gr -i ./tests/harmonic/SH_2pi_REF.json -n 1000
./bin/gr -i ./tests/harmonic/SH_2pi.json     -n 1000

fi

# Analyze
./bin/fitharmonic -r SH_2pi_REF -i SH_2pi -l 4 -f HE -m 40,0.4,1.5

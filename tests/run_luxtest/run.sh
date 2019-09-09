#!/bin/sh
#
#
# LUXpdf and inclusive lepton pair production via gamma-gamma fusion
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
./bin/gr -i ./tests/run_luxtest/LUX.json -w true -n 50000

fi

# Scale factor 3 x for three lepton flavors (we generate a sample only for one flavor)
SCALE=3

# Analyze
./bin/analyze \
-i LUX \
-g 11 \
-n 2 \
-l "l^{+}l^{-} / 13 TeV" \
-M "95, 0.0, 5000" \
-Y "95,-1.5, 1.5" \
-P "95, 0.0, 2.0" \
-u fb \
-S $SCALE #-X 1000

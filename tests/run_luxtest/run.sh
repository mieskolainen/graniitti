#!/bin/sh
#
#
# LUXpdf and inclusive lepton pair production via gamma-gamma fusion
#
# Run with: source ./tests/run_xxx/run.sh

# Generate
./bin/gr -i ./tests/run_luxtest/LUX.json -w true -n 50000

# Scale factor 3 x for three lepton flavors (we generate a sample only for one flavor)
SCALE=3

# Analyze
./bin/analyze -i LUX \
-g 11 \
-n 2 \
-l 'l^{+}l^{-} / 13 TeV ' \
-M 5000 -Y 1.5 -P 2.0 \
-u fb \
-S $SCALE #-X 1000

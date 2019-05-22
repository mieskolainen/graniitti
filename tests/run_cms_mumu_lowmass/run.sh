#!/bin/sh
#
# Simulation and analysis of low mass mumu
#
# Run with: source ./tests/run_xxx/run.sh

POMLOOP=false

# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=1.0

read -p "Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Simulate
make -j4 && ./bin/gr -i ./tests/processes/CMS_mumu_lowmass.json -n 1000000 -w true

fi

# Analyze
make -j4 && ./bin/analyze \
-i CMS_mumu_lowmass \
-g 13 \
-n 2 \
-l '#gamma#gamma #rightarrow #mu^{+}#mu^{-}' \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.1 GeV' \
-M 95,0.0,3.0 \
-Y 95,-3.0,3.0 \
-P 95,0.0,2.0 \
-u nb \
-S $S2


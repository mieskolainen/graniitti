#!/bin/sh
#
# Simulation and analysis with ALICE cuts
#
# Run with: source ./tests/run_xxx/run.sh

POMLOOP=false

# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
#S2=1.0
S2=0.18

read -p "Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
make -j4 && ./bin/gr -i ./tests/processes/ALICE7_2pi.json -n 100000 -l $POMLOOP -w true

fi

# Analyze
make -j4 && ./bin/analyze \
-i ALICE7_2pi,ALICE7_2pi.csv \
-g 211,211 \
-n 2,2 \
-l 'GRANIITTI','ALICE (arb.norm)' \
-t '#pi^{+}#pi^{-}, |#eta| < 0.9 #wedge p_{T} > 0.1 GeV, #sqrt{s} = 7 TeV' \
-M 2.5 -Y 1.25 -P 2.0 \
-u ub \
-S $S2,1.0

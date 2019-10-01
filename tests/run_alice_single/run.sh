#!/bin/sh
#
# Simulation pi+pi- and analysis with ALICE cuts
#
# Run with: source ./tests/run_xxx/run.sh

POMLOOP=false

# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.15

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Not set
if [ -z "$EVENTS" ]
then
	EVENTS=10000
fi

# Generate
./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[RES+CON]<F> -> pi+ pi-" -e 7000 -n $EVENTS -l $POMLOOP -w true -o "ALICE7_2pi" -f "hepmc3"

fi


# Analyze
./bin/analyze \
-i "ALICE7_2pi, ALICE7_2pi.csv" \
-g "211, 211" \
-n "2, 2" \
-l "GRANIITTI, ALICE (arb.norm)" \
-t "#sqrt{s} = 7 TeV, #pi^{+}#pi^{-}, |#eta| < 0.9 #wedge p_{T} > 0.15 GeV" \
-M "95, 0.0, 2.5" \
-Y "95, -1.25, 1.25" \
-P "95, 0.0, 2.0" \
-u ub \
-S "$S2, 1.0"

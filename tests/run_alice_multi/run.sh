#!/bin/sh
#
# Simulation pi+pi-, K+K-, ppbar and analysis with ALICE cuts
#
# Run with: source ./tests/run_xxx/run.sh

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
./bin/gr -i ./tests/processes/ALICE_2pi.json -p "PP[RES+CON]<F> -> pi+ pi-" -w true -l false -n $EVENTS
./bin/gr -i ./tests/processes/ALICE_2pi.json -p "PP[RES+CON]<F> -> K+ K-"   -w true -l false -n $EVENTS
./bin/gr -i ./tests/processes/ALICE_2pi.json -p "PP[RES+CON]<F> -> p+ p-"   -w true -l false -n $EVENTS

fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.15

./bin/analyze \
-i "ALICE_2pi, ALICE_2K, ALICE_ppbar" \
-g "211, 321, 2212" \
-n "2, 2, 2" \
-l "#pi^{+}#pi^{-}, #it{K}^{+}#it{K}^{-}, #it{p}#bar{#it{p}}" \
-M "95, 0.0, 3.0" \
-Y "95,-1.5, 1.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 0.9, p_{T} > 0.15 GeV' \
-S "$S2, $S2, $S2" #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

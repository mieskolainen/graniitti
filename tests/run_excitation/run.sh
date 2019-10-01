#!/bin/sh
#
# Simulation with forward proton excitation (0, 1 or 2)
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
./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[CON]<F> -> pi+ pi-" \
-w true -l false -n $EVENTS -s 0 -o "2pi_excite_0" -f "hepmc3"

./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[CON]<F> -> pi+ pi-" \
-w true -l false -n $EVENTS -s 1 -o "2pi_excite_1" -f "hepmc3"

./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[CON]<F> -> pi+ pi-" \
-w true -l false -n $EVENTS -s 2 -o "2pi_excite_2" -f "hepmc3"

fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.15

./bin/analyze \
-i "2pi_excite_0, 2pi_excite_1, 2pi_excite_2" \
-g "211, 211, 211" \
-n "2, 2, 2" \
-l "p + #pi^{+}#pi^{-} + p, X + #pi^{+}#pi^{-} + p, X + #pi^{+}#pi^{-} + Y" \
-M "95, 0.0, 3.0" \
-Y "95,-1.5, 1.5" \
-P "95, 0.0, 3.0" \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 0.9, p_{T} > 0.15 GeV' \
-S "$S2, $S2, $S2" -X 100000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

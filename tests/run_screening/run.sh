#!/bin/sh
#
# Simulation K+K- with screening loop off/on
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
./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[CON]<F> -> K+ K-" \
-p "PP[CON]<F> -> K+ K-" -w true -l false -n $EVENTS -o "continuum" -f "hepmc3"

./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[CON]<F> -> K+ K-" \
-p "PP[CON]<F> -> K+ K-" -w true -l true  -n $EVENTS -o "continuum_screened" -f "hepmc3"

fi

# Analyze

./bin/analyze \
-i "continuum, continuum_screened" \
-g "321, 321" \
-n "2, 2" \
-l "K^{+}K^{-}, K^{+}K^{-} (S^{2})" \
-M "95, 0.0, 3.0" \
-Y "95, -2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t "#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV"

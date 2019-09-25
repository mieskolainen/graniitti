#!/bin/sh
#
# Simulation of 2pi, 4pi, 6pi with different phase space constructions.
# With these cuts, <F> and <C> should match one to one!
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Not set
if [ -z "$EVENTS" ]
then
	EVENTS=100000
fi

# Generate
./bin/gr -i ./tests/run_compare_PS/card.json -p "PP[CON]<F> -> pi+ pi-" -w true -l false -n $EVENTS -o "F2" -f "hepmc3"
./bin/gr -i ./tests/run_compare_PS/card.json -p "PP[CON]<C> -> pi+ pi-" -w true -l false -n $EVENTS -o "C2" -f "hepmc3"

./bin/gr -i ./tests/run_compare_PS/card.json -p "PP[CON]<F> -> pi+ pi- pi+ pi-" -w true -l false -n $EVENTS -o "F4" -f "hepmc3"
./bin/gr -i ./tests/run_compare_PS/card.json -p "PP[CON]<C> -> pi+ pi- pi+ pi-" -w true -l false -n $EVENTS -o "C4" -f "hepmc3"

./bin/gr -i ./tests/run_compare_PS/card.json -p "PP[CON]<F> -> pi+ pi- pi+ pi- pi+ pi-" -w true -l false -n $EVENTS -o "F6" -f "hepmc3"
#./bin/gr -i ./tests/run_compare_PS/card.json -p "PP[CON]<C> -> pi+ pi- pi+ pi- pi+ pi-" -w true -l false -n $EVENTS -o "C6" -f "hepmc3"


fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.15

./bin/analyze \
-i "F2, C2, F4, C4, F6" \
-g "211, 211, 211, 211, 211" \
-n "2, 2, 4, 4, 6" \
-l "2#pi <F>, 2#pi <C>, 4#pi <F>, 4#pi <C>, 6#pi <F>" \
-M "95, 0.0, 4.0" \
-Y "95,-1.5, 1.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 0.9, p_{T} > 0.1 GeV' \
-S "$S2, $S2, $S2, $S2, $S2" #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

#!/bin/sh
#
# Simulation pi+pi- with ATLAS roman pot cuts
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

# Full
./bin/gr -i ./tests/run_atlas_multi/ATLAS.json \
	-w true -l false -n $EVENTS -p "PP[RES+CON]<F> -> pi+ pi-" -o "ATLAS"     -f "hepmc3"

# DeltaPhi < 0
./bin/gr -i ./tests/run_atlas_multi/ATLAS_NEG.json \
	-w true -l false -n $EVENTS -p "PP[RES+CON]<F> -> pi+ pi-" -o "ATLAS_NEG" -f "hepmc3"

# DeltaPhi > 0
./bin/gr -i ./tests/run_atlas_multi/ATLAS_POS.json \
	-w true -l false -n $EVENTS -p "PP[RES+CON]<F> -> pi+ pi-" -o "ATLAS_POS" -f "hepmc3"

fi
# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.15

./bin/analyze \
-i "ATLAS, ATLAS_NEG, ATLAS_POS" \
-g "211,211,211" \
-n "2,2,2" \
-l "#pi^{+}#pi^{-}, #pi^{+}#pi^{-} [#delta#phi < #pi/2],  #pi^{+}#pi^{-} [#delta#phi > #pi/2]" \
-M "95, 0.0, 3.0" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.2 GeV, |t| > 0.03 GeV^{2}' \
-S "$S2, $S2, $S2" #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

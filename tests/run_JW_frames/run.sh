#!/bin/sh
#
# Simulation and analysis with different rest frame definitions for the spin polarization density
#
# Basic check: sample generated with CM should look as it should in CM frame plots,
#              same for HX in HX frame, and GJ in GJ frame (all identical)
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

./bin/gr -i ./tests/processes/CMS13.json -p "PP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:0, JZ1:1, JZ2:0} @FRAME:CM" \
-w true -l false -n $EVENTS -o "f2_JZ0_CM" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json -p "PP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:0, JZ1:1, JZ2:0} @FRAME:HX" \
-w true -l false -n $EVENTS -o "f2_JZ1_HX" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json -p "PP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:0, JZ1:1, JZ2:0} @FRAME:GJ" \
-w true -l false -n $EVENTS -o "f2_JZ2_GJ" -f "hepmc3"

fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.15

./bin/analyze \
-i "f2_JZ0_CM, f2_JZ1_HX, f2_JZ2_GJ" \
-g "211, 211, 211" \
-n "2, 2, 2" \
-l "f_{2}: #lambda=#pm1 (CM), f_{2}: #lambda=#pm1 (HX), f_{2}: #lambda=#pm1 (GJ)" \
-M "95, 0.5, 1.6" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV' \
-S "$S2, $S2, $S2" #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

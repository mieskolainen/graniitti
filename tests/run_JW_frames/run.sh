#!/bin/sh
#
# Simulation and analysis with different rest frame definitions for the spin polarization density
# Basic check: sample generated with SR should look as it should in SR frame plots, same for HE in HE frame, GJ in GJ frame
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
N=100000

./bin/gr -i ./tests/processes/CMS13.json -p "PP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:0, JZ1:1, JZ2:0} @FRAME:SR" \
-w true -l false -n $N -o "f2_JZ0_SR" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json -p "PP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:0, JZ1:1, JZ2:0} @FRAME:HE" \
-w true -l false -n $N -o "f2_JZ1_HE" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json -p "PP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:0, JZ1:1, JZ2:0} @FRAME:GJ" \
-w true -l false -n $N -o "f2_JZ2_GJ" -f "hepmc3"

fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.15

./bin/analyze \
-i "f2_JZ0_SR, f2_JZ1_HE, f2_JZ2_GJ" \
-g "211, 211, 211" \
-n "2, 2, 2" \
-l "f_{2}: #lambda=#pm1 (SR), f_{2}: #lambda=#pm1 (HE), f_{2}: #lambda=#pm1 (GJ)" \
-M "95, 0.5, 1.6" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV' \
-S "$S2, $S2, $S2" #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

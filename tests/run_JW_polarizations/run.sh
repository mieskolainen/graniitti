#!/bin/sh
#
# Simulation and analysis with different tensor Pomeron spin-0 couplings
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
N=100000

./bin/gr -i ./tests/processes/CMS13.json -p "PP[RES]<F> -> pi+ pi- @RES{f0_980:1}  @FRAME:SR" \
-w true -l false -n $N -o "f0_980" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json -p "yP[RES]<F> -> pi+ pi- @RES{rho_770:1} @R[rho_770]{JZ0:1, JZ1:0} @FRAME:SR" \
-w true -l false -n $N -o "rho_JZ0" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json -p "yP[RES]<F> -> pi+ pi- @RES{rho_770:1} @R[rho_770]{JZ0:0, JZ1:1} @FRAME:SR" \
-w true -l false -n $N -o "rho_JZ1" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json -p "PP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:1, JZ1:0, JZ2:0} @FRAME:SR" \
-w true -l false -n $N -o "f2_JZ0" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json -p "PP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:0, JZ1:1, JZ2:0} @FRAME:SR" \
-w true -l false -n $N -o "f2_JZ1" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json -p "PP[RES]<F> -> pi+ pi- @RES{f2_1270:1} @R[f2_1270]{JZ0:0, JZ1:0, JZ2:1} @FRAME:SR" \
-w true -l false -n $N -o "f2_JZ2" -f "hepmc3"

fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.18
S2_photo=0.7

./bin/analyze \
-i "f0_980, rho_JZ0, rho_JZ1, f2_JZ0, f2_JZ1, f2_JZ2" \
-g "211, 211, 211, 211, 211, 211" \
-n "2, 2, 2, 2, 2, 2" \
-l "f_{0}, #rho: #lambda=0, #rho: #lambda=#pm1, f_{2}: #lambda=0, f_{2}: #lambda=#pm1, f_{2}: #lambda=#pm2" \
-M "95, 0.5, 1.6" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV' \
-S "$S2, $S2_photo, $S2_photo, $S2, $S2, $S2" #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

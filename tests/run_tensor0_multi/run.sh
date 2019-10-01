#!/bin/sh
#
# Simulation and analysis with different tensor Pomeron (scalar, pseudoscalar) couplings
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

./bin/gr -i ./tests/processes/CMSCUTS.json \
-p "PP[RESTENSOR]<F> -> pi+ pi- @RES{f0_980:1}  @R[f0_980]{g0:1.0, g1:0.0}"  -w true -l false -n $EVENTS -o "f0_0" -f "hepmc3"

./bin/gr -i ./tests/processes/CMSCUTS.json \
-p "PP[RESTENSOR]<F> -> pi+ pi- @RES{f0_980:1}  @R[f0_980]{g0:0.0, g1:1.0}"  -w true -l false -n $EVENTS -o "f0_1" -f "hepmc3"

./bin/gr -i ./tests/processes/CMSCUTS.json \
-p "PP[RESTENSOR]<F> -> 22 22   @RES{eta:1}     @R[eta]{g0:1.0, g1:0.0}"     -w true -l false -n $EVENTS -o "eta_0" -f "hepmc3"

./bin/gr -i ./tests/processes/CMSCUTS.json \
-p "PP[RESTENSOR]<F> -> 22 22-  @RES{eta:1}     @R[eta]{g0:0.0, g1:1.0}"     -w true -l false -n $EVENTS -o "eta_1" -f "hepmc3"

./bin/gr -i ./tests/processes/CMSCUTS.json \
-p "yP[RESTENSOR]<F> -> pi+ pi- @RES{rho_770:1}" -w true -l false -n $EVENTS -o "rho" -f "hepmc3"

fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.15
S2photo=0.7

./bin/analyze \
-i "f0_0, f0_1, eta_0, eta_1, rho" \
-g "211, 211, 22, 22, 211" \
-n "2, 2, 2, 2, 2" \
-l "f_{0}: <00>, f_{0}: <22>, #eta: <11>, #eta: <33>, #rho" \
-M "95, 0.3, 1.3" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV' \
-S "$S2, $S2, $S2, $S2, $S2photo" #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg


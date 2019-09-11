#!/bin/sh
#
# Simulation Odderon-Pomeron and Gamma-Pomeron -> phi(1020)
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Not set
if [ -z "$EVENTS" ]
then
	EVENTS=1000000
fi

# Generate

./bin/gr -i ./tests/processes/CMS13.json \
-p "yP[RES]<F> -> K+ K- @RES{phi_1020:1}" -w true -l false -n $EVENTS -o "photo_phi"        -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13POTS.json \
-p "yP[RES]<F> -> K+ K- @RES{phi_1020:1}" -w true -l false -n $EVENTS -o "photo_phi_pots"   -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json \
-p "OP[RES]<F> -> K+ K- @RES{phi_1020:1}" -w true -l false -n $EVENTS -o "odderon_phi"      -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13POTS.json \
-p "OP[RES]<F> -> K+ K- @RES{phi_1020:1}" -w true -l false -n $EVENTS -o "odderon_phi_pots" -f "hepmc3"

fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.15
S2_photo=0.7

./bin/analyze \
-i "photo_phi, photo_phi_pots, odderon_phi, odderon_phi_pots" \
-g "321, 321, 321, 321" \
-n "2, 2, 2, 2" \
-l "#gammaP #rightarrow #phi, #gammaP #rightarrow #phi (RP), OP #rightarrow #phi, OP #rightarrow #phi (RP)"  \
-M "95, 0.99, 1.05" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u nb \
-t '#phi #rightarrow #it{K}^{+}#it{K}^{-} | #sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV (RP: |t_{1,2}| > 0.05 GeV^{2})' \
-S "$S2_photo, $S2_photo, $S2, $S2" #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

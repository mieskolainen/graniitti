#!/bin/sh
#
# Simulation and analysis with ALICE cuts
#
# Run with: source ./tests/run_xxx/run.sh

read -p "alicerun: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

N=100000

# Generate
./bin/gr -i ./tests/processes/ALICE_2pi.json -w true -l false -n $N
./bin/gr -i ./tests/processes/ALICE_2K.json -w true -l false -n $N
./bin/gr -i ./tests/processes/ALICE_ppbar.json -w true -l false -n $N

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

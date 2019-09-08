#!/bin/sh
#
# Simulation with screening loop off/on
#
# Run with: source ./tests/run_xxx/run.sh

read -p "cmsrun: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

N=100000

# Generate
./bin/gr -i ./tests/processes/ALICE_2K.json -p "PP[CON]<F> -> K+ K-" -w true -l false -n $N -o "continuum" -f "hepmc3"
./bin/gr -i ./tests/processes/ALICE_2K.json -p "PP[CON]<F> -> K+ K-" -w true -l true  -n $N -o "continuum_screened" -f "hepmc3"

fi

# Analyze

./bin/analyze \
-i "continuum, continuum_screened" \
-g "321, 321" \
-n "2, 2" \
-l "K^{+}K^{-}, K^{+}K^{-} (SCREENED)" \
-M "95, 0.0, 3.0" \
-Y "95, -2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t "#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV"

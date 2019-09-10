#!/bin/sh
#
# Simulation and analysis of Durham QCD MMbar continuum pi+pi-, K+K-, etaeta, eta'eta'
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

N=1000000

# Generate
./bin/gr -i ./tests/processes/gg2MMbar.json -p "gg[CON]<F> -> pi+ pi-"               -w true -l false -n $N -o "gg2pipi"     -f "hepmc3"
./bin/gr -i ./tests/processes/gg2MMbar.json -p "gg[CON]<F> -> K+ K-"                 -w true -l false -n $N -o "gg2KK"       -f "hepmc3"
./bin/gr -i ./tests/processes/gg2MMbar.json -p "gg[CON]<F> -> eta0 eta0"             -w true -l false -n $N -o "gg2etaeta"   -f "hepmc3"
./bin/gr -i ./tests/processes/gg2MMbar.json -p "gg[CON]<F> -> eta'(958)0 eta'(958)0" -w true -l false -n $N -o "gg2etapetap" -f "hepmc3"

fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.05

./bin/analyze \
-i "gg2pipi, gg2KK, gg2etaeta, gg2etapetap" \
-g "211, 321, 221, 331" \
-n "2, 2, 2, 2" \
-l "#pi^{+}#pi^{-}, #it{K}^{+}#it{K}^{-}, #eta#eta, #eta'#eta'" \
-M "95, 3.0, 15.0" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 6.0" \
-u nb \
-t '#sqrt{s} = 7 TeV, |#eta| < 1.8, E_{T} > 2 GeV (MSTW2008lo68cl)' \
-S "$S2, $S2, $S2, $S2" #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

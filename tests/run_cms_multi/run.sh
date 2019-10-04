#!/bin/sh
#
# Simulation pi+pi-, K+K-, ppbar and analysis with CMS cuts
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

# Generate
./bin/gr -i ./tests/processes/CMSCUTS.json -p "PP[RES+CON]<F> -> pi+ pi-" -w true -l false -n $EVENTS -o "CMS19_2pi"   -f "hepmc3"
./bin/gr -i ./tests/processes/CMSCUTS.json -p "PP[RES+CON]<F> -> K+ K-"   -w true -l false -n $EVENTS -o "CMS19_2K"    -f "hepmc3"
./bin/gr -i ./tests/processes/CMSCUTS.json -p "PP[CON]<F> -> p+ p-"   -w true -l false -n $EVENTS -o "CMS19_ppbar" -f "hepmc3"

fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.15

./bin/analyze \
-i "CMS19_2pi, CMS19_2K, CMS19_ppbar" \
-g "211, 321, 2212" \
-n "2, 2, 2" \
-l "#pi^{+}#pi^{-}, #it{K}^{+}#it{K}^{-}, #it{p}#bar{#it{p}}" \
-M "95, 0.0, 3.0" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t "#sqrt{s} = 13 TeV, |#eta| < 2.4, p_{T} > 0.2 GeV" \
-S "$S2, $S2, $S2" #-X 1000

#!/bin/sh

# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.18

read -p "Generate events (or only analyze)? [y/n]" -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
make -j4 && ./bin/gr -i ./tests/processes/ALICE7_2pi.json -n 100000 -l false -w true

fi

# Analyze
make -j4 ROOT=TRUE && ./bin/analyze -i ALICE7_2pi,ALICE7_2pi.csv -g 211,211 -n 2,2 -l \
	'GRANIITTI','ALICE (arb.norm)' -t '#pi^{+}#pi^{-}, |#eta| < 0.9, p_{T} > 0.1 GeV, #sqrt{s} = 7 TeV' \
	-M 2.5 -Y 1.25 -P 0.4 -u ub -S $S2,1.0

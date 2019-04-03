#!/bin/sh

# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.2

read -p "Generate events (or only analyze)? [y/n]" -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
make -j4 && ./bin/gr -i ./tests/processes/CDF14_2pi.json -n 300000 -l false -w true

fi

# Analyze
make -j4 ROOT=TRUE && ./bin/analyze -i CDF14_2pi -g 211 -n 2 -l \
	'#pi^{+}#pi^{-}' -t '#sqrt{s} = 1.96 TeV, |#eta| < 1.3, p_{T} > 0.4 GeV, |Y_{X}| < 1.0' \
	-M 2.5 -Y 1.5 -P 1.5 -u ub -S $S2

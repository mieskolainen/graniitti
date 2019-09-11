#!/bin/sh
#
# Simulation yy -> mu+mu- (low mass and high mass domains) and analysis
#
# Run with: source ./tests/run_xxx/run.sh

POMLOOP=false

# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=1.0

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Not set
if [ -z "$EVENTS" ]
then
	EVENTS=10000
fi

# Simulate mu+mu-
./bin/gr -i ./tests/run_cms_mumu_low_high/CMS_mumu_lowmass.json \
-p "yy[CON]<F> -> mu+ mu-" -n $EVENTS -w true -o "yy_EPA_lo" -f "hepmc3"

./bin/gr -i ./tests/run_cms_mumu_low_high/CMS_mumu_lowmass.json \
-p "yy[QED]<F> -> mu+ mu-" -n $EVENTS -w true -o "yy_QED_lo" -f "hepmc3"

./bin/gr -i ./tests/run_cms_mumu_low_high/CMS_mumu_himass.json  \
-p "yy[CON]<F> -> mu+ mu-" -n $EVENTS -w true -o "yy_EPA_hi" -f "hepmc3"

./bin/gr -i ./tests/run_cms_mumu_low_high/CMS_mumu_himass.json  \
-p "yy[QED]<F> -> mu+ mu-" -n $EVENTS -w true -o "yy_QED_hi" -f "hepmc3"

fi

# Analyze lowmass
./bin/analyze \
-i "yy_EPA_lo, yy_QED_lo" \
-g "13, 13" \
-n "2, 2" \
-l "#gamma#gamma #rightarrow #mu^{+}#mu^{-} (kt-EPA), #gamma#gamma #rightarrow #mu^{+}#mu^{-} (QED)" \
-t "#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.1 GeV, 0.5 < M < 5 GeV" \
-M "95, 0.0, 6.0" \
-Y "95,-3.0, 3.0" \
-P "95, 0.0, 2.0" \
-u nb \
-S "$S2, $S2"

# Analyze highmass
./bin/analyze \
-i "yy_EPA_hi, yy_QED_hi" \
-g "13, 13" \
-n "2, 2" \
-l "#gamma#gamma #rightarrow #mu^{+}#mu^{-} (kt-EPA), #gamma#gamma #rightarrow #mu^{+}#mu^{-} (QED)" \
-t "#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.1 GeV, M > 5 GeV" \
-M "95, 4.0, 20.0" \
-Y "95,-3.0, 3.0" \
-P "95, 0.0, 2.0" \
-u nb \
-S "$S2, $S2"

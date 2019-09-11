#!/bin/sh
#
# Simulation and analysis of Durham QCD chic0
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

# Generate, operator &> decouples the decay phase space integral volume, thus we can apply branching ratio manually

./bin/gr -i ./tests/processes/gg2chic0.json -p "gg[chic(0)]<F> &> pi+ pi-" -w true -l false -n $EVENTS \
-o "chic0_MSTW2008lo68cl" -f "hepmc3" -q "MSTW2008lo68cl"

./bin/gr -i ./tests/processes/gg2chic0.json -p "gg[chic(0)]<F> &> pi+ pi-" -w true -l false -n $EVENTS \
-o "chic0_CT10nlo"        -f "hepmc3" -q "CT10nlo"

./bin/gr -i ./tests/processes/gg2chic0.json -p "gg[chic(0)]<F> &> pi+ pi-" -w true -l false -n $EVENTS \
-o "chic0_MMHT2014lo68cl" -f "hepmc3" -q "MMHT2014lo68cl"

fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
# 0.05 (integrated screening) x 8E-3 (PDG chic0 -> pipi)
SCALE=0.0003

./bin/analyze \
-i "chic0_MSTW2008lo68cl, chic0_CT10nlo, chic0_MMHT2014lo68cl" \
-g "211, 211, 211" \
-n "2, 2, 2" \
-l "#chi_{c0} (MSTW2008lo68cl), #chi_{c0} (CT10nlo), #chi_{c0} (MMHT2014lo68cl)" \
-M "95, 3.2, 3.6" \
-Y "95,-1.0, 1.0" \
-P "95, 0.0, 3.0" \
-u pb \
-t '#sqrt{s} = 13 TeV, |#eta| < 0.9, p_{T} > 0.1 GeV | #chi_{c0}#rightarrow#pi^{+}#pi^{-} (BR~(2/3) #times 8E-3, S^{2}~0.05)' \
-S "$SCALE, $SCALE, $SCALE" #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

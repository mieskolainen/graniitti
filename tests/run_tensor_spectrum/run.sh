#!/bin/sh
#
# Simulation pi+pi-, K+K-, ppbar and analysis with Tensor Pomeron
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
-p "PP[RES+CONTENSOR]<F> -> pi+ pi- @RES{f0_500:1,rho_770:1,f0_980:1,f2_1270:1,f0_1500:1,f2_1525:1,f0_1710:1}" \
-w true -l false -n $EVENTS -o "tensor_pipi" -f "hepmc3"

./bin/gr -i ./tests/processes/CMSCUTS.json \
-p "PP[RES+CONTENSOR]<F> -> K+ K-   @RES{f0_980:1,phi_1020:1,f2_1270:1,f0_1500:1,f2_1525:1,f0_1710:1}" \
-w true -l false -n $EVENTS -o "tensor_KK" -f "hepmc3"

./bin/gr -i ./tests/processes/CMSCUTS.json \
-p "PP[RES+CONTENSOR]<F> -> p+ p- " \
-w true -l false -n $EVENTS -o "tensor_ppbar" -f "hepmc3"

fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.15

./bin/analyze \
-i "tensor_pipi,tensor_KK,tensor_ppbar" \
-g "211, 321, 2212" \
-n "2, 2, 2" \
-l "#pi^{+}#pi^{-} (TP), K^{+}K^{-} (TP), p#bar{p} (TP)" \
-M "95, 0.0, 2.5" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t "#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV" \
-S "$S2, $S2, $S2" #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

#!/bin/sh
#
# Simulation of rhorho, phiphi and with the same with tensor pomeron
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

N=100000

./bin/gr -i ./tests/processes/CMS13.json -p "PP[CON]<F> -> rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}" \
-w true -l false -n $N -o "rhorho" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json -p "PP[CON]<F> -> phi(1020)0 > {K+ K-} phi(1020)0 > {K+ K-}" \
-w true -l false -n $N -o "phiphi" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json -p "PP[CONTENSOR]<F> -> rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}" \
-w true -l false -n $N -o "rhorho_tensor" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json -p "PP[CONTENSOR]<F> -> phi(1020)0 > {K+ K-} phi(1020)0 > {K+ K-}" \
-w true -l false -n $N -o "phiphi_tensor" -f "hepmc3"

fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.15

./bin/analyze \
-i "rhorho, phiphi, rhorho_tensor, phiphi_tensor" \
-g "113, 333, 113, 333" \
-n "2, 2, 2, 2" \
-l "#rho^{0}#rho^{0} #rightarrow #pi^{+}#pi^{-}, #phi#phi #rightarrow K^{+}K^{-}, #rho^{0}#rho^{0} #rightarrow #pi^{+}#pi^{-} (TP), #phi#phi #rightarrow K^{+}K^{-} (TP)" \
-M "95, 0.0, 4.0" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 2.0" \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV' \
-S "$S2, $S2, $S2, $S2" #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

#!/bin/sh
#
# Simulation of kt-epa, DZ-epa and Durham flux only processes (matrix element = 1)
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Not set
if [ -z "$EVENTS" ]
then
	EVENTS=100000
fi

# Generate
./bin/gr -i ./tests/processes/fluxcuts.json -p "gg[FLUX]<F> -> 22 22" \
-w true -l false -n $EVENTS -o "pure_gg_flux" -f "hepmc3"

./bin/gr -i ./tests/processes/fluxcuts.json -p "yy[FLUX]<F> -> 22 22" \
-w true -l false -n $EVENTS -o "pure_yy_flux" -f "hepmc3"

./bin/gr -i ./tests/processes/fluxcuts.json -p "yy_DZ[FLUX]<P> -> 22 22" \
-w true -l false -n $EVENTS -o "pure_yy_dz_flux" -f "hepmc3"

fi

# Analyze
# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=1.0

./bin/analyze \
-i "pure_gg_flux, pure_yy_flux, pure_yy_dz_flux" \
-g "22, 22, 22" \
-n "2, 2, 2" \
-l "Durham flux, kt-EPA flux , DZ flux" \
-M "95, 50.0, 900.0" \
-Y "95,-2.5, 2.5" \
-P "95, 0.0, 3.0" \
-u nb \
-R false \
-t '|A|^{2} = 1.0, #sqrt{s} = 14 TeV, Y_{X} = [-2.5, 2.5], no screening' \
-S "$S2, $S2, $S2"
#-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

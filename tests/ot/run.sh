#!/bin/s
#
# Optimal Transport tests
#

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Generate
./bin/gr -i gencard/test.json -p "PP[CON]<C> -> pi+ pi-"     -o continuum1    -f hepmc3 -n 5800 -r 222
./bin/gr -i gencard/test.json -p "PP[CON]<C> -> pi+ pi-"     -o continuum2    -f hepmc3 -n 5900 -r 333
./bin/gr -i gencard/test.json -p "PP[RES+CON]<C> -> pi+ pi-" -o fullspectrum1 -f hepmc3 -n 6100 -r 555
./bin/gr -i gencard/test.json -p "PP[RES+CON]<C> -> pi+ pi-" -o fullspectrum2 -f hepmc3 -n 6000 -r 777
fi

# Analyze
LAMBDA=0.01
ITER=1500

#./bin/ot -i ./output/fullspectrum1.hepmc3,./output/fullspectrum1.hepmc3 -a $LAMBDA -r $ITER
#./bin/ot -i ./output/fullspectrum1.hepmc3,./output/fullspectrum2.hepmc3 -a $LAMBDA -r $ITER

./bin/ot -i ./output/fullspectrum1.hepmc3,./output/continuum1.hepmc3    -a $LAMBDA -r $ITER
#./bin/ot -i ./output/fullspectrum1.hepmc3,./output/continuum2.hepmc3    -a $LAMBDA -r $ITER


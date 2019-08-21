#!/bin/sh

# weighted events for speed

N=100000

./bin/gr -i ./tests/processes/CMS13.json -p "PP[RESTENSOR]<F> -> pi+ pi- @RES{f0_980:1}  @R[f0_980]{g0:1.0, g1:0.0}"  -w true -l false -n $N -o "f0_0" -f "hepmc3"
./bin/gr -i ./tests/processes/CMS13.json -p "PP[RESTENSOR]<F> -> pi+ pi- @RES{f0_980:1}  @R[f0_980]{g0:0.0, g1:1.0}"  -w true -l false -n $N -o "f0_1" -f "hepmc3"
./bin/gr -i ./tests/processes/CMS13.json -p "PP[RESTENSOR]<F> -> 22 22   @RES{eta:1}     @R[eta]{g0:1.0, g1:0.0}"     -w true -l false -n $N -o "eta_0" -f "hepmc3"
./bin/gr -i ./tests/processes/CMS13.json -p "PP[RESTENSOR]<F> -> 22 22-  @RES{eta:1}     @R[eta]{g0:0.0, g1:1.0}"     -w true -l false -n $N -o "eta_1" -f "hepmc3"
./bin/gr -i ./tests/processes/CMS13.json -p "yP[RESTENSOR]<F> -> pi+ pi- @RES{rho_770:1}" -w true -l false -n $N -o "rho" -f "hepmc3"

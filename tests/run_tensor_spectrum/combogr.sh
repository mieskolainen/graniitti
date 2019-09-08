#!/bin/sh

# weighted events for speed

N=100000

./bin/gr -i ./tests/processes/CMS13.json \
-p "PP[RES+CONTENSOR]<F> -> pi+ pi- @RES{f0_500:1,rho_770:1,f0_980:1,f2_1270:1,f0_1500:1,f2_1525:1,f0_1710:1}" \
-w true -l false -n $N -o "tensor_pipi" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json \
-p "PP[RES+CONTENSOR]<F> -> K+ K-   @RES{f0_980:1,phi_1020:1,f2_1270:1,f0_1500:1,f2_1525:1,f0_1710:1}" \
-w true -l false -n $N -o "tensor_KK" -f "hepmc3"

./bin/gr -i ./tests/processes/CMS13.json \
-p "PP[RES+CONTENSOR]<F> -> p+ p- " \
-w true -l false -n $N -o "tensor_ppbar" -f "hepmc3"

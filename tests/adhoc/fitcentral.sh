#!/bin/sh
#
# Soft central production fit example
#
# ** Note that if you fit several final state simultaneously,
#    all .json files should contain the same set of resonances under fitcard!
#
# Run with: source /tests/fitcentral.sh

make -j4 && ./bin/fitcentral \
-i 0 \ # input
-c 1 \ # continuum form
-a f0_500,rho_770,f0_980,phi_1020,f0_1500,f2_1525,f0_1710 \  # fix amplitude 
-p f0_500,rho_770,f0_980,phi_1020,f0_1500,f2_1525,f0_1710    # fix phase

# add 'continuum' to fix continuum parameters

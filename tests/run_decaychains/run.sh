#!/bin/sh
#
# Simulation of different decay chains
#
# Run with: source ./tests/run_xxx/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

# Not set
#if [ -z "$EVENTS" ]
#then
EVENTS=10
#fi

# Generate
./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[RES]<F> -> rho(770)0 > {pi0 > {22 22} 22} rho(770)0 > {pi0 pi+} @RES{f2_2200:1}" -w true -l false -n $EVENTS -h 0
./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[RES]<F> -> rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-} @RES{f2_2200:1}" -w true -l false -n $EVENTS -h 0
./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[RES]<F> -> rho(770)0 > {pi+ pi-} rho(770)0 @RES{f2_2200:1}" -w true -l false -n $EVENTS -h 0
./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[RES]<F> -> 22 22 @RES{eta:1}" -w true -l false -n $EVENTS -h 0
./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[RES]<F> -> 22 22 @RES{f0_980:1}" -w true -l false -n $EVENTS -h 0
./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[RES]<F> -> pi+ pi- @RES{f2_2200:1}" -w true -l false -n $EVENTS -h 0
./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[RES]<F> -> 22 22 @RES{f2_2200:1}" -w true -l false -n $EVENTS -h 0
./bin/gr -i ./tests/processes/ALICECUTS.json -p "PP[RES]<F> -> pi0 > {22 22} pi0 @RES{f2_2200:1}" -w true -l false -n $EVENTS -h 0


fi

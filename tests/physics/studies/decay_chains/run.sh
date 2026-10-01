#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize 10
#
# Simulation of different decay chains
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/decay_chains/run.sh

if ! study_generates
then
	study_skip "decay chain study has no separate analysis stage"
fi

if study_generates
then

# Generate
# Use &> for illustrative decays without physical decay input in DECAYS.json
study_gr -i ./tests/physics/studies/decay_chains/gencard_ALICECUTS.json -p "GP[RES]<F> &> pi0 > {22 22} pi0 > {22 22} @RES{f2_2150:1}" -h 0 -o decay_chain_1
study_gr -i ./tests/physics/studies/decay_chains/gencard_ALICECUTS.json -p "GP[RES]<F> -> rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-} @RES{f2_2150:1}" -h 0 -o decay_chain_2
study_gr -i ./tests/physics/studies/decay_chains/gencard_ALICECUTS.json -p "GP[RES]<F> -> rho(770)0 > {pi+ pi-} rho(770)0 @RES{f2_2150:1}" -h 0 -o decay_chain_3
study_gr -i ./tests/physics/studies/decay_chains/gencard_ALICECUTS.json -p "GP[RES]<F> -> 22 22 @RES{eta:1}" -h 0 -o decay_chain_4
study_gr -i ./tests/physics/studies/decay_chains/gencard_ALICECUTS.json -p "GP[RES]<F> &> 22 22 @RES{f0_980:1}" -h 0 -o decay_chain_5
study_gr -i ./tests/physics/studies/decay_chains/gencard_ALICECUTS.json -p "GP[RES]<F> -> pi+ pi- @RES{f2_2150:1}" -h 0 -o decay_chain_6
study_gr -i ./tests/physics/studies/decay_chains/gencard_ALICECUTS.json -p "GP[RES]<F> &> 22 22 @RES{f2_2150:1}" -h 0 -o decay_chain_7
study_gr -i ./tests/physics/studies/decay_chains/gencard_ALICECUTS.json -p "GP[RES]<F> &> pi0 > {22 22} pi0 @RES{f2_2150:1}" -h 0 -o decay_chain_8


fi

#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize 1000000
#
# Simulation and spherical harmonic expansion with ALICE cuts
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/alice_harmonic/run.sh

HARMONIC_INPUTS=(
	input/SH_2pi_J0_ALICE_fullsim.hepmc3
	input/SH_2pi_ALICE_fullsim.hepmc3
	input/ALICE_exclusive_pipi_data.hepmc3
)
study_require_opt_in RUN_HARMONIC "ALICE harmonic study"
study_require_files "ALICE harmonic study" "${HARMONIC_INPUTS[@]}"

if study_generates
then

# Generate
study_gr -i ./tests/physics/studies/alice_harmonic/gencard_SH_2pi_J0_ALICE.json
study_gr -i ./tests/physics/studies/alice_harmonic/gencard_SH_2pi_ALICE.json

fi

# Run the central only closure measurement
study_fitharmonic --card tests/physics/studies/alice_harmonic/central_measurement.json

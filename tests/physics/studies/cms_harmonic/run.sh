#!/usr/bin/env bash
set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize 1000000
#
# Simulation pi+pi- and spherical harmonic expansion with CMS cuts
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/cms_harmonic/run.sh

HARMONIC_INPUTS=(
	input/SH_2pi_J0_CMS_fullsim.hepmc3
	input/SH_2pi_CMS_fullsim.hepmc3
	input/CMS_exclusive_pipi_data.hepmc3
	input/SH_2pi_J0_CMS_tagged_fullsim.hepmc3
	input/SH_2pi_CMS_tagged_fullsim.hepmc3
	input/CMS_tagged_exclusive_pipi_data.hepmc3
)
study_require_opt_in RUN_HARMONIC "CMS harmonic study"
study_require_files "CMS harmonic study" "${HARMONIC_INPUTS[@]}"

if study_generates
then

# Generate
study_gr -i ./tests/physics/studies/cms_harmonic/gencard_SH_2pi_J0_CMS.json
study_gr -i ./tests/physics/studies/cms_harmonic/gencard_SH_2pi_CMS.json

fi

# Run central only and proton tagged closure measurements
study_fitharmonic --card tests/physics/studies/cms_harmonic/central_measurement.json
study_fitharmonic --card tests/physics/studies/cms_harmonic/tagged_measurement.json

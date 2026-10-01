#!/usr/bin/env bash
#
# Generate simulations and analysis (ROOT based) plots
#
# Run with: bash tests/physics/studies/run_plots.sh
#
# (c) 2017-2025 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

set -euo pipefail

source tests/physics/studies/common.sh
study_initialize 10000000

# Run one study with the shared action and event count
run_study() {
    STUDY_ACTION="${STUDY_ACTION}" NEVENTS="${NEVENTS}" \
      WEIGHTED="${WEIGHTED}" LOOPSCREEN="${LOOPSCREEN}" bash "$1"
}

cmake -S . -B build -DWITH_TEST=ON && cmake --build build -j4
run_study ./tests/physics/studies/atlas_multi/run.sh

run_study ./tests/physics/studies/excitation/run.sh

run_study ./tests/physics/studies/screening/run.sh
run_study ./tests/physics/studies/alice_multi/run.sh
run_study ./tests/physics/studies/jw_polarization/run.sh
run_study ./tests/physics/studies/jw_frames/run.sh


# Tensor Pomeron
#run_study ./tests/physics/studies/tensor0/run.sh
run_study ./tests/physics/studies/tensor2/run.sh
#run_study ./tests/physics/studies/tensor_spectrum/run.sh


# Spherical harmonic expansion requires independently propagated detector and data inputs
if [[ "${RUN_HARMONIC:-0}" == "1" ]]
then
	run_study ./tests/physics/studies/cms_harmonic/run.sh
	#run_study ./tests/physics/studies/alice_harmonic/run.sh
fi

echo "[run_plots.sh: done]"

if [[ "${COPY_PLOTS:-0}" == "1" ]]
then
	bash ./tests/physics/studies/copy_plots.sh
fi

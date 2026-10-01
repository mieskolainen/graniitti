#!/usr/bin/env bash
#
# Run all physics studies with a run.sh launcher
#
# Run with: bash ./tests/physics/studies/run_all.sh > run_all.out
#
# (c) 2017-2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>

set -euo pipefail

source tests/physics/studies/common.sh
study_initialize 5000

# Run one discovered study with the selected action
run_study() {
    STUDY_ACTION="${STUDY_ACTION}" NEVENTS="${NEVENTS}" \
      WEIGHTED="${WEIGHTED}" LOOPSCREEN="${LOOPSCREEN}" bash "$1"
}

cmake -S . -B build -DWITH_TEST=ON && cmake --build build -j4

# Discover every canonical study launcher recursively
mapfile -t STUDIES < <(find tests/physics/studies -mindepth 2 -type f -name run.sh -print | sort)
for STUDY in "${STUDIES[@]}"
do
	run_study "${STUDY}"
done

echo "[run_all.sh: done]"

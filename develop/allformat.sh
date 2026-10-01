#!/usr/bin/env bash
# Format repository C++ sources except generated MG5 amplitudes
#
# Run with: bash develop/allformat.sh

set -euo pipefail

mapfile -d '' files < <(find ./src ./include/Graniitti -type f \
    \( -name '*.cc' -o -name '*.h' \) \
    ! -path '*/Amplitude/MG5/*' ! -path '*._old*' -print0 | sort -z)
wait "$!"
bash ./develop/format.sh "${files[@]}"

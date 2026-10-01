#!/usr/bin/env bash
# Print the build configuration, MG5 converter and C++ sources as PDFs
#
# Run with: bash develop/code2pdf.sh

set -euo pipefail

options=(--verbose --color=1 -C -fCourier8 -o - '--header=$N|$%/$=|')
enscript "${options[@]}" CMakeLists.txt | ps2pdf - code_cmake.pdf

mapfile -d '' converter < <(find ./develop/MG2GRA -name '*.py' -type f ! -path '*._old*' -print0 | sort -z)
wait "$!"
enscript -Epython "${options[@]}" "${converter[@]}" | ps2pdf - code_MG2GRA.pdf

mapfile -d '' sources < <(find ./src -name '*.cc' -type f ! -path '*._old*' -print0 | sort -z)
wait "$!"
enscript -Ecpp "${options[@]}" "${sources[@]}" | ps2pdf - code_cc.pdf

mapfile -d '' headers < <(find ./include -name '*.h' -type f ! -path '*._old*' -print0 | sort -z)
wait "$!"
enscript -Ecpp "${options[@]}" "${headers[@]}" | ps2pdf - code_h.pdf

#!/usr/bin/env bash
#
# Central production at different energies, screening / no screening, excitation (El,SD,DD)
#
# Run with: NEVENTS=10000 bash ./tests/physics/studies/cross_section_scan/run.sh

set -euo pipefail
STUDY_LAUNCHER_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${STUDY_LAUNCHER_DIR}/../common.sh"
study_initialize

P=./tests/physics/studies/cross_section_scan
A="$P/analysis"

if study_generates
then

E=1.995262E+01,7.170601E+01,2.576980E+02,9.261187E+02,3.328298E+03,1.196128E+04,4.298662E+04,1.544859E+05,5.551936E+05,1.995262E+06
E="${ENERGIES:-${E}}"

./bin/xscan -i $P/gencard_CEP_EL.json,$P/gencard_CEP_SD.json,$P/gencard_CEP_DD.json -e $E -l false
cp scan.csv scan_screening_false.csv
./bin/xscan -i $P/gencard_CEP_EL.json,$P/gencard_CEP_SD.json,$P/gencard_CEP_DD.json -e $E -l true
cp scan.csv scan_screening_true.csv

# Retain both scan tables for analysis without repeating integration
cp scan_screening_false.csv scan_screening_true.csv "$A/"
else
    study_require_outputs "cross section analysis" "$A/scan_screening_false.csv" "$A/scan_screening_true.csv"
fi

python "$A/analyze.py"

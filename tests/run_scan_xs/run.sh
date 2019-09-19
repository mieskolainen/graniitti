#!/bin/sh
#
# Central production at different energies, screening / no screening, excitation (El,SD,DD)
#
# Run with: source ./tests/run_scan_xs/run.sh

read -p "run: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then

P=./tests/run_scan_xs

#E=62,500,546,1800,7000,8000,13000,60000
E=1.995262E+01,7.170601E+01,2.576980E+02,9.261187E+02,3.328298E+03,1.196128E+04,4.298662E+04,1.544859E+05,5.551936E+05,1.995262E+06

./bin/xscan -i $P/CEP_EL.json,$P/CEP_SD.json,$P/CEP_DD.json -e $E -l false
cp scan.csv scan_screening_false.csv
./bin/xscan -i $P/CEP_EL.json,$P/CEP_SD.json,$P/CEP_DD.json -e $E -l true
cp scan.csv scan_screening_true.csv

# Copy
cp scan_screening_* ./tests/analysis_scan_xs/

fi

#!/bin/bash
#
# LHAPDF 6.x compilation and install script
#
# Run first: source setenv.sh
# 
# mikael.mieskolainen@cern.ch (2026)

set -euo pipefail

#wget http://www.hepforge.org/archive/lhapdf/LHAPDF-6.5.5.tar.gz
tar -xf LHAPDF-6.5.5.tar.gz
cd LHAPDF-6.5.5

# Remove old
rm ${GRANIITTI_IO_PATH}/LHAPDF -f -r

# Compile and install new
./configure --prefix=${GRANIITTI_IO_PATH}/LHAPDF
make -j4
make install
cd ..
sleep 3
rm LHAPDF-6.5.5 -f -r

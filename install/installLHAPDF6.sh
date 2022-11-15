#!/bin/bash
#
# LHAPDF 6.x compilation and install script
#
# Ubuntu requirements:
#   sudo apt install cmake g++ python3-dev
#
# Run with:
#
#   INSTALLPATH=$HOME/local
#   source installLHAPDF6.sh
#

#wget http://www.hepforge.org/archive/lhapdf/LHAPDF-6.5.3.tar.gz
tar -xf LHAPDF-6.5.3.tar.gz
cd LHAPDF-6.5.3

# Remove old
rm ${INSTALLPATH}/LHAPDF -f -r

# Compile and install new
./configure --prefix=${INSTALLPATH}/LHAPDF
make -j4
make install
cd ..
sleep 3
rm LHAPDF-6.5.3 -f -r


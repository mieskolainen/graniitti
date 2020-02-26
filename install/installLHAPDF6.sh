#!/bin/bash
#
# LHAPDF 6.x compilation and install script
#
# Ubuntu requirements:
# sudo apt-get python-dev
#
# Run with:
# INSTALLPATH=$HOME/local
# source installLHAPDF6.sh
#
#wget http://www.hepforge.org/archive/lhapdf/LHAPDF-6.2.3.tar.gz

tar -xf LHAPDF-6.2.3.tar.gz
cd LHAPDF-6.2.3
./configure --prefix=${INSTALLPATH}/LHAPDF
make -j4
make install
cd ..
wait 2
rm LHAPDF-6.2.3 -f -r

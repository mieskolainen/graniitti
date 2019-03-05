#!/bin/bash
#
# Download and install some LHAPDF 6.x PDFSETS
#
# Run with:
# INSTALLPATH=$HOME/local
# source instalLPDFSET.sh

array=(CT10nlo MMHT2014lo68cl MSTW2008lo68cl)

for i in "${array[@]}"
do
	wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.2/$i.tar.gz -O- | tar xz -C ${INSTALLPATH}/LHAPDF/share/LHAPDF
done

#!/bin/bash
#
# Download and install some LHAPDF 6.x PDFSETS
#
# Run first: source setenv.sh
# 
# mikael.mieskolainen@cern.ch (2026)

set -euo pipefail

array=(CT10nlo MMHT2014lo68cl MSHT20lo_as130 MSTW2008lo68cl NNPDF31_lo_as_0118 GKG18_DPDF_FitB_LO GKG18_DPDF_FitB_NLO LUXqed17_plus_PDF4LHC15_nnlo_100)
#array=(LUXqed17_plus_PDF4LHC15_nnlo_100)

for i in "${array[@]}"
do
  wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/$i.tar.gz
  tar -xf $i.tar.gz -C ${GRANIITTI_IO_PATH}/LHAPDF/share/LHAPDF
  rm $i.tar.gz
done

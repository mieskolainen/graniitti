#!/bin/sh

make -j4 && ./bin/analyze -i CMS_13_2pi,CMS_13_2K -g 211,321 -n 2,2 -l \
'#pi^{+}#pi^{-}','#it{K}^{+}#it{K}^{-} / 13 TeV' \
-M 3.0 -Y 2.5 -P 2.0 -u ub -S 0.18 #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

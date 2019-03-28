#!/bin/sh

# Integrated screening factor (for speed, set 1 if Pomeron loop is on)
S2=0.18

make -j4 && ./bin/analyze -i CMS19_2pi,CMS19_2K -g 211,321 -n 2,2 -l \
'#pi^{+}#pi^{-}','#it{K}^{+}#it{K}^{-} / 13 TeV' \
-M 3.0 -Y 2.5 -P 2.0 -u ub -S $S2 #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

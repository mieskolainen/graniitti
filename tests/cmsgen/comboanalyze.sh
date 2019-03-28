#!/bin/sh

# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.18

make -j4 && ./bin/analyze -i CMS19_2pi,CMS19_2K -g 211,321 -n 2,2 -l \
'#pi^{+}#pi^{-}','#it{K}^{+}#it{K}^{-}' \
-M 3.0 -Y 2.5 -P 2.0 -u ub -t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.2 GeV' -S $S2 #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

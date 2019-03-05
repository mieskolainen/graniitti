#!/bin/sh

make -j4 && ./bin/analyze -i ALICE2pi,ALICE2K,ALICEppbar -g 211,321,2212 -n 2,2,2 -l \
'#pi^{+}#pi^{-}','#it{K}^{+}#it{K}^{-}','#it{p}#bar{#it{p}} / 13 TeV ' \
-M 3.0 -Y 1.5 -P 2.0 -u ub -S 0.18 #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

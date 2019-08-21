#!/bin/sh

# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.18

./bin/analyze \
-i tensor_pipi,tensor_KK \
-g 211,321 \
-n 2,2 \
-l '#pi^{+}#pi^{-}','K^{+}K^{-}' \
-M 95,0.0,2.5 \
-Y 95,-2.5,2.5 \
-P 95,0.0,2.0 \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV' \
-S $S2,$S2 #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

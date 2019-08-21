#!/bin/sh

# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.18

./bin/analyze \
-i f0_0,f0_1,eta_0,eta_1 \
-g 211,211,22,22 \
-n 2,2,2,2 \
-l 'f0: <02>','f0: <22>','eta: <11>','eta: <33>' \
-M 95,0.3,1.3 \
-Y 95,-2.5,2.5 \
-P 95,0.0,2.0 \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV' \
-S $S2,$S2,$S2,$S2 #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

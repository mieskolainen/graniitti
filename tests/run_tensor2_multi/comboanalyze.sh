#!/bin/sh

# Hard-coded integrated screening factor (for speed, set 1 if Pomeron loop was on)
S2=0.18

./bin/analyze \
-i f2_0,f2_1,f2_2,f2_3,f2_4,f2_5,f2_6 \
-g 211,211,211,211,211,211,211 \
-n 2,2,2,2,2,2,2 \
-l 'f2: <02>','f2: <20> - <22>','f2: <20> + <22>','f2: <24>','f2: <42>','f2: <44>','f2: <64>' \
-M 95,0.75,1.75 \
-Y 95,-2.5,2.5 \
-P 95,0.0,2.0 \
-u ub \
-t '#sqrt{s} = 13 TeV, |#eta| < 2.5, p_{T} > 0.15 GeV' \
-S $S2,$S2,$S2,$S2,$S2,$S2,$S2 #-X 1000

#convert -density 600 -trim .h1_S_M_logy.pdf -quality 100 output.jpg

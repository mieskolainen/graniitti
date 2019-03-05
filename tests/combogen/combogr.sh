#!/bin/sh

make -j4 \
&& ./bin/gr -i ./tests/experiment/ALICE_2pi.json -w true -l false -n 4000000 \
&& ./bin/gr -i ./tests/experiment/ALICE_2K.json -w true -l false -n 320000 \
&& ./bin/gr -i ./tests/experiment/ALICE_ppbar.json -w true -l false -n 32000


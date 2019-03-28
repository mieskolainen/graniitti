#!/bin/sh

make -j4 \
&& ./bin/gr -i ./tests/processes/ALICE_2pi.json -w true -l false -n 40000 \
&& ./bin/gr -i ./tests/processes/ALICE_2K.json -w true -l false -n 3200 \
&& ./bin/gr -i ./tests/processes/ALICE_ppbar.json -w true -l false -n 320


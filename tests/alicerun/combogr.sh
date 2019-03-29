#!/bin/sh

# weighted events for speed

make -j4 \
&& ./bin/gr -i ./tests/processes/ALICE_2pi.json -w true -l false -n 1000000 \
&& ./bin/gr -i ./tests/processes/ALICE_2K.json -w true -l false -n 50000 \
&& ./bin/gr -i ./tests/processes/ALICE_ppbar.json -w true -l false -n 50000


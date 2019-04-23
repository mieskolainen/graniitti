#!/bin/sh

# weighted events for speed

make -j4 \
&& ./bin/gr -i ./tests/processes/CMS19_2pi.json -w true -l false -n 100000 \
&& ./bin/gr -i ./tests/processes/CMS19_2K.json -w true -l false -n 100000 \
&& ./bin/gr -i ./tests/processes/CMS19_ppbar.json -w true -l false -n 100000

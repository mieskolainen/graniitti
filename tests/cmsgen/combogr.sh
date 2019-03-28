#!/bin/sh

make -j4 \
&& ./bin/gr -i ./tests/processes/CMS19_2pi.json -w true -l false -n 40000 \
&& ./bin/gr -i ./tests/processes/CMS19_2K.json -w true -l false -n 3200


#!/bin/sh
#
# Run with: source all.sh

yes yy_ll      | python MG2GRA.py
yes yy_ww      | python MG2GRA.py
yes yy_ww_evev | python MG2GRA.py
yes gg_gg      | python MG2GRA.py
yes gg_ggg     | python MG2GRA.py

#!/bin/sh
#
# Soft central production fit example
#
# Requires python tools (see FAQ / README)

ray start --head
python python/icetune --tuneset default

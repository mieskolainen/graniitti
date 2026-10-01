#!/bin/sh
# build on lxplus (-DMARCH=x86-64-v3)

cmake -S . -B build -DMARCH=x86-64-v3
cmake --build build -j8

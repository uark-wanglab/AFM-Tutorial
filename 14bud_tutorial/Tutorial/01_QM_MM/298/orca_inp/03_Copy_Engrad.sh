#!/usr/bin/bash

CWD=$(pwd)

rm final*/*densities
rm final*/*gbw

for item in "$@"; do
    cp -a $item/${item}.orca.engrad ../genref
done

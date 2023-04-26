#!/bin/bash

CWD=$(pwd)


for item in "$@"; do
    cd $item && sbatch runme.pinnacle ${item}.orca.inp
    cd $CWD
done
